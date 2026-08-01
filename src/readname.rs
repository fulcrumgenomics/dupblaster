//! Locating the sequencing unit and imaging tile within a read name.
//!
//! Duplicates that arise **on the flowcell** — ExAmp/cluster duplicates — are
//! copies of one template imaged in one place, so they share a tile. Duplicates
//! that arise **in the library** (PCR copies, or two genuinely distinct
//! molecules at one locus) are independent and so land on tiles independently.
//! That asymmetry is the entire basis of the sequencing-vs-library
//! decomposition, and the read name is the only place the input records it.
//!
//! Both tokens are treated as **opaque bytes and never parsed as numbers**. The
//! decomposition needs identity, not geometry: two templates are on one tile iff
//! their tile tokens are equal. That keeps extraction immune to instruments
//! widening tile numbers or changing their encoding, and it is why no pixel
//! distance appears anywhere in this module — fixed pixel radii do not
//! generalize across runs, while tile identity does.
//!
//! **The format is chosen by the user, never guessed.** Read-name layouts differ
//! between platforms in ways that can silently mis-parse rather than fail (the
//! pre-CASAVA-1.8 Illumina layout puts a y coordinate where the modern one puts
//! a tile), so auto-detection would risk producing a confident wrong number. A
//! name the chosen format cannot parse is an error, not something to work around.

use std::str::FromStr;

use anyhow::{Context, Result, anyhow, bail};
use memchr::memchr_iter;
use regex::bytes::Regex;

/// Colons bounding the fields a [`ReadNameFormat::ColonDelimited`] name needs.
///
/// Fields are 0-based and field `i` is bounded by colons `i-1` and `i`, so
/// flowcell/lane/tile (2, 3, 4) need colons 1 through 4. Requiring a *sixth*
/// colon on top of those enforces the seven-field minimum that separates this
/// layout from the five-field pre-CASAVA-1.8 one, whose field 4 is a y
/// coordinate rather than a tile.
const REQUIRED_COLONS: usize = 6;

/// Spelling that introduces a user-supplied regex in a format specification.
const REGEX_PREFIX: &str = "regex:";

/// Where a template was imaged, as opaque tokens borrowed from its read name.
///
/// The unit is kept separate from the tile because **tiles are only comparable
/// within a unit**: tile `1101` of one flowcell and tile `1101` of another are
/// different physical places, and conflating them is exactly the bug that makes
/// `samtools markdup` call cross-flowcell duplicates optical (samtools#1996) —
/// measured at 1.98% of duplicates on real data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ImagingLocation<'a> {
    /// The sequencing unit: what identifies one independently-loaded surface,
    /// i.e. flowcell plus lane.
    pub unit: &'a [u8],
    /// The imaging tile within `unit`.
    pub tile: &'a [u8],
}

/// How to locate the sequencing unit and tile within a read name.
///
/// Parsed from a user-supplied specification: `illumina`, `element`, or
/// `regex:<pattern>` — see [`FromStr`].
#[derive(Clone, Debug)]
pub enum ReadNameFormat {
    /// Seven or more colon-delimited fields —
    /// `instrument:run:flowcell:lane:tile:x:y` — with the unit taken as the
    /// contiguous `flowcell:lane` span and the tile as field 4 (0-based).
    ///
    /// Illumina (CASAVA 1.8+, bcl2fastq, BCL Convert) and Element AVITI
    /// (`P2-05:Dual-Index-Sim:FC-63067cc:1:10201:0000:0015`) share this layout
    /// exactly, so they are one extractor reachable under either spelling.
    /// Singular G4 documents lane/tile/X/Y and is very likely the same, but is
    /// unverified against a real file. Fields beyond the seventh are ignored, so
    /// names carrying an appended UMI still parse.
    ColonDelimited,
    /// A user-supplied regex, for platforms without a preset. The escape hatch
    /// for undelimited layouts such as MGI DNBSEQ
    /// (`F350009384L1C001R0010008170`), which no field-splitting rule can
    /// handle. Precedent: Picard's `READ_NAME_REGEX`.
    Custom(CustomReadNameRegex),
}

impl Default for ReadNameFormat {
    /// [`Self::ColonDelimited`] — the layout of every Illumina instrument since
    /// CASAVA 1.8, and of Element AVITI, so the overwhelmingly common case needs
    /// no flag. A read name this layout cannot parse fails the run: see
    /// [`Self::parse_error`], and `--sequencing-duplicate-detection off` for the way out.
    fn default() -> Self {
        Self::ColonDelimited
    }
}

impl ReadNameFormat {
    /// Extract the imaging location from `name`, or `None` if `name` does not
    /// fit this format. A `None` here is a hard error at the call site: the user
    /// named the format, so a name that does not match means the wrong format
    /// was chosen or the data is not what it claims to be.
    #[inline]
    pub(crate) fn extract<'a>(&self, name: &'a [u8]) -> Option<ImagingLocation<'a>> {
        match self {
            Self::ColonDelimited => colon_delimited_location(name),
            Self::Custom(custom) => custom.extract(name),
        }
    }

    /// Build the error for a read name this format could not parse, quoting the
    /// name back alongside what the format expected.
    pub(crate) fn parse_error(&self, name: &[u8]) -> anyhow::Error {
        anyhow!(
            "read name {:?} does not match the {} read-name format",
            String::from_utf8_lossy(name),
            self.describe()
        )
    }

    /// How to refer to this format in a message, in the user's own terms.
    fn describe(&self) -> String {
        match self {
            Self::ColonDelimited => {
                "illumina/element (instrument:run:flowcell:lane:tile:x:y)".to_string()
            }
            Self::Custom(custom) => format!("{REGEX_PREFIX}{}", custom.regex.as_str()),
        }
    }
}

impl FromStr for ReadNameFormat {
    type Err = anyhow::Error;

    fn from_str(spec: &str) -> Result<Self> {
        match spec {
            "illumina" | "element" => Ok(Self::ColonDelimited),
            _ => match spec.strip_prefix(REGEX_PREFIX) {
                Some(pattern) => CustomReadNameRegex::new(pattern).map(Self::Custom),
                None => bail!(
                    "unknown read-name format {spec:?}; expected `illumina`, `element`, or \
                     `{REGEX_PREFIX}<pattern>`"
                ),
            },
        }
    }
}

/// A validated user-supplied read-name regex with its capture groups resolved.
///
/// The group indices are resolved once at construction so a missing `su` or
/// `tile` group is reported when the option is parsed rather than on the first
/// read, and so matching costs no name lookups per template.
#[derive(Clone, Debug)]
pub struct CustomReadNameRegex {
    /// Matched against the raw read-name bytes, so the captured tokens borrow
    /// from the name rather than being copied through a `str`.
    regex: Regex,
    /// Capture-group index of the sequencing unit.
    unit_group: usize,
    /// Capture-group index of the tile.
    tile_group: usize,
}

impl CustomReadNameRegex {
    /// Compile `pattern`, requiring both named capture groups.
    fn new(pattern: &str) -> Result<Self> {
        let regex =
            Regex::new(pattern).with_context(|| format!("invalid read-name regex {pattern:?}"))?;
        let unit_group = named_group(&regex, "su")?;
        let tile_group = named_group(&regex, "tile")?;
        Ok(Self { regex, unit_group, tile_group })
    }

    /// Match `name` and return the captured tokens, which borrow from `name`.
    #[inline]
    fn extract<'a>(&self, name: &'a [u8]) -> Option<ImagingLocation<'a>> {
        let captures = self.regex.captures(name)?;
        let unit = captures.get(self.unit_group)?.as_bytes();
        let tile = captures.get(self.tile_group)?.as_bytes();
        if unit.is_empty() || tile.is_empty() {
            return None;
        }
        Some(ImagingLocation { unit, tile })
    }
}

/// Locate the `flowcell:lane` unit and the tile in a colon-delimited read name.
///
/// The unit is returned as the single contiguous `flowcell:lane` span rather
/// than two slices, so it interns as one token and prints as one readable
/// column in the per-sequencing-unit report.
fn colon_delimited_location(name: &[u8]) -> Option<ImagingLocation<'_>> {
    let mut colon = [0usize; REQUIRED_COLONS];
    let mut found = 0;
    for pos in memchr_iter(b':', name) {
        colon[found] = pos;
        found += 1;
        if found == REQUIRED_COLONS {
            break;
        }
    }
    if found < REQUIRED_COLONS {
        return None;
    }

    // Reject an empty flowcell, lane, or tile. The tokens are compared only for
    // equality, so an empty one would silently merge places that are distinct.
    let flowcell = name.get(colon[1] + 1..colon[2])?;
    let lane = name.get(colon[2] + 1..colon[3])?;
    let tile = name.get(colon[3] + 1..colon[4])?;
    if flowcell.is_empty() || lane.is_empty() || tile.is_empty() {
        return None;
    }

    Some(ImagingLocation { unit: name.get(colon[1] + 1..colon[3])?, tile })
}

/// Index of the capture group named `name` in `regex`.
fn named_group(regex: &Regex, name: &str) -> Result<usize> {
    regex.capture_names().position(|group| group == Some(name)).ok_or_else(|| {
        anyhow!("read-name regex {:?} has no (?<{name}>...) capture group", regex.as_str())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A CASAVA 1.8+ Illumina read name from the liquidcell test sample.
    const ILLUMINA: &[u8] = b"A00354:1305:H72CFDSXF:2:1101:1027:1986";
    /// An Element AVITI read name — the same seven-field layout.
    const ELEMENT: &[u8] = b"P2-05:Dual-Index-Sim:FC-63067cc:1:10201:0000:0015";

    fn extract(name: &[u8]) -> Option<ImagingLocation<'_>> {
        ReadNameFormat::ColonDelimited.extract(name)
    }

    fn format(spec: &str) -> ReadNameFormat {
        spec.parse().expect("spec parses")
    }

    #[test]
    fn illumina_name_yields_flowcell_lane_unit_and_tile() {
        let loc = extract(ILLUMINA).expect("parses");
        assert_eq!(loc.unit, b"H72CFDSXF:2");
        assert_eq!(loc.tile, b"1101");
    }

    #[test]
    fn element_aviti_name_parses_with_the_same_layout() {
        let loc = extract(ELEMENT).expect("parses");
        assert_eq!(loc.unit, b"FC-63067cc:1");
        assert_eq!(loc.tile, b"10201");
    }

    #[test]
    fn trailing_fields_beyond_the_seventh_are_ignored() {
        let loc = extract(b"A00354:1305:H72CFDSXF:2:1101:1027:1986:ACGTACGT").expect("parses");
        assert_eq!(loc.unit, b"H72CFDSXF:2");
        assert_eq!(loc.tile, b"1101");
    }

    #[test]
    fn mate_suffix_does_not_disturb_the_tile() {
        let loc = extract(b"A00354:1305:H72CFDSXF:2:1101:1027:1986/1").expect("parses");
        assert_eq!(loc.tile, b"1101");
    }

    #[test]
    fn legacy_five_field_illumina_name_is_rejected() {
        // Field 4 here is the y coordinate, not a tile; parsing it as one would
        // give nearly every read its own "tile" and silently destroy the estimate.
        assert!(extract(b"HWUSI-EAS100R:6:73:941:1973#0/1").is_none());
    }

    #[test]
    fn name_with_no_colons_is_rejected() {
        assert!(extract(b"SRR1234567.1").is_none());
    }

    #[test]
    fn name_with_empty_tile_field_is_rejected() {
        assert!(extract(b"A00354:1305:H72CFDSXF:2::1027:1986").is_none());
    }

    #[test]
    fn name_with_empty_flowcell_field_is_rejected() {
        assert!(extract(b"A00354:1305::2:1101:1027:1986").is_none());
    }

    #[test]
    fn name_with_empty_lane_field_is_rejected() {
        assert!(extract(b"A00354:1305:H72CFDSXF::1101:1027:1986").is_none());
    }

    #[test]
    fn tile_tokens_are_compared_verbatim_without_numeric_normalization() {
        let padded = extract(b"A:1:FC:1:0001:0:0").expect("parses");
        let bare = extract(b"A:1:FC:1:1:0:0").expect("parses");
        assert_ne!(padded.tile, bare.tile);
    }

    #[test]
    fn same_tile_number_on_two_flowcells_yields_different_units() {
        // The samtools#1996 failure mode: matching tile *numbers* across
        // flowcells must not read as the same physical tile.
        let a = extract(b"A00354:1305:H72CFDSXF:2:1101:1027:1986").expect("parses");
        let b = extract(b"A00354:1305:22T3L2LT4:2:1101:1027:1986").expect("parses");
        assert_eq!(a.tile, b.tile);
        assert_ne!(a.unit, b.unit);
    }

    #[test]
    fn different_lanes_of_one_flowcell_are_different_units() {
        let a = extract(b"A00354:1305:H72CFDSXF:1:1101:1027:1986").expect("parses");
        let b = extract(b"A00354:1305:H72CFDSXF:2:1101:1027:1986").expect("parses");
        assert_ne!(a.unit, b.unit);
    }

    #[test]
    fn illumina_and_element_specs_select_the_same_extractor() {
        for spec in ["illumina", "element"] {
            let loc = format(spec).extract(ILLUMINA).expect("parses");
            assert_eq!(loc.unit, b"H72CFDSXF:2");
        }
    }

    #[test]
    fn unknown_format_spec_is_rejected_and_lists_the_alternatives() {
        let err = "novaseq".parse::<ReadNameFormat>().expect_err("must be rejected").to_string();
        assert!(err.contains("illumina"), "{err}");
        assert!(err.contains("regex:"), "{err}");
    }

    #[test]
    fn custom_regex_extracts_tokens_from_an_undelimited_name() {
        // MGI DNBSEQ packs its fields with no delimiter, so no field-splitting
        // rule can reach them; a regex can. Fixed widths here rather than `\d+`
        // both to pin the tokens exactly and because MGI's true field widths are
        // unverified against a real file — the point is the mechanism.
        let mgi = format(r"regex:^(?<su>F\w+L\d)(?<tile>C\d{3}R\d{3})");
        let loc = mgi.extract(b"F350009384L1C001R0010008170").expect("parses");
        assert_eq!(loc.unit, b"F350009384L1");
        assert_eq!(loc.tile, b"C001R001");
    }

    #[test]
    fn custom_regex_need_not_be_anchored_to_the_whole_name() {
        // `captures` searches rather than requiring a full match, so a pattern
        // may pick tokens out of the middle of a name.
        let loc = format(r"regex:lane(?<su>\d+)_tile(?<tile>\d+)")
            .extract(b"run7_lane3_tile21_x0_y0")
            .expect("parses");
        assert_eq!(loc.unit, b"3");
        assert_eq!(loc.tile, b"21");
    }

    #[test]
    fn custom_regex_without_an_su_group_is_rejected() {
        let err = r"regex:^(?<tile>\d+)".parse::<ReadNameFormat>().expect_err("must be rejected");
        assert!(err.to_string().contains("(?<su>"), "{err}");
    }

    #[test]
    fn custom_regex_without_a_tile_group_is_rejected() {
        let err = r"regex:^(?<su>\d+)".parse::<ReadNameFormat>().expect_err("must be rejected");
        assert!(err.to_string().contains("(?<tile>"), "{err}");
    }

    #[test]
    fn syntactically_invalid_custom_regex_is_rejected() {
        assert!(r"regex:^(?<su>[".parse::<ReadNameFormat>().is_err());
    }

    #[test]
    fn custom_regex_that_does_not_match_yields_no_location() {
        let mgi = format(r"regex:^(?<su>F\w+L\d)(?<tile>C\d+R\d+)");
        assert!(mgi.extract(ILLUMINA).is_none());
    }

    #[test]
    fn parse_error_names_the_read_name_and_the_chosen_format() {
        let err = format("illumina").parse_error(b"SRR1234567.1").to_string();
        assert!(err.contains("SRR1234567.1"), "{err}");
        assert!(err.contains("instrument:run:flowcell:lane:tile:x:y"), "{err}");
    }

    #[test]
    fn parse_error_for_a_custom_format_quotes_the_pattern() {
        let err = format(r"regex:^(?<su>F\w+L\d)(?<tile>C\d+R\d+)")
            .parse_error(b"SRR1234567.1")
            .to_string();
        assert!(err.contains("(?<su>"), "{err}");
    }
}
