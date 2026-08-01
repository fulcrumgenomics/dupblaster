//! Run-summary metrics emitted as a per-library TSV (one row per library) when
//! `--metrics-prefix PREFIX` is set on the command line.
//!
//! Rows are plain `#[derive(Serialize)]` structs written through
//! [`fgoxide::io::DelimFile`] (the same serde-driven path riker uses), so the
//! struct fields are the single source of truth for the columns — there is no
//! separate column list to keep in sync. Fixed-precision floats use a
//! `serialize_with` helper.
//!
//! Column layout is informed by Picard's [`DuplicationMetrics`][picard] but
//! adapted to our template-level data model and Riker-style `frac_` naming.
//! Library-size estimation is the standard Lander-Waterman bisection ported
//! from Picard's `estimateLibrarySize` (40 bisection steps, expanding upper
//! bound until bracketed).
//!
//! [picard]: https://github.com/broadinstitute/picard/blob/main/src/main/java/picard/sam/DuplicationMetrics.java

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;
use serde::Serialize;

use crate::DUPBLASTER_BUILD;
use crate::dedup::{LibraryStats, Stats};
use crate::tiles::{Decomposition, SequencingUnitStats};

/// Summary metrics for one library within a dupblaster run.
///
/// One row in the emitted TSV. In single-library mode there is one row; in
/// library-aware mode (`>1` distinct `@RG LB:`) there is one row per library.
/// Optional numeric fields (currently just `estimated_library_size`) render as
/// an empty cell when `None` (serde serializes `Option::None` to an empty
/// field). Field declaration order is the TSV column order.
#[derive(Debug, Clone, Serialize)]
pub struct Metrics {
    /// Sample name — `--sample` if set, else comma-joined unique `@RG SM:` values.
    pub sample: String,
    /// Library name for this row — an `@RG LB:` value, `"Unknown Library"`, or
    /// `"All Reads"` when library splitting is off.
    pub library: String,
    /// `CARGO_PKG_VERSION` at build time.
    pub dupblaster_version: &'static str,
    /// Total templates processed.
    pub total_templates: u64,
    /// Templates flagged as duplicate (rollup of pair + orphan dups).
    pub duplicate_templates: u64,
    /// Read-level duplicate fraction in [0, 1], matching Picard's
    /// `PERCENT_DUPLICATION` formula:
    /// `(orphan_dups + 2*pair_dups) / (orphan_reads + 2*pair_reads)`.
    /// Emitted with 6 decimal places.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_duplicates: f64,
    /// Templates where both reads of the pair are mapped.
    pub mapped_pairs: u64,
    /// Templates where both reads are unmapped — never dup-checked.
    pub unmapped_pairs: u64,
    /// Duplicates among `mapped_pairs`.
    pub duplicate_pairs: u64,
    /// Duplicate pairs that are copies made on the flowcell (cluster/ExAmp
    /// duplicates), counted as `Σ (k − tiles)` over duplicate groups: a tile
    /// holding `c` members of a group seeded one molecule and copied the rest.
    ///
    /// Uncorrected, so it slightly over-states the count wherever a group's
    /// molecules collide on a tile by chance — negligible on WGS, material for
    /// groups of hundreds. The companion `<PREFIX>.sequencing-units.tsv` sums to
    /// exactly this figure.
    ///
    /// `None` — an empty cell — whenever the split was not computed:
    /// `--sequencing-duplicate-detection off` was passed, the library has no both-ends-mapped
    /// pairs, or it sits on a single tile and so carries no information at all.
    /// The per-sequencing-unit file's `tiles` column distinguishes the last case.
    pub raw_sequencing_duplicate_pairs: Option<u64>,
    /// `raw_sequencing_duplicate_pairs` corrected for tiles that collide by chance,
    /// by inferring how many independent molecules the observed tile count
    /// implies. This is the figure to use.
    ///
    /// Deliberately *not* named `READ_PAIR_OPTICAL_DUPLICATES`: this is a
    /// threshold-free, tile-identity definition rather than Picard's pixel-radius
    /// one, and the two should not invite direct comparison.
    pub corrected_sequencing_duplicate_pairs: Option<u64>,
    /// Duplicate pairs from independent molecules: PCR copies, plus genuinely
    /// distinct molecules that happen to share a locus. The residual of
    /// `duplicate_pairs − corrected_sequencing_duplicate_pairs`, so the two sum exactly.
    pub library_duplicate_pairs: Option<u64>,
    /// `duplicate_pairs / mapped_pairs`, in [0, 1] — the pair-level duplicate
    /// rate, as distinct from the read-level `frac_duplicates` above.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_duplicate_pairs: f64,
    /// `corrected_sequencing_duplicate_pairs / mapped_pairs`, in [0, 1]: what
    /// share of the library's pairs were lost to the flowcell.
    ///
    /// Over `mapped_pairs` rather than over `duplicate_pairs`, so it shares a
    /// denominator with `frac_duplicate_pairs` and is always ≤ it — the two are
    /// directly comparable, and one is a component of the other. Two `frac_`
    /// columns in one row measuring against different denominators would be a trap.
    /// The "what share of my duplication was optical" reading is still available as
    /// the ratio of the two, or from the counts.
    #[serde(serialize_with = "serialize_opt_f64_6dp")]
    pub frac_sequencing_duplicate_pairs: Option<f64>,
    /// Lander-Waterman estimate of the library's distinct molecules, with
    /// sequencing duplicates removed from the observed total — Picard's
    /// convention for `ESTIMATED_LIBRARY_SIZE`, which subtracts them from `n` but
    /// not from the unique count. Flowcell duplicates are not evidence a library
    /// is exhausted, so counting them as saturation understates it badly.
    ///
    /// Empty when it cannot be estimated: no pairs, no duplicate pairs, or — a
    /// degenerate case — when *every* duplicate is a sequencing duplicate, which
    /// leaves the observed total equal to the unique count and nothing for the
    /// estimator to work from.
    pub estimated_library_size: Option<u64>,
    /// Templates where exactly one read is mapped (mate unmapped or absent).
    pub mapped_orphans: u64,
    /// Duplicates among `mapped_orphans`.
    pub duplicate_orphans: u64,
    /// Single-record templates whose only record is unmapped: an unmapped
    /// single-end read (always), or a paired primary whose mate is absent and
    /// itself unmapped (only under `--ignore-unmated`).
    pub unmapped_orphans: u64,
    /// Templates whose paired flag was set but the mate was missing (only
    /// under `--ignore-unmated`).
    pub unmated_templates: u64,
}

/// Serialize an `f64` with 6 decimal places (fixed precision for the duplicate
/// fraction). Used via `#[serde(serialize_with = "serialize_f64_6dp")]`. Shared
/// with [`crate::complexity`] so both TSV outputs render fractions identically.
pub(crate) fn serialize_f64_6dp<S: serde::Serializer>(
    value: &f64,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    serializer.serialize_str(&format!("{value:.6}"))
}

/// Serialize an optional `f64` with 6 decimal places, `None` as an empty cell.
/// Plain `Option<u64>` already renders `None` as empty; a float needs this because
/// the fixed-precision formatting has to be skipped rather than applied to a zero.
fn serialize_opt_f64_6dp<S: serde::Serializer>(
    value: &Option<f64>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    match value {
        Some(value) => serializer.serialize_str(&format!("{value:.6}")),
        None => serializer.serialize_str(""),
    }
}

/// A count as a fraction of `mapped_pairs`, `0.0` when there are none.
///
/// Shared by both pair-level `frac_` columns so their denominator cannot drift
/// apart: they are meant to be directly comparable, with one a component of the
/// other.
fn fraction_of_mapped_pairs(count: u64, mapped_pairs: u64) -> f64 {
    if mapped_pairs == 0 { 0.0 } else { count as f64 / mapped_pairs as f64 }
}

impl Metrics {
    /// Build one TSV row from a single library's [`LibraryStats`] plus the
    /// already-resolved `sample` name. The library name comes from
    /// `library_stats.name`.
    pub fn from_library_stats(
        library_stats: &LibraryStats,
        sample: &str,
        decomposition: Option<&Decomposition>,
    ) -> Self {
        let mapped_orphans = library_stats.mapped_orphan_id_count;
        let mapped_pairs = library_stats.both_mapped_id_count;
        let duplicate_orphans = library_stats.orphan_dup_count;
        let duplicate_pairs = library_stats.both_mapped_dup_count;
        let denom = mapped_orphans + 2 * mapped_pairs;
        let frac_duplicates = if denom == 0 {
            0.0
        } else {
            (duplicate_orphans + 2 * duplicate_pairs) as f64 / denom as f64
        };
        // `duplicate_pairs` should never exceed `mapped_pairs` by construction
        // in `mark_dups`, but use saturating_sub defensively — if a future
        // refactor introduces a bug, we'd rather emit a nonsensical library
        // size than panic on u64 underflow.
        // A library on fewer than two tiles carries no information: every
        // duplicate is on "the" tile whether it was clustered or not. Report the
        // tile evidence so the reason is visible, but not counts that cannot mean
        // anything. The same blanking covers a library with no pairs at all.
        let estimable = decomposition.filter(|split| split.tile_count > 1 && mapped_pairs > 0);
        // Picard's convention: sequencing duplicates come off the observed total,
        // not off the unique count.
        let estimated_library_size = match estimable {
            Some(split) => estimate_library_size(
                mapped_pairs.saturating_sub(split.corrected_sequencing_duplicates),
                mapped_pairs.saturating_sub(duplicate_pairs),
            ),
            None => {
                estimate_library_size(mapped_pairs, mapped_pairs.saturating_sub(duplicate_pairs))
            }
        };
        Self {
            sample: sample.to_string(),
            library: library_stats.name.clone(),
            dupblaster_version: DUPBLASTER_BUILD,
            total_templates: library_stats.id_count,
            duplicate_templates: library_stats.dup_count,
            frac_duplicates,
            mapped_pairs,
            unmapped_pairs: library_stats.both_unmapped_id_count,
            duplicate_pairs,
            raw_sequencing_duplicate_pairs: estimable.map(|s| s.raw_sequencing_duplicates),
            corrected_sequencing_duplicate_pairs: estimable
                .map(|s| s.corrected_sequencing_duplicates),
            library_duplicate_pairs: estimable.map(|s| s.library_duplicates),
            frac_duplicate_pairs: fraction_of_mapped_pairs(duplicate_pairs, mapped_pairs),
            // Same denominator as `frac_duplicate_pairs`, so the two compose.
            frac_sequencing_duplicate_pairs: estimable.map(|split| {
                fraction_of_mapped_pairs(split.corrected_sequencing_duplicates, mapped_pairs)
            }),
            estimated_library_size,
            mapped_orphans,
            duplicate_orphans,
            unmapped_orphans: library_stats.unmapped_orphan_id_count,
            unmated_templates: library_stats.unmated_count,
        }
    }

    /// Build the run-summary rows from end-of-run [`Stats`]: one row per library
    /// that processed at least one template (skipping an empty "Unknown
    /// Library" catch-all). Falls back to a single run-wide "All Reads" row
    /// when no library saw any data. `header` and `sample_override` resolve the
    /// shared `sample` column.
    pub fn rows_from_stats(
        stats: &Stats,
        header: &Header,
        sample_override: Option<&str>,
        decomposition: Option<&[Decomposition]>,
    ) -> Vec<Metrics> {
        let sample = resolve_sample(header, sample_override);
        // Enumerate before filtering: `decomposition` is indexed by library
        // bucket, but rows skip libraries that saw no templates, so the row
        // position is not the bucket index.
        let mut rows: Vec<Metrics> = stats
            .libraries
            .iter()
            .enumerate()
            .filter(|(_, ls)| ls.id_count > 0)
            .map(|(lib, ls)| {
                Metrics::from_library_stats(ls, &sample, decomposition.and_then(|d| d.get(lib)))
            })
            .collect();
        if rows.is_empty() {
            // No library saw data, so there is nothing to decompose either.
            rows.push(Metrics::from_library_stats(&stats.totals(), &sample, None));
        }
        rows
    }
}

/// One row of the per-sequencing-unit QC table written alongside the run summary when
/// `--read-name-format` is on.
///
/// Its own file because the granularity is different: a sequencing unit is a
/// flowcell-and-lane, and one library spans many of them. The variation this
/// exposes is large enough to matter — three flowcells of a single library
/// measured 26.8%, 13.1% and 2.6% sequencing duplicates as a share of their own
/// templates, a 9× spread invisible in the per-library row.
#[derive(Debug, Clone, Serialize)]
pub struct SequencingUnitMetrics {
    /// Sample name, matching the run-summary rows.
    pub sample: String,
    /// Library this unit's reads belong to.
    pub library: String,
    /// The unit as it appeared in the read names, e.g. `H72CFDSXF:2`.
    pub sequencing_unit: String,
    /// Templates observed on this unit.
    pub templates: u64,
    /// Distinct imaging tiles seen on this unit.
    pub tiles: usize,
    /// Sequencing duplicates on this unit's own tiles: each tile contributes
    /// `members - 1` of whatever duplicate group it holds part of. A group
    /// straddling two units splits across them exactly, so this column sums over
    /// all units to `raw_sequencing_duplicate_pairs` in the per-library file.
    pub sequencing_duplicate_pairs: u64,
    /// `sequencing_duplicates / templates`, in [0, 1] — the loading-density
    /// signal, comparable across units.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_sequencing_duplicate_pairs: f64,
}

impl SequencingUnitMetrics {
    /// Build the per-unit rows, resolving library buckets to their names.
    ///
    /// A library that produced no sequencing units — no both-ends-mapped pairs —
    /// still gets a row, with zeros and an empty `sequencing_unit`. Emitting
    /// nothing would leave a header-less empty file that every header-driven
    /// reader chokes on, and a visible zero row says "measured, found nothing"
    /// where an absent file says nothing at all.
    pub fn rows(units: &[SequencingUnitStats], stats: &Stats, sample: &str) -> Vec<Self> {
        let name_of = |library: u32| {
            stats.libraries.get(library as usize).map_or_else(String::new, |ls| ls.name.clone())
        };
        let mut rows: Vec<Self> = units
            .iter()
            .map(|unit| Self {
                sample: sample.to_string(),
                library: name_of(unit.library),
                sequencing_unit: unit.unit.clone(),
                templates: unit.templates,
                tiles: unit.tiles,
                sequencing_duplicate_pairs: unit.sequencing_duplicates,
                frac_sequencing_duplicate_pairs: if unit.templates == 0 {
                    0.0
                } else {
                    unit.sequencing_duplicates as f64 / unit.templates as f64
                },
            })
            .collect();
        for (library, library_stats) in stats.libraries.iter().enumerate() {
            let library = library as u32;
            if library_stats.id_count == 0 || units.iter().any(|unit| unit.library == library) {
                continue;
            }
            rows.push(Self {
                sample: sample.to_string(),
                library: name_of(library),
                sequencing_unit: String::new(),
                templates: 0,
                tiles: 0,
                sequencing_duplicate_pairs: 0,
                frac_sequencing_duplicate_pairs: 0.0,
            });
        }
        rows
    }
}

/// Path of the run-summary table: `<PREFIX>.duplicate-metrics.tsv`.
pub fn duplicate_metrics_path(prefix: &Path) -> PathBuf {
    suffixed(prefix, ".duplicate-metrics.tsv")
}

/// Path of the per-sequencing-unit table: `<PREFIX>.sequencing-units.tsv`. A
/// different granularity from [`duplicate_metrics_path`] — one flowcell-and-lane
/// per row rather than one library — so it gets its own file.
pub fn sequencing_units_path(prefix: &Path) -> PathBuf {
    suffixed(prefix, ".sequencing-units.tsv")
}

/// Append `suffix` to `prefix` as raw bytes rather than via `with_extension`, so
/// a prefix that already contains dots keeps them: `sampleA.hg38` yields
/// `sampleA.hg38.duplicate-metrics.tsv`, not `sampleA.duplicate-metrics.tsv`.
fn suffixed(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// Write the per-sequencing-unit rows to `path` as a tab-separated file.
pub fn write_unit_rows_to_path(rows: &[SequencingUnitMetrics], path: &Path) -> Result<()> {
    DelimFile::default()
        .write_tsv(path, rows.iter())
        .with_context(|| format!("writing sequencing-unit TSV to {}", path.display()))
}

/// Write the metrics rows to `path` as a tab-separated file: a header row of
/// serde field names followed by one row per [`Metrics`]. Serialization is
/// handled by [`fgoxide::io::DelimFile`], so the [`Metrics`] field set is the
/// single source of truth for the columns.
pub fn write_rows_to_path(rows: &[Metrics], path: &Path) -> Result<()> {
    DelimFile::default()
        .write_tsv(path, rows.iter())
        .with_context(|| format!("writing duplicate-metrics TSV to {}", path.display()))
}

/// Resolve the `sample` value: explicit override wins, else comma-join the
/// unique `@RG SM:` tags from the header, else empty string.
pub fn resolve_sample(header: &Header, sample_override: Option<&str>) -> String {
    if let Some(s) = sample_override {
        return s.to_string();
    }
    let mut samples: BTreeSet<String> = BTreeSet::new();
    for (_id, map) in header.read_groups() {
        if let Some(sm) = map.other_fields().get(&rg_tag::SAMPLE) {
            let s = sm.to_string();
            if !s.is_empty() {
                samples.insert(s);
            }
        }
    }
    samples.into_iter().collect::<Vec<_>>().join(",")
}

/// Estimate the library size from observed read pairs using the
/// Lander-Waterman formula. Ports Picard's `estimateLibrarySize`:
/// finds `x` such that `c/x - 1 + exp(-n/x) = 0`, where `n` is total pairs
/// observed and `c` is unique pairs observed. Returns `None` if there are
/// no pairs or no duplicates (in which case the library is effectively
/// "infinite" given the observation).
pub fn estimate_library_size(read_pairs: u64, unique_read_pairs: u64) -> Option<u64> {
    if read_pairs == 0 || unique_read_pairs >= read_pairs {
        return None;
    }
    let n = read_pairs as f64;
    let c = unique_read_pairs as f64;
    // Multipliers of `c` for the bisection bounds. Picard uses [1, 100]
    // initially and expands the upper bound by 10× until f(M*c) > 0.
    let mut lo = 1.0_f64;
    let mut hi = 100.0_f64;
    // Faithful port of Picard's lower-bound sanity check. It is in fact
    // unreachable for valid inputs: with `lo == 1`, `f_lw(c, c, n) =
    // exp(-n/c)`, which is strictly positive for any finite `n, c > 0` (and
    // `read_pairs > 0` is guaranteed above). Kept for bit-for-bit parity
    // with Picard rather than removed.
    if f_lw(lo * c, c, n) < 0.0 {
        return None;
    }
    while f_lw(hi * c, c, n) > 0.0 {
        hi *= 10.0;
        if !hi.is_finite() {
            return None;
        }
    }
    for _ in 0..40 {
        let mid = (lo + hi) / 2.0;
        let v = f_lw(mid * c, c, n);
        if v == 0.0 {
            break;
        } else if v > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let est = c * (lo + hi) / 2.0;
    if est.is_finite() && est >= 0.0 { Some(est as u64) } else { None }
}

/// The Lander-Waterman residual used by the bisection: `c/x - 1 + exp(-n/x)`.
fn f_lw(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write `rows` to a temp TSV via the real [`write_rows_to_path`] writer and
    /// return the file contents, so tests assert on the actual serialized output.
    fn write_rows_to_string(rows: &[Metrics]) -> String {
        let tmp = tempfile::NamedTempFile::new().expect("temp file");
        write_rows_to_path(rows, tmp.path()).expect("write rows");
        std::fs::read_to_string(tmp.path()).expect("read back")
    }

    #[test]
    fn library_size_returns_none_when_no_pairs() {
        assert_eq!(estimate_library_size(0, 0), None);
    }

    #[test]
    fn library_size_returns_none_when_no_dups() {
        // All pairs unique → library is effectively infinite.
        assert_eq!(estimate_library_size(1000, 1000), None);
    }

    #[test]
    fn library_size_is_sensible_at_50pct_dup() {
        // At 50% dup rate over 1M pairs, the Picard formula yields a
        // library size in the small-millions ballpark. We don't need
        // an exact value here — just that it's finite and meaningfully
        // larger than the observed unique count.
        let est = estimate_library_size(1_000_000, 500_000).expect("estimable");
        assert!(est > 500_000, "library size {est} should exceed unique pairs");
        assert!(est < 100_000_000, "library size {est} should be finite");
    }

    #[test]
    fn library_size_handles_extreme_low_dup_rate_without_panicking() {
        // 10M pairs with a single duplicate is an extreme low-dup regime
        // that stresses the bisection's upper-bound expansion. The result
        // must be either a sensible finite estimate or None (Picard's
        // give-up), never a panic or a value below the unique count.
        if let Some(est) = estimate_library_size(10_000_000, 9_999_999) {
            assert!(est >= 9_999_999, "library size {est} should be >= the unique count");
        }
    }

    #[test]
    fn library_size_grows_as_dup_rate_drops() {
        let est_high_dup = estimate_library_size(1_000_000, 200_000).unwrap();
        let est_low_dup = estimate_library_size(1_000_000, 900_000).unwrap();
        assert!(
            est_low_dup > est_high_dup,
            "lower dup rate ({est_low_dup}) should imply larger library than higher dup rate ({est_high_dup})"
        );
    }

    #[test]
    fn frac_duplicates_uses_picard_read_level_formula() {
        // 100 mapped pairs (200 mapped reads), 30 pair dups (60 dup reads).
        // 50 mapped orphans (50 reads), 10 orphan dups (10 reads).
        // Expected: (60 + 10) / (200 + 50) = 70 / 250 = 0.28.
        let ls = LibraryStats {
            name: "lib1".to_string(),
            id_count: 200,
            dup_count: 40,
            both_mapped_id_count: 100,
            both_mapped_dup_count: 30,
            mapped_orphan_id_count: 50,
            orphan_dup_count: 10,
            ..Default::default()
        };
        let m = Metrics::from_library_stats(&ls, "", None);
        assert!((m.frac_duplicates - 0.28).abs() < 1e-9, "got {}", m.frac_duplicates);
    }

    #[test]
    fn frac_duplicates_is_zero_when_no_mapped_data() {
        let ls = LibraryStats { id_count: 5, both_unmapped_id_count: 5, ..Default::default() };
        let m = Metrics::from_library_stats(&ls, "", None);
        assert_eq!(m.frac_duplicates, 0.0);
    }

    #[test]
    fn sample_override_wins_over_header() {
        let header = Header::default();
        let s = resolve_sample(&header, Some("forced"));
        assert_eq!(s, "forced");
    }

    #[test]
    fn sample_empty_when_no_override_and_no_read_groups() {
        let header = Header::default();
        let s = resolve_sample(&header, None);
        assert_eq!(s, "");
    }

    #[test]
    fn tsv_header_and_value_have_same_column_count() {
        let ls = LibraryStats {
            name: "lib1".to_string(),
            id_count: 10,
            both_mapped_id_count: 5,
            both_mapped_dup_count: 1,
            ..Default::default()
        };
        let m = Metrics::from_library_stats(&ls, "test", None);
        let text = write_rows_to_string(&[m]);
        let mut lines = text.lines();
        let hdr_cols = lines.next().unwrap().split('\t').count();
        let val_cols = lines.next().unwrap().split('\t').count();
        // Equality matters most for the trailing columns, which are `None` here:
        // a writer that dropped empty cells would silently shorten the row.
        assert_eq!(hdr_cols, val_cols);
        assert_eq!(hdr_cols, 19, "expected 19 metric columns");
    }

    #[test]
    fn rows_from_stats_emits_one_row_per_nonempty_library() {
        // Bucket 0 (Unknown Library) saw no data and must be skipped; the two
        // real libraries each get a row, with their own counts.
        let stats = Stats {
            libraries: vec![
                LibraryStats { name: "Unknown Library".to_string(), ..Default::default() },
                LibraryStats {
                    name: "libA".to_string(),
                    id_count: 3,
                    both_mapped_id_count: 3,
                    both_mapped_dup_count: 1,
                    ..Default::default()
                },
                LibraryStats {
                    name: "libB".to_string(),
                    id_count: 2,
                    both_mapped_id_count: 2,
                    ..Default::default()
                },
            ],
            clamped_template_count: 0,
        };
        let rows = Metrics::rows_from_stats(&stats, &Header::default(), None, None);
        let names: Vec<&str> = rows.iter().map(|m| m.library.as_str()).collect();
        assert_eq!(names, ["libA", "libB"]);
        assert_eq!(rows[0].mapped_pairs, 3);
        assert_eq!(rows[0].duplicate_pairs, 1);
        assert_eq!(rows[1].mapped_pairs, 2);
        assert_eq!(rows[1].duplicate_pairs, 0);
    }

    #[test]
    fn unestimable_library_size_renders_as_empty_cell() {
        let ls = LibraryStats {
            id_count: 10,
            both_mapped_id_count: 10,
            both_mapped_dup_count: 0, // no dups → library size None
            ..Default::default()
        };
        let m = Metrics::from_library_stats(&ls, "", None);
        assert!(m.estimated_library_size.is_none());
        let text = write_rows_to_string(&[m]);
        // Look the column up by name rather than by position: an earlier version
        // asserted on the *last* field, which only worked while the library-size
        // column happened to be last and passed by accident once it moved.
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let values: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let index = header
            .iter()
            .position(|column| *column == "estimated_library_size")
            .expect("column present");
        assert_eq!(values[index], "", "library-size cell should be empty");
    }
}
