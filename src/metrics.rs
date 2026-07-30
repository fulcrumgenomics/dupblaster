//! Run-summary metrics emitted as a per-library TSV (one row per library) when
//! `--stats PATH` is set on the command line.
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
    /// Duplicates among `mapped_pairs`.
    pub duplicate_pairs: u64,
    /// Templates where exactly one read is mapped (mate unmapped or absent).
    pub mapped_orphans: u64,
    /// Duplicates among `mapped_orphans`.
    pub duplicate_orphans: u64,
    /// Single-record templates whose only record is unmapped: an unmapped
    /// single-end read (always), or a paired primary whose mate is absent and
    /// itself unmapped (only under `--ignore-unmated`).
    pub unmapped_orphans: u64,
    /// Templates where both reads are unmapped — never dup-checked.
    pub unmapped_pairs: u64,
    /// Templates whose paired flag was set but the mate was missing (only
    /// under `--ignore-unmated`).
    pub unmated_templates: u64,
    /// `None` when no duplicate pairs were observed or there were no pairs
    /// to estimate from — written as an empty cell.
    pub estimated_library_size: Option<u64>,
    /// Duplicate pairs that are copies of one molecule made on the flowcell
    /// (cluster/ExAmp duplicates), counted as `E[tiles | k] − tiles` per group
    /// and so already corrected for tiles that collide by chance.
    ///
    /// Deliberately *not* named `READ_PAIR_OPTICAL_DUPLICATES`: this is a
    /// threshold-free, tile-identity definition rather than Picard's pixel-radius
    /// one, and the two should not invite direct comparison.
    ///
    /// `None` — an empty cell — whenever the split was not computed:
    /// `--read-name-format` was not given, the library has no pairs, or it sits
    /// on a single tile and therefore carries no information. `tile_count` and
    /// `tile_collision_rate` say which.
    pub sequencing_duplicates: Option<u64>,
    /// Duplicate pairs from independent molecules: PCR copies, plus genuinely
    /// distinct molecules that happen to share a locus. The residual of
    /// `duplicate_pairs − sequencing_duplicates`, so the two always sum exactly.
    pub library_duplicates: Option<u64>,
    /// `sequencing_duplicates / duplicate_pairs`, in [0, 1].
    #[serde(serialize_with = "serialize_opt_f64_6dp")]
    pub frac_sequencing_duplicates: Option<f64>,
    /// Uncorrected `Σ (k − tiles)`. The gap to `sequencing_duplicates` is how
    /// much of the naive count was chance tile collision: negligible on WGS,
    /// over 20% for groups above 100 members.
    pub naive_sequencing_duplicates: Option<u64>,
    /// Distinct imaging tiles seen for this library. A misconfigured
    /// `--read-name-format` shows up here as an implausibly large number.
    pub tile_count: Option<usize>,
    /// `q = Σ w_t²`: the chance two unrelated templates of this library share a
    /// tile. The validity indicator for the split — at `q` near 1 there is only
    /// one tile in play and no duplicate can be attributed, which is why it is
    /// reported even when the counts are blank.
    #[serde(serialize_with = "serialize_opt_f64_6dp")]
    pub tile_collision_rate: Option<f64>,
    /// Library size re-estimated with sequencing duplicates removed from the
    /// observed total, following Picard's convention of subtracting optical
    /// duplicates from `n` but not from the unique count.
    ///
    /// Flowcell duplicates say nothing about how many distinct molecules the
    /// library held, so counting them as evidence of saturation makes a library
    /// look far smaller than it is — worth 3.1x on one 30x WGS sample
    /// (73.0M → 226.4M). Reported next to the uncorrected
    /// `estimated_library_size` rather than replacing it, so the plain column
    /// stays comparable across runs whether or not the split was computed.
    pub estimated_library_size_corrected: Option<u64>,
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
        let estimated_library_size =
            estimate_library_size(mapped_pairs, mapped_pairs.saturating_sub(duplicate_pairs));
        // A library on fewer than two tiles carries no information: every
        // duplicate is on "the" tile whether it was clustered or not. Report the
        // tile evidence so the reason is visible, but not counts that cannot mean
        // anything. The same blanking covers a library with no pairs at all.
        let estimable = decomposition.filter(|split| split.tile_count > 1 && mapped_pairs > 0);
        Self {
            sample: sample.to_string(),
            library: library_stats.name.clone(),
            dupblaster_version: DUPBLASTER_BUILD,
            total_templates: library_stats.id_count,
            duplicate_templates: library_stats.dup_count,
            frac_duplicates,
            mapped_pairs,
            duplicate_pairs,
            mapped_orphans,
            duplicate_orphans,
            unmapped_orphans: library_stats.unmapped_orphan_id_count,
            unmapped_pairs: library_stats.both_unmapped_id_count,
            unmated_templates: library_stats.unmated_count,
            estimated_library_size,
            sequencing_duplicates: estimable.map(|split| split.sequencing_duplicates),
            library_duplicates: estimable.map(|split| split.library_duplicates),
            frac_sequencing_duplicates: estimable.and_then(|split| {
                (split.duplicate_pairs > 0)
                    .then(|| split.sequencing_duplicates as f64 / split.duplicate_pairs as f64)
            }),
            naive_sequencing_duplicates: estimable.map(|split| split.naive_sequencing_duplicates),
            tile_count: decomposition.map(|split| split.tile_count),
            tile_collision_rate: decomposition.map(|split| split.tile_collision_rate),
            estimated_library_size_corrected: estimable.and_then(|split| {
                estimate_library_size(
                    mapped_pairs.saturating_sub(split.sequencing_duplicates),
                    mapped_pairs.saturating_sub(duplicate_pairs),
                )
            }),
        }
    }

    /// Build the `--stats` rows from end-of-run [`Stats`]: one row per library
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

/// One row of the per-sequencing-unit QC table written alongside `--stats` when
/// `--read-name-format` is on.
///
/// Its own file because the granularity is different: a sequencing unit is a
/// flowcell-and-lane, and one library spans many of them. The variation this
/// exposes is large enough to matter — three flowcells of a single library
/// measured 26.8%, 13.1% and 2.6% sequencing duplicates as a share of their own
/// templates, a 9× spread invisible in the per-library row.
#[derive(Debug, Clone, Serialize)]
pub struct SequencingUnitMetrics {
    /// Sample name, matching the `--stats` rows.
    pub sample: String,
    /// Library this unit's reads belong to.
    pub library: String,
    /// The unit as it appeared in the read names, e.g. `H72CFDSXF:2`.
    pub sequencing_unit: String,
    /// Templates observed on this unit.
    pub templates: u64,
    /// Distinct imaging tiles seen on this unit.
    pub tiles: usize,
    /// Sequencing duplicates credited to this unit. A duplicate group straddling
    /// two units is credited whole to whichever holds most of its members — a
    /// cluster duplicate physically happened on one flowcell, so splitting it
    /// would be meaningless.
    pub sequencing_duplicates: u64,
    /// `sequencing_duplicates / templates`, in [0, 1] — the loading-density
    /// signal, comparable across units.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_sequencing_duplicates: f64,
}

impl SequencingUnitMetrics {
    /// Build the per-unit rows, resolving library buckets to their names.
    pub fn rows(units: &[SequencingUnitStats], stats: &Stats, sample: &str) -> Vec<Self> {
        units
            .iter()
            .map(|unit| Self {
                sample: sample.to_string(),
                library: stats
                    .libraries
                    .get(unit.library as usize)
                    .map_or_else(String::new, |ls| ls.name.clone()),
                sequencing_unit: unit.unit.clone(),
                templates: unit.templates,
                tiles: unit.tiles,
                sequencing_duplicates: unit.sequencing_duplicates,
                frac_sequencing_duplicates: if unit.templates == 0 {
                    0.0
                } else {
                    unit.sequencing_duplicates as f64 / unit.templates as f64
                },
            })
            .collect()
    }
}

/// Path of the per-sequencing-unit table, derived from the `--stats` path by
/// inserting `.sequencing-units` before the extension. Derived rather than given
/// its own flag so enabling the decomposition needs one option, not two.
pub fn sequencing_units_path(stats_path: &Path) -> PathBuf {
    let extension = stats_path.extension().and_then(|ext| ext.to_str()).unwrap_or("tsv");
    stats_path.with_extension(format!("sequencing-units.{extension}"))
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
/// single source of truth for the columns. A `.gz`/`.bgz` suffix on `path`
/// transparently gzip-compresses the output.
pub fn write_rows_to_path(rows: &[Metrics], path: &Path) -> Result<()> {
    DelimFile::default()
        .write_tsv(path, rows.iter())
        .with_context(|| format!("writing stats TSV to {}", path.display()))
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
        assert_eq!(hdr_cols, 21, "expected 21 metric columns");
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
        // The last column on the value row should be empty (rendered as
        // an empty cell after the final tab).
        // Assert the *last field* is empty rather than relying on a trailing
        // tab — the latter passes by accident regardless of which column is
        // last.
        let value_line = text.lines().nth(1).unwrap();
        let last_field = value_line.rsplit('\t').next().unwrap();
        assert_eq!(last_field, "", "library-size cell should be empty, line: {value_line:?}");
    }
}
