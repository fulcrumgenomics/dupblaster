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
use std::path::Path;

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;
use serde::Serialize;

use crate::DUPBLASTER_BUILD;
use crate::dedup::{LibraryStats, Stats};

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
    /// Single-record templates whose only record is unmapped (only under
    /// `--ignore-unmated`).
    pub unmapped_orphans: u64,
    /// Templates where both reads are unmapped — never dup-checked.
    pub unmapped_pairs: u64,
    /// Templates whose paired flag was set but the mate was missing (only
    /// under `--ignore-unmated`).
    pub unmated_templates: u64,
    /// `None` when no duplicate pairs were observed or there were no pairs
    /// to estimate from — written as an empty cell.
    pub estimated_library_size: Option<u64>,
}

/// Serialize an `f64` with 6 decimal places (fixed precision for the duplicate
/// fraction). Used via `#[serde(serialize_with = "serialize_f64_6dp")]`.
fn serialize_f64_6dp<S: serde::Serializer>(value: &f64, serializer: S) -> Result<S::Ok, S::Error> {
    serializer.serialize_str(&format!("{value:.6}"))
}

impl Metrics {
    /// Build one TSV row from a single library's [`LibraryStats`] plus the
    /// already-resolved `sample` name. The library name comes from
    /// `library_stats.name`.
    pub fn from_library_stats(library_stats: &LibraryStats, sample: &str) -> Self {
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
    ) -> Vec<Metrics> {
        let sample = resolve_sample(header, sample_override);
        let mut rows: Vec<Metrics> = stats
            .libraries
            .iter()
            .filter(|ls| ls.id_count > 0)
            .map(|ls| Metrics::from_library_stats(ls, &sample))
            .collect();
        if rows.is_empty() {
            rows.push(Metrics::from_library_stats(&stats.totals(), &sample));
        }
        rows
    }
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
        let m = Metrics::from_library_stats(&ls, "");
        assert!((m.frac_duplicates - 0.28).abs() < 1e-9, "got {}", m.frac_duplicates);
    }

    #[test]
    fn frac_duplicates_is_zero_when_no_mapped_data() {
        let ls = LibraryStats { id_count: 5, both_unmapped_id_count: 5, ..Default::default() };
        let m = Metrics::from_library_stats(&ls, "");
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
        let m = Metrics::from_library_stats(&ls, "test");
        let text = write_rows_to_string(&[m]);
        let mut lines = text.lines();
        let hdr_cols = lines.next().unwrap().split('\t').count();
        let val_cols = lines.next().unwrap().split('\t').count();
        assert_eq!(hdr_cols, val_cols);
        assert_eq!(hdr_cols, 14, "expected 14 metric columns");
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
        let rows = Metrics::rows_from_stats(&stats, &Header::default(), None);
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
        let m = Metrics::from_library_stats(&ls, "");
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
