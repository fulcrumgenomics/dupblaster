//! Duplication-complexity metrics: the sampled ladder here, the group-size
//! spectrum in [`crate::counts`].
//!
//! This module produces the **duplication-sampled ladder**: periodic snapshots
//! taken every N templates as the input streams by, each carrying a cumulative
//! set (`total`, `unique`, `duplicates`, `frac_duplicates`) and a matching
//! per-window set (`window_total`, `window_unique`, `window_duplicates`,
//! `window_frac_duplicates`) covering the templates added since the previous
//! snapshot. Plotted as duplication rate vs. depth; the window columns are
//! usually the more legible view.
//!
//! **The ladder's shape depends on the order templates arrive in.** On
//! coordinate- or queryname-sorted multi-flowcell input the curve is dominated
//! by contiguous per-read-group blocks (flowcell/ExAmp duplicates cluster within
//! a read group and land adjacent under sorting), so it is an order-dependent
//! *diagnostic*, not a complexity estimator: use the order-independent count
//! histogram in [`crate::counts`] (η_k) for library complexity. The ladder is
//! cleanest on single-lane or homogeneous-lane input; that same order dependence
//! is what a later flowcell-vs-library decomposition can exploit.
//!
//! ## One category per library
//!
//! Each library is reported on exactly one category (matching the histogram in
//! [`crate::counts`]): `pairs` (both-ends-mapped templates) if it has any
//! pairs, otherwise `single_end` (mapped orphan/single-end templates). Both
//! categories are sampled on **their own** running count — `pairs` every N
//! both-mapped templates, `single_end` every N mapped orphans — so each curve
//! is evenly spaced in its own depth even under picard-exact, where pairs
//! stream in pass 1 and fragments in pass 2. At finalize the non-reported
//! category's rows are dropped.
//!
//! Everything here is gated behind `--duplication-sampled`; when the flag is off
//! no recorder is constructed and the hot path is byte-for-byte unchanged.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use serde::Serialize;

use crate::dedup::{CATEGORY_PAIRS, CATEGORY_SINGLE_END, LibraryStats, Stats};
use crate::metrics::serialize_f64_6dp;

/// Snapshot bookkeeping for one `(library, category)` curve: the next running
/// count at which to emit a snapshot, and the counters at the previous snapshot
/// (needed for the marginal per-window delta). The cumulative totals
/// themselves are not stored here — they live in [`LibraryStats`], which
/// [`LadderRecorder::observe`] reads.
#[derive(Clone)]
struct CategoryLadder {
    /// Next running count (in this category) at which to take a snapshot.
    next: u64,
    /// Running count at the previous snapshot (0 before the first).
    last_total: u64,
    /// `duplicates` at the previous snapshot (0 before the first).
    last_dup: u64,
}

impl CategoryLadder {
    fn new(interval: u64) -> Self {
        Self { next: interval, last_total: 0, last_dup: 0 }
    }
}

/// Per-library ladder state: independent snapshot bookkeeping for the library's
/// `pairs` and `single_end` curves.
#[derive(Clone)]
struct LibraryLadder {
    pairs: CategoryLadder,
    single_end: CategoryLadder,
}

impl LibraryLadder {
    fn new(interval: u64) -> Self {
        Self { pairs: CategoryLadder::new(interval), single_end: CategoryLadder::new(interval) }
    }
}

/// Accumulates ladder snapshots as templates stream by. One recorder per run,
/// with per-library, per-category snapshot bookkeeping; constructed only when
/// `--duplication-sampled` is on.
pub struct LadderRecorder {
    /// Snapshot cadence in templates (the `--sampling-interval` value).
    interval: u64,
    /// Per-library snapshot bookkeeping, indexed by library bucket.
    libs: Vec<LibraryLadder>,
    /// Resolved sample name, copied into every row.
    sample: String,
    /// Accumulated rows, filtered to the reported category and sorted in
    /// [`Self::finalize`].
    rows: Vec<DuplicateLadderRow>,
}

impl LadderRecorder {
    /// Build a recorder for `num_libs` library buckets snapshotting every
    /// `interval` templates (per category), stamping `sample` into every row.
    pub fn new(num_libs: usize, interval: u64, sample: String) -> Self {
        Self {
            interval,
            libs: vec![LibraryLadder::new(interval); num_libs],
            sample,
            rows: Vec::new(),
        }
    }

    /// Record one processed template for library `lib`, given that library's
    /// current running counters. Emits a `pairs` snapshot each time the
    /// library's both-mapped count crosses a multiple of `interval`, and a
    /// `single_end` snapshot each time its mapped-orphan count does. Cheap
    /// enough for the hot path: two comparisons per template.
    pub fn observe(&mut self, lib: u32, stats: &LibraryStats) {
        let i = lib as usize;
        Self::snapshot(
            &mut self.rows,
            &self.sample,
            self.interval,
            &mut self.libs[i].pairs,
            &stats.name,
            CATEGORY_PAIRS,
            stats.both_mapped_id_count,
            stats.both_mapped_dup_count,
        );
        Self::snapshot(
            &mut self.rows,
            &self.sample,
            self.interval,
            &mut self.libs[i].single_end,
            &stats.name,
            CATEGORY_SINGLE_END,
            stats.mapped_orphan_id_count,
            stats.orphan_dup_count,
        );
    }

    /// Emit snapshots for one category up to the current `total`, advancing its
    /// threshold each time. Each row records the *current* `(total, duplicates)`
    /// against the previous snapshot's counters (the marginal delta). A `while`
    /// (not `if`) so a jump larger than one interval can't skip a snapshot.
    #[expect(
        clippy::too_many_arguments,
        reason = "disjoint field borrows of the recorder plus the category's counters; \
                  bundling them would just reintroduce the borrow conflict this signature avoids"
    )]
    fn snapshot(
        rows: &mut Vec<DuplicateLadderRow>,
        sample: &str,
        interval: u64,
        cat: &mut CategoryLadder,
        library: &str,
        category: &'static str,
        total: u64,
        duplicates: u64,
    ) {
        while total >= cat.next {
            rows.push(build_ladder_row(
                sample,
                library,
                category,
                total,
                duplicates,
                cat.last_total,
                cat.last_dup,
            ));
            cat.last_total = total;
            cat.last_dup = duplicates;
            cat.next += interval;
        }
    }

    /// After the stream ends: add a final snapshot for each library's reported
    /// category (so the ladder ends at the true total), drop the non-reported
    /// category's rows, and sort for a stable, plottable file order.
    pub fn finalize(&mut self, stats: &Stats) {
        for (i, ls) in stats.libraries.iter().enumerate() {
            if ls.id_count == 0 {
                continue;
            }
            if ls.both_mapped_id_count > 0 {
                let cat = &self.libs[i].pairs;
                if cat.last_total != ls.both_mapped_id_count {
                    self.rows.push(build_ladder_row(
                        &self.sample,
                        &ls.name,
                        CATEGORY_PAIRS,
                        ls.both_mapped_id_count,
                        ls.both_mapped_dup_count,
                        cat.last_total,
                        cat.last_dup,
                    ));
                }
            } else {
                let cat = &self.libs[i].single_end;
                if cat.last_total != ls.mapped_orphan_id_count {
                    self.rows.push(build_ladder_row(
                        &self.sample,
                        &ls.name,
                        CATEGORY_SINGLE_END,
                        ls.mapped_orphan_id_count,
                        ls.orphan_dup_count,
                        cat.last_total,
                        cat.last_dup,
                    ));
                }
            }
        }
        // Keep only the reported category per library: `pairs` if it has any
        // pairs, else `single_end`. (Under picard-exact a paired library still
        // accrued `single_end` snapshots in pass 2 for its orphans; drop them.)
        let reported: HashMap<&str, &'static str> =
            stats.libraries.iter().map(|ls| (ls.name.as_str(), ls.reported_category())).collect();
        self.rows.retain(|r| reported.get(r.library.as_str()).copied() == Some(r.category));
        self.rows.sort_by(|a, b| {
            (a.library.as_str(), a.category, a.total).cmp(&(
                b.library.as_str(),
                b.category,
                b.total,
            ))
        });
    }

    /// The accumulated rows (call [`Self::finalize`] first).
    pub fn rows(&self) -> &[DuplicateLadderRow] {
        &self.rows
    }
}

/// Build one ladder row, deriving the cumulative (`unique`, `frac_duplicates`)
/// and per-window (`window_*`) fields. `prev_total`/`prev_duplicates` are the
/// previous snapshot's counters for this `(library, category)` series (both 0 for
/// the first snapshot). `duplicates` is `<= total` and monotonic across snapshots
/// by construction (a template is a duplicate only if the same molecule was
/// already present).
fn build_ladder_row(
    sample: &str,
    library: &str,
    category: &'static str,
    total: u64,
    duplicates: u64,
    prev_total: u64,
    prev_duplicates: u64,
) -> DuplicateLadderRow {
    let frac = if total == 0 { 0.0 } else { duplicates as f64 / total as f64 };
    let window_total = total - prev_total;
    let window_duplicates = duplicates - prev_duplicates;
    let window_frac =
        if window_total == 0 { 0.0 } else { window_duplicates as f64 / window_total as f64 };
    DuplicateLadderRow {
        sample: sample.to_string(),
        library: library.to_string(),
        category,
        total,
        unique: total - duplicates,
        duplicates,
        frac_duplicates: frac,
        window_total,
        window_unique: window_total - window_duplicates,
        window_duplicates,
        window_frac_duplicates: window_frac,
    }
}

/// One row of the duplication-sampled ladder: a `(library, category)` observed at
/// a cumulative depth of `total` templates. Each snapshot carries a matching
/// **cumulative** set (`total`/`unique`/`duplicates`/`frac_duplicates`) and
/// **window** set (`window_*`, covering just the templates added since the
/// previous snapshot), so the file is directly plottable either way without
/// post-processing.
///
/// Field declaration order is the TSV column order; `sample` is first by
/// convention (matches [`crate::metrics::Metrics`]).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DuplicateLadderRow {
    /// Sample name — `--sample` if set, else comma-joined unique `@RG SM:`.
    pub sample: String,
    /// Library name for this row.
    pub library: String,
    /// `pairs` for a library with any both-mapped pairs, else `single_end`.
    pub category: &'static str,
    /// Cumulative templates observed in this category at this snapshot (x-axis).
    pub total: u64,
    /// Distinct molecules observed so far = `total - duplicates`.
    pub unique: u64,
    /// Templates in this category flagged as duplicates so far.
    pub duplicates: u64,
    /// Cumulative `duplicates / total` in `[0, 1]`, 6 dp (0 when `total == 0`).
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_duplicates: f64,
    /// Templates added *in this window* (since the previous snapshot); equals
    /// `total` for the first snapshot.
    pub window_total: u64,
    /// Distinct molecules in this window = `window_total - window_duplicates`.
    pub window_unique: u64,
    /// Duplicates flagged in this window = `duplicates` minus the previous
    /// snapshot's.
    pub window_duplicates: u64,
    /// Marginal duplicate fraction over this window = `window_duplicates /
    /// window_total` in `[0, 1]`, 6 dp. Usually the more legible view than the
    /// cumulative `frac_duplicates`.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub window_frac_duplicates: f64,
}

/// `foo/sampleA` → `foo/sampleA.duplication-sampled.tsv`.
pub fn ladder_path(prefix: &Path) -> std::path::PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(".duplication-sampled.tsv");
    std::path::PathBuf::from(name)
}

/// Write the ladder rows to `path` as a TSV (header + one row each), via
/// [`fgoxide::io::DelimFile`].
pub fn write_ladder_rows(rows: &[DuplicateLadderRow], path: &Path) -> Result<()> {
    DelimFile::default()
        .write_tsv(path, rows.iter())
        .with_context(|| format!("writing duplication-sampled TSV to {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `LibraryStats` with the counters the ladder reads.
    fn ls(
        name: &str,
        both_mapped: u64,
        both_mapped_dup: u64,
        orphan: u64,
        orphan_dup: u64,
    ) -> LibraryStats {
        LibraryStats {
            name: name.to_string(),
            id_count: both_mapped + orphan,
            dup_count: both_mapped_dup + orphan_dup,
            both_mapped_id_count: both_mapped,
            both_mapped_dup_count: both_mapped_dup,
            mapped_orphan_id_count: orphan,
            orphan_dup_count: orphan_dup,
            ..Default::default()
        }
    }

    fn row<'a>(
        rows: &'a [DuplicateLadderRow],
        library: &str,
        category: &str,
        total: u64,
    ) -> &'a DuplicateLadderRow {
        rows.iter()
            .find(|r| r.library == library && r.category == category && r.total == total)
            .unwrap_or_else(|| panic!("no {library}/{category}/total={total} row"))
    }

    #[test]
    fn pairs_curve_snapshots_every_interval_and_derives_fields() {
        let mut rec = LadderRecorder::new(1, 10, "s".to_string());
        for i in 1..=25u64 {
            rec.observe(0, &ls("lib", i, i / 2, 0, 0)); // ~50% pair dup
        }
        let at10 = row(rec.rows_unfiltered(), "lib", "pairs", 10);
        assert_eq!(at10.duplicates, 5);
        assert_eq!(at10.unique, 5);
        assert!((at10.frac_duplicates - 0.5).abs() < 1e-9);
        assert!(rec.rows_unfiltered().iter().any(|r| r.category == "pairs" && r.total == 20));
    }

    #[test]
    fn paired_library_reports_only_pairs_and_ends_at_true_total() {
        let mut rec = LadderRecorder::new(1, 10, "s".to_string());
        // A PE library with a few orphans mixed in.
        for i in 1..=23u64 {
            rec.observe(0, &ls("lib", i, 3, i / 5, 0));
        }
        let stats = Stats { libraries: vec![ls("lib", 23, 3, 4, 0)], clamped_template_count: 0 };
        rec.finalize(&stats);
        // Only pairs rows survive; ends at the true both-mapped total (23).
        assert!(rec.rows().iter().all(|r| r.category == "pairs"));
        let totals: Vec<u64> = rec.rows().iter().map(|r| r.total).collect();
        assert_eq!(totals, vec![10, 20, 23]);
    }

    #[test]
    fn single_end_only_library_reports_single_end() {
        let mut rec = LadderRecorder::new(1, 10, "s".to_string());
        for i in 1..=15u64 {
            rec.observe(0, &ls("lib", 0, 0, i, i / 3));
        }
        let stats = Stats { libraries: vec![ls("lib", 0, 0, 15, 5)], clamped_template_count: 0 };
        rec.finalize(&stats);
        assert!(rec.rows().iter().all(|r| r.category == "single_end"));
        let totals: Vec<u64> = rec.rows().iter().map(|r| r.total).collect();
        assert_eq!(totals, vec![10, 15]);
    }

    #[test]
    fn no_duplicate_final_row_on_exact_boundary() {
        let mut rec = LadderRecorder::new(1, 10, "s".to_string());
        for i in 1..=20u64 {
            rec.observe(0, &ls("lib", i, 0, 0, 0));
        }
        let stats = Stats { libraries: vec![ls("lib", 20, 0, 0, 0)], clamped_template_count: 0 };
        rec.finalize(&stats);
        let n_at_20 = rec.rows().iter().filter(|r| r.total == 20).count();
        assert_eq!(n_at_20, 1);
    }

    #[test]
    fn tsv_has_sample_first_and_expected_columns() {
        let mut rec = LadderRecorder::new(1, 5, "NA12878".to_string());
        rec.observe(0, &ls("lib", 5, 1, 0, 0));
        rec.finalize(&Stats { libraries: vec![ls("lib", 5, 1, 0, 0)], clamped_template_count: 0 });
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_ladder_rows(rec.rows(), tmp.path()).unwrap();
        let text = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(
            text.lines().next().unwrap(),
            "sample\tlibrary\tcategory\ttotal\tunique\tduplicates\tfrac_duplicates\t\
             window_total\twindow_unique\twindow_duplicates\twindow_frac_duplicates"
        );
        let first_val = text.lines().nth(1).unwrap();
        assert!(first_val.starts_with("NA12878\t"), "sample must be column 1: {first_val}");
        assert!(first_val.contains("0.200000"), "frac at 6dp: {first_val}");
    }

    #[test]
    fn window_columns_are_marginal_deltas_between_snapshots() {
        let mut rec = LadderRecorder::new(1, 10, "s".to_string());
        // Cumulative duplicates at the snapshot points: dup[10]=2, dup[20]=7,
        // dup[23]=9. Values between snapshots are never read.
        for i in 1..=23u64 {
            let dup = if i <= 10 {
                2
            } else if i <= 20 {
                7
            } else {
                9
            };
            rec.observe(0, &ls("lib", i, dup, 0, 0));
        }
        rec.finalize(&Stats { libraries: vec![ls("lib", 23, 9, 0, 0)], clamped_template_count: 0 });

        // First snapshot: window == cumulative (previous is zero).
        let at10 = row(rec.rows(), "lib", "pairs", 10);
        assert_eq!((at10.window_total, at10.window_duplicates, at10.window_unique), (10, 2, 8));
        assert!((at10.window_frac_duplicates - 0.2).abs() < 1e-9);
        // Full window: 7 - 2 = 5 dups over 10 templates.
        let at20 = row(rec.rows(), "lib", "pairs", 20);
        assert_eq!((at20.window_total, at20.window_duplicates, at20.window_unique), (10, 5, 5));
        assert!((at20.window_frac_duplicates - 0.5).abs() < 1e-9);
        // Partial final window: 9 - 7 = 2 dups over 23 - 20 = 3 templates.
        let at23 = row(rec.rows(), "lib", "pairs", 23);
        assert_eq!((at23.window_total, at23.window_duplicates, at23.window_unique), (3, 2, 1));
        assert!((at23.window_frac_duplicates - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn ladder_path_appends_suffix() {
        assert_eq!(
            ladder_path(Path::new("out/sampleA")),
            Path::new("out/sampleA.duplication-sampled.tsv")
        );
    }

    // Test-only accessor for rows before finalize's filtering.
    impl LadderRecorder {
        fn rows_unfiltered(&self) -> &[DuplicateLadderRow] {
            &self.rows
        }
    }
}
