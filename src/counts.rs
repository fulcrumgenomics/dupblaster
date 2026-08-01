//! Per-library occurrence counting for the duplicate group-size histogram
//! (η_k), the `--duplication-spectrum` output.
//!
//! ## One category per library
//!
//! A library is reported on its **paired** subset if it has any both-mapped
//! pairs, and on its **single-end** signatures only if it is *solely*
//! single-end. This is both the clean statistics (a pair signature pins two
//! endpoints, so it's a far better complexity estimator than a coincidental-
//! collision-prone single coordinate) and the thing that makes the feature
//! correct under every single-end strategy: single-end counts are only ever
//! *emitted* for a library with no pairs, and "no pairs" means no pair-ends were
//! ever folded into the fragment keyspace — so the dedup dup-verdict is a
//! faithful "have I seen this fragment before" signal in that case, even under
//! the Picard strategies (whose fragment table otherwise holds pair-ends for
//! the "fragments don't beat pairs" rule).
//!
//! To avoid spending single-end memory on a library that turns out to be
//! paired, [`CountsMap`] flips a `has_pairs` flag on the first pair and **drops
//! the single-end table** at that moment (freeing any orphans seen before the
//! first pair) and skips all further single-end inserts.
//!
//! ## Side-table + keying
//!
//! Counts live in a **count side-table** ([`CountTable`]) that sits *alongside*
//! the main dedup [`crate::sig::DupTable`], mirroring its cell partitioning and
//! reusing the very same signature [`crate::sig::Slot`]s the dedup path
//! computes — so a counted signature groups exactly as the dedup verdict does,
//! by construction, without recomputing anything. Each cell is a struct-of-
//! arrays [`crate::countset::CountSet`] (u64 pair keys / u32 single-end keys).
//! Only signatures observed **≥2×** are stored (singletons are recovered by
//! subtraction, `η₁ = distinct − |side table|`), keeping it small. Per-signature
//! counts saturate at `u16::MAX`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use serde::Serialize;

use crate::countset::{CountKey, CountSet};
use crate::dedup::Stats;
use crate::sig::{FragmentSlot, PairSlot, Slot, stride_for};

/// A partitioned, struct-of-arrays count side-table that shadows the cell layout
/// of the dedup [`crate::sig::DupTable`]: one [`CountSet`] per partition cell,
/// indexed by a signature's [`Slot::off`]. Only signatures seen ≥2× are entered.
struct CountTable<K> {
    cells: Vec<CountSet<K>>,
}

impl<K: CountKey> CountTable<K> {
    /// A 2D pair side-table (`stride²` cells), matching [`crate::sig::PairDupTable`].
    fn new_pair(bin_count: u32) -> Self {
        let stride = stride_for(bin_count) as usize;
        Self::with_cells(stride.saturating_mul(stride))
    }

    /// A 1D single-end side-table (`stride` cells), matching the fragment/orphan keying.
    fn new_single_end(bin_count: u32) -> Self {
        Self::with_cells(stride_for(bin_count) as usize)
    }

    fn with_cells(n: usize) -> Self {
        Self { cells: (0..n).map(|_| CountSet::new()).collect() }
    }

    /// Record one repeat observation at `slot` (create its count at 2, else bump).
    #[inline]
    fn bump(&mut self, slot: Slot<K>) {
        self.cells[slot.off].bump(slot.sig);
    }

    /// Every held signature's occurrence count, across every cell.
    fn counts(&self) -> impl Iterator<Item = u16> + '_ {
        self.cells.iter().flat_map(|c| c.counts())
    }
}

/// Per-library occurrence counter. Tracks the distinct-signature count and a
/// count side-table (signatures seen ≥2×) for both-mapped pairs and — only while
/// the library has seen no pairs — single-end/orphan fragments. Constructed only
/// when `--duplication-spectrum` is on. Fed pre-computed [`Slot`]s from the dedup
/// path, so it holds no coordinate or strand logic of its own.
pub struct CountsMap {
    /// Bin count, retained to rebuild the single-end table when it is dropped.
    bin_count: u32,
    /// Set once the library sees its first pair; freezes single-end counting.
    has_pairs: bool,
    /// Distinct pair signatures observed (incremented on each first sighting).
    pair_distinct: u64,
    /// Pair signatures seen ≥2×, keyed by the pair [`Slot`] (2D-partitioned).
    pairs: CountTable<u64>,
    /// Distinct single-end/orphan signatures observed (only meaningful while
    /// `!has_pairs`).
    se_distinct: u64,
    /// Single-end signatures seen ≥2×, keyed by the fragment [`Slot`] (1D).
    single_end: CountTable<u32>,
}

impl CountsMap {
    /// Build an empty counter sized for `bin_count` partition cells (matching the
    /// dedup tables' stride).
    pub fn new(bin_count: u32) -> Self {
        Self {
            bin_count,
            has_pairs: false,
            pair_distinct: 0,
            pairs: CountTable::new_pair(bin_count),
            se_distinct: 0,
            single_end: CountTable::new_single_end(bin_count),
        }
    }

    /// Observe one both-mapped pair template via its precomputed pair [`Slot`]
    /// and the dedup verdict `is_dup`.
    #[inline]
    pub fn observe_pair(&mut self, slot: PairSlot, is_dup: bool) {
        if !self.has_pairs {
            // First pair: this library is paired, so its single-end counts will
            // never be emitted. Drop the table (freeing any pre-first-pair
            // orphans) and stop counting single-end from here on.
            self.has_pairs = true;
            self.single_end = CountTable::new_single_end(self.bin_count);
            self.se_distinct = 0;
        }
        if !is_dup {
            self.pair_distinct += 1;
            return;
        }
        self.pairs.bump(slot);
    }

    /// Observe one mapped single-end / orphan template via its precomputed
    /// fragment [`Slot`]. A no-op once the library has seen a pair.
    #[inline]
    pub fn observe_single_end(&mut self, slot: FragmentSlot, is_dup: bool) {
        if self.has_pairs {
            return;
        }
        if !is_dup {
            self.se_distinct += 1;
            return;
        }
        self.single_end.bump(slot);
    }

    /// η_k for pairs as an ascending `count → n_molecules` map.
    fn pair_histogram(&self) -> BTreeMap<u32, u64> {
        histogram(self.pair_distinct, self.pairs.counts())
    }

    /// η_k for single-end/orphan signatures.
    fn se_histogram(&self) -> BTreeMap<u32, u64> {
        histogram(self.se_distinct, self.single_end.counts())
    }
}

/// Fold a distinct-count and the side-table counts into an ascending
/// `count → n_molecules` histogram. `η₁ = distinct − |side table|`; every side
/// entry contributes to its own `count` bucket.
fn histogram(distinct: u64, side: impl Iterator<Item = u16>) -> BTreeMap<u32, u64> {
    let mut hist: BTreeMap<u32, u64> = BTreeMap::new();
    let mut side_len = 0u64;
    for count in side {
        *hist.entry(count as u32).or_insert(0) += 1;
        side_len += 1;
    }
    // `side_len <= distinct` by construction (every repeated signature was also
    // counted once as distinct on its first sighting); saturating for safety.
    let singletons = distinct.saturating_sub(side_len);
    if singletons > 0 {
        hist.insert(1, singletons);
    }
    hist
}

/// One row of the group-size histogram: for a `(library, category)`,
/// `n_molecules` distinct molecules were each observed `n_observations` times.
/// `sample` is first by convention.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CountHistogramRow {
    /// Sample name — `--sample` if set, else comma-joined unique `@RG SM:`.
    pub sample: String,
    /// Library name for this row.
    pub library: String,
    /// `pairs` for a library with any both-mapped pairs, else `single_end`.
    pub category: &'static str,
    /// Occurrence count `k`: how many times each distinct molecule was observed.
    /// Counts saturate at 65,535 (`u16::MAX`), so a molecule seen more often is
    /// reported at `k = 65535`; see the module docs.
    pub n_observations: u32,
    /// How many distinct molecules were observed exactly `n_observations` times.
    pub n_molecules: u64,
}

/// Build the histogram rows: one category per library — `pairs` if the library
/// has any both-mapped pairs, otherwise `single_end`. Libraries that saw no
/// data are skipped.
pub fn histogram_rows(counts: &[CountsMap], stats: &Stats, sample: &str) -> Vec<CountHistogramRow> {
    let mut rows = Vec::new();
    for (i, ls) in stats.libraries.iter().enumerate() {
        if ls.id_count == 0 {
            continue;
        }
        let cm = &counts[i];
        let category = ls.reported_category();
        let hist =
            if ls.both_mapped_id_count > 0 { cm.pair_histogram() } else { cm.se_histogram() };
        for (n_observations, n_molecules) in hist {
            rows.push(CountHistogramRow {
                sample: sample.to_string(),
                library: ls.name.clone(),
                category,
                n_observations,
                n_molecules,
            });
        }
    }
    rows
}

/// `foo/sampleA` → `foo/sampleA.duplication-spectrum.tsv`.
pub fn histogram_path(prefix: &Path) -> PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(".duplication-spectrum.tsv");
    PathBuf::from(name)
}

/// Write the histogram rows to `path` as a TSV via [`fgoxide::io::DelimFile`].
pub fn write_histogram_rows(rows: &[CountHistogramRow], path: &Path) -> Result<()> {
    DelimFile::default()
        .write_tsv(path, rows.iter())
        .with_context(|| format!("writing duplication-spectrum TSV to {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dedup::{LibraryStats, Stats};

    // Tests exercise the counting logic directly with hand-built slots; the
    // signature/partition math (and its strand handling) is tested where it
    // lives, in `sig.rs`. `off = 0` puts every signature in one cell, so a slot
    // is a duplicate iff its `sig` repeats.
    fn cm() -> CountsMap {
        CountsMap::new(0)
    }
    fn pair(sig: u64) -> PairSlot {
        Slot { off: 0, sig }
    }
    fn se(sig: u32) -> FragmentSlot {
        Slot { off: 0, sig }
    }

    #[test]
    fn repeated_pair_accumulates_full_count() {
        let mut c = cm();
        c.observe_pair(pair(100), false);
        c.observe_pair(pair(100), true);
        c.observe_pair(pair(100), true);
        assert_eq!(c.pair_distinct, 1);
        assert_eq!(c.pair_histogram().get(&3), Some(&1), "one signature seen 3×");
    }

    #[test]
    fn pair_histogram_recovers_singletons_by_subtraction() {
        let mut c = cm();
        c.observe_pair(pair(100), false);
        c.observe_pair(pair(100), true);
        c.observe_pair(pair(100), true);
        c.observe_pair(pair(500), false);
        c.observe_pair(pair(900), false);
        let h = c.pair_histogram();
        assert_eq!(h.get(&1), Some(&2), "two singletons");
        assert_eq!(h.get(&3), Some(&1), "one triple");
    }

    #[test]
    fn single_end_counted_when_library_has_no_pairs() {
        let mut c = cm();
        c.observe_single_end(se(100), false);
        c.observe_single_end(se(100), true);
        c.observe_single_end(se(500), false);
        let h = c.se_histogram();
        assert_eq!(h.get(&1), Some(&1));
        assert_eq!(h.get(&2), Some(&1));
    }

    #[test]
    fn first_pair_drops_single_end_counts_and_freezes_them() {
        let mut c = cm();
        // Orphans arrive before any pair.
        c.observe_single_end(se(100), false);
        c.observe_single_end(se(100), true);
        assert_eq!(c.single_end.counts().count(), 1);
        // First pair: single-end table is dropped and se_distinct reset.
        c.observe_pair(pair(200), false);
        assert!(c.has_pairs);
        assert_eq!(
            c.single_end.counts().count(),
            0,
            "single-end table must be freed on first pair"
        );
        assert_eq!(c.se_distinct, 0);
        // Later orphans are ignored.
        c.observe_single_end(se(700), false);
        c.observe_single_end(se(700), true);
        assert_eq!(c.single_end.counts().count(), 0);
        assert_eq!(c.se_distinct, 0);
        assert!(c.se_histogram().is_empty());
    }

    #[test]
    fn counts_saturate_at_u16_max() {
        let mut c = cm();
        c.observe_single_end(se(100), false);
        for _ in 0..70_000 {
            c.observe_single_end(se(100), true);
        }
        assert_eq!(c.se_histogram().get(&(u16::MAX as u32)), Some(&1));
    }

    #[test]
    fn histogram_rows_report_pairs_for_a_paired_library() {
        let mut counts = vec![cm()];
        // An orphan seen before the first pair — must not surface (library is paired).
        counts[0].observe_single_end(se(9), false);
        counts[0].observe_pair(pair(1), false);
        counts[0].observe_pair(pair(1), true);
        counts[0].observe_pair(pair(3), false);
        let stats = Stats {
            libraries: vec![LibraryStats {
                name: "libA".to_string(),
                id_count: 4,
                both_mapped_id_count: 3,
                both_mapped_dup_count: 1,
                mapped_orphan_id_count: 1,
                ..Default::default()
            }],
            clamped_template_count: 0,
        };
        let rows = histogram_rows(&counts, &stats, "NA12878");
        assert!(rows.iter().all(|r| r.category == "pairs"));
        assert!(rows.iter().all(|r| r.sample == "NA12878"));
        assert_eq!(rows.iter().find(|r| r.n_observations == 1).unwrap().n_molecules, 1);
        assert_eq!(rows.iter().find(|r| r.n_observations == 2).unwrap().n_molecules, 1);
    }

    #[test]
    fn histogram_rows_report_single_end_for_a_se_only_library() {
        let mut counts = vec![cm()];
        counts[0].observe_single_end(se(1), false);
        counts[0].observe_single_end(se(1), true);
        counts[0].observe_single_end(se(5), false);
        let stats = Stats {
            libraries: vec![LibraryStats {
                name: "se".to_string(),
                id_count: 3,
                mapped_orphan_id_count: 3,
                orphan_dup_count: 1,
                ..Default::default()
            }],
            clamped_template_count: 0,
        };
        let rows = histogram_rows(&counts, &stats, "");
        assert!(rows.iter().all(|r| r.category == "single_end"));
        assert_eq!(rows.iter().find(|r| r.n_observations == 1).unwrap().n_molecules, 1);
        assert_eq!(rows.iter().find(|r| r.n_observations == 2).unwrap().n_molecules, 1);
    }

    #[test]
    fn empty_library_is_skipped() {
        let counts = vec![cm()];
        let stats = Stats {
            libraries: vec![LibraryStats { name: "empty".to_string(), ..Default::default() }],
            clamped_template_count: 0,
        };
        assert!(histogram_rows(&counts, &stats, "").is_empty());
    }

    #[test]
    fn tsv_header_is_sample_first_with_expected_columns() {
        let mut counts = vec![cm()];
        counts[0].observe_single_end(se(1), false);
        let stats = Stats {
            libraries: vec![LibraryStats {
                name: "lib".to_string(),
                id_count: 1,
                mapped_orphan_id_count: 1,
                ..Default::default()
            }],
            clamped_template_count: 0,
        };
        let rows = histogram_rows(&counts, &stats, "s");
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_histogram_rows(&rows, tmp.path()).unwrap();
        let text = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(
            text.lines().next().unwrap(),
            "sample\tlibrary\tcategory\tn_observations\tn_molecules"
        );
    }

    #[test]
    fn histogram_path_appends_suffix() {
        assert_eq!(histogram_path(Path::new("out/x")), Path::new("out/x.duplication-spectrum.tsv"));
    }
}
