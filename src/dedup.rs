//! Per-block record processing — the heart of dupblaster.
//!
//! Each "block" is the run of consecutive records sharing a QNAME (dupblaster
//! requires query-grouped input). For each block we:
//!
//! 1. Identify the primary first-of-pair and primary second-of-pair (or a
//!    singleton primary if the read isn't paired).
//! 2. Compute a strand-aware 5'-aligned signature for the pair, look it up
//!    in a hash set, and if it's a repeat mark *every* record in the block
//!    as duplicate (FLAG 0x400).
//! 3. Optionally add MC/MQ mate tags to all paired records.
//! 4. Write each record to the output, either flagged-dup or, with
//!    `--removeDups`, omitted entirely.

use std::collections::HashMap;

use anyhow::{Result, bail};
use fgumi_raw_bam::RawRecord;
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

use crate::cigar::CigarInfo;
use crate::raw_writer::RawBamWriter;
use crate::sig::{
    BinIndex, FragmentDupTable, MethylationMode, PairDupTable, SingleEndStrategy,
    five_prime_aligned_pos, orphan_pos_override,
};

// ── SAM flag constants ──────────────────────────────────────────────────────

/// SAM FLAG 0x1: template has multiple segments (i.e. paired-end).
pub const FLAG_PAIRED: u16 = 0x1;
/// SAM FLAG 0x4: segment unmapped.
pub const FLAG_UNMAPPED: u16 = 0x4;
/// SAM FLAG 0x8: mate of the next segment is unmapped.
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
/// SAM FLAG 0x10: this segment is on the reverse strand.
pub const FLAG_REVERSE: u16 = 0x10;
/// SAM FLAG 0x40: first segment in the template (R1 in PE sequencing).
pub const FLAG_FIRST_SEGMENT: u16 = 0x40;
/// SAM FLAG 0x80: last segment in the template (R2 in PE sequencing).
pub const FLAG_LAST_SEGMENT: u16 = 0x80;
/// SAM FLAG 0x100: secondary alignment.
pub const FLAG_SECONDARY: u16 = 0x100;
/// SAM FLAG 0x400: PCR/optical duplicate (the bit we set).
pub const FLAG_DUPLICATE: u16 = 0x400;
/// SAM FLAG 0x800: supplementary (chimeric) alignment.
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

// ── Public types ───────────────────────────────────────────────────────────

/// Subset of CLI options the processor actually uses at runtime.
#[derive(Debug, Clone)]
pub struct ProcessorOptions {
    /// Remove dup records from output instead of just flagging them.
    pub remove_dups: bool,
    /// Add `MC` (mate CIGAR) and `MQ` (mate MAPQ) tags to paired records.
    pub add_mate_tags: bool,
    /// Don't fail on unmated alignments — see `--ignore-unmated`.
    pub ignore_unmated: bool,
    /// Width of the synthetic-genome padding around each contig. A read
    /// whose 5' clip extends more than this past a contig edge has its
    /// coordinate clamped into the contig (counted as a processing artifact).
    pub max_read_length: i32,
    /// How to key single-end / orphan reads in the dup table. See
    /// [`SingleEndStrategy`] for the design discussion.
    pub single_end_strategy: SingleEndStrategy,
    /// Methylation-aware pair keying. `None` = standard WGS (coordinate-
    /// canonical) keying. `Some(Directional)` keys pairs in template order so
    /// the two original strands of a bisulfite/enzymatic-conversion fragment at
    /// one locus stay distinct. See [`MethylationMode`]. Only the pair path is
    /// affected — the single-end/orphan path is already strand-aware, so it
    /// keeps OT/OB orphans separate in every mode.
    pub methylation_mode: Option<MethylationMode>,
}

/// Display name used for reads whose library can't be determined (no `RG`
/// tag, an `RG` absent from the header, or an `@RG` line with no `LB`).
/// Matches Picard MarkDuplicates' bucket name.
const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// Display name for the single combined bucket when library splitting is
/// turned off (`--library-unaware`) over a header that *does* declare
/// multiple libraries — signalling the merge to anyone reading `--stats`.
const ALL_READS: &str = "All Reads";

/// Maps each read to a dense library bucket, built once from the header.
///
/// Duplicates are only ever called *within* a library, so the library is part
/// of the dedup key (Picard's model; samblaster, by contrast, ignores it).
/// dupblaster realizes "library bits in the key" as an outer partition: the
/// bucket returned here selects which per-library dedup table a template is
/// checked against (see [`RecordProcessor`]).
///
/// Libraries are deduplicated by their `LB` value, **not** by `RG:ID` — two
/// read groups sharing one `LB` are one library. Bucket 0 is the catch-all
/// "unknown library". In single-library mode (`--library-unaware`, or a header
/// with ≤1 distinct `LB`) there is exactly one bucket and `lookup` is
/// never consulted — the whole mechanism is a no-op with no per-read RG scan.
pub struct LibraryIndex {
    /// `RG:ID` bytes → library bucket. Empty in single-library mode.
    rg_to_lib: HashMap<Box<[u8]>, u32>,
    /// Display name per bucket; `names.len()` is the bucket count.
    names: Vec<String>,
}

impl LibraryIndex {
    /// Build the index from the SAM header. `disabled` forces single-library
    /// mode regardless of the header (the `--library-unaware` escape hatch).
    pub fn from_header(header: &Header, disabled: bool) -> Self {
        // Collect (RG:ID, LB) pairs and the distinct LB set. BTreeSet keeps the
        // library order deterministic across runs (sorted by LB), so bucket
        // indices — and therefore the `--stats` row order — are stable.
        let mut rg_lb: Vec<(Vec<u8>, String)> = Vec::new();
        let mut distinct: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        for (id, map) in header.read_groups() {
            if let Some(lb) = map.other_fields().get(&rg_tag::LIBRARY) {
                let lb = lb.to_string();
                if !lb.is_empty() {
                    distinct.insert(lb.clone());
                    rg_lb.push((id.as_slice().to_vec(), lb));
                }
            }
        }

        // Single-table mode: splitting disabled, or the data can't be separated
        // by library because the header declares 0 or 1 distinct LB.
        if disabled || distinct.len() <= 1 {
            let name = match distinct.len() {
                1 => distinct.into_iter().next().expect("one library present"),
                0 => UNKNOWN_LIBRARY.to_string(),
                _ => ALL_READS.to_string(),
            };
            return Self { rg_to_lib: HashMap::new(), names: vec![name] };
        }

        // Library-aware: bucket 0 is the "unknown library" catch-all; buckets
        // 1..=N are the distinct libraries in sorted order.
        let mut names = Vec::with_capacity(distinct.len() + 1);
        names.push(UNKNOWN_LIBRARY.to_string());
        let mut lb_to_idx: HashMap<String, u32> = HashMap::new();
        for lb in distinct {
            let idx = names.len() as u32;
            lb_to_idx.insert(lb.clone(), idx);
            names.push(lb);
        }
        let rg_to_lib =
            rg_lb.into_iter().map(|(id, lb)| (id.into_boxed_slice(), lb_to_idx[&lb])).collect();
        Self { rg_to_lib, names }
    }

    /// Number of library buckets (≥ 1). `1` means single-library mode.
    pub fn num_libs(&self) -> u32 {
        self.names.len() as u32
    }

    /// Display name for a bucket, for `--stats` rows.
    pub fn name(&self, idx: u32) -> &str {
        &self.names[idx as usize]
    }

    /// Resolve an `RG:ID` to its library bucket; bucket 0 (unknown) when the
    /// read group isn't mapped to a library.
    fn lookup(&self, rg: &[u8]) -> u32 {
        self.rg_to_lib.get(rg).copied().unwrap_or(0)
    }
}

/// Per-library, template-level counters. In single-library mode there is one
/// of these; in library-aware mode there is one per [`LibraryIndex`] bucket.
/// Every counter is a count of QNAME blocks (templates), not alignments.
#[derive(Debug, Default, Clone)]
pub struct LibraryStats {
    /// Library display name (from [`LibraryIndex::name`]).
    pub name: String,
    /// Total templates attributed to this library.
    pub id_count: u64,
    /// Templates marked as duplicates (rollup of orphan + both-mapped dups).
    pub dup_count: u64,
    /// Templates whose primary reads are both unmapped — never dup-checked.
    pub both_unmapped_id_count: u64,
    /// Templates with a single unmapped primary: an unmapped single-end read
    /// (always), or a paired primary whose mate is absent and itself unmapped
    /// (only under `--ignore-unmated`).
    pub unmapped_orphan_id_count: u64,
    /// Templates with one mapped + one unmapped/absent primary read.
    pub mapped_orphan_id_count: u64,
    /// Duplicates among `mapped_orphan_id_count`.
    pub orphan_dup_count: u64,
    /// Templates with both primary reads mapped.
    pub both_mapped_id_count: u64,
    /// Duplicates among `both_mapped_id_count`.
    pub both_mapped_dup_count: u64,
    /// Paired records whose mate was missing (only under `--ignore-unmated`).
    pub unmated_count: u64,
}

impl LibraryStats {
    /// Add `other`'s counters into `self` (used to roll per-library counters up
    /// into a run-wide total). The name is left untouched.
    fn accumulate(&mut self, other: &LibraryStats) {
        self.id_count += other.id_count;
        self.dup_count += other.dup_count;
        self.both_unmapped_id_count += other.both_unmapped_id_count;
        self.unmapped_orphan_id_count += other.unmapped_orphan_id_count;
        self.mapped_orphan_id_count += other.mapped_orphan_id_count;
        self.orphan_dup_count += other.orphan_dup_count;
        self.both_mapped_id_count += other.both_mapped_id_count;
        self.both_mapped_dup_count += other.both_mapped_dup_count;
        self.unmated_count += other.unmated_count;
    }
}

/// Run-wide statistics — printed at exit and exposed via `--stats` TSV.
///
/// Template-level counters live in [`LibraryStats`], one per library bucket
/// (just one in single-library mode). The run-wide totals consumed by the
/// stderr summary are recovered by summing them with [`Self::totals`].
#[derive(Debug, Default, Clone)]
pub struct Stats {
    /// Per-library counters, indexed by [`LibraryIndex`] bucket.
    pub libraries: Vec<LibraryStats>,
    /// Templates with at least one read whose 5' coordinate had to be
    /// clamped into its contig because its clipping extended more than
    /// `--max-read-length` past a contig edge. A processing artifact (not a
    /// library metric): counted once per template, surfaced via stderr/log
    /// and the exit code, never written to `--stats`. Run-wide (not split by
    /// library) because it's a property of the run configuration.
    pub clamped_template_count: u64,
}

impl Stats {
    /// Build empty per-library counters named from `library_index`.
    pub fn new(library_index: &LibraryIndex) -> Self {
        let libraries = (0..library_index.num_libs())
            .map(|i| LibraryStats { name: library_index.name(i).to_string(), ..Default::default() })
            .collect();
        Self { libraries, clamped_template_count: 0 }
    }

    /// Run-wide totals across all libraries (name = "All Reads"). Used for the
    /// stderr summary and the end-of-run footer.
    pub fn totals(&self) -> LibraryStats {
        let mut total = LibraryStats { name: ALL_READS.to_string(), ..Default::default() };
        for lib in &self.libraries {
            total.accumulate(lib);
        }
        total
    }
}

/// Main per-block driver. Holds the partitioned dedup tables.
///
/// Every record in a block is either all dup-flagged or all unflagged, so we
/// only need a single `bool` of per-block state — no per-record scratch.
///
/// The dedup state is partitioned by library: `dups[lib]` / `frag_dups[lib]`
/// hold the signatures seen for library bucket `lib` (see [`LibraryIndex`]).
/// Each table is allocated lazily on first use so a header that declares many
/// libraries but only populates a few pays only for the ones it touches, and
/// the single-library common case is one table — identical to pre-library
/// behavior. Cell pre-sizing is scaled down by library count
/// (`scaled_partition_cap`) so the empty-table baseline grows ~√L rather
/// than linearly with the library count.
pub struct RecordProcessor {
    /// Bin layout for the synthetic super-contig — translates a `(tid, pos)` to
    /// `(bin_num, bin_pos)` for signature computation.
    bins: BinIndex,
    /// `bins.bin_count()` cached at construction — sizes a per-library table at
    /// its first allocation. Cached so the hot-path lazy-table accessors don't
    /// recompute it per template.
    bin_count: u32,
    /// Per-library pair tables, indexed by library bucket; `None` until the
    /// library's first pair is seen.
    dups: Vec<Option<PairDupTable>>,
    /// Per-library fragment tables holding the 5' position of every paired
    /// read end, so a later orphan / single-end read at the same
    /// `(bin_num, bin_pos, strand)` is detected as a duplicate of the pair via
    /// [`FragmentDupTable::check_or_insert`]. Two strategies populate it, at
    /// different times:
    ///
    /// * [`SingleEndStrategy::PicardApprox`]: filled *concurrently* with the
    ///   pair pass — each pair-end is inserted as its pair is seen.
    ///   Approximate, because an orphan that precedes its pair in the stream
    ///   isn't yet registered.
    /// * [`SingleEndStrategy::PicardExact`]: built *after* the pair pass by
    ///   [`Self::finalize_fragment_table`] draining each library's pair table,
    ///   so it's complete before any deferred fragment is checked.
    ///
    /// Stays all-`None` under `StrandAware` / `SamblasterLegacy` (those key
    /// orphans on the pair table's `s1 = 0` row). Uses `u32`-keyed cells (vs
    /// the pair table's `u64`) and a 1D `stride`-cell layout (vs 2D `stride²`),
    /// keeping the memory overhead bounded.
    frag_dups: Vec<Option<FragmentDupTable>>,
    /// Resolves each template's read group to a library bucket. `num_libs() ==
    /// 1` means single-library mode, where [`Self::resolve_library`] short-
    /// circuits to bucket 0 without ever scanning a record's `RG` tag.
    library_index: LibraryIndex,
    /// Per-cell initial capacity for lazily-built tables, scaled down by the
    /// library count so the baseline doesn't grow linearly. See
    /// [`scaled_partition_cap`].
    partition_cap: usize,
    /// One-entry cache for [`Self::resolve_library`]: the last `RG` bytes seen
    /// and the bucket they resolved to. Real inputs have long runs of one read
    /// group, so this collapses the per-template lookup to a `memcmp` almost
    /// always. Empty until the first RG is resolved.
    last_rg: Vec<u8>,
    /// Library bucket corresponding to `last_rg`. Updated in tandem.
    last_lib_idx: u32,
    /// Configuration options in effect for this run — remove-dups, mate-tags,
    /// methylation mode, etc.
    opts: ProcessorOptions,
    /// `add_mate_tags` writes mate CIGAR into here so we don't allocate a
    /// fresh `String` per pair. Reused across blocks.
    mate_cigar_scratch: Vec<u8>,
}

impl RecordProcessor {
    /// Build from a flat slice of reference-sequence lengths (in BAM order,
    /// 1-indexed — `ref_lengths[0]` is contig at `tid=0`). `library_index`
    /// determines how many per-library dedup tables there are (one in
    /// single-library mode).
    pub fn from_ref_lengths(
        ref_lengths: &[i32],
        opts: ProcessorOptions,
        min_bin_count: u32,
        library_index: LibraryIndex,
    ) -> Self {
        let bins = BinIndex::from_ref_lengths(ref_lengths, opts.max_read_length, min_bin_count);
        let num_libs = library_index.num_libs() as usize;
        let partition_cap = scaled_partition_cap(library_index.num_libs());
        // Cache `bin_count` once: it's needed only to size a table at first
        // allocation, but the lazy-table accessors are on the hot path, so
        // recomputing it per call would be wasted arithmetic on every template.
        let bin_count = bins.bin_count();
        // Tables are allocated lazily (on first use), so every library starts
        // `None`; the single-library case allocates exactly one table on its
        // first template, matching the pre-library memory profile.
        let dups = (0..num_libs).map(|_| None).collect();
        let frag_dups = (0..num_libs).map(|_| None).collect();
        Self {
            bins,
            bin_count,
            dups,
            frag_dups,
            library_index,
            partition_cap,
            last_rg: Vec::new(),
            last_lib_idx: 0,
            opts,
            mate_cigar_scratch: Vec::with_capacity(64),
        }
    }

    /// Lazily get (allocating on first use) the pair table for `lib`.
    fn pair_table(&mut self, lib: u32) -> &mut PairDupTable {
        let (bin_count, cap) = (self.bin_count, self.partition_cap);
        self.dups[lib as usize].get_or_insert_with(|| PairDupTable::new_pair(bin_count, cap))
    }

    /// Lazily get (allocating on first use) the fragment table for `lib`.
    fn frag_table(&mut self, lib: u32) -> &mut FragmentDupTable {
        let (bin_count, cap) = (self.bin_count, self.partition_cap);
        self.frag_dups[lib as usize]
            .get_or_insert_with(|| FragmentDupTable::new_single_end(bin_count, cap))
    }

    /// Resolve a template's library bucket from the read group on its first
    /// record. Short-circuits to bucket 0 in single-library mode (no RG scan).
    /// A record with *no* `RG` tag maps to the unknown bucket (0) directly,
    /// without touching the one-entry cache; an `RG` that's present but maps to
    /// the unknown bucket (e.g. absent from the header) still caches normally so
    /// repeats stay on the fast path.
    fn resolve_library(&mut self, block: &[RawRecord]) -> u32 {
        if self.library_index.num_libs() == 1 {
            return 0;
        }
        // All records of a QNAME block share a read group in practice, so the
        // first record is a faithful (and cheapest) source.
        match block[0].tags().find_string(b"RG") {
            Some(rg) => {
                if !self.last_rg.is_empty() && self.last_rg == rg {
                    return self.last_lib_idx;
                }
                let idx = self.library_index.lookup(rg);
                self.last_rg.clear();
                self.last_rg.extend_from_slice(rg);
                self.last_lib_idx = idx;
                idx
            }
            None => 0,
        }
    }

    /// Expose the chosen `bin_shift` for diagnostic logging.
    pub fn bin_shift(&self) -> u32 {
        self.bins.bin_shift()
    }

    /// Expose `bin_count` for diagnostic logging.
    pub fn bin_count(&self) -> u32 {
        self.bin_count
    }

    /// Process one block of records sharing a QNAME and emit them to the output.
    pub fn process_block(
        &mut self,
        block: &mut [RawRecord],
        stats: &mut Stats,
        out: &mut RawBamWriter,
    ) -> Result<()> {
        let lib = self.resolve_library(block);
        let is_dup = self.mark_dups(block, lib, stats)?;
        self.emit(block, is_dup, out)?;
        Ok(())
    }

    /// Classify `block`, update `stats`, and return `true` if the template is a
    /// duplicate. Does not write output — the caller (e.g. [`Self::process_block`])
    /// passes the result to [`Self::emit`].
    fn mark_dups(&mut self, block: &mut [RawRecord], lib: u32, stats: &mut Stats) -> Result<bool> {
        let class = self.classify_block(block)?;
        let li = lib as usize;
        stats.libraries[li].id_count += 1;
        Ok(match class {
            BlockClass::BothUnmapped => {
                stats.libraries[li].both_unmapped_id_count += 1;
                false
            }
            BlockClass::UnmappedOrphan => {
                stats.libraries[li].unmapped_orphan_id_count += 1;
                false
            }
            BlockClass::Unmated => {
                stats.libraries[li].unmated_count += 1;
                false
            }
            BlockClass::Pair { first, second } => {
                if self.opts.add_mate_tags {
                    self.add_mate_tags_pair(block, first, second);
                }
                stats.libraries[li].both_mapped_id_count += 1;
                let dup = self.check_pair_signature(block, first, second, lib, stats)?;
                if dup {
                    stats.libraries[li].dup_count += 1;
                    stats.libraries[li].both_mapped_dup_count += 1;
                }
                dup
            }
            BlockClass::Fragment { mapped, mate } => {
                if self.opts.add_mate_tags
                    && let Some(m) = mate
                {
                    self.add_mate_tags_pair(block, mapped, m);
                }
                stats.libraries[li].mapped_orphan_id_count += 1;
                let dup = self.check_fragment_signature(block, mapped, lib, stats)?;
                if dup {
                    stats.libraries[li].dup_count += 1;
                    stats.libraries[li].orphan_dup_count += 1;
                }
                dup
            }
        })
    }

    /// Identify the primary alignments in a QNAME block and classify it,
    /// without mutating any state (no stats, no dup tables, no mate tags).
    ///
    /// This is the single source of truth for "what kind of template is
    /// this?": it is shared by the single-pass [`Self::mark_dups`] and by the
    /// picard-exact two-pass driver ([`Self::process_block_phase1`] /
    /// [`Self::process_fragment_block`]), so the two paths can never disagree
    /// about which blocks are pairs vs fragments. Bails on the same
    /// non-query-grouped / broken-block conditions the single-pass path did.
    fn classify_block(&self, block: &[RawRecord]) -> Result<BlockClass> {
        let mut first: Option<usize> = None;
        let mut second: Option<usize> = None;
        // Read each record's flag word exactly once per iteration. The
        // compiler can usually CSE multiple `rec.flags()` calls but
        // hoisting it explicitly removes the dependency on inlining
        // decisions and matches what we'd write in C.
        for (i, rec) in block.iter().enumerate() {
            let f = rec.flags();
            // Primary = neither secondary nor supplementary. One mask + cmp.
            if f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
                continue;
            }
            if f & FLAG_PAIRED == 0 {
                // Single-end primary: there is no first/second distinction, so
                // park it in `second`. The lone-primary branch below retrieves
                // it via `first.or(second)` regardless of which slot is used.
                second = Some(i);
            } else if f & FLAG_FIRST_SEGMENT != 0 {
                first = Some(i);
            } else if f & FLAG_LAST_SEGMENT != 0 {
                second = Some(i);
            }
        }

        if first.is_none() && second.is_none() {
            // The block contains records but none of them are primary
            // alignments — i.e. it has only secondary and/or supplementary
            // records. This nearly always means the input is not
            // query-grouped (the primary alignment for this QNAME is
            // somewhere else in the stream, typically because the file is
            // coordinate-sorted). Fail loudly regardless of --ignore-unmated;
            // this condition is a configuration bug, not a tolerable edge case.
            bail!("{}", non_query_grouped_message(block));
        }

        if first.is_none() || second.is_none() {
            let only_idx = first.or(second).expect("at least one primary present");
            let only_flags = block[only_idx].flags();

            // Single-end (unpaired) read: there is no mate, so none of the
            // broken-pair / not-query-grouped conditions apply — the template
            // is already complete. An unmapped SE read simply passes through
            // untouched (counted as an unmapped orphan, never dup-checked); a
            // mapped SE read enters fragment dedup. This must come before the
            // paired-orphan handling below, which gates a *missing mate* behind
            // --ignore-unmated; an SE read never had a mate to miss.
            if !has(only_flags, FLAG_PAIRED) {
                if has(only_flags, FLAG_UNMAPPED) {
                    return Ok(BlockClass::UnmappedOrphan);
                }
                return Ok(BlockClass::Fragment { mapped: only_idx, mate: None });
            }

            // From here the read IS paired but its mate is absent from the
            // block. Check UNMAPPED first: a paired+unmapped lone primary is an
            // unmapped orphan, not an "unmated" record. (Both conditions hold
            // but the more specific one wins.)
            if has(only_flags, FLAG_UNMAPPED) {
                if self.opts.ignore_unmated {
                    return Ok(BlockClass::UnmappedOrphan);
                }
                bail!("{}", broken_block_message(block));
            }
            // Paired primary, mate flag says mate is *mapped*, but the mate
            // is missing from this block → input is not properly QNAME-grouped.
            if !has(only_flags, FLAG_MATE_UNMAPPED) {
                if self.opts.ignore_unmated {
                    return Ok(BlockClass::Unmated);
                }
                bail!("{}", broken_block_message(block));
            }
            // A paired primary that is mapped, with its mate flagged unmapped
            // but absent from the block (e.g. unmapped mate filtered out
            // upstream). A lone mapped orphan with no mate to tag.
            return Ok(BlockClass::Fragment { mapped: only_idx, mate: None });
        }

        // Both primaries present.
        let f = first.expect("both primaries present");
        let s = second.expect("both primaries present");
        let ff = block[f].flags();
        let sf = block[s].flags();
        if has(ff, FLAG_UNMAPPED) && has(sf, FLAG_UNMAPPED) {
            return Ok(BlockClass::BothUnmapped);
        }
        // One end unmapped → the mapped end is a fragment; the unmapped mate
        // is present in-block, so it's available for mate tagging.
        if has(ff, FLAG_UNMAPPED) {
            return Ok(BlockClass::Fragment { mapped: s, mate: Some(f) });
        }
        if has(sf, FLAG_UNMAPPED) {
            return Ok(BlockClass::Fragment { mapped: f, mate: Some(s) });
        }
        Ok(BlockClass::Pair { first: f, second: s })
    }

    /// Compute the doubly-mapped pair signature, record it, and return
    /// whether the pair is a duplicate. Under `PicardApprox` this also
    /// registers both pair ends in the concurrent fragment table; under
    /// `PicardExact` the fragment table is built later by draining
    /// [`Self::finalize_fragment_table`], so there's nothing to register here.
    fn check_pair_signature(
        &mut self,
        block: &[RawRecord],
        first: usize,
        second: usize,
        lib: u32,
        stats: &mut Stats,
    ) -> Result<bool> {
        let f_derived = self.derive_for_dup(&block[first], false)?;
        let s_derived = self.derive_for_dup(&block[second], false)?;
        if f_derived.clamped || s_derived.clamped {
            self.note_clamp(stats);
        }
        // Directional methylation mode keys in template order (first-of-pair →
        // slot A) so the two original strands (OT/OB) at a locus produce
        // distinct keys — coordinate canonicalization would merge them. The
        // per-end strand bits in the cell index carry orientation either way,
        // so dropping the swap is the whole change. Standard WGS keying stays
        // coordinate-canonical (Picard/samblaster-exact). The `matches!` (not
        // `is_some`) is deliberate: a future non-directional/PBAT mode must NOT
        // take this branch — it needs canonical positions plus a strand bit.
        let template_order =
            matches!(self.opts.methylation_mode, Some(MethylationMode::Directional));
        let (a, b) = if !template_order && need_swap(&f_derived, &s_derived) {
            (&s_derived, &f_derived)
        } else {
            (&f_derived, &s_derived)
        };
        let is_dup = self.pair_table(lib).check_dm(
            a.bin_num,
            a.bin_pos,
            a.is_reverse,
            b.bin_num,
            b.bin_pos,
            b.is_reverse,
        );
        // PicardApprox side effect: register both ends of the pair in this
        // library's fragment table so later orphans at those coordinates are
        // marked as duplicates of this pair ("fragments don't beat pairs",
        // streaming approximation). PicardExact builds its fragment tables
        // later by draining the pair tables; the other strategies have none.
        if self.opts.single_end_strategy == SingleEndStrategy::PicardApprox {
            // Strand-aware keying always; FragmentDupTable doesn't expose the
            // SamblasterLegacy strand-drop variant.
            let frag = self.frag_table(lib);
            let _ =
                frag.check_or_insert(f_derived.bin_num, f_derived.bin_pos, f_derived.is_reverse);
            let _ =
                frag.check_or_insert(s_derived.bin_num, s_derived.bin_pos, s_derived.is_reverse);
        }
        Ok(is_dup)
    }

    /// Compute the single-end / orphan signature for the mapped fragment end
    /// and return whether it's a duplicate, against library `lib`'s tables.
    /// `PicardApprox` / `PicardExact` consult the fragment table (which also
    /// holds every paired read end, so an orphan collides with a pair sharing
    /// its 5' coordinate); `StrandAware` / `SamblasterLegacy` use the orphan
    /// row of the pair table.
    fn check_fragment_signature(
        &mut self,
        block: &[RawRecord],
        mapped: usize,
        lib: u32,
        stats: &mut Stats,
    ) -> Result<bool> {
        let d = self.derive_for_dup(&block[mapped], true)?;
        if d.clamped {
            self.note_clamp(stats);
        }
        Ok(match self.opts.single_end_strategy {
            SingleEndStrategy::PicardApprox | SingleEndStrategy::PicardExact => {
                self.frag_table(lib).check_or_insert(d.bin_num, d.bin_pos, d.is_reverse)
            }
            strategy @ (SingleEndStrategy::StrandAware | SingleEndStrategy::SamblasterLegacy) => {
                self.pair_table(lib).check_orphan(d.bin_num, d.bin_pos, d.is_reverse, strategy)
            }
        })
    }

    // ── picard-exact two-pass driver ────────────────────────────────────────

    /// Phase 1 of [`SingleEndStrategy::PicardExact`]: emit pairs and
    /// unmapped/unmated blocks straight to `out` exactly as the single-pass
    /// path would, but *defer* every mapped-fragment block by writing its
    /// whole record group to `temp` (an uncompressed BAM). Deferred blocks are
    /// neither dup-marked nor counted here — that happens in phase 2, after
    /// the fragment table is complete. Mate tags are added to fragment blocks
    /// before buffering so the buffered copy matches single-pass output.
    pub fn process_block_phase1(
        &mut self,
        block: &mut [RawRecord],
        stats: &mut Stats,
        out: &mut RawBamWriter,
        temp: &mut RawBamWriter,
    ) -> Result<()> {
        let lib = self.resolve_library(block);
        let li = lib as usize;
        match self.classify_block(block)? {
            BlockClass::Fragment { mapped, mate } => {
                if self.opts.add_mate_tags
                    && let Some(m) = mate
                {
                    self.add_mate_tags_pair(block, mapped, m);
                }
                // Deferred to phase 2 (re-read from the temp BAM); the RG tag
                // round-trips, so library attribution happens there.
                for rec in block.iter() {
                    temp.write_record(rec)?;
                }
                Ok(())
            }
            BlockClass::BothUnmapped => {
                stats.libraries[li].id_count += 1;
                stats.libraries[li].both_unmapped_id_count += 1;
                self.emit(block, false, out)
            }
            BlockClass::UnmappedOrphan => {
                stats.libraries[li].id_count += 1;
                stats.libraries[li].unmapped_orphan_id_count += 1;
                self.emit(block, false, out)
            }
            BlockClass::Unmated => {
                stats.libraries[li].id_count += 1;
                stats.libraries[li].unmated_count += 1;
                self.emit(block, false, out)
            }
            BlockClass::Pair { first, second } => {
                stats.libraries[li].id_count += 1;
                if self.opts.add_mate_tags {
                    self.add_mate_tags_pair(block, first, second);
                }
                stats.libraries[li].both_mapped_id_count += 1;
                let dup = self.check_pair_signature(block, first, second, lib, stats)?;
                if dup {
                    stats.libraries[li].dup_count += 1;
                    stats.libraries[li].both_mapped_dup_count += 1;
                }
                self.emit(block, dup, out)
            }
        }
    }

    /// Transition between phases 1 and 2 of [`SingleEndStrategy::PicardExact`]:
    /// consume each library's completed pair table into its fragment table,
    /// which then holds the 5' position of every paired read end for that
    /// library. Call exactly once, after phase 1 and before phase 2.
    pub fn finalize_fragment_table(&mut self) {
        debug_assert!(
            self.frag_dups.iter().all(Option::is_none),
            "fragment tables must be built exactly once"
        );
        let cap = self.partition_cap;
        for lib in 0..self.dups.len() {
            if let Some(pair) = self.dups[lib].as_mut() {
                self.frag_dups[lib] = Some(pair.drain_into_fragment_table(cap));
            }
        }
    }

    /// Phase 2 of [`SingleEndStrategy::PicardExact`]: process one buffered
    /// fragment block (re-read from the temp BAM) against the now-complete
    /// fragment table, marking it duplicate if its 5' coordinate matches any
    /// pair end or earlier fragment, then emit it. Fragment-related stats are
    /// accrued here. Requires [`Self::finalize_fragment_table`] to have run.
    pub fn process_fragment_block(
        &mut self,
        block: &mut [RawRecord],
        stats: &mut Stats,
        out: &mut RawBamWriter,
    ) -> Result<()> {
        let lib = self.resolve_library(block);
        let li = lib as usize;
        stats.libraries[li].id_count += 1;
        // The temp BAM holds only blocks classified as Fragment in phase 1,
        // and a BAM round-trip preserves flags, so re-classification yields
        // Fragment again. Anything else is a bug in the buffering logic.
        let dup = match self.classify_block(block)? {
            BlockClass::Fragment { mapped, .. } => {
                stats.libraries[li].mapped_orphan_id_count += 1;
                let d = self.check_fragment_signature(block, mapped, lib, stats)?;
                if d {
                    stats.libraries[li].dup_count += 1;
                    stats.libraries[li].orphan_dup_count += 1;
                }
                d
            }
            other => bail!(
                "picard-exact phase 2 encountered a non-fragment block ({}); \
                 this is an internal error in fragment buffering",
                other.describe()
            ),
        };
        self.emit(block, dup, out)
    }

    /// Add `MC`/`MQ` mate tags in both directions for a block whose two
    /// primary ends are at indices `a` and `b`.
    fn add_mate_tags_pair(&mut self, block: &mut [RawRecord], a: usize, b: usize) {
        self.add_mate_tags(block, a, b);
        self.add_mate_tags(block, b, a);
    }

    /// Compute the [`DerivedAlignment`] (bin, position, strand, clamp flag) for
    /// a single record. `orphan` is `true` for single-end / mapped-orphan records,
    /// which may have their position overridden under `SamblasterLegacy`.
    fn derive_for_dup(&self, rec: &RawRecord, orphan: bool) -> Result<DerivedAlignment> {
        let info = CigarInfo::from_cigar_ops(rec.cigar_ops_iter());
        let rapos = rapos_of(rec);
        let is_reverse = has(rec.flags(), FLAG_REVERSE);
        let mut pos =
            five_prime_aligned_pos(rapos, info.sclip, info.eclip, info.ra_len, is_reverse);
        // Under the SamblasterLegacy strategy, override the reverse-strand
        // 5' coord with the leftmost-aligned coord so fwd/rev orphans
        // whose alignments share a leftmost position collide.
        // StrandAware and PicardApprox keep the regular strand-aware 5'
        // position.
        if orphan
            && is_reverse
            && self.opts.single_end_strategy == SingleEndStrategy::SamblasterLegacy
        {
            pos = orphan_pos_override(rapos, info.sclip);
        }
        let seq_num = seq_num_of(rec);
        // Guard against a malformed BAM whose record references a contig id
        // beyond the header's @SQ list; without this the `seq_offs` lookup
        // in `bin_for` would panic with an opaque index-out-of-bounds.
        if seq_num >= self.bins.num_seqs() {
            bail!(
                "record {} references reference id {} but the header declares only {} \
                 reference sequence(s) — input BAM is malformed",
                String::from_utf8_lossy(rec.read_name()),
                rec.ref_id(),
                self.bins.num_seqs() - 1,
            );
        }
        // A clip extending more than --max-read-length past a contig edge is
        // clamped into the contig; `clamped` is folded into a per-template
        // counter by the caller (see `note_clamp`).
        let (bin_num, bin_pos, clamped) = self.bins.bin_for(seq_num, pos);
        Ok(DerivedAlignment { bin_num, bin_pos, is_reverse, seq_num, pos, clamped })
    }

    /// Record that a template had at least one end's coordinate clamped at a
    /// contig edge: warn once (the first time, mid-stream) and bump the
    /// per-template counter. Call at most once per template.
    fn note_clamp(&self, stats: &mut Stats) {
        if stats.clamped_template_count == 0 {
            log::warn!(
                "A read's 5' coordinate was clamped to its contig: its clipping extends more \
                 than --max-read-length ({}) bases past a contig edge, so duplicate marking may \
                 be imprecise for such reads. Re-run with a larger --max-read-length to avoid it.",
                self.opts.max_read_length,
            );
        }
        stats.clamped_template_count += 1;
    }

    /// Write `MC` and `MQ` tags onto every record in `block` that belongs to
    /// the same end as `target_first` (matched by `FLAG_FIRST_SEGMENT`), using
    /// the alignment at index `mate` as the source. Skips if the mate is unmapped.
    fn add_mate_tags(&mut self, block: &mut [RawRecord], target_first: usize, mate: usize) {
        let mate_flags = block[mate].flags();
        if has(mate_flags, FLAG_UNMAPPED) {
            return;
        }
        // Render the mate's CIGAR into our scratch Vec once, instead of
        // calling `cigar_to_string()` (which allocates N+1 Strings per call
        // — one per op via `to_string()`).
        self.mate_cigar_scratch.clear();
        write_cigar_text(block[mate].cigar_ops_iter(), &mut self.mate_cigar_scratch);
        let mq = block[mate].mapq();
        let target_first_bit = has(block[target_first].flags(), FLAG_FIRST_SEGMENT);
        for rec in block.iter_mut() {
            if has(rec.flags(), FLAG_FIRST_SEGMENT) != target_first_bit {
                continue;
            }
            let has_mc = rec.tags().find_string(b"MC").is_some();
            // `find_int` matches any of c/C/s/S/i/I subtypes. We must use it
            // (not `find_uint8`) because fgumi's `append_int` emits `c` for
            // MAPQ values 0..=127 — virtually all real MAPQ scores. Using
            // `find_uint8` (only matches `C`) would miss those and append a
            // duplicate `MQ` tag on every run-through.
            let has_mq = rec.tags().find_int(b"MQ").is_some();
            let mut editor = rec.tags_editor();
            if !has_mc {
                editor.append_string(b"MC", &self.mate_cigar_scratch);
            }
            if !has_mq {
                editor.append_int(b"MQ", i32::from(mq));
            }
        }
    }

    /// Write every record in `block` to `out`, overwriting each record's
    /// `FLAG_DUPLICATE` bit to match `is_dup`. If `remove_dups` is set and the
    /// block is a duplicate, the records are dropped entirely.
    fn emit(
        &mut self,
        block: &mut [RawRecord],
        is_dup: bool,
        out: &mut RawBamWriter,
    ) -> Result<()> {
        for rec in block.iter_mut() {
            // Always overwrite the duplicate flag with this run's decision.
            // Picard MarkDuplicates does the same — any pre-existing
            // FLAG_DUPLICATE bit (e.g. legacy markings on 1000G archive
            // BAMs from earlier Picard runs) is cleared before we set our
            // own. Without this, re-running dup-marking on already-marked
            // BAMs would union the two passes' markings.
            let f = rec.flags() & !FLAG_DUPLICATE;
            rec.set_flags(if is_dup { f | FLAG_DUPLICATE } else { f });
            if !(self.opts.remove_dups && is_dup) {
                out.write_record(rec)?;
            }
        }
        Ok(())
    }
}

/// Classification of a QNAME block by its primary alignments, produced by
/// [`RecordProcessor::classify_block`]. Carries the record indices the
/// downstream marking logic needs, but no mutable state.
#[derive(Debug, Clone, Copy)]
enum BlockClass {
    /// Both primary ends are unmapped — never dup-checked.
    BothUnmapped,
    /// A lone unmapped primary: an unmapped single-end read (always), or a
    /// paired primary whose mate is absent and itself unmapped (under
    /// `--ignore-unmated`). Emitted untouched, never dup-checked.
    UnmappedOrphan,
    /// `--ignore-unmated`: a paired primary whose mate should be mapped but is
    /// absent from the block (unmated / not properly query-grouped).
    Unmated,
    /// Both primary ends are mapped — a full pair.
    Pair {
        /// Block index of the primary first-of-pair record.
        first: usize,
        /// Block index of the primary second-of-pair record.
        second: usize,
    },
    /// Exactly one mapped end participates in fragment dedup: a single-end
    /// read or a mapped orphan. `mapped` is its index; `mate` is `Some(idx)`
    /// when the unmapped mate is present in the block (used only for mate
    /// tagging), or `None` for a lone primary.
    Fragment {
        /// Block index of the mapped primary (the end that enters the dup table).
        mapped: usize,
        /// Block index of the unmapped mate, if it is present in this block.
        mate: Option<usize>,
    },
}

impl BlockClass {
    /// Short human-readable label, used only in internal-error messages.
    fn describe(&self) -> &'static str {
        match self {
            BlockClass::BothUnmapped => "both-unmapped",
            BlockClass::UnmappedOrphan => "unmapped-orphan",
            BlockClass::Unmated => "unmated",
            BlockClass::Pair { .. } => "pair",
            BlockClass::Fragment { .. } => "fragment",
        }
    }
}

/// Dedup-relevant properties derived from a single alignment record, produced
/// by [`RecordProcessor::derive_for_dup`]. Carries just enough info to
/// compute and look up the duplicate signature without keeping the full record
/// in scope.
#[derive(Debug, Clone, Copy)]
struct DerivedAlignment {
    /// Which bin (partition cell dimension) this alignment's super-contig
    /// position falls into.
    bin_num: u32,
    /// Position within the bin, packed as the low 32 bits of the pair sig.
    bin_pos: u32,
    /// Whether the alignment is on the reverse strand.
    is_reverse: bool,
    /// 1-based sequence number in the super-contig (`tid + 1`; 0 = unmapped).
    seq_num: usize,
    /// Strand-aware 5'-aligned reference position (before binning).
    pos: i32,
    /// True if [`BinIndex::bin_for`] had to clamp this end's coordinate into
    /// its contig (clip exceeded `--max-read-length` past a contig edge).
    clamped: bool,
}

/// Return `true` if all bits in `bit` are set in `flags`. Used throughout
/// `dedup.rs` as a concise alternative to the inline `flags & bit != 0` pattern.
#[inline]
fn has(flags: u16, bit: u16) -> bool {
    flags & bit != 0
}

/// Convert htslib's 0-based pos (and -1-for-unmapped) to samblaster's
/// 1-based rapos.
#[inline]
fn rapos_of(rec: &RawRecord) -> i32 {
    let p = rec.pos();
    if p < 0 { 0 } else { p + 1 }
}

/// Convert ref_id (-1 = unmapped) to samblaster's 1-based seq_num where 0 is
/// reserved for `"*"`.
#[inline]
fn seq_num_of(rec: &RawRecord) -> usize {
    let tid = rec.ref_id();
    if tid < 0 { 0 } else { (tid + 1) as usize }
}

/// Should we swap `a` and `b` to canonicalize the pair before signature
/// computation? The desired order is `(pos asc, seq_num asc, forward strand
/// first)`. Equivalent lexicographic compare on the canonical tuple — we
/// want the "smaller" alignment in slot a, so swap iff a > b.
fn need_swap(a: &DerivedAlignment, b: &DerivedAlignment) -> bool {
    (a.pos, a.seq_num, a.is_reverse) > (b.pos, b.seq_num, b.is_reverse)
}

/// Base per-cell initial capacity for a single-library run (the value tuned
/// by profiling — see [`crate::sig::DupTable::with_n_cells`]).
const BASE_PARTITION_CAP: usize = 64;

/// Floor on the scaled per-cell capacity so pre-sizing still buys something at
/// large library counts.
const MIN_PARTITION_CAP: usize = 8;

/// Per-cell initial capacity for the dedup tables, scaled down for multi-
/// library runs.
///
/// Partitioning the dedup state into one table per library multiplies the
/// empty-cell baseline (`stride²` cells × cap) by the number of tables. To keep
/// that baseline from growing linearly with the library count, we divide the
/// base capacity by `ceil(√num_libs)`. The square root (rather than dividing by
/// `num_libs` outright) is deliberate: it leaves each cell large enough to
/// absorb a *skewed-dominant* library — where one library holds most of the
/// reads — without excessive rehashing, while still keeping the total baseline
/// growth sub-linear (≈√L). At `num_libs == 1` the divisor is 1, so the cap is
/// the unscaled base and single-library runs are byte-for-byte unchanged.
fn scaled_partition_cap(num_libs: u32) -> usize {
    let divisor = isqrt_ceil(num_libs.max(1)) as usize;
    (BASE_PARTITION_CAP / divisor.max(1)).max(MIN_PARTITION_CAP)
}

/// `ceil(sqrt(n))` via integer square root (no float rounding surprises at
/// perfect squares).
fn isqrt_ceil(n: u32) -> u32 {
    if n == 0 {
        return 0;
    }
    let s = n.isqrt();
    if s * s == n { s } else { s + 1 }
}

/// Write a CIGAR string (e.g. `"5S45M5S"`) into `out` from packed BAM ops.
/// Avoids the per-op `String` allocation that `cigar_to_string()` performs.
fn write_cigar_text<I: IntoIterator<Item = u32>>(ops: I, out: &mut Vec<u8>) {
    use std::io::Write as _;
    for word in ops {
        let len = word >> 4;
        let code = (word & 0xf) as usize;
        // Same op-code table as `parse_cigar_ops` in sam_reader.rs.
        const OPS: [u8; 9] = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];
        // itoa would be faster, but `write!` into a Vec<u8> uses no
        // intermediate heap allocation — the lookup-table-based digit
        // path in core::fmt fills directly into our buffer.
        let _ = write!(out, "{len}");
        // Valid BAM only uses op codes 0..=8; a higher code means corrupt
        // input. `'?'` keeps us from panicking on it, but flag it in debug.
        debug_assert!(code < OPS.len(), "invalid CIGAR op code {code}");
        out.push(*OPS.get(code).unwrap_or(&b'?'));
    }
    if out.is_empty() {
        out.push(b'*');
    }
}

/// Error shown when a paired read is missing its mate or appears without a
/// proper primary alignment, indicating the input is not query-grouped.
fn broken_block_message(block: &[RawRecord]) -> String {
    let qname = block.first().map(|r| r.read_name().to_vec()).unwrap_or_default();
    let qname_str = String::from_utf8_lossy(&qname);
    format!(
        "Can't find first and/or second of pair in sam block of length {} for id: {}\n\
         dupblaster: Are you sure the input is sorted by read ids?",
        block.len(),
        qname_str
    )
}

/// Error when a QNAME block contains records but no primary alignment.
/// This is the signature of non-query-grouped input — typically a
/// coordinate-sorted BAM where the primary lives at a different position
/// from the secondary/supplementary alignments we just collected.
fn non_query_grouped_message(block: &[RawRecord]) -> String {
    let qname = block.first().map(|r| r.read_name().to_vec()).unwrap_or_default();
    let qname_str = String::from_utf8_lossy(&qname);
    format!(
        "QNAME {} appeared with {} secondary/supplementary record(s) but no \
         primary — the primary must be in a different part of the stream. \
         This is almost certainly because the input is not query-grouped \
         (e.g. coordinate-sorted). Re-sort with `samtools sort -n` or \
         `mako sort --queryname` and re-run.",
        qname_str,
        block.len(),
    )
}

#[cfg(test)]
mod tests {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReadGroup;

    use super::*;

    /// Build a header with the given `(RG:ID, LB)` read groups (LB optional).
    fn header_with_read_groups(rgs: &[(&str, Option<&str>)]) -> Header {
        let mut header = Header::default();
        for (id, lb) in rgs {
            let mut rg = Map::<ReadGroup>::default();
            if let Some(lb) = lb {
                rg.other_fields_mut().insert(rg_tag::LIBRARY, (*lb).into());
            }
            header.read_groups_mut().insert((*id).into(), rg);
        }
        header
    }

    #[test]
    fn isqrt_ceil_rounds_up_to_the_next_integer_root() {
        assert_eq!(isqrt_ceil(0), 0);
        assert_eq!(isqrt_ceil(1), 1);
        assert_eq!(isqrt_ceil(2), 2);
        assert_eq!(isqrt_ceil(4), 2);
        assert_eq!(isqrt_ceil(5), 3);
        assert_eq!(isqrt_ceil(9), 3);
        assert_eq!(isqrt_ceil(10), 4);
    }

    #[test]
    fn scaled_partition_cap_is_unscaled_for_a_single_library() {
        // The common case must be byte-for-byte identical to pre-library runs.
        assert_eq!(scaled_partition_cap(1), BASE_PARTITION_CAP);
    }

    #[test]
    fn scaled_partition_cap_shrinks_by_ceil_sqrt_of_library_count() {
        assert_eq!(scaled_partition_cap(4), BASE_PARTITION_CAP / 2); // ceil(√4)=2
        assert_eq!(scaled_partition_cap(5), BASE_PARTITION_CAP / 3); // ceil(√5)=3
        assert_eq!(scaled_partition_cap(9), BASE_PARTITION_CAP / 3); // ceil(√9)=3
    }

    #[test]
    fn scaled_partition_cap_never_drops_below_the_floor() {
        // A pathologically large library count still pre-sizes each cell a bit.
        assert_eq!(scaled_partition_cap(10_000), MIN_PARTITION_CAP);
    }

    #[test]
    fn library_index_is_single_bucket_with_no_libraries() {
        let idx = LibraryIndex::from_header(&header_with_read_groups(&[]), false);
        assert_eq!(idx.num_libs(), 1);
        assert_eq!(idx.name(0), UNKNOWN_LIBRARY);
    }

    #[test]
    fn library_index_is_single_bucket_with_one_library() {
        let idx =
            LibraryIndex::from_header(&header_with_read_groups(&[("A", Some("lib1"))]), false);
        assert_eq!(idx.num_libs(), 1);
        assert_eq!(idx.name(0), "lib1");
    }

    #[test]
    fn library_index_assigns_a_bucket_per_distinct_library() {
        let header = header_with_read_groups(&[("A", Some("lib1")), ("B", Some("lib2"))]);
        let idx = LibraryIndex::from_header(&header, false);
        // Bucket 0 is the unknown catch-all; lib1/lib2 take 1 and 2 (sorted).
        assert_eq!(idx.num_libs(), 3);
        assert_eq!(idx.name(0), UNKNOWN_LIBRARY);
        assert_ne!(idx.lookup(b"A"), idx.lookup(b"B"));
        assert_ne!(idx.lookup(b"A"), 0);
        assert_ne!(idx.lookup(b"B"), 0);
    }

    #[test]
    fn library_index_dedups_read_groups_by_library() {
        // Two read groups, one LB → one library → single-bucket mode.
        let header = header_with_read_groups(&[("A", Some("lib1")), ("B", Some("lib1"))]);
        let idx = LibraryIndex::from_header(&header, false);
        assert_eq!(idx.num_libs(), 1);
    }

    #[test]
    fn library_index_unknown_read_group_maps_to_bucket_zero() {
        let header = header_with_read_groups(&[("A", Some("lib1")), ("B", Some("lib2"))]);
        let idx = LibraryIndex::from_header(&header, false);
        // An RG that isn't in the header, or has no LB, falls to the unknown bucket.
        assert_eq!(idx.lookup(b"Z"), 0);
    }

    #[test]
    fn library_index_disabled_collapses_to_one_all_reads_bucket() {
        let header = header_with_read_groups(&[("A", Some("lib1")), ("B", Some("lib2"))]);
        let idx = LibraryIndex::from_header(&header, true);
        assert_eq!(idx.num_libs(), 1);
        assert_eq!(idx.name(0), ALL_READS);
    }
}
