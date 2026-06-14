//! Signature computation and the partitioned duplicate hash table.
//!
//! This module replicates samblaster's binning/padding strategy:
//!
//! * Every reference contig is padded by `2 * max_read_length` on both sides so
//!   that clipping never produces a negative or out-of-bounds genomic offset.
//! * A "super-contig" coordinate is `seq_off[seq_num] + pad_pos(pos)`.
//! * The super-contig is chopped into 2^27-position bins (a `bin_num`,
//!   `bin_pos` pair).
//!
//! The dedup state is a **2D-indexed array of small hash sets**, mirroring
//! the C++ `state->sigs` layout. The (bin_num, strand) tuple of each side of
//! a pair becomes an array index; the (bin_pos1, bin_pos2) pair becomes the
//! `u64` signature stored in that bucket. Compared to a single
//! `HashSet<(u32,u32,u32,u32)>` this halves the per-entry footprint (8 B
//! payload vs 16 B) and matches the C++ algorithm bit-for-bit.

use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hash, Hasher};

/// Fast hasher specialized for `u64` keys (our signature type).
///
/// We hash with a single `wrapping_mul` by the 64-bit golden-ratio
/// constant. For a u64 input that arrives via a single `write_u64`
/// call (which is what `u64`'s `Hash` impl does), this is
/// algebraically equivalent to the inner step of FxHash but skips the
/// `rotate_left + xor` (those collapse to no-ops on the all-zero
/// initial state). Crucially the multiply mixes input entropy into the
/// *high* bits of the output: hashbrown uses the top 7 bits as its H2
/// fingerprint byte, and an identity hash would collide every H2 slot
/// whenever the sig's high bits are zero. That's the common case here: the
/// sig packs `(bin_pos1 << 32) | bin_pos2` with each `bin_pos < 2^bin_shift`,
/// and `bin_shift` is well below 32 for any real genome (see
/// [`pick_bin_shift`]), so the top bits of both halves are zero. The
/// golden-ratio multiply spreads that entropy across the whole word.
#[derive(Default)]
pub struct U64Hasher(u64);

impl Hasher for U64Hasher {
    #[inline(always)]
    fn write_u64(&mut self, n: u64) {
        // 64-bit fractional part of the golden ratio — same constant
        // used by `splitmix64`, `fxhash`, and many other "one-shot"
        // hashers for u64 keys.
        self.0 = n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, _bytes: &[u8]) {
        // We only ever hash u64 keys; if this is called the caller is
        // wrong about the key type.
        unreachable!("U64Hasher only supports u64 keys");
    }
}

/// Fast hasher specialized for `u32` keys.
///
/// Up-casts to `u64` and multiplies by the **64-bit** golden-ratio
/// constant (the same as [`U64Hasher`]). The 64-bit multiply is
/// essential: hashbrown reads the H2 fingerprint from the *top* 7 bits
/// of the hash, so the multiplier must propagate input entropy into
/// the high half of the output. An early implementation used a 32-bit
/// multiplier and a `(u32) as u64` cast — the high 32 bits were
/// always zero, every entry shared the same H2 fingerprint, and
/// hashbrown's group probe degraded to a full per-cell linear scan.
/// That regressed wall time by ~70 % on the fragment table; the
/// 64-bit constant restores it. The internal state is `u64` so the
/// product is preserved end-to-end.
#[derive(Default)]
pub struct U32Hasher(u64);

impl Hasher for U32Hasher {
    #[inline(always)]
    fn write_u32(&mut self, n: u32) {
        self.0 = (n as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, _bytes: &[u8]) {
        unreachable!("U32Hasher only supports u32 keys");
    }
}

/// Maximum `bin_shift` we can use: `bin_pos1` and `bin_pos2` are packed
/// 32-and-32 into the u64 sig, so `bin_pos < 2^32`.
const MAX_BIN_SHIFT: u32 = 32;

/// Pick `bin_shift` dynamically given the genome super-contig length.
///
/// Returns the largest shift such that:
///   1. The sig fits in u64 (`bin_shift ≤ 32`).
///   2. The resulting `bin_count = ceil(total_len / 2^bin_shift)` is at
///      least `min_bin_count` (the per-side floor).
///
/// `min_bin_count` controls the partition cell count
/// (`cells ≈ (bin_count+1)² × 4`). Smaller floor → fewer, larger cells
/// (less per-cell overhead but bigger resize peaks). Larger floor → more,
/// smaller cells (peaks stay tiny, slightly more bookkeeping at low N).
pub fn pick_bin_shift(total_len: u64, min_bin_count: u32) -> u32 {
    let min = min_bin_count.max(1) as u64;
    for shift in (0..=MAX_BIN_SHIFT).rev() {
        let bin_count = if shift >= 64 { 0 } else { (total_len + (1u64 << shift) - 1) >> shift };
        if bin_count >= min {
            return shift;
        }
    }
    0
}

/// Per-genome bin lookup table — built once from the SAM header.
pub struct BinIndex {
    /// `seq_offs[seq_num]` — start of contig `seq_num` in the super-contig.
    seq_offs: Vec<u64>,
    /// Total padded length of the synthetic super-contig (for `bin_count`).
    total_len: u64,
    /// Padding added to each side of every contig (`--max-read-length`).
    /// Stored so [`Self::bin_for`] can pad a raw alignment position correctly.
    max_read_length: i32,
    /// Bits per `bin_pos`, chosen at startup so the sig fits in u64 while
    /// keeping `bin_count` ≥ the requested floor.
    bin_shift: u32,
    /// Cached `(1 << bin_shift) - 1`.
    bin_mask: u64,
}

impl BinIndex {
    /// Build from a flat slice of contig lengths in BAM `tid` order
    /// (`ref_lengths[0]` is `tid=0`). Sequence number 0 is reserved for `"*"`;
    /// contigs are mapped to seq numbers 1..=N in the order they appear.
    /// `min_bin_count` is the floor on bins per side (controls partition
    /// density — see [`pick_bin_shift`]).
    pub fn from_ref_lengths(ref_lengths: &[i32], max_read_length: i32, min_bin_count: u32) -> Self {
        let mut seq_offs = Vec::with_capacity(ref_lengths.len() + 1);
        let mut total_len: u64 = 0;
        // "*" entry — pad_length(0) = 2 * max_read_length, then "+1" terminator.
        seq_offs.push(0);
        total_len += pad_length_for(0, max_read_length) as u64 + 1;
        for &len in ref_lengths {
            let seq_off = total_len;
            seq_offs.push(seq_off);
            total_len += pad_length_for(len, max_read_length) as u64 + 1;
        }
        let bin_shift = pick_bin_shift(total_len, min_bin_count);
        let bin_mask = if bin_shift >= 64 { u64::MAX } else { (1u64 << bin_shift) - 1 };
        Self { seq_offs, total_len, max_read_length, bin_shift, bin_mask }
    }

    /// Bits used by `bin_pos` in the packed sig.
    pub fn bin_shift(&self) -> u32 {
        self.bin_shift
    }

    /// `bin_count` per side. With the dynamic `bin_shift`, this is
    /// `ceil(total_len / 2^bin_shift)`.
    pub fn bin_count(&self) -> u32 {
        if self.bin_shift >= 64 {
            1
        } else {
            ((self.total_len + self.bin_mask) >> self.bin_shift) as u32
        }
    }

    /// Number of sequence slots: one per `@SQ` contig plus the reserved
    /// slot 0 for the unmapped (`*`) reference. A record's `seq_num`
    /// (`tid + 1`) must be `< num_seqs()`.
    pub fn num_seqs(&self) -> usize {
        self.seq_offs.len()
    }

    /// Compute `(bin_num, bin_pos, was_clamped)` for an alignment's 5'
    /// aligned position.
    ///
    /// dupblaster lays every contig end-to-end on one synthetic axis with
    /// `max_read_length` of padding on each side of each contig (see
    /// [`Self::from_ref_lengths`]). A read whose 5' position clips more than
    /// `max_read_length` past a contig edge would otherwise land in a
    /// *neighbouring* contig's bins (or underflow the axis), producing
    /// spurious cross-contig duplicate collisions. We clamp the padded
    /// coordinate into this contig's own block so that can't happen.
    ///
    /// `was_clamped` is `true` only when the clamp actually moved the
    /// coordinate — i.e. the clip exceeded the padding. The caller counts
    /// those (per template) so the run can warn that `--max-read-length` is
    /// too small for the data and that edge duplicate-marking may be
    /// imprecise. Reads clipped within the padding are represented exactly
    /// and never flagged.
    pub fn bin_for(&self, seq_num: usize, pos: i32) -> (u32, u32, bool) {
        let padded = self.seq_offs[seq_num] as i64 + pad_pos(pos, self.max_read_length) as i64;
        // This contig occupies the half-open block `[seq_offs[seq_num],
        // next_off)`; clamp into it (inclusive of the final slot) so the
        // coordinate can never reach the next contig's block or go negative.
        let lo = self.seq_offs[seq_num] as i64;
        let next_off = self.seq_offs.get(seq_num + 1).copied().unwrap_or(self.total_len);
        let hi = next_off as i64 - 1;
        let clamped = padded.clamp(lo, hi);
        let was_clamped = clamped != padded;
        let combined = clamped as u64;
        let bin_num = (combined >> self.bin_shift) as u32;
        let bin_pos = (combined & self.bin_mask) as u32;
        (bin_num, bin_pos, was_clamped)
    }
}

/// `pad_length(n) = n + 2 * max_read_length`.
pub fn pad_length_for(n: i32, max_read_length: i32) -> i32 {
    n + 2 * max_read_length
}

/// `pad_pos(p) = p + max_read_length`.
pub fn pad_pos(p: i32, max_read_length: i32) -> i32 {
    p + max_read_length
}

/// Strand-aware 5' aligned reference position from CIGAR analysis.
pub fn five_prime_aligned_pos(
    rapos: i32,
    sclip: i32,
    eclip: i32,
    ra_len: i32,
    is_reverse: bool,
) -> i32 {
    if is_reverse { rapos + ra_len + eclip - 1 } else { rapos - sclip }
}

/// samblaster-legacy override for reverse-strand orphans: replace the
/// strand-aware 5' coordinate with the leftmost-aligned coordinate so
/// fwd/rev orphans at the same locus collide. Behavior introduced in
/// samblaster v0.1.23 (March 2020). Used only when the single-end
/// strategy is [`SingleEndStrategy::SamblasterLegacy`].
pub fn orphan_pos_override(rapos: i32, sclip: i32) -> i32 {
    rapos - sclip
}

/// How dupblaster keys single-end / orphan reads.
///
/// See the README's "single-end strategies" section for the algorithmic
/// background and the reason this is configurable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SingleEndStrategy {
    /// Strand-aware 5'-aligned coordinate (the default). Two single-end
    /// or orphan reads are duplicates only if their reads' strand-aware
    /// 5' positions match. Picard MarkDuplicates uses an equivalent key
    /// in its `fragSort` collection.
    StrandAware,
    /// Strand-aware keying plus a Picard-style cross-check: each end of
    /// a fully-mapped pair is also registered in a fragment-level table
    /// so that subsequent orphans / single-end reads at those positions
    /// are marked as duplicates of the pair (approximating Picard's
    /// "fragments don't beat pairs" rule). Approximate in a single
    /// streaming pass — an orphan that arrives *before* its
    /// corresponding pair will pass through as non-dup. ~2x memory of
    /// `StrandAware` because of the extra fragment table.
    PicardApprox,
    /// Exact Picard fragment semantics via a deferred two-pass scheme:
    /// fully-mapped/unmapped pairs stream straight to the output while
    /// every mapped-orphan / single-end "fragment" block is buffered to a
    /// temporary uncompressed BAM. After the pair pass completes, the pair
    /// table is *consumed* into a fragment table (via
    /// [`PairDupTable::drain_into_fragment_table`]) holding the 5' position
    /// of every paired read end; the buffered fragments are then re-read
    /// and checked against it. Because the fragment table is complete
    /// before any fragment is checked, the "fragments never beat pairs"
    /// rule holds exactly and independent of input order — unlike
    /// `PicardApprox`. Costs a temporary on-disk copy of the fragment
    /// blocks (a small fraction of paired data) and emits those fragments
    /// at the *end* of the output stream rather than in input order.
    /// Intended for paired data only.
    PicardExact,
    /// samblaster v0.1.23+ behavior: leftmost-aligned reference
    /// coordinate with strand bit dropped, so forward and reverse
    /// orphans whose alignments share a leftmost coord collide.
    /// Discouraged for short-read PE data — see the README's
    /// "single-end strategies" section for the discussion.
    SamblasterLegacy,
}

/// Whether (and how) to key duplicates in a methylation-aware way.
///
/// Standard WGS keying canonicalizes each pair to a coordinate-ordered
/// signature (leftmost end first). That merges the two original strands
/// (OT/CTOT and OB/CTOB) of a bisulfite- or enzymatic-conversion fragment at
/// the same locus — wrong for methylation data, where those strands carry
/// independent methylation and must be counted separately. Methylation mode
/// keys the pair in *template order* instead, so opposite-strand fragments
/// separate while same-strand PCR copies still collapse.
///
/// Only the directional case is implemented (see the variant). Non-directional
/// (PBAT) libraries are deliberately out of scope: their first-of-pair-to-
/// strand relationship is not fixed, so correct keying needs per-read strand
/// tags *plus* a canonicalize-within-strand key — a separate, larger piece of
/// work. The enum exists (rather than a bool) so that variant can slot in later
/// without changing the option's type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethylationMode {
    /// Directional preps (WGBS / EM-seq / TAPS), where sequencing adapters are
    /// ligated to the intact double-stranded fragment *before* conversion, so
    /// first-of-pair is locked to a consistent end of the original strand. The
    /// pair is keyed in template (first-of-pair → second-of-pair) order rather
    /// than coordinate-canonical order. The per-end strand bits already in the
    /// partition cell index then encode orientation for free, so an OT key
    /// differs from the OB key at the same locus while PCR copies of one strand
    /// still collide.
    Directional,
}

// ── DupTable: partitioned dedup hash sets ──────────────────────────────────

/// Partitioned hash table holding dedup signatures.
///
/// Generic over the signature type `S` and its specialized hasher `H`.
/// Two concrete instantiations are used in dupblaster, via the
/// [`PairDupTable`] and [`FragmentDupTable`] type aliases:
///
/// * **Pair table** ([`PairDupTable`]): `S = u64` because pair signatures
///   pack `(bin_pos1 << 32) | bin_pos2`. Built with [`Self::new_pair`]
///   which allocates a 2D `stride²` cell array indexed by
///   `off = s1 * stride + s2`; each cell stores u64 sigs for that
///   `(bin_num, strand) × (bin_num, strand)` cohort. The orphan path
///   on this table uses the `s1 = 0` row.
/// * **Fragment table** ([`FragmentDupTable`]): `S = u32` because each
///   entry is a single `bin_pos` (≤ 27 bits). Built with
///   [`Self::new_single_end`] which allocates a 1D `stride`-length cell
///   array indexed by `off = s = bin_num*2 + strand`. Used only by
///   [`SingleEndStrategy::PicardApprox`].
///
/// ## Why this partitioning level?
///
/// The pair-sig u64 has ~10 free bits (5 above each 27-bit `bin_pos`),
/// so it's tempting to pack `bin_num` into them and partition by strand
/// only — 5 cells instead of ~2200. Measured: same wall time, **+60 MB
/// memory**. The reason is hashbrown's power-of-2 growth: collapsing
/// entries into a single hot cell rounds the table up to the next
/// power-of-2 (~8M slots × 9 B = 72 MB for the hot FR cell), while
/// the 2200-cell scheme has many small tables each rounding up to a
/// much smaller power-of-2 (total ~60 MB). Many small powers-of-2 sum
/// to less than one big power-of-2 — so this layout is the memory
/// sweet spot.
pub struct DupTable<S, H>
where
    S: Copy + Eq + Hash,
    H: Hasher + Default,
{
    /// Cell-row width = `(bin_count + 1) * 2`. Used by [`Self::check_dm`]
    /// to compute the 2D offset `s1 * stride + s2`; not consulted by
    /// fragment-table operations (which only use the `s1 = 0` row).
    stride: u32,
    /// The flat array of partition cells, each a `HashSet` of signature values.
    /// Indexed by a 1D offset derived from `(s1 * stride + s2)` for pair tables
    /// or just `s` for fragment tables.
    sets: Vec<HashSet<S, BuildHasherDefault<H>>>,
}

/// Pair-key dedup table (`u64` sigs, 2D-partitioned). Used for the
/// fully-mapped pair-signature index and — under `StrandAware` /
/// `SamblasterLegacy` — also for the orphan-on-pair-table single-end
/// keying.
pub type PairDupTable = DupTable<u64, U64Hasher>;

/// Fragment-key dedup table (`u32` sigs, 1D-partitioned). Stores each end of
/// every fully-mapped PE pair plus every orphan, enabling the orphan ↔ pair-end
/// cross-check. Allocated under [`SingleEndStrategy::PicardApprox`] (populated
/// concurrently with the pair pass) and [`SingleEndStrategy::PicardExact`]
/// (populated after the pair pass by draining the pair table; see
/// [`PairDupTable::drain_into_fragment_table`]). Not used by `StrandAware` /
/// `SamblasterLegacy`, which key orphans on the pair table's `s1 = 0` row.
pub type FragmentDupTable = DupTable<u32, U32Hasher>;

impl<S, H> DupTable<S, H>
where
    S: Copy + Eq + Hash,
    H: Hasher + Default,
{
    /// Build a 2D-partitioned pair table (`stride²` cells). `bin_count`
    /// is `total_padded_len >> BIN_SHIFT`. Because positions can land
    /// exactly on a bin boundary (so `bin_num` reaches `bin_count`,
    /// inclusive) and `s = bin_num*2 + strand`, we size the stride for
    /// `s ∈ [0, (bin_count+1)*2)` and the array as `stride²`.
    ///
    /// `cell_cap` is the initial bucket capacity pre-allocated for each
    /// partition cell (see `with_n_cells`).
    pub fn new_pair(bin_count: u32, cell_cap: usize) -> Self {
        let stride = (bin_count + 1) * 2;
        let n = (stride as usize).saturating_mul(stride as usize);
        Self::with_n_cells(stride, n, cell_cap)
    }

    /// Build a 1D-partitioned single-end/fragment table (`stride` cells).
    /// Used when only the `s1 = 0` row of a 2D table would be populated,
    /// e.g. the fragment table under `PicardApprox`. Saves the
    /// `stride²` − `stride` empty cells the 2D layout would otherwise
    /// allocate.
    pub fn new_single_end(bin_count: u32, cell_cap: usize) -> Self {
        let stride = (bin_count + 1) * 2;
        let n = stride as usize;
        Self::with_n_cells(stride, n, cell_cap)
    }

    /// `cell_cap` is the initial bucket capacity for each partition cell.
    /// Pre-sizing absorbs the early 2× rehash steps in the hot partitions
    /// where most entries land, at a small upfront memory cost. Profile
    /// showed `hashbrown::reserve_rehash` at ~3% of total wall before
    /// pre-sizing; after, it disappears from the flamegraph. The caller
    /// scales it down by library count (see `scaled_partition_cap` in
    /// `dedup.rs`) so the empty-table baseline doesn't grow linearly when
    /// the dedup state is partitioned per library.
    fn with_n_cells(stride: u32, n: usize, cell_cap: usize) -> Self {
        let mut sets = Vec::with_capacity(n);
        sets.resize_with(n, || HashSet::with_capacity_and_hasher(cell_cap, Default::default()));
        Self { stride, sets }
    }
}

// ── pair-table-specific operations (u64 sigs) ──────────────────────────────

impl PairDupTable {
    /// Insert and check the doubly-mapped (DM) signature. Returns `true`
    /// if this signature was already present (i.e. the pair is a
    /// duplicate). Caller must ensure the table was built with
    /// [`Self::new_pair`].
    #[inline]
    pub fn check_dm(
        &mut self,
        bin_num1: u32,
        bin_pos1: u32,
        rev1: bool,
        bin_num2: u32,
        bin_pos2: u32,
        rev2: bool,
    ) -> bool {
        let s1 = bin_num1 * 2 + (rev1 as u32);
        let s2 = bin_num2 * 2 + (rev2 as u32);
        let off = (s1 * self.stride + s2) as usize;
        let sig = ((bin_pos1 as u64) << 32) | (bin_pos2 as u64);
        !self.sets[off].insert(sig)
    }

    /// Insert and check the orphan / single-end signature on the pair
    /// table (used by `StrandAware` and `SamblasterLegacy`, which don't
    /// allocate a separate fragment table). The first side is the
    /// (placeholder) unmapped read — its `s1` and `bin_pos` contribute
    /// zero.
    ///
    /// Only `StrandAware` and `SamblasterLegacy` reach this method —
    /// `PicardApprox` / `PicardExact` key orphans on their own
    /// [`FragmentDupTable`] instead, never here.
    ///
    /// * [`SingleEndStrategy::StrandAware`]: `s2 = bin_num*2 + (rev as u32)`
    ///   — strand-aware, so a forward orphan and a reverse orphan at the
    ///   same 5'-aligned position land in different cells and are *not*
    ///   duplicates of each other.
    /// * [`SingleEndStrategy::SamblasterLegacy`]: `s2 = bin_num*2` —
    ///   strand bit is dropped, so any two orphans whose alignments
    ///   share a leftmost-aligned coordinate (the upstream caller is
    ///   expected to pass that coordinate via `bin_pos2`) collide,
    ///   regardless of strand.
    #[inline]
    pub fn check_orphan(
        &mut self,
        bin_num2: u32,
        bin_pos2: u32,
        rev2: bool,
        strategy: SingleEndStrategy,
    ) -> bool {
        let s2 = match strategy {
            SingleEndStrategy::StrandAware => bin_num2 * 2 + (rev2 as u32),
            SingleEndStrategy::SamblasterLegacy => bin_num2 * 2,
            // The picard strategies route orphans through a FragmentDupTable
            // (see `check_fragment_signature` in `dedup.rs`), never here;
            // panic loudly if a future refactor wires them in by mistake.
            SingleEndStrategy::PicardApprox | SingleEndStrategy::PicardExact => {
                unreachable!("check_orphan is only reached from StrandAware/SamblasterLegacy")
            }
        };
        let off = s2 as usize;
        let sig = bin_pos2 as u64;
        !self.sets[off].insert(sig)
    }

    /// Consume every stored pair signature, inserting *both* 5'-aligned ends
    /// of each pair into a freshly-built [`FragmentDupTable`], freeing each
    /// pair cell as it is drained so peak memory stays near a single table
    /// rather than holding the full pair and fragment tables at once.
    ///
    /// Used by [`SingleEndStrategy::PicardExact`]: after the streaming pair
    /// pass finishes, the returned table holds the 5' position of every
    /// paired read end, so a deferred orphan / single-end read sharing one of
    /// those 5' coordinates is marked as a duplicate of the pair ("fragments
    /// never beat pairs"), exactly and independent of input order.
    ///
    /// Reconstruction is lossless because the pair table stores a pair's full
    /// signature: the cell index encodes each end's `(bin_num, strand)` (via
    /// `s = bin_num*2 + strand`) and the `u64` sig packs the two `bin_pos`
    /// halves. Pair canonicalization (the `need_swap`-style ordering
    /// applied before insertion) is irrelevant here because each end is
    /// re-inserted independently.
    ///
    /// `cell_cap` is forwarded to the new fragment table's per-cell
    /// pre-sizing (see `with_n_cells`).
    pub fn drain_into_fragment_table(&mut self, cell_cap: usize) -> FragmentDupTable {
        let stride = self.stride as usize;
        // `new_single_end(bin_count)` rebuilds `stride = (bin_count + 1) * 2`,
        // so this reproduces the pair table's stride exactly.
        let bin_count = self.stride / 2 - 1;
        let mut frag = FragmentDupTable::new_single_end(bin_count, cell_cap);
        for off in 0..self.sets.len() {
            // Recover each end's (bin_num, strand) from the 2D cell index.
            let s1 = (off / stride) as u32;
            let s2 = (off % stride) as u32;
            let (bin_num1, rev1) = (s1 >> 1, s1 & 1 == 1);
            let (bin_num2, rev2) = (s2 >> 1, s2 & 1 == 1);
            // Move the cell out so its backing allocation is dropped as soon
            // as we've drained it — this is what keeps the peak near one table.
            let cell = std::mem::take(&mut self.sets[off]);
            for sig in cell {
                let bin_pos1 = (sig >> 32) as u32;
                let bin_pos2 = (sig & 0xFFFF_FFFF) as u32;
                frag.check_or_insert(bin_num1, bin_pos1, rev1);
                frag.check_or_insert(bin_num2, bin_pos2, rev2);
            }
        }
        frag
    }
}

// ── fragment-table-specific operations (u32 sigs) ──────────────────────────

impl FragmentDupTable {
    /// Insert `(bin_num, bin_pos, rev)` into the fragment table; return
    /// `true` if the entry was already present. Strand-aware keying only
    /// (the table is allocated only under `PicardApprox`, which doesn't
    /// use the legacy strand-drop variant). Caller must ensure the
    /// table was built with [`DupTable::new_single_end`].
    ///
    /// This is the only operation the table supports — there is no
    /// "register-without-checking" variant because the underlying
    /// `HashSet::insert` already conveys both pieces of information in
    /// one call; callers that don't care about the return value
    /// (paired-end-end inserts in `RecordProcessor`) simply ignore it.
    #[inline]
    pub fn check_or_insert(&mut self, bin_num: u32, bin_pos: u32, rev: bool) -> bool {
        let s = bin_num * 2 + (rev as u32);
        let off = s as usize;
        !self.sets[off].insert(bin_pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin_for_does_not_clamp_positions_within_the_padding() {
        // Positions in [-max_read_length, len+max_read_length] are inside the
        // contig's padded block and must be represented exactly (no clamp).
        let bins = BinIndex::from_ref_lengths(&[1_000_000], 1000, 32);
        assert!(!bins.bin_for(1, 0).2);
        assert!(!bins.bin_for(1, 999_999).2);
        assert!(!bins.bin_for(1, -1000).2); // exactly at the start padding edge
    }

    #[test]
    fn bin_for_clamps_positions_beyond_the_padding() {
        // A 5' position pushed more than max_read_length before the contig
        // start (e.g. a long read with thousands of bases of leading clip)
        // would bleed into a neighbouring contig's bins; it is clamped into
        // this contig's block and flagged so the caller can warn.
        let bins = BinIndex::from_ref_lengths(&[1_000_000], 1000, 32);
        let (_, _, clamped_start) = bins.bin_for(1, -10_000);
        assert!(clamped_start, "a clip far past the start should clamp");
        // Likewise far past the contig end.
        let (_, _, clamped_end) = bins.bin_for(1, 1_010_000);
        assert!(clamped_end, "a clip far past the end should clamp");
    }

    #[test]
    fn num_seqs_counts_contigs_plus_unmapped_slot() {
        let bins = BinIndex::from_ref_lengths(&[100, 200, 300], 1000, 32);
        assert_eq!(bins.num_seqs(), 4); // 3 contigs + the reserved `*` slot
    }

    #[test]
    fn dm_repeat_detected_as_duplicate() {
        let mut t = PairDupTable::new_pair(8, 64);
        assert!(!t.check_dm(1, 100, false, 2, 200, true));
        assert!(t.check_dm(1, 100, false, 2, 200, true));
    }

    #[test]
    fn dm_strand_differs_no_collision() {
        let mut t = PairDupTable::new_pair(8, 64);
        assert!(!t.check_dm(1, 100, false, 2, 200, false));
        // Different strand on the second alignment → different cell.
        assert!(!t.check_dm(1, 100, false, 2, 200, true));
    }

    #[test]
    fn strand_aware_orphan_repeat_detected() {
        let mut t = PairDupTable::new_pair(8, 64);
        // Two reverse-strand orphans with the same strand-aware 5' coord
        // (the caller already passed `rev2=true` and the 5'-aligned bin_pos).
        assert!(!t.check_orphan(3, 555, true, SingleEndStrategy::StrandAware));
        assert!(t.check_orphan(3, 555, true, SingleEndStrategy::StrandAware));
    }

    #[test]
    fn strand_aware_orphan_fwd_rev_do_not_collide() {
        let mut t = PairDupTable::new_pair(8, 64);
        // Forward orphan at bin (3, 555) and reverse orphan at the same
        // (3, 555) — under StrandAware they go to different cells and
        // do NOT collide.
        assert!(!t.check_orphan(3, 555, false, SingleEndStrategy::StrandAware));
        assert!(!t.check_orphan(3, 555, true, SingleEndStrategy::StrandAware));
    }

    #[test]
    fn samblaster_legacy_orphan_fwd_rev_collide() {
        let mut t = PairDupTable::new_pair(8, 64);
        // Same locus — strand is intentionally not part of the key for
        // SamblasterLegacy (the caller is expected to pass the leftmost-
        // aligned coord via bin_pos2, regardless of strand).
        assert!(!t.check_orphan(3, 555, false, SingleEndStrategy::SamblasterLegacy));
        assert!(t.check_orphan(3, 555, true, SingleEndStrategy::SamblasterLegacy));
    }

    #[test]
    fn fragment_dup_table_repeat_detected() {
        let mut t = FragmentDupTable::new_single_end(8, 64);
        assert!(!t.check_or_insert(3, 555, true));
        assert!(t.check_or_insert(3, 555, true));
    }

    #[test]
    fn fragment_dup_table_strand_aware() {
        let mut t = FragmentDupTable::new_single_end(8, 64);
        // Fwd and rev at the same (bin_num, bin_pos) go to different
        // cells under strand-aware keying.
        assert!(!t.check_or_insert(3, 555, false));
        assert!(!t.check_or_insert(3, 555, true));
    }

    #[test]
    fn drain_recovers_both_pair_ends_into_fragment_table() {
        let mut pairs = PairDupTable::new_pair(8, 64);
        // Insert one pair: end A = (bin 1, pos 100, fwd), end B = (bin 2, pos 200, rev).
        assert!(!pairs.check_dm(1, 100, false, 2, 200, true));
        let mut frag = pairs.drain_into_fragment_table(64);
        // Both ends of the pair are now present in the fragment table: a
        // re-insert at either coordinate reports a collision.
        assert!(frag.check_or_insert(1, 100, false), "pair end A should be present");
        assert!(frag.check_or_insert(2, 200, true), "pair end B should be present");
        // A coordinate that was never part of any pair is absent.
        assert!(!frag.check_or_insert(3, 300, false), "unrelated coord should be absent");
    }

    #[test]
    fn drain_is_strand_aware() {
        let mut pairs = PairDupTable::new_pair(8, 64);
        // Pair end at (bin 4, pos 555, forward).
        assert!(!pairs.check_dm(4, 555, false, 5, 600, false));
        let mut frag = pairs.drain_into_fragment_table(64);
        // The forward end collides; the reverse-strand coordinate at the same
        // bin/pos is a distinct key and does not.
        assert!(!frag.check_or_insert(4, 555, true), "reverse coord must not collide with fwd end");
        assert!(frag.check_or_insert(4, 555, false), "forward end should be present");
    }

    #[test]
    fn fragment_dup_table_register_and_query() {
        let mut t = FragmentDupTable::new_single_end(8, 64);
        // First insert returns false ("not previously seen").
        assert!(!t.check_or_insert(3, 555, false));
        // A subsequent query at the same coord returns true.
        assert!(t.check_or_insert(3, 555, false));
    }
}
