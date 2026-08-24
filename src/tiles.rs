//! Interning imaging locations, and the tile statistics derived from them.
//!
//! Every template's `(library, sequencing unit, tile)` triple is interned to a
//! small serial ID. Duplicate groups are later reconstructed from those IDs to
//! split duplicates into a sequencing component (copies of one molecule imaged
//! in one place) and a library component (independent molecules), so this
//! dictionary is the bridge between a read name and the decomposition.
//!
//! **IDs are serial and reversible, not hashed.** That buys four things a hash
//! would not: it cannot collide, so a tile's share of reads is measured on real
//! tiles rather than on hash buckets; it is deterministic, with no seeded-hasher
//! caveat on the reported metrics; it reverses, so the per-sequencing-unit report
//! can print real flowcell and lane names instead of opaque numbers; and it
//! self-diagnoses a misconfigured extractor, because pointing the tile field at
//! an x coordinate makes the cardinality explode instead of silently producing a
//! plausible wrong answer.
//!
//! IDs are assigned in first-seen order, so the *order* of assignment depends on
//! the input's order — but the equality relations between triples do not, and the
//! metrics depend only on those. Integer counts are therefore exactly
//! order-invariant. The floating-point figures (`q`, and the chance-corrected
//! count) sum in ID order and so can differ in their last bits between orderings
//! of the same file; they are stable to far more precision than the six decimals
//! reported.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use tempfile::TempDir;

use crate::readname::{ImagingLocation, ReadNameFormat};
use crate::sig::{PairSlot, stride_for};

/// Default number of files the spill is split across.
///
/// One file per partition cell is not possible — a run reports on the order of
/// 9,216 cells, against soft descriptor limits of 256 (macOS) and 1,024 (Linux)
/// — yet every record of a duplicate group must land in one file for the group to
/// be reassembled. Hashing the group key into a small K gives both, for 64
/// descriptors. See [`TileSpiller::bucket_of`].
pub(crate) const DEFAULT_SPILL_BUCKETS: u32 = 64;

/// 64-bit fractional part of the golden ratio, the same constant
/// [`crate::sig::U64Hasher`] uses to spread low-entropy signature bits.
const GOLDEN_RATIO_64: u64 = 0x9E37_79B9_7F4A_7C15;

/// Descriptors held back from the bucket count for everything else the run
/// needs: input, output, the picard-exact orphan temp, the metrics and plot
/// files, plus slack for whatever the runtime holds open.
const RESERVED_DESCRIPTORS: u64 = 32;

/// Write buffer per bucket. At the default bucket count this is 4 MB of buffers
/// in total, which keeps each bucket's writes large enough to stay sequential.
const BUCKET_BUFFER_BYTES: usize = 64 * 1024;

/// Read buffer for a bucket in the post-pass, a whole number of records so no
/// record ever straddles two reads.
const BUCKET_READ_RECORDS: usize = 64 * 1024;

/// Bisection steps used to invert `E[n_tiles | m]`. Well past `f64` resolution
/// for the magnitudes involved; the loop breaks early once the bracket is tight,
/// so this is only a backstop.
const BISECTION_STEPS: usize = 60;

/// Ceiling on the molecule count the inverse searches to, so a group occupying
/// nearly every tile cannot send the bracket search unbounded.
const INVERSE_SEARCH_CEILING: f64 = 1e12;

/// Distinct triples above which the extractor is very likely misconfigured.
///
/// Real geometry does not come close: a NovaSeq S4 run has 704 tiles per lane,
/// and three flowcells of merged data measured 6,334 distinct triples. A million
/// is also where the dictionary itself starts to cost real memory (~50 MB).
const CARDINALITY_WARN: usize = 1_000_000;

/// Distinct triples treated as pathological rather than merely suspicious.
///
/// At this point the dictionary alone would be roughly 800 MB, which no real
/// flowcell geometry can justify — it means the extractor is pointed at
/// something per-read, such as an x/y coordinate.
const CARDINALITY_LIMIT: usize = 16_777_216;

/// One interned `(library, sequencing unit, tile)` triple.
#[derive(Clone, Debug)]
pub(crate) struct TileEntry {
    /// Library bucket, as assigned by [`crate::dedup::LibraryIndex`]. Part of the
    /// triple because duplicates are only ever called within a library, so tiles
    /// of different libraries are never compared.
    pub library: u32,
    /// Sequencing-unit token, verbatim from the read name, so the per-unit report
    /// prints real flowcell and lane names rather than opaque numbers.
    ///
    /// The tile token is deliberately *not* kept beside it: only tile *identity*
    /// is ever needed, and the packed dictionary key already carries it. No output
    /// names an individual tile.
    pub unit: Box<[u8]>,
    /// Templates observed at this triple, the numerator of the tile's share.
    pub templates: u64,
}

/// Tile aggregates for one library, all derived from the same single pass.
#[derive(Clone, Debug, Default)]
pub(crate) struct LibraryTiles {
    /// Templates observed for this library across all of its tiles.
    pub templates: u64,
    /// Distinct tiles seen for this library.
    pub tiles: usize,
    /// `q = Σ w_t²`: the chance two unrelated templates share a tile. The
    /// validity indicator for the whole split — a library on one tile has `q = 1`
    /// and carries no information, since every duplicate is on "the" tile whether
    /// it was clustered or not.
    pub collision_rate: f64,
    /// The shares `w_t` themselves, for the chance model.
    pub shares: Vec<f64>,
}

/// Interns imaging locations and accumulates per-tile template counts.
pub(crate) struct TileDictionary {
    /// How to pull the unit and tile out of a read name; chosen by the user.
    format: ReadNameFormat,
    /// Packed triple key → serial ID. See [`Self::pack_key`] for the encoding.
    ///
    /// FxHash rather than std's default SipHash: this map is probed once per
    /// template on input whose neighbours don't share a tile, and the keys are
    /// ~22 bytes, short enough that SipHash's per-word mixing dominates the
    /// probe. Fixed-key hashing forgoes SipHash's HashDoS resistance, which
    /// matches the threat model the dedup tables' fixed multiply hashers
    /// ([`crate::sig::U64Hasher`]) already accept — collisions cost probes,
    /// never correctness, and the cardinality bail bounds the dictionary.
    ids: HashMap<Box<[u8]>, u32, rustc_hash::FxBuildHasher>,
    /// Interned triples, indexed by ID.
    entries: Vec<TileEntry>,
    /// Reused buffer for the lookup key, so probing allocates nothing.
    key: Vec<u8>,
    /// The previous template's triple and ID, reused in place.
    memo: Memo,
    /// Whether [`CARDINALITY_WARN`] has already been reported.
    warned: bool,
}

impl TileDictionary {
    /// Create an empty dictionary that extracts read names using `format`.
    pub(crate) fn new(format: ReadNameFormat) -> Self {
        Self {
            format,
            ids: HashMap::default(),
            entries: Vec::new(),
            key: Vec::new(),
            memo: Memo::default(),
            warned: false,
        }
    }

    /// Record one template of `library` whose read name is `name`, returning the
    /// ID of its triple.
    ///
    /// Call exactly once per template — not once per read — since both mates
    /// share a QNAME and therefore a tile. Errors if `name` does not match the
    /// chosen format: the user named that format, so a name that does not fit it
    /// means the wrong one was chosen or the data is not what it claims to be.
    pub(crate) fn observe(&mut self, library: u32, name: &[u8]) -> Result<u32> {
        let location = self.format.extract(name).ok_or_else(|| self.format.parse_error(name))?;
        let id = self.intern(library, location)?;
        self.entries[id as usize].templates += 1;
        Ok(id)
    }

    /// Resolve `location` to its ID, assigning a fresh one if it is new.
    fn intern(&mut self, library: u32, location: ImagingLocation<'_>) -> Result<u32> {
        self.pack_key(library, location)?;
        // Consecutive templates of a name-sorted file almost always share a
        // tile, so this one-entry memo of the previous packed key skips the
        // hash probe on the common case — one length-checked memcmp against
        // the key just packed.
        if self.memo.key == self.key {
            return Ok(self.memo.id);
        }

        let id = match self.ids.get(self.key.as_slice()) {
            Some(&id) => id,
            None => self.insert(library, location)?,
        };
        // Swap rather than copy: the old memo buffer becomes the next
        // pack_key scratch, so storing the memo moves no bytes at all.
        std::mem::swap(&mut self.memo.key, &mut self.key);
        self.memo.id = id;
        Ok(id)
    }

    /// Build the packed lookup key for a triple into [`Self::key`].
    ///
    /// The unit is length-prefixed rather than separated by a delimiter byte:
    /// the tokens are opaque, so a custom regex could capture any byte at all,
    /// and any delimiter we picked could appear inside a token and make
    /// `("AB", "C")` collide with `("A", "BC")`.
    fn pack_key(&mut self, library: u32, location: ImagingLocation<'_>) -> Result<()> {
        let Ok(unit_len) = u16::try_from(location.unit.len()) else {
            bail!(
                "sequencing-unit token is {} bytes, which cannot be a read-name field \
                 (SAM limits a QNAME to 254 bytes) — check the --read-name-format pattern",
                location.unit.len()
            );
        };
        self.key.clear();
        self.key.extend_from_slice(&library.to_le_bytes());
        self.key.extend_from_slice(&unit_len.to_le_bytes());
        self.key.extend_from_slice(location.unit);
        self.key.extend_from_slice(location.tile);
        Ok(())
    }

    /// Assign a fresh ID to a triple not yet in the dictionary.
    fn insert(&mut self, library: u32, location: ImagingLocation<'_>) -> Result<u32> {
        if self.entries.len() >= CARDINALITY_LIMIT {
            bail!(
                "more than {CARDINALITY_LIMIT} distinct (library, sequencing unit, tile) \
                 triples: the --read-name-format is almost certainly extracting a per-read \
                 field such as an x/y coordinate rather than a tile"
            );
        }
        let id = self.entries.len() as u32;
        self.ids.insert(self.key.clone().into_boxed_slice(), id);
        self.entries.push(TileEntry { library, unit: location.unit.into(), templates: 0 });
        if !self.warned && self.entries.len() >= CARDINALITY_WARN {
            self.warned = true;
            log::warn!(
                "{} distinct (library, sequencing unit, tile) triples seen, far more than any \
                 real flowcell geometry — check that --read-name-format names the tile field \
                 and not an x/y coordinate.",
                self.entries.len()
            );
        }
        Ok(id)
    }

    /// Every interned triple, indexed by ID.
    pub(crate) fn entries(&self) -> &[TileEntry] {
        &self.entries
    }

    /// Per-library tile aggregates, in one pass over the dictionary.
    ///
    /// One pass rather than a scan per statistic per library: the post-pass needs
    /// template counts, tile counts and collision rates for every library, and
    /// computing them separately made it `O(entries x libraries)`.
    pub(crate) fn library_tiles(&self, num_libs: u32) -> Vec<LibraryTiles> {
        let mut stats = vec![LibraryTiles::default(); num_libs as usize];
        for entry in &self.entries {
            if let Some(library) = stats.get_mut(entry.library as usize) {
                library.templates += entry.templates;
                library.tiles += 1;
            }
        }
        // `q` needs the totals, so shares can only be formed on a second look.
        for entry in &self.entries {
            if let Some(library) = stats.get_mut(entry.library as usize)
                && library.templates > 0
            {
                let share = entry.templates as f64 / library.templates as f64;
                library.collision_rate += share * share;
                library.shares.push(share);
            }
        }
        stats
    }
}

/// One spilled pair-template: its dedup signature plus the tile it was imaged on.
///
/// **`off` cannot be dropped as recoverable from `sig`.** They carry different
/// information: `off` is the 2D partition cell (`bin_num * 2 + strand` per end,
/// into `stride²` cells) and `sig` is the pair of *within-bin* positions. Two
/// templates are the same duplicate signature only if both agree.
///
/// The **on-disk** width is fixed at [`SPILL_RECORD_BYTES`] by [`Self::to_bytes`],
/// independent of how this struct happens to be laid out in memory.
///
/// The field order still matters, but for memory rather than disk: the post-pass
/// loads a whole bucket into a `Vec<SpillRecord>` to sort it, so a padding-free
/// 16-byte element keeps that buffer at ~83 MB per bucket for a 333M-template
/// file instead of ~125 MB. Putting the `u64` first is what achieves that —
/// declaring `off` first under `#[repr(C)]` would pad to 24 (4 bytes before `sig`
/// to align it, 4 more at the tail to keep the size an 8-multiple). There is no
/// `repr(C)` here, so rustc reorders fields itself and would reach 16 either way;
/// the declaration order and the assertion below simply make it explicit rather
/// than dependent on a layout rustc leaves unspecified.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct SpillRecord {
    /// Within-cell signature: `(bin_pos1 << 32) | bin_pos2`.
    sig: u64,
    /// Partition cell index.
    off: u32,
    /// Interned `(library, sequencing unit, tile)` ID.
    id: u32,
}

/// Bytes one [`SpillRecord`] occupies on disk.
const SPILL_RECORD_BYTES: usize = 16;

const _: () = assert!(
    size_of::<SpillRecord>() == SPILL_RECORD_BYTES,
    "SpillRecord must stay padding-free so a bucket's sort buffer stays 16 B/record"
);

const _: () = assert!(
    BUCKET_BUFFER_BYTES.is_multiple_of(SPILL_RECORD_BYTES),
    "observe_pair flushes a bucket on len == BUCKET_BUFFER_BYTES exactly, so the buffer must \
     hold a whole number of records"
);

impl SpillRecord {
    /// Encode to its on-disk form.
    ///
    /// Explicit little-endian conversion rather than a cast of the struct's own
    /// bytes: it needs no `unsafe`, no alignment assumptions about the read
    /// buffer, and the cost is invisible next to sorting the bucket.
    #[inline]
    fn to_bytes(self) -> [u8; SPILL_RECORD_BYTES] {
        let mut bytes = [0u8; SPILL_RECORD_BYTES];
        bytes[..8].copy_from_slice(&self.sig.to_le_bytes());
        bytes[8..12].copy_from_slice(&self.off.to_le_bytes());
        bytes[12..].copy_from_slice(&self.id.to_le_bytes());
        bytes
    }

    /// Decode from its on-disk form.
    #[inline]
    fn from_bytes(bytes: &[u8; SPILL_RECORD_BYTES]) -> Self {
        let sig = u64::from_le_bytes(bytes[..8].try_into().expect("8 bytes"));
        let off = u32::from_le_bytes(bytes[8..12].try_into().expect("4 bytes"));
        let id = u32::from_le_bytes(bytes[12..16].try_into().expect("4 bytes"));
        Self { sig, off, id }
    }
}

/// The sequencing/library split for one library.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Decomposition {
    /// Duplicate pair-templates the spill accounted for: `Σ (k − 1)`. Should
    /// equal the run's `duplicate_pairs` for this library.
    pub duplicate_pairs: u64,
    /// Chance-corrected sequencing duplicates.
    pub corrected_sequencing_duplicates: u64,
    /// The residual `duplicate_pairs − sequencing_duplicates`, so the two always
    /// sum to the total exactly.
    pub library_duplicates: u64,
    /// Uncorrected `Σ (k − n_tiles)`. Kept because the gap between this and the
    /// corrected figure is the diagnostic for large-group data: negligible on
    /// WGS, material wherever groups are large. Also the figure the
    /// per-sequencing-unit table sums to exactly.
    pub raw_sequencing_duplicates: u64,
    /// `q = Σ w_t²`, the chance two unrelated templates share a tile.
    pub tile_collision_rate: f64,
    /// Distinct tiles seen for this library.
    pub tile_count: usize,
}

/// Per-sequencing-unit rollup: the QC view that exposes loading differences
/// between flowcells and lanes.
///
/// Worth its own granularity because the variation inside one sample is large.
/// Three flowcells of one library measured sequencing-duplicate rates of 26.8%,
/// 13.1% and 2.6% of their own templates — a 9× spread that a per-library number
/// averages away completely.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SequencingUnitStats {
    /// Library bucket this unit's reads belong to.
    pub library: u32,
    /// The unit token as it appeared in the read names, e.g. `H72CFDSXF:2`.
    pub unit: String,
    /// Templates observed on this unit.
    pub templates: u64,
    /// Distinct tiles seen on this unit.
    pub tiles: usize,
    /// Sequencing duplicates on this unit's own tiles.
    ///
    /// Exact, with no attribution heuristic: each tile of a duplicate group seeded
    /// one molecule and copied the rest, so it contributes `members - 1` to
    /// whichever unit holds it. A group straddling two units therefore splits
    /// across them by construction, and summing this column over every unit gives
    /// back `raw_sequencing_duplicates` precisely.
    pub sequencing_duplicates: u64,
}

/// Streams every pair-template's `(signature, tile)` to disk, then reassembles
/// duplicate groups from it to split duplicates into sequencing and library
/// components.
///
/// A group's tiles cannot be counted in one streaming pass: the moment a
/// template is recognised as a duplicate, the first member of its group has
/// already gone past, and holding every signature's first tile in memory is
/// exactly the per-signature RAM this avoids. So every pair-template is spilled
/// and the groups are rebuilt afterwards. Measured alternatives — a per-tile
/// cache (needs tile-contiguous input, and never sees a whole group, so it can
/// neither apply the chance correction nor break results down per sequencing
/// unit) and keeping only each signature's first tile (87.3% accurate for 250 MB
/// of RAM) — were both dominated by this.
pub(crate) struct TileSpiller {
    /// Interns triples and accumulates the tile shares the correction needs.
    dictionary: TileDictionary,
    /// One append-only buffered file per bucket. Emptied by
    /// [`TileSpiller::finish_spill`], which must run before the buckets are read.
    buckets: Vec<SpillBucket>,
    /// Owns the bucket files; removes them when dropped.
    dir: TempDir,
    /// `buckets.len()`, cached as a `u32` for the hot-path modulo.
    bucket_count: u32,
    /// Records appended so far, for reporting how much temp space was used.
    spilled: u64,
    /// zstd level the bucket files were written with; `None` means raw records.
    /// Retained only so the post-pass knows *whether* to decode — zstd recovers the
    /// level itself from each frame header.
    level: Option<i32>,
}

impl TileSpiller {
    /// Create a spiller writing bucket files under `tmp_dir` (the system temp
    /// directory when `None`).
    ///
    /// `buckets` is clamped to what the process's descriptor limit allows.
    /// `bin_count` is validated here so the hot path can narrow a cell index to
    /// `u32` without a per-template check.
    pub(crate) fn new(
        format: ReadNameFormat,
        bin_count: u32,
        buckets: u32,
        tmp_dir: Option<&Path>,
        level: Option<i32>,
    ) -> Result<Self> {
        let stride = u64::from(stride_for(bin_count));
        if stride * stride > u64::from(u32::MAX) {
            bail!(
                "partition cell count {} exceeds what a spill record can address; \
                 lower --min-bins",
                stride * stride
            );
        }

        let bucket_count = clamp_buckets(buckets);
        let dir = match tmp_dir {
            Some(dir) => TempDir::new_in(dir),
            None => TempDir::new(),
        }
        .context("creating temp directory for the duplicate-decomposition spill")?;

        let mut writers = Vec::with_capacity(bucket_count as usize);
        for bucket in 0..bucket_count {
            let path = dir.path().join(format!("spill-{bucket:04}"));
            let sink = SpillSink::create(&path, level)?;
            writers.push(SpillBucket { buf: Vec::with_capacity(BUCKET_BUFFER_BYTES), sink });
        }

        Ok(Self {
            dictionary: TileDictionary::new(format),
            buckets: writers,
            dir,
            bucket_count,
            spilled: 0,
            level,
        })
    }

    /// Record one both-ends-mapped template: intern its tile and append the
    /// spill record.
    ///
    /// Only pairs reach here, which is deliberate — single-end and orphan
    /// signatures are too noisy to decompose and library size is not estimated
    /// from them. It also means the tile shares behind the chance correction are
    /// measured over exactly the population that forms the groups.
    ///
    /// A failed write is a hard error rather than a silently dropped metric: if
    /// the temp volume is full there is a good chance the output volume is too,
    /// and the user should hear about it while they can still act.
    pub(crate) fn observe_pair(&mut self, library: u32, name: &[u8], slot: PairSlot) -> Result<()> {
        // The split is on by default, so a read name it cannot parse stops a run
        // the user may not have realised was doing this work at all. The error has
        // to name both ways forward.
        let id = self.dictionary.observe(library, name).context(
            "cannot split sequencing from library duplicates. Pass --read-name-format to name \
             this platform's read-name layout, or --sequencing-duplicate-detection off to skip the split",
        )?;
        // Narrowing is checked in `new` via the partition cell count.
        let record = SpillRecord { sig: slot.sig, off: slot.off as u32, id };
        let bucket_idx = self.bucket_of(record);
        let bucket = &mut self.buckets[bucket_idx];
        // An open-coded buffer rather than a `BufWriter`: the 16-byte append is
        // a fixed-size store the compiler inlines, where `BufWriter::write_all`
        // paid a `memmove` call per record.
        if bucket.buf.len() == BUCKET_BUFFER_BYTES {
            bucket.sink.write_all(&bucket.buf).context(
                "writing to the duplicate-decomposition spill (is the temp volume full?)",
            )?;
            bucket.buf.clear();
        }
        bucket.buf.extend_from_slice(&record.to_bytes());
        self.spilled += 1;
        Ok(())
    }

    /// Rebuild duplicate groups from the spill and decompose them.
    ///
    /// Call only after the output BAM has been closed: this reads back gigabytes
    /// and can take tens of seconds, and dupblaster sits in pipelines where a
    /// downstream sort is blocked on its stdout.
    pub(crate) fn decompose(mut self, num_libs: u32) -> Result<DecompositionResult> {
        self.finish_spill(false)?;

        let mut walker = GroupWalker::new(&self.dictionary, num_libs);
        let mut records: Vec<SpillRecord> = Vec::new();
        for bucket in 0..self.bucket_count {
            let path = self.dir.path().join(format!("spill-{bucket:04}"));
            read_bucket(&path, &mut records, self.level)?;
            sort_bucket(&mut records, &walker.library_of, num_libs);
            walker.walk(&records);
        }
        Ok(walker.finish())
    }

    /// Which bucket a record belongs in.
    ///
    /// Hashes the whole group key rather than taking `off` modulo the bucket
    /// count. Either keeps a group intact — its members share `off` *and* `sig` —
    /// but `off` alone distributes terribly, because it encodes the bin number and
    /// a file spanning few bins has only a handful of distinct values. A
    /// single-chromosome extract put the entire spill in one bucket, raising peak
    /// RSS by 125 MB; hashing `sig` too spreads evenly no matter how few bins the
    /// input covers, since `sig` varies per locus.
    #[inline]
    fn bucket_of(&self, record: SpillRecord) -> usize {
        let mixed = (record.sig ^ u64::from(record.off).wrapping_mul(GOLDEN_RATIO_64))
            .wrapping_mul(GOLDEN_RATIO_64);
        // The high half, because a multiply leaves the most entropy there.
        ((mixed >> 32) % u64::from(self.bucket_count)) as usize
    }

    /// Close every bucket stream, optionally totalling the bytes they occupy on
    /// disk.
    ///
    /// Separate from [`Self::decompose`] because a compressed bucket cannot be
    /// decoded until its frame epilogue is written — a correctness requirement, not
    /// just a hook for reporting the size. `measure` is opt-in because totalling
    /// costs one `stat` per bucket and the answer is already known when nothing
    /// compressed the records.
    ///
    /// Safe to call twice, which is what lets [`Self::decompose`] guarantee the
    /// streams are closed without knowing whether the caller already did it. It
    /// accepts no further records afterwards.
    pub(crate) fn finish_spill(&mut self, measure: bool) -> Result<Option<u64>> {
        for (bucket, mut writer) in std::mem::take(&mut self.buckets).into_iter().enumerate() {
            writer.sink.write_all(&writer.buf).with_context(|| {
                format!("flushing spill bucket {bucket} (is the temp volume full?)")
            })?;
            writer.sink.finish().with_context(|| {
                format!("closing spill bucket {bucket} (is the temp volume full?)")
            })?;
        }
        if !measure {
            return Ok(None);
        }
        let mut on_disk = 0;
        for bucket in 0..self.bucket_count {
            let path = self.dir.path().join(format!("spill-{bucket:04}"));
            on_disk += std::fs::metadata(&path)
                .with_context(|| format!("sizing spill bucket {}", path.display()))?
                .len();
        }
        Ok(Some(on_disk))
    }

    /// Logical bytes handed to the spill: `records × 16`, before compression.
    pub(crate) fn spilled_bytes(&self) -> u64 {
        self.spilled * SPILL_RECORD_BYTES as u64
    }
}

/// Everything the decomposition produces.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DecompositionResult {
    /// One entry per library bucket, parallel to the run-summary library rows.
    pub libraries: Vec<Decomposition>,
    /// One entry per distinct `(library, sequencing unit)`.
    pub units: Vec<SequencingUnitStats>,
}

/// Walks sorted spill buckets, accumulating the decomposition as it goes.
struct GroupWalker {
    /// Library per triple ID. One indirection per sort comparison, over a table
    /// small enough to stay in cache — groups must not span libraries, and only
    /// the ID knows which library a record belongs to.
    library_of: Vec<u32>,
    /// Sequencing-unit index per triple ID.
    unit_of: Vec<u32>,
    /// `E[n_tiles | k]` per library.
    models: Vec<ChanceModel>,
    /// Running tallies per library.
    totals: Vec<GroupTotals>,
    /// Tile aggregates per library, computed once up front.
    library_tiles: Vec<LibraryTiles>,
    /// Per-unit rollup, indexed by the values in `unit_of`.
    units: Vec<SequencingUnitStats>,
}

impl GroupWalker {
    /// Build the lookup tables and empty accumulators.
    fn new(dictionary: &TileDictionary, num_libs: u32) -> Self {
        let mut unit_index: HashMap<(u32, &[u8]), u32> = HashMap::new();
        let mut units: Vec<SequencingUnitStats> = Vec::new();
        let mut unit_of = Vec::with_capacity(dictionary.entries().len());
        for entry in dictionary.entries() {
            let unit = *unit_index.entry((entry.library, &entry.unit)).or_insert_with(|| {
                units.push(SequencingUnitStats {
                    library: entry.library,
                    unit: String::from_utf8_lossy(&entry.unit).into_owned(),
                    templates: 0,
                    tiles: 0,
                    sequencing_duplicates: 0,
                });
                units.len() as u32 - 1
            });
            units[unit as usize].templates += entry.templates;
            units[unit as usize].tiles += 1;
            unit_of.push(unit);
        }
        let library_tiles = dictionary.library_tiles(num_libs);
        Self {
            library_of: dictionary.entries().iter().map(|entry| entry.library).collect(),
            unit_of,
            models: library_tiles
                .iter()
                .map(|library| ChanceModel::new(library.shares.clone()))
                .collect(),
            totals: vec![GroupTotals::default(); num_libs as usize],
            library_tiles,
            units,
        }
    }

    /// Accumulate every duplicate group in one sorted bucket.
    ///
    /// Records arrive sorted by `(library, off, sig, id)`, so a duplicate group is
    /// a maximal run of equal `(library, off, sig)`, and within it a tile is a run
    /// of equal `id`.
    fn walk(&mut self, records: &[SpillRecord]) {
        let mut start = 0;
        while start < records.len() {
            let key = self.group_key(&records[start]);
            let mut end = start + 1;
            while end < records.len() && self.group_key(&records[end]) == key {
                end += 1;
            }
            let group = &records[start..end];
            start = end;

            // A signature seen once is not a duplicate group. Skipping it also
            // keeps the inverse solve off the path of every unique template.
            let k = group.len() as u64;
            if k < 2 {
                continue;
            }

            // A tile holding `members` of the group seeded one molecule and copied
            // the rest, so it contributed exactly `members - 1` sequencing
            // duplicates. Summed over the group's tiles that is `k - tiles`, the
            // group's own total — so the per-unit column is exact integer
            // arithmetic that reconciles with the per-library raw figure, with no
            // attribution heuristic and nothing to round.
            let mut tiles = 0u64;
            let mut raw_sequencing = 0u64;
            let mut run = 0;
            while run < group.len() {
                let id = group[run].id;
                let mut run_end = run + 1;
                while run_end < group.len() && group[run_end].id == id {
                    run_end += 1;
                }
                let members = (run_end - run) as u64;
                tiles += 1;
                raw_sequencing += members - 1;
                let unit = self.unit_of[id as usize] as usize;
                self.units[unit].sequencing_duplicates += members - 1;
                run = run_end;
            }

            let library = key.0 as usize;
            // Chance-correct by asking how many *independent* molecules the
            // observed tile count implies, rather than by subtracting the
            // collisions expected of all `k` members. Only the independent
            // molecules can collide, so the latter over-corrects — by 21% for a
            // 3,102-member group on real tile shares.
            let independent = self.models[library].independent_molecules(tiles).min(k as f64);
            self.totals[library].duplicates += k - 1;
            self.totals[library].raw_sequencing += raw_sequencing;
            self.totals[library].corrected_sequencing += k as f64 - independent;
        }
    }

    /// The `(library, cell, signature)` a record groups under.
    #[inline]
    fn group_key(&self, record: &SpillRecord) -> (u32, u32, u64) {
        (self.library_of[record.id as usize], record.off, record.sig)
    }

    /// Resolve the accumulators into the reported result.
    fn finish(self) -> DecompositionResult {
        let libraries = self
            .totals
            .iter()
            .zip(&self.library_tiles)
            .map(|(totals, tiles)| totals.finish(tiles))
            .collect();
        let mut units = self.units;
        // Interning order is input-order dependent; sort so the emitted table is
        // not.
        units.sort_by(|a, b| a.library.cmp(&b.library).then_with(|| a.unit.cmp(&b.unit)));
        DecompositionResult { libraries, units }
    }
}

/// Running per-library tallies while groups are walked.
#[derive(Clone, Copy, Debug, Default)]
struct GroupTotals {
    /// `Σ (k − 1)` over groups.
    duplicates: u64,
    /// `Σ (k − n_tiles)` over groups: the uncorrected count.
    raw_sequencing: u64,
    /// `Σ (E[n_tiles | k] − n_tiles)` over groups, accumulated as a float
    /// because the per-group correction is fractional. Individual terms may go
    /// slightly negative when a group's tiles happen to spread more than chance
    /// predicts; that is real information and is not clamped away per group.
    corrected_sequencing: f64,
}

impl GroupTotals {
    /// Resolve the tallies into the reported [`Decomposition`].
    ///
    /// The sequencing count is rounded and bounded by the duplicate total, and
    /// the library count is then taken as the residual, so the two always sum to
    /// the total exactly however the float accumulation landed.
    fn finish(&self, tiles: &LibraryTiles) -> Decomposition {
        let sequencing = (self.corrected_sequencing.round().max(0.0) as u64).min(self.duplicates);
        Decomposition {
            duplicate_pairs: self.duplicates,
            corrected_sequencing_duplicates: sequencing,
            library_duplicates: self.duplicates - sequencing,
            raw_sequencing_duplicates: self.raw_sequencing,
            tile_collision_rate: tiles.collision_rate,
            tile_count: tiles.tiles,
        }
    }
}

/// Infers how many independent molecules a group's observed tile count implies.
struct ChanceModel {
    /// Tile shares `w_t` for this library.
    shares: Vec<f64>,
    /// Independent-molecule count per distinct *observed tile count*, cached.
    ///
    /// The key is the tile count alone, not the group size: the inverse of
    /// `E[n_tiles | m]` depends only on how many tiles were seen and on the
    /// shares. Nearly every group occupies one or two tiles, so this collapses to
    /// a handful of solves per run.
    implied: HashMap<u64, f64>,
}

impl ChanceModel {
    /// Build the model from `library`'s tile shares.
    fn new(shares: Vec<f64>) -> Self {
        Self { shares, implied: HashMap::new() }
    }

    /// The number of independent molecules that `observed` distinct tiles implies.
    ///
    /// Some of a group's molecules land on the same tile by chance, so `observed`
    /// under-states the molecule count and the difference has to be credited to the
    /// library rather than to the flowcell. This inverts
    /// `E[n_tiles | m] = Σ_t (1 − (1 − w_t)^m)` at `observed`, which is what
    /// "how many independent molecules would look like this?" means.
    ///
    /// Inverting is not the same as subtracting `E[n_tiles | k] − observed`, the
    /// obvious-looking form: only the *independent* molecules can collide, and
    /// there are fewer of them than `k`, so subtracting the collisions expected of
    /// all `k` members over-corrects. On real tile shares that costs 0.2% of the
    /// count at `k = 20`, 3.8% at `k = 500` and 21% at `k = 3102` — precisely the
    /// large-group regime the correction exists for.
    ///
    /// `E` is strictly increasing in `m`, so bisection is enough; the caller bounds
    /// the result by the group size, since a group cannot hold more molecules than
    /// members.
    fn independent_molecules(&mut self, observed: u64) -> f64 {
        // Fewer than two tiles carries no information: `E[n_tiles | m]` is then
        // constant, so the inverse is not unique and bisection would return 1,
        // claiming every duplicate was made on the flowcell. Every template had to
        // land on the one tile whether it was clustered or not. Reporting an
        // unbounded molecule count lets the caller's `min(k)` resolve it to "as
        // many molecules as members", i.e. no sequencing duplicates at all.
        if self.shares.len() < 2 {
            return f64::INFINITY;
        }
        let target = observed as f64;
        if let Some(&implied) = self.implied.get(&observed) {
            return implied;
        }
        // `E[n | 1] = Σ w_t = 1`, and `E[n | m] < m` for more than one tile, so
        // the answer is at least `observed`. Grow the upper bound until it brackets.
        let mut low = target;
        let mut high = (target * 2.0).max(2.0);
        while self.expected_tiles(high) < target && high < INVERSE_SEARCH_CEILING {
            high *= 2.0;
        }
        for _ in 0..BISECTION_STEPS {
            let mid = 0.5 * (low + high);
            if self.expected_tiles(mid) < target {
                low = mid;
            } else {
                high = mid;
            }
            if high - low <= 1e-9 * high {
                break;
            }
        }
        let implied = 0.5 * (low + high);
        self.implied.insert(observed, implied);
        implied
    }

    /// `E[n_tiles | m] = Σ_t (1 − (1 − w_t)^m)`, for real-valued `m`.
    fn expected_tiles(&self, molecules: f64) -> f64 {
        self.shares.iter().map(|share| 1.0 - (1.0 - share).powf(molecules)).sum()
    }
}

/// The previous template's packed key, cached to skip the hash probe.
///
/// Holds the packed key (see [`TileDictionary::pack_key`]) rather than the
/// triple's fields: matching is then one length-checked memcmp against the key
/// the current template just packed, and storing is a buffer swap that moves no
/// bytes — the old memo buffer becomes the next template's pack scratch.
#[derive(Default)]
struct Memo {
    /// Packed key of the previously interned triple; empty until the first
    /// intern. A packed key is never empty (it starts with a 4-byte library id
    /// and a 2-byte length prefix), so a fresh memo cannot match a real key.
    key: Vec<u8>,
    id: u32,
}

/// One spill bucket: its pending-record buffer and the file it flushes to.
///
/// The buffer is managed by hand rather than through a `BufWriter` so the
/// per-record 16-byte append inlines to fixed-size stores; records are appended
/// by [`TileSpiller::observe_pair`] and the buffer is drained a full 64 KB at a
/// time.
struct SpillBucket {
    /// Records not yet written to `sink`; capacity [`BUCKET_BUFFER_BYTES`].
    buf: Vec<u8>,
    /// The bucket file, optionally behind zstd.
    sink: SpillSink,
}

/// Write half of a spill bucket.
///
/// An enum rather than `Box<dyn Write>` because finishing a zstd stream consumes
/// the encoder to write its frame epilogue, which a trait object cannot express.
/// Dispatch cost is irrelevant: the buffer in [`SpillBucket`] means a variant is
/// selected once per 64 KB, not once per record.
enum SpillSink {
    Raw(File),
    /// Boxed because the encoder is far larger than a `File`, and an unboxed
    /// variant would inflate every element of the writer vector.
    Zstd(Box<zstd::stream::write::Encoder<'static, File>>),
}

impl SpillSink {
    /// Create the sink for one bucket file, compressing at `level` when given.
    fn create(path: &Path, level: Option<i32>) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("creating spill bucket {}", path.display()))?;
        match level {
            None => Ok(Self::Raw(file)),
            Some(level) => {
                let encoder = zstd::stream::write::Encoder::new(file, level)
                    .with_context(|| format!("starting zstd for {}", path.display()))?;
                Ok(Self::Zstd(Box::new(encoder)))
            }
        }
    }

    /// Close the stream, writing the zstd frame epilogue if there is one.
    ///
    /// Must happen before the bucket is read back: without the epilogue the
    /// decoder rejects the frame as truncated.
    fn finish(self) -> Result<()> {
        match self {
            Self::Raw(mut file) => file.flush().context("flushing a spill bucket")?,
            Self::Zstd(encoder) => {
                encoder.finish().context("finishing a zstd spill bucket")?;
            }
        }
        Ok(())
    }
}

impl Write for SpillSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Raw(file) => file.write(buf),
            Self::Zstd(encoder) => encoder.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Raw(file) => file.flush(),
            Self::Zstd(encoder) => encoder.flush(),
        }
    }
}

/// Read half of a spill bucket, mirroring [`SpillSink`].
enum SpillSource {
    Raw(File),
    Zstd(Box<zstd::stream::read::Decoder<'static, std::io::BufReader<File>>>),
}

impl SpillSource {
    /// Open one bucket file for the post-pass. `level` only selects the decoder;
    /// zstd reads its own parameters from the frame header.
    fn open(path: &Path, level: Option<i32>) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("opening spill bucket {}", path.display()))?;
        match level {
            None => Ok(Self::Raw(file)),
            Some(_) => {
                let decoder = zstd::stream::read::Decoder::new(file)
                    .with_context(|| format!("starting zstd decode for {}", path.display()))?;
                Ok(Self::Zstd(Box::new(decoder)))
            }
        }
    }
}

impl Read for SpillSource {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Raw(file) => file.read(buf),
            Self::Zstd(decoder) => decoder.read(buf),
        }
    }
}

/// Sort one bucket into the `(library, off, sig, id)` order [`GroupWalker::walk`]
/// reads: a duplicate group is a run of equal `(library, off, sig)`, and sorting
/// by ID last puts a group's records for one tile together, so its distinct tiles
/// are countable in a single pass with no scratch.
///
/// A single-library run — the overwhelmingly common case — holds the library term
/// constant, leaving exactly `(off, sig, id)`: 32 + 64 + 32 bits, which pack into
/// one `u128` whose numeric order *is* the tuple's lexicographic order. That
/// compares in a couple of instructions rather than walking a four-field tuple,
/// and drops the `library_of` lookup that would otherwise run twice per
/// comparison. A whole-genome bucket holds millions of records and the sort costs
/// `O(n log n)` comparisons, so the narrower key is worth the branch.
fn sort_bucket(records: &mut [SpillRecord], library_of: &[u32], num_libs: u32) {
    if num_libs == 1 {
        records.sort_unstable_by_key(packed_sort_key);
    } else {
        records.sort_unstable_by_key(|record| {
            (library_of[record.id as usize], packed_sort_key(record))
        });
    }
}

/// The `(off, sig, id)` portion of the sort key, packed into one `u128` whose
/// numeric order is the tuple's lexicographic order. Both arms of
/// [`sort_bucket`] key through this, so the field order is defined exactly once.
#[inline]
fn packed_sort_key(record: &SpillRecord) -> u128 {
    (u128::from(record.off) << 96) | (u128::from(record.sig) << 32) | u128::from(record.id)
}

/// Read every record of a bucket file into `records`, replacing its contents.
fn read_bucket(path: &Path, records: &mut Vec<SpillRecord>, level: Option<i32>) -> Result<()> {
    records.clear();
    let mut source = SpillSource::open(path, level)?;
    let mut buffer = vec![0u8; BUCKET_READ_RECORDS * SPILL_RECORD_BYTES];
    loop {
        let filled = fill_buffer(&mut source, &mut buffer)
            .with_context(|| format!("reading spill bucket {}", path.display()))?;
        if filled == 0 {
            break;
        }
        if filled % SPILL_RECORD_BYTES != 0 {
            bail!(
                "spill bucket {} ends mid-record ({filled} bytes is not a multiple of \
                 {SPILL_RECORD_BYTES}) — the temp volume may have filled",
                path.display()
            );
        }
        // The bail above guarantees `filled` is a whole number of records, so
        // the `as_chunks` remainder is empty.
        let (chunks, _) = buffer[..filled].as_chunks::<SPILL_RECORD_BYTES>();
        records.extend(chunks.iter().map(SpillRecord::from_bytes));
        if filled < buffer.len() {
            break;
        }
    }
    Ok(())
}

/// Fill `buffer` from `source`, returning how many bytes were read. Short only at
/// end of input, so a record is never split across two calls.
///
/// The loop is load-bearing under compression: a decoder returns whatever the
/// current frame block holds, so short reads are routine rather than exceptional.
fn fill_buffer(source: &mut impl Read, buffer: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        match source.read(&mut buffer[filled..]) {
            Ok(0) => break,
            Ok(read) => filled += read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

/// Clamp a requested bucket count to what this process's descriptor limit allows.
fn clamp_buckets(requested: u32) -> u32 {
    let requested = requested.max(1);
    let Some(limit) = open_file_limit() else {
        return requested;
    };
    let usable = u32::try_from(limit.saturating_sub(RESERVED_DESCRIPTORS)).unwrap_or(u32::MAX);
    let clamped = usable.min(requested).max(1);
    if clamped < requested {
        log::warn!(
            "reducing duplicate-decomposition spill buckets from {requested} to {clamped}: the \
             open-file limit is {limit}. Raise it with `ulimit -n` for larger buckets."
        );
    }
    clamped
}

/// The process's soft limit on open file descriptors.
///
/// Read rather than assumed: the usual defaults (256 on macOS, 1,024 on Linux)
/// bracket the bucket count closely enough that guessing would either waste
/// buckets or exhaust descriptors, and a dev machine can report 1,048,576.
fn open_file_limit() -> Option<u64> {
    let mut limit: libc::rlimit = unsafe { std::mem::zeroed() };
    // SAFETY: `getrlimit` writes only into `limit`, which is a valid, fully
    // initialized `rlimit` for the duration of the call, and reports failure
    // through its return value rather than leaving `limit` untouched.
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut limit) } == 0 {
        Some(limit.rlim_cur as u64)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A dictionary over the colon-delimited (Illumina/Element) layout.
    fn dictionary() -> TileDictionary {
        TileDictionary::new("illumina".parse().expect("illumina is a valid format"))
    }

    /// Read name on `flowcell`/`lane`/`tile`, with the other fields fixed.
    fn name(flowcell: &str, lane: u32, tile: u32) -> Vec<u8> {
        format!("A00354:1305:{flowcell}:{lane}:{tile}:1027:1986").into_bytes()
    }

    #[test]
    fn one_tile_interns_to_one_id() {
        let mut dict = dictionary();
        let first = dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        let second = dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        assert_eq!(first, second);
        assert_eq!(dict.entries().len(), 1);
    }

    #[test]
    fn distinct_tiles_intern_to_distinct_ids() {
        let mut dict = dictionary();
        let a = dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        let b = dict.observe(0, &name("FC", 1, 1102)).expect("parses");
        assert_ne!(a, b);
        assert_eq!(dict.entries().len(), 2);
    }

    #[test]
    fn same_tile_number_on_two_flowcells_interns_separately() {
        let mut dict = dictionary();
        let a = dict.observe(0, &name("H72CFDSXF", 2, 1101)).expect("parses");
        let b = dict.observe(0, &name("22T3L2LT4", 2, 1101)).expect("parses");
        assert_ne!(a, b, "matching tile numbers across flowcells are different places");
    }

    #[test]
    fn same_tile_in_two_libraries_interns_separately() {
        let mut dict = dictionary();
        let a = dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        let b = dict.observe(1, &name("FC", 1, 1101)).expect("parses");
        assert_ne!(a, b, "duplicates are only called within a library");
        assert_eq!(dict.entries()[a as usize].library, 0);
        assert_eq!(dict.entries()[b as usize].library, 1);
    }

    #[test]
    fn ids_are_assigned_serially_from_zero() {
        let mut dict = dictionary();
        for (expected, tile) in [1101, 1102, 1103].into_iter().enumerate() {
            let id = dict.observe(0, &name("FC", 1, tile)).expect("parses");
            assert_eq!(id as usize, expected);
        }
    }

    #[test]
    fn ids_reverse_to_the_original_unit_and_tile_names() {
        let mut dict = dictionary();
        let id = dict.observe(0, &name("H72CFDSXF", 2, 1101)).expect("parses");
        let entry = &dict.entries()[id as usize];
        assert_eq!(&*entry.unit, b"H72CFDSXF:2");
    }

    #[test]
    fn templates_are_counted_per_tile() {
        let mut dict = dictionary();
        for _ in 0..3 {
            dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        }
        dict.observe(0, &name("FC", 1, 1102)).expect("parses");
        assert_eq!(dict.entries()[0].templates, 3);
        assert_eq!(dict.entries()[1].templates, 1);
        assert_eq!(dict.library_tiles(1)[0].templates, 4);
    }

    #[test]
    fn the_memo_does_not_change_which_id_is_returned() {
        // Alternating tiles defeats the one-entry memo on every template; the
        // same input with runs of a tile hits it on nearly every template. Both
        // must intern identically.
        let mut alternating = dictionary();
        let mut runs = dictionary();
        for _ in 0..4 {
            alternating.observe(0, &name("FC", 1, 1101)).expect("parses");
            alternating.observe(0, &name("FC", 1, 1102)).expect("parses");
        }
        for tile in [1101, 1102] {
            for _ in 0..4 {
                runs.observe(0, &name("FC", 1, tile)).expect("parses");
            }
        }
        assert_eq!(alternating.entries().len(), runs.entries().len());
        assert_eq!(alternating.library_tiles(1)[0].templates, runs.library_tiles(1)[0].templates);
        assert_eq!(
            alternating.library_tiles(1)[0].collision_rate,
            runs.library_tiles(1)[0].collision_rate
        );
    }

    #[test]
    fn the_single_library_sort_matches_the_general_one() {
        // Ids stay small so the general path's `library_of` table can cover them;
        // one library throughout, which is the case the fast path assumes.
        let records = vec![
            SpillRecord { sig: 7, off: 2, id: 5 },
            SpillRecord { sig: 7, off: 1, id: 9 },
            SpillRecord { sig: 3, off: 2, id: 1 },
            SpillRecord { sig: 7, off: 2, id: 2 },
            SpillRecord { sig: 3, off: 2, id: 4 },
            SpillRecord { sig: u64::MAX, off: 0, id: 0 },
        ];
        let library_of = vec![0u32; 10];
        let mut packed = records.clone();
        let mut tupled = records;
        sort_bucket(&mut packed, &library_of, 1);
        sort_bucket(&mut tupled, &library_of, 2);
        assert_eq!(packed, tupled);
    }

    #[test]
    fn the_single_library_sort_orders_by_off_then_sig_then_id() {
        // Includes the extremes of every field, which the packed key must carry
        // without the shifts colliding. `library_of` is unread on this path.
        let mut records = vec![
            SpillRecord { sig: 7, off: 2, id: 5 },
            SpillRecord { sig: 7, off: 1, id: 9 },
            SpillRecord { sig: 3, off: 2, id: 1 },
            SpillRecord { sig: 7, off: 2, id: 2 },
            SpillRecord { sig: 3, off: 2, id: 4 },
            SpillRecord { sig: u64::MAX, off: 0, id: 0 },
            SpillRecord { sig: 0, off: u32::MAX, id: u32::MAX },
        ];
        sort_bucket(&mut records, &[], 1);
        let keys: Vec<(u32, u64, u32)> = records.iter().map(|r| (r.off, r.sig, r.id)).collect();
        assert_eq!(
            keys,
            vec![
                (0, u64::MAX, 0),
                (1, 7, 9),
                (2, 3, 1),
                (2, 3, 4),
                (2, 7, 2),
                (2, 7, 5),
                (u32::MAX, 0, u32::MAX),
            ]
        );
    }

    #[test]
    fn the_multi_library_sort_groups_by_library_first() {
        // Ids 0 and 2 are library 1; id 1 is library 0. Library must outrank
        // `off`, so the lone library-0 record sorts ahead of both others despite
        // having the largest `off`.
        let library_of = vec![1, 0, 1];
        let mut records = vec![
            SpillRecord { sig: 1, off: 1, id: 0 },
            SpillRecord { sig: 1, off: 9, id: 1 },
            SpillRecord { sig: 1, off: 2, id: 2 },
        ];
        sort_bucket(&mut records, &library_of, 2);
        assert_eq!(records.iter().map(|r| r.id).collect::<Vec<_>>(), vec![1, 0, 2]);
    }

    #[test]
    fn a_library_on_one_tile_has_a_collision_rate_of_one() {
        let mut dict = dictionary();
        for _ in 0..10 {
            dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        }
        assert_eq!(dict.library_tiles(1)[0].collision_rate, 1.0, "one tile carries no information");
    }

    #[test]
    fn evenly_used_tiles_have_a_collision_rate_of_one_over_n() {
        let mut dict = dictionary();
        for tile in 0..4 {
            dict.observe(0, &name("FC", 1, tile)).expect("parses");
        }
        let q = dict.library_tiles(1)[0].collision_rate;
        assert!((q - 0.25).abs() < 1e-12, "q = {q}");
    }

    #[test]
    fn a_skewed_tile_distribution_raises_the_collision_rate() {
        let mut even = dictionary();
        let mut skewed = dictionary();
        for tile in 0..4 {
            for _ in 0..4 {
                even.observe(0, &name("FC", 1, tile)).expect("parses");
            }
        }
        for _ in 0..13 {
            skewed.observe(0, &name("FC", 1, 0)).expect("parses");
        }
        for tile in 1..4 {
            skewed.observe(0, &name("FC", 1, tile)).expect("parses");
        }
        assert!(
            skewed.library_tiles(1)[0].collision_rate > even.library_tiles(1)[0].collision_rate,
            "{:?} should exceed {:?}",
            skewed.library_tiles(1)[0].collision_rate,
            even.library_tiles(1)[0].collision_rate
        );
    }

    #[test]
    fn collision_rate_is_computed_within_a_library_not_across_libraries() {
        // Each library sits on one tile, so each has q = 1 despite there being
        // two tiles in the dictionary overall.
        let mut dict = dictionary();
        dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        dict.observe(1, &name("FC", 1, 1102)).expect("parses");
        assert_eq!(dict.library_tiles(1)[0].collision_rate, 1.0);
        assert_eq!(dict.library_tiles(2)[1].collision_rate, 1.0);
    }

    #[test]
    fn a_library_with_no_templates_has_no_tiles_and_no_shares() {
        let stats = &dictionary().library_tiles(1)[0];
        assert_eq!(stats.templates, 0);
        assert_eq!(stats.tiles, 0);
        assert!(stats.shares.is_empty());
    }

    #[test]
    fn tiles_are_counted_per_library() {
        let mut dict = dictionary();
        dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        dict.observe(0, &name("FC", 1, 1102)).expect("parses");
        dict.observe(1, &name("FC", 1, 1103)).expect("parses");
        assert_eq!(dict.library_tiles(2)[0].tiles, 2);
        assert_eq!(dict.library_tiles(2)[1].tiles, 1);
    }

    #[test]
    fn an_unparseable_read_name_is_an_error_naming_the_name() {
        let err = dictionary().observe(0, b"SRR1234567.1").expect_err("must not parse");
        assert!(err.to_string().contains("SRR1234567.1"), "{err}");
    }

    /// Bins used by the spiller tests; only the cell-count validation reads it.
    const TEST_BIN_COUNT: u32 = 47;

    /// A spiller over the colon-delimited layout, writing under a temp dir.
    fn spiller() -> TileSpiller {
        TileSpiller::new(
            "illumina".parse().expect("illumina is a valid format"),
            TEST_BIN_COUNT,
            DEFAULT_SPILL_BUCKETS,
            None,
            None,
        )
        .expect("spiller opens")
    }

    /// Observe one pair-template of library 0 at `(off, sig)` on `tile`.
    fn observe(spiller: &mut TileSpiller, off: usize, sig: u64, tile: u32) {
        spiller
            .observe_pair(0, &name("FC", 1, tile), PairSlot { off, sig })
            .expect("observation succeeds");
    }

    /// Decompose a single-library spiller.
    fn decompose(spiller: TileSpiller) -> Decomposition {
        spiller.decompose(1).expect("decomposition succeeds").libraries[0]
    }

    /// Tiles enough that a chance collision between independent molecules is
    /// negligible — the regime real WGS runs in (6,334 tiles were observed on a
    /// three-flowcell sample).
    const DIVERSE_TILES: u32 = 2000;

    /// Give the library a realistic spread of tiles by observing one singleton
    /// template on each of `tiles` tiles.
    ///
    /// Tile shares drive the chance correction, so a test asserting a *corrected*
    /// count has to establish them: with only one tile in the data every
    /// duplicate must land on it, so there is nothing to distinguish clustering
    /// from coincidence and the correct corrected answer is zero. Each
    /// observation gets a unique signature so none of them form groups.
    fn spread_over_tiles(spiller: &mut TileSpiller, tiles: u32) {
        for tile in 0..tiles {
            observe(spiller, 0, 1_000_000 + u64::from(tile), tile);
        }
    }

    #[test]
    fn a_spill_record_round_trips_through_its_on_disk_form() {
        let record = SpillRecord { sig: 0xDEAD_BEEF_1234_5678, off: 9215, id: 6333 };
        assert_eq!(SpillRecord::from_bytes(&record.to_bytes()), record);
    }

    #[test]
    fn a_spill_record_is_sixteen_bytes() {
        // Guards the field order: declaring `off` before `sig` pads it to 24.
        assert_eq!(size_of::<SpillRecord>(), 16);
    }

    #[test]
    fn a_group_entirely_on_one_tile_is_all_sequencing_duplicates() {
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, DIVERSE_TILES);
        for _ in 0..4 {
            observe(&mut spiller, 7, 99, 1101);
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 3);
        assert_eq!(result.raw_sequencing_duplicates, 3);
        assert_eq!(result.corrected_sequencing_duplicates, 3);
        assert_eq!(result.library_duplicates, 0);
    }

    #[test]
    fn a_single_tile_library_reports_no_sequencing_duplicates() {
        // Everything the run saw is on one tile, so `q == 1` and a same-tile
        // duplicate is no evidence of clustering at all. The naive rule would
        // credit every duplicate to sequencing; the correction refuses to, which
        // is why `q` has to be reported alongside the split.
        let mut spiller = spiller();
        for _ in 0..4 {
            observe(&mut spiller, 7, 99, 1101);
        }
        let result = decompose(spiller);
        assert_eq!(result.tile_collision_rate, 1.0);
        assert_eq!(result.raw_sequencing_duplicates, 3);
        assert_eq!(result.corrected_sequencing_duplicates, 0);
    }

    #[test]
    fn a_group_with_one_member_per_tile_is_all_library_duplicates() {
        let mut spiller = spiller();
        for tile in 0..4 {
            observe(&mut spiller, 7, 99, tile);
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 3);
        assert_eq!(result.raw_sequencing_duplicates, 0);
        assert_eq!(result.corrected_sequencing_duplicates, 0);
        assert_eq!(result.library_duplicates, 3);
    }

    #[test]
    fn two_members_on_each_of_two_tiles_is_two_sequencing_and_one_library() {
        // The {A,A,B,B} case. A rule that star-pairs every duplicate to one
        // original scores this 1 (it never compares the two B's to each other);
        // a rule of "has any same-tile partner" scores it 3. The truth is 2.
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, DIVERSE_TILES);
        for tile in [1101, 1102] {
            for _ in 0..2 {
                observe(&mut spiller, 7, 99, tile);
            }
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 3);
        assert_eq!(result.raw_sequencing_duplicates, 2);
        assert_eq!(result.corrected_sequencing_duplicates, 2);
        assert_eq!(result.library_duplicates, 1);
    }

    #[test]
    fn sequencing_and_library_duplicates_always_sum_to_the_duplicate_total() {
        let mut spiller = spiller();
        // A deliberately mixed spill: a singleton, a same-tile pair, a
        // cross-tile pair, and a lopsided larger group.
        observe(&mut spiller, 1, 10, 1101);
        for _ in 0..2 {
            observe(&mut spiller, 2, 20, 1101);
        }
        observe(&mut spiller, 3, 30, 1101);
        observe(&mut spiller, 3, 30, 1102);
        for tile in [1101, 1101, 1101, 1102, 1103] {
            observe(&mut spiller, 4, 40, tile);
        }
        let result = decompose(spiller);
        assert_eq!(
            result.corrected_sequencing_duplicates + result.library_duplicates,
            result.duplicate_pairs
        );
    }

    #[test]
    fn a_signature_seen_once_contributes_no_duplicates() {
        let mut spiller = spiller();
        for sig in 0..8 {
            observe(&mut spiller, 1, sig, 1101);
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 0);
        assert_eq!(result.corrected_sequencing_duplicates, 0);
        assert_eq!(result.library_duplicates, 0);
    }

    #[test]
    fn signatures_differing_only_in_the_partition_cell_are_different_groups() {
        // `off` is not recoverable from `sig`, so a record that dropped it would
        // merge these two into one four-member group.
        let mut spiller = spiller();
        for off in [1, 2] {
            for _ in 0..2 {
                observe(&mut spiller, off, 99, 1101);
            }
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 2, "two groups of two, not one group of four");
    }

    #[test]
    fn groups_do_not_span_libraries() {
        let mut spiller = spiller();
        for library in 0..2 {
            spiller
                .observe_pair(library, &name("FC", 1, 1101), PairSlot { off: 7, sig: 99 })
                .expect("observation succeeds");
        }
        let results = spiller.decompose(2).expect("decomposition succeeds").libraries;
        assert_eq!(results[0].duplicate_pairs, 0, "one template each is not a duplicate");
        assert_eq!(results[1].duplicate_pairs, 0);
    }

    #[test]
    fn the_chance_correction_removes_a_split_that_chance_alone_explains() {
        // Four evenly-loaded tiles and an eight-member group spread two-per-tile.
        // The raw rule scores 4 sequencing duplicates (each tile's members beyond
        // the first), but eight independent molecules would be expected to cover
        // all four tiles anyway — so occupying four tiles is no evidence of
        // clustering, and the corrected count should collapse to nothing.
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, 4);
        for tile in 0..4 {
            for _ in 0..2 {
                observe(&mut spiller, 1, 99, tile);
            }
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 7);
        assert_eq!(result.raw_sequencing_duplicates, 4, "two per tile over four tiles");
        assert_eq!(result.corrected_sequencing_duplicates, 0, "chance explains all of it");
        assert_eq!(result.library_duplicates, 7);
    }

    #[test]
    fn a_group_on_one_of_many_tiles_keeps_its_full_raw_count() {
        // The inverse must not shave the count when the evidence is unambiguous:
        // with plenty of tiles available, all eight members on one tile implies
        // one molecule and seven copies, exactly as the raw rule says. The older
        // subtract-E[n|k] form under-counted here, by 21% at k=3102.
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, 4);
        for _ in 0..8 {
            observe(&mut spiller, 1, 99, 0);
        }
        let result = decompose(spiller);
        assert_eq!(result.raw_sequencing_duplicates, 7);
        assert_eq!(result.corrected_sequencing_duplicates, 7);
    }

    #[test]
    fn the_chance_correction_is_negligible_when_tiles_are_many() {
        // The same eight-member same-tile group, but spread over enough tiles
        // that a chance collision is unlikely: the correction should leave the
        // naive count essentially alone. This is the WGS regime.
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, DIVERSE_TILES);
        for _ in 0..8 {
            observe(&mut spiller, 1, 99, 0);
        }
        let result = decompose(spiller);
        assert_eq!(result.raw_sequencing_duplicates, 7);
        assert_eq!(result.corrected_sequencing_duplicates, 7);
    }

    #[test]
    fn the_reported_split_does_not_depend_on_the_order_templates_arrive_in() {
        // Serial IDs are assigned in first-seen order, so the bytes on disk
        // differ between these two runs while every reported number must not.
        let observations: Vec<(usize, u64, u32)> = vec![
            (1, 10, 1101),
            (1, 10, 1101),
            (1, 10, 1102),
            (2, 20, 1103),
            (2, 20, 1103),
            (3, 30, 1101),
        ];
        let mut forward = spiller();
        for &(off, sig, tile) in &observations {
            observe(&mut forward, off, sig, tile);
        }
        let mut reversed = spiller();
        for &(off, sig, tile) in observations.iter().rev() {
            observe(&mut reversed, off, sig, tile);
        }
        assert_eq!(decompose(forward), decompose(reversed));
    }

    #[test]
    fn the_collision_rate_and_tile_count_are_reported_with_the_split() {
        let mut spiller = spiller();
        for tile in 0..4 {
            observe(&mut spiller, 1, u64::from(tile), tile);
        }
        let result = decompose(spiller);
        assert_eq!(result.tile_count, 4);
        assert!((result.tile_collision_rate - 0.25).abs() < 1e-12);
    }

    #[test]
    fn an_unparseable_read_name_fails_the_spill_rather_than_being_skipped() {
        let mut spiller = spiller();
        let err = spiller
            .observe_pair(0, b"SRR1234567.1", PairSlot { off: 1, sig: 1 })
            .expect_err("must not be silently skipped");
        // `{:#}` walks the context chain; the offending name is on the inner error
        // and the guidance about how to proceed on the outer one.
        let message = format!("{err:#}");
        assert!(message.contains("SRR1234567.1"), "{message}");
        assert!(message.contains("--sequencing-duplicate-detection off"), "{message}");
    }

    #[test]
    fn bucket_count_is_clamped_to_the_descriptor_limit() {
        // Whatever this machine's limit, an absurd request must come back
        // smaller than asked and never zero.
        let clamped = clamp_buckets(u32::MAX);
        assert!(clamped >= 1);
        assert!(clamped < u32::MAX);
    }

    #[test]
    fn a_bucket_count_of_zero_still_yields_one_bucket() {
        assert_eq!(clamp_buckets(0), 1);
    }

    #[test]
    fn every_bucket_count_reaches_the_same_answer() {
        // Bucketing is only a way to keep descriptors bounded; it must not
        // change which records meet each other.
        let observations: Vec<(usize, u64, u32)> =
            (0..64).map(|i| (i as usize % 9, i % 5, (i % 3) as u32)).collect();
        let mut results = Vec::new();
        for buckets in [1, 2, 16, 64] {
            let mut spiller = TileSpiller::new(
                "illumina".parse().expect("valid format"),
                TEST_BIN_COUNT,
                buckets,
                None,
                None,
            )
            .expect("spiller opens");
            for &(off, sig, tile) in &observations {
                observe(&mut spiller, off, sig, tile);
            }
            results.push(decompose(spiller));
        }
        assert!(
            results.windows(2).all(|pair| pair[0] == pair[1]),
            "bucket count changed the answer: {results:?}"
        );
    }

    /// Records per bucket needed to refill the write buffer several times over, so
    /// each compressed stream takes many writes rather than one.
    ///
    /// A real run always works this way; spreading the same count over the default
    /// 64 buckets would leave each one under a single buffer's worth, so every
    /// encoder would see exactly one write and the multi-write path would go
    /// untested.
    const RECORDS_PER_STREAM: u64 = 4 * (BUCKET_BUFFER_BYTES / SPILL_RECORD_BYTES) as u64;
    /// Bucket count for round-trip tests, small so each stream stays busy.
    const BUSY_STREAM_BUCKETS: u32 = 4;

    /// A spiller over the default bucket count, compressing at `level`.
    fn spiller_with_level(level: Option<i32>) -> TileSpiller {
        TileSpiller::new(
            "illumina".parse().expect("valid format"),
            TEST_BIN_COUNT,
            DEFAULT_SPILL_BUCKETS,
            None,
            level,
        )
        .expect("spiller opens")
    }

    /// Round-trip a busy spill at `level` and return the decomposition.
    fn round_trip_at_level(level: Option<i32>) -> Decomposition {
        let mut spiller = TileSpiller::new(
            "illumina".parse().expect("valid format"),
            TEST_BIN_COUNT,
            BUSY_STREAM_BUCKETS,
            None,
            level,
        )
        .expect("spiller opens");
        for i in 0..RECORDS_PER_STREAM * u64::from(BUSY_STREAM_BUCKETS) {
            observe(&mut spiller, i as usize % 9, i % 977, (i % 31) as u32);
        }
        decompose(spiller)
    }

    #[test]
    fn a_fast_tier_round_trips_every_record() {
        // Negative levels select a different internal strategy, not just different
        // parameters, so they need their own coverage.
        assert_eq!(round_trip_at_level(Some(-5)), round_trip_at_level(None));
    }

    #[test]
    fn the_default_level_round_trips_every_record() {
        assert_eq!(round_trip_at_level(Some(3)), round_trip_at_level(None));
    }

    #[test]
    fn the_highest_accepted_level_round_trips_every_record() {
        assert_eq!(round_trip_at_level(Some(9)), round_trip_at_level(None));
    }

    #[test]
    fn zstd_shrinks_a_busy_spill_below_its_logical_size() {
        let mut spiller = TileSpiller::new(
            "illumina".parse().expect("valid format"),
            TEST_BIN_COUNT,
            BUSY_STREAM_BUCKETS,
            None,
            Some(1),
        )
        .expect("spiller opens");
        for i in 0..RECORDS_PER_STREAM * u64::from(BUSY_STREAM_BUCKETS) {
            observe(&mut spiller, i as usize % 9, i % 977, (i % 31) as u32);
        }
        let logical = spiller.spilled_bytes();
        let on_disk = spiller.finish_spill(true).expect("spill closes").expect("size measured");
        assert!(on_disk < logical, "zstd grew the spill: {on_disk} vs {logical} logical");
    }

    #[test]
    fn a_nearly_empty_spill_can_come_out_larger_than_its_logical_size() {
        // Every bucket gets a complete frame whether or not it saw a record, so on a
        // tiny input the per-bucket framing outweighs anything there was to
        // compress. Pinned because the documented guidance warns about it.
        let mut spiller = spiller_with_level(Some(3));
        observe(&mut spiller, 0, 0, 1101);
        let logical = spiller.spilled_bytes();
        let on_disk = spiller.finish_spill(true).expect("spill closes").expect("size measured");
        assert!(
            on_disk > logical,
            "expected framing to dominate one record: {on_disk} vs {logical} logical"
        );
    }

    #[test]
    fn tokens_are_packed_unambiguously_so_a_split_cannot_collide() {
        // "AB"/"C" and "A"/"BC" concatenate to the same bytes; the length prefix
        // is what keeps them distinct.
        let mut dict =
            TileDictionary::new(r"regex:^(?<su>\w+)-(?<tile>\w+)$".parse().expect("valid regex"));
        let a = dict.observe(0, b"AB-C").expect("parses");
        let b = dict.observe(0, b"A-BC").expect("parses");
        assert_ne!(a, b);
        assert_eq!(dict.entries().len(), 2);
    }
}
