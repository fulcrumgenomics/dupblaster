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
//! the input's order — but the equality relations between triples do not, and
//! every metric depends only on those. The reported numbers are therefore
//! order-invariant even though the ID values are not.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
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

/// Rough upper bound on the spill's size as a fraction of the input's, used only
/// for an advisory warning.
///
/// The spill is 16 bytes per pair against roughly 141 bytes per template of a
/// well-compressed input BAM, so about 8.5%; the margin covers tighter
/// compression. It is a poor bound in the other direction — an uncompressed BAM
/// runs several times larger per template, so this over-states the need badly —
/// which is exactly why it only ever warns. See [`warn_if_space_looks_short`].
const SPILL_SIZE_FRACTION_OF_INPUT: f64 = 0.15;

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

/// Interns imaging locations and accumulates per-tile template counts.
pub(crate) struct TileDictionary {
    /// How to pull the unit and tile out of a read name; chosen by the user.
    format: ReadNameFormat,
    /// Packed triple key → serial ID. See [`Self::pack_key`] for the encoding.
    ids: HashMap<Box<[u8]>, u32>,
    /// Interned triples, indexed by ID.
    entries: Vec<TileEntry>,
    /// Reused buffer for the lookup key, so probing allocates nothing.
    key: Vec<u8>,
    /// The previous template's triple and ID.
    memo: Option<Memo>,
    /// Whether [`CARDINALITY_WARN`] has already been reported.
    warned: bool,
}

impl TileDictionary {
    /// Create an empty dictionary that extracts read names using `format`.
    pub(crate) fn new(format: ReadNameFormat) -> Self {
        Self {
            format,
            ids: HashMap::new(),
            entries: Vec::new(),
            key: Vec::new(),
            memo: None,
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
        // Consecutive templates of a name-sorted file almost always share a
        // tile, so this one-entry memo skips both the key build and the hash
        // probe on the common case. It degrades to nothing on coordinate-sorted
        // input, where the map carries the load.
        if let Some(memo) = &self.memo
            && memo.matches(library, location)
        {
            return Ok(memo.id);
        }

        self.pack_key(library, location)?;
        let id = match self.ids.get(self.key.as_slice()) {
            Some(&id) => id,
            None => self.insert(library, location)?,
        };
        self.memo = Some(Memo::new(library, location, id));
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

    /// The chance that two unrelated templates of `library` land on one tile:
    /// `q = Σ w_t²`, over that library's tile shares `w_t`.
    ///
    /// This is the validity indicator for the whole decomposition. A single-tile
    /// or single-unit library has `q ≈ 1` and carries no information — every
    /// duplicate looks like a sequencing duplicate — so the estimate must be
    /// suppressed rather than reported. It is also the basis of the chance
    /// correction for large groups, which collide on tiles by accident.
    ///
    /// `None` when the library has no templates.
    pub(crate) fn collision_rate(&self, library: u32) -> Option<f64> {
        let total = self.template_count(library);
        if total == 0 {
            return None;
        }
        let total = total as f64;
        Some(
            self.entries
                .iter()
                .filter(|entry| entry.library == library)
                .map(|entry| (entry.templates as f64 / total).powi(2))
                .sum(),
        )
    }

    /// Templates observed for `library` across all of its tiles.
    pub(crate) fn template_count(&self, library: u32) -> u64 {
        self.entries
            .iter()
            .filter(|entry| entry.library == library)
            .map(|entry| entry.templates)
            .sum()
    }

    /// Distinct tiles seen for `library`.
    pub(crate) fn tile_count(&self, library: u32) -> usize {
        self.entries.iter().filter(|entry| entry.library == library).count()
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
    fn from_bytes(bytes: &[u8]) -> Self {
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
    pub sequencing_duplicates: u64,
    /// The residual `duplicate_pairs − sequencing_duplicates`, so the two always
    /// sum to the total exactly.
    pub library_duplicates: u64,
    /// Uncorrected `Σ (k − n_tiles)`. Kept because the gap between this and the
    /// corrected figure is the diagnostic for large-group data: negligible on
    /// WGS, over 20% for groups with `k > 100`.
    pub naive_sequencing_duplicates: u64,
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
    /// Sequencing duplicates attributed to this unit.
    ///
    /// A duplicate group can straddle units, and there is no principled way to
    /// split a cluster duplicate between two flowcells — it physically happened
    /// on one. Each group's sequencing duplicates are therefore credited whole to
    /// the unit holding most of its members, which is the same attribution the
    /// reference implementation uses. Exact per-unit quantities (`templates`,
    /// `tiles`) carry no such caveat.
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
    /// One append-only buffered file per bucket.
    buckets: Vec<BufWriter<File>>,
    /// Owns the bucket files; removes them when dropped.
    dir: TempDir,
    /// `buckets.len()`, cached as a `u32` for the hot-path modulo.
    bucket_count: u32,
    /// Records appended so far, for reporting how much temp space was used.
    spilled: u64,
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
        input_bytes: Option<u64>,
    ) -> Result<Self> {
        let stride = u64::from(stride_for(bin_count));
        if stride * stride > u64::from(u32::MAX) {
            bail!(
                "partition cell count {} exceeds what a spill record can address; \
                 lower --min-bins-per-side",
                stride * stride
            );
        }

        let bucket_count = clamp_buckets(buckets);
        let dir = match tmp_dir {
            Some(dir) => TempDir::new_in(dir),
            None => TempDir::new(),
        }
        .context("creating temp directory for the duplicate-decomposition spill")?;
        warn_if_space_looks_short(dir.path(), input_bytes);

        let mut writers = Vec::with_capacity(bucket_count as usize);
        for bucket in 0..bucket_count {
            let path = dir.path().join(format!("spill-{bucket:04}"));
            let file = File::create(&path)
                .with_context(|| format!("creating spill bucket {}", path.display()))?;
            writers.push(BufWriter::with_capacity(BUCKET_BUFFER_BYTES, file));
        }

        Ok(Self {
            dictionary: TileDictionary::new(format),
            buckets: writers,
            dir,
            bucket_count,
            spilled: 0,
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
        let id = self.dictionary.observe(library, name)?;
        // Narrowing is checked in `new` via the partition cell count.
        let record = SpillRecord { sig: slot.sig, off: slot.off as u32, id };
        let bucket = self.bucket_of(record);
        self.buckets[bucket]
            .write_all(&record.to_bytes())
            .context("writing to the duplicate-decomposition spill (is the temp volume full?)")?;
        self.spilled += 1;
        Ok(())
    }

    /// Rebuild duplicate groups from the spill and decompose them.
    ///
    /// Call only after the output BAM has been closed: this reads back gigabytes
    /// and can take tens of seconds, and dupblaster sits in pipelines where a
    /// downstream sort is blocked on its stdout.
    pub(crate) fn decompose(mut self, num_libs: u32) -> Result<DecompositionResult> {
        for (bucket, writer) in self.buckets.iter_mut().enumerate() {
            writer.flush().with_context(|| {
                format!("flushing spill bucket {bucket} (is the temp volume full?)")
            })?;
        }

        let mut walker = GroupWalker::new(&self.dictionary, num_libs);
        let mut records: Vec<SpillRecord> = Vec::new();
        for bucket in 0..self.bucket_count {
            let path = self.dir.path().join(format!("spill-{bucket:04}"));
            read_bucket(&path, &mut records)?;
            // Sorting by ID last puts a group's records for one tile together, so
            // its distinct tiles are countable in a single pass with no scratch.
            records.sort_unstable_by_key(|record| {
                (walker.library_of[record.id as usize], record.off, record.sig, record.id)
            });
            walker.walk(&records);
        }
        Ok(walker.finish(&self.dictionary))
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

    /// Bytes written to the spill, for the run's resource reporting.
    pub(crate) fn spilled_bytes(&self) -> u64 {
        self.spilled * SPILL_RECORD_BYTES as u64
    }
}

/// Everything the decomposition produces.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DecompositionResult {
    /// One entry per library bucket, parallel to the `--stats` library rows.
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
    /// Per-unit rollup, indexed by the values in `unit_of`.
    units: Vec<SequencingUnitStats>,
    /// Unrounded sequencing duplicates per unit, parallel to `units`. Kept as a
    /// float for the same reason the per-library total is, so the two agree.
    unit_sequencing: Vec<f64>,
    /// Reused `(unit, members)` tally for the group in hand. Groups span very few
    /// units, so a short linear scan beats a map.
    tally: Vec<(u32, u64)>,
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
        Self {
            library_of: dictionary.entries().iter().map(|entry| entry.library).collect(),
            unit_of,
            models: (0..num_libs).map(|lib| ChanceModel::new(dictionary, lib)).collect(),
            totals: vec![GroupTotals::default(); num_libs as usize],
            unit_sequencing: vec![0.0; units.len()],
            units,
            tally: Vec::new(),
        }
    }

    /// Accumulate every duplicate group in one sorted bucket.
    ///
    /// Records arrive sorted by `(library, off, sig, id)`, so a duplicate group is
    /// a maximal run of equal `(library, off, sig)`, and within it the distinct
    /// tiles are the runs of equal `id`.
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

            // A signature seen once is not a duplicate group. Skipping it is not
            // only an optimization: it avoids one model lookup for each of the
            // hundreds of millions of unique templates, and avoids accumulating
            // the float residue of `E[n_tiles | 1] − 1`, which is zero only up to
            // the rounding of a sum of shares.
            let k = group.len() as u64;
            if k < 2 {
                continue;
            }
            let tiles =
                1 + group.windows(2).filter(|adjacent| adjacent[0].id != adjacent[1].id).count()
                    as u64;
            let library = key.0 as usize;
            let sequencing = self.models[library].expected_tiles(k) - tiles as f64;
            self.totals[library].duplicates += k - 1;
            self.totals[library].naive_sequencing += k - tiles;
            self.totals[library].corrected_sequencing += sequencing;

            // Accumulate the unrounded value, exactly as the per-library total
            // does, and round once in `finish`. Rounding per group instead would
            // stop the per-unit column summing to the per-library figure —
            // ten groups worth 0.76 each are 8 duplicates, not 10.
            if let Some(unit) = self.majority_unit(group) {
                self.unit_sequencing[unit as usize] += sequencing;
            }
        }
    }

    /// The `(library, cell, signature)` a record groups under.
    #[inline]
    fn group_key(&self, record: &SpillRecord) -> (u32, u32, u64) {
        (self.library_of[record.id as usize], record.off, record.sig)
    }

    /// The sequencing unit holding most of `group`'s members.
    ///
    /// Ties break on the lexicographically smallest unit name rather than on the
    /// unit index, because indices are assigned in first-seen order and would
    /// make the attribution depend on the input's order.
    fn majority_unit(&mut self, group: &[SpillRecord]) -> Option<u32> {
        self.tally.clear();
        for record in group {
            let unit = self.unit_of[record.id as usize];
            match self.tally.iter_mut().find(|(tallied, _)| *tallied == unit) {
                Some((_, members)) => *members += 1,
                None => self.tally.push((unit, 1)),
            }
        }
        let units = &self.units;
        self.tally
            .iter()
            .max_by(|a, b| {
                a.1.cmp(&b.1).then_with(|| units[b.0 as usize].unit.cmp(&units[a.0 as usize].unit))
            })
            .map(|(unit, _)| *unit)
    }

    /// Resolve the accumulators into the reported result.
    fn finish(self, dictionary: &TileDictionary) -> DecompositionResult {
        let libraries = self
            .totals
            .iter()
            .enumerate()
            .map(|(lib, totals)| totals.finish(dictionary, lib as u32))
            .collect();
        let mut units = self.units;
        for (unit, sequencing) in units.iter_mut().zip(&self.unit_sequencing) {
            unit.sequencing_duplicates = sequencing.round().max(0.0) as u64;
        }
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
    /// `Σ (k − n_tiles)` over groups.
    naive_sequencing: u64,
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
    fn finish(self, dictionary: &TileDictionary, library: u32) -> Decomposition {
        let sequencing = (self.corrected_sequencing.round().max(0.0) as u64).min(self.duplicates);
        Decomposition {
            duplicate_pairs: self.duplicates,
            sequencing_duplicates: sequencing,
            library_duplicates: self.duplicates - sequencing,
            naive_sequencing_duplicates: self.naive_sequencing,
            tile_collision_rate: dictionary.collision_rate(library).unwrap_or(0.0),
            tile_count: dictionary.tile_count(library),
        }
    }
}

/// `E[n_tiles | k]` for one library, memoized per group size.
struct ChanceModel {
    /// Tile shares `w_t` for this library.
    shares: Vec<f64>,
    /// `E[n_tiles | k]`, cached per distinct `k`. Group sizes repeat heavily —
    /// nearly every group is small — so this collapses the cost to one
    /// evaluation per distinct size.
    expected: HashMap<u64, f64>,
}

impl ChanceModel {
    /// Build the model from `library`'s tile shares.
    fn new(dictionary: &TileDictionary, library: u32) -> Self {
        let total = dictionary.template_count(library);
        let shares = if total == 0 {
            Vec::new()
        } else {
            dictionary
                .entries()
                .iter()
                .filter(|entry| entry.library == library)
                .map(|entry| entry.templates as f64 / total as f64)
                .collect()
        };
        Self { shares, expected: HashMap::new() }
    }

    /// Distinct tiles expected of `k` templates drawn independently in
    /// proportion to the tile shares: `Σ_t (1 − (1 − w_t)^k)`.
    ///
    /// This is the baseline the observed tile count is measured against. A group
    /// of `k` *independent* library molecules does not occupy `k` tiles — some
    /// collide by chance — so crediting `k − n_tiles` to sequencing duplication
    /// over-counts. The over-count is negligible for small groups and reaches
    /// 21% for groups above 100 members, the regime RNA-seq lives in.
    fn expected_tiles(&mut self, k: u64) -> f64 {
        if let Some(&expected) = self.expected.get(&k) {
            return expected;
        }
        let expected: f64 =
            self.shares.iter().map(|share| 1.0 - (1.0 - share).powf(k as f64)).sum();
        self.expected.insert(k, expected);
        expected
    }
}

/// The previous template's triple, cached to skip the hash probe.
struct Memo {
    library: u32,
    unit: Box<[u8]>,
    tile: Box<[u8]>,
    id: u32,
}

impl Memo {
    fn new(library: u32, location: ImagingLocation<'_>, id: u32) -> Self {
        Self { library, unit: location.unit.into(), tile: location.tile.into(), id }
    }

    /// Whether this memo is for exactly `library` and `location`.
    #[inline]
    fn matches(&self, library: u32, location: ImagingLocation<'_>) -> bool {
        self.library == library && &*self.unit == location.unit && &*self.tile == location.tile
    }
}

/// Read every record of a bucket file into `records`, replacing its contents.
fn read_bucket(path: &Path, records: &mut Vec<SpillRecord>) -> Result<()> {
    records.clear();
    let mut file =
        File::open(path).with_context(|| format!("opening spill bucket {}", path.display()))?;
    let mut buffer = vec![0u8; BUCKET_READ_RECORDS * SPILL_RECORD_BYTES];
    loop {
        let filled = fill_buffer(&mut file, &mut buffer)
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
        records
            .extend(buffer[..filled].chunks_exact(SPILL_RECORD_BYTES).map(SpillRecord::from_bytes));
        if filled < buffer.len() {
            break;
        }
    }
    Ok(())
}

/// Fill `buffer` from `file`, returning how many bytes were read. Short only at
/// end of file, so a record is never split across two calls.
fn fill_buffer(file: &mut File, buffer: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        match file.read(&mut buffer[filled..]) {
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

/// Warn before the main pass if `dir` looks too small to hold the spill.
///
/// Deliberately a warning, not an error. The estimate is only as good as the
/// input's bytes-per-template, which swings by around 4x between an uncompressed
/// BAM and a well-compressed one — so sizing a *hard* failure from it would abort
/// runs that had ample room, which is far worse than not checking. The precise
/// protection is the hard error on a failed spill write; this only buys an early
/// heads-up in the case that error comes too late to be cheap, namely temp and
/// output on different volumes.
///
/// Silent when the input size is unknown (a stream) or the volume cannot be
/// interrogated.
fn warn_if_space_looks_short(dir: &Path, input_bytes: Option<u64>) {
    let (Some(input_bytes), Some(available)) = (input_bytes, available_bytes(dir)) else {
        return;
    };
    let needed = (input_bytes as f64 * SPILL_SIZE_FRACTION_OF_INPUT) as u64;
    if available < needed {
        log::warn!(
            "{} has {:.1} GiB free, which may be too little for the duplicate-decomposition \
             spill of a {:.1} GiB input (a rough upper bound is {:.1} GiB, less for \
             uncompressed input). If it runs out, point --tmp-dir at a larger volume or drop \
             --read-name-format.",
            dir.display(),
            gibibytes(available),
            gibibytes(input_bytes),
            gibibytes(needed),
        );
    }
}

/// Bytes available to this user on the volume holding `dir`.
fn available_bytes(dir: &Path) -> Option<u64> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;

    let path = CString::new(dir.as_os_str().as_bytes()).ok()?;
    let mut stat: libc::statvfs = unsafe { std::mem::zeroed() };
    // SAFETY: `path` is a valid NUL-terminated C string that outlives the call,
    // and `statvfs` only writes into `stat`, reporting failure via its return
    // value rather than leaving `stat` partly written.
    if unsafe { libc::statvfs(path.as_ptr(), &mut stat) } != 0 {
        return None;
    }
    Some(stat.f_bavail as u64 * stat.f_frsize as u64)
}

/// Render a byte count in GiB, for messages about disk space.
fn gibibytes(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0 * 1024.0)
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
        assert_eq!(dict.template_count(0), 4);
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
        assert_eq!(alternating.template_count(0), runs.template_count(0));
        assert_eq!(alternating.collision_rate(0), runs.collision_rate(0));
    }

    #[test]
    fn a_library_on_one_tile_has_a_collision_rate_of_one() {
        let mut dict = dictionary();
        for _ in 0..10 {
            dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        }
        assert_eq!(dict.collision_rate(0), Some(1.0), "one tile carries no information");
    }

    #[test]
    fn evenly_used_tiles_have_a_collision_rate_of_one_over_n() {
        let mut dict = dictionary();
        for tile in 0..4 {
            dict.observe(0, &name("FC", 1, tile)).expect("parses");
        }
        let q = dict.collision_rate(0).expect("has templates");
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
            skewed.collision_rate(0) > even.collision_rate(0),
            "{:?} should exceed {:?}",
            skewed.collision_rate(0),
            even.collision_rate(0)
        );
    }

    #[test]
    fn collision_rate_is_computed_within_a_library_not_across_libraries() {
        // Each library sits on one tile, so each has q = 1 despite there being
        // two tiles in the dictionary overall.
        let mut dict = dictionary();
        dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        dict.observe(1, &name("FC", 1, 1102)).expect("parses");
        assert_eq!(dict.collision_rate(0), Some(1.0));
        assert_eq!(dict.collision_rate(1), Some(1.0));
    }

    #[test]
    fn collision_rate_of_an_unseen_library_is_undefined() {
        assert_eq!(dictionary().collision_rate(0), None);
    }

    #[test]
    fn tiles_are_counted_per_library() {
        let mut dict = dictionary();
        dict.observe(0, &name("FC", 1, 1101)).expect("parses");
        dict.observe(0, &name("FC", 1, 1102)).expect("parses");
        dict.observe(1, &name("FC", 1, 1103)).expect("parses");
        assert_eq!(dict.tile_count(0), 2);
        assert_eq!(dict.tile_count(1), 1);
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
        assert_eq!(result.naive_sequencing_duplicates, 3);
        assert_eq!(result.sequencing_duplicates, 3);
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
        assert_eq!(result.naive_sequencing_duplicates, 3);
        assert_eq!(result.sequencing_duplicates, 0);
    }

    #[test]
    fn a_group_with_one_member_per_tile_is_all_library_duplicates() {
        let mut spiller = spiller();
        for tile in 0..4 {
            observe(&mut spiller, 7, 99, tile);
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 3);
        assert_eq!(result.naive_sequencing_duplicates, 0);
        assert_eq!(result.sequencing_duplicates, 0);
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
        assert_eq!(result.naive_sequencing_duplicates, 2);
        assert_eq!(result.sequencing_duplicates, 2);
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
            result.sequencing_duplicates + result.library_duplicates,
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
        assert_eq!(result.sequencing_duplicates, 0);
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
    fn the_chance_correction_reduces_sequencing_duplicates_when_tiles_are_few() {
        // Four evenly-loaded tiles, and one eight-member group all on tile 0.
        // Eight independent molecules would already be expected to cover only
        // ~3.6 of the four tiles, so crediting all seven duplicates to
        // sequencing over-counts badly.
        let mut spiller = spiller();
        spread_over_tiles(&mut spiller, 4);
        for _ in 0..8 {
            observe(&mut spiller, 1, 99, 0);
        }
        let result = decompose(spiller);
        assert_eq!(result.duplicate_pairs, 7);
        assert_eq!(result.naive_sequencing_duplicates, 7);
        assert!(
            result.sequencing_duplicates < result.naive_sequencing_duplicates,
            "correction should fire: got {}",
            result.sequencing_duplicates
        );
        assert_eq!(
            result.sequencing_duplicates + result.library_duplicates,
            result.duplicate_pairs
        );
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
        assert_eq!(result.naive_sequencing_duplicates, 7);
        assert_eq!(result.sequencing_duplicates, 7);
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
        assert!(err.to_string().contains("SRR1234567.1"), "{err}");
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
    fn a_temp_volume_that_looks_too_small_warns_rather_than_aborting() {
        // The estimate is only as good as the input's bytes-per-template, which
        // swings ~4x with compression, so an over-eager hard failure would abort
        // runs that had ample room. The authoritative check is the write itself.
        let dir = TempDir::new().expect("temp dir");
        warn_if_space_looks_short(dir.path(), Some(u64::MAX));
        warn_if_space_looks_short(dir.path(), Some(1024));
        warn_if_space_looks_short(dir.path(), None);
    }

    #[test]
    fn a_spiller_opens_even_for_an_input_larger_than_the_volume() {
        // Same point at the level that matters: an implausible input size must not
        // stop the spiller from being created.
        assert!(
            TileSpiller::new(
                "illumina".parse().expect("valid format"),
                TEST_BIN_COUNT,
                DEFAULT_SPILL_BUCKETS,
                None,
                Some(u64::MAX),
            )
            .is_ok()
        );
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
