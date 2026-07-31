# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Duplicates are now split into a sequencing and a library component, on by default.** A **sequencing** duplicate is a copy made on the flowcell (cluster/ExAmp); a **library** duplicate is an independent molecule — a PCR copy, or one that happens to share a locus. They call for opposite responses (loading density versus library complexity) and a single duplicate rate cannot tell them apart. Two duplicates imaged on the same tile are copies of one molecule; on different tiles they are independent.
  - Adds `raw_sequencing_duplicate_pairs`, `corrected_sequencing_duplicate_pairs`, `library_duplicate_pairs`, `frac_duplicate_pairs`, and `frac_sequencing_duplicate_pairs` to `--stats`, and writes a per-sequencing-unit table beside it (`<STATS>.sequencing-units.tsv`) with per-flowcell/lane template counts, tile counts and sequencing-duplicate rates. That table's counts sum *exactly* to `raw_sequencing_duplicate_pairs`.
  - Both-ends-mapped pairs only. **`--no-sequencing-dups`** turns it off; it only ever *classifies* duplicates, so disabling it never changes which reads are marked.
  - **Threshold-free**: tile *identity* only, never a pixel radius. Same-tile displacement distributions differ radically between runs (40.7% of same-tile pairs within 20 px on one sample against 2.8% on another), so no fixed radius generalizes — `samtools markdup -d 2500` over-called one of our two labelled samples and under-called the other. It also handles coincidental duplicates with no special case, which is what makes the metric meaningful on RNA-seq and amplicon data.
  - Corrected for tiles that collide **by chance**, by inferring how many independent molecules a group's observed tile count implies, rather than by subtracting the collisions expected of all its members, which over-corrects. A library on a single tile carries no information at all — every duplicate is on "the" tile whether it was clustered or not — so its counts are left blank rather than reported as zero.
  - Validated against read-name ground truth on a 333M-template EM-seq sample: 62.8% of chr20 duplicate pairs sequencing, 64.9% genome-wide, with a tile count exact to the tile (6,334). The per-sequencing-unit column sums to `raw_sequencing_duplicate_pairs` exactly (53,070,208) rather than approximately.
  - Costs **+7.2% wall time** (359.7 s → 385.6 s on 333M templates), no change in peak memory, and 16 bytes of temporary disk per pair. Only about a third of that (27 ns per pair) falls in the main pass; the rest is a post-pass that runs after the output stream is closed, so it does not hold up a downstream process — the overhead a pipeline feels is roughly **+2.5%**.
- `--read-name-format <FORMAT>`: the read-name layout the split reads. `illumina` (the default — `instrument:run:flowcell:lane:tile:x:y`, CASAVA 1.8+, bcl2fastq, BCL Convert), `element` (Element AVITI, the same layout), or `regex:<PATTERN>` with `(?<su>…)` and `(?<tile>…)` capture groups for platforms without a preset. The layout is never guessed beyond that default, because mis-parsing one would silently produce a confident wrong number.
- Short flags `-l` for `--compression-level` and `-m` for `--add-mate-tags`.
- dupblaster now credits Fulcrum Genomics and links to its repository — two lines on stderr at startup (suppressed by `--quiet`), and a banner at the top of `-h`/`--help`.

### Changed

- **`--help` is far more compact.** Every option's description was rewritten to say the same thing in less prose, and the detail that earns its keep moved into the long form: `-h` is now a one-to-two-line-per-option summary, with the full text under `--help`. Per-option "see the README" pointers are gone — the repository link in the banner is the single pointer to full documentation.
- `--check-crc`, `--no-check-crc`, `--read-buffer-mb`, `--write-buffer-mb`, and `--max-read-length` are grouped under an **"Advanced tuning (rarely needed)"** heading in `--help`, separating the knobs whose defaults suit essentially every run from the options a caller actually chooses between. No behavior change.
- **`estimated_library_size` now has sequencing duplicates removed from the observed total**, following Picard's `ESTIMATED_LIBRARY_SIZE` convention (subtracted from `n`, not from the unique count). Flowcell duplicates are not evidence a library is exhausted, so counting them as saturation understates it — worth 2.3x on one 30x WGS sample. Under `--no-sequencing-dups` the column falls back to the previous uncorrected value.

### Fixed

- A run that aborts part-way now leaves its partial BAM **without** the BGZF EOF marker, so `samtools quickcheck` and other readers correctly report it as truncated. Previously the writer's `Drop` emitted that marker on the way out, so a file missing records passed every integrity check it had. The partial output is still left on disk — discarding it would be its own surprise — it just no longer claims to be complete. Applies to every aborting run, independently of `--no-sequencing-dups`.

### Upgrading from 0.2.0

Two consequences of the split being on by default. Both are one flag to resolve, and `--no-sequencing-dups` restores 0.2.0 behaviour for everything the split touches. (The partial-BAM change under **Fixed** applies regardless of that flag.)

- **Temporary disk is now used on every run** — 16 bytes per both-ends-mapped pair under `--tmp-dir` (`$TMPDIR` by default): ~5 GB for a 30x human genome, ~50 GB at 300x. dupblaster deliberately does not try to predict the requirement, since with a streamed input it cannot know the shape of what is coming; a spill write that fails is a hard error rather than a silently dropped metric. Point `--tmp-dir` at a volume with room, or pass `--no-sequencing-dups`.
- **Read names that are not Illumina/Element-shaped now fail the run.** MGI, Ultima, pre-CASAVA-1.8 Illumina, and anything whose names were rewritten (SRA accessions, simulated data) need either `--no-sequencing-dups` or `--read-name-format regex:PATTERN`. The error message names both.
## [0.2.0] - 2026-07-28

### Added

- `--complexity-metrics <PREFIX>`: opt-in, off-by-default per-library duplication-complexity QC (no extra cost when unset). Writes a **duplicate-rate ladder** (`.duplication-sampled.tsv`/`.pdf`) — duplicate rate vs. sequencing depth, sampled every `--complexity-interval` templates, with matching cumulative and per-window columns — and a **group-size histogram** η_k (`.duplication-spectrum.tsv`/`.pdf`) — how many molecules were seen exactly *k* times. Each library is reported on one category: `pairs` if it has any both-ends-mapped pairs, else `single_end` (the cleaner estimator, and what keeps the metrics correct under every `--single-end-strategy`). PDF plots render via [kuva](https://crates.io/crates/kuva); see the README for details.
- `--complexity-interval <N>`: snapshot cadence (in templates) for the complexity ladder. Default 1,000,000.

### Fixed

- Fixed a deadlock where dupblaster could hang indefinitely — instead of exiting with an error — when the process downstream of it in a pipe died (e.g. `… | dupblaster --output - | sorter` where `sorter` exits on a full disk). After the broken-pipe error was surfaced once, the threaded writer would accept further writes into a ring buffer its already-exited IO thread could never drain, then park forever waiting to be woken. The re-entrant write that the `bgzf` writer performs while being dropped (flushing buffered blocks and the BGZF EOF marker) triggered exactly this. The threaded writer now rejects all writes and flushes once its IO thread has reported an error, so the failure propagates and the process exits non-zero.

## [0.1.1] - 2026-06-24

### Fixed

- Single-end **unmapped** reads (SAM flag `0x4` with no paired bit) no longer
  abort the run with "Can't find first and/or second of pair". They now pass
  through untouched — counted as unmapped orphans and never dup-checked, the
  same as fully-unmapped pairs. Mapped single-end reads were already handled
  correctly; only unmapped single-end reads were affected.

## [0.1.0] - 2026-06-14

Initial release of dupblaster — a fast, streaming duplicate marker for
**query-grouped** SAM/BAM, inspired by samblaster (Faust & Hall, 2014) and
Picard MarkDuplicates.

### Added

- **Streaming duplicate marking** in a single pass over query-grouped input,
  using a strand-aware 5'-aligned signature. Marks duplicates by default;
  `--remove-dups` drops them instead.
- **Library-aware marking, on by default** (matching Picard MarkDuplicates):
  duplicates are called only *within* a library. Library membership comes from
  each read's `RG:Z` tag mapped through the header's `@RG ... LB:` field — read
  groups sharing an `LB` are one library, and reads with no resolvable library
  share an "Unknown Library" bucket. It activates only when the header declares
  more than one distinct `LB`, so single-library runs are byte-for-byte
  identical to single-table mode (no per-read RG scan). The dedup state is
  partitioned into one lazily-allocated table per library, with per-cell
  pre-sizing scaled by `ceil(√library_count)` so the empty-table memory baseline
  grows ~√L rather than linearly. `--library-unaware` forces the single-table,
  library-agnostic behavior (samblaster's behavior).
- **Single-end / orphan strategies** via `--single-end-strategy`:
  - `strand-aware` (default) — a forward orphan and a reverse orphan at the same
    5'-aligned position are *not* duplicates, matching Picard's `fragSort`.
  - `picard-approx` — a fragment-level table registers each end of every
    fully-mapped pair, so later orphans / single-end reads at those positions
    are marked; approximates Picard's "fragments don't beat pairs" rule in one
    streaming pass (order-sensitive: an orphan arriving before its pair passes
    through as non-dup).
  - `picard-exact` — an exact, order-independent implementation: orphans /
    single-end reads are buffered to a temporary uncompressed BAM (`--tmp-dir`,
    default `$TMPDIR`) and marked against a fragment table after the pair pass,
    so an orphan is marked regardless of stream order. Buffered fragments are
    emitted at the end of the output. Matches Picard's fragment dup counts and
    partitions, not its choice of representative read.
  - `samblaster-legacy` — samblaster v0.1.23+'s leftmost-aligned, strand-dropped
    key (not recommended for short-read PE data).
- **`--methylation-mode directional`** — methylation-aware marking for
  directional bisulfite / enzymatic libraries (WGBS, EM-seq, TAPS); off by
  default. Keys each pair in template order (first-of-pair → second-of-pair)
  rather than coordinate-canonically, keeping the two original strands (OT/OB)
  of a fragment distinct while genuine same-strand PCR copies still collapse.
  Works across pair orientations and cross-contig chimeras, and composes with
  `--single-end-strategy`, `--remove-dups`, and `--add-mate-tags`.
  Non-directional / PBAT libraries are out of scope (`--methylation-mode pbat`
  is rejected, not silently mis-handled).
- **Marks based on the current run only** — any pre-existing `FLAG_DUPLICATE`
  bit on input is cleared before marking, matching Picard MarkDuplicates and
  samtools markdup (samblaster instead ORs the old and new flags). Re-running on
  an already-marked BAM reflects only the current pass.
- **SAM and BAM input**, auto-detected from the first byte. **Output is always
  BAM**; `-o` / `--output` must be `-` (stdout) or end in `.bam` (any other
  extension is rejected at startup).
- **`--compression-level <0-12>`** for BGZF output compression. Default 0
  (uncompressed BGZF, same as `samtools view -u`), since most dupblaster
  pipelines pipe into a sort that recompresses; the valid range is delegated to
  `bgzf::CompressionLevel`.
- **`--add-mate-tags`** adds MC (mate CIGAR) / MQ (mate MAPQ) tags;
  **`--ignore-unmated`** tolerates unmated records.
- **Non-query-grouped input detection** — a QNAME block containing only
  secondary / supplementary alignments (no primary) aborts with a clear error
  pointing at probable coordinate-sorted input, rather than silently skipping
  the block (not suppressed by `--ignore-unmated`).
- **`--stats <PATH>`** writes a wide run-summary TSV — one row per library —
  with sample, template / duplicate counts, Picard-style `frac_duplicates`, and
  a Lander-Waterman `estimated_library_size`. A `.gz` / `.bgz` suffix
  transparently gzip-compresses the file; the path must be a real file (not `-`,
  which would interleave with the BAM stream). **`--sample <NAME>`** overrides
  the sample column (otherwise derived from `@RG SM:` tags, comma-joined).
- **Threaded IO** — dedicated read and write threads with lock-free ring buffers
  (sized by `--read-buffer-mb` / `--write-buffer-mb`) decouple the worker from
  kernel-pipe blips in bursty pipelines (`bwa mem | dupblaster | samtools sort`).
- **`--check-crc` / `--no-check-crc`** control BGZF CRC32 verification on input
  (default: on for files, off for stdin).
- **`--max-read-length`** (default 1000) controls synthetic-genome padding.
- **`@PG` provenance** — auto-chains via `PP:` to the existing chain leaf
  (re-running on its own output disambiguates the ID), and validates the input
  `PP:` chain up-front, failing cleanly on a dangling reference rather than
  panicking mid-run.
- **End-of-run resource footer** (wall time, user / system CPU, max RSS) on
  Unix; suppressed by `--quiet`.
- **Reproducible benchmark pipeline** (`benchmark-pipeline/`, Snakemake + pixi):
  downloads an NYGC 1000G high-coverage CRAM (bwa-mem to GRCh38, ~30× WGS),
  subsamples and query-groups it, and times the dup-marking tools in the suite
  (dupblaster's modes, samblaster, Picard MarkDuplicates, samtools markdup,
  dupsifter). The `bench-compare` tool co-streams each tool's output against
  Picard's `kf`-tagged output to produce set-equivalence, orphan-discordance,
  and supplementary-flag-inheritance TSVs.

[Unreleased]: https://github.com/fulcrumgenomics/dupblaster/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/fulcrumgenomics/dupblaster/releases/tag/v0.2.0
[0.1.1]: https://github.com/fulcrumgenomics/dupblaster/releases/tag/v0.1.1
[0.1.0]: https://github.com/fulcrumgenomics/dupblaster/releases/tag/v0.1.0
