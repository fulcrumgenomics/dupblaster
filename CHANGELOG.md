# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `--read-name-format <FORMAT>`: opt-in, off-by-default decomposition of duplicates into a **sequencing** component (copies made on the flowcell — cluster/ExAmp duplicates) and a **library** component (PCR copies, plus genuinely distinct molecules that happen to share a locus). Two duplicates imaged on the same tile are copies of one molecule; on different tiles they are independent. Adds `sequencing_duplicates`, `library_duplicates`, `frac_sequencing_duplicates`, `naive_sequencing_duplicates`, `tile_count` and `tile_collision_rate` to `--stats`, and writes a per-sequencing-unit table beside it (`<STATS>.sequencing-units.tsv`) with per-flowcell/lane template counts, tile counts and sequencing-duplicate rates. Both-ends-mapped pairs only. Validated against read-name ground truth on a 333M-template EM-seq sample: 62.78% of duplicate pairs sequencing, against 62.73% from an independent reference implementation.
  - Accepted formats: `illumina` (`instrument:run:flowcell:lane:tile:x:y` — CASAVA 1.8+, bcl2fastq, BCL Convert), `element` (Element AVITI, the same layout), and `regex:<PATTERN>` with `(?<su>…)` and `(?<tile>…)` capture groups for platforms without a preset. The layout is never guessed — a read name the chosen format cannot parse is an error, because mis-parsing one would silently produce a confident wrong number.
  - **Threshold-free**: uses tile *identity* only, never a pixel radius. Same-tile displacement distributions differ radically between runs (40.7% of same-tile pairs within 20 px on one sample versus 2.8% on another), so a fixed radius cannot generalize — `samtools markdup -d 2500` over-called one of our two test samples and under-called the other. It also handles coincidental duplicates correctly with no special case, which is what makes the metric meaningful for RNA-seq and amplicon data.
  - Duplicate groups are corrected for tiles that collide **by chance**, reported alongside as `tile_collision_rate` (`q = Σ w_t²`). The correction is immaterial on WGS (0.09 pp) but removes over 20% of the naive count for groups above 100 members, so it is load-bearing wherever large groups are normal. A library on a single tile carries no information at all (`q = 1`); its counts are left blank rather than reported as zero, with `tile_count` and `tile_collision_rate` showing why.
  - Costs ~16 bytes of temporary disk per pair (about 5 GB for a 30× human genome, under `--tmp-dir`) and roughly 25–30% wall time, which is why it is opt-in rather than silently on. A pre-flight check fails at startup if the temp volume plainly cannot hold the spill; a write failure mid-run is a hard error rather than a silently dropped metric. The temp spill is read back only after the output BAM is closed, so a downstream process is never left blocked on our stdout.
- `estimated_library_size_corrected` in `--stats`: the Lander-Waterman library-size estimate re-run with sequencing duplicates removed from the observed total (Picard's convention — subtracted from `n`, not from the unique count). Flowcell duplicates are not evidence a library is exhausted, so counting them as saturation understates it; worth 3.1× on one 30× WGS sample (73.0M → 226.4M). Reported alongside the uncorrected `estimated_library_size`, which is unchanged, so the plain column stays comparable across runs.

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
