# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/fulcrumgenomics/dupblaster/compare/HEAD
