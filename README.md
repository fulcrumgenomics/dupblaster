[![Build](https://github.com/fulcrumgenomics/dupblaster/actions/workflows/check.yml/badge.svg)](https://github.com/fulcrumgenomics/dupblaster/actions/workflows/check.yml)
[![Version at crates.io](https://img.shields.io/crates/v/dupblaster)](https://crates.io/crates/dupblaster)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/dupblaster.svg?label=bioconda)](https://bioconda.github.io/recipes/dupblaster/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/dupblaster/blob/main/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21445780-blue.svg)](https://doi.org/10.5281/zenodo.21445780)

# dupblaster

A modern, performance-forward successor to [samblaster][samblaster] for
marking and removing PCR duplicates in query-grouped SAM/BAM files —
streaming, BAM-native, threaded IO, and tuned to disappear into the
aligner pipeline.

[samblaster]: https://github.com/GregoryFaust/samblaster

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/dupblaster/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/dupblaster/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/dupblaster/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn
more about how we can power your bioinformatics with dupblaster and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Highlights

- **Coordinate-based duplicate marking, no coordinate sort required.** dupblaster
  reads alignments directly from the aligner's output and marks duplicates in a
  single streaming pass over query-grouped data. This is the same approach
  pioneered by [samblaster][samblaster] (Faust & Hall, [*Bioinformatics* 2014][faust-hall]).
- **CPU-lean in the hot path.** dupblaster minimizes per-record work — SIMD SAM
  parsing, minimal allocations, and a hand-tuned partitioned hash table for the
  coordinate index — so a `bwa-mem | dupblaster | samtools sort` pipeline spends
  its cycles on alignment and sorting rather than on dup-marking.
- **Dedicated IO threads and buffers on both sides of the pipeline.** Reading and
  writing run on their own threads with lock-free ring buffers between them and
  the worker, so dupblaster rarely blocks on IO: a brief stall in `samtools sort`
  (e.g. flushing a sort chunk to disk) does not back-pressure the aligner.
- **BAM-native I/O, no SAM text adapter step.** Input is auto-detected (SAM or
  BAM); output is always BAM (uncompressed by default; see
  `--compression-level`). Pair dupblaster with
  [bwa-mem3][bwa-mem3]'s `--bam=0` flag to skip SAM-text encoding
  end-to-end through the aligner pipeline.
- **Library-aware by default, like Picard MarkDuplicates.** When the header
  declares more than one library (`@RG ... LB:`), duplicates are called only
  *within* a library; single-library inputs are unaffected (and unchanged in
  speed/memory). `--library-unaware` forces samblaster's library-agnostic
  behavior. See [§ Library awareness](#library-awareness).
- **A per-library stats TSV, ready for QC pipelines.** `--stats <PATH>` writes a
  wide TSV — one row per library — with sample, template/duplicate counts,
  Picard-style `frac_duplicates`, and a Lander-Waterman library-size estimate.
- **Opt-in library-complexity QC.** `--complexity-metrics <PREFIX>` adds a duplicate-rate-vs-depth ladder and a group-size histogram (η_k) per library — TSVs plus ready-made PDF plots — to answer "how complex is this library, and would sequencing deeper pay off?". Off by default, and free when unset. See [§ Complexity metrics](#complexity-metrics---complexity-metrics).
- **Sequencing vs. library duplicates, with no pixel threshold.** dupblaster splits duplicates into those made on the flowcell (optical/ExAmp — "the flowcell was loaded too densely") and those from independent molecules ("the library was over-amplified"), which call for opposite responses. It uses imaging-tile *identity* rather than a fixed pixel radius, because same-tile displacement distributions differ radically between runs, and it corrects for tiles that collide by chance — so it stays honest on RNA-seq and amplicon data too. **On by default** at ~7% wall time, of which only a third is felt by a downstream process; `--no-sequencing-dups` turns it off. See [§ Sequencing vs. library duplicates](#sequencing-vs-library-duplicates---no-sequencing-dups).
- **Modern, gnu-style CLI.** `--remove-dups`, `--add-mate-tags`,
  `--ignore-unmated`, `--max-read-length`, `--stats`, … no camelCase flags.

dupblaster is also the fastest option in the suite: on a compute-bound 8× WGS
benchmark it marks duplicates ~14× faster than Picard MarkDuplicates and
`samtools markdup` (run single-threaded) on x86, and ~21–25× faster on
Graviton4, at a fraction of the memory.

<p align="center">
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/dupblaster/main/docs/img/benchmark-walltime.png"
       alt="Duplicate-marking runtime by tool and CPU architecture — 8× WGS, compute-bound, log scale; dupblaster is fastest on both x86 and Graviton4."
       width="760">
</p>

See [§ Benchmarks](#benchmarks) for the full per-architecture tables and
methodology, and [§ Functional equivalence](#functional-equivalence) for
concordance with Picard MarkDuplicates.

[faust-hall]: https://doi.org/10.1093/bioinformatics/btu314
[bwa-mem3]: https://github.com/fg-labs/bwa-mem3

**Jump to:** [Install](#install) · [Quick start](#quick-start) · [Recipes](#recipes) · [Input assumptions](#important-assumptions) · [CLI summary](#cli-summary) · [Stats](#stats-output---stats) · [Complexity metrics](#complexity-metrics---complexity-metrics) · [Algorithm](#algorithm-sketch) · [Benchmarks](#benchmarks) · [Limitations](#limitations)

## Install

```sh
# From crates.io (recommended for Rust users):
cargo install dupblaster

# Via bioconda (recommended for genomics pipelines):
conda install -c bioconda dupblaster

# From source:
git clone https://github.com/fulcrumgenomics/dupblaster.git
cd dupblaster
cargo build --release
# binary at target/release/dupblaster
```

## Quick start

Drop dupblaster into the standard align → mark-dups → sort pipeline,
directly after the aligner:

```sh
# With bwa-mem3 (recommended — emits uncompressed BAM with --bam=0,
# skipping the SAM-text round trip entirely):
bwa-mem3 mem --bam=0 -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --stats sample.dupblaster.tsv -o - \
    | mako sort -o sample.bam -

# Or with samtools sort:
bwa-mem3 mem --bam=0 -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --stats sample.dupblaster.tsv -o - \
    | samtools sort -@ 4 -o sample.bam -

# Or with classic bwa-mem (SAM output; dupblaster auto-detects):
bwa mem -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --stats sample.dupblaster.tsv -o sample.dups.bam
```

The pipeline above pairs dupblaster with bwa-mem3[^bwa-mem3] and mako[^mako].

[^bwa-mem3]: [bwa-mem3](https://github.com/fg-labs/bwa-mem3) is Fulcrum Genomics'
    bwa-mem successor; its `--bam=0` flag emits uncompressed BAM directly, so the
    aligner → dupblaster → sorter pipeline skips SAM-text encoding end to end.
[^mako]: [mako](https://github.com/fg-labs/mako) is Fulcrum Genomics' fast SAM/BAM
    sorter, a drop-in replacement for `samtools sort` for the common cases.

## Recipes

All of these read query-grouped SAM/BAM and write BAM; the flags compose
freely. Multi-library inputs need no flag — dupblaster splits on `@RG LB:`
automatically (see [§ Library awareness](#library-awareness)).

```sh
# Remove duplicates instead of flagging them (leaner BAM out):
bwa-mem3 mem --bam=0 -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --remove-dups -o - \
    | mako sort -o sample.bam -

# Add MC (mate CIGAR) and MQ (mate MAPQ) tags, which some downstream callers
# and UMI tools expect:
bwa-mem3 mem --bam=0 -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --add-mate-tags -o - \
    | mako sort -o sample.bam -

# Exact, order-independent orphan handling (Picard's "fragments don't beat
# pairs"). Orphans are emitted at the end of the stream, so sort downstream:
bwa-mem3 mem --bam=0 -t 8 ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --single-end-strategy picard-exact --tmp-dir /scratch -o - \
    | mako sort -o sample.bam -

# Bisulfite / EM-seq / TAPS (directional preps): keep the two original strands
# (OT/OB) of each fragment distinct so methylation isn't lost to dup-collapsing.
# Use any bisulfite aligner that emits query-grouped output (e.g. bwa-meth):
bwa-meth.py --reference ref.fa r1.fq.gz r2.fq.gz \
    | dupblaster --methylation-mode directional -o - \
    | mako sort -o sample.bam -
```

See [§ Methylation mode](#methylation-mode) for what `directional` does and
why non-directional / PBAT libraries are out of scope.

## Important assumptions

Two assumptions that, if violated, either fail the run or produce wrong answers.

### Input must be **query-grouped**

Every record for a given QNAME must appear in one contiguous run, with no other
QNAME's records interleaved; the order of QNAMEs relative to each other doesn't
matter. In SAM/BAM terms this is `@HD SO:unsorted GO:query` (equivalently
`grouporder=queryname`), which `bwa-mem`, `bwa-mem3`, `bwa-mem2`, and `bowtie2`
emit naturally.

dupblaster makes a best effort to catch coordinate-sorted input and fail loudly:

- **Paired-end:** the first QNAME's block holds only one mate-half, so dupblaster
  aborts with an "unmated record" error on record one.
- **Chimeric or multi-mapped reads** (essentially all modern WGS): a block
  eventually contains only secondary/supplementary alignments — its primary sits
  at a different coordinate — triggering a "QNAME … but no primary" error.
- **Undetectable:** pure single-end data with no secondary/supplementary
  alignments. Every block legitimately holds one primary, so the ordering is
  invisible and coordinate-sorted input silently produces wrong dup calls.

To dedupe a coordinate-sorted BAM, either re-sort to query-grouped
(`mako sort --queryname` or `samtools sort -n`), or use a coordinate-sort-aware
tool (`samtools markdup`, Picard `MarkDuplicates`).

### Output is always BAM (uncompressed by default)

dupblaster writes only BAM. The `-o`/`--output` path must be `-` (stdout) or end
in `.bam`; any other extension is rejected at startup rather than producing a
misnamed file.

Output defaults to uncompressed BGZF (level 0, "stored" blocks, same as
`samtools view -u`): the dominant pipeline pipes into a sort step that
recompresses anyway, so skipping the round-trip saves wall time. Use
`--compression-level <0-12>` when writing to durable storage or to a sink that
won't recompress.

## CLI summary

The most common flags:

| Flag | Purpose |
|---|---|
| `-i, --input <PATH>` | Input SAM/BAM file (default: stdin). |
| `-o, --output <PATH>` | Output BAM file or `-` for stdout. Must end in `.bam`. |
| `-r, --remove-dups` | Drop duplicate reads from output instead of just flagging them. |
| `-m, --add-mate-tags` | Add MC (mate CIGAR) and MQ (mate MAPQ) tags to all paired records. |
| `--ignore-unmated` | Don't abort if a primary record's mate is missing. |
| `-l, --compression-level <N>` | BGZF compression level for output (0-12). Default 0 = uncompressed. |
| `--single-end-strategy <NAME>` | How to key single-end / orphan reads. `strand-aware` (default), `picard-exact`, `picard-approx`, or `samblaster-legacy`. See [§ Single-end / orphan handling](#single-end--orphan-handling). |
| `--methylation-mode <MODE>` | Methylation-aware keying for bisulfite / enzymatic-conversion data. Off by default. `directional` keeps the two original strands (OT/OB) of a fragment distinct. See [§ Methylation mode](#methylation-mode). |
| `--tmp-dir <DIR>` | Directory for dupblaster's temporary files, deleted on exit (default: `$TMPDIR`). **The sequencing-vs-library duplicate split writes here on every run** — 16 bytes per pair, so ~5 GB for a 30x human genome and ~50 GB at 300x. Point this at a volume with room, or pass `--no-sequencing-dups`. Also used by `--single-end-strategy picard-exact` for its orphan buffer. |
| `--tmp-compression-level <LEVEL>` | Compress dupblaster's temporary files with zstd, `-7` to `9` (default: uncompressed) — both the duplicate-split spill and the `picard-exact` orphan buffer. On a whole-genome run, levels `-5` to `1` cut the spill by roughly a third to a half for a few percent more runtime; the orphan buffer compresses better still. The saving is smaller on smaller inputs, and memory grows with the level (see below). Negative levels are zstd's fast tiers. Storage only — duplicate flags and metrics are identical either way. |
| `--library-unaware` | Disable library-aware marking; use one dedup table across all reads (samblaster behavior). No effect when the header has ≤1 library. See [§ Library awareness](#library-awareness). |
| `--stats <PATH>` | Write a per-library TSV of run-summary metrics (one row per library). |
| `--sample <NAME>` | Override the `sample` column in `--stats` (and `--complexity-metrics`) output. |
| `--complexity-metrics <PREFIX>` | Write per-library duplication-complexity QC: a duplicate-rate ladder and a group-size histogram. Off by default. See [§ Complexity metrics](#complexity-metrics---complexity-metrics). |
| `--complexity-interval <N>` | Snapshot cadence (in templates) for the complexity ladder. Default 1,000,000. |
| `--no-sequencing-dups` | Don't split duplicates into sequencing (on-flowcell) and library components. The split is on by default; it classifies duplicates rather than changing which reads are marked. See [§ Sequencing vs. library duplicates](#sequencing-vs-library-duplicates---no-sequencing-dups). |
| `--read-name-format <FORMAT>` | Read-name layout the split takes the sequencing unit and tile from: `illumina` (default), `element`, or `regex:PATTERN`. A name the layout cannot parse is a hard error. |

Run `dupblaster --help` for the full list. The IO ring buffers
(`--read-buffer-mb`, `--write-buffer-mb`), BGZF CRC verification
(`--check-crc` / `--no-check-crc`), and `--max-read-length` sit under an
"Advanced tuning" heading there — their defaults suit essentially every run.

## Stats output (`--stats`)

`--stats <PATH>` writes a TSV with one column per metric and **one row per
library** — in library-aware runs, one row for each `@RG LB:` that saw data
(the catch-all "Unknown Library" row appears only if it did); single-library
runs emit a single row. This is the format we use in our own QC pipelines —
easy to concatenate across many samples and load into pandas / data.table /
DuckDB. Give `--stats` a `.gz` (or `.bgz`) suffix to gzip-compress the file.

| Column | Meaning |
|---|---|
| `sample` | From `--sample`, or comma-joined `@RG SM:` values, or empty. |
| `library` | Library name (`@RG LB:`), `Unknown Library`, or `All Reads` under `--library-unaware`. |
| `dupblaster_version` | Version of dupblaster that produced this row. |
| `total_templates` | Templates (QNAMEs) seen with a usable primary. |
| `duplicate_templates` | Templates marked as duplicates of another template. |
| `frac_duplicates` | Picard-style read-level fraction: `(orphan_dups + 2*pair_dups) / (orphan_reads + 2*pair_reads)`. |
| `mapped_pairs` | Templates with both reads mapped. |
| `duplicate_pairs` | Mapped pairs marked duplicate. |
| `mapped_orphans` | Templates with exactly one read mapped. |
| `duplicate_orphans` | Mapped orphans marked duplicate. |
| `unmapped_orphans` | Templates with one read present and unmapped (no mapped mate). |
| `unmapped_pairs` | Templates with both reads unmapped. |
| `unmated_templates` | Templates with a stray half (skipped unless `--ignore-unmated`). |
| `estimated_library_size` | Lander-Waterman estimate of the library's distinct molecules, with sequencing duplicates removed from the observed total (Picard's `ESTIMATED_LIBRARY_SIZE` convention). Empty when not estimable — including the degenerate case where *every* duplicate is a sequencing duplicate, which leaves no resampling to infer from. |
| `raw_sequencing_duplicate_pairs` | Duplicate pairs made on the flowcell, uncorrected: `Σ (tile members − 1)` over duplicate groups. The per-sequencing-unit table sums to exactly this. Empty under `--no-sequencing-dups`, or when a library has no pairs or sits on a single tile. |
| `corrected_sequencing_duplicate_pairs` | The same count corrected for tiles that collide by chance. **This is the figure to use.** |
| `library_duplicate_pairs` | The residual `duplicate_pairs − corrected_sequencing_duplicates`; the two sum exactly. |
| `frac_duplicate_pairs` | `duplicate_pairs / mapped_pairs` — the pair-level rate, as distinct from the read-level `frac_duplicates`. |
| `frac_sequencing_duplicate_pairs` | `corrected_sequencing_duplicates / duplicate_pairs`. |

## Sequencing vs. library duplicates (`--no-sequencing-dups`)

Not all duplicates mean the same thing. A **library duplicate** is a second copy of a molecule that existed before sequencing — a PCR product, or two genuinely distinct molecules that happen to share a locus — and it tells you the library was over-amplified or under-complex. A **sequencing duplicate** is made *on the flowcell*, when one cluster is read as two (optical duplicates on unpatterned flowcells, ExAmp duplicates on patterned ones). It tells you the flowcell was loaded too densely and says nothing at all about the library.

They call for opposite responses, and a single duplicate rate cannot distinguish them. dupblaster does, using the one place the input records where a read was imaged: its name.

```bash
# On by default — nothing to enable.
dupblaster -i in.bam -o out.bam --stats sample.dupblaster.tsv
```

The split is **on by default**, assumes the Illumina read-name layout (see [§ Choosing the format](#choosing-the-format)), and covers **both-ends-mapped pairs only** — single-end and orphan reads are not decomposed. Pass `--no-sequencing-dups` to turn it off; it only ever *classifies* duplicates, so turning it off never changes which reads get marked. Measured on a 333M-template (50 GB) name-sorted BAM:

| | without | with | cost |
|---|---|---|---|
| wall time | 359.7 s | 385.6 s | **+7.2%** |
| user CPU | 358.6 s | 381.4 s | +6.3% |
| peak RSS | 5,419 MB | 5,386 MB | **none** |
| temp disk | — | 4.95 GiB | 16 B per pair |

So the cost is a few percent of runtime and no extra memory, against 16 bytes of temp disk per pair (about 5 GB for a 30x human genome — see `--tmp-dir`).

If that temp space is the binding constraint rather than CPU, **`--tmp-compression-level`** compresses the spill with zstd, without changing a single reported number. On the same 333M-template sample, against an uncompressed control:

| level | spill | ratio | wall | Δ wall |
|---|---|---|---|---|
| *(omitted)* | 4.95 GiB | 1.00x | 391.2 s | — |
| `-5` | 2.95 GiB | 1.68x | 391.2 s | +1.4%¹ |
| `-1` | 2.56 GiB | 1.94x | 407.0 s | +4.0% |
| `1` | 2.44 GiB | 2.03x | 405.9 s | +3.7% |

¹ Measured in a separate sweep against its own control. Repeated uncompressed runs varied by 1.7% (384.7 / 386.0 / 391.2 s), so treat differences of a few seconds between levels as noise — the ratios reproduce exactly, the timings do not separate the levels.

The same flag also compresses the orphan buffer that `--single-end-strategy picard-exact` writes. That one compresses better — **3.0x at `-5` to 4.4x at `1`**, measured on query-grouped BAM records, because it holds whole records rather than dense integer keys — and it costs the main pass nothing, since it compresses on the IO thread that file already runs on. On well-paired data it is a rounding error either way (orphans were 0.03% of reads on the sample above); it matters for single-end or orphan-heavy input, where that buffer is the bulk of the temporary footprint.

Two things to know before reaching for a high level. The spill is split across 64 bucket files that compress independently, so **smaller inputs compress less** — a whole-genome run gets the ratios above, a small panel may get nothing. And **memory grows with the level**, because every spill bucket holds its own compression context: negligible at level 1, but hundreds of megabytes by level 9. Levels above 9 are rejected for that reason. If a run is memory-constrained as well as disk-constrained, stay at or below 1.

**In a pipeline, most of that cost is invisible.** The work splits into a per-pair part during the main pass and a post-pass that reassembles the groups, and the post-pass runs only *after* the output stream is closed:

| phase | cost | blocks a downstream process? |
|---|---|---|
| main pass — one dictionary lookup and a 16-byte append per pair (27 ns/pair) | +9.1 s | yes |
| post-pass — read back 4.95 GiB, sort each bucket, walk groups | +16.9 s | **no** |

Those two account for the whole 25.9 s above.

In `bwa-mem … \| dupblaster … \| samtools sort`, the sort sees end-of-input as soon as the main pass ends and proceeds while dupblaster is still decomposing. The overhead the pipeline actually feels is the first row alone — about **+2.5%**.

### How it works

Two duplicates imaged on the same tile are copies of one molecule. On different tiles they are independent, because one cluster cannot span two tiles. So for a duplicate group of `k` templates spread over `n` distinct tiles:

```
sequencing_duplicates = k - n      one molecule seeds each tile; the rest are copies
library_duplicates    = n - 1      each extra tile is an independent molecule
                        ─────
total                 = k - 1      the duplicates dupblaster already reports
```

The two always sum to `duplicate_pairs` exactly. That is the **raw** rule, reported as `raw_sequencing_duplicate_pairs`.

#### Correcting for chance

Independent molecules sometimes land on the same tile by accident, so the tile count `n` under-states how many molecules a group really held, and the raw rule credits the difference to the flowcell. The reported `corrected_sequencing_duplicate_pairs` therefore does not use `k - n`. It asks instead **how many independent molecules would produce the `n` tiles we saw**, by inverting

```
E[n | m] = Σ_t (1 - (1 - w_t)^m)
```

for `m` (where `w_t` is tile `t`'s share of the library's templates), and reports `k - m`, bounded by `k`.

Inverting is not the same as subtracting `E[n | k] - n`, which is the form that first suggests itself. Only the *independent* molecules can collide, and there are fewer of them than `k`, so subtracting the collisions expected of all `k` members over-corrects — by 0.2% of the count at `k = 20`, 3.8% at `k = 500` and **21% at `k = 3102`** on real tile shares. The gap between the raw and corrected columns is a useful diagnostic in its own right: it is 0.05 pp on our WGS sample and grows with group size, so it is where RNA-seq and amplicon data will differ most from WGS.

A library on fewer than two tiles is a special case: `E[n | m]` is then constant, every duplicate is on "the" tile whether it was clustered or not, and no attribution is possible. Both counts are left blank rather than reported as zero.

### No pixel radius

Picard and `samtools markdup` call a duplicate optical when two reads fall within a fixed pixel distance. dupblaster uses tile **identity** only, and reads no x/y coordinates at all, because same-tile displacement distributions differ radically between runs. On one of our two labelled samples 40.7% of same-tile pairs were within 20 px and 2.5% were beyond 2500 px; on the other, 2.8% were within 20 px and **10.7%** were beyond 2500 px. A single radius cannot serve both — `markdup -d 2500` over-called one and under-called the other:

| | sample A | sample B |
|---|---|---|
| truth (tile identity) | 62.8% of dups | 71.5% |
| `markdup -d 2500` | 59.4% | 65.0% |
| `markdup -d 100` | 33.2% | 10.3% |

Working from identity also handles coincidental duplicates with no special case, which is what makes the metric meaningful for RNA-seq and amplicon data where such duplicates are the norm. Two distinct molecules at one locus are imaged independently, so they usually occupy two tiles and read as one library duplicate and zero sequencing. They *can* collide on one tile, with probability `q = Σ w_t²`, which is exactly what the correction above exists to account for rather than something the raw rule gets right on its own.

### Choosing the format

The layout defaults to `illumina`, and dupblaster never *guesses* beyond that — read-name formats differ between platforms in ways that **mis-parse rather than fail** (the pre-CASAVA-1.8 Illumina layout puts a y coordinate exactly where the modern one puts a tile), and a confident wrong number is worse than no number.

**A read name the layout cannot parse is a hard error.** dupblaster will not quietly drop a metric that is switched on, so a platform without a preset needs one of two things from you: `--read-name-format regex:PATTERN` to describe its layout, or `--no-sequencing-dups` to skip the split. The error names both.

That applies to MGI, Ultima, pre-CASAVA-1.8 Illumina, and anything whose names were rewritten (SRA accessions, for instance) — see [§ Upgrading](#upgrading-from-020).

| Value | Layout |
|---|---|
| `illumina` | `instrument:run:flowcell:lane:tile:x:y` — CASAVA 1.8+, bcl2fastq, BCL Convert. Extra trailing fields (an appended UMI) are ignored. |
| `element` | Element AVITI, which uses the same seven-field layout. |
| `regex:PATTERN` | A pattern with `(?<su>…)` and `(?<tile>…)` named capture groups, for platforms without a preset — including undelimited layouts such as MGI DNBSEQ, e.g. `regex:^(?<su>F\w+L\d)(?<tile>C\d{3}R\d{3})`. |

The sequencing unit and tile are treated as **opaque tokens** and never parsed as numbers, so an instrument widening its tile numbering cannot break extraction. The unit is flowcell *and* lane, so a tile number reused across flowcells or lanes is correctly seen as a different physical place — a mistake that causes `samtools markdup` to call cross-flowcell duplicates optical ([samtools#1996](https://github.com/samtools/samtools/issues/1996)), which we measured at 1.98% of duplicates on real data.

### When it can't be estimated

The split is left **blank** rather than zero whenever it cannot mean anything, so a missing number is visibly a missing number rather than a measurement. That happens when `--no-sequencing-dups` was passed, when a library has no both-ends-mapped pairs, or when a library sits on a single tile.

That last case is worth understanding: with one tile, every duplicate is on "the" tile whether it was clustered or not, so nothing can be attributed. The per-sequencing-unit file's `tiles` column is what distinguishes it from the other reasons a cell is blank. A misconfigured `--read-name-format` shows up there too, as an implausible tile count — past a million distinct tiles dupblaster warns, and past 16,777,216 it fails.

### Per-sequencing-unit table

With `--stats <PATH>`, a companion `<PATH>.sequencing-units.tsv` breaks the same numbers down per flowcell-and-lane, with real flowcell names. This is where the metric earns its keep — sequencing-duplicate rate varies enormously *within* one sample, so a per-library average hides it. Three flowcells of one library measured 26.0%, 12.6% and **2.4%** of their own templates, and one of those flowcells ranged from 3.9% to 19.1% across its four lanes.

| Column | Meaning |
|---|---|
| `sample`, `library` | As in `--stats`. |
| `sequencing_unit` | Flowcell and lane, verbatim from the read names, e.g. `H72CFDSXF:2`. |
| `templates` | Templates observed on this unit. |
| `tiles` | Distinct tiles seen on this unit. |
| `sequencing_duplicate_pairs` | Sequencing duplicates on this unit's own tiles — each tile contributes `members - 1` of any group it holds part of, so a group straddling two units splits across them exactly. Summed over all units this equals `raw_sequencing_duplicate_pairs`. |
| `frac_sequencing_duplicate_pairs` | `sequencing_duplicates / templates` — the loading-density signal, comparable across units. |

### Corrected library size

Flowcell duplicates are not evidence that a library is exhausted, so counting them as saturation makes it look smaller than it is. `estimated_library_size` therefore runs the Lander-Waterman solve with sequencing duplicates removed from the observed total — Picard's `ESTIMATED_LIBRARY_SIZE` convention, subtracted from `n` and not from the unique count. On one 30x WGS sample that raised the estimate 2.3x. There is deliberately only one such column: under `--no-sequencing-dups`, or wherever the split is not estimable, it falls back to the uncorrected value rather than appearing beside a second one.

### Upgrading from 0.2.0

Two things changed for existing pipelines, both because the split is now on by default:

1. **dupblaster reads temporary disk on every run** — 16 bytes per both-ends-mapped pair under `--tmp-dir` (`$TMPDIR` by default), so ~5 GB for a 30x human genome and ~50 GB at 300x. It does not try to predict the requirement, because with a streamed input it cannot know the shape of what is coming; a spill write that fails is a hard error rather than a silently dropped metric. Where that space is tight, `--tmp-compression-level` trades a little CPU for a substantially smaller spill.
2. **Read names that are not Illumina/Element-shaped now fail the run.** If your BAMs come from MGI or Ultima, predate CASAVA 1.8, or had their names rewritten (SRA accessions, simulated data), add `--no-sequencing-dups` — or describe the layout with `--read-name-format regex:PATTERN` to get the metric.

Both are single-flag fixes; `--no-sequencing-dups` restores 0.2.0 behaviour for everything the split touches.

One change applies regardless of that flag: a run that aborts part-way now leaves its partial BAM **without** the BGZF EOF marker, so `samtools quickcheck` and other readers report it as truncated. Previously such a file looked complete while missing records.

### Known limitation

A cluster duplicate that straddles a tile boundary reads as a library duplicate, because the two copies are on different tiles. We measured this at 1.4-1.9% of duplicates on two samples. The direction is known (it under-estimates sequencing duplicates), and fixing it would need flowcell geometry that Illumina does not publish.

## Complexity metrics (`--complexity-metrics`)

`--complexity-metrics <PREFIX>` opts into extra per-library QC output aimed at the questions "how complex is this library?" and "would sequencing deeper pay off?". It writes **four files** — a TSV and a PDF plot for each of a *duplication-sampled ladder* and a *duplication spectrum (η_k)*:

- `<PREFIX>.duplication-sampled.tsv` / `.pdf`
- `<PREFIX>.duplication-spectrum.tsv` / `.pdf`

It is **off by default**; when unset, dupblaster allocates nothing extra and its speed and memory are unchanged. It works under every `--single-end-strategy`, and `--complexity-interval <N>` sets the ladder's snapshot cadence (default 1,000,000 templates).

### Pairs vs. single-end: one category per library

Each library is reported on **exactly one** category:

- **`pairs`** — the both-ends-mapped subset — if the library has *any* pairs. A pair signature pins *two* genomic endpoints, so coincidental "same coordinates by chance" collisions are rare and the numbers are a clean library-complexity signal. A library's mapped-orphan / single-end reads are a small, noisy minority (a single-end signature pins only one endpoint), so for a paired library they are intentionally left out.
- **`single_end`** — the mapped single-end / orphan signatures — only if the library is *solely* single-end (no pairs at all). This covers single-end runs (e.g. much RNA-seq) where the single-end signatures are all you have.

This split is also what keeps the metrics correct and consistent across every single-end strategy: single-end numbers are only ever reported for a library with no pairs, which is exactly the case where the underlying keying is unambiguous. The `category` column names which subset a row describes.

### `<PREFIX>.duplication-sampled.tsv` — duplication rate vs. depth

Snapshots taken every `--complexity-interval` templates (plus a final row at the true total). Each row carries a **cumulative** set (`total`/`unique`/`duplicates`/`frac_duplicates`) and a matching **window** set (`window_*`) covering just the templates added since the previous snapshot — plot cumulative `frac_duplicates` vs. `total` for the running rate, or `window_frac_duplicates` vs. `total` for the marginal view (usually the more legible one).

| Column | Meaning |
|---|---|
| `sample` | As in `--stats`. |
| `library` | Library name. |
| `category` | `pairs` or `single_end` (see above). |
| `total` | Cumulative templates observed in this category at this snapshot. |
| `unique` | Distinct molecules so far (`total − duplicates`). |
| `duplicates` | Templates flagged as duplicates so far. |
| `frac_duplicates` | Cumulative `duplicates / total`. |
| `window_total` | Templates added *in this window* (since the previous snapshot; equals `total` for the first row). |
| `window_unique` | Distinct molecules in this window (`window_total − window_duplicates`). |
| `window_duplicates` | Duplicates flagged in this window. |
| `window_frac_duplicates` | Marginal rate: `window_duplicates ÷ window_total`. |

> **The ladder's shape depends on input order.** It is cleanest on single-lane or homogeneous-lane input, where it reads as a genuine saturation curve. On coordinate- or queryname-sorted **multi-flowcell** input the curve is dominated by contiguous per-read-group blocks — flowcell/ExAmp (optical) duplicates cluster within a read group and land adjacent under sorting — so the ladder is an order-dependent *diagnostic* (a flowcell-heterogeneity fingerprint), **not** a complexity estimator. For order-independent library complexity, use the duplication spectrum (η_k) below.

### `<PREFIX>.duplication-sampled.pdf` — ladder plot

A ready-made PDF of the **marginal** (per-window) duplication rate vs. depth, one line per library on a shared fraction-duplicated axis so libraries compare directly (a legend appears when there is more than one). The cumulative rate is left to the TSV; the marginal is the more legible view, and the same order-dependence caveat above applies.

### `<PREFIX>.duplication-spectrum.tsv` — group-size histogram (η_k)

For each distinct molecule, how many times it was observed; tallied into "how many molecules were seen exactly once, exactly twice, …". This is the same shape as RSeQC's duplication plot and the input to library-complexity estimators. Internally it is kept cheap with a side table that stores counts only for the duplicate signatures seen ≥2× (singletons are recovered by subtraction), and that table is dropped the moment a library is known to be paired.

> **Per-molecule counts saturate at 65,535** (held as `u16` to keep the side table compact). A molecule observed more than that is reported at `n_observations = 65535`, so for an extremely amplified molecule — e.g. a single 5′ position in very deep RNA-seq — the reads/pairs figures can slightly undercount. The distinct-molecule counts (`n_molecules`) are never affected.

| Column | Meaning |
|---|---|
| `sample` | As in `--stats`. |
| `library` | Library name. |
| `category` | `pairs` or `single_end` (see above). |
| `n_observations` | Occurrence count *k* — how many times each distinct molecule was seen (capped at 65,535; see above). |
| `n_molecules` | How many distinct molecules were observed exactly `n_observations` times. |

### `<PREFIX>.duplication-spectrum.pdf` — η_k plot

A ready-made PDF of the spectrum: the fraction of **molecules** and of **reads/pairs** (weighting each family by its size *k*) observed exactly *k* times, as points on a log-y / linear-x percent plot. The x-axis is trimmed to where the spectrum is dense (extended while ≥50% of occurrence bins in `[1..k]` are populated, capped at *k*=500); molecules and reads beyond that are summarized in the x-axis label rather than clipped. Divergence between the two series — reads riding above molecules in the tail — flags high-copy families (PCR/optical jackpots, or highly-expressed transcripts in RNA) that carry a disproportionate share of the reads. With more than one library the plot is a **faceted grid**, one panel per library on the single page.

## Algorithm sketch

dupblaster uses the same coordinate-based dedup approach as samblaster
(see Faust & Hall, [*Bioinformatics* 2014][faust-hall] for the full
algorithm): for each template, compute a key from the 5'-aligned
positions and strands of the primary first-of-pair and second-of-pair
(or the unpaired primary), and call a template a duplicate if its key
has been seen before. Secondary and supplementary alignments inherit
the duplicate flag from their primary.

The coordinate index is a partitioned hash table sized to the
genome's contig list. The per-contig position cap is **2^31 − 1 bp
(~2.15 Gb)** — this is the SAM/BAM format's own limit, not a dupblaster
constraint. Every realistic reference genome fits well under that:
human chr1 is 0.25 Gb, the largest plant chromosomes (wheat) are
under 1 Gb, axolotl tops out at around 3 Gb chromosomes (which would
hit the cap — file an issue if this affects you). The number of
contigs is unbounded.

Comparisons account for soft-clipping at the 5' end, so two reads
with different clipping but the same true alignment start are
correctly identified as duplicates.

### Single-end / orphan handling

`--single-end-strategy` selects how single-end reads and orphans (a
mapped read whose mate is unmapped) are keyed in the dedup table.
Four strategies are supported:

* **`strand-aware`** (default). The dedup key is the strand-specific
  5'-aligned position. A forward orphan and a reverse orphan at the
  same 5' coord are *not* duplicates. This matches Picard
  MarkDuplicates' `fragSort` keying at the single-end level and is
  the right answer for short-read PE data.
* **`picard-approx`**. Strand-aware key plus a Picard-style cross-
  check: each end of every fully-mapped PE pair is also registered
  in a fragment-level table, so a later orphan / single-end read at
  the same 5' coord is marked as a duplicate of the pair. This
  approximates Picard's "fragments don't beat pairs" rule in a
  single streaming pass — approximate because an orphan that arrives
  *before* its corresponding pair passes through as non-dup. Roughly
  doubles dupblaster's memory footprint at run time. Recommended when
  you want Picard-equivalent dup partitions and have the memory
  headroom.
* **`picard-exact`**. The same "fragments don't beat pairs" rule as
  `picard-approx`, but **exact and order-independent**. dupblaster runs
  two passes: fully-mapped and unmapped pairs stream straight to the
  output, while every mapped-orphan / single-end read is buffered to a
  temporary BAM (see `--tmp-dir`, and `--tmp-compression-level` to compress
  it). After the pair pass, the
  pair table is consumed into a fragment table holding the 5' position of
  every paired read end, and the buffered fragments are re-read and marked
  against it — so an orphan is marked a duplicate of a pair regardless of
  which came first in the stream, unlike `picard-approx`. Two trade-offs:
  (1) buffered fragments are emitted at the **end** of the output stream
  rather than in input order (re-sort downstream if order matters — most
  pipelines already do), and (2) it writes a temporary on-disk copy of the
  fragment reads. **Intended for paired data**, where orphans are a small
  fraction; on single-end-only libraries it would buffer the entire input.
  Matches Picard's fragment *dup counts / partitions* exactly; it does not
  attempt to reproduce Picard's choice of *which* read in a duplicate set
  is the representative.
* **`samblaster-legacy`**. samblaster's post-v0.1.25 (Feb 2020)
  behavior: leftmost-aligned reference coordinate with the strand
  bit dropped, so a forward orphan and a reverse orphan at the same
  leftmost-aligned position collide. Produces false positives on
  short-read PE data (two distinct molecules sharing only their
  leftmost-aligned coord on opposite strands are marked dup).
  Provided for byte-compatibility with samblaster output on long-
  read singleton workflows where it was originally validated; **not
  recommended otherwise**.

### Methylation mode

`--methylation-mode directional` adapts duplicate marking for bisulfite
and enzymatic-conversion libraries (WGBS, EM-seq, TAPS). In these data
the two strands of a fragment — the original-top (OT/CTOT) and
original-bottom (OB/CTOB) — carry **independent** methylation and must be
counted as *separate* molecules, not duplicates of each other.

Standard WGS keying canonicalizes each pair to a coordinate-ordered
signature (leftmost end first). At a given locus the OT and OB fragments
occupy the same two coordinates, so canonicalization gives them the same
key and the second is wrongly marked a duplicate. Directional mode keys
the pair in **template order** instead — first-of-pair into slot A,
second-of-pair into slot B, with no coordinate swap. Because directional
preps ligate adapters to the intact double-stranded fragment *before*
conversion, first-of-pair is locked to a consistent end of the original
strand, so:

- the OT and OB pairs produce *different* keys (their first-of-pair reads
  sit on opposite ends and strands) and are kept distinct, while
- genuine PCR copies of one strand reproduce the same template-order
  geometry and still collapse.

This holds across all pair orientations (FR/RF/FF/RR) and cross-contig
chimeras: orientation is simply part of the key. The single-end / orphan
path is unchanged — it is already strand-aware, so OT/OB orphans stay
separate in every mode.

**Scope.** Only **directional** libraries are supported. Non-directional
/ **PBAT** libraries (where the first-of-pair-to-strand relationship is
not fixed) are out of scope: correct keying there needs per-read strand
tags *plus* a canonicalize-within-strand key, which is a separate piece
of work. `--methylation-mode pbat` is intentionally rejected rather than
silently mis-handled. The flag is also independent of (and composes with)
`--single-end-strategy`, `--remove-dups`, and `--add-mate-tags`.

### Library awareness

Duplicates are a property of a *library* (a PCR-amplified pool), not of the
genome: two reads at the same coordinates from **different** libraries are
independent observations, not copies of one molecule. Picard MarkDuplicates
keys on the library; samblaster ignores it. dupblaster follows Picard **by
default**.

Library membership comes from each read's `RG:Z` tag, mapped through the
header's `@RG ... LB:` field. Read groups that share one `LB` are one library;
reads with no `RG`, an `RG` absent from the header, or an `@RG` line with no
`LB` share a single **"Unknown Library"** bucket. The dedup state is then
partitioned into one independent table per library.

This activates only when the header declares **more than one** distinct `LB`.
With zero or one library, there's nothing to separate, so dupblaster runs in
single-table mode — identical results, and identical speed and memory, to
before (no per-read `RG` scan). Pass `--library-unaware` to force single-table
mode even with a multi-library header (samblaster's behavior).

**Memory.** Each library's table is allocated lazily (only libraries that
actually appear cost anything), and per-cell pre-sizing is scaled down by
`ceil(√library_count)`. The empty-table baseline therefore grows roughly with
the *square root* of the library count rather than linearly, and the stored
signatures — the part that scales with data — are essentially conserved when a
fixed amount of sequencing is split across libraries (identical fragments
rarely recur across independent preps). In practice a 3-library run and a
1-library run over the same data use within a few percent of the same memory.

## Benchmarks

### Overview

Benchmarked on two AWS instance classes — `r8i.4xlarge` (Intel x86) and
`r8g.4xlarge` (Graviton4), each 16 vCPU / 123 GiB RAM. The reproducible
Snakemake pipeline in
[`benchmark-pipeline/`](benchmark-pipeline/) drives the runs and
[`bench-compare/`](bench-compare/) evaluates functional equivalence against
Picard MarkDuplicates. The input is 8× WGS downsampled from NYGC 1000 Genomes
sample HG03953 (bwa-mem aligned, retaining supplementary alignments): a 67 GB
query-grouped BAM, or 83 GB as SAM text.

The streaming tools here — dupblaster and samblaster — process uncompressed data
faster than all but the fastest local SSDs can supply or absorb it. They are
built for a different deployment: a single pipe between an aligner and a sorter,
where the surrounding stages keep the stream moving and storage is never the
bottleneck. Benchmarking them against a disk would measure the disk, not the
tool, and would not reflect how they are actually run.

We model that no-I/O-bottleneck case directly. Each tool's input is pre-warmed
into the page cache so the timed read runs at RAM speed, and each tool writes its
output to a FIFO drained to `/dev/null` (the writer thread still serializes every
record, so its CPU is counted). Format conversion and any sort a tool requires
happen outside the timed window.

The numbers below are dupblaster 0.1.0, a single replicate, taken 2026-06-13; the
raw collated reports are committed in
[`benchmark-pipeline/published-results/`](benchmark-pipeline/published-results/).

### Performance

**x86 — `r8i.4xlarge` (Intel, 16 vCPU)**

| Tool | Runtime (s) | CPU (s) | RSS (MB) | × fastest |
|------|------------:|--------:|---------:|----------:|
| dupblaster 0.1.0 (BAM) | 60.0 | 68.6 | 1248 | 1.0× |
| dupblaster 0.1.0 (picard-approx) | 66.1 | 91.0 | 2314 | 1.1× |
| dupblaster 0.1.0 (picard-exact) | 60.8 | 71.2 | 1466 | 1.0× |
| dupblaster 0.1.0 (SAM) | 93.0 | 125.3 | 1229 | 1.5× |
| samblaster 0.1.26 | 123.4 | 123.4 | 1425 | 2.1× |
| dupsifter 1.3.0 | 203.2 | 203.1 | 3114 | 3.4× |
| Picard 3.4.0 | 825.2 | 1250.1 | 8578 | 13.8× |
| samtools markdup 1.23 | 808.2 | 795.2 | 218 | 13.5× |
| samtools markdup 1.23 (`-S`) | 991.9 | 948.1 | 228 | 16.5× |
| samtools markdup 1.23 (`-m s -S`) | 986.5 | 955.6 | 228 | 16.5× |

**Graviton4 — `r8g.4xlarge` (Arm, 16 vCPU)**

| Tool | Runtime (s) | CPU (s) | RSS (MB) | × fastest |
|------|------------:|--------:|---------:|----------:|
| dupblaster 0.1.0 (BAM) | 40.2 | 60.0 | 1236 | 1.0× |
| dupblaster 0.1.0 (picard-approx) | 67.0 | 96.1 | 2327 | 1.7× |
| dupblaster 0.1.0 (picard-exact) | 44.6 | 62.9 | 1438 | 1.1× |
| dupblaster 0.1.0 (SAM) | 98.1 | 115.0 | 1228 | 2.4× |
| samblaster 0.1.26 | 151.1 | 151.1 | 1424 | 3.8× |
| dupsifter 1.3.0 | 217.0 | 216.9 | 3113 | 5.4× |
| Picard 3.4.0 | 1006.1 | 1504.1 | 8523 | 25.0× |
| samtools markdup 1.23 | 865.2 | 864.0 | 217 | 21.5× |
| samtools markdup 1.23 (`-S`) | 1014.2 | 1012.8 | 227 | 25.2× |
| samtools markdup 1.23 (`-m s -S`) | 1024.3 | 1017.7 | 227 | 25.5× |

`samtools markdup`'s timed window excludes the `fixmate -m` + coordinate-sort
prep it requires; the other tools consume the query-grouped input directly.

dupblaster is ~14× faster than Picard and single-threaded samtools markdup on
x86 and ~21–25× faster on Graviton4. The margin is larger on Graviton4:
dupblaster's lighter per-record work scales across the cores (40 s vs 60 s on
x86), while Picard's heavier per-read cost does not (1006 s vs 825 s). Resident
memory differs similarly — dupblaster holds ~1.2 GB versus Picard's ~8.5 GB.

dupblaster exposes two orthogonal options. **Input format** affects only speed: BAM
is the native on-disk shape, while the SAM path must text-parse the stream
(~1.5–2.4× the wall time, identical marking). **Single-end strategy** affects only
orphan handling (see [§ Single-end / orphan handling](#single-end--orphan-handling)):
the default `strand-aware` is fastest,
`picard-approx` adds a streaming pair-end cross-check, and `picard-exact`
buffers orphans for an order-independent second pass. All three mark paired
reads identically; pick a strategy for orphan fidelity, a format for speed.

### Functional equivalence

| Tool | PE concordance | SE/orphan concordance | Supp. marked |
|------|---------------:|----------------------:|:---------------:|
| Picard 3.4.0 | reference | reference | only on query-grouped input |
| dupblaster 0.1.0 (strand-aware) | 100% | 75.5% | yes |
| dupblaster 0.1.0 (picard-approx) | 100% | 95.9% | yes |
| dupblaster 0.1.0 (picard-exact) | 100% | 100% | yes |
| samblaster 0.1.26 | 100% | 75.6% | yes |
| dupsifter 1.3.0 | 86.9% | 75.5% | yes |
| samtools markdup 1.23 | 99.9% | 99.8% | no |
| samtools markdup 1.23 (`-S`) | 99.9% | 99.8% | yes |
| samtools markdup 1.23 (`-m s -S`) | 100% | 99.8% | yes |

Concordance is *set-equivalence*, not per-read agreement. We tag every Picard
primary with a canonical fragment key (unclipped 5′ position + strand), so a
duplicate set is all templates sharing a key — a pair keys on the sorted keys
of both ends, an orphan on its one mapped end. A tool is concordant on a set
when it marks the **same number** of templates in that set as duplicates as
Picard does: same group, same count, regardless of which template each tool
elects as the representative. PE counts both-ends-mapped sets; SE/orphan counts
one-end-mapped sets. "Supp. marked" is whether the duplicate flag propagates
from a primary to its supplementary alignments (the benchmark data has
supplementary but no secondary alignments).

### Key differences

- **samblaster / dupblaster (strand-aware)** — identical to Picard on every
  paired set. The orphan gap is one effect: Picard additionally cross-checks an
  orphan against the 5′ positions of mapped pairs and marks it a duplicate of a
  pair; the strand-aware key does not. dupblaster's `picard-approx` and
  `picard-exact` strategies add that cross-check (approx in one streaming pass;
  exact order-independently), raising orphan concordance to 95.9% / 100%.

- **[dupsifter][dupsifter]** — paired concordance is only 86.9% because its signature is
  strand-of-origin-aware (inherited from its WGBS design) and does not collapse
  FR/RF orientation. For a fragment captured from both strands at the same
  coordinates, Picard counts one duplicate set, but dupsifter keeps one
  representative per orientation — marking exactly one fewer duplicate in every
  affected set. `-W` (WGS mode) disables only the bisulfite-strand inference,
  not this orientation split, so the gap persists on plain WGS.

- **samtools markdup** — paired concordance is 99.9% in its default `-m t`
  (template) mode; the few disagreements are FF/RR (same-strand) and
  inter-chromosomal pair geometries, where template mode folds R1/R2
  (first/second-in-template) identity into the key and so diverges from
  Picard. Running `-m s` (sequence mode) keys on the unclipped 5′ ends
  without that distinction — coordinate-canonical like Picard and dupblaster —
  and closes the gap to **100%** paired concordance (the `-m s -S` row), at no
  meaningful runtime cost. It cross-checks orphans like Picard (99.8%). By
  default it does *not* mark the secondary/supplementary alignments of a
  duplicate; `-S` enables that, at a runtime cost. This matters for short-read
  structural-variant callers, which read supplementary alignments as breakpoint
  evidence: when the supplementary records of a dup-marked template go
  unflagged, they count as independent observations of what is really one
  PCR-amplified event, inflating breakpoint support. dupblaster (and samblaster)
  propagate the primary's flag to supplementary alignments by default, avoiding
  this; with `samtools markdup` pass `-S` (or re-flag downstream).

### Recommendations

- **Streaming, non-UMI data — use dupblaster.** It is the fastest option and
  reads directly from the aligner with no coordinate sort. Paired marking is
  identical to Picard, and `picard-exact` also matches Picard on orphans.
- **Input not query-grouped — use [Picard][picard-md].** It coordinate-sorts
  internally and adds optical-duplicate detection and per-read duplicate-set
  tags.
- **Coordinate-sorted input that must stay a streaming pass — use
  [samtools markdup][samtools-markdup].** It marks duplicates in a single pass
  over coordinate-sorted BAM without re-grouping. Note its FF/RR and
  inter-chromosomal differences from Picard, and pass `-S` for
  secondary/supplementary marking.
- **UMI libraries — use [fgumi][fgumi].** It distinguishes true PCR duplicates
  from coincidental position collisions; positional dedup cannot.

[picard-md]: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
[samtools-markdup]: http://www.htslib.org/doc/samtools-markdup.html
[dupsifter]: https://github.com/huishenlab/dupsifter
[fgumi]: https://github.com/fulcrumgenomics/fgumi

## Limitations

- No optical / sequencing duplicate detection.
- No UMI awareness — use [fgumi][fgumi] for UMI-aware dedup.
- Methylation mode (`--methylation-mode`) supports **directional**
  libraries only (WGBS / EM-seq / TAPS); non-directional / PBAT is not
  supported.

## Differences from samblaster (C++)

For users coming from samblaster, the high-level changes are:

- **Modern CLI:** GNU-style `--kebab-case` flags (e.g. `--remove-dups`,
  not `--removeDups`).
- **BAM-native:** input is SAM *or* BAM; output is always uncompressed
  BAM. No `samtools view` adapter step.
- **No SV-extraction flags** (`-d`, `-s`, `-u`, `-a`, `-e`): dropped.
- **Library-aware:** duplicates are called within a library (`@RG LB:`)
  by default, like Picard; samblaster is library-agnostic.
  `--library-unaware` restores the samblaster behavior.
- **Larger genomes:** no hardcoded position cap.
- **Threaded IO:** dedicated read and write threads with ring buffers,
  so the worker doesn't block on pipe stalls.
- **Per-library stats TSV:** structured `--stats` output for QC pipelines,
  one row per library.
- **No header preservation guarantee:** `@PG` records are auto-chained
  via `PP:` (samblaster does not chain).
- **Idempotent dup flag:** dupblaster *overwrites* `FLAG_DUPLICATE` on
  every output record with the current run's decision (matching Picard
  and samtools markdup); samblaster *ORs* it in, preserving prior
  markings. Overwrite is idempotent across re-runs (`bwa mem | dupblaster`
  and `bwa mem | dupblaster | dupblaster` produce identical output); OR
  is not.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, the test
suite, and release flow.

## License

MIT. See [LICENSE](LICENSE).

## Citing dupblaster

Every release is archived on Zenodo. Cite the concept DOI [10.5281/zenodo.21445780](https://doi.org/10.5281/zenodo.21445780), which always resolves to the latest release, or the version-specific DOI from that record if you need to pin the exact version you ran (`dupblaster --version` reports it). Please also cite the samblaster paper, whose algorithm dupblaster adapts — see [§ Acknowledgements](#acknowledgements).

## Acknowledgements

dupblaster is inspired by, and adapts the coordinate-based duplicate
detection algorithm of, [samblaster][samblaster] by Greg Faust and Ira
Hall. If you use dupblaster in published work, please cite the
samblaster paper:

> Faust, G.G. and Hall, I.M., *SAMBLASTER: fast duplicate marking and
> structural variant read extraction*, **Bioinformatics**
> 30(17): 2503-2505 (2014).
> [doi:10.1093/bioinformatics/btu314][faust-hall]
