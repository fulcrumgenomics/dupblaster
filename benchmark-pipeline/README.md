# dupblaster benchmark pipeline

A Snakemake + pixi pipeline that downloads an NYGC 1000G high-coverage
CRAM (bwa-mem aligned, GRCh38, ~30×), subsamples it to ~8× effective
coverage, strips any prior duplicate flags, name-sorts it, then times
ten duplicate-marking tools. Comparison is done by the Rust
`bench-compare` tool (a workspace sibling of dupblaster) using Picard's
`kf` fragment-key tag as the oracle key — no template-coordinate sort is
needed.

Designed to run on a single host with a large local disk. Tested on
Amazon Linux 2023 (x86_64 and aarch64) and macOS arm64.

## Why NYGC instead of the original phase-3 BAMs

The 1000G phase-3 release BAMs (`.mapped.bam`) were aligned with
**bwa-aln**, which does not emit supplementary alignments — there are zero
records with `FLAG 0x800` set, so the supplementary inheritance code path
is never exercised. The 2019 NYGC re-alignments use **bwa-mem**, which
produces supplementary alignments for split / chimeric reads, and are more
representative of modern aligned data.

## Why we don't use Picard RevertSam

The source CRAM ships with Picard `MarkDuplicates` already run — every
record carries either a set `FLAG_DUPLICATE` (0x400) bit or not,
reflecting the production centre's call at sequencing time. Without
cleanup, tools that OR their decision into the existing flag (rather than
overwriting it) would conflate historical and current markings.

The obvious cleanup is `Picard RevertSam`, but it has a fatal flaw for
this benchmark: by design `RevertSam` **drops every secondary and
supplementary record**, which defeats the reason we picked the NYGC
bwa-mem alignments in the first place. Instead the prep pipeline does:

```
samtools view --remove-flags 0x400   # strip prior dup bits on all records
  | samtools sort -n                 # lexicographic name sort
```

This preserves all records (primary, secondary, and supplementary) and
produces a fresh name-sorted BAM that every tool can consume on equal
footing.

## Tool matrix

Ten tools are benchmarked (`TOOLS` list in the Snakefile):

| Tool                       | Mode / flag                           | Input prep               |
|----------------------------|---------------------------------------|--------------------------|
| `samblaster`               | default                               | name-sorted SAM text     |
| `dupblaster-sam`           | default                               | name-sorted SAM text     |
| `dupblaster-bam`           | default                               | name-sorted BAM          |
| `dupblaster-picard-approx` | `--single-end-strategy picard-approx` | name-sorted BAM          |
| `dupblaster-picard-exact`  | `--single-end-strategy picard-exact`  | name-sorted BAM          |
| `picard`                   | `MarkDuplicates ASSUME_SORT_ORDER=queryname` | name-sorted BAM  |
| `samtools-markdup`         | default (no `-S`)                     | coord-sorted BAM + fixmate |
| `samtools-markdup-S`       | `-S` (propagate to secondary/supp)    | coord-sorted BAM + fixmate |
| `samtools-markdup-seq-S`   | `-m s -S` (sequence mode + propagate) | coord-sorted BAM + fixmate |
| `dupsifter`                | `-W` (WGS-only mode)                  | name-sorted BAM          |

The three dupblaster modes (`-sam`, `-bam`, `-picard-approx`) run the
same binary; only the input format and `--single-end-strategy` differ.
`dupblaster-picard-exact` buffers orphan / single-end reads to a temp
BAM and re-marks them in a second pass after the pair pass — its timed
window includes both passes and the temp-BAM write, and its output is
name-sorted by the untimed conversion step.

`samtools-markdup` requires `fixmate -m | sort` upstream (for `MC` / `ms`
mate tags and coordinate order); that prep lives outside the timed window,
analogous to how a coord-sorted input is prepped before other
coord-order tools.

`samtools-markdup-seq-S` is the same as `samtools-markdup-S` but adds
`-m s` (sequence mode). markdup's default `-m t` (template mode) folds
R1/R2 (first/second-in-template) identity into the pair key for same-strand
pairs, which diverges from Picard and dupblaster (both coordinate-canonical
and R1/R2-agnostic) on FF/RR and inter-chromosomal geometries. Sequence
mode drops that distinction and reaches 100% paired-set concordance with
Picard (vs 99.87% for `-m t`), at no meaningful runtime cost — so it is the
apples-to-apples samtools mode for this comparison.

`dupsifter` is the huishenlab samblaster-lineage tool, run in `-W`
(WGS-only) mode because the bench data is unconverted WGS.

## Timing design

Each timed rule:

1. Acquires the whole `bench` Snakemake resource pool (size 100), so no
   two timed rules ever run concurrently — wall-time isn't contaminated by
   another rule sharing CPU.
2. Runs `env time -v` over the tool only, writing the tool's native,
   **uncompressed** output (L0 BGZF or SAM text).
3. After `env time -v` exits (outside the timed window), converts the raw
   output to the shared comparison artifact: a name-sorted L1 BGZF BAM
   (`out.nsort.bam`). Order-preserving tools recompress with
   `samtools view -b --output-fmt-option level=1`; coordinate-order tools
   (`samtools-markdup`, `samtools-markdup-S`, `samtools-markdup-seq-S`) and
   `dupblaster-picard-exact`
   (whose orphans are appended at the end) are name-sorted with
   `samtools sort -n -l 1`.

Keeping recompression out of the timed window means output-format cost
never penalizes a tool, while every comparison still gets uniformly
compressed, uniformly ordered input.

## Comparison design

`bench-compare` co-streams Picard's `kf`-tagged `out.nsort.bam` against
each partner tool's `out.nsort.bam` in one pass. Picard runs with
`TAG_DUPLICATE_KEY=true` (a patched Picard build — see the Picard jar
note below), which emits a canonical `kf` fragment-key tag on every
primary mapped record. `bench-compare` reads that tag directly off the
raw BAM bytes without a SAM-text round-trip.

From one pass the tool produces:

- **set-equivalence concordance** (`setcmp.tsv`): for each Picard
  duplicate set (keyed by the sorted combination of the ends' `kf` tags),
  how many templates did Picard mark vs how many did the partner mark?
  Reported by template category (`pe_both_mapped`, `orphan`, `single_end`,
  …). Sets of size 1 (Picard singletons) are tallied separately.
- **orphan-discordance triage** (`orphan_triage.tsv`): four-bucket
  classification of per-orphan disagreements —
  `cat1_cross_table` (orphan's position coincides with a pair end),
  `cat2_tiebreaker` (tools picked different representatives at the same
  position with no cross-table reason), `cat3_picard_only_other`,
  `cat4_partner_only_other`.
- **dup-flag inheritance** (`inheritance.tsv`): per tool, the fraction of
  QNAME groups that have inconsistent dup flags across primary + secondary +
  supplementary records.

A per-comparison human-readable summary is written to
`results/runs/<sample>/<partner>/rep1/compare.txt`.

The comparison is done only on rep 1 (dup-marking is deterministic). All
9 tools except `picard` are compared against Picard.

> **Picard jar (bring your own):** `PICARD_JAR` defaults to
> `benchmark-pipeline/picard.jar` (gitignored — drop a jar there), or override
> with `--config picard_jar=/path/to/picard.jar`. Two cases:
> - **Timing benchmark** (`./run.sh results/bench.tsv`) runs Picard *vanilla*,
>   so a **stock** Picard release works.
> - **kf-concordance comparison** (`run_picard` + `bench-compare`) needs a build
>   with the opt-in `TAG_DUPLICATE_KEY=true` argument that emits the `kf`
>   fragment-key tag `bench-compare` keys on. That lives on the
>   [`tf_tag_duplicate_key`](https://github.com/broadinstitute/picard) Picard
>   branch and is not yet upstreamed; you must build it yourself for the
>   concordance numbers (the timing numbers don't need it).

## Quick start

```bash
cd benchmark-pipeline
./install.sh
./run.sh --dry-run     # preview the DAG
./run.sh               # full run; first invocation downloads ~17 GB
```

`install.sh` fetches `pixi` if it isn't on PATH, materializes the pixi
environment (snakemake, samtools, samblaster, dupsifter, openjdk, GNU
time, aria2c), and builds dupblaster + bench-compare from the parent
crate in release mode.

Override knobs at the CLI:

```bash
./run.sh --config replicates=3          # default is 1; bump to measure variance
./run.sh --config samples=HG03953       # run only this sample
./run.sh --config dupblaster_bin=/abs/path/to/dupblaster
./run.sh --config bench_compare_bin=/abs/path/to/bench-compare
./run.sh --config picard_jar=/abs/path/to/picard.jar
./run.sh --config picard_heap=12g       # default is 8g
./run.sh --cores 8
```

Anything after `--` is forwarded verbatim to Snakemake:

```bash
./run.sh -- --rerun-triggers mtime input
```

## Requirements

- **pixi** — `install.sh` fetches it automatically if missing.
- **cargo** — needed to build dupblaster and bench-compare from the parent
  crate. Install via [rustup](https://rustup.rs) if missing; or skip the
  build step with `./install.sh --skip-build` and supply pre-built
  binaries via `--config dupblaster_bin=…`.
- **Patched Picard jar** — see the Picard jar note above.
- **Disk** — ~200 GB free per sample at peak working set. For HG03953:
  reference (~3 GB) + CRAM (~14 GB, persistent) + subsampled BAM (~2 GB,
  temp) + name-sorted BAM (~2 GB, persistent) + per-tool raw + nsort outputs
  + temp sort scratch. Local NVMe materially improves wall-time.
- **Memory** — 16 GiB minimum, 24 GiB recommended. Picard runs with an
  8 GiB heap (`picard_heap=8g` by default); the JVM adds ~1–2 GiB
  overhead. The bench-pool resource model prevents timed rules from running
  concurrently, but 16 GiB leaves limited headroom for OS page cache.
- **Bandwidth** — first run downloads ~17 GB (CRAM + reference) from
  AWS S3 us-east-1 (materially faster than EBI FTP for these files).

> **Reproducibility note: `pixi.lock` is authoritative.**
> The committed `pixi.lock` resolves every dependency to an exact hash.
> `pixi install` (called by `install.sh`) honors the lockfile and reproduces
> the published environment byte-for-byte. Do not run `pixi update` — that
> re-resolves the dependency graph and will drift you off the published
> numbers.

## Samples

The sample dictionary lives in the `Snakefile`. The default single sample
is **HG03953** (NYGC re-alignment, STU population, ~29.8× CRAM,
`ERR3242904`). The subsample fraction is `0.27`, producing ≈ 8× effective
coverage (0.27 of ~30×). The seed is fixed at `42` for reproducibility.

Add more samples by appending entries to `SAMPLES` in the Snakefile:

```python
"NAMEXX": {
    "population":     "STU",
    "run_id":         "ERR3242904",
    "coverage_x":     29.8,
    "cram_md5":       "...",        # from NYGC sequence.index col MD5
    "subsample_frac": 0.27,         # fraction of templates to keep
    "subsample_seed": 42,           # samtools view --subsample-seed
},
```

To run only a subset of configured samples:

```bash
./run.sh --config samples=HG03953
```

The GRCh38 + decoy + HLA reference is shared across all samples and cached
after the first run.

## Job graph

```
download_reference ──> index_reference ─┐
                                         │
download_cram(N) ───────────────────────> prep_subsampled_bam(N)
                                                  │
                                           prep_namesorted_bam(N)  [persistent]
                                          /        |         \
               prep_namesorted_sam(N) <--/         |          \--> prep_fixmate_coord_bam(N)
                       |                           |                        |
               run_samblaster           run_dupblaster_bam        run_samtools_markdup
               run_dupblaster_sam       run_dupblaster_picard_approx
                                        run_dupblaster_picard_exact
                                        run_picard
                                        run_dupsifter

every run_* emits out.nsort.bam + time.txt
    │
    ├── collate_bench  ────────────────> results/bench.tsv
    └── compare  ──────────────────────> results/setcmp.tsv
                                         results/orphan_triage.tsv
                                         results/inheritance.tsv
```

Prep intermediates marked `temp()` in the Snakefile are deleted by
Snakemake once the last downstream consumer completes.
`results/prepped/<sample>.qnsort.bam` is **not** `temp()` — it is the
canonical cached input for the whole bench; persisting it means re-runs
after partial failures never re-download or re-subsample the (large) CRAM.

## Outputs

Persistent outputs from a complete run:

```
results/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}
results/download/<sample>.final.cram
results/prepped/<sample>.qnsort.bam

results/runs/<sample>/<tool>/rep<N>/time.txt        # GNU time -v output
results/runs/<sample>/<tool>/rep<N>/out.nsort.bam   # name-sorted L1 comparison artifact
results/runs/<sample>/<tool>/rep<N>/compare.txt     # bench-compare human-readable summary
results/runs/<sample>/picard/rep<N>/metrics.txt     # Picard DuplicationMetrics
results/runs/<sample>/dupsifter/rep<N>/dupsifter_stats.txt

results/bench.tsv          # wall_s, user_s, sys_s, cpu_percent, max_rss_kb, exit_status
results/setcmp.tsv         # set-equivalence concordance vs Picard (kf-keyed)
results/orphan_triage.tsv  # four-bucket orphan-discordance triage
results/inheritance.tsv    # dup-flag consistency across primary+secondary+supplementary
```

`bench.tsv` columns: `sample`, `tool`, `rep`, `wall_s`, `user_s`, `sys_s`,
`cpu_percent`, `max_rss_kb`, `exit_status`.

`setcmp.tsv` columns: `sample`, `partner`, `pe_sets`, `pe_concordant`,
`pe_concord_pct`, `orphan_sets`, `orphan_concordant`, `orphan_concord_pct`,
`picard_singletons_partner_marked`.

`orphan_triage.tsv` columns: `sample`, `partner`, `total_orphans`,
`picard_marked`, `partner_marked`, `concordant_dup`, `concordant_nondup`,
`discordant`, `cat1_cross_table`, `cat2_tiebreaker`, `cat3_picard_only_other`,
`cat4_partner_only_other`.

`inheritance.tsv` columns: `sample`, `tool`, `total_qnames`,
`groups_with_supp`, `consistent`, `inconsistent`, `inconsistency_pct`.

## Adding a new tool

1. Add the tool name to the `TOOLS` list near the top of the Snakefile.
2. Add a `run_<tool-name>` rule that:
   - Takes the appropriate prepped input (`prep_namesorted_bam`,
     `prep_namesorted_sam`, or `prep_fixmate_coord_bam`).
   - Declares `resources: bench = 100` to serialize it with the other
     timed rules.
   - Runs `env time -v -o {output.time_txt} <tool ...>` writing native
     uncompressed output.
   - Converts the raw output to `out.nsort.bam` (name-sorted L1 BAM) in
     the same rule, outside the `env time -v` block.
   - `run_dupblaster_bam` is a good template for a single-binary tool;
     `run_samtools_markdup` shows how to handle a coordinate-order tool
     that needs a post-hoc name sort.
3. If the tool needs an input format that none of the existing prep rules
   produces, add a `prep_<format>` rule with `temp()` output.

`collate_bench` and `compare` derive their tool lists from `TOOL_NAMES`
and `[t for t in TOOL_NAMES if t != "picard"]` respectively, so the new
tool is picked up automatically.

## Adding a new sample

Append an entry to the `SAMPLES` dict in the Snakefile (schema above).
Source `run_id` and `cram_md5` from the NYGC 1000G high-coverage sequence
index:

```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
```

The reference download is shared and already cached after the first run.

## Caveats

- **Hardware sensitivity.** Wall-clock numbers depend heavily on the host:
  storage (local NVMe vs workstation SSD — slower IO widens the spread
  between BAM-native and SAM tools) and CPU (Picard's JVM shows 1.5–2× more
  variance on older cores). Treat the top-level README's Benchmarks numbers
  as representative, not absolute.
- **Subsampled coverage.** The default 0.27 subsample fraction on a ~30×
  CRAM produces ≈ 8× effective coverage. Concordance numbers are stable
  across subsample fractions; runtime and RSS scale with template count.
- **Comparison uses rep 1 only.** Dup-marking is deterministic; running
  multiple replicates (`--config replicates=N`) measures runtime variance,
  not marking variance.
- **Picard and samblaster both pass `--ignore-unmated` / `VALIDATION_STRINGENCY=LENIENT`.**
  The NYGC CRAMs contain occasional genuine orphan reads (mates that were
  filtered downstream) that survived the prep because we don't use
  RevertSam. The flag turns hard errors into counters; `0 unmated` in
  samblaster's stats means the data was clean.
- **macOS clamps GNU `time` precision to seconds** on wall time, but the
  timer subsystem value (reported by the `time` package from conda-forge)
  should be sub-second accurate.
