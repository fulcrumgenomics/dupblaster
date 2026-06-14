//! dupblaster — the command-line application.
//!
//! This binary is the whole tool: CLI parsing ([`Args`]), end-to-end run
//! orchestration ([`run`]), and the end-of-run resource-usage footer. The
//! reusable engine and IO live in `dupblaster_lib` (`dedup`, `sig`, the
//! readers/writers, `metrics`, …) — this file wires them together and talks to
//! the user.
//!
//! Long-form flags follow GNU style (`--kebab-case`) — diverging from C++
//! samblaster's `--camelCase` convention. Short flags (`-i`, `-o`, `-r`, `-q`)
//! match upstream. The SV-extraction flags from upstream (`-d`, `-s`, `-u`,
//! `-a`, `-e`, `-M`, `--maxSplitCount`, …) have been intentionally dropped —
//! see README for the rationale.

mod cigar;
mod dedup;
mod io_threading;
mod metrics;
mod raw_reader;
mod raw_writer;
mod sam_reader;
mod sig;

use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use bgzf::CompressionLevel;
use clap::{Parser, ValueEnum};
use fgumi_raw_bam::RawRecord;
use noodles_sam::Header;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::Program;
use noodles_sam::header::record::value::map::program::tag as program_tag;

use crate::dedup::{LibraryIndex, ProcessorOptions, RecordProcessor, Stats};
use crate::io_threading::ThreadedReader;
use crate::metrics::{Metrics, write_rows_to_path};
use crate::raw_reader::RawBamReader;
use crate::raw_writer::RawBamWriter;
use crate::sam_reader::SamReader;
use crate::sig::{MethylationMode, SingleEndStrategy};

/// Crate-level build identifier shown in `--version` and the `@PG VN:` tag.
const DUPBLASTER_BUILD: &str = env!("CARGO_PKG_VERSION");

/// Global allocator — mimalloc outperforms the system allocator on the
/// workload of many small, short-lived hash-set entries that dominate
/// dupblaster's hot path.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

// ── Command-line interface ──────────────────────────────────────────────────

/// CLI value-enum mirror of [`SingleEndStrategy`]. The mirror exists so that
/// the kebab-case CLI spellings and per-variant help text live next to the rest
/// of the CLI definition, while the algorithmic enum stays in the library near
/// the code that consumes it.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum SingleEndStrategyCli {
    /// (Default) Strand-aware 5'-aligned coordinate key. A forward
    /// orphan and a reverse orphan at the same 5' position are NOT
    /// duplicates of each other. Picard's `fragSort` uses an equivalent
    /// key.
    #[value(name = "strand-aware")]
    StrandAware,
    /// Strand-aware keying plus a Picard-style cross-check: each end of
    /// a fully-mapped pair is also registered in a fragment-level
    /// table, so subsequent orphans / single-end reads at those
    /// positions are marked as duplicates of the pair. Approximate
    /// (order-sensitive in a streaming pass — an orphan that arrives
    /// before its corresponding pair will pass through as non-dup).
    /// Uses ~2x the dup-table memory of `strand-aware`.
    #[value(name = "picard-approx")]
    PicardApprox,
    /// Exact Picard fragment semantics. Pairs stream straight through;
    /// mapped-orphan / single-end reads are buffered to a temporary
    /// uncompressed BAM, then re-processed after the pair pass against a
    /// fragment table holding every paired read end — so "fragments never
    /// beat pairs" holds exactly, independent of input order (unlike
    /// `picard-approx`). Costs a temporary on-disk copy of the fragment
    /// blocks and emits those fragments at the *end* of the output stream
    /// (not in input order). Intended for paired data, where orphans are a
    /// small fraction. See `--tmp-dir`.
    #[value(name = "picard-exact")]
    PicardExact,
    /// samblaster v0.1.23+ legacy behavior: leftmost-aligned coordinate
    /// with the strand bit dropped, so a forward orphan and a reverse
    /// orphan whose alignments share a leftmost position collide. NOT
    /// recommended for short-read PE data — see the README's "single-end
    /// strategies" section for the full discussion. Provided only for
    /// byte-compatibility with samblaster output on long-read singleton
    /// workflows.
    #[value(name = "samblaster-legacy")]
    SamblasterLegacy,
}

impl SingleEndStrategyCli {
    /// Map to the algorithm-level [`SingleEndStrategy`] consumed by the engine.
    /// (An inherent method rather than a `From` impl: the orphan rule forbids
    /// `impl From<LocalCliEnum> for ForeignLibEnum` in the binary crate.)
    fn to_strategy(self) -> SingleEndStrategy {
        match self {
            Self::StrandAware => SingleEndStrategy::StrandAware,
            Self::PicardApprox => SingleEndStrategy::PicardApprox,
            Self::PicardExact => SingleEndStrategy::PicardExact,
            Self::SamblasterLegacy => SingleEndStrategy::SamblasterLegacy,
        }
    }
}

/// CLI value-enum mirror of [`MethylationMode`]. The methylation field is an
/// `Option` on [`Args`] (`None` = off = standard WGS keying), so this enum only
/// names the *modes*. Only `directional` exists today; a non-directional / PBAT
/// variant is intentionally deferred (see [`MethylationMode`]).
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum MethylationModeCli {
    /// Directional bisulfite / enzymatic-conversion prep (WGBS, EM-seq, TAPS).
    /// Keys each pair in template (first-of-pair → second-of-pair) order
    /// instead of coordinate-canonical order, so the two original strands
    /// (OT/OB) of a fragment at one locus are kept distinct while same-strand
    /// PCR copies still collapse. Correct for any prep that ligates adapters
    /// before conversion. Does NOT handle non-directional / PBAT libraries.
    #[value(name = "directional")]
    Directional,
}

impl MethylationModeCli {
    /// Map to the algorithm-level [`MethylationMode`] consumed by the engine.
    fn to_mode(self) -> MethylationMode {
        match self {
            Self::Directional => MethylationMode::Directional,
        }
    }
}

/// Parse a `--compression-level` argument by delegating validation to
/// `bgzf::CompressionLevel::new`, so the bgzf crate stays the single
/// source of truth for the valid range (0..=12 in practice; 0 = stored,
/// 1-9 standard zlib, 10-12 libdeflate's extra-strong tiers).
fn parse_compression_level(s: &str) -> Result<CompressionLevel, String> {
    let n: u8 = s.parse().map_err(|e| format!("not a u8: {e}"))?;
    CompressionLevel::new(n).map_err(|e| format!("{e}"))
}

/// Short one-line `about` text shown by `--help`.
const SHORT_ABOUT: &str = "Mark or remove duplicate reads in a query-grouped SAM/BAM file.";
/// Full `about` text shown by `--help --help` (long form).
const LONG_ABOUT: &str = "Mark or remove duplicates in a query-grouped SAM/BAM file.\n\
Input must be query-grouped (all alignments for the same QNAME adjacent, \
typically straight from the aligner); output is always BAM (uncompressed by \
default — see --compression-level).";

/// Parsed command-line arguments — see `--help` for per-flag descriptions.
#[derive(Parser, Debug, Clone)]
#[command(name = "dupblaster", disable_version_flag = true, about = SHORT_ABOUT, long_about = LONG_ABOUT)]
pub struct Args {
    /// Input SAM/BAM file [stdin].
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output BAM file [stdout]. The path must end in `.bam` — output is
    /// always BAM, never SAM. Use `-` for stdout (no extension check).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// BGZF compression level for the output BAM (0-9 typical, up to 12
    /// for libdeflate's strongest tier). Level 0 (the default) produces
    /// "stored" BGZF blocks — same as `samtools view -u` and what most
    /// dupblaster pipelines want, since the downstream sort recompresses.
    /// Use a non-zero value when piping to a non-recompressing sink.
    #[arg(long = "compression-level",
          default_value = "0",
          value_parser = parse_compression_level)]
    pub compression_level: CompressionLevel,

    /// Remove duplicate reads from the output instead of just flagging them.
    #[arg(short = 'r', long = "remove-dups")]
    pub remove_dups: bool,

    /// Add MC (mate CIGAR) and MQ (mate MAPQ) tags to all paired SAM output.
    #[arg(long = "add-mate-tags")]
    pub add_mate_tags: bool,

    /// Suppress abort on unmated alignments.
    #[arg(long = "ignore-unmated")]
    pub ignore_unmated: bool,

    /// Maximum expected read length (in bp). Affects only the synthetic-
    /// genome padding used by the duplicate-detection index and a "reads
    /// longer than this" warning counter; no algorithmic effect at typical
    /// short-read sizes. Bump if running on long reads. The upper bound
    /// keeps `2 * max_read_length` (the per-contig padding) well within
    /// `i32`, so the super-contig length arithmetic can't overflow.
    #[arg(long = "max-read-length",
          default_value_t = 1000,
          value_parser = clap::value_parser!(i32).range(1..=10_000_000))]
    pub max_read_length: i32,

    /// Strategy for keying single-end / orphan reads in the dedup hash
    /// table. The default (`strand-aware`) matches Picard's `fragSort`
    /// behavior at the single-end-key level. See the README's "single-end
    /// strategies" section for a discussion of the options.
    #[arg(long = "single-end-strategy",
          value_enum,
          default_value_t = SingleEndStrategyCli::StrandAware)]
    pub single_end_strategy: SingleEndStrategyCli,

    /// Disable library-aware duplicate marking. By default, when the header
    /// declares more than one library (distinct `@RG LB:` values), duplicates
    /// are only called *within* a library — matching Picard MarkDuplicates, and
    /// the `--stats` TSV reports one row per library. This flag forces a single
    /// combined dedup table across all reads (samblaster's library-agnostic
    /// behavior). No effect when the header has ≤1 library, where the two modes
    /// are identical.
    #[arg(long = "library-unaware")]
    pub library_unaware: bool,

    /// Methylation-aware duplicate marking for bisulfite / enzymatic-conversion
    /// data. Omitted by default (standard WGS keying, which stays
    /// Picard-exact-capable). With `directional`, pairs are keyed in template
    /// order so the two original strands (OT/OB) of a fragment at the same
    /// locus are kept distinct — the correct behavior for WGBS / EM-seq / TAPS,
    /// where collapsing them would discard independent methylation information.
    /// Non-directional / PBAT libraries are not supported. See the README's
    /// "methylation mode" section.
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<MethylationModeCli>,

    /// Directory for the temporary uncompressed BAM that `--single-end-strategy
    /// picard-exact` uses to buffer orphan / single-end reads between its two
    /// passes. Defaults to the system temp dir (`$TMPDIR`). The file is
    /// deleted when dupblaster exits. No effect under other strategies.
    #[arg(long = "tmp-dir")]
    pub tmp_dir: Option<PathBuf>,

    /// Write a per-library TSV summary of run metrics to PATH (one row per
    /// library; a `.gz`/`.bgz` suffix gzip-compresses it). Schema: see the
    /// README. Columns include sample, library, template/dup counts,
    /// `frac_duplicates`, and a Picard-style `estimated_library_size`.
    #[arg(long = "stats")]
    pub stats: Option<PathBuf>,

    /// Sample name written to the `sample` column of `--stats` output. If
    /// omitted, dupblaster comma-joins the unique `@RG SM:` values from the
    /// input header (empty if none are present).
    #[arg(long = "sample")]
    pub sample: Option<String>,

    /// Output fewer statistics.
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Minimum bins per side for the partitioned dedup hash table.
    /// Controls cell count via `cells ≈ (bins+1)² × 4`. Lower values
    /// give fewer/larger cells (lower steady-state RSS at high
    /// coverage, but bigger resize-peak spikes — single hashbrown
    /// resizes can momentarily double RSS); higher values give
    /// more/smaller cells (resize peaks stay tiny; slight overhead
    /// at low coverage).
    ///
    /// Default 32 was picked from a sweep of values 1..=128 — it sits
    /// just past the "U-curve" memory minimum at low coverage while
    /// keeping resize peaks small at all coverages.
    ///
    /// Capped at 8192: `bin_count` tracks ~2×`min_bins`, and the pair
    /// table's flat index is `s1 * stride + s2` with `stride =
    /// (bin_count+1)*2`, computed in `u32`. The cap keeps `stride²` (the
    /// largest index) safely below `u32::MAX`, so the index can't wrap.
    #[arg(long = "min-bins",
          default_value_t = 32,
          hide = true,
          value_parser = clap::value_parser!(u32).range(1..=8192))]
    pub min_bins: u32,

    /// Verify BGZF CRC32 on input. Default: on for file input, off for stdin
    /// (where the producer is assumed trusted, e.g. piped from bwa-mem).
    /// Mutually exclusive with `--no-check-crc`.
    #[arg(long = "check-crc", conflicts_with = "no_check_crc")]
    pub check_crc: bool,

    /// Skip BGZF CRC32 verification on input regardless of source.
    /// Mutually exclusive with `--check-crc`.
    #[arg(long = "no-check-crc", conflicts_with = "check_crc")]
    pub no_check_crc: bool,

    /// Size (MB) of the user-space ring buffer between the input IO thread
    /// and the worker. Big enough to absorb upstream bursts (bwa-mem dumps
    /// ~250 MB of reads in tight bursts every few seconds). Default 16 MB
    /// is plenty when the worker drain rate exceeds the producer's burst
    /// rate — increase if your producer is unusually fast or bursty.
    #[arg(long = "read-buffer-mb", default_value_t = 16, value_parser = clap::value_parser!(u32).range(1..=4096))]
    pub read_buffer_mb: u32,

    /// Size (MB) of the user-space ring buffer between the worker and the
    /// output IO thread. Larger than the read buffer by default because
    /// downstream sorters (samtools sort, mako) periodically pause input
    /// for 1-3 s while flushing a sort chunk to disk. Default 64 MB
    /// absorbs ~2-5 s of downstream blocking at typical bwa-mem-limited
    /// output rates; bump to 256+ for slower downstreams.
    #[arg(long = "write-buffer-mb", default_value_t = 64, value_parser = clap::value_parser!(u32).range(1..=4096))]
    pub write_buffer_mb: u32,

    /// Print the dupblaster version to stdout and exit.
    #[arg(short = 'V', long = "version")]
    pub show_version: bool,
}

impl Args {
    /// Resolved CRC-verify setting. Explicit `--check-crc` / `--no-check-crc`
    /// override the default; otherwise it's on for file input (data may have
    /// been transferred/archived) and off for stdin (the producer is assumed
    /// trusted, e.g. piped fresh from bwa-mem).
    pub fn effective_check_crc(&self) -> bool {
        if self.check_crc {
            return true;
        }
        if self.no_check_crc {
            return false;
        }
        matches!(self.input.as_deref(), Some(p) if p.to_string_lossy() != "-")
    }

    /// Validate caller-facing invariants that clap can't express directly.
    /// Currently: `-o` must be `-` (stdout) or end in `.bam`, since output
    /// is always BAM and a non-`.bam` extension would silently produce a
    /// confusingly-named file.
    pub fn validate(&self) -> Result<()> {
        if let Some(p) = &self.output {
            let s = p.to_string_lossy();
            if s != "-" && !s.ends_with(".bam") {
                bail!(
                    "output path {} must end in `.bam` (dupblaster only writes BAM); \
                     use `-` to send BAM to stdout",
                    p.display()
                );
            }
        }
        Ok(())
    }
}

// ── Binary entry point ──────────────────────────────────────────────────────

fn main() -> ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    let args = Args::parse();
    match run(args) {
        // `run` chooses the exit code: SUCCESS normally, FAILURE when the run
        // completed but coordinate clamping made the result suspect (the
        // warning is printed inside `run`). A returned `Err` is a genuine
        // mid-run failure.
        Ok(code) => code,
        Err(err) => {
            eprintln!("dupblaster: {err:#}");
            eprintln!("dupblaster: Premature exit (return code 1).");
            ExitCode::FAILURE
        }
    }
}

// ── Run orchestration ───────────────────────────────────────────────────────

/// Run dupblaster end to end: detect the input format, wire reader → processor
/// → writer, drive the QNAME-block loop, and emit the summary, metrics, and
/// resource footer. Returns the process exit code.
fn run(args: Args) -> Result<ExitCode> {
    if args.show_version {
        // Version goes to stdout (not stderr) so `dupblaster --version | head`
        // and tooling that captures it work as expected.
        println!("dupblaster {DUPBLASTER_BUILD}");
        return Ok(ExitCode::SUCCESS);
    }
    args.validate()?;
    let started = StartedRun::now();
    if !args.quiet {
        eprintln!("dupblaster: Version {DUPBLASTER_BUILD}");
    }

    // Open input. The first byte tells us SAM vs BAM: 0x1f is the BGZF
    // gzip magic; '@' (0x40) is a SAM header line. Anything else is an
    // error.
    //
    // Input goes through a dedicated IO read thread + ring buffer
    // (see [`ThreadedReader`]) so the worker never blocks on the kernel pipe.
    let raw_source: Box<dyn std::io::Read + Send> = match args.input.as_deref() {
        Some(p) if p.to_string_lossy() != "-" => {
            let f = File::open(p).with_context(|| format!("opening {} for read", p.display()))?;
            Box::new(f)
        }
        _ => Box::new(std::io::stdin()),
    };
    let read_buf_bytes = (args.read_buffer_mb as usize).saturating_mul(1024 * 1024);
    let mut reader_box: Box<dyn BufRead> =
        Box::new(ThreadedReader::new(raw_source, read_buf_bytes));
    let input_name =
        args.input.as_deref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdin".into());
    let input_format = detect_format(&mut reader_box)?;
    let check_crc = args.effective_check_crc();
    let mut reader: Reader = match input_format {
        Format::Bam => {
            if !args.quiet {
                eprintln!(
                    "dupblaster: Reading BAM from {input_name} (CRC verify: {}).",
                    if check_crc { "on" } else { "off" }
                );
            }
            Reader::Bam(RawBamReader::new(reader_box, check_crc))
        }
        Format::Sam => {
            if !args.quiet {
                eprintln!("dupblaster: Reading SAM from {input_name}.");
            }
            Reader::Sam(SamReader::new(reader_box))
        }
    };

    let mut header = reader.read_header()?;
    if header.reference_sequences().is_empty() {
        bail!("Input has no @SQ reference sequences. Exiting.");
    }
    if !args.quiet {
        eprintln!(
            "dupblaster: Loaded {} header sequence entries.",
            header.reference_sequences().len()
        );
    }

    let ref_lengths: Vec<i32> = header
        .reference_sequences()
        .iter()
        .map(|(name, m)| {
            // BAM stores contig lengths as i32; bail on overflow rather than
            // silently wrapping a >2 Gb contig into a negative value.
            i32::try_from(usize::from(m.length())).map_err(|_| {
                anyhow::anyhow!(
                    "Contig {} length {} exceeds i32::MAX",
                    String::from_utf8_lossy(name),
                    usize::from(m.length()),
                )
            })
        })
        .collect::<Result<Vec<i32>>>()?;

    append_dupblaster_pg(&mut header)?;

    let write_buf_bytes = (args.write_buffer_mb as usize).saturating_mul(1024 * 1024);
    let mut out = RawBamWriter::open(
        args.output.as_deref(),
        &header,
        write_buf_bytes,
        args.compression_level,
    )?;

    let opts = processor_options(&args);
    // Library-aware duplicate marking: when the header declares >1 distinct
    // `@RG LB:` value (and `--library-unaware` wasn't passed), dedup state is
    // partitioned per library so duplicates are only called within a library.
    let library_index = LibraryIndex::from_header(&header, args.library_unaware);
    let num_libs = library_index.num_libs();
    let mut stats = Stats::new(&library_index);
    let mut processor =
        RecordProcessor::from_ref_lengths(&ref_lengths, opts, args.min_bins, library_index);
    if !args.quiet {
        eprintln!(
            "dupblaster: bin_shift={} bin_count={} (min-bins={}), partition cells ≈ {}",
            processor.bin_shift(),
            processor.bin_count(),
            args.min_bins,
            ((processor.bin_count() + 1) * 2).pow(2),
        );
        if num_libs > 1 {
            eprintln!(
                "dupblaster: library-aware mode — {num_libs} library buckets (incl. 'Unknown \
                 Library'); duplicates are called within a library and reported per library."
            );
        }
        if args.methylation_mode.is_some() {
            eprintln!(
                "dupblaster: methylation mode = directional — pairs keyed in template order so \
                 opposite-strand (OT/OB) fragments at one locus are kept distinct."
            );
        }
    }

    // dupblaster requires query-grouped input: records are read in maximal
    // runs sharing a QNAME ("blocks") and each block is processed as a unit.
    // `for_each_block` owns that grouping loop. picard-exact takes a separate
    // two-pass driver (pairs stream out, fragments are buffered then
    // re-processed); every other strategy is a single streaming pass.
    if args.single_end_strategy.to_strategy() == SingleEndStrategy::PicardExact {
        run_picard_exact(&mut reader, &mut processor, &mut stats, &mut out, &header, &args)?;
    } else {
        let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
        for_each_block(&mut reader, &mut pool, |block| {
            processor.process_block(block, &mut stats, &mut out)
        })
        .context("processing record block")?;
    }

    print_run_stats(&stats, &args);

    if let Some(stats_path) = args.stats.as_deref() {
        let rows = Metrics::rows_from_stats(&stats, &header, args.sample.as_deref());
        write_rows_to_path(&rows, stats_path).context("writing --stats TSV")?;
    }

    out.finish().context("finishing main output")?;

    // Footer (wall/CPU/RSS) is emitted *after* output finalization so its
    // numbers include any flush/close time.
    report(&started, stats.totals().id_count, args.quiet);

    // Output is fully written; coordinate clamping at contig edges is a
    // correctness compromise (potential false duplicates among the clamped
    // reads), so we exit non-zero to flag it even though the run completed.
    // The actionable warning was printed by `print_run_stats`.
    if stats.clamped_template_count > 0 {
        return Ok(ExitCode::FAILURE);
    }
    Ok(ExitCode::SUCCESS)
}

/// Build the engine's [`ProcessorOptions`] from parsed CLI [`Args`]. Lives here
/// (rather than as a `From` impl on `ProcessorOptions`) because `Args` is a
/// binary-crate type the library can't reference.
fn processor_options(args: &Args) -> ProcessorOptions {
    ProcessorOptions {
        remove_dups: args.remove_dups,
        add_mate_tags: args.add_mate_tags,
        ignore_unmated: args.ignore_unmated,
        max_read_length: args.max_read_length,
        single_end_strategy: args.single_end_strategy.to_strategy(),
        methylation_mode: args.methylation_mode.map(MethylationModeCli::to_mode),
    }
}

/// Input format auto-detected from the first byte.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    /// SAM text input (first byte is `@` = 0x40).
    Sam,
    /// BAM binary input (first byte is BGZF magic 0x1f).
    Bam,
}

/// Peek the first byte of `reader` (without consuming it) to determine
/// whether the input is SAM (`@`) or BAM (BGZF magic `0x1f`).
fn detect_format<R: BufRead>(reader: &mut R) -> Result<Format> {
    let head = reader.fill_buf().context("reading input to detect format")?;
    if head.is_empty() {
        bail!("Empty input: no SAM/BAM header detected");
    }
    match head[0] {
        0x1f => Ok(Format::Bam),
        b'@' => Ok(Format::Sam),
        b => bail!(
            "Input doesn't look like SAM or BAM (first byte 0x{:02x}); expected '@' or 0x1f",
            b
        ),
    }
}

/// Unified reader enum so the main loop is format-agnostic.
enum Reader {
    /// BAM binary reader (BGZF-decompressed via [`RawBamReader`]).
    Bam(RawBamReader<Box<dyn BufRead>>),
    /// SAM text reader (line-parsed via [`SamReader`]).
    Sam(SamReader<Box<dyn BufRead>>),
}

impl Reader {
    /// Read and parse the SAM/BAM header, leaving the reader positioned at the
    /// first alignment record.
    fn read_header(&mut self) -> Result<Header> {
        match self {
            Reader::Bam(r) => r.read_header().context("reading BAM header"),
            Reader::Sam(r) => r.read_header().context("reading SAM header"),
        }
    }

    /// Read one alignment record into `rec`. Returns `Ok(true)` on success,
    /// `Ok(false)` at EOF.
    fn read_record(&mut self, rec: &mut RawRecord) -> std::io::Result<bool> {
        match self {
            Reader::Bam(r) => r.read_record(rec),
            Reader::Sam(r) => r.read_record(rec),
        }
    }
}

/// Drive the QNAME-block grouping loop over `reader`, invoking `on_block`
/// for each maximal run of records sharing a QNAME.
///
/// `pool` is a growing `Vec<RawRecord>` whose allocations are reused across
/// blocks (and across calls — picard-exact reuses the same loop for its
/// second pass). On a QNAME change we hand the active slice to `on_block`,
/// swap the just-read record into slot 0, and continue.
fn for_each_block(
    reader: &mut Reader,
    pool: &mut Vec<RawRecord>,
    mut on_block: impl FnMut(&mut [RawRecord]) -> Result<()>,
) -> Result<()> {
    if pool.is_empty() {
        pool.push(RawRecord::new());
    }
    let mut block_len: usize = 0;
    let mut current_qname: Vec<u8> = Vec::new();
    loop {
        if pool.len() == block_len {
            pool.push(RawRecord::new());
        }
        let read_idx = block_len;
        let got = reader.read_record(&mut pool[read_idx]).context("reading record")?;
        if !got {
            break;
        }
        let new_qname_differs =
            block_len > 0 && pool[read_idx].read_name() != current_qname.as_slice();
        if new_qname_differs {
            on_block(&mut pool[..block_len])?;
            if read_idx != 0 {
                pool.swap(0, read_idx);
            }
            block_len = 1;
            current_qname.clear();
            current_qname.extend_from_slice(pool[0].read_name());
        } else {
            if block_len == 0 {
                current_qname.clear();
                current_qname.extend_from_slice(pool[0].read_name());
            }
            block_len += 1;
        }
    }
    if block_len > 0 {
        on_block(&mut pool[..block_len])?;
    }
    Ok(())
}

/// Two-pass driver for [`SingleEndStrategy::PicardExact`].
///
/// **Pass 1** streams pairs and unmapped/unmated blocks straight to `out`
/// while buffering every mapped-orphan / single-end block to a temporary
/// uncompressed BAM. **Transition** drains the completed pair table into a
/// fragment table holding every paired read end. **Pass 2** re-reads the
/// buffered blocks and marks each fragment that collides with a pair end (or
/// an earlier fragment), emitting them after the pairs.
///
/// The temp file is created via [`create_temp_bam`] and unlinked when the
/// `NamedTempFile` drops at the end of this function — after both the writer
/// and reader handles are closed.
fn run_picard_exact(
    reader: &mut Reader,
    processor: &mut RecordProcessor,
    stats: &mut Stats,
    out: &mut RawBamWriter,
    header: &Header,
    args: &Args,
) -> Result<()> {
    // A modest ring is plenty: fragments are a small fraction of paired data,
    // and this writer targets local disk rather than a bursty downstream.
    const TEMP_RING_BYTES: usize = 8 * 1024 * 1024;
    let stored = CompressionLevel::new(0).expect("compression level 0 is valid");

    let temp = create_temp_bam(args.tmp_dir.as_deref())?;
    let temp_path = temp.path().to_path_buf();
    if !args.quiet {
        eprintln!(
            "dupblaster: picard-exact mode — buffering orphan/single-end reads to {}",
            temp_path.display()
        );
    }

    // Pass 1: pairs out, fragments deferred to the temp BAM.
    {
        let mut temp_writer = RawBamWriter::open(Some(&temp_path), header, TEMP_RING_BYTES, stored)
            .context("opening temp orphan BAM")?;
        let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
        for_each_block(reader, &mut pool, |block| {
            processor.process_block_phase1(block, stats, out, &mut temp_writer)
        })
        .context("processing record block (picard-exact pass 1)")?;
        temp_writer.finish().context("finishing temp orphan BAM")?;
    }

    // Transition: build the fragment table from every pair end.
    processor.finalize_fragment_table();

    // Pass 2: re-read the buffered fragments and mark them against the
    // now-complete fragment table. CRC verification is off — we wrote this
    // file ourselves moments ago.
    let file = File::open(&temp_path).context("reopening temp orphan BAM for pass 2")?;
    let buffered: Box<dyn BufRead> = Box::new(std::io::BufReader::new(file));
    let mut temp_reader = Reader::Bam(RawBamReader::new(buffered, false));
    temp_reader.read_header().context("reading temp orphan BAM header")?;
    let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
    for_each_block(&mut temp_reader, &mut pool, |block| {
        processor.process_fragment_block(block, stats, out)
    })
    .context("processing fragment block (picard-exact pass 2)")?;

    Ok(())
}

/// Create the temporary uncompressed BAM used by picard-exact to buffer
/// fragment blocks. Honors `dir` (`--tmp-dir`) when given, else the system
/// temp dir (`$TMPDIR`). The returned handle unlinks the file on drop.
fn create_temp_bam(dir: Option<&std::path::Path>) -> Result<tempfile::NamedTempFile> {
    let mut builder = tempfile::Builder::new();
    builder.prefix("dupblaster-orphans-").suffix(".bam");
    let file = match dir {
        Some(d) => builder.tempfile_in(d),
        None => builder.tempfile(),
    }
    .context("creating temp file for picard-exact orphan buffering")?;
    Ok(file)
}

/// Append dupblaster's @PG line to the noodles header.
///
/// The `CL:` tag records the actual `std::env::args()` (with `argv[0]`
/// normalized to its basename), so the recorded command-line includes
/// every flag the user passed instead of a hand-curated subset.
///
/// noodles auto-chains via `PP:`: for every existing chain-leaf @PG, our
/// new record is appended as its child. If an `ID:DUPBLASTER` record
/// already exists (e.g. the input was already processed by dupblaster),
/// noodles disambiguates by appending `-{previous_id}` to our ID.
fn append_dupblaster_pg(header: &mut Header) -> Result<()> {
    // Defensive validation: every existing @PG with a PP: tag must point
    // at an ID that exists in the header. noodles' `programs.add(...)`
    // walks the existing chain to find leaves and panics with "no entry
    // found for key" if it encounters a broken PP reference. Convert
    // that into a clean error here so malformed input headers fail
    // cleanly instead of crashing.
    let programs = header.programs().as_ref();
    let known_ids: std::collections::HashSet<&[u8]> =
        programs.keys().map(|k| k.as_slice()).collect();
    for (id, map) in programs.iter() {
        if let Some(pp) = map.other_fields().get(&program_tag::PREVIOUS_PROGRAM_ID) {
            let pp_bytes = pp.as_ref();
            if !known_ids.contains(pp_bytes) {
                bail!(
                    "input header @PG ID:{} has PP:{} but no @PG with that ID exists. \
                     This is a malformed SAM header. Strip the broken PP tag or rewrite \
                     the @PG chain (e.g. via `samtools reheader`) before re-running.",
                    String::from_utf8_lossy(id.as_slice()),
                    String::from_utf8_lossy(pp_bytes),
                );
            }
        }
    }

    let cl = command_line_for_pg();
    let mut map = Map::<Program>::default();
    map.other_fields_mut().insert(program_tag::VERSION, DUPBLASTER_BUILD.into());
    map.other_fields_mut().insert(program_tag::COMMAND_LINE, cl.into());
    header.programs_mut().add("DUPBLASTER", map).context("appending @PG DUPBLASTER record")?;
    Ok(())
}

/// Build the `@PG CL:` string from `std::env::args`. Replaces `argv[0]`
/// with its file-name component so the recorded command isn't tied to the
/// absolute install path (e.g. `/nix/store/.../bin/dupblaster` → `dupblaster`).
/// Args are space-joined; paths with spaces aren't quoted (matches the
/// upstream samblaster convention).
fn command_line_for_pg() -> String {
    let mut args = std::env::args();
    let prog = args
        .next()
        .map(|a| {
            std::path::Path::new(&a)
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or(a)
        })
        .unwrap_or_else(|| "dupblaster".to_string());
    let rest: Vec<String> = args.collect();
    if rest.is_empty() { prog } else { format!("{prog} {}", rest.join(" ")) }
}

/// Print the end-of-run duplicate summary to stderr, including the coordinate-
/// clamping warning when any templates were affected.
fn print_run_stats(stats: &Stats, args: &Args) {
    // The stderr summary is run-wide; per-library breakdowns go to `--stats`.
    let totals = stats.totals();
    if totals.id_count == 0 {
        eprintln!("dupblaster: No reads processed.");
        return;
    }
    if stats.clamped_template_count > 0 {
        let pct = 100.0 * stats.clamped_template_count as f64 / totals.id_count as f64;
        eprintln!(
            "dupblaster: WARNING: {} of {} ({:.3}%) templates had a read whose 5' coordinate \
             was clamped to its contig because its clipping extends more than \
             --max-read-length({}) bases past a contig edge.",
            stats.clamped_template_count, totals.id_count, pct, args.max_read_length
        );
        eprintln!(
            "dupblaster: Duplicate marking may be imprecise for those templates; re-run with a \
             larger --max-read-length to eliminate the clamping. (Exiting non-zero.)"
        );
    }
    if args.ignore_unmated {
        let pct = 100.0 * totals.unmated_count as f64 / totals.id_count as f64;
        eprintln!(
            "dupblaster: Found {:>10} of {:>10} ({:5.3}%) total read ids are marked paired yet are unmated.",
            totals.unmated_count, totals.id_count, pct
        );
        if totals.unmated_count > 0 {
            eprintln!(
                "dupblaster: Please double check that input file is query-grouped (QNAME grouped)."
            );
        }
    }
    let verb = if args.remove_dups { "Removed" } else { "Marked " };
    let pct = 100.0 * totals.dup_count as f64 / totals.id_count as f64;
    eprintln!(
        "dupblaster: {} {:>10} of {:>10} ({:5.3}%) total read ids as duplicates.",
        verb, totals.dup_count, totals.id_count, pct
    );
}

// ── End-of-run resource-usage footer ────────────────────────────────────────
//
// Mirrors the C++ samblaster footer that prints memory + timing after the
// dup-marking summary. Suppressed by `--quiet`. CPU and max-RSS reads use
// `getrusage(RUSAGE_SELF)` and are Unix-only — on other platforms only wall
// time is reported.

/// Captures a starting wall-clock timestamp; pair with [`report`] at end.
struct StartedRun {
    /// Snapshot of `Instant::now()` at process start.
    wall_start: Instant,
}

impl StartedRun {
    /// Snapshot the current wall-clock timestamp.
    fn now() -> Self {
        Self { wall_start: Instant::now() }
    }
}

/// Print a single-line resource-usage footer to stderr. No-op if `quiet`.
fn report(started: &StartedRun, n_templates: u64, quiet: bool) {
    if quiet {
        return;
    }
    let wall = started.wall_start.elapsed().as_secs_f64();
    let stderr = std::io::stderr();
    let mut stderr = stderr.lock();

    #[cfg(unix)]
    if let Some(ru) = read_rusage() {
        let user = ru.user_secs;
        let sys = ru.sys_secs;
        let rss_mb = ru.max_rss_bytes as f64 / (1024.0 * 1024.0);
        let _ = writeln!(
            stderr,
            "dupblaster: Processed {n_templates} templates in {wall:.2}s wall, \
             {user:.2}s user CPU, {sys:.2}s system CPU, max RSS {rss_mb:.1} MB.",
        );
        return;
    }

    let _ = writeln!(stderr, "dupblaster: Processed {n_templates} templates in {wall:.2}s wall.");
}

/// Normalized resource-usage snapshot from `getrusage(RUSAGE_SELF)`.
#[cfg(unix)]
struct Rusage {
    /// User-mode CPU time in seconds.
    user_secs: f64,
    /// Kernel-mode CPU time in seconds.
    sys_secs: f64,
    /// Peak resident-set size in bytes (normalized from OS-native units).
    max_rss_bytes: u64,
}

/// Snapshot the current process's user/sys CPU and max RSS via
/// `getrusage(RUSAGE_SELF)`. Returns `None` if the syscall fails.
///
/// `ru_maxrss` is in **bytes** on macOS and **kilobytes** on Linux — we
/// normalize to bytes here so the caller doesn't need to care.
#[cfg(unix)]
fn read_rusage() -> Option<Rusage> {
    // SAFETY: `getrusage` writes a `rusage` struct that we just allocated.
    // Both the syscall number and the struct layout come from `libc`.
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) };
    if rc != 0 {
        return None;
    }
    let user_secs = ru.ru_utime.tv_sec as f64 + ru.ru_utime.tv_usec as f64 * 1e-6;
    let sys_secs = ru.ru_stime.tv_sec as f64 + ru.ru_stime.tv_usec as f64 * 1e-6;
    let max_rss = ru.ru_maxrss as u64;
    // macOS reports bytes; Linux + BSD report kilobytes.
    #[cfg(target_os = "macos")]
    let max_rss_bytes = max_rss;
    #[cfg(not(target_os = "macos"))]
    let max_rss_bytes = max_rss.saturating_mul(1024);
    Some(Rusage { user_secs, sys_secs, max_rss_bytes })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args_with_output(out: Option<&str>) -> Args {
        Args {
            input: None,
            output: out.map(PathBuf::from),
            remove_dups: false,
            add_mate_tags: false,
            ignore_unmated: false,
            max_read_length: 1000,
            library_unaware: false,
            single_end_strategy: SingleEndStrategyCli::StrandAware,
            methylation_mode: None,
            tmp_dir: None,
            stats: None,
            sample: None,
            quiet: true,
            min_bins: 32,
            check_crc: false,
            no_check_crc: false,
            read_buffer_mb: 16,
            write_buffer_mb: 64,
            compression_level: CompressionLevel::new(0).expect("level 0 is valid"),
            show_version: false,
        }
    }

    #[test]
    fn validate_accepts_bam_extension() {
        assert!(args_with_output(Some("out.bam")).validate().is_ok());
        assert!(args_with_output(Some("/tmp/sample-A.bam")).validate().is_ok());
    }

    #[test]
    fn validate_accepts_stdout_dash() {
        assert!(args_with_output(Some("-")).validate().is_ok());
    }

    #[test]
    fn validate_accepts_no_output_specified() {
        assert!(args_with_output(None).validate().is_ok());
    }

    #[test]
    fn validate_rejects_non_bam_extension() {
        let err = args_with_output(Some("out.sam")).validate().unwrap_err();
        assert!(err.to_string().contains("must end in `.bam`"));
    }

    #[test]
    fn validate_rejects_missing_extension() {
        let err = args_with_output(Some("out")).validate().unwrap_err();
        assert!(err.to_string().contains("must end in `.bam`"));
    }

    #[test]
    fn max_read_length_rejects_non_positive() {
        assert!(Args::try_parse_from(["dupblaster", "--max-read-length", "0"]).is_err());
        assert!(Args::try_parse_from(["dupblaster", "--max-read-length", "-5"]).is_err());
    }

    #[test]
    fn max_read_length_rejects_overflowing_value() {
        // Above the 10M cap that keeps 2*max_read_length within i32.
        assert!(Args::try_parse_from(["dupblaster", "--max-read-length", "20000000"]).is_err());
    }

    #[test]
    fn max_read_length_accepts_in_range() {
        let args = Args::try_parse_from(["dupblaster", "--max-read-length", "150000"]).unwrap();
        assert_eq!(args.max_read_length, 150_000);
    }

    #[test]
    fn min_bins_rejects_zero_and_oversized() {
        assert!(Args::try_parse_from(["dupblaster", "--min-bins", "0"]).is_err());
        assert!(Args::try_parse_from(["dupblaster", "--min-bins", "100000"]).is_err());
    }

    #[test]
    fn min_bins_accepts_in_range() {
        let args = Args::try_parse_from(["dupblaster", "--min-bins", "8192"]).unwrap();
        assert_eq!(args.min_bins, 8192);
    }

    #[test]
    fn effective_check_crc_true_when_flag_set() {
        let mut a = args_with_output(None);
        a.check_crc = true;
        assert!(a.effective_check_crc());
    }

    #[test]
    fn effective_check_crc_false_when_no_check_flag_set() {
        let mut a = args_with_output(None);
        a.no_check_crc = true;
        assert!(!a.effective_check_crc());
    }

    #[test]
    fn effective_check_crc_defaults_on_for_file_input() {
        let mut a = args_with_output(None);
        a.input = Some(PathBuf::from("reads.bam"));
        assert!(a.effective_check_crc());
    }

    #[test]
    fn effective_check_crc_defaults_off_for_stdin() {
        // input == None means stdin.
        assert!(!args_with_output(None).effective_check_crc());
    }

    #[test]
    fn effective_check_crc_defaults_off_for_stdin_dash() {
        let mut a = args_with_output(None);
        a.input = Some(PathBuf::from("-"));
        assert!(!a.effective_check_crc());
    }

    #[test]
    fn check_crc_and_no_check_crc_are_mutually_exclusive() {
        assert!(Args::try_parse_from(["dupblaster", "--check-crc", "--no-check-crc"]).is_err());
    }

    #[test]
    fn methylation_mode_defaults_to_none() {
        let args = Args::try_parse_from(["dupblaster"]).unwrap();
        assert_eq!(args.methylation_mode, None);
    }

    #[test]
    fn methylation_mode_parses_directional() {
        let args =
            Args::try_parse_from(["dupblaster", "--methylation-mode", "directional"]).unwrap();
        assert_eq!(args.methylation_mode, Some(MethylationModeCli::Directional));
    }

    #[test]
    fn methylation_mode_rejects_unknown_value() {
        // `pbat` is intentionally not implemented yet — it must be rejected at
        // parse time rather than silently treated as directional.
        assert!(Args::try_parse_from(["dupblaster", "--methylation-mode", "pbat"]).is_err());
    }
}
