//! dupblaster — the command-line application.
//!
//! This binary is the whole tool: CLI parsing ([`Args`]), end-to-end run
//! orchestration ([`run`]), and the end-of-run resource-usage footer. The
//! reusable engine and IO live in `dupblaster_lib` (`dedup`, `sig`, the
//! readers/writers, `metrics`, …) — this file wires them together and talks to
//! the user.
//!
//! Long-form flags follow GNU style (`--kebab-case`) — diverging from C++
//! samblaster's `--camelCase` convention. Short flags `-i`, `-o`, `-r`, and `-q`
//! match upstream; `-l` and `-m` are dupblaster additions. The SV-extraction
//! flags from upstream (`-d`, `-s`, `-u`, `-a`, `-e`, `-M`, `--maxSplitCount`,
//! …) have been intentionally dropped — see README for the rationale.

mod cigar;
mod complexity;
mod counts;
mod countset;
mod dedup;
mod metrics;
mod plots;
mod raw_reader;
mod raw_writer;
mod readname;
mod sam_reader;
mod sig;
mod tiles;

use std::fs::File;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
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
use rawb_io::ReadAhead;

use crate::complexity::{LadderRecorder, ladder_path, write_ladder_rows};
use crate::counts::{histogram_path, histogram_rows, write_histogram_rows};
use crate::dedup::{LibraryIndex, ProcessorOptions, RecordProcessor, Stats};
use crate::metrics::{
    Metrics, SequencingUnitMetrics, resolve_sample, sequencing_units_path, write_rows_to_path,
    write_unit_rows_to_path,
};
use crate::raw_reader::RawBamReader;
use crate::raw_writer::RawBamWriter;
use crate::readname::ReadNameFormat;
use crate::sam_reader::SamReader;
use crate::sig::{MethylationMode, SingleEndStrategy};
use crate::tiles::{DEFAULT_SPILL_BUCKETS, TileSpiller};

/// Crate-level build identifier shown in `--version` and the `@PG VN:` tag.
const DUPBLASTER_BUILD: &str = env!("CARGO_PKG_VERSION");

/// Upstream repository, shown in the `--help` banner and the startup banner.
/// A macro rather than a `const` so the one definition can be `concat!`-ed into
/// the compile-time `--help` text as well as printed at runtime.
macro_rules! repo_url {
    () => {
        "https://github.com/fulcrumgenomics/dupblaster"
    };
}

/// Attribution banner heading both `-h` and `--help`. It carries the pointer to
/// the full documentation, so no individual option's help has to.
macro_rules! help_banner {
    () => {
        concat!("dupblaster by Fulcrum Genomics\n", repo_url!(), "\n\n")
    };
}

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
    /// Identical to the position key used in paired-end data.
    #[value(name = "strand-aware")]
    StrandAware,
    /// Matches Picard MarkDuplicates exactly by buffering orphans to a temp BAM
    /// (see `--tmp-dir`) and processing after pairs. Emits orphans at the end of
    /// the output rather than in input order.
    #[value(name = "picard-exact")]
    PicardExact,
    /// Strand-aware, plus each end of a mapped pair is registered in a fragment
    /// table so later orphans are marked as dups of pairs. Order-sensitive in one
    /// pass, but an orphan preceding its pair passes as non-dup. Approximately
    /// doubles memory usage.
    #[value(name = "picard-approx")]
    PicardApprox,
    /// samblaster v0.1.23+ behavior: leftmost coordinate with the strand bit
    /// dropped, so forward and reverse orphans sharing a leftmost position
    /// collide. NOT recommended for short-read PE data; provided for
    /// byte-compatibility on long-read singleton workflows.
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
            Self::PicardExact => SingleEndStrategy::PicardExact,
            Self::PicardApprox => SingleEndStrategy::PicardApprox,
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
    /// Keys pairs in template order rather than coordinate-canonical order, so a
    /// fragment's two original strands (OT/OB) at one locus stay distinct while
    /// same-strand PCR copies still collapse. Does NOT handle non-directional /
    /// PBAT.
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

/// Bounds on `--tmp-compression-level`, deliberately far narrower than the
/// range libzstd advertises.
///
/// The ceiling is a memory limit, not a taste judgement. The duplicate-split spill
/// keeps one compression context per bucket alive for the whole main pass, so the
/// cost is `DEFAULT_SPILL_BUCKETS ×` one context — under 100 MB at level 1, but
/// tens of gigabytes at zstd's maximum, which would OOM exactly the constrained
/// runs this flag exists to rescue. Levels past 9 also barely shrink this data.
///
/// The floor is where the fast tiers stop giving anything back: below roughly -7
/// the temporary files stop getting smaller while still costing CPU.
const MIN_TMP_COMPRESSION_LEVEL: i32 = -7;
const MAX_TMP_COMPRESSION_LEVEL: i32 = 9;

/// Parse and validate a `--tmp-compression-level` value.
///
/// Rejects 0 rather than forwarding it: zstd reads level 0 as "use the default"
/// (3), while `--compression-level` in this same CLI reads 0 as "stored". One tool
/// cannot have 0 mean both things, so here it means neither and the error says so.
fn parse_tmp_compression_level(s: &str) -> Result<i32, String> {
    let level: i32 = s.parse().map_err(|_| format!("not an integer: {s}"))?;
    // libzstd stays the outer authority, in case a future release narrows below
    // the band above.
    let supported = zstd::compression_level_range();
    let (min, max) = (
        MIN_TMP_COMPRESSION_LEVEL.max(*supported.start()),
        MAX_TMP_COMPRESSION_LEVEL.min(*supported.end()),
    );
    if level == 0 {
        return Err("0 is ambiguous here: omit --tmp-compression-level entirely to leave \
                    temporary files uncompressed, or pass 3 for zstd's default level"
            .to_string());
    }
    if !(min..=max).contains(&level) {
        return Err(format!("must be between {min} and {max}"));
    }
    Ok(level)
}

/// `--help` heading for the IO and indexing knobs whose defaults suit
/// essentially every run — grouped away from the options a caller chooses
/// between.
const TUNING_HEADING: &str = "Advanced tuning (rarely needed)";

/// Short `about` text shown by `-h`.
const SHORT_ABOUT: &str =
    concat!(help_banner!(), "Mark or remove duplicate reads in a query-grouped SAM/BAM file.");
/// Full `about` text shown by `--help` (long form).
const LONG_ABOUT: &str = concat!(
    help_banner!(),
    "Mark or remove duplicates in a query-grouped SAM/BAM file.\n\
     Input must be query-grouped (same-QNAME alignments adjacent, as straight from \
     an aligner); output is always BAM, uncompressed unless -l says otherwise.",
);

/// Parsed command-line arguments — see `--help` for per-flag descriptions.
#[derive(Parser, Debug, Clone)]
#[command(name = "dupblaster", disable_version_flag = true, about = SHORT_ABOUT, long_about = LONG_ABOUT)]
pub struct Args {
    /// Input SAM/BAM file [default: stdin].
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output BAM file [default: stdout]. Must end in `.bam`, or be `-` for
    /// stdout.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// BGZF compression level for the output BAM, 0-12.
    #[arg(short = 'l',
          long = "compression-level",
          default_value = "0",
          value_parser = parse_compression_level)]
    pub compression_level: CompressionLevel,

    /// Remove duplicates from the output instead of flagging them.
    #[arg(short = 'r', long = "remove-dups")]
    pub remove_dups: bool,

    /// Add MC (mate CIGAR) and MQ (mate MAPQ) tags to all paired output.
    #[arg(short = 'm', long = "add-mate-tags")]
    pub add_mate_tags: bool,

    /// Suppress abort on unmated alignments.
    #[arg(long = "ignore-unmated")]
    pub ignore_unmated: bool,

    /// How to key single-end / orphan reads.
    #[arg(long = "single-end-strategy",
          value_enum,
          default_value_t = SingleEndStrategyCli::StrandAware)]
    pub single_end_strategy: SingleEndStrategyCli,

    /// Identify duplicates across libraries instead of within libraries.
    ///
    /// Libraries are identified from `@RG LB:` entries, and by default
    /// de-duplication is restricted to within a library, with statistics reported
    /// per-library. Has no effect when the header declares ≤1 library.
    #[arg(long = "library-unaware")]
    pub library_unaware: bool,

    /// Methylation-aware pair keying for strand-aware libraries; off by default.
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<MethylationModeCli>,

    /// Directory for temporary files [default: the system temp dir, `$TMPDIR`].
    #[arg(long = "tmp-dir")]
    pub tmp_dir: Option<PathBuf>,

    /// Write a TSV of per-library metrics to PATH.
    #[arg(long = "stats")]
    pub stats: Option<PathBuf>,

    /// Sample name for the `sample` column of the TSV outputs. Defaults to the
    /// unique `@RG SM:` values, comma-joined.
    #[arg(long = "sample")]
    pub sample: Option<String>,

    /// Write per-library duplication-complexity QC with this path prefix. Off by
    /// default, and costs nothing when unset.
    ///
    /// Writes `<PREFIX>.duplication-{sampled,spectrum}.{tsv,pdf}` containing
    /// i) duplication stats sampled every `--complexity-interval`, and ii)
    /// a family size histogram of how many molecules were seen how many times.
    #[arg(long = "complexity-metrics")]
    pub complexity_metrics: Option<PathBuf>,

    /// Snapshot cadence for the `--complexity-metrics` sampling, in templates.
    #[arg(long = "complexity-interval",
          default_value_t = 1_000_000,
          value_parser = clap::value_parser!(u64).range(1..))]
    pub complexity_interval: u64,

    /// Don't split duplicates into sequencing ("optical") and library components.
    ///
    /// Sequencing vs. library duplicate classification is on by default. Disabling
    /// i) does not affect the output BAM, ii) reduces tmp file usage by
    /// approximately 10-20MiB per 1m read pairs, iii) reduces CPU usage by 5-10%.
    #[arg(long = "no-sequencing-dups")]
    pub no_sequencing_dups: bool,

    /// Read/query name format in the BAM file - used to extract flowcell/lane/tile
    /// or equivalent for sequencing duplicate classification [default: illumina].
    ///
    /// `illumina`: instrument:run:flowcell:lane:tile:x:y (CASAVA 1.8+, bcl2fastq,
    /// BCL Convert).
    ///
    /// `element`: matches the illumina format, for Element AVITI.
    ///
    /// `regex:PATTERN`: PATTERN needs a `(?<su>...)` capture group for the
    /// sequencing unit (e.g. flowcell+lane) and a `(?<tile>...)` group for the
    /// tile. The groups are found by name, not by position, and a pattern missing
    /// either is rejected up front. Both tokens are compared as opaque text, so
    /// they need no particular shape.
    ///
    /// `illumina` written as a regex, to adapt for another platform:
    ///
    /// `--read-name-format 'regex:^[^:]+:[^:]+:(?<su>[^:]+:[^:]+):(?<tile>[^:]+):[^:]+:[^:]+'`
    ///
    /// The unit spans two fields because a tile number only identifies a place
    /// within one lane of one flowcell.
    #[arg(long = "read-name-format", value_name = "FORMAT")]
    pub read_name_format: Option<ReadNameFormat>,

    /// Print only the duplicate summary on stderr, not the progress and
    /// resource-usage lines.
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Print the dupblaster version to stdout and exit.
    #[arg(short = 'V', long = "version")]
    pub show_version: bool,

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

    /// Verify BGZF CRC32 on input. Defaults to on for files and off for stdin.
    #[arg(long = "check-crc", conflicts_with = "no_check_crc", help_heading = TUNING_HEADING)]
    pub check_crc: bool,

    /// Skip BGZF CRC32 verification whatever the input source.
    #[arg(long = "no-check-crc", conflicts_with = "check_crc", help_heading = TUNING_HEADING)]
    pub no_check_crc: bool,

    /// Size in MB of the ring buffer between the input IO thread and the worker.
    ///
    /// The default absorbs typical aligner bursts (bwa-mem emits ~250 MB every
    /// few seconds); raise it for a faster or burstier producer.
    #[arg(long = "read-buffer-mb", default_value_t = 16, help_heading = TUNING_HEADING,
          value_parser = clap::value_parser!(u32).range(1..=4096))]
    pub read_buffer_mb: u32,

    /// Size in MB of the ring buffer between the worker and the output IO thread.
    ///
    /// Bigger than the read buffer because downstream sorters pause input for
    /// 1-3 s while flushing a chunk to disk; raise it to 256+ for slower sinks.
    #[arg(long = "write-buffer-mb", default_value_t = 64, help_heading = TUNING_HEADING,
          value_parser = clap::value_parser!(u32).range(1..=4096))]
    pub write_buffer_mb: u32,

    /// Maximum expected read length in bp; raise it for long reads.
    ///
    /// Sets the per-contig padding in the duplicate-detection index and the
    /// threshold for the "reads longer than this" warning; no algorithmic effect
    /// at short-read sizes.
    // The 10 Mbp ceiling keeps `2 * max_read_length` (the per-contig padding)
    // well inside `i32`, so the super-contig length arithmetic can't overflow.
    #[arg(long = "max-read-length", default_value_t = 1000, help_heading = TUNING_HEADING,
          value_parser = clap::value_parser!(i32).range(1..=10_000_000))]
    pub max_read_length: i32,

    /// Ztsd compression level used to compress temporary files (-7 to 9).
    ///
    /// Unspecified uses no compression.  Levels -5 to 1 tend to offer
    /// meaningful (1.5-2x) reduction in temp usage, for marginal (2-7%)
    /// more CPU usage on realistic sized data.
    #[arg(long = "tmp-compression-level", value_name = "LEVEL",
          help_heading = TUNING_HEADING, allow_negative_numbers = true,
          value_parser = parse_tmp_compression_level)]
    pub tmp_compression_level: Option<i32>,
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
/// Fail fast if the directory that would hold `path` (an output file, or the
/// prefix the `--complexity-metrics` files derive from) is missing or not
/// writable — so a bad location errors up front instead of after a full
/// processing pass (those files are written only at the end). Probes by creating
/// and removing a temp file in that directory.
fn ensure_dir_writable(path: &Path, flag: &str) -> Result<()> {
    let dir = match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    };
    if !dir.is_dir() {
        bail!("{flag}: output directory {} does not exist", dir.display());
    }
    let probe = dir.join(format!(".dupblaster-write-probe.{}", std::process::id()));
    std::fs::File::create(&probe)
        .and_then(|_| std::fs::remove_file(&probe))
        .with_context(|| format!("{flag}: cannot write to output directory {}", dir.display()))?;
    Ok(())
}

fn run(args: Args) -> Result<ExitCode> {
    if args.show_version {
        // Version goes to stdout (not stderr) so `dupblaster --version | head`
        // and tooling that captures it work as expected.
        println!("dupblaster {DUPBLASTER_BUILD}");
        return Ok(ExitCode::SUCCESS);
    }
    args.validate()?;
    // Fail fast if the optional QC outputs (written only at the very end) can't
    // be created, rather than after a full processing pass. The main `-o` output
    // is validated separately when its writer is opened, below.
    if let Some(stats) = args.stats.as_deref() {
        ensure_dir_writable(stats, "--stats")?;
    }
    if let Some(prefix) = args.complexity_metrics.as_deref() {
        // All four complexity files share the prefix's directory.
        ensure_dir_writable(prefix, "--complexity-metrics")?;
    }
    let started = StartedRun::now();
    if !args.quiet {
        eprintln!("dupblaster: Version {DUPBLASTER_BUILD} by Fulcrum Genomics");
        eprintln!("dupblaster: {}", repo_url!());
    }

    // Open input. The first byte tells us SAM vs BAM: 0x1f is the BGZF
    // gzip magic; '@' (0x40) is a SAM header line. Anything else is an
    // error.
    //
    // Input goes through a dedicated IO read thread + ring buffer
    // (a `rawb_io::ReadAhead`) so the worker never blocks on the kernel pipe.
    let raw_source: Box<dyn std::io::Read + Send> = match args.input.as_deref() {
        Some(p) if p.to_string_lossy() != "-" => {
            let f = File::open(p).with_context(|| format!("opening {} for read", p.display()))?;
            Box::new(f)
        }
        _ => Box::new(std::io::stdin()),
    };
    let read_buf_bytes = (args.read_buffer_mb as usize).saturating_mul(1024 * 1024);
    let mut reader_box: Box<dyn BufRead> =
        Box::new(ReadAhead::with_thread_name(raw_source, read_buf_bytes, "dupblaster"));
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
    // Attached after construction because the spiller is sized against
    // `bin_count`, which only exists once the bins are built.
    if !args.no_sequencing_dups {
        let spiller = TileSpiller::new(
            args.read_name_format.clone().unwrap_or_default(),
            processor.bin_count(),
            DEFAULT_SPILL_BUCKETS,
            args.tmp_dir.as_deref(),
            args.tmp_compression_level,
        )
        .context("preparing the sequencing-vs-library duplicate spill")?;
        processor.attach_tile_spiller(spiller);
    }
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

    // Optional duplication-complexity metrics. The recorder and per-library
    // counters are constructed only when `--complexity-metrics` is set, so the
    // default path allocates nothing and its performance/RSS are unchanged.
    // Works under every single-end strategy: each library is reported on its
    // paired subset when it has any pairs (a strategy-independent key), and on
    // single-end signatures only when it has none — the case where the fragment
    // keyspace is guaranteed free of pair-ends and so is clean under Picard too.
    let mut ladder = args.complexity_metrics.as_ref().map(|_| {
        LadderRecorder::new(
            num_libs as usize,
            args.complexity_interval,
            resolve_sample(&header, args.sample.as_deref()),
        )
    });

    // dupblaster requires query-grouped input: records are read in maximal
    // runs sharing a QNAME ("blocks") and each block is processed as a unit.
    // `for_each_block` owns that grouping loop. picard-exact takes a separate
    // two-pass driver (pairs stream out, fragments are buffered then
    // re-processed); every other strategy is a single streaming pass.
    let processed = if args.single_end_strategy.to_strategy() == SingleEndStrategy::PicardExact {
        run_picard_exact(
            &mut reader,
            &mut processor,
            &mut stats,
            &mut out,
            &header,
            &args,
            &mut ladder,
        )
    } else {
        let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
        for_each_block(&mut reader, &mut pool, |block| {
            let lib = processor.process_block(block, &mut stats, &mut out)?;
            // One comparison per template in the common case; the recorder
            // snapshots only when a library crosses an interval boundary.
            if let Some(rec) = ladder.as_mut() {
                rec.observe(lib, &stats.libraries[lib as usize]);
            }
            Ok(())
        })
        .context("processing record block")
    };
    // A run that died part-way has written a prefix of its output. Leave that
    // prefix on disk — silently discarding a user's partial output would be its
    // own surprise — but abandon the writer rather than finishing it, so the file
    // lacks the BGZF EOF marker and every reader can see it is incomplete.
    // Dropping the writer instead would emit that marker and make the truncated
    // file pass `samtools quickcheck`.
    if let Err(error) = processed {
        out.abandon();
        return Err(error);
    }

    print_run_stats(&stats, &args);

    // Finalize the primary BAM output *before* writing the optional QC side
    // files (`--stats`, `--complexity-metrics`). Their rows/counts are already
    // fully in memory, so a hiccup writing them (a bad prefix path, a full disk,
    // a kuva render error) must not abort the run before `finish()` has flushed
    // the BGZF EOF block and surfaced any main-output write error — an error
    // that relying on `Drop` alone would silently swallow.
    out.finish().context("finishing main output")?;

    // Only now that the output stream is closed may the decomposition read its
    // spill back: it touches gigabytes and can take tens of seconds, and
    // dupblaster sits in pipelines where a downstream sort or caller is blocked
    // on our stdout. This ordering is a hard requirement, not a preference.
    let decomposition = match processor.take_tile_spiller() {
        Some(mut spiller) => {
            // Closes the bucket streams, which the post-pass requires before it can
            // decode them, and yields the on-disk total while the files are still
            // there. Only measured when compression is on: statting every bucket is
            // pointless work on the default path, where the answer is the logical
            // size by construction.
            let measure = !args.quiet && args.tmp_compression_level.is_some();
            let on_disk =
                spiller.finish_spill(measure).context("closing the duplicate-split spill")?;
            if !args.quiet {
                let logical = spiller.spilled_bytes();
                let mib = |bytes: u64| bytes as f64 / (1024.0 * 1024.0);
                // A ratio needs something to divide; an empty spill has none. Note
                // the figure can legitimately come out below 1.0 on a tiny input,
                // where per-bucket framing outweighs what there was to compress.
                let compressed = match (on_disk, logical) {
                    (Some(on_disk), logical) if logical > 0 => format!(
                        " ({:.1} MiB on disk, {:.2}x)",
                        mib(on_disk),
                        logical as f64 / on_disk.max(1) as f64
                    ),
                    _ => String::new(),
                };
                eprintln!(
                    "dupblaster: decomposing duplicates from {:.1} MiB of spilled tile \
                     records{compressed}",
                    mib(logical),
                );
            }
            Some(
                spiller
                    .decompose(num_libs)
                    .context("decomposing sequencing vs library duplicates")?,
            )
        }
        None => None,
    };

    if let Some(stats_path) = args.stats.as_deref() {
        let rows = Metrics::rows_from_stats(
            &stats,
            &header,
            args.sample.as_deref(),
            decomposition.as_ref().map(|result| result.libraries.as_slice()),
        );
        write_rows_to_path(&rows, stats_path).context("writing --stats TSV")?;

        // The per-unit table is a different granularity (one flowcell-and-lane
        // per row), so it gets its own file beside the per-library one.
        if let Some(result) = decomposition.as_ref() {
            let sample = resolve_sample(&header, args.sample.as_deref());
            let rows = SequencingUnitMetrics::rows(&result.units, &stats, &sample);
            write_unit_rows_to_path(&rows, &sequencing_units_path(stats_path))
                .context("writing sequencing-unit TSV")?;
        }
    }

    // `--complexity-metrics <PREFIX>`, shared by the two QC blocks below.
    let prefix = args.complexity_metrics.as_deref();

    // Emit the duplicate-rate ladder (only when `--complexity-metrics` was set,
    // which is the only case `ladder` is `Some`). `finalize` appends the final
    // per-library snapshot and sorts rows for a stable file order.
    if let (Some(mut rec), Some(prefix)) = (ladder, prefix) {
        rec.finalize(&stats);
        let path = ladder_path(prefix);
        write_ladder_rows(rec.rows(), &path).context("writing duplication-sampled TSV")?;
        plots::write_duplicate_ladder_pdf(rec.rows(), &plots::ladder_plot_path(prefix))
            .context("writing duplication-sampled plot")?;
    }

    // Group-size histogram (η_k), from the per-library occurrence counts. Only
    // populated when `--complexity-metrics` was set (so `counts()` is `Some`).
    if let (Some(prefix), Some(counts)) = (prefix, processor.counts()) {
        let sample = resolve_sample(&header, args.sample.as_deref());
        let rows = histogram_rows(counts, &stats, &sample);
        write_histogram_rows(&rows, &histogram_path(prefix))
            .context("writing duplication-spectrum TSV")?;
        plots::write_count_histogram_pdfs(&rows, prefix)
            .context("writing duplication-spectrum plot")?;
    }

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
        collect_counts: args.complexity_metrics.is_some(),
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
/// BAM. **Transition** drains the completed pair table into a
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
    ladder: &mut Option<LadderRecorder>,
) -> Result<()> {
    // A modest ring is plenty: fragments are a small fraction of paired data,
    // and this writer targets local disk rather than a bursty downstream.
    const TEMP_RING_BYTES: usize = 8 * 1024 * 1024;

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
        let mut temp_writer = RawBamWriter::open_temp(
            &temp_path,
            header,
            TEMP_RING_BYTES,
            args.tmp_compression_level,
        )
        .context("opening temp orphan BAM")?;
        let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
        for_each_block(reader, &mut pool, |block| {
            let lib = processor.process_block_phase1(block, stats, out, &mut temp_writer)?;
            if let Some(rec) = ladder.as_mut() {
                rec.observe(lib, &stats.libraries[lib as usize]);
            }
            Ok(())
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
    // Unlike the write side, this runs on the worker thread — there is no
    // read-ahead for the temp file — so decompression is not free here.
    let buffered: Box<dyn BufRead> = match args.tmp_compression_level {
        None => Box::new(std::io::BufReader::new(file)),
        Some(_) => Box::new(std::io::BufReader::new(
            zstd::stream::read::Decoder::new(file).context("starting zstd decode of temp BAM")?,
        )),
    };
    let mut temp_reader = Reader::Bam(RawBamReader::new(buffered, false));
    temp_reader.read_header().context("reading temp orphan BAM header")?;
    let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
    for_each_block(&mut temp_reader, &mut pool, |block| {
        let lib = processor.process_fragment_block(block, stats, out)?;
        if let Some(rec) = ladder.as_mut() {
            rec.observe(lib, &stats.libraries[lib as usize]);
        }
        Ok(())
    })
    .context("processing fragment block (picard-exact pass 2)")?;

    Ok(())
}

/// Create the temporary BAM used by picard-exact to buffer fragment blocks.
/// Honors `dir` (`--tmp-dir`) when given, else the system temp dir (`$TMPDIR`).
/// The returned handle unlinks the file on drop.
///
/// Whether its contents are compressed is [`RawBamWriter::open_temp`]'s business,
/// driven by `--tmp-compression-level`.
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
            complexity_metrics: None,
            complexity_interval: 1_000_000,
            no_sequencing_dups: true,
            read_name_format: None,
            tmp_compression_level: None,
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

    #[test]
    fn short_l_sets_compression_level() {
        let args = Args::try_parse_from(["dupblaster", "-l", "6"]).unwrap();
        assert_eq!(u8::from(args.compression_level), 6);
    }

    #[test]
    fn short_m_sets_add_mate_tags() {
        let args = Args::try_parse_from(["dupblaster", "-m"]).unwrap();
        assert!(args.add_mate_tags);
    }
}
