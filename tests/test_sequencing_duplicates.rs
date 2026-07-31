//! Integration tests for `--read-name-format`: splitting duplicates into
//! sequencing (on-flowcell) and library components end to end through the CLI.

mod helpers;

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use helpers::*;

/// An Illumina read name on `flowcell`/`lane`/`tile`, made unique by `cluster`.
///
/// The x/y fields have to vary even though nothing here reads them: two distinct
/// molecules imaged on one tile are two templates with two QNAMEs, and reusing a
/// name would make dupblaster see one template with four records instead of four
/// duplicate templates. Exactly the same is true on a real flowcell.
fn read_name(flowcell: &str, lane: u32, tile: u32, cluster: u32) -> String {
    format!("A00354:1305:{flowcell}:{lane}:{tile}:{cluster}:1986")
}

/// Run dupblaster with the decomposition on.
fn run(input: &Path, stats: &Path, output: &Path, extra: &[&str]) -> std::process::Output {
    let mut cmd = Command::new(rust_binary());
    cmd.args(["-i"])
        .arg(input)
        .args(["-o"])
        .arg(output)
        .args(["--stats"])
        .arg(stats)
        .args(["--read-name-format", "illumina"])
        .args(extra);
    cmd.output().expect("dupblaster ran")
}

/// Run and require success.
fn run_ok(input: &Path, stats: &Path, output: &Path, extra: &[&str]) {
    let out = run(input, stats, output, extra);
    assert!(out.status.success(), "dupblaster failed: {}", String::from_utf8_lossy(&out.stderr));
}

/// Parse a two-line TSV into a column → value map.
fn parse_row(path: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(path).expect("read TSV");
    let mut lines = text.lines();
    let cols: Vec<&str> = lines.next().expect("header line").split('\t').collect();
    let vals: Vec<&str> = lines.next().expect("value line").split('\t').collect();
    assert_eq!(cols.len(), vals.len(), "column count mismatch");
    cols.into_iter().map(String::from).zip(vals.into_iter().map(String::from)).collect()
}

/// Parse every data row of a TSV into column → value maps.
fn parse_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).expect("read TSV");
    let mut lines = text.lines();
    let cols: Vec<&str> = lines.next().expect("header line").split('\t').collect();
    lines
        .map(|line| {
            let vals: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), vals.len(), "column count mismatch");
            cols.iter().map(|c| c.to_string()).zip(vals.into_iter().map(String::from)).collect()
        })
        .collect()
}

/// Accumulates SAM input, handing out a fresh cluster coordinate per template.
struct Run {
    builder: SamBuilder,
    /// Next x coordinate, so every template gets a distinct QNAME.
    cluster: u32,
}

impl Run {
    fn new() -> Self {
        Self { builder: SamBuilder::new().sq("chr1", 10_000_000), cluster: 0 }
    }

    /// Add a both-ends-mapped pair on `flowcell`/`lane`/`tile`, aligned at `pos`.
    /// Two pairs sharing a `pos` are duplicates of each other.
    fn pair(&mut self, flowcell: &str, lane: u32, tile: u32, pos: u32) -> &mut Self {
        self.pair_rg(flowcell, lane, tile, pos, None)
    }

    /// As [`Self::pair`], optionally tagged with a read group.
    fn pair_rg(
        &mut self,
        flowcell: &str,
        lane: u32,
        tile: u32,
        pos: u32,
        rg: Option<&str>,
    ) -> &mut Self {
        self.cluster += 1;
        let name = read_name(flowcell, lane, tile, self.cluster);
        let builder = std::mem::replace(&mut self.builder, SamBuilder::new());
        self.builder = match rg {
            Some(rg) => builder
                .rec_simple_rg(&name, 99, "chr1", pos, "50M", "=", pos + 100, 150, rg)
                .rec_simple_rg(&name, 147, "chr1", pos + 100, "50M", "=", pos, -150, rg),
            None => builder
                .rec_simple(&name, 99, "chr1", pos, "50M", "=", pos + 100, 150)
                .rec_simple(&name, 147, "chr1", pos + 100, "50M", "=", pos, -150),
        };
        self
    }

    /// Add a mapped single-end read, which `picard-exact` defers to its temp BAM.
    fn single_end(&mut self, flowcell: &str, lane: u32, tile: u32, pos: u32) -> &mut Self {
        self.cluster += 1;
        let name = read_name(flowcell, lane, tile, self.cluster);
        let builder = std::mem::replace(&mut self.builder, SamBuilder::new());
        self.builder = builder.rec_simple(&name, 0, "chr1", pos, "50M", "*", 0, 0);
        self
    }

    /// Put one non-duplicate pair on each of `tiles` tiles, giving the library a
    /// realistic tile spread.
    ///
    /// The chance correction needs this: with only a handful of tiles in play, a
    /// same-tile duplicate is no evidence of clustering at all, and the corrected
    /// count is correctly zero.
    fn spread_over_tiles(&mut self, tiles: u32) -> &mut Self {
        for tile in 0..tiles {
            self.pair("FC", 1, tile, 1_000_000 + 1_000 * tile);
        }
        self
    }

    fn write_to(&self, path: &Path) {
        self.builder.write_to(path);
    }
}

#[test]
fn duplicates_on_one_tile_are_reported_as_sequencing_duplicates() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    // Four distinct molecules at one locus, all imaged on tile 1101.
    for _ in 0..4 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "3");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "3");
    assert_eq!(row["library_duplicate_pairs"], "0");
}

/// Build an input mixing a same-tile group with a distinct-tile group, run it with
/// `extra`, and return the metrics row. The mix means a level that corrupted the
/// spill would have to corrupt it consistently to still produce these numbers.
fn split_with_flags(extra: &[&str]) -> HashMap<String, String> {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut input = Run::new();
    input.spread_over_tiles(200);
    for _ in 0..4 {
        input.pair("FC", 1, 1101, 500_000);
    }
    for tile in [2101, 2102, 2103] {
        input.pair("FC", 1, tile, 900_000);
    }
    input.write_to(&env.input);
    run_ok(&env.input, &stats, &out, extra);
    parse_row(&stats)
}

/// The metrics the mixed input must report at any compression level.
fn assert_expected_split(row: &HashMap<String, String>) {
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "3");
    assert_eq!(row["library_duplicate_pairs"], "2");
}

// The unit tests cover the round trip inside TileSpiller; these cover the wiring,
// so a level that never reaches the spiller cannot pass unnoticed.

#[test]
fn an_uncompressed_spill_reports_the_expected_split() {
    assert_expected_split(&split_with_flags(&[]));
}

#[test]
fn a_fast_tier_spill_reports_the_expected_split() {
    assert_expected_split(&split_with_flags(&["--tmp-compression-level", "-5"]));
}

#[test]
fn a_compressed_spill_reports_the_expected_split() {
    assert_expected_split(&split_with_flags(&["--tmp-compression-level", "1"]));
}

#[test]
fn the_highest_accepted_level_reports_the_expected_split() {
    assert_expected_split(&split_with_flags(&["--tmp-compression-level", "9"]));
}

#[test]
fn an_uncompressed_run_reports_no_on_disk_size() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut input = Run::new();
    input.spread_over_tiles(200);
    for _ in 0..4 {
        input.pair("FC", 1, 1101, 500_000);
    }
    input.write_to(&env.input);

    let result = run(&env.input, &stats, &out, &[]);
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("spilled tile records"), "expected the spill line: {stderr}");
    assert!(
        !stderr.contains("on disk"),
        "an uncompressed run should not report an on-disk size: {stderr}"
    );
}

#[test]
fn a_compressed_run_reports_the_on_disk_size_and_a_ratio() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut input = Run::new();
    input.spread_over_tiles(200);
    for _ in 0..4 {
        input.pair("FC", 1, 1101, 500_000);
    }
    input.write_to(&env.input);

    let result = run(&env.input, &stats, &out, &["--tmp-compression-level", "1"]);
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("on disk"), "expected an on-disk size: {stderr}");
    assert!(stderr.contains("x)"), "expected a ratio: {stderr}");
}

#[test]
fn a_compression_level_above_the_accepted_range_is_rejected() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut input = Run::new();
    input.pair("FC", 1, 1101, 500_000);
    input.write_to(&env.input);

    let result = run(&env.input, &stats, &out, &["--tmp-compression-level", "22"]);
    assert!(!result.status.success(), "level 22 should be rejected");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("must be between"), "expected a range error, got: {stderr}");
}

#[test]
fn a_compression_level_of_zero_is_rejected_as_ambiguous() {
    // zstd reads 0 as "the default level" while `--compression-level 0` in the same
    // CLI means "stored", so 0 here would silently turn compression on.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut input = Run::new();
    input.pair("FC", 1, 1101, 500_000);
    input.write_to(&env.input);

    let result = run(&env.input, &stats, &out, &["--tmp-compression-level", "0"]);
    assert!(!result.status.success(), "level 0 should be rejected");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("ambiguous"), "expected the ambiguity error, got: {stderr}");
}

/// Both temporary files — the tile spill and the `picard-exact` orphan buffer —
/// are compressed by the one flag, so this is the only path where two zstd
/// streams are live in a single run. A level that reached one but not the other
/// would still pass every other test.
#[test]
fn compressing_both_temp_files_at_once_reports_the_same_numbers() {
    fn run_both(extra: &[&str]) -> HashMap<String, String> {
        let env = TestEnv::new();
        let stats = env._tmp.path().join("stats.tsv");
        let out = env._tmp.path().join("out.bam");
        let mut input = Run::new();
        input.spread_over_tiles(200);
        for _ in 0..4 {
            input.pair("FC", 1, 1101, 500_000);
        }
        for tile in [2101, 2102, 2103] {
            input.pair("FC", 1, tile, 900_000);
        }
        // Single-end reads route through the picard-exact temp BAM; duplicates
        // among them prove that buffer round-tripped rather than merely existed.
        for _ in 0..3 {
            input.single_end("FC", 1, 1101, 700_000);
        }
        input.write_to(&env.input);
        let mut args = vec!["--single-end-strategy", "picard-exact"];
        args.extend_from_slice(extra);
        run_ok(&env.input, &stats, &out, &args);
        parse_row(&stats)
    }

    let plain = run_both(&[]);
    assert_eq!(run_both(&["--tmp-compression-level", "1"]), plain);
    assert_eq!(run_both(&["--tmp-compression-level", "-5"]), plain);
    assert_expected_split(&plain);
    assert_eq!(plain["duplicate_orphans"], "2", "the temp BAM's own duplicates");
}

#[test]
fn duplicates_on_distinct_tiles_are_reported_as_library_duplicates() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for tile in [2101, 2102, 2103, 2104] {
        run.pair("FC", 1, tile, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "3");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "0");
    assert_eq!(row["library_duplicate_pairs"], "3");
}

#[test]
fn two_copies_on_each_of_two_tiles_yields_two_sequencing_and_one_library() {
    // The {A,A,B,B} case: two molecules, each duplicated once on its own tile.
    // A rule that star-pairs every duplicate to one original scores this 1; a
    // rule of "has any same-tile partner" scores it 3. The answer is 2.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for tile in [2101, 2102] {
        for _ in 0..2 {
            run.pair("FC", 1, tile, 500_000);
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "3");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "2");
    assert_eq!(row["library_duplicate_pairs"], "1");
}

#[test]
fn the_same_tile_number_on_two_flowcells_is_not_a_sequencing_duplicate() {
    // The samtools#1996 failure mode, measured at 1.98% of duplicates on real
    // data: matching tile *numbers* across flowcells are different places.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for flowcell in ["H72CFDSXF", "22T3L2LT4"] {
        run.pair(flowcell, 2, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "1");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "0");
    assert_eq!(row["library_duplicate_pairs"], "1");
}

#[test]
fn the_same_tile_number_in_two_lanes_is_not_a_sequencing_duplicate() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for lane in [3, 4] {
        run.pair("FC", lane, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "1");
    assert_eq!(row["library_duplicate_pairs"], "1");
}

#[test]
fn sequencing_and_library_duplicates_sum_to_the_duplicate_pair_total() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    // A same-tile triple, a cross-tile pair, and a lopsided larger group.
    for _ in 0..3 {
        run.pair("FC", 1, 3101, 400_000);
    }
    for tile in [3102, 3103] {
        run.pair("FC", 1, tile, 500_000);
    }
    for tile in [3104, 3104, 3104, 3105] {
        run.pair("FC", 1, tile, 600_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    let total: u64 = row["duplicate_pairs"].parse().unwrap();
    let sequencing: u64 = row["corrected_sequencing_duplicate_pairs"].parse().unwrap();
    let library: u64 = row["library_duplicate_pairs"].parse().unwrap();
    assert_eq!(total, 6, "3 + 1 + 2 duplicates across the three groups");
    assert_eq!(sequencing + library, total);
}

#[test]
fn a_single_tile_library_blanks_the_split() {
    // Everything the run saw is on one tile, so a same-tile duplicate is exactly
    // what coincidence predicts and no attribution is possible. The counts must be
    // blank rather than zero or total — a measured zero would be a claim.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    for _ in 0..3 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "2", "duplicates are still marked and counted");
    assert_eq!(row["raw_sequencing_duplicate_pairs"], "");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "");
    assert_eq!(row["library_duplicate_pairs"], "");
    assert_eq!(row["frac_sequencing_duplicate_pairs"], "");
}

#[test]
fn both_pair_fractions_share_a_denominator_so_one_bounds_the_other() {
    // `frac_sequencing_duplicates` is over `mapped_pairs`, not over
    // `duplicate_pairs`, so it is directly comparable with `frac_duplicate_pairs`
    // and can never exceed it.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for group in 0..10 {
        for _ in 0..3 {
            run.pair("FC", 1, 1101, 500_000 + 1_000 * group);
        }
    }
    for group in 0..10 {
        for tile in [2101, 2102] {
            run.pair("FC", 1, tile, 700_000 + 1_000 * group);
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    let mapped: f64 = row["mapped_pairs"].parse().unwrap();
    let dups: f64 = row["duplicate_pairs"].parse().unwrap();
    let sequencing: f64 = row["corrected_sequencing_duplicate_pairs"].parse().unwrap();
    let frac_dups: f64 = row["frac_duplicate_pairs"].parse().unwrap();
    let frac_seq: f64 = row["frac_sequencing_duplicate_pairs"].parse().unwrap();

    assert!((frac_dups - dups / mapped).abs() < 1e-6, "frac_duplicate_pairs over mapped_pairs");
    assert!((frac_seq - sequencing / mapped).abs() < 1e-6, "frac_sequencing over mapped_pairs");
    assert!(frac_seq > 0.0 && frac_seq <= frac_dups, "{frac_seq} must not exceed {frac_dups}");
}

#[test]
fn removing_sequencing_duplicates_raises_the_library_size_estimate() {
    // Flowcell duplicates say nothing about how many distinct molecules the
    // library held, so counting them as saturation makes it look smaller. With the
    // split off the estimate falls back to the uncorrected form, so these two runs
    // isolate exactly that effect.
    //
    // The data mixes same-tile and cross-tile groups on purpose: a library whose
    // duplicates are *entirely* sequencing duplicates has shown no evidence of
    // resampling at all, so its corrected size is legitimately not estimable.
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for group in 0..20 {
        for _ in 0..3 {
            run.pair("FC", 1, 1101, 500_000 + 1_000 * group);
        }
    }
    for group in 0..20 {
        for tile in [2101, 2102] {
            run.pair("FC", 1, tile, 700_000 + 1_000 * group);
        }
    }
    run.write_to(&env.input);

    let corrected_stats = env._tmp.path().join("on.tsv");
    run_ok(&env.input, &corrected_stats, &out, &[]);
    let plain_stats = env._tmp.path().join("off.tsv");
    run_ok(&env.input, &plain_stats, &out, &["--no-sequencing-dups"]);

    let corrected: u64 = parse_row(&corrected_stats)["estimated_library_size"].parse().unwrap();
    let plain: u64 = parse_row(&plain_stats)["estimated_library_size"].parse().unwrap();
    assert!(
        corrected > plain,
        "corrected {corrected} should exceed uncorrected {plain} once flowcell duplicates \
         stop counting as saturation"
    );
}

#[test]
fn the_library_size_estimate_is_blank_when_every_duplicate_is_a_sequencing_duplicate() {
    // The degenerate case of the correction: removing the sequencing duplicates
    // leaves the observed total equal to the unique count, so there is no
    // resampling for Lander-Waterman to work from and no size to report.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for group in 0..20 {
        for _ in 0..3 {
            run.pair("FC", 1, 1101, 500_000 + 1_000 * group);
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(
        row["corrected_sequencing_duplicate_pairs"], row["duplicate_pairs"],
        "every duplicate here is a flowcell duplicate"
    );
    assert_eq!(row["estimated_library_size"], "", "nothing left to estimate from");
}

#[test]
fn the_library_size_estimate_falls_back_when_the_split_is_not_estimable() {
    // One tile means the split cannot be computed, but the library-size estimate
    // still has everything it needs — it just cannot subtract a correction. It must
    // fall back to the uncorrected value rather than going blank.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    for _ in 0..3 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "", "split not estimable");
    assert!(!row["estimated_library_size"].is_empty(), "but the estimate still stands");
}

#[test]
fn a_read_name_the_chosen_format_cannot_parse_is_a_hard_error() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("SRR1234567.1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("SRR1234567.1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = run(&env.input, &stats, &out, &[]);

    assert!(!result.status.success(), "an unparseable read name must not be skipped silently");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("SRR1234567.1"), "error should quote the read name: {stderr}");
}

#[test]
fn an_unknown_read_name_format_is_rejected_before_any_work() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--read-name-format", "novaseq"])
        .output()
        .expect("dupblaster ran");

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("illumina"), "should list the valid formats: {stderr}");
}

#[test]
fn a_custom_regex_format_extracts_the_unit_and_tile() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    // Read names in nobody's standard layout, reachable only via a pattern.
    let mut builder = SamBuilder::new().sq("chr1", 1_000_000);
    for tile in 0..8 {
        let name = format!("run7_lane1_tile{tile}_x{tile}_y0");
        let pos = 10_000 + 1_000 * tile;
        builder = builder
            .rec_simple(&name, 99, "chr1", pos, "50M", "=", pos + 100, 150)
            .rec_simple(&name, 147, "chr1", pos + 100, "50M", "=", pos, -150);
    }
    builder.write_to(&env.input);

    let result = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--stats"])
        .arg(&stats)
        .args(["--read-name-format", r"regex:lane(?<su>\d+)_tile(?<tile>\d+)"])
        .output()
        .expect("dupblaster ran");
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));

    let row = parse_row(&stats);
    assert_eq!(
        row["corrected_sequencing_duplicate_pairs"], "0",
        "eight distinct tiles, no clustering"
    );
}

#[test]
fn a_regex_without_the_required_capture_groups_is_rejected() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--read-name-format", r"regex:^(?<su>\d+)"])
        .output()
        .expect("dupblaster ran");

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("(?<tile>"), "should name the missing group: {stderr}");
}

#[test]
fn the_per_sequencing_unit_table_is_written_beside_the_stats_file() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    // Two flowcells, two lanes each, so four distinct sequencing units.
    let mut pos = 10_000;
    for (flowcell, lane) in [("H72CFDSXF", 1), ("H72CFDSXF", 2), ("22T3L2LT4", 1), ("22T3L2LT4", 2)]
    {
        for tile in 0..4 {
            run.pair(flowcell, lane, tile, pos);
            pos += 1_000;
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let rows = parse_rows(&env._tmp.path().join("stats.sequencing-units.tsv"));
    let names: Vec<&str> = rows.iter().map(|r| r["sequencing_unit"].as_str()).collect();
    assert_eq!(
        names,
        ["22T3L2LT4:1", "22T3L2LT4:2", "H72CFDSXF:1", "H72CFDSXF:2"],
        "units are named from the read names and ordered deterministically"
    );
    for row in &rows {
        assert_eq!(row["templates"], "4");
        assert_eq!(row["tiles"], "4");
    }
}

#[test]
fn the_per_sequencing_unit_table_attributes_sequencing_duplicates_to_its_flowcell() {
    // One flowcell is loaded too densely, the other is not. The per-unit table is
    // what makes that visible — a per-library average hides it entirely, and the
    // spread across flowcells of one real sample was 9x.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    let mut pos = 10_000;
    for flowcell in ["DENSE", "SPARSE"] {
        for tile in 0..100 {
            run.pair(flowcell, 1, tile, pos);
            pos += 1_000;
        }
    }
    // Only DENSE has same-tile duplicate groups.
    for group in 0..5 {
        for _ in 0..3 {
            run.pair("DENSE", 1, 7, 5_000_000 + 1_000 * group);
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let rows = parse_rows(&env._tmp.path().join("stats.sequencing-units.tsv"));
    let by_unit: HashMap<&str, &HashMap<String, String>> =
        rows.iter().map(|r| (r["sequencing_unit"].as_str(), r)).collect();
    let dense: u64 = by_unit["DENSE:1"]["sequencing_duplicate_pairs"].parse().unwrap();
    let sparse: u64 = by_unit["SPARSE:1"]["sequencing_duplicate_pairs"].parse().unwrap();
    assert_eq!(dense, 10, "five groups of three, two sequencing duplicates each");
    assert_eq!(sparse, 0, "the sparse flowcell has no same-tile groups");
}

#[test]
fn the_per_unit_column_sums_exactly_to_the_raw_per_library_figure() {
    // Each tile contributes `members - 1`, so both figures are one integer sum
    // viewed two ways: no rounding, no attribution heuristic, no tolerance.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    let mut pos = 10_000;
    for flowcell in ["FCA", "FCB"] {
        for tile in 0..6 {
            run.pair(flowcell, 1, tile, pos);
            pos += 1_000;
        }
    }
    for flowcell in ["FCA", "FCB"] {
        for group in 0..10 {
            for _ in 0..3 {
                run.pair(flowcell, 1, 2, 5_000_000 + 1_000 * group);
            }
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let raw: u64 = parse_row(&stats)["raw_sequencing_duplicate_pairs"].parse().unwrap();
    let per_unit: u64 = parse_rows(&env._tmp.path().join("stats.sequencing-units.tsv"))
        .iter()
        .map(|r| r["sequencing_duplicate_pairs"].parse::<u64>().unwrap())
        .sum();
    assert!(raw > 0, "the test data should produce some sequencing duplicates");
    assert_eq!(per_unit, raw, "per-unit must sum exactly, not approximately");
}

#[test]
fn the_per_unit_column_reconciles_across_many_heterogeneous_units() {
    // Twenty units with differing tile spreads and group shapes, so per-unit
    // residues cannot cancel by symmetry. Exactness must not depend on scale.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    let mut pos = 10_000;
    let flowcells: Vec<String> = (0..20).map(|i| format!("FC{i}")).collect();
    for (index, flowcell) in flowcells.iter().enumerate() {
        for tile in 0..(4 + index as u32 % 9) {
            run.pair(flowcell, 1, tile, pos);
            pos += 1_000;
        }
    }
    for (index, flowcell) in flowcells.iter().enumerate() {
        for group in 0..(1 + index as u32 % 5) {
            for _ in 0..(2 + index as u32 % 5) {
                run.pair(flowcell, 1, 0, 5_000_000 + 1_000 * group);
            }
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let raw: u64 = parse_row(&stats)["raw_sequencing_duplicate_pairs"].parse().unwrap();
    let per_unit: u64 = parse_rows(&env._tmp.path().join("stats.sequencing-units.tsv"))
        .iter()
        .map(|r| r["sequencing_duplicate_pairs"].parse::<u64>().unwrap())
        .sum();
    assert!(raw > 0);
    assert_eq!(per_unit, raw, "20 units must still reconcile exactly");
}

/// Run without `--read-name-format`, so the layout defaults to Illumina.
fn run_default(input: &Path, stats: &Path, output: &Path, extra: &[&str]) -> std::process::Output {
    Command::new(rust_binary())
        .args(["-i"])
        .arg(input)
        .args(["-o"])
        .arg(output)
        .args(["--stats"])
        .arg(stats)
        .args(extra)
        .output()
        .expect("dupblaster ran")
}

#[test]
fn the_split_runs_by_default_with_no_flags_at_all() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for _ in 0..4 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    let result = run_default(&env.input, &stats, &out, &[]);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));

    let row = parse_row(&stats);
    assert_eq!(
        row["corrected_sequencing_duplicate_pairs"], "3",
        "the split should be on by default"
    );
    assert_eq!(row["library_duplicate_pairs"], "0");
    assert!(env._tmp.path().join("stats.sequencing-units.tsv").exists());
}

#[test]
fn no_sequencing_dups_turns_the_split_off() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(200);
    for _ in 0..4 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    let result = run_default(&env.input, &stats, &out, &["--no-sequencing-dups"]);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));

    let row = parse_row(&stats);
    assert_eq!(row["duplicate_pairs"], "3", "duplicates are still marked and counted");
    assert_eq!(row["corrected_sequencing_duplicate_pairs"], "", "but not classified");
    assert_eq!(row["library_duplicate_pairs"], "");
    assert!(!env._tmp.path().join("stats.sequencing-units.tsv").exists());
}

#[test]
fn unparseable_read_names_fail_the_run_even_with_no_flags() {
    // The split is on by default, so this is what a user of a platform without a
    // preset meets first. It must fail rather than quietly drop the metric, and the
    // message has to name both ways forward.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("SRR1234567.1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("SRR1234567.1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = run_default(&env.input, &stats, &out, &[]);

    assert!(!result.status.success(), "must not silently skip the metric");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("SRR1234567.1"), "should quote the name: {stderr}");
    assert!(stderr.contains("--read-name-format"), "should offer the layout flag: {stderr}");
    assert!(stderr.contains("--no-sequencing-dups"), "should offer the opt-out: {stderr}");
}

#[test]
fn the_split_can_be_skipped_on_a_platform_without_a_preset() {
    // The escape hatch the failure above points at.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("F350009384L1C001R0010000000", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("F350009384L1C001R0010000000", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = run_default(&env.input, &stats, &out, &["--no-sequencing-dups"]);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    assert_eq!(parse_row(&stats)["duplicate_pairs"], "0");
}

#[test]
fn naming_the_format_explicitly_makes_an_unparseable_name_fatal() {
    // Same input as above, but the user asserted the layout — so silently
    // dropping the metric would hide their mistake.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("SRR1234567.1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("SRR1234567.1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let result = run(&env.input, &stats, &out, &[]);
    assert!(!result.status.success(), "an explicitly named format must fail loudly");
}

#[test]
fn read_names_that_stop_parsing_partway_are_fatal_even_by_default() {
    // Parseable names then unparseable ones means two platforms' data got merged.
    // Reporting a split computed from the parseable prefix would silently describe
    // part of the file, so this fails however the format was chosen.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(20);
    run.builder = std::mem::replace(&mut run.builder, SamBuilder::new())
        .rec_simple("SRR1234567.1", 99, "chr1", 900_000, "50M", "=", 900_100, 150)
        .rec_simple("SRR1234567.1", 147, "chr1", 900_100, "50M", "=", 900_000, -150);
    run.write_to(&env.input);
    let result = run_default(&env.input, &stats, &out, &[]);
    assert!(!result.status.success(), "a mid-file format change must not pass silently");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("SRR1234567"), "{stderr}");
}

#[test]
fn the_split_is_unchanged_under_every_single_end_strategy() {
    // Only pairs are spilled, so the orphan-keying strategy must not touch the
    // split. `picard-exact` is the one that matters: it runs a second pass over
    // buffered orphan blocks, and a pair counted in both passes would break the
    // seq + lib == duplicate_pairs identity.
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(100);
    for _ in 0..3 {
        run.pair("FC", 1, 1101, 500_000);
    }
    // A mapped orphan, so the deferred-fragment path actually has work to do.
    run.builder = std::mem::replace(&mut run.builder, SamBuilder::new()).rec_simple(
        &read_name("FC", 1, 1101, 90_001),
        73,
        "chr1",
        700_000,
        "50M",
        "*",
        0,
        0,
    );
    run.write_to(&env.input);

    let mut seen: Vec<String> = Vec::new();
    for strategy in ["strand-aware", "picard-approx", "picard-exact", "samblaster-legacy"] {
        let stats = env._tmp.path().join(format!("stats-{strategy}.tsv"));
        run_ok(&env.input, &stats, &out, &["--single-end-strategy", strategy, "--ignore-unmated"]);
        let row = parse_row(&stats);
        let total: u64 = row["duplicate_pairs"].parse().unwrap();
        let sequencing: u64 = row["corrected_sequencing_duplicate_pairs"].parse().unwrap();
        let library: u64 = row["library_duplicate_pairs"].parse().unwrap();
        assert_eq!(sequencing + library, total, "identity broken under {strategy}");
        seen.push(format!("{sequencing}/{library}/{total}"));
    }
    assert!(
        seen.windows(2).all(|w| w[0] == w[1]),
        "orphan strategy changed the pair split: {seen:?}"
    );
}

#[test]
fn the_output_bam_is_complete_and_readable_with_the_decomposition_on() {
    // The decomposition runs a post-pass after the output stream is closed; that
    // must not truncate or disturb the BAM itself.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(20);
    for _ in 0..3 {
        run.pair("FC", 1, 1101, 500_000);
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let records = read_records(&out);
    assert_eq!(records.len(), 46, "20 spread pairs plus 3 at one locus, two records each");
}

#[test]
fn the_decomposition_is_reported_per_library() {
    // Two libraries on one flowcell: lib1's duplicates are same-tile and lib2's
    // are cross-tile, so the rows must not be averaged together.
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.builder = std::mem::replace(&mut run.builder, SamBuilder::new())
        .rg("rg1", "sample", Some("lib1"))
        .rg("rg2", "sample", Some("lib2"));
    let mut pos = 10_000;
    for (rg, dup_tiles) in [("rg1", [1101, 1101, 1101]), ("rg2", [2101, 2102, 2103])] {
        for tile in 0..50 {
            run.pair_rg("FC", 1, 5_000 + tile, pos, Some(rg));
            pos += 1_000;
        }
        for tile in dup_tiles {
            run.pair_rg("FC", 1, tile, 5_000_000, Some(rg));
        }
    }
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let rows = parse_rows(&stats);
    let by_library: HashMap<&str, &HashMap<String, String>> =
        rows.iter().map(|r| (r["library"].as_str(), r)).collect();
    assert_eq!(by_library["lib1"]["corrected_sequencing_duplicate_pairs"], "2");
    assert_eq!(by_library["lib1"]["library_duplicate_pairs"], "0");
    assert_eq!(by_library["lib2"]["corrected_sequencing_duplicate_pairs"], "0");
    assert_eq!(by_library["lib2"]["library_duplicate_pairs"], "2");
}
