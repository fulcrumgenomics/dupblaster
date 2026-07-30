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
    assert_eq!(row["sequencing_duplicates"], "3");
    assert_eq!(row["library_duplicates"], "0");
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
    assert_eq!(row["sequencing_duplicates"], "0");
    assert_eq!(row["library_duplicates"], "3");
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
    assert_eq!(row["sequencing_duplicates"], "2");
    assert_eq!(row["library_duplicates"], "1");
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
    assert_eq!(row["sequencing_duplicates"], "0");
    assert_eq!(row["library_duplicates"], "1");
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
    assert_eq!(row["library_duplicates"], "1");
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
    let sequencing: u64 = row["sequencing_duplicates"].parse().unwrap();
    let library: u64 = row["library_duplicates"].parse().unwrap();
    assert_eq!(total, 6, "3 + 1 + 2 duplicates across the three groups");
    assert_eq!(sequencing + library, total);
}

#[test]
fn a_single_tile_library_blanks_the_split_but_still_reports_the_collision_rate() {
    // Everything on one tile means q == 1 and no duplicate can be attributed:
    // a same-tile duplicate is exactly what coincidence predicts. The counts
    // must be blank rather than zero or total, with the tile evidence present so
    // the reason is visible instead of looking like a bug.
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
    assert_eq!(row["tile_count"], "1");
    assert_eq!(row["tile_collision_rate"], "1.000000");
    assert_eq!(row["sequencing_duplicates"], "");
    assert_eq!(row["library_duplicates"], "");
}

#[test]
fn the_tile_collision_rate_is_reported_for_an_ordinary_library() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(4);
    run.write_to(&env.input);
    run_ok(&env.input, &stats, &out, &[]);

    let row = parse_row(&stats);
    assert_eq!(row["tile_count"], "4");
    assert_eq!(row["tile_collision_rate"], "0.250000", "four evenly-used tiles give q = 1/4");
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
    assert_eq!(row["tile_count"], "8", "each of the eight tiles should be distinct");
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
    let dense: u64 = by_unit["DENSE:1"]["sequencing_duplicates"].parse().unwrap();
    let sparse: u64 = by_unit["SPARSE:1"]["sequencing_duplicates"].parse().unwrap();
    assert_eq!(dense, 10, "five groups of three, two sequencing duplicates each");
    assert_eq!(sparse, 0, "the sparse flowcell has no same-tile groups");
}

#[test]
fn no_sequencing_unit_table_is_written_without_the_read_name_format() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
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
        .args(["--stats"])
        .arg(&stats)
        .output()
        .expect("dupblaster ran");
    assert!(result.status.success());
    assert!(!env._tmp.path().join("stats.sequencing-units.tsv").exists());
}

#[test]
fn the_split_is_identical_under_every_spill_bucket_count() {
    // Bucketing exists only to bound descriptor use; it must not change which
    // records meet each other in a group.
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    let mut run = Run::new();
    run.spread_over_tiles(50);
    for group in 0..10 {
        for tile in [4101, 4101, 4102] {
            run.pair("FC", 1, tile, 5_000_000 + 1_000 * group);
        }
    }
    run.write_to(&env.input);

    let mut seen: Vec<String> = Vec::new();
    for buckets in ["1", "3", "64"] {
        let stats = env._tmp.path().join(format!("stats-{buckets}.tsv"));
        run_ok(&env.input, &stats, &out, &["--spill-buckets", buckets]);
        let row = parse_row(&stats);
        let total: u64 = row["duplicate_pairs"].parse().unwrap();
        let sequencing: u64 = row["sequencing_duplicates"].parse().unwrap();
        let library: u64 = row["library_duplicates"].parse().unwrap();
        assert_eq!(total, 20, "ten groups of three contribute two duplicates each");
        assert_eq!(sequencing + library, total);
        seen.push(format!("{sequencing}/{library}"));
    }
    assert!(seen.windows(2).all(|w| w[0] == w[1]), "bucket count changed the answer: {seen:?}");
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
    assert_eq!(by_library["lib1"]["sequencing_duplicates"], "2");
    assert_eq!(by_library["lib1"]["library_duplicates"], "0");
    assert_eq!(by_library["lib2"]["sequencing_duplicates"], "0");
    assert_eq!(by_library["lib2"]["library_duplicates"], "2");
}
