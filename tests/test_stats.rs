//! Integration tests for the `--stats` TSV output.

mod helpers;

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use helpers::*;

fn run_with_stats(input: &Path, stats: &Path, output: &Path, extra: &[&str]) {
    let mut cmd = Command::new(rust_binary());
    cmd.args(["-i"]).arg(input).args(["-o"]).arg(output).args(["--stats"]).arg(stats).args(extra);
    let out = cmd.output().expect("rust dupblaster ran");
    assert!(
        out.status.success(),
        "rust dupblaster failed: {}",
        String::from_utf8_lossy(&out.stderr),
    );
}

/// Parse the two-line TSV into a column → value map.
fn parse_stats_tsv(path: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(path).expect("read stats");
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let values = lines.next().expect("value line");
    let cols: Vec<&str> = header.split('\t').collect();
    let vals: Vec<&str> = values.split('\t').collect();
    assert_eq!(cols.len(), vals.len(), "column count mismatch in stats TSV");
    cols.into_iter().map(String::from).zip(vals.into_iter().map(String::from)).collect()
}

#[test]
fn stats_tsv_reports_correct_counts_for_simple_dup_input() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    // 3 mapped pairs, one of which is a duplicate of another. Expected:
    //   total_templates = 3, mapped_pairs = 3, duplicate_pairs = 1.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r3", 99, "chr1", 500, "50M", "=", 600, 150)
        .rec_simple("r3", 147, "chr1", 600, "50M", "=", 500, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &[]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["total_templates"], "3");
    assert_eq!(m["mapped_pairs"], "3");
    assert_eq!(m["duplicate_pairs"], "1");
    assert_eq!(m["duplicate_templates"], "1");
    // Picard formula: (0 + 2*1) / (0 + 2*3) = 2/6 = 0.3333...
    // Tolerance matches the `{:.6}` render precision.
    let frac: f64 = m["frac_duplicates"].parse().unwrap();
    assert!((frac - 1.0 / 3.0).abs() < 1e-6, "frac_duplicates was {frac}");
    assert!(!m["estimated_library_size"].is_empty(), "library size should be set");
    assert_eq!(m["dupblaster_version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(m["sample"], "");
}

#[test]
fn stats_tsv_uses_sample_override() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &["--sample", "NA12878"]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["sample"], "NA12878");
}

#[test]
fn stats_tsv_pulls_sample_from_read_group_sm() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("rg1", "SAMPLE_A", Some("libA"))
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &[]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["sample"], "SAMPLE_A");
}

#[test]
fn stats_tsv_comma_joins_multiple_sm_values() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("rg1", "SAMPLE_A", None)
        .rg("rg2", "SAMPLE_B", None)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &[]);
    let m = parse_stats_tsv(&stats);
    // BTreeSet sort order: alphabetical.
    assert_eq!(m["sample"], "SAMPLE_A,SAMPLE_B");
}

/// Regression for the unmapped-vs-unmated stats accounting bug. A block
/// containing only one primary record that is both paired and unmapped
/// should land in `unmapped_orphans`, NOT `unmated_templates`.
#[test]
fn paired_unmapped_singleton_counts_as_unmapped_orphan_not_unmated() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    // FLAG 5 = paired (1) + unmapped (4). Single primary in its qname block.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // r1: a normal mapped pair so total_templates > 0.
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        // r2: paired + first + unmapped, only one record in the block.
        // FLAG 69 = 1 (paired) + 4 (unmapped) + 64 (first).
        .rec_simple("r2", 69, "*", 0, "*", "*", 0, 0)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &["--ignore-unmated"]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["unmapped_orphans"], "1", "row: {m:?}");
    assert_eq!(m["unmated_templates"], "0", "row: {m:?}");
}

#[test]
fn stats_tsv_has_all_expected_columns() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &[]);
    let text = std::fs::read_to_string(&stats).unwrap();
    let header = text.lines().next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let expected = [
        "sample",
        "library",
        "dupblaster_version",
        "total_templates",
        "duplicate_templates",
        "frac_duplicates",
        "mapped_pairs",
        "duplicate_pairs",
        "mapped_orphans",
        "duplicate_orphans",
        "unmapped_orphans",
        "unmapped_pairs",
        "unmated_templates",
        "estimated_library_size",
    ];
    assert_eq!(cols, expected);
}

/// `--remove-dups` drops duplicate records from the output but the `--stats`
/// counts still reflect every template examined (Picard semantics): the
/// duplicate is counted, not hidden.
#[test]
fn stats_counts_unaffected_by_remove_dups() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    run_with_stats(&env.input, &stats, &out, &["--remove-dups"]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["total_templates"], "2");
    assert_eq!(m["duplicate_pairs"], "1", "the duplicate is still counted");
    assert_eq!(m["duplicate_templates"], "1");
    // The duplicate pair's two records are removed; only the kept pair's
    // two records remain in the output.
    let records = read_records(&out);
    assert_eq!(records.len(), 2, "duplicate pair should be removed from output");
}

/// Empty input (header only, zero records) must not crash the `--stats`
/// writer. Risk: `frac_duplicates`'s read-level fraction has 0 in the
/// denominator when no reads were observed; that needs to land as `0.0`
/// (or empty) rather than NaN / panic. `estimated_library_size` should
/// be empty since there are no pairs to estimate from.
#[test]
fn stats_tsv_handles_empty_input_without_crashing() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new().sq("chr1", 1_000_000).write_to(&env.input); // header, no records
    run_with_stats(&env.input, &stats, &out, &[]);
    let m = parse_stats_tsv(&stats);
    assert_eq!(m["total_templates"], "0");
    assert_eq!(m["duplicate_templates"], "0");
    assert_eq!(m["mapped_pairs"], "0");
    assert_eq!(m["duplicate_pairs"], "0");
    // frac_duplicates: well-formed number (not NaN, not inf).
    let f: f64 = m["frac_duplicates"].parse().expect("frac_duplicates is a number");
    assert!(f.is_finite(), "frac_duplicates was non-finite: {}", m["frac_duplicates"]);
    assert_eq!(f, 0.0, "frac_duplicates should be 0 with no reads");
    // estimated_library_size: empty cell since there are no pairs.
    assert_eq!(m["estimated_library_size"], "");
}
