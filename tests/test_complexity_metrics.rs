//! Integration tests for the complexity metrics (duplicate-rate ladder +
//! group-size histogram), including coverage across all single-end strategies.

mod helpers;

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use helpers::*;

/// Run dupblaster with complexity metrics left at their default (on) and return
/// the parsed rows of the file with the given suffix.
///
/// Builds its own command rather than using [`dupblaster`], which switches
/// complexity metrics off — this file is the one that must exercise the default.
/// The sequencing split is off because these inputs use short synthetic read names
/// carrying no flowcell or tile.
fn run_and_read(
    input: &Path,
    prefix: &Path,
    suffix: &str,
    extra: &[&str],
) -> Vec<HashMap<String, String>> {
    let out = prefix.with_extension("out.bam");
    let status = Command::new(rust_binary())
        .args(["-i"])
        .arg(input)
        .args(["-o"])
        .arg(&out)
        .args(["--metrics-prefix"])
        .arg(prefix)
        .args(["--sequencing-dups", "off"])
        .args(extra)
        .output()
        .expect("dupblaster ran");
    assert!(
        status.status.success(),
        "dupblaster failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    let mut path = prefix.as_os_str().to_owned();
    path.push(suffix);
    parse_rows(Path::new(&path))
}

fn run_ladder(input: &Path, prefix: &Path, extra: &[&str]) -> Vec<HashMap<String, String>> {
    run_and_read(input, prefix, ".duplication-sampled.tsv", extra)
}

fn run_histogram(input: &Path, prefix: &Path, extra: &[&str]) -> Vec<HashMap<String, String>> {
    run_and_read(input, prefix, ".duplication-spectrum.tsv", extra)
}

/// Parse a multi-row TSV into a vector of column→value maps.
fn parse_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = text.lines();
    let cols: Vec<String> = lines.next().expect("header").split('\t').map(String::from).collect();
    lines
        .map(|line| {
            let vals: Vec<&str> = line.split('\t').collect();
            assert_eq!(vals.len(), cols.len(), "ragged row: {line:?}");
            cols.iter().cloned().zip(vals.iter().map(|s| s.to_string())).collect()
        })
        .collect()
}

/// Find the row for a given (library, category, total).
fn find<'a>(
    rows: &'a [HashMap<String, String>],
    library: &str,
    category: &str,
    total: &str,
) -> &'a HashMap<String, String> {
    rows.iter()
        .find(|r| r["library"] == library && r["category"] == category && r["total"] == total)
        .unwrap_or_else(|| panic!("no {library}/{category}/total={total} row"))
}

/// `n_molecules` for a (category, n_observations); 0 if absent.
fn n_sigs(rows: &[HashMap<String, String>], category: &str, k: &str) -> u64 {
    rows.iter()
        .find(|r| r["category"] == category && r["n_observations"] == k)
        .map(|r| r["n_molecules"].parse().unwrap())
        .unwrap_or(0)
}

/// A PE input: one pair coordinate observed 3× (a triple) plus two unique pairs.
fn write_pe_triple_plus_two_singletons(path: &Path) {
    SamBuilder::new()
        .sq("chr1", 2_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 300, 250)
        .rec_simple("r1", 147, "chr1", 300, "50M", "=", 100, -250)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 300, 250)
        .rec_simple("r2", 147, "chr1", 300, "50M", "=", 100, -250)
        .rec_simple("r3", 99, "chr1", 100, "50M", "=", 300, 250)
        .rec_simple("r3", 147, "chr1", 300, "50M", "=", 100, -250)
        .rec_simple("r4", 99, "chr1", 500, "50M", "=", 700, 250)
        .rec_simple("r4", 147, "chr1", 700, "50M", "=", 500, -250)
        .rec_simple("r5", 99, "chr1", 900, "50M", "=", 1100, 250)
        .rec_simple("r5", 147, "chr1", 1100, "50M", "=", 900, -250)
        .write_to(path);
}

/// Three single-end reads: two share a 5' position, one unique.
fn write_se_two_plus_one(path: &Path) {
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .record("s1", 0, "chr1", 100, 60, "50M", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .record("s2", 0, "chr1", 100, 60, "50M", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .record("s3", 0, "chr1", 500, 60, "50M", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(path);
}

/// One library holding both both-mapped pairs and mapped orphans, the orphans
/// arriving both *before* and *after* the first pair. `r1`/`r2` share a pair
/// coordinate (one duplicate pair); `r3` is a unique pair. `s1`/`s2` are duplicate
/// orphans seen before any pair; `s3` is an orphan after it. All one (default)
/// library, so it is reported on `pairs` and the orphans must never surface.
fn write_mixed_pairs_and_orphans(path: &Path) {
    let seq = "A".repeat(50);
    let qual = "I".repeat(50);
    SamBuilder::new()
        .sq("chr1", 2_000_000)
        // Orphans before any pair (s1, s2 share a 5' coord → duplicate orphans).
        .record("s1", 0, "chr1", 100, 60, "50M", "*", 0, 0, &seq, &qual)
        .record("s2", 0, "chr1", 100, 60, "50M", "*", 0, 0, &seq, &qual)
        // First pair — flips `has_pairs`, dropping the orphan counts above.
        .rec_simple("r1", 99, "chr1", 500, "50M", "=", 700, 250)
        .rec_simple("r1", 147, "chr1", 700, "50M", "=", 500, -250)
        // Duplicate pair (same coordinates as r1).
        .rec_simple("r2", 99, "chr1", 500, "50M", "=", 700, 250)
        .rec_simple("r2", 147, "chr1", 700, "50M", "=", 500, -250)
        // Unique pair.
        .rec_simple("r3", 99, "chr1", 900, "50M", "=", 1100, 250)
        .rec_simple("r3", 147, "chr1", 1100, "50M", "=", 900, -250)
        // Orphan after the first pair (must be ignored by the counter).
        .record("s3", 0, "chr1", 1500, 60, "50M", "*", 0, 0, &seq, &qual)
        .write_to(path);
}

/// Two read groups → two libraries (libA, libB), each with a duplicate pair, so
/// the histogram spans multiple libraries and the plot must facet into one file.
fn write_two_library_pairs(path: &Path) {
    SamBuilder::new()
        .sq("chr1", 2_000_000)
        .rg("rg1", "S", Some("libA"))
        .rg("rg2", "S", Some("libB"))
        // libA: a duplicate pair (a1/a2) plus a unique pair (a3).
        .rec_simple_rg("a1", 99, "chr1", 100, "50M", "=", 300, 250, "rg1")
        .rec_simple_rg("a1", 147, "chr1", 300, "50M", "=", 100, -250, "rg1")
        .rec_simple_rg("a2", 99, "chr1", 100, "50M", "=", 300, 250, "rg1")
        .rec_simple_rg("a2", 147, "chr1", 300, "50M", "=", 100, -250, "rg1")
        .rec_simple_rg("a3", 99, "chr1", 500, "50M", "=", 700, 250, "rg1")
        .rec_simple_rg("a3", 147, "chr1", 700, "50M", "=", 500, -250, "rg1")
        // libB: a duplicate pair (b1/b2).
        .rec_simple_rg("b1", 99, "chr1", 900, "50M", "=", 1100, 250, "rg2")
        .rec_simple_rg("b1", 147, "chr1", 1100, "50M", "=", 900, -250, "rg2")
        .rec_simple_rg("b2", 99, "chr1", 900, "50M", "=", 1100, 250, "rg2")
        .rec_simple_rg("b2", 147, "chr1", 1100, "50M", "=", 900, -250, "rg2")
        .write_to(path);
}

// ── ladder ────────────────────────────────────────────────────────────────

#[test]
fn ladder_reports_pairs_only_for_pe_input_sample_first() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_pe_triple_plus_two_singletons(&env.input);
    let rows = run_ladder(&env.input, &prefix, &["--sample", "NA12878"]);
    assert!(!rows.is_empty());
    assert_eq!(rows[0]["sample"], "NA12878");
    // One category only — pairs — never `all` or `single_end`.
    assert!(rows.iter().all(|r| r["category"] == "pairs"), "PE input → pairs only");
    let lib = rows[0]["library"].clone();
    let final_row = find(&rows, &lib, "pairs", "5");
    assert_eq!(final_row["duplicates"], "2"); // r2, r3
    assert_eq!(final_row["unique"], "3");
}

#[test]
fn ladder_reports_single_end_for_se_only_input() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_se_two_plus_one(&env.input);
    let rows = run_ladder(&env.input, &prefix, &[]);
    assert!(rows.iter().all(|r| r["category"] == "single_end"), "SE-only → single_end");
    let lib = rows[0]["library"].clone();
    let final_row = find(&rows, &lib, "single_end", "3");
    assert_eq!(final_row["duplicates"], "1");
    assert_eq!(final_row["unique"], "2");
}

#[test]
fn ladder_pairs_curve_snapshots_at_intervals() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    let mut b = SamBuilder::new().sq("chr1", 5_000_000);
    for i in 0..5u32 {
        let p = 100 + i * 1000;
        b = b.rec_simple(&format!("r{i}"), 99, "chr1", p, "50M", "=", p + 200, 250).rec_simple(
            &format!("r{i}"),
            147,
            "chr1",
            p + 200,
            "50M",
            "=",
            p,
            -250,
        );
    }
    b.write_to(&env.input);
    let rows = run_ladder(&env.input, &prefix, &["--complexity-interval", "2"]);
    let lib = rows[0]["library"].clone();
    let totals: Vec<&str> = rows
        .iter()
        .filter(|r| r["library"] == lib && r["category"] == "pairs")
        .map(|r| r["total"].as_str())
        .collect();
    assert_eq!(totals, ["2", "4", "5"], "snapshots at 2, 4, and a final at 5");
}

// ── histogram ────────────────────────────────────────────────────────────

#[test]
fn histogram_counts_pair_group_sizes() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_pe_triple_plus_two_singletons(&env.input);
    let rows = run_histogram(&env.input, &prefix, &[]);
    assert_eq!(n_sigs(&rows, "pairs", "1"), 2, "two singletons");
    assert_eq!(n_sigs(&rows, "pairs", "3"), 1, "one triple");
    assert!(rows.iter().all(|r| r["category"] == "pairs"));
}

#[test]
fn histogram_single_end_only_input() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_se_two_plus_one(&env.input);
    let rows = run_histogram(&env.input, &prefix, &["--sample", "RNA1"]);
    assert_eq!(n_sigs(&rows, "single_end", "1"), 1);
    assert_eq!(n_sigs(&rows, "single_end", "2"), 1);
    assert!(rows.iter().all(|r| r["category"] == "single_end"));
    assert!(rows.iter().all(|r| r["sample"] == "RNA1"));
}

#[test]
fn mixed_library_reports_only_pairs_and_ignores_orphans() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_mixed_pairs_and_orphans(&env.input);

    // Ladder: pairs only, ending at the 3 both-mapped templates (r2 the dup).
    let ladder = run_ladder(&env.input, &prefix, &[]);
    assert!(ladder.iter().all(|r| r["category"] == "pairs"), "mixed library → pairs only");
    let lib = ladder[0]["library"].clone();
    let final_row = find(&ladder, &lib, "pairs", "3");
    assert_eq!(final_row["duplicates"], "1"); // r2
    assert_eq!(final_row["unique"], "2"); // r1, r3

    // Histogram: the r1/r2 pair is one signature seen 2×, r3 a singleton; the
    // orphan signatures (s1/s2/s3) must not surface in any category.
    let hist = run_histogram(&env.input, &prefix, &[]);
    assert!(hist.iter().all(|r| r["category"] == "pairs"), "no single_end rows for a paired lib");
    assert_eq!(n_sigs(&hist, "pairs", "1"), 1, "r3 singleton");
    assert_eq!(n_sigs(&hist, "pairs", "2"), 1, "r1/r2 duplicate pair");
    assert_eq!(n_sigs(&hist, "single_end", "1"), 0, "orphans must not surface");
    assert_eq!(n_sigs(&hist, "single_end", "2"), 0, "duplicate orphan must not surface");
}

// ── plots ──────────────────────────────────────────────────────────────────

#[test]
fn writes_pdf_plots_for_both_metrics() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_pe_triple_plus_two_singletons(&env.input);
    // `run_ladder` runs the binary, which writes all four outputs; then check
    // the two PDFs exist and are real PDFs (a single-library run, so the
    // histogram plot is the unqualified `<prefix>.duplication-spectrum.pdf`).
    let _ = run_ladder(&env.input, &prefix, &[]);
    for suffix in [".duplication-sampled.pdf", ".duplication-spectrum.pdf"] {
        let mut p = prefix.as_os_str().to_owned();
        p.push(suffix);
        let bytes =
            std::fs::read(Path::new(&p)).unwrap_or_else(|e| panic!("reading {suffix}: {e}"));
        assert!(bytes.starts_with(b"%PDF"), "{suffix} is not a PDF");
    }
}

#[test]
fn count_histogram_multi_library_writes_one_faceted_pdf() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_two_library_pairs(&env.input);
    // The TSV carries both libraries...
    let hist = run_histogram(&env.input, &prefix, &[]);
    assert!(hist.iter().any(|r| r["library"] == "libA"), "libA present");
    assert!(hist.iter().any(|r| r["library"] == "libB"), "libB present");
    // ...and the plot is a single faceted file, never one per library.
    let path = |suffix: &str| {
        let mut p = prefix.as_os_str().to_owned();
        p.push(suffix);
        std::path::PathBuf::from(p)
    };
    assert!(std::fs::read(path(".duplication-spectrum.pdf")).unwrap().starts_with(b"%PDF"));
    assert!(!path(".duplication-spectrum.libA.pdf").exists(), "no per-library plot file");
    assert!(!path(".duplication-spectrum.libB.pdf").exists(), "no per-library plot file");
}

#[test]
fn unwritable_metrics_dir_fails_early_with_clear_error() {
    // The metrics files are written only at the very end, so a bad prefix must be
    // caught up front rather than after a full processing pass.
    let env = TestEnv::new();
    write_pe_triple_plus_two_singletons(&env.input);
    let out = env._tmp.path().join("out.bam");
    let bad_prefix = env._tmp.path().join("no-such-dir").join("pfx");
    let res = dupblaster()
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--metrics-prefix"])
        .arg(&bad_prefix)
        .output()
        .expect("dupblaster ran");
    assert!(!res.status.success(), "must fail on an unwritable --metrics-prefix directory");
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(
        err.contains("--metrics-prefix") && err.contains("does not exist"),
        "error should name the flag and the missing directory: {err}"
    );
}

// ── coverage across every single-end strategy ──────────────────────────────

/// The pairs histogram is keyed identically (`check_dm`) regardless of the
/// single-end strategy, so it must be byte-identical across all four.
#[test]
fn pairs_histogram_is_identical_across_all_single_end_strategies() {
    let env = TestEnv::new();
    write_pe_triple_plus_two_singletons(&env.input);
    let mut reference: Option<Vec<(String, String)>> = None;
    for strat in ["strand-aware", "picard-approx", "picard-exact", "samblaster-legacy"] {
        let prefix = env._tmp.path().join(format!("cm-{strat}"));
        let rows = run_histogram(&env.input, &prefix, &["--single-end-strategy", strat]);
        // Reduce to (n_observations, n_molecules) for the pairs category.
        let mut pairs: Vec<(String, String)> = rows
            .iter()
            .filter(|r| r["category"] == "pairs")
            .map(|r| (r["n_observations"].clone(), r["n_molecules"].clone()))
            .collect();
        pairs.sort();
        assert!(!pairs.is_empty(), "{strat}: expected pairs rows");
        match &reference {
            None => reference = Some(pairs),
            Some(reference) => assert_eq!(&pairs, reference, "{strat} pairs histogram differs"),
        }
    }
}

/// picard-exact (two-pass, reorders output) must still produce both files, with
/// pairs data — it used to be rejected.
#[test]
fn picard_exact_produces_both_files() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_pe_triple_plus_two_singletons(&env.input);
    let hist = run_histogram(&env.input, &prefix, &["--single-end-strategy", "picard-exact"]);
    assert_eq!(n_sigs(&hist, "pairs", "3"), 1);
    // The ladder file was written too, with the pairs curve.
    let mut ladder = prefix.as_os_str().to_owned();
    ladder.push(".duplication-sampled.tsv");
    let ladder_rows = parse_rows(Path::new(&ladder));
    assert!(ladder_rows.iter().any(|r| r["category"] == "pairs"));
}

/// A single-end-only library under a Picard strategy: no pairs means the
/// fragment keyspace has no pair-ends, so the single-end histogram is clean.
#[test]
fn se_only_under_picard_approx_is_clean() {
    let env = TestEnv::new();
    let prefix = env._tmp.path().join("cm");
    write_se_two_plus_one(&env.input);
    let rows = run_histogram(&env.input, &prefix, &["--single-end-strategy", "picard-approx"]);
    assert_eq!(n_sigs(&rows, "single_end", "1"), 1);
    assert_eq!(n_sigs(&rows, "single_end", "2"), 1);
}
