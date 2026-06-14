//! Integration tests for library-aware duplicate marking.
//!
//! Duplicates are called only *within* a library (matching Picard
//! MarkDuplicates). Library membership comes from each read's `RG:Z` tag mapped
//! through the header's `@RG ... LB:` field. Library awareness is on by default
//! whenever the header declares more than one distinct `LB`; `--library-unaware`
//! forces the old single-table (library-agnostic) behavior.

mod helpers;

use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use helpers::*;

/// Map each qname to whether any of its records carries the duplicate flag.
fn dup_by_qname(flags: &[(String, u16)]) -> HashMap<String, bool> {
    let mut m: HashMap<String, bool> = HashMap::new();
    for (q, f) in flags {
        let e = m.entry(q.clone()).or_insert(false);
        *e |= f & FLAG_DUPLICATE != 0;
    }
    m
}

/// Two pairs at the *same* coordinates but in *different* libraries are not
/// duplicates of each other — the whole point of library-aware marking.
#[test]
fn cross_library_pairs_at_same_locus_are_not_duplicates() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        // Pair p1 in lib1.
        .rec_simple_rg("p1", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        // Pair p2 in lib2 at identical coordinates.
        .rec_simple_rg("p2", 99, "chr1", 100, "50M", "=", 200, 150, "B")
        .rec_simple_rg("p2", 147, "chr1", 200, "50M", "=", 100, -150, "B")
        .write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &[]);
    let dup = dup_by_qname(&flags);
    assert!(!dup["p1"], "p1 should not be a duplicate");
    assert!(!dup["p2"], "p2 is a different library — not a duplicate of p1");
}

/// Two read groups sharing one `LB` are one library, so identical-coordinate
/// pairs across them still collapse to a duplicate.
#[test]
fn same_library_pairs_across_read_groups_are_duplicates() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Two RGs, same LB → same library.
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib1"))
        .rec_simple_rg("p1", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .rec_simple_rg("p2", 99, "chr1", 100, "50M", "=", 200, 150, "B")
        .rec_simple_rg("p2", 147, "chr1", 200, "50M", "=", 100, -150, "B")
        .write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &[]);
    let dup = dup_by_qname(&flags);
    assert!(!dup["p1"], "first pair seen is the original");
    assert!(dup["p2"], "same library, same coords → duplicate");
}

/// `--library-unaware` restores the single-table behavior: cross-library pairs
/// at the same locus collapse to a duplicate again.
#[test]
fn library_unaware_flag_marks_cross_library_dups() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        .rec_simple_rg("p1", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .rec_simple_rg("p2", 99, "chr1", 100, "50M", "=", 200, 150, "B")
        .rec_simple_rg("p2", 147, "chr1", 200, "50M", "=", 100, -150, "B")
        .write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &["--library-unaware"]);
    let dup = dup_by_qname(&flags);
    assert!(dup["p2"], "with --library-unaware, different libraries collapse");
}

/// Reads with no `RG` tag share the single "Unknown Library" bucket and so
/// still dedup against each other, even when other libraries are present.
#[test]
fn reads_without_read_group_share_unknown_library() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        // One pair carries a real RG (keeps the header multi-library so
        // library-aware mode is active), the other two carry none.
        .rec_simple_rg("anchored", 99, "chr1", 700, "50M", "=", 800, 150, "A")
        .rec_simple_rg("anchored", 147, "chr1", 800, "50M", "=", 700, -150, "A")
        .rec_simple("u1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("u1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("u2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("u2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &[]);
    let dup = dup_by_qname(&flags);
    assert!(!dup["u1"], "first RG-less pair is the original");
    assert!(dup["u2"], "second RG-less pair dedups within Unknown Library");
}

/// The "Unknown Library" bucket is distinct from every named library: an
/// RG-less pair at the same locus as a named-library pair is NOT a duplicate of
/// it. The named library also still dedups correctly across the interleaved
/// no-RG block (exercises library resolution / the per-block RG lookup).
#[test]
fn unknown_library_is_distinct_from_named_libraries() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        // lib1 pair, then an RG-less pair, then another lib1 pair — all at the
        // same coordinates.
        .rec_simple_rg("p1", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .rec_simple("pN", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("pN", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple_rg("p2", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p2", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &[]);
    let dup = dup_by_qname(&flags);
    assert!(!dup["p1"], "first lib1 pair is the original");
    assert!(!dup["pN"], "RG-less pair is the Unknown Library — not a dup of lib1");
    assert!(dup["p2"], "second lib1 pair dedups within lib1 across the no-RG block");
}

/// Library-aware marking composes with the picard-exact two-pass strategy: the
/// per-library fragment tables (drained from per-library pair tables) keep
/// orphan dedup within a library. An orphan sharing a pair end's 5' position is
/// a duplicate only when it's in the *same* library as that pair.
#[test]
fn picard_exact_orphans_dedup_within_library_only() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    let unmapped_mate = |b: SamBuilder, q: &str| {
        // R2 unmapped (flag 133 = paired | unmapped | last); its mate (R1) is
        // the mapped orphan end. No RG needed here — library is resolved from
        // the block's first record (the mapped R1, which carries RG).
        b.record(q, 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
    };
    let mut b = SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        // lib1 pair with its R1 5' end at chr1:100.
        .rec_simple_rg("pairA", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("pairA", 147, "chr1", 200, "50M", "=", 100, -150, "A");
    // lib1 orphan at chr1:100 — same library as pairA → dup of pairA's end.
    b = b.rec_simple_rg("orphanA", 73, "chr1", 100, "50M", "=", 100, 0, "A");
    b = unmapped_mate(b, "orphanA");
    // lib2 orphan at chr1:100 — different library, no lib2 pair/fragment there.
    b = b.rec_simple_rg("orphanB", 73, "chr1", 100, "50M", "=", 100, 0, "B");
    b = unmapped_mate(b, "orphanB");
    b.write_to(&env.input);

    let flags = run_and_extract_flags(&env.input, &out, &["--single-end-strategy", "picard-exact"]);
    let dup = dup_by_qname(&flags);
    assert!(!dup["pairA"], "the pair is the original");
    assert!(dup["orphanA"], "lib1 orphan at the pair's 5' end is a duplicate of the pair");
    assert!(!dup["orphanB"], "lib2 orphan at the same locus is a different library — not a dup");
}

/// Parse a multi-row stats TSV into one column→value map per data row.
fn parse_stats_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).expect("read stats");
    let mut lines = text.lines();
    let cols: Vec<String> = lines.next().expect("header").split('\t').map(String::from).collect();
    lines
        .map(|line| {
            let vals: Vec<&str> = line.split('\t').collect();
            assert_eq!(vals.len(), cols.len(), "row/column count mismatch");
            cols.iter().cloned().zip(vals.iter().map(|s| s.to_string())).collect()
        })
        .collect()
}

/// `--stats` emits one row per library, each with its own counts and dup rate,
/// in sorted-LB order. The empty "Unknown Library" bucket is omitted.
#[test]
fn stats_reports_one_row_per_library() {
    let env = TestEnv::new();
    let stats = env._tmp.path().join("stats.tsv");
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rg("A", "s1", Some("lib1"))
        .rg("B", "s1", Some("lib2"))
        // lib1: 3 pairs, one a duplicate of another (p1 == p1dup), plus unique p3.
        .rec_simple_rg("p1", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .rec_simple_rg("p1dup", 99, "chr1", 100, "50M", "=", 200, 150, "A")
        .rec_simple_rg("p1dup", 147, "chr1", 200, "50M", "=", 100, -150, "A")
        .rec_simple_rg("p3", 99, "chr1", 500, "50M", "=", 600, 150, "A")
        .rec_simple_rg("p3", 147, "chr1", 600, "50M", "=", 500, -150, "A")
        // lib2: a single unique pair at the same coords as lib1's p1 — not a dup.
        .rec_simple_rg("p4", 99, "chr1", 100, "50M", "=", 200, 150, "B")
        .rec_simple_rg("p4", 147, "chr1", 200, "50M", "=", 100, -150, "B")
        .write_to(&env.input);

    let run = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--stats"])
        .arg(&stats)
        .output()
        .expect("dupblaster ran");
    assert!(run.status.success(), "dupblaster failed: {}", String::from_utf8_lossy(&run.stderr));

    let rows = parse_stats_rows(&stats);
    let libs: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(libs, ["lib1", "lib2"], "one row per library, sorted by LB");

    let lib1 = &rows[0];
    assert_eq!(lib1["total_templates"], "3");
    assert_eq!(lib1["mapped_pairs"], "3");
    assert_eq!(lib1["duplicate_pairs"], "1");

    let lib2 = &rows[1];
    assert_eq!(lib2["total_templates"], "1");
    assert_eq!(lib2["mapped_pairs"], "1");
    assert_eq!(lib2["duplicate_pairs"], "0");
}
