//! Verify that dupblaster's @PG record chains via PP: to the previous
//! program in the header — and that re-running dupblaster on its own
//! output disambiguates the ID instead of producing a duplicate.

mod helpers;

use std::path::Path;
use std::process::Command;

use helpers::*;

/// Run the Rust binary writing BAM, then return its `@PG` header lines
/// (read in-process via noodles).
fn run_and_get_pg_lines(sam_in: &Path, bam_out: &Path) -> Vec<String> {
    let out = Command::new(rust_binary())
        .args(["-i"])
        .arg(sam_in)
        .args(["-o"])
        .arg(bam_out)
        .args(["--quiet"])
        .output()
        .expect("rust dupblaster ran");
    assert!(
        out.status.success(),
        "rust dupblaster failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    read_pg_lines(bam_out)
}

/// Pick out the value of a SAM tag (e.g. "ID" or "PP") from a tab-separated
/// `@PG` line, or return None.
fn tag<'a>(pg_line: &'a str, tag: &str) -> Option<&'a str> {
    pg_line.split('\t').find_map(|f| f.strip_prefix(&format!("{tag}:")))
}

#[test]
fn pg_record_chains_to_previous_leaf_program() {
    let env = TestEnv::new();
    let mut sb = SamBuilder::new().sq("chr1", 1_000_000);
    sb.header.push_str("@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n");
    sb.header.push_str("@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.18\n");
    sb = sb
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150);
    sb.write_to(&env.input);

    let bam_out = env._tmp.path().join("out.bam");
    let pg_lines = run_and_get_pg_lines(&env.input, &bam_out);

    // We expect three @PG lines: bwa, samtools (chained to bwa),
    // DUPBLASTER (chained to samtools — the existing chain leaf).
    assert_eq!(pg_lines.len(), 3, "got {pg_lines:?}");
    let dupblaster_line = pg_lines
        .iter()
        .find(|l| tag(l, "ID") == Some("DUPBLASTER"))
        .expect("DUPBLASTER PG present");
    assert_eq!(
        tag(dupblaster_line, "PP"),
        Some("samtools"),
        "expected PP:samtools, got {dupblaster_line}"
    );
}

#[test]
fn pg_record_has_no_pp_when_input_header_has_no_programs() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);

    let bam_out = env._tmp.path().join("out.bam");
    let pg_lines = run_and_get_pg_lines(&env.input, &bam_out);

    assert_eq!(pg_lines.len(), 1, "got {pg_lines:?}");
    assert_eq!(tag(&pg_lines[0], "ID"), Some("DUPBLASTER"));
    // No PP because we're the root.
    assert_eq!(tag(&pg_lines[0], "PP"), None);
}

#[test]
fn rerunning_dupblaster_disambiguates_pg_id() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);

    let pass1 = env._tmp.path().join("pass1.bam");
    let _ = run_and_get_pg_lines(&env.input, &pass1);
    // Second pass: feed pass1.bam back in. The existing header already has
    // @PG ID:DUPBLASTER from pass 1, so noodles must disambiguate.
    let pass2 = env._tmp.path().join("pass2.bam");
    let out = Command::new(rust_binary())
        .args(["-i"])
        .arg(&pass1)
        .args(["-o"])
        .arg(&pass2)
        .args(["--quiet"])
        .output()
        .unwrap();
    assert!(out.status.success(), "second pass failed: {}", String::from_utf8_lossy(&out.stderr));
    let pg_ids: Vec<String> =
        read_pg_lines(&pass2).iter().filter_map(|l| tag(l, "ID").map(String::from)).collect();
    // Two distinct DUPBLASTER-related IDs, and they should not duplicate.
    let unique: std::collections::HashSet<_> = pg_ids.iter().collect();
    assert_eq!(unique.len(), pg_ids.len(), "duplicate @PG ID in {pg_ids:?}");
    // Each ID starts with "DUPBLASTER".
    for id in &pg_ids {
        assert!(id.starts_with("DUPBLASTER"), "unexpected PG ID: {id}");
    }
}

#[test]
fn input_with_broken_pp_pointer_fails_with_clear_error() {
    // Input header contains a @PG whose PP points at an ID that doesn't
    // exist in the header. dupblaster's defensive validation should
    // produce a clean non-zero exit with a message naming the broken
    // chain. Without this check, noodles' programs.add() panics in the
    // middle of header processing with "no entry found for key".
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    let mut sb = SamBuilder::new().sq("chr1", 1_000_000);
    sb.header.push_str("@PG\tID:orphaned\tPN:orphaned\tPP:nonexistent\tVN:1.0\n");
    sb = sb
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150);
    sb.write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(!r.status.success(), "expected non-zero exit on broken PP chain");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("PP:nonexistent") || stderr.contains("malformed SAM header"),
        "stderr should explain the broken PP, got: {stderr}"
    );
}

#[test]
fn pg_record_chains_to_one_leaf_when_input_has_two_independent_chains() {
    // Two independent @PG chains (each a root with no PP). dupblaster
    // delegates leaf detection to noodles; this test pins the outcome:
    // the run completes, our @PG is present, and its PP points at one of
    // the two existing leaves (whichever noodles picks).
    let env = TestEnv::new();
    let mut sb = SamBuilder::new().sq("chr1", 1_000_000);
    sb.header.push_str("@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n");
    sb.header.push_str("@PG\tID:trimmer\tPN:trimmer\tVN:0.4.2\n");
    sb = sb
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150);
    sb.write_to(&env.input);

    let bam_out = env._tmp.path().join("out.bam");
    let pg_lines = run_and_get_pg_lines(&env.input, &bam_out);
    let dup = pg_lines
        .iter()
        .find(|l| tag(l, "ID").is_some_and(|id| id.starts_with("DUPBLASTER")))
        .expect("dupblaster @PG present");
    let pp = tag(dup, "PP").expect("dupblaster @PG has PP when input has existing leaves");
    assert!(
        pp == "bwa" || pp == "trimmer",
        "dupblaster PP should point at one of the existing leaves, got PP:{pp}"
    );
}
