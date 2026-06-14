//! Tests for dupblaster's handling of unmated paired records — a primary
//! record whose mate is missing from the input stream entirely.

mod helpers;

use std::process::Command;

use helpers::*;

#[test]
fn single_paired_primary_with_no_mate_aborts_without_ignore_unmated() {
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    // One record only: FLAG=99 (paired + proper_pair + mate_rev + first).
    // No mate, no second record. dupblaster should abort with a non-zero
    // exit code and a "broken block" / read-id sortedness error.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(!r.status.success(), "expected dupblaster to fail on unmated input");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("first and/or second of pair") || stderr.contains("sorted by read ids"),
        "stderr should explain unmated/sortedness, got: {stderr}"
    );
}

#[test]
fn single_paired_primary_with_no_mate_succeeds_under_ignore_unmated() {
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    // Mix one normal pair (so total_templates > 0) with one unmated record.
    // Under --ignore-unmated, the run should succeed; the unmated template
    // is counted in the stats but not crashed on.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 500, "50M", "=", 600, 150)
        .write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .args(["--ignore-unmated"])
        .output()
        .unwrap();
    assert!(
        r.status.success(),
        "rust failed under --ignore-unmated: {}",
        String::from_utf8_lossy(&r.stderr)
    );
}
