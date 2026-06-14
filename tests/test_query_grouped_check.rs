//! Tests for dupblaster's detection of non-query-grouped input.
//!
//! When a QNAME block contains only secondary/supplementary alignments
//! (no primary), the primary must live elsewhere in the stream — the
//! signature of non-query-grouped (typically coordinate-sorted) input.
//! dupblaster should fail loudly with a message that names the cause.

mod helpers;
use std::process::Command;

use helpers::*;

/// SAM FLAG constants used in the test fixtures below.
const PAIRED: u16 = 0x1;
const FIRST_SEGMENT: u16 = 0x40;
const LAST_SEGMENT: u16 = 0x80;
const SUPPLEMENTARY: u16 = 0x800;
const SECONDARY: u16 = 0x100;

#[test]
fn fails_when_block_has_only_supplementary_records() {
    // A QNAME block containing nothing but a supplementary alignment.
    // In query-grouped input the primary would be in the same block;
    // its absence means the input is not query-grouped.
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", PAIRED | FIRST_SEGMENT | SUPPLEMENTARY, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(!r.status.success(), "should have failed on non-query-grouped input");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("not query-grouped"), "stderr was: {stderr}");
}

#[test]
fn fails_when_block_has_only_secondary_records() {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", PAIRED | LAST_SEGMENT | SECONDARY, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(!r.status.success(), "should have failed on non-query-grouped input");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("not query-grouped"), "stderr was: {stderr}");
}

#[test]
fn fails_even_with_ignore_unmated() {
    // --ignore-unmated must NOT silence the non-query-grouped detector;
    // this is a configuration-error condition, not a tolerable edge case.
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", PAIRED | FIRST_SEGMENT | SUPPLEMENTARY, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--ignore-unmated"])
        .output()
        .unwrap();
    assert!(!r.status.success(), "--ignore-unmated must not suppress the non-query-grouped check");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("not query-grouped"), "stderr was: {stderr}");
}
