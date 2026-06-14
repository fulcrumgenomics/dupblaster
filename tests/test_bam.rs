//! BAM input/output integration tests.
//!
//! These tests verify that:
//! * SAM input → BAM round-trips: same dedup decisions.
//! * BAM input automatically produces BGZF-framed BAM output.
//! * BAM on stdin works end-to-end (the production aligner-on-stdin path).
//!
//! BAM is built and decoded in-process via noodles (`helpers::sam_to_bam`,
//! `helpers::read_recs_and_header`); no external `samtools` is required.

mod helpers;
use std::process::Command;

use helpers::*;

#[test]
fn bam_input_produces_bam_output_with_matching_dups() {
    let env = TestEnv::new();
    // Build a SAM with two duplicate pairs.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_in = env._tmp.path().join("in.bam");
    let bam_out = env._tmp.path().join("rust.out.bam");
    sam_to_bam(&env.input, &bam_in);

    let out = Command::new(rust_binary())
        .args(["-i"])
        .arg(&bam_in)
        .args(["-o"])
        .arg(&bam_out)
        .output()
        .unwrap();
    assert!(out.status.success(), "rust failed: {}", String::from_utf8_lossy(&out.stderr));

    // Decode the BAM output back to records for inspection.
    let (_header, records) = read_recs_and_header(&bam_out);
    assert_eq!(records.len(), 4);
    // r1 should keep its original flags (99, 147); r2 should have +0x400.
    let flags: Vec<u16> = records.iter().map(|r| u16::from(r.flags())).collect();
    assert_eq!(flags, vec![99, 147, 99 + 0x400, 147 + 0x400]);
}

#[test]
fn bam_output_uses_bgzf_framing() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_in = env._tmp.path().join("in.bam");
    let bam_out = env._tmp.path().join("rust.out.bam");
    sam_to_bam(&env.input, &bam_in);
    Command::new(rust_binary())
        .args(["-i"])
        .arg(&bam_in)
        .args(["-o"])
        .arg(&bam_out)
        .output()
        .unwrap();
    let bytes = std::fs::read(&bam_out).unwrap();
    // BGZF magic.
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04], "output is not BGZF-framed BAM");
}

#[test]
fn bam_input_routes_through_pipeline_via_stdin() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_in = env._tmp.path().join("in.bam");
    sam_to_bam(&env.input, &bam_in);

    // Pipe the BAM into the rust binary on stdin; output to a file.
    let bam_bytes = std::fs::read(&bam_in).unwrap();
    let bam_out = env._tmp.path().join("rust.stdin.bam");
    let mut child = Command::new(rust_binary())
        .args(["-o"])
        .arg(&bam_out)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    use std::io::Write as _;
    child.stdin.as_mut().unwrap().write_all(&bam_bytes).unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "rust failed on stdin BAM: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    // Output should be BAM (matching input format detected from the bytes).
    let bytes = std::fs::read(&bam_out).unwrap();
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04], "stdin BAM produced non-BAM output");
}
