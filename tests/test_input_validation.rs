//! Tests for dupblaster's handling of reads whose 5' clip extends past a
//! contig edge by more than `--max-read-length`. Such reads can't be
//! represented exactly on dupblaster's synthetic single-axis genome, so their
//! coordinate is clamped into the contig. Clamping is a correctness
//! compromise (potential false duplicates among clamped reads), so dupblaster
//! still writes its output, warns, and exits non-zero — and bumping
//! `--max-read-length` past the clip makes the run clean again.

mod helpers;

use std::process::Command;

use helpers::*;

/// A read mapped near a contig start with thousands of bases of leading
/// hard-clip clips past the default (1000 bp) padding. dupblaster clamps it,
/// warns, and exits non-zero — but still writes the output.
#[test]
fn clip_beyond_padding_clamps_warns_and_exits_nonzero() {
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    // A proper pair whose first end has a 5000-base leading hard clip at
    // chr1:1. forward 5' pos = 1 - 5000 ≪ 0, clipping ~4000 bp past the
    // 1000 bp padding.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 1, "5000H50M", "=", 200, 249)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 1, -249)
        .write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(!r.status.success(), "clamping should make dupblaster exit non-zero");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("clamped"), "stderr should mention clamping, got: {stderr}");
    assert!(
        stderr.contains("max-read-length"),
        "stderr should point at --max-read-length, got: {stderr}"
    );
    // Output is still written despite the non-zero exit: the pair round-trips.
    let records = read_records(&out_bam);
    assert_eq!(records.len(), 2, "output should still be written when clamping");
}

/// The same input runs clean (exit 0, no clamp warning) once
/// `--max-read-length` is bumped past the clip.
#[test]
fn bumping_max_read_length_avoids_clamping() {
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 1, "5000H50M", "=", 200, 249)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 1, -249)
        .write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .args(["--max-read-length", "6000"])
        .output()
        .unwrap();
    assert!(
        r.status.success(),
        "no clamping expected with a larger --max-read-length: {}",
        String::from_utf8_lossy(&r.stderr)
    );
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(!stderr.contains("clamped"), "no clamp warning expected, got: {stderr}");
    let records = read_records(&out_bam);
    assert_eq!(records.len(), 2);
}

/// A normal pair well inside the contig neither clamps nor exits non-zero.
#[test]
fn normal_reads_do_not_clamp() {
    let env = TestEnv::new();
    let out_bam = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 5000, "50M", "=", 5200, 250)
        .rec_simple("r1", 147, "chr1", 5200, "50M", "=", 5000, -250)
        .write_to(&env.input);

    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(r.status.success(), "stderr: {}", String::from_utf8_lossy(&r.stderr));
    assert!(!String::from_utf8_lossy(&r.stderr).contains("clamped"));
}
