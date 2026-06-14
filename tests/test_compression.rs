//! `--compression-level` integration tests.
//!
//! Default (level 0) produces uncompressed BGZF; non-zero levels produce
//! truly compressed BGZF blocks. Either way the output must remain a
//! valid BAM that a standard reader (noodles) can decode.

mod helpers;
use std::process::Command;

use helpers::*;

fn build_input(env: &TestEnv) {
    // Repeat the same record many times so a compressor has something to
    // crunch — at level 0 vs level 6 we want a clear size delta.
    let mut sb = SamBuilder::new().sq("chr1", 1_000_000);
    for i in 0..500 {
        let qname = format!("r{i}");
        sb = sb.rec_simple(&qname, 99, "chr1", 100 + i, "50M", "=", 200 + i, 150).rec_simple(
            &qname,
            147,
            "chr1",
            200 + i,
            "50M",
            "=",
            100 + i,
            -150,
        );
    }
    sb.write_to(&env.input);
}

#[test]
fn default_output_is_uncompressed_bgzf() {
    let env = TestEnv::new();
    build_input(&env);
    let out = env._tmp.path().join("out.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));

    let bytes = std::fs::read(&out).unwrap();
    // Still BGZF-framed even though "uncompressed" — that's the standard
    // samtools-view-u format.
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04], "expected BGZF framing");
}

#[test]
fn level_6_produces_smaller_output_than_default() {
    let env = TestEnv::new();
    build_input(&env);

    let out0 = env._tmp.path().join("out0.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out0)
        .output()
        .unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let size_0 = std::fs::metadata(&out0).unwrap().len();

    let out6 = env._tmp.path().join("out6.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out6)
        .args(["--compression-level", "6"])
        .output()
        .unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let size_6 = std::fs::metadata(&out6).unwrap().len();

    assert!(
        size_6 < size_0,
        "level 6 output ({size_6} B) should be smaller than level 0 ({size_0} B)"
    );
}

#[test]
fn compressed_output_round_trips() {
    // Level-6 output must remain a valid BAM whose records decode back intact.
    let env = TestEnv::new();
    build_input(&env);
    let out = env._tmp.path().join("out.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--compression-level", "6"])
        .output()
        .unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));

    let records = read_records(&out);
    assert_eq!(records.len(), 1000, "should round-trip 1000 records");
}

#[test]
fn rejects_out_of_range_compression_level() {
    let env = TestEnv::new();
    build_input(&env);
    let out = env._tmp.path().join("out.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(["--compression-level", "13"])
        .output()
        .unwrap();
    assert!(!r.status.success(), "level 13 should be rejected by the value parser");
}
