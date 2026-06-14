//! Tests for the end-of-run resource footer (wall/CPU/RSS) and its
//! suppression under `--quiet`.

mod helpers;

use std::process::Command;

use helpers::*;

/// Build a tiny one-pair input and run dupblaster with `extra` flags,
/// returning the captured process output.
fn run(extra: &[&str]) -> std::process::Output {
    let env = TestEnv::new();
    let out = env._tmp.path().join("out.bam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out)
        .args(extra)
        .output()
        .expect("ran")
}

#[test]
fn footer_is_printed_without_quiet() {
    let r = run(&[]);
    assert!(r.status.success(), "stderr: {}", String::from_utf8_lossy(&r.stderr));
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("Processed") && stderr.contains("templates in"),
        "expected a resource footer, got: {stderr}"
    );
    // getrusage is available on the unix CI/dev targets, so the RSS field
    // should be present.
    assert!(stderr.contains("max RSS"), "expected RSS in the footer, got: {stderr}");
}

#[test]
fn footer_is_suppressed_with_quiet() {
    let r = run(&["--quiet"]);
    assert!(r.status.success(), "stderr: {}", String::from_utf8_lossy(&r.stderr));
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(!stderr.contains("Processed"), "footer should be suppressed under --quiet: {stderr}");
}
