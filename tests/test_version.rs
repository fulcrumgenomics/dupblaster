//! Tests for the `--version` / `-V` flag: it prints the version to stdout
//! (so it can be piped/captured) and exits 0 without requiring any input.

mod helpers;

use std::process::Command;

use helpers::*;

#[test]
fn version_long_flag_prints_to_stdout_and_exits_zero() {
    let r = Command::new(rust_binary()).arg("--version").output().unwrap();
    assert!(r.status.success(), "--version should exit 0");
    let stdout = String::from_utf8_lossy(&r.stdout);
    assert!(stdout.starts_with("dupblaster"), "stdout: {stdout:?}");
    assert!(stdout.contains(env!("CARGO_PKG_VERSION")), "stdout: {stdout:?}");
    assert!(
        r.stderr.is_empty(),
        "stderr should be empty: {:?}",
        String::from_utf8_lossy(&r.stderr)
    );
}

#[test]
fn version_short_flag_matches_long_flag() {
    let long = Command::new(rust_binary()).arg("--version").output().unwrap();
    let short = Command::new(rust_binary()).arg("-V").output().unwrap();
    assert!(short.status.success(), "-V should exit 0");
    assert_eq!(short.stdout, long.stdout, "-V and --version should print identically");
}
