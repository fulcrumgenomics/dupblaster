//! Tests for the startup attribution banner: dupblaster names Fulcrum Genomics
//! and links to its repository on stderr, keeps that off stdout (which carries
//! BAM), and suppresses it under `--quiet`.

mod helpers;

use helpers::*;
use tempfile::tempdir;

/// One duplicate pair's worth of input — enough to drive a full run.
fn tiny_input(dir: &std::path::Path) -> std::path::PathBuf {
    let sam = dir.join("in.sam");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("q1", 99, "chr1", 1000, "10M", "=", 1200, 210)
        .rec_simple("q1", 147, "chr1", 1200, "10M", "=", 1000, -210)
        .write_to(&sam);
    sam
}

fn stderr_of_run(extra_args: &[&str]) -> String {
    let dir = tempdir().unwrap();
    let sam = tiny_input(dir.path());
    let out = dir.path().join("out.bam");
    run_and_capture(&sam, &out, extra_args).stderr
}

#[test]
fn startup_banner_credits_fulcrum_genomics() {
    assert!(stderr_of_run(&[]).contains("by Fulcrum Genomics"));
}

#[test]
fn startup_banner_links_to_the_repository() {
    assert!(stderr_of_run(&[]).contains("https://github.com/fulcrumgenomics/dupblaster"));
}

#[test]
fn startup_banner_reports_the_version() {
    assert!(stderr_of_run(&[]).contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn quiet_suppresses_the_startup_banner() {
    let stderr = stderr_of_run(&["--quiet"]);
    assert!(!stderr.contains("Fulcrum Genomics"), "stderr: {stderr:?}");
    assert!(!stderr.contains("github.com"), "stderr: {stderr:?}");
}

#[test]
fn quiet_still_reports_the_duplicate_summary() {
    // `--quiet` trims the progress lines, not the result the run exists to
    // produce — the reworded flag help promises exactly this.
    assert!(stderr_of_run(&["--quiet"]).contains("as duplicates"));
}

#[test]
fn piped_stdout_is_a_valid_bam_stream() {
    // With `-o -` the banner shares a process with the BAM stream. Parse stdout
    // as BAM rather than scanning it for banner text: that way *any* stray write
    // fails the test, not just this one string.
    let dir = tempdir().unwrap();
    let sam = tiny_input(dir.path());
    // The metrics prefix goes in the temp dir, not beside `-`, which would derive
    // a path in the working directory.
    let r = dupblaster()
        .arg("-i")
        .arg(&sam)
        .arg("-o")
        .arg("-")
        .arg("--metrics-prefix")
        .arg(metrics_prefix_in(dir.path()))
        .output()
        .unwrap();
    assert!(r.status.success(), "stderr: {}", String::from_utf8_lossy(&r.stderr));

    let piped = dir.path().join("piped.bam");
    std::fs::write(&piped, &r.stdout).expect("capture piped BAM");
    let (header, records) = read_recs_and_header(&piped);
    assert!(!header.reference_sequences().is_empty(), "piped BAM lost its @SQ header");
    let names: Vec<String> =
        records.iter().filter_map(|rec| rec.name().map(|n| n.to_string())).collect();
    assert_eq!(names, ["q1", "q1"], "piped BAM should carry both records intact");

    // ...and the attribution still reached the user, on stderr.
    assert!(String::from_utf8_lossy(&r.stderr).contains("Fulcrum Genomics"));
}

#[test]
fn help_banner_names_fulcrum_genomics_and_the_repository() {
    let r = dupblaster().arg("--help").output().unwrap();
    let stdout = String::from_utf8_lossy(&r.stdout);
    assert!(stdout.contains("dupblaster by Fulcrum Genomics"), "stdout: {stdout:?}");
    assert!(stdout.contains("https://github.com/fulcrumgenomics/dupblaster"), "stdout: {stdout:?}");
}

#[test]
fn help_no_longer_points_at_the_readme() {
    // The banner URL is the single pointer to full documentation; per-option
    // help must stand on its own.
    for flag in ["-h", "--help"] {
        let r = dupblaster().arg(flag).output().unwrap();
        let stdout = String::from_utf8_lossy(&r.stdout);
        assert!(!stdout.contains("README"), "{flag} still mentions the README: {stdout:?}");
    }
}
