//! --add-mate-tags behavioral tests.

mod helpers;
use std::process::Command;

use helpers::*;

/// Running --add-mate-tags on output that already has MC/MQ must NOT
/// append duplicate tags. The bug this covers: fgumi's append_int emits
/// the BAM `c` subtype for MAPQ values 0..=127, but find_uint8 only
/// matched `C` — so the dup-check missed the existing tag.
///
/// We count *physical* tag occurrences via the lazy BAM reader: `RecordBuf`
/// would de-dupe the fields and hide the regression.
#[test]
fn add_mate_tags_is_idempotent_across_two_runs() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "5S45M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "30M20S", "=", 100, -150)
        .write_to(&env.input);
    // First pass: add tags.
    let pass1 = env._tmp.path().join("pass1.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&pass1)
        .args(["--add-mate-tags"])
        .output()
        .unwrap();
    assert!(r.status.success(), "pass1: {}", String::from_utf8_lossy(&r.stderr));
    // Second pass: feed pass1 back in with --add-mate-tags.
    let pass2 = env._tmp.path().join("pass2.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&pass1)
        .args(["-o"])
        .arg(&pass2)
        .args(["--add-mate-tags"])
        .output()
        .unwrap();
    assert!(r.status.success(), "pass2: {}", String::from_utf8_lossy(&r.stderr));

    // Each tag key must appear exactly once per record (no duplicate append).
    for (i, n) in count_tag_occurrences(&pass2, *b"MC").iter().enumerate() {
        assert_eq!(*n, 1, "record {i} has {n} MC tags, expected 1");
    }
    for (i, n) in count_tag_occurrences(&pass2, *b"MQ").iter().enumerate() {
        assert_eq!(*n, 1, "record {i} has {n} MQ tags, expected 1");
    }
}

#[test]
fn add_mate_tags_writes_mate_cigar_and_mate_mapq() {
    // Build a pair with distinct CIGARs and MAPQs on each end so that
    // mate-tag values are unambiguously identifiable.
    //   R1: CIGAR=10M5S, MAPQ=55 → R2 should get MC:Z:10M5S, MQ:55.
    //   R2: CIGAR=50M,   MAPQ=33 → R1 should get MC:Z:50M,   MQ:33.
    let env = TestEnv::new();
    let mut sb = SamBuilder::new().sq("chr1", 1_000_000);
    sb.records.push(format!(
        "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}",
        qname = "r1",
        flag = 99,
        rname = "chr1",
        pos = 100,
        mapq = 55,
        cigar = "10M5S",
        rnext = "=",
        pnext = 200,
        tlen = 150,
        seq = "A".repeat(15),
        qual = "I".repeat(15),
    ));
    sb.records.push(format!(
        "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}",
        qname = "r1",
        flag = 147,
        rname = "chr1",
        pos = 200,
        mapq = 33,
        cigar = "50M",
        rnext = "=",
        pnext = 100,
        tlen = -150,
        seq = "A".repeat(50),
        qual = "I".repeat(50),
    ));
    sb.write_to(&env.input);

    let out_bam = env._tmp.path().join("out.bam");
    let r = Command::new(rust_binary())
        .args(["-i"])
        .arg(&env.input)
        .args(["-o"])
        .arg(&out_bam)
        .args(["--add-mate-tags"])
        .output()
        .unwrap();
    assert!(r.status.success(), "rust failed: {}", String::from_utf8_lossy(&r.stderr));

    let (_header, recs) = read_recs_and_header(&out_bam);
    assert_eq!(recs.len(), 2);

    // R1: MC should equal R2's CIGAR (50M), MQ should equal R2's MAPQ (33).
    assert_eq!(tag_string(&recs[0], *b"MC").as_deref(), Some("50M"), "R1 MC");
    assert_eq!(tag_int(&recs[0], *b"MQ"), Some(33), "R1 MQ");

    // R2: MC should equal R1's CIGAR (10M5S), MQ should equal R1's MAPQ (55).
    assert_eq!(tag_string(&recs[1], *b"MC").as_deref(), Some("10M5S"), "R2 MC");
    assert_eq!(tag_int(&recs[1], *b"MQ"), Some(55), "R2 MQ");
}
