//! Verify that when a primary read pair is detected as a duplicate, the
//! FLAG 0x400 (PCR/optical duplicate) bit is set on *every* record in the
//! block — primary first, primary second, secondary alignments, and
//! supplementary alignments alike.
//!
//! These tests check absolute flag values in the output BAM (decoded via
//! samtools) rather than comparing fingerprints with C++ samblaster, so
//! they document our intent independently and run without a C++ binary.

mod helpers;

use helpers::*;

#[test]
fn duplicate_block_propagates_flag_to_supplementary() {
    let env = TestEnv::new();
    // r2 is a duplicate of r1 and has a supplementary alignment.
    // FLAG 2113 = 0x1|0x40|0x800 (paired + first + supplementary).
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .record(
            "r2",
            2113,
            "chr1",
            500,
            60,
            "25S25M",
            "=",
            200,
            0,
            &"A".repeat(50),
            &"I".repeat(50),
        )
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    // r1 records: no dup flag.
    assert_eq!(recs[0], ("r1".into(), 99));
    assert_eq!(recs[1], ("r1".into(), 147));
    // r2 records (all three): dup flag set.
    assert_eq!(recs[2], ("r2".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("r2".into(), 147 | FLAG_DUPLICATE));
    assert_eq!(recs[4], ("r2".into(), 2113 | FLAG_DUPLICATE));
}

#[test]
fn duplicate_block_propagates_flag_to_secondary() {
    let env = TestEnv::new();
    // r2 is a duplicate of r1 and has a secondary alignment.
    // FLAG 355 = 0x1|0x2|0x40|0x100|0x20 (paired + propPair + first + secondary + mateRev).
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 355, "chr1", 500, "50M", "=", 100, 0)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("r1".into(), 99));
    assert_eq!(recs[1], ("r1".into(), 147));
    assert_eq!(recs[2], ("r2".into(), 99 | FLAG_DUPLICATE));
    // Secondary picks up the dup flag too.
    assert_eq!(recs[3], ("r2".into(), 355 | FLAG_DUPLICATE));
    assert_eq!(recs[4], ("r2".into(), 147 | FLAG_DUPLICATE));
}

#[test]
fn duplicate_block_propagates_flag_to_both_secondary_and_supplementary() {
    let env = TestEnv::new();
    // r2 is a duplicate of r1 with BOTH a secondary AND a supplementary.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 355, "chr1", 500, "50M", "=", 100, 0) // secondary
        .record(
            "r2",
            2113,
            "chr1",
            700,
            60,
            "25S25M",
            "=",
            200,
            0,
            &"A".repeat(50),
            &"I".repeat(50),
        ) // supplementary
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("r1".into(), 99));
    assert_eq!(recs[1], ("r1".into(), 147));
    // All four r2 records flagged dup.
    for (i, want_flag) in [99u16, 355, 2113, 147].into_iter().enumerate() {
        assert_eq!(
            recs[2 + i],
            ("r2".into(), want_flag | FLAG_DUPLICATE),
            "r2 record {} (input flag {}) did not get FLAG_DUPLICATE",
            i,
            want_flag,
        );
    }
}

#[test]
fn non_duplicate_block_with_secondary_and_supplementary_keeps_flags_unchanged() {
    let env = TestEnv::new();
    // Two distinct read pairs, each with a secondary and a supplementary —
    // none are duplicates of each other, so no record should gain 0x400.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r1", 355, "chr1", 800, "50M", "=", 100, 0)
        .record(
            "r1",
            2113,
            "chr1",
            900,
            60,
            "25S25M",
            "=",
            200,
            0,
            &"A".repeat(50),
            &"I".repeat(50),
        )
        .rec_simple("r2", 99, "chr1", 1000, "50M", "=", 1100, 150)
        .rec_simple("r2", 147, "chr1", 1100, "50M", "=", 1000, -150)
        .rec_simple("r2", 355, "chr1", 1500, "50M", "=", 1000, 0)
        .record(
            "r2",
            2113,
            "chr1",
            1600,
            60,
            "25S25M",
            "=",
            1100,
            0,
            &"A".repeat(50),
            &"I".repeat(50),
        )
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    // No record should have 0x400 set.
    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{} (flag {}) unexpectedly got FLAG_DUPLICATE",
            qname,
            flag,
        );
    }
    // Sanity: we have all 8 records.
    assert_eq!(recs.len(), 8);
}
