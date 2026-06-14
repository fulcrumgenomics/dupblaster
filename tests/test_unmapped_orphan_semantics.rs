//! Generative tests for dupblaster's handling of unmapped reads and
//! orphan pairs (one mate unmapped, the other mapped).

mod helpers;

use helpers::*;

#[test]
fn both_unmapped_pairs_are_never_marked_duplicate() {
    let env = TestEnv::new();
    // Two pairs where every read is unmapped (PAIRED|UNMAPPED|MATE_UNMAPPED).
    // R1: 0x1|0x4|0x8|0x40 = 77. R2: 0x1|0x4|0x8|0x80 = 141.
    // Unmapped pairs never participate in dup-marking; both must stay clean.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .record("r1", 77, "*", 0, 0, "*", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .record("r1", 141, "*", 0, 0, "*", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .record("r2", 77, "*", 0, 0, "*", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .record("r2", 141, "*", 0, 0, "*", "*", 0, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs.len(), 4);
    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) is fully unmapped and must not be marked dup"
        );
    }
}

#[test]
fn two_orphan_pairs_with_same_mapped_position_are_dups() {
    let env = TestEnv::new();
    // Pair A: R1 mapped at chr1:100 forward, R2 unmapped (mate of R1).
    // Pair B: same shape, same coordinates — must be detected as a dup.
    //
    // R1 flag = PAIRED|MATE_UNMAPPED|FIRST_SEGMENT = 0x1|0x8|0x40 = 73.
    // R2 flag = PAIRED|UNMAPPED|LAST_SEGMENT       = 0x1|0x4|0x80 = 133.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rA", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        .rec_simple("rB", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rB", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    // rA is the kept representative; rB is the duplicate. The dup flag
    // propagates to every record in rB's block (including the unmapped
    // mate), matching the "primary dup → mark all" rule.
    assert_eq!(recs[0], ("rA".into(), 73));
    assert_eq!(recs[1], ("rA".into(), 133));
    assert_eq!(recs[2], ("rB".into(), 73 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 133 | FLAG_DUPLICATE));
}

#[test]
fn orphan_pairs_with_different_mapped_positions_are_not_dups() {
    let env = TestEnv::new();
    // Two orphan pairs where the mapped read's 5' position differs.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rA", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        .rec_simple("rB", 73, "chr1", 500, "50M", "=", 500, 0)
        .record("rB", 133, "chr1", 500, 0, "*", "=", 500, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) has a distinct mapped 5' position and must not be marked dup"
        );
    }
}

#[test]
fn two_reverse_strand_single_end_orphans_at_same_position_are_dups() {
    // Single-end reads (FLAG_PAIRED = 0) on the reverse strand. The
    // signature code path for orphans is distinct from paired-end, and
    // the strand bit shouldn't break dup detection between two
    // identically-mapped reverse-strand reads. (Forward-vs-reverse
    // orphan equivalence is a separate question tracked in task #56.)
    let env = TestEnv::new();
    // FLAG 16 = REVERSE only (no PAIRED bit). 50M @ chr1:100.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), FLAG_REVERSE));
    assert_eq!(recs[1], ("rB".into(), FLAG_REVERSE | FLAG_DUPLICATE));
}

#[test]
fn two_reverse_strand_single_end_orphans_at_different_positions_are_not_dups() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 500, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) has a distinct 5' position and must not be marked dup"
        );
    }
}
