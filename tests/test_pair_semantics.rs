//! Generative tests for dupblaster's pair-level duplicate semantics:
//! ordering invariance, three-deep duplicate sets, --remove-dups, and
//! cross-chromosome signature uniqueness.

mod helpers;

use helpers::*;

#[test]
fn pair_order_swap_does_not_change_duplicate_decision() {
    let env = TestEnv::new();
    // Pair A: R1 first, then R2. Pair B: R2 first, then R1. Same
    // coordinates and strands — the signature must canonicalize so
    // both pairs are detected as duplicates regardless of arrival order.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        // rB written R2 first, then R1.
        .rec_simple("rB", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("rB", 99, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    // Output order matches input order: rB's R2 came first, then R1.
    assert_eq!(recs[2], ("rB".into(), 147 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 99 | FLAG_DUPLICATE));
}

#[test]
fn three_identical_pairs_only_first_kept_as_non_duplicate() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r3", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r3", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    // r1 is the representative (no dup flag); r2 and r3 are duplicates.
    assert_eq!(recs[0], ("r1".into(), 99));
    assert_eq!(recs[1], ("r1".into(), 147));
    assert_eq!(recs[2], ("r2".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("r2".into(), 147 | FLAG_DUPLICATE));
    assert_eq!(recs[4], ("r3".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[5], ("r3".into(), 147 | FLAG_DUPLICATE));
}

#[test]
fn remove_dups_drops_duplicate_records_from_output() {
    let env = TestEnv::new();
    // Three identical pairs. With -r, only r1's two records should remain.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("r3", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("r3", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["-r"]);

    // Only the kept pair survives. No record carries FLAG_DUPLICATE
    // because the dup records have been removed.
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0], ("r1".into(), 99));
    assert_eq!(recs[1], ("r1".into(), 147));
}

#[test]
fn cross_chromosome_pairs_at_distinct_positions_are_not_dups() {
    let env = TestEnv::new();
    // Two distinct cross-chrom pairs. The two-key signature includes the
    // chromosome of each end, so different (chr,pos) combinations must
    // not collide.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .sq("chr2", 1_000_000)
        // Pair A: R1 on chr1:100 fwd, R2 on chr2:200 rev.
        .rec_simple("rA", 97, "chr1", 100, "50M", "chr2", 200, 0)
        .rec_simple("rA", 145, "chr2", 200, "50M", "chr1", 100, 0)
        // Pair B: R1 on chr1:300 fwd, R2 on chr2:400 rev — different positions.
        .rec_simple("rB", 97, "chr1", 300, "50M", "chr2", 400, 0)
        .rec_simple("rB", 145, "chr2", 400, "50M", "chr1", 300, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) should not be marked duplicate"
        );
    }
}

#[test]
fn cross_chromosome_pairs_at_identical_positions_are_dups() {
    let env = TestEnv::new();
    // Same as the "distinct positions" test but pair B mirrors A's
    // coordinates exactly — must be detected as a dup.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .sq("chr2", 1_000_000)
        .rec_simple("rA", 97, "chr1", 100, "50M", "chr2", 200, 0)
        .rec_simple("rA", 145, "chr2", 200, "50M", "chr1", 100, 0)
        .rec_simple("rB", 97, "chr1", 100, "50M", "chr2", 200, 0)
        .rec_simple("rB", 145, "chr2", 200, "50M", "chr1", 100, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), 97));
    assert_eq!(recs[1], ("rA".into(), 145));
    assert_eq!(recs[2], ("rB".into(), 97 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 145 | FLAG_DUPLICATE));
}
