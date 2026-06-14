//! Tests for the `--single-end-strategy` flag.
//!
//! The two strategies wired up here:
//!
//!   * `strand-aware` (default): forward-strand and reverse-strand
//!     orphans at the same 5' position are NOT duplicates of each
//!     other.
//!   * `samblaster-legacy`: forward and reverse orphans whose
//!     alignments share a *leftmost-aligned* reference coordinate
//!     collide regardless of strand (samblaster v0.1.23+ behavior).
//!
//! The third value `picard-approx` is implemented in a later commit
//! and tested separately.

mod helpers;

use helpers::*;

/// Forward-strand single-end at chr1:100, reverse-strand single-end
/// at chr1:151 (a 50bp 50M alignment starting at 151 has 5'-end at
/// 151 + 50 - 1 = 200; its leftmost-aligned coord is 151). Under
/// `samblaster-legacy` we *also* override the rev pos to leftmost, so
/// the forward read at leftmost=100 and the reverse read at
/// leftmost=151 do NOT collide. So this test exercises something
/// different: two reverse-strand orphans at the same leftmost (=
/// same rapos+sclip) — should collide under both strategies.
#[test]
fn samblaster_legacy_two_rev_orphans_same_leftmost_collide() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(
        &env.input,
        &bam_out,
        &["--single-end-strategy", "samblaster-legacy"],
    );
    assert_eq!(recs[0], ("rA".into(), FLAG_REVERSE));
    assert_eq!(recs[1], ("rB".into(), FLAG_REVERSE | FLAG_DUPLICATE));
}

/// Under `samblaster-legacy`, a forward orphan whose alignment starts
/// at chr1:100 (leftmost = 100) and a reverse orphan whose alignment
/// starts at chr1:100 (leftmost = 100 by `rapos - sclip`, with the
/// post-2020 leftmost override) collide as duplicates — *even though*
/// their 5'-aligned positions are different (forward 5' = 100,
/// reverse 5' = 149 for a 50M alignment).
///
/// This is the test that pins the legacy false-positive behavior.
#[test]
fn samblaster_legacy_fwd_and_rev_orphan_same_leftmost_collide() {
    let env = TestEnv::new();
    // FLAG 0 (forward, no PAIRED) and FLAG 16 (reverse) both at chr1:100
    // with the same 50M CIGAR. Leftmost = 100 for both.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(
        &env.input,
        &bam_out,
        &["--single-end-strategy", "samblaster-legacy"],
    );
    assert_eq!(recs[0], ("rA".into(), 0));
    // Under legacy, the rev orphan at the same leftmost as rA IS marked dup.
    assert_eq!(recs[1], ("rB".into(), FLAG_REVERSE | FLAG_DUPLICATE));
}

/// Same input as the legacy test above, but under the default
/// `strand-aware` strategy: fwd and rev orphans at the same leftmost
/// do NOT collide because their strand-aware 5' positions are
/// different.
#[test]
fn strand_aware_fwd_and_rev_orphan_same_leftmost_do_not_collide() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    // strand-aware is the default; pass it explicitly for documentation.
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "strand-aware"]);
    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) must NOT be marked dup under strand-aware"
        );
    }
}

/// Strand-aware default: two reverse-strand single-end reads at the
/// same chr1:100 still collide (same strand, same 5'-aligned position).
/// This is the case that should keep working regardless of strategy.
#[test]
fn strand_aware_two_rev_orphans_same_position_collide() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", FLAG_REVERSE, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "strand-aware"]);
    assert_eq!(recs[0], ("rA".into(), FLAG_REVERSE));
    assert_eq!(recs[1], ("rB".into(), FLAG_REVERSE | FLAG_DUPLICATE));
}

/// Strand-aware default: two forward-strand single-end reads at the
/// same chr1:100 collide. Sanity check that the default doesn't break
/// the common case.
#[test]
fn strand_aware_two_fwd_orphans_same_position_collide() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);
    assert_eq!(recs[0], ("rA".into(), 0));
    assert_eq!(recs[1], ("rB".into(), FLAG_DUPLICATE));
}

/// picard-approx: when a fully-mapped PE template arrives BEFORE an
/// orphan at the same 5'-coord, the orphan is marked as a duplicate
/// of the pair (the canonical "fragments don't beat pairs" case).
///
/// Under strand-aware, the same input does NOT mark the orphan as a
/// dup — there is no cross-check between the pair table and the
/// orphan table.
#[test]
fn picard_approx_pe_then_orphan_at_same_coord_marks_orphan_dup() {
    let env = TestEnv::new();
    // Pair rA: R1 at chr1:100 (fwd) + R2 at chr1:200 (rev). Standard
    // proper-pair flags. Then orphan rB: a single forward-strand read
    // at chr1:100 with mate unmapped.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        // rB orphan: R1 mapped at chr1:100 forward, R2 unmapped.
        // FLAG 0x1 | 0x8 | 0x40 = 73 (paired + mate-unmapped + first).
        .rec_simple("rB", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rB", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");

    // Under picard-approx, the orphan rB collides with rA's R1 in the
    // fragment table → marked dup.
    let recs_picard =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-approx"]);
    assert_eq!(recs_picard[0], ("rA".into(), 99));
    assert_eq!(recs_picard[1], ("rA".into(), 147));
    assert_eq!(recs_picard[2], ("rB".into(), 73 | FLAG_DUPLICATE));
    assert_eq!(recs_picard[3], ("rB".into(), 133 | FLAG_DUPLICATE));

    // Under strand-aware, the same input does NOT cross-check; rB
    // passes through as non-dup.
    let bam_out2 = env._tmp.path().join("out2.bam");
    let recs_strand =
        run_and_extract_flags(&env.input, &bam_out2, &["--single-end-strategy", "strand-aware"]);
    for (qname, flag) in &recs_strand {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "strand-aware: {qname} (flag {flag}) should NOT be dup-marked (no cross-check)"
        );
    }
}

/// picard-approx is order-sensitive: an orphan that arrives BEFORE its
/// corresponding pair is NOT marked as a duplicate. This is the
/// documented approximation — Picard's coord-sorted design lets it
/// catch both orderings; our streaming design can only catch the
/// pair-first ordering.
#[test]
fn picard_approx_orphan_then_pe_at_same_coord_does_not_mark_orphan_dup() {
    let env = TestEnv::new();
    // Same setup as the previous test but with the orphan FIRST.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Orphan rB first.
        .rec_simple("rB", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rB", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        // Pair rA second.
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-approx"]);
    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "picard-approx (orphan-first): {qname} (flag {flag}) should NOT be dup-marked"
        );
    }
}

/// Picard-approx must still correctly mark plain pair-pair duplicates
/// the same way strand-aware does — the cross-check is additive, not
/// a replacement for the pair table.
#[test]
fn picard_approx_still_marks_pair_pair_duplicates() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("rB", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rB", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-approx"]);
    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    assert_eq!(recs[2], ("rB".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 147 | FLAG_DUPLICATE));
}
