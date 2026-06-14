//! Tests for `--single-end-strategy picard-exact`.
//!
//! picard-exact is the two-pass mode: fully-mapped/unmapped pairs stream
//! straight to the output, while mapped-orphan / single-end "fragment" reads
//! are buffered to a temp BAM and re-processed after the pair pass against a
//! fragment table holding every paired read end. Two properties distinguish
//! it from `picard-approx`:
//!
//!   1. **Order-independence.** An orphan that arrives *before* its matching
//!      pair is still marked as a duplicate of that pair (picard-approx, being
//!      a single streaming pass, misses this ordering).
//!   2. **Output reordering.** Buffered fragments are emitted at the *end* of
//!      the stream, after every pair — not in input order.

mod helpers;

use helpers::*;

/// Look up the flag emitted for the i-th occurrence (0-based) of `qname` in
/// `recs`, which are `(qname, flag)` pairs in output order.
fn nth_flag(recs: &[(String, u16)], qname: &str, n: usize) -> u16 {
    recs.iter()
        .filter(|(q, _)| q == qname)
        .nth(n)
        .unwrap_or_else(|| panic!("no occurrence #{n} of {qname} in output: {recs:?}"))
        .1
}

/// The headline property: an orphan that arrives BEFORE its corresponding
/// pair is still marked as a duplicate of that pair. This is exactly the
/// input `picard-approx` fails to mark (see test_single_end_strategy.rs);
/// picard-exact gets it right because the fragment table is complete before
/// any fragment is checked.
///
/// Also pins the output reordering: the pair (rA) is emitted before the
/// orphan (rB) even though rB came first in the input.
#[test]
fn picard_exact_orphan_before_pair_marks_orphan_dup() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Orphan rB FIRST: R1 mapped at chr1:100 forward, R2 unmapped.
        // FLAG 73 = paired | mate-unmapped | first.
        .rec_simple("rB", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rB", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        // Pair rA SECOND: R1 at chr1:100 (fwd), R2 at chr1:200 (rev).
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);

    // Output reordering: the pair flushes in pass 1, the orphan in pass 2.
    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    // The orphan (both of its records) is marked dup of rA's R1 end.
    assert_eq!(recs[2], ("rB".into(), 73 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 133 | FLAG_DUPLICATE));
}

/// The pair-first ordering must also mark the orphan dup (picard-approx gets
/// this one too; picard-exact must not regress it).
#[test]
fn picard_exact_pair_before_orphan_marks_orphan_dup() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Pair rA first.
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("rA", 147, "chr1", 200, "50M", "=", 100, -150)
        // Orphan rB second.
        .rec_simple("rB", 73, "chr1", 100, "50M", "=", 100, 0)
        .record("rB", 133, "chr1", 100, 0, "*", "=", 100, 0, &"A".repeat(50), &"I".repeat(50))
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);
    assert_eq!(nth_flag(&recs, "rA", 0), 99);
    assert_eq!(nth_flag(&recs, "rB", 0), 73 | FLAG_DUPLICATE);
    assert_eq!(nth_flag(&recs, "rB", 1), 133 | FLAG_DUPLICATE);
}

/// picard-exact must still mark plain pair-pair duplicates the same way the
/// other strategies do — the fragment table is additive, not a replacement
/// for the pair table.
#[test]
fn picard_exact_still_marks_pair_pair_duplicates() {
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
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);
    // No fragments → no reordering; pairs stream in input order.
    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    assert_eq!(recs[2], ("rB".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 147 | FLAG_DUPLICATE));
}

/// Two single-end orphans at the same 5' position with NO pair there: the
/// first is kept, the second is a duplicate (count = group size − 1, matching
/// Picard's fragment handling — the keeper's identity is arbitrary, the
/// count is not).
#[test]
fn picard_exact_two_orphans_no_pair_one_marked_dup() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);
    let dup_count = recs.iter().filter(|(_, f)| f & FLAG_DUPLICATE != 0).count();
    assert_eq!(dup_count, 1, "exactly one of two identical orphans should be dup: {recs:?}");
}

/// A lone single-end orphan with no pair and no other fragment at its
/// position is never a duplicate.
#[test]
fn picard_exact_lone_orphan_not_dup() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rSolo", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);
    assert_eq!(recs[0], ("rSolo".into(), 0));
}

/// A forward orphan does NOT collide with a reverse-strand pair end at the
/// same leftmost coordinate: fragment keying is strand-aware (the pair end's
/// 5' position is on the opposite strand), so the orphan survives.
#[test]
fn picard_exact_fragment_keying_is_strand_aware() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Pair rA: R1 reverse at chr1:100, R2 forward at chr1:400.
        // FLAG R1 = paired|proper|reverse|first = 1|2|16|64 = 83.
        // FLAG R2 = paired|proper|mate-reverse|last = 1|2|32|128 = 163.
        .rec_simple("rA", 83, "chr1", 100, "50M", "=", 400, -350)
        .rec_simple("rA", 163, "chr1", 400, "50M", "=", 100, 350)
        // Orphan rB: forward single-end at chr1:100 — same leftmost as rA's
        // R1, but opposite strand, so a different fragment key.
        .rec_simple("rB", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["--single-end-strategy", "picard-exact"]);
    assert_eq!(
        nth_flag(&recs, "rB", 0) & FLAG_DUPLICATE,
        0,
        "fwd orphan must not match rev pair end"
    );
}

/// `--tmp-dir` is honored: pointing it at a writable directory succeeds and
/// produces the same result as the default temp dir.
#[test]
fn picard_exact_respects_tmp_dir_flag() {
    let env = TestEnv::new();
    let tmp_dir = env._tmp.path().join("scratch");
    std::fs::create_dir(&tmp_dir).expect("create scratch tmp dir");
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("rB", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(
        &env.input,
        &bam_out,
        &["--single-end-strategy", "picard-exact", "--tmp-dir", tmp_dir.to_str().unwrap()],
    );
    let dup_count = recs.iter().filter(|(_, f)| f & FLAG_DUPLICATE != 0).count();
    assert_eq!(dup_count, 1, "result must be unchanged when using a custom --tmp-dir");
}
