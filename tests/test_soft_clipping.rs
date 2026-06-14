//! Generative tests for dupblaster's soft-clip handling.
//!
//! Duplicate detection works on a read's **5'-aligned reference position**,
//! not its raw POS. For a forward-strand read with leading soft-clip
//! `KS<rest>`, the 5' fragment position is `pos - K`. For a reverse-strand
//! read with trailing soft-clip `<rest>KS`, the 5' fragment position is
//! `pos + ref_aligned_len + K - 1`. Two reads whose POS differs but whose
//! 5'-aligned position is identical must be detected as duplicates.

mod helpers;

use helpers::*;

#[test]
fn forward_leading_softclip_is_compensated_in_signature() {
    let env = TestEnv::new();
    // Pair A: R1 = 50M @100 (5'=100), R2 = 50M @151 reverse (5'=200).
    // Pair B: R1 = 5S45M @105 (5'=105-5=100, same), R2 same as A.
    // Pair B must be detected as a duplicate of A.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 151, 100)
        .rec_simple("rA", 147, "chr1", 151, "50M", "=", 100, -100)
        .rec_simple("rB", 99, "chr1", 105, "5S45M", "=", 151, 95)
        .rec_simple("rB", 147, "chr1", 151, "50M", "=", 105, -95)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    // rB's primary records inherit the dup flag.
    assert_eq!(recs[2], ("rB".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 147 | FLAG_DUPLICATE));
}

#[test]
fn reverse_trailing_softclip_is_compensated_in_signature() {
    let env = TestEnv::new();
    // Pair A: R1 = 50M @100 (fwd, 5'=100), R2 = 50M @151 (rev, 5'=200).
    // Pair B: R1 same as A, R2 = 45M5S @151 (rev, 5'=151+45+5-1=200).
    // Both pairs share the same {fwd 5'=100, rev 5'=200} signature.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 151, 100)
        .rec_simple("rA", 147, "chr1", 151, "50M", "=", 100, -100)
        .rec_simple("rB", 99, "chr1", 100, "50M", "=", 151, 100)
        .rec_simple("rB", 147, "chr1", 151, "45M5S", "=", 100, -95)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    assert_eq!(recs[2], ("rB".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 147 | FLAG_DUPLICATE));
}

#[test]
fn different_softclip_lengths_yielding_different_five_prime_are_not_dups() {
    let env = TestEnv::new();
    // Pair A: R1 = 50M @100 (5'=100). Pair B: R1 = 10S40M @100 (5'=90).
    // Same POS, different 5'-aligned position → not duplicates.
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "50M", "=", 151, 100)
        .rec_simple("rA", 147, "chr1", 151, "50M", "=", 100, -100)
        .rec_simple("rB", 99, "chr1", 100, "10S40M", "=", 151, 90)
        .rec_simple("rB", 147, "chr1", 151, "50M", "=", 100, -90)
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
fn softclip_pushing_five_prime_before_contig_start_does_not_panic() {
    let env = TestEnv::new();
    // R1 starts at POS=5 with a 10-base leading soft-clip → 5' fragment
    // position is 5-10 = -5, i.e. before the contig start. dupblaster
    // must process this without panic (the partitioned hash table pads
    // by --max-read-length so the binned position stays non-negative).
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 5, "10S40M", "=", 55, 90)
        .rec_simple("rA", 147, "chr1", 55, "50M", "=", 5, -90)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    // Single pair, no duplicates expected — just verifying the run
    // completes and produces output.
    assert_eq!(recs.len(), 2);
    for (_, flag) in &recs {
        assert_eq!(flag & FLAG_DUPLICATE, 0);
    }
}

#[test]
fn leading_hard_clip_counts_toward_five_prime_clip_like_soft_clip() {
    // Picard / samblaster convention: BOTH leading-H and leading-S are
    // treated as 5' clipping when computing the fragment 5' position.
    // Verify equivalence: 15S35M and 5H10S35M both have 15 bp of leading
    // clip → same 5' position → must be marked as duplicates.
    //
    // R1 query lengths differ between the two records (50 vs 45 bp)
    // because H is not query-consuming, but the signatures still match.
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "15S35M", "=", 151, 85)
        .rec_simple("rA", 147, "chr1", 151, "50M", "=", 100, -85)
        .rec_simple("rB", 99, "chr1", 100, "5H10S35M", "=", 151, 85)
        .rec_simple("rB", 147, "chr1", 151, "50M", "=", 100, -85)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("rA".into(), 99));
    assert_eq!(recs[1], ("rA".into(), 147));
    assert_eq!(recs[2], ("rB".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rB".into(), 147 | FLAG_DUPLICATE));
}

#[test]
fn leading_hard_clip_changes_signature_relative_to_soft_clip_only() {
    // Inverse: a record with H+S leading clip should NOT collide with
    // an otherwise-identical record that has only S leading clip — the
    // total leading clip differs, so the 5' positions differ.
    //   rA: 10S40M @100 → 5' = 100 - 10 = 90
    //   rB: 5H10S40M @100 → 5' = 100 - 15 = 85
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("rA", 99, "chr1", 100, "10S40M", "=", 151, 90)
        .rec_simple("rA", 147, "chr1", 151, "50M", "=", 100, -90)
        .rec_simple("rB", 99, "chr1", 100, "5H10S40M", "=", 151, 90)
        .rec_simple("rB", 147, "chr1", 151, "50M", "=", 100, -90)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) should NOT be marked dup: total leading clips differ"
        );
    }
}
