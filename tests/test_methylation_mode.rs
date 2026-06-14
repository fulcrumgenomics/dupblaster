//! Directional methylation-mode (`--methylation-mode directional`) semantics.
//!
//! The defining behavior: two fragments from the **same genomic locus but
//! opposite original strands** (OT/CTOT vs OB/CTOB) must NOT be collapsed as
//! duplicates, because they carry independent methylation — while genuine PCR
//! copies of the *same* strand still must. Standard WGS keying (the default)
//! coordinate-canonicalizes the pair and therefore merges the two strands;
//! directional mode keys in template (first-of-pair → second-of-pair) order,
//! which keeps them distinct.
//!
//! Flag legend for the OT/OB pairs used throughout (all on a single fragment
//! whose left read spans `pos..pos+49` and right read `pos+100..pos+149`):
//!
//! * **OT** (Watson-derived): R1 = left, forward (flag 99); R2 = right,
//!   reverse (flag 147). This is the ordinary `FR` "innie" pair.
//! * **OB** (Crick-derived, *same locus*): R1 = right, reverse (flag 83);
//!   R2 = left, forward (flag 163). Same two alignments as OT, but the
//!   first-of-pair / second-of-pair roles — and strands — are swapped.
//!
//! Coordinate canonicalization maps OT and OB to the *same* key (leftmost end
//! first); template order keeps them apart.

mod helpers;

use helpers::*;

/// Without the flag (standard WGS keying), an OT pair and an OB pair at the
/// same locus canonicalize to the same signature, so the second is a duplicate.
/// This is the behavior directional mode deliberately changes.
#[test]
fn wgs_mode_collapses_opposite_strand_pairs_at_same_locus() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // OT: R1 left-fwd (99), R2 right-rev (147).
        .rec_simple("ot", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot", 147, "chr1", 200, "50M", "=", 100, -150)
        // OB: R1 right-rev (83), R2 left-fwd (163) — same two alignments.
        .rec_simple("ob", 83, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 163, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    assert_eq!(recs[0], ("ot".into(), 99));
    assert_eq!(recs[1], ("ot".into(), 147));
    // OB collapses onto OT under coordinate-canonical (WGS) keying.
    assert_eq!(recs[2], ("ob".into(), 83 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("ob".into(), 163 | FLAG_DUPLICATE));
}

/// The cornerstone: in directional mode the *same* input keeps OT and OB
/// distinct — neither is flagged a duplicate of the other.
#[test]
fn directional_mode_keeps_opposite_strand_pairs_distinct() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("ot", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 83, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 163, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    for (qname, flag) in &recs {
        assert_eq!(
            flag & FLAG_DUPLICATE,
            0,
            "{qname} (flag {flag}) must not be a duplicate — OT and OB are different strands"
        );
    }
}

/// Directional mode must still collapse genuine PCR copies of the *same*
/// strand: two identical OT pairs → the second is a duplicate.
#[test]
fn directional_mode_collapses_same_strand_pcr_duplicates() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("ot1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ot2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot2", 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0], ("ot1".into(), 99));
    assert_eq!(recs[1], ("ot1".into(), 147));
    assert_eq!(recs[2], ("ot2".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("ot2".into(), 147 | FLAG_DUPLICATE));
}

/// Per-strand independence: OT, OT, OB at one locus → exactly the second OT is
/// a duplicate; the OB is untouched.
#[test]
fn directional_mode_three_deep_keeps_other_strand() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("ot1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ot2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 83, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 163, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0, "first OT is the representative");
    assert_eq!(recs[2], ("ot2".into(), 99 | FLAG_DUPLICATE), "second OT is a duplicate");
    assert_eq!(recs[3], ("ot2".into(), 147 | FLAG_DUPLICATE));
    assert_eq!(recs[4].1 & FLAG_DUPLICATE, 0, "OB is a different strand, not a duplicate");
    assert_eq!(recs[5].1 & FLAG_DUPLICATE, 0);
}

// ── orientation spectrum: all four template-order orientations dedup stably ──

/// `FF` (both forward) same-strand pairs: identical copies still collapse.
#[test]
fn directional_mode_collapses_identical_ff_pairs() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // R1 fwd (PAIRED|FIRST=67), R2 fwd (PAIRED|LAST=131).
        .rec_simple("ff1", 67, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("ff1", 131, "chr1", 200, "50M", "=", 100, 0)
        .rec_simple("ff2", 67, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("ff2", 131, "chr1", 200, "50M", "=", 100, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0);
    assert_eq!(recs[2], ("ff2".into(), 67 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("ff2".into(), 131 | FLAG_DUPLICATE));
}

/// `RR` (both reverse) same-strand pairs: identical copies still collapse.
#[test]
fn directional_mode_collapses_identical_rr_pairs() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // R1 rev (PAIRED|REVERSE|MATE_REVERSE|FIRST=115),
        // R2 rev (PAIRED|REVERSE|MATE_REVERSE|LAST=179).
        .rec_simple("rr1", 115, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("rr1", 179, "chr1", 200, "50M", "=", 100, 0)
        .rec_simple("rr2", 115, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("rr2", 179, "chr1", 200, "50M", "=", 100, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0);
    assert_eq!(recs[2], ("rr2".into(), 115 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rr2".into(), 179 | FLAG_DUPLICATE));
}

/// `RF` "outie" pairs (R1 reverse on the left, R2 forward on the right):
/// identical copies still collapse, confirming the outie geometry keys stably.
#[test]
fn directional_mode_collapses_identical_rf_outie_pairs() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // R1 rev-left (PAIRED|REVERSE|FIRST=81), R2 fwd-right (PAIRED|MATE_REVERSE|LAST=161).
        .rec_simple("rf1", 81, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("rf1", 161, "chr1", 200, "50M", "=", 100, 0)
        .rec_simple("rf2", 81, "chr1", 100, "50M", "=", 200, 0)
        .rec_simple("rf2", 161, "chr1", 200, "50M", "=", 100, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0);
    assert_eq!(recs[2], ("rf2".into(), 81 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("rf2".into(), 161 | FLAG_DUPLICATE));
}

/// Cross-contig (chimeric) pairs: a strand-swapped chimeric pair at the same
/// two loci is kept distinct (different strands), while an exact copy of the
/// first chimera still collapses.
#[test]
fn directional_mode_keeps_strand_swapped_cross_contig_pairs_distinct() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .sq("chr2", 1_000_000)
        // Chimera A ("OT"-like): R1 chr1:100 fwd (97), R2 chr2:200 rev (145).
        .rec_simple("a1", 97, "chr1", 100, "50M", "chr2", 200, 0)
        .rec_simple("a1", 145, "chr2", 200, "50M", "chr1", 100, 0)
        // Exact copy of A → duplicate.
        .rec_simple("a2", 97, "chr1", 100, "50M", "chr2", 200, 0)
        .rec_simple("a2", 145, "chr2", 200, "50M", "chr1", 100, 0)
        // Chimera B ("OB"-like): R1 chr2:200 rev (81), R2 chr1:100 fwd (161) —
        // same two alignments, first/second roles and strands swapped.
        .rec_simple("b", 81, "chr2", 200, "50M", "chr1", 100, 0)
        .rec_simple("b", 161, "chr1", 100, "50M", "chr2", 200, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0, "chimera A is the representative");
    assert_eq!(recs[2], ("a2".into(), 97 | FLAG_DUPLICATE), "exact copy of A is a duplicate");
    assert_eq!(recs[3], ("a2".into(), 145 | FLAG_DUPLICATE));
    assert_eq!(recs[4].1 & FLAG_DUPLICATE, 0, "strand-swapped chimera B is kept distinct");
    assert_eq!(recs[5].1 & FLAG_DUPLICATE, 0);
}

// ── composition with the rest of the pipeline ────────────────────────────────

/// A supplementary alignment in a duplicate block inherits the dup flag in
/// directional mode, just as it does for standard keying.
#[test]
fn directional_mode_propagates_dup_flag_to_supplementary() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("ot1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot1", 147, "chr1", 200, "50M", "=", 100, -150)
        // Duplicate OT pair, with a supplementary alignment on R1
        // (PAIRED|MATE_REVERSE|FIRST|SUPPLEMENTARY = 99 | 0x800 = 2147).
        .rec_simple("ot2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ot2", 99 | FLAG_SUPPLEMENTARY, "chr1", 500, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0].1 & FLAG_DUPLICATE, 0, "first OT is the representative");
    // Every record of the duplicate block — including the supplementary.
    assert_eq!(recs[2], ("ot2".into(), 99 | FLAG_DUPLICATE));
    assert_eq!(recs[3], ("ot2".into(), 147 | FLAG_DUPLICATE));
    assert_eq!(recs[4], ("ot2".into(), 99 | FLAG_SUPPLEMENTARY | FLAG_DUPLICATE));
}

/// `--remove-dups` in directional mode drops same-strand copies but keeps both
/// original strands (OT and OB both survive; the OT copy is removed).
#[test]
fn directional_mode_remove_dups_keeps_both_strands() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple("ot1", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot1", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ot2", 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple("ot2", 147, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 83, "chr1", 200, "50M", "=", 100, -150)
        .rec_simple("ob", 163, "chr1", 100, "50M", "=", 200, 150)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs =
        run_and_extract_flags(&env.input, &bam_out, &["-r", "--methylation-mode", "directional"]);

    // ot2 (the same-strand copy) is removed; ot1 and the OB pair remain.
    let names: Vec<&str> = recs.iter().map(|(q, _)| q.as_str()).collect();
    assert_eq!(names, vec!["ot1", "ot1", "ob", "ob"]);
    for (qname, flag) in &recs {
        assert_eq!(flag & FLAG_DUPLICATE, 0, "{qname} survivors carry no dup flag after removal");
    }
}

/// Directional mode composes with single-end / orphan handling, which is
/// untouched by the mode (the fragment path is already strand-aware): two
/// identical single-end reads still dedup.
#[test]
fn directional_mode_still_marks_single_end_duplicates() {
    let env = TestEnv::new();
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        // Unpaired (single-end) forward reads at the same 5' position.
        .rec_simple("se1", 0, "chr1", 100, "50M", "*", 0, 0)
        .rec_simple("se2", 0, "chr1", 100, "50M", "*", 0, 0)
        .write_to(&env.input);
    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &["--methylation-mode", "directional"]);

    assert_eq!(recs[0], ("se1".into(), 0));
    assert_eq!(recs[1], ("se2".into(), FLAG_DUPLICATE));
}
