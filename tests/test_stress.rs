//! Large-input stress test — 11K templates, 22K records, runs in
//! under a second in debug. Verifies exact dup counts on a known
//! workload size that's representative enough to catch off-by-one
//! errors in block boundary handling.

mod helpers;

use helpers::*;

/// Build 10K distinct pairs at evenly-spaced positions, then append 1K
/// pairs that duplicate the first 1K. We expect exactly 1000 duplicate
/// templates marked, and 11000 total templates processed.
#[test]
fn ten_thousand_pairs_with_thousand_dups() {
    let env = TestEnv::new();

    // Build the SAM programmatically. 10000 unique pairs at positions
    // spaced 1000 bp apart; then 1000 duplicates of the first 1000.
    let mut sb = SamBuilder::new().sq("chr1", 20_000_000);
    for i in 0..10_000u32 {
        let pos = 100 + i * 1000;
        let mate = pos + 200;
        let qname = format!("u{i:05}");
        sb = sb
            .rec_simple(&qname, 99, "chr1", pos, "50M", "=", mate, 250)
            .rec_simple(&qname, 147, "chr1", mate, "50M", "=", pos, -250);
    }
    for i in 0..1_000u32 {
        let pos = 100 + i * 1000;
        let mate = pos + 200;
        let qname = format!("d{i:05}");
        sb = sb
            .rec_simple(&qname, 99, "chr1", pos, "50M", "=", mate, 250)
            .rec_simple(&qname, 147, "chr1", mate, "50M", "=", pos, -250);
    }
    sb.write_to(&env.input);

    let bam_out = env._tmp.path().join("out.bam");
    let recs = run_and_extract_flags(&env.input, &bam_out, &[]);

    let total = recs.len();
    let dup_count = recs.iter().filter(|(_, f)| f & FLAG_DUPLICATE != 0).count();
    // 11000 templates * 2 reads/template = 22000 records.
    assert_eq!(total, 22_000, "expected 22000 records, got {total}");
    // 1000 dup templates * 2 reads = 2000 records flagged.
    assert_eq!(dup_count, 2_000, "expected 2000 dup-flagged records, got {dup_count}");
}
