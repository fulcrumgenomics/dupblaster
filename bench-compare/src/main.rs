//! Single-pass duplicate-marking comparison: Picard (oracle) vs one partner.
//!
//! Rust port of `workflow/scripts/compare_against_picard.py` with
//! byte-identical **TSV** output (the human-readable stdout report is not
//! line-identical — it labels itself `(rust)` and formats numbers differently).
//! It joins two QNAME-sorted BAMs by name, reading each record's qname, flag,
//! and `kf` tag directly off the raw BAM bytes (no SAM-text round-trip), and
//! computes:
//!   - per-category duplicate-set concordance (Picard's kf groups vs partner
//!     dup-mark counts), and
//!   - the four-bucket orphan-discordance triage (cross-table / tiebreaker /
//!     picard-only / partner-only).
//!
//! Picard's canonical key is read off the `kf` tag rather than re-derived. See
//! the Python script for the full rationale.
//!
//! Both inputs must be QNAME-sorted (`samtools sort -n`) with the same record
//! set in the same order; we co-stream primaries and verify QNAMEs match.

use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use bstr::{BString, ByteSlice};
use clap::Parser;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::data::field::Value;
use rustc_hash::{FxHashMap, FxHashSet};

/// The buffered, BGZF-decoded byte stream underneath the BAM reader — the
/// concrete inner type of [`bam::io::Reader`] that `open_bam` produces.
type BamSource = bgzf::io::Reader<BufReader<File>>;

/// Global allocator — mimalloc for faster hash-table operations on large BAMs.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// SAM FLAG bit: template has multiple segments (paired-end).
const FLAG_PAIRED: u16 = 0x1;
/// SAM FLAG bit: segment is unmapped.
const FLAG_UNMAPPED: u16 = 0x4;
/// SAM FLAG bit: PCR or optical duplicate.
const FLAG_DUP: u16 = 0x400;
/// SAM FLAG bit: not primary alignment.
const FLAG_SECONDARY: u16 = 0x100;
/// SAM FLAG bit: supplementary alignment.
const FLAG_SUPPLEMENTARY: u16 = 0x800;
/// Combined mask for records that are neither primary nor unmapped — used to skip
/// sec/supp records during the per-template comparison.
const FLAG_SEC_SUPP: u16 = FLAG_SECONDARY | FLAG_SUPPLEMENTARY;

/// The two-byte BAM tag name for Picard's duplicate key field.
const KF: &[u8; 2] = b"kf";

/// Command-line arguments for the bench-compare binary.
#[derive(Parser)]
#[command(about = "Compare a partner tool's duplicate marking against Picard via the kf tag")]
struct Args {
    /// QNAME-sorted Picard output carrying kf tags.
    #[arg(long)]
    picard_bam: PathBuf,
    /// QNAME-sorted partner output.
    #[arg(long)]
    partner_bam: PathBuf,
    /// Label for the partner tool (written to TSV rows).
    #[arg(long)]
    partner_label: String,
    /// Sample name (written to TSV rows).
    #[arg(long, default_value = "")]
    sample: String,
    /// Append a set-equivalence summary row to this TSV.
    #[arg(long)]
    setcmp_tsv: Option<PathBuf>,
    /// Append an orphan-triage summary row to this TSV.
    #[arg(long)]
    orphan_tsv: Option<PathBuf>,
    /// Append per-tool dup-flag inheritance rows (picard + partner) to this TSV.
    #[arg(long)]
    inheritance_tsv: Option<PathBuf>,
    /// Debug: dump every discordant `pe_both_mapped` set (where Picard's and the
    /// partner's dup-mark counts differ) as `group_key, size, picard_marked,
    /// partner_marked`. The group_key embeds each end's `refidx:unclipped5':strand`,
    /// so the geometric pattern of disagreement is visible. Off by default.
    #[arg(long)]
    dump_pe_discordant: Option<PathBuf>,
}

/// Category names mirror compare_against_picard.py's template_category().
fn template_category(flags: &[u16]) -> &'static str {
    match flags.len() {
        0 => "no_primary",
        1 => {
            let f = flags[0];
            if f & FLAG_UNMAPPED != 0 {
                "unmapped"
            } else if f & FLAG_PAIRED == 0 {
                "single_end"
            } else {
                "orphan"
            }
        }
        2 => {
            let a_unm = flags[0] & FLAG_UNMAPPED != 0;
            let b_unm = flags[1] & FLAG_UNMAPPED != 0;
            if a_unm && b_unm {
                "both_unmapped"
            } else if a_unm || b_unm {
                "orphan"
            } else {
                "pe_both_mapped"
            }
        }
        // rare; matches Python's f"multi_primary_{n}" only in spirit — these
        // never participate in pe/orphan rows so the exact label is moot.
        _ => "multi_primary",
    }
}

/// A template's duplicate-set key from its primary records' kf tags: the sorted
/// combination of the ends (faithful proxy for Picard's canonical pair key) or
/// the single fragment key for an orphan/single-end read.
fn group_key(kfs: &[BString]) -> BString {
    if kfs.len() == 1 {
        return kfs[0].clone();
    }
    let mut refs: Vec<&BString> = kfs.iter().collect();
    refs.sort_unstable();
    let mut out = Vec::with_capacity(refs.iter().map(|k| k.len() + 1).sum());
    for (i, k) in refs.iter().enumerate() {
        if i > 0 {
            out.push(b'|');
        }
        out.extend_from_slice(k.as_bytes());
    }
    BString::from(out)
}

/// Accumulated statistics for a single duplicate set (group of templates sharing
/// the same canonical pair key as reported by Picard's `kf` tag).
struct SetEntry {
    /// Total number of templates in this set.
    size: u32,
    /// Number of templates Picard flagged as duplicate within this set.
    picard_marked: u32,
    /// Number of templates the partner tool flagged as duplicate within this set.
    partner_marked: u32,
    /// Template category (e.g. `"pe_both_mapped"`), or `"mixed"` when the set spans
    /// multiple categories (rare but possible with pathological data).
    category: &'static str,
}

/// Collects per-template dup-marking data from both sides of the comparison and
/// later summarizes it into set-concordance and orphan-triage statistics.
#[derive(Default)]
struct Aggregator {
    /// All duplicate sets seen so far, keyed by their canonical group key.
    sets: FxHashMap<BString, SetEntry>,
    /// The per-end `kf` keys that appear in `pe_both_mapped` templates — used to
    /// classify orphan discordances as cross-table (cat1) vs. other.
    pair_end_kfs: FxHashSet<BString>,
    /// Orphan templates grouped by their single mapped-end `kf` key; each entry
    /// records `(picard_marked, partner_marked)` for that orphan.
    orphan_recs: FxHashMap<BString, Vec<(bool, bool)>>,
}

impl Aggregator {
    /// Record one template's dup-marking outcome from both sides.
    ///
    /// `flags` contains the SAM flags of all primary records for this template;
    /// `kfs` holds the `kf` tag values from those primaries (absent on unmapped
    /// ends); `picard_dup`/`partner_dup` indicate whether either tool flagged the
    /// template as a duplicate.
    fn add_template(
        &mut self,
        flags: &[u16],
        kfs: &[BString],
        picard_dup: bool,
        partner_dup: bool,
    ) {
        let cat = template_category(flags);
        if !kfs.is_empty() {
            let gkey = group_key(kfs);
            match self.sets.get_mut(&gkey) {
                Some(e) => {
                    e.size += 1;
                    e.picard_marked += picard_dup as u32;
                    e.partner_marked += partner_dup as u32;
                    if e.category != cat {
                        e.category = "mixed";
                    }
                }
                None => {
                    self.sets.insert(
                        gkey,
                        SetEntry {
                            size: 1,
                            picard_marked: picard_dup as u32,
                            partner_marked: partner_dup as u32,
                            category: cat,
                        },
                    );
                }
            }
        }
        if cat == "pe_both_mapped" {
            for k in kfs {
                self.pair_end_kfs.insert(k.clone());
            }
        } else if cat == "orphan" && !kfs.is_empty() {
            // Picard writes kf only on mapped primaries, so an orphan template
            // has exactly one kf — the mapped end's fragment key.
            self.orphan_recs.entry(kfs[0].clone()).or_default().push((picard_dup, partner_dup));
        }
    }

    /// Roll up duplicate-set concordance by template category, splitting out
    /// singletons (sets of size 1 can't disagree on set membership).
    fn summarize_sets(&self) -> SetSummary {
        let mut cat_total: FxHashMap<&'static str, u64> = FxHashMap::default();
        let mut cat_concordant: FxHashMap<&'static str, u64> = FxHashMap::default();
        let mut singleton_partner_marked = 0u64;
        let mut singleton_picard_marked = 0u64;
        for e in self.sets.values() {
            if e.size < 2 {
                singleton_picard_marked += e.picard_marked as u64;
                singleton_partner_marked += e.partner_marked as u64;
                continue;
            }
            *cat_total.entry(e.category).or_default() += 1;
            if e.picard_marked == e.partner_marked {
                *cat_concordant.entry(e.category).or_default() += 1;
            }
        }
        SetSummary { cat_total, cat_concordant, singleton_partner_marked, singleton_picard_marked }
    }

    /// Triage orphan-record dup-flag disagreements into the explanatory
    /// buckets (cross-table, tie-breaker, picard-only, partner-only).
    fn triage_orphans(&self) -> OrphanTriage {
        let mut t = OrphanTriage::default();
        for (kf, recs) in &self.orphan_recs {
            let n_orphans = recs.len() as u64;
            let picard_count = recs.iter().filter(|(p, _)| *p).count() as u64;
            let partner_count = recs.iter().filter(|(_, o)| *o).count() as u64;
            let has_pair = self.pair_end_kfs.contains(kf);
            t.total_orphans += n_orphans;
            t.picard_marked += picard_count;
            t.partner_marked += partner_count;

            let is_tiebreaker_pos = n_orphans > 1
                && picard_count == partner_count
                && picard_count > 0
                && picard_count < n_orphans
                && !has_pair;

            for &(picard_marked, partner_marked) in recs {
                if picard_marked == partner_marked {
                    if picard_marked {
                        t.concordant_dup += 1;
                    } else {
                        t.concordant_nondup += 1;
                    }
                    continue;
                }
                if picard_marked && !partner_marked {
                    if has_pair {
                        t.cat1_cross_table += 1;
                    } else if is_tiebreaker_pos {
                        t.cat2_tiebreaker += 1;
                    } else {
                        t.cat3_picard_only_other += 1;
                    }
                } else if is_tiebreaker_pos {
                    t.cat2_tiebreaker += 1;
                } else {
                    t.cat4_partner_only_other += 1;
                }
            }
        }
        t.discordant = t.cat1_cross_table
            + t.cat2_tiebreaker
            + t.cat3_picard_only_other
            + t.cat4_partner_only_other;
        t
    }
}

/// Ordered list of set categories for the stdout report. Categories not in
/// this list still appear but are appended after these in an unspecified order.
const SET_CATS: [&str; 6] =
    ["pe_both_mapped", "orphan", "single_end", "unmapped", "both_unmapped", "mixed"];

/// Rolled-up set-equivalence concordance statistics produced by
/// [`Aggregator::summarize_sets`].
struct SetSummary {
    /// Total number of duplicate sets (size >= 2) per template category.
    cat_total: FxHashMap<&'static str, u64>,
    /// Number of concordant sets (Picard and partner marked the same count)
    /// per template category.
    cat_concordant: FxHashMap<&'static str, u64>,
    /// Templates that are Picard singletons (no set) but the partner marked them
    /// as duplicates — indicates extra calls by the partner.
    singleton_partner_marked: u64,
    /// Templates that are Picard singletons but Picard itself marked them; should
    /// always be zero and is reported as a sanity check.
    singleton_picard_marked: u64,
}

/// Four-bucket triage of orphan dup-flag disagreements, produced by
/// [`Aggregator::triage_orphans`]. Mirrors the categories defined in
/// `compare_against_picard.py`.
#[derive(Default)]
struct OrphanTriage {
    /// Total number of orphan primary records examined.
    total_orphans: u64,
    /// Orphans Picard flagged as duplicate.
    picard_marked: u64,
    /// Orphans the partner tool flagged as duplicate.
    partner_marked: u64,
    /// Orphans both tools agreed are duplicate.
    concordant_dup: u64,
    /// Orphans both tools agreed are not duplicate.
    concordant_nondup: u64,
    /// Total discordant orphans (sum of the four cat* fields).
    discordant: u64,
    /// Cat 1: Picard marks an orphan whose mapped-end key also appears in a
    /// paired set — the partner missed a cross-table duplicate.
    cat1_cross_table: u64,
    /// Cat 2: discordance explained by tiebreaking (multiple orphans share the
    /// same key and each tool picks a different one to keep).
    cat2_tiebreaker: u64,
    /// Cat 3: Picard marks the orphan for a reason other than cross-table or
    /// tiebreaking, and the partner does not.
    cat3_picard_only_other: u64,
    /// Cat 4: the partner marks the orphan but Picard does not, for a reason
    /// other than tiebreaking (e.g. samblaster strand-drop edge case).
    cat4_partner_only_other: u64,
}

/// Open a BAM at `path` and consume its header, returning a ready-to-read
/// noodles [`bam::io::Reader`].
fn open_bam(path: &Path) -> Result<bam::io::Reader<BamSource>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut reader = bam::io::Reader::new(BufReader::with_capacity(1 << 20, file));
    reader.read_header().context("reading BAM header")?;
    Ok(reader)
}

/// The fields of one record we keep while grouping: its FLAG and (for primary
/// mapped records) its kf tag. Duplicate/secondary/supplementary status is
/// derived from the flag; sec/supp records carry no kf.
struct RecView {
    /// SAM FLAG field for this record.
    flag: u16,
    /// Picard's `kf` duplicate-key tag, present only on mapped primary records.
    kf: Option<BString>,
}

/// The read name of a BAM record as raw bytes (empty if absent).
fn name_bytes(rec: &bam::Record) -> &[u8] {
    rec.name().map_or(&[][..], |n| n.as_ref())
}

/// Extract the fields of interest from a BAM record into a [`RecView`].
fn extract(rec: &bam::Record) -> RecView {
    let kf = match rec.data().get(KF) {
        Some(Ok(Value::String(s))) => Some(s.to_owned()),
        _ => None,
    };
    RecView { flag: u16::from(rec.flags()), kf }
}

/// Streams a QNAME-grouped BAM one whole group at a time — all records of a read
/// name (primary, secondary, and supplementary) together. Requires the input to
/// be QNAME-grouped (e.g. `samtools sort -n`); a group boundary is a change in
/// read name. Reads one record ahead to detect that boundary.
struct GroupReader<R: Read> {
    /// The underlying BAM reader.
    reader: bam::io::Reader<R>,
    /// Scratch record reused for each `read_record` call to avoid per-read allocation.
    rec: bam::Record,
    /// The first record of the next group, already read to detect the boundary.
    pending: Option<(Vec<u8>, RecView)>,
    /// Set to `true` after the very first read so we can distinguish "never started"
    /// from "EOF reached" when `pending` is `None`.
    started: bool,
}

impl<R: Read> GroupReader<R> {
    /// Construct a `GroupReader` from an already-header-consumed [`bam::io::Reader`].
    fn new(reader: bam::io::Reader<R>) -> Self {
        Self { reader, rec: bam::Record::default(), pending: None, started: false }
    }

    /// Return the next QNAME group as `(qname, records)`, or `None` at EOF.
    /// noodles `read_record` returns `0` at EOF.
    fn next_group(&mut self) -> Result<Option<(Vec<u8>, Vec<RecView>)>> {
        let (qname, mut recs) = match self.pending.take() {
            Some((qn, rv)) => (qn, vec![rv]),
            None => {
                // pending is None only before the first read, or once EOF is reached.
                if self.started {
                    return Ok(None);
                }
                self.started = true;
                if self.reader.read_record(&mut self.rec).context("reading record")? == 0 {
                    return Ok(None);
                }
                (name_bytes(&self.rec).to_vec(), vec![extract(&self.rec)])
            }
        };
        loop {
            if self.reader.read_record(&mut self.rec).context("reading record")? == 0 {
                return Ok(Some((qname, recs))); // EOF: final group complete, pending stays None
            }
            if name_bytes(&self.rec) == qname.as_slice() {
                recs.push(extract(&self.rec));
            } else {
                self.pending = Some((name_bytes(&self.rec).to_vec(), extract(&self.rec)));
                return Ok(Some((qname, recs)));
            }
        }
    }
}

/// Append a tab-separated `row` to the TSV at `path`, writing the `header` line
/// first only if the file is new or empty. Creates the file if it doesn't exist.
fn append_row(path: &Path, header: &[&str], row: &[String]) -> Result<()> {
    let existed = path.metadata().map(|m| m.len() > 0).unwrap_or(false);
    let mut fh = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .with_context(|| format!("opening {} for append", path.display()))?;
    if !existed {
        writeln!(fh, "{}", header.join("\t"))?;
    }
    writeln!(fh, "{}", row.join("\t"))?;
    Ok(())
}

/// Format `num/den * 100` as a six-decimal string, or `"0"` when `den == 0`.
fn pct(num: u64, den: u64) -> String {
    if den == 0 { "0".to_string() } else { format!("{:.6}", 100.0 * num as f64 / den as f64) }
}

/// Per-tool dup-flag inheritance: a QNAME group is "consistent" iff all its
/// records (primary + secondary + supplementary) share the same duplicate flag.
/// Mirrors the old check_inheritance.py; computed per side as a byproduct of the
/// comparison pass. `inconsistency_pct` is reported as inconsistent/with_supp.
#[derive(Default)]
struct Inheritance {
    /// Total QNAME groups (templates) observed.
    total_qnames: u64,
    /// Groups that contain at least one supplementary record — only these can be
    /// inconsistent, so `inconsistency_pct` is reported relative to this count.
    groups_with_supp: u64,
    /// Groups where every record shares the same duplicate flag value.
    consistent: u64,
    /// Groups where primary and sec/supp records disagree on the duplicate flag.
    inconsistent: u64,
}

impl Inheritance {
    /// Record the inheritance consistency of one QNAME group's full record set.
    fn observe(&mut self, recs: &[RecView]) {
        self.total_qnames += 1;
        let mut has_supp = false;
        let mut seen_dup = false;
        let mut seen_nondup = false;
        for r in recs {
            if r.flag & FLAG_SEC_SUPP != 0 {
                has_supp = true;
            }
            if r.flag & FLAG_DUP != 0 {
                seen_dup = true;
            } else {
                seen_nondup = true;
            }
        }
        if has_supp {
            self.groups_with_supp += 1;
        }
        if seen_dup && seen_nondup {
            self.inconsistent += 1;
        } else {
            self.consistent += 1;
        }
    }
}

/// Everything one comparison pass produces.
struct Comparison {
    /// Per-template dup-marking outcomes accumulated during the pass.
    agg: Aggregator,
    /// Dup-flag inheritance statistics for the Picard side.
    picard_inh: Inheritance,
    /// Dup-flag inheritance statistics for the partner tool side.
    partner_inh: Inheritance,
    /// Total number of primary records seen (sec/supp excluded).
    n_primary: u64,
    /// Total number of QNAME groups (templates) processed.
    n_groups: u64,
}

/// Co-stream two QNAME-grouped BAMs one QNAME group at a time, matching groups by
/// name. For each group we (a) record per-side dup-flag inheritance over all its
/// records, and (b) feed the group's *primary* records into the Aggregator for
/// the set/orphan comparison. Group-level lockstep tolerates the two tools
/// differing in secondary/supplementary content; only the QNAME *sequence* must
/// match (both inputs QNAME-sorted the same way). Errors if the sequences diverge.
fn compare<R: Read>(
    picard: &mut GroupReader<R>,
    partner: &mut GroupReader<R>,
) -> Result<Comparison> {
    let mut c = Comparison {
        agg: Aggregator::default(),
        picard_inh: Inheritance::default(),
        partner_inh: Inheritance::default(),
        n_primary: 0,
        n_groups: 0,
    };
    // Reused per-group scratch for the comparison (primary records only).
    let mut flags: Vec<u16> = Vec::new();
    let mut kfs: Vec<BString> = Vec::new();

    loop {
        match (picard.next_group()?, partner.next_group()?) {
            (None, None) => break,
            (Some((pq, precs)), Some((oq, orecs))) => {
                if pq != oq {
                    bail!(
                        "QNAME mismatch at group {}: picard={} partner={}",
                        c.n_groups,
                        pq.as_bstr(),
                        oq.as_bstr()
                    );
                }
                c.n_groups += 1;
                c.picard_inh.observe(&precs);
                c.partner_inh.observe(&orecs);

                // Comparison: primary records only.
                flags.clear();
                kfs.clear();
                let mut picard_dup = false;
                for r in &precs {
                    if r.flag & FLAG_SEC_SUPP != 0 {
                        continue;
                    }
                    flags.push(r.flag);
                    if let Some(kf) = &r.kf {
                        kfs.push(kf.clone());
                    }
                    if r.flag & FLAG_DUP != 0 {
                        picard_dup = true;
                    }
                    c.n_primary += 1;
                }
                let mut partner_dup = false;
                for r in &orecs {
                    if r.flag & FLAG_SEC_SUPP == 0 && r.flag & FLAG_DUP != 0 {
                        partner_dup = true;
                    }
                }
                c.agg.add_template(&flags, &kfs, picard_dup, partner_dup);

                if c.n_groups.is_multiple_of(20_000_000) {
                    eprintln!("  ... {} groups, {} primary recs", c.n_groups, c.n_primary);
                }
            }
            (Some((pq, _)), None) => {
                bail!("EOF mismatch: partner ended, picard still at group {}", pq.as_bstr())
            }
            (None, Some((oq, _))) => {
                bail!("EOF mismatch: picard ended, partner still at group {}", oq.as_bstr())
            }
        }
    }
    Ok(c)
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut picard = GroupReader::new(open_bam(&args.picard_bam)?);
    let mut partner = GroupReader::new(open_bam(&args.partner_bam)?);
    let cmp = compare(&mut picard, &mut partner)?;

    if cmp.agg.sets.is_empty() {
        bail!("no kf tags found on the Picard BAM — was it produced with TAG_DUPLICATE_KEY=true?");
    }

    let sets = cmp.agg.summarize_sets();
    let orphans = cmp.agg.triage_orphans();
    let label = &args.partner_label;
    let n_recs = cmp.n_primary;
    let n_templates = cmp.n_groups;

    // ----- report (stdout) -----
    println!("\n=== Picard vs {label}: set-equivalence concordance (rust) ===");
    println!("primary records: {n_recs:>14}");
    println!("templates:       {n_templates:>14}");
    println!(
        "\n{:16}  {:>10}  {:>12}  {:>12}  {:>11}",
        "category", "sets", "concordant", "discordant", "concord_pct"
    );
    let mut seen: Vec<&str> =
        SET_CATS.iter().copied().filter(|c| sets.cat_total.contains_key(c)).collect();
    for c in sets.cat_total.keys() {
        if !seen.contains(c) {
            seen.push(c);
        }
    }
    let (mut gt, mut gc) = (0u64, 0u64);
    for cat in &seen {
        let t = *sets.cat_total.get(cat).unwrap_or(&0);
        let c = *sets.cat_concordant.get(cat).unwrap_or(&0);
        println!("{:16}  {:>10}  {:>12}  {:>12}  {:>10}%", cat, t, c, t - c, pct(c, t));
        gt += t;
        gc += c;
    }
    if gt > 0 {
        println!("{:16}  {:>10}  {:>12}  {:>12}  {:>10}%", "TOTAL", gt, gc, gt - gc, pct(gc, gt));
    }
    println!(
        "\nPicard singletons, {label} dup-marked (extra calls): {}",
        sets.singleton_partner_marked
    );
    println!(
        "Picard singletons, Picard dup-marked (should be 0):  {}",
        sets.singleton_picard_marked
    );

    println!("\n=== Picard vs {label}: orphan discordance triage (rust) ===");
    for (k, v) in [
        ("total_orphans", orphans.total_orphans),
        ("picard_marked", orphans.picard_marked),
        ("partner_marked", orphans.partner_marked),
        ("concordant_dup", orphans.concordant_dup),
        ("concordant_nondup", orphans.concordant_nondup),
        ("discordant", orphans.discordant),
        ("cat1_cross_table", orphans.cat1_cross_table),
        ("cat2_tiebreaker", orphans.cat2_tiebreaker),
        ("cat3_picard_only_other", orphans.cat3_picard_only_other),
        ("cat4_partner_only_other", orphans.cat4_partner_only_other),
    ] {
        println!("  {k:24} {v:>12}");
    }

    // ----- TSV rows (byte-identical to compare_against_picard.py) -----
    if let Some(p) = &args.setcmp_tsv {
        let pe_t = *sets.cat_total.get("pe_both_mapped").unwrap_or(&0);
        let pe_c = *sets.cat_concordant.get("pe_both_mapped").unwrap_or(&0);
        let or_t = *sets.cat_total.get("orphan").unwrap_or(&0);
        let or_c = *sets.cat_concordant.get("orphan").unwrap_or(&0);
        append_row(
            p,
            &[
                "sample",
                "partner",
                "pe_sets",
                "pe_concordant",
                "pe_concord_pct",
                "orphan_sets",
                "orphan_concordant",
                "orphan_concord_pct",
                "picard_singletons_partner_marked",
            ],
            &[
                args.sample.clone(),
                label.clone(),
                pe_t.to_string(),
                pe_c.to_string(),
                pct(pe_c, pe_t),
                or_t.to_string(),
                or_c.to_string(),
                pct(or_c, or_t),
                sets.singleton_partner_marked.to_string(),
            ],
        )?;
    }
    if let Some(p) = &args.orphan_tsv {
        append_row(
            p,
            &[
                "sample",
                "partner",
                "total_orphans",
                "picard_marked",
                "partner_marked",
                "concordant_dup",
                "concordant_nondup",
                "discordant",
                "cat1_cross_table",
                "cat2_tiebreaker",
                "cat3_picard_only_other",
                "cat4_partner_only_other",
            ],
            &[
                args.sample.clone(),
                label.clone(),
                orphans.total_orphans.to_string(),
                orphans.picard_marked.to_string(),
                orphans.partner_marked.to_string(),
                orphans.concordant_dup.to_string(),
                orphans.concordant_nondup.to_string(),
                orphans.discordant.to_string(),
                orphans.cat1_cross_table.to_string(),
                orphans.cat2_tiebreaker.to_string(),
                orphans.cat3_picard_only_other.to_string(),
                orphans.cat4_partner_only_other.to_string(),
            ],
        )?;
    }

    // ----- dup-flag inheritance (folds in check_inheritance.py) -----
    let inh_row = |tool: &str, inh: &Inheritance| {
        println!(
            "inheritance {tool:16} total={} with_supp={} consistent={} inconsistent={}",
            inh.total_qnames, inh.groups_with_supp, inh.consistent, inh.inconsistent
        );
    };
    inh_row("picard", &cmp.picard_inh);
    inh_row(label, &cmp.partner_inh);
    if let Some(p) = &args.inheritance_tsv {
        // One row for the tool under test (the partner). Picard is the oracle and
        // its inheritance is identical across every comparison, so we don't repeat
        // it per-partner; it is shown on stdout above for visibility.
        append_row(
            p,
            &[
                "sample",
                "tool",
                "total_qnames",
                "groups_with_supp",
                "consistent",
                "inconsistent",
                "inconsistency_pct",
            ],
            &[
                args.sample.clone(),
                label.clone(),
                cmp.partner_inh.total_qnames.to_string(),
                cmp.partner_inh.groups_with_supp.to_string(),
                cmp.partner_inh.consistent.to_string(),
                cmp.partner_inh.inconsistent.to_string(),
                pct(cmp.partner_inh.inconsistent, cmp.partner_inh.groups_with_supp),
            ],
        )?;
    }

    if let Some(p) = &args.dump_pe_discordant {
        let mut fh =
            std::fs::File::create(p).with_context(|| format!("creating {}", p.display()))?;
        writeln!(fh, "group_key\tsize\tpicard_marked\tpartner_marked")?;
        for (gkey, e) in &cmp.agg.sets {
            if e.size >= 2 && e.category == "pe_both_mapped" && e.picard_marked != e.partner_marked
            {
                writeln!(
                    fh,
                    "{}\t{}\t{}\t{}",
                    gkey.as_bstr(),
                    e.size,
                    e.picard_marked,
                    e.partner_marked
                )?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bs(s: &str) -> BString {
        BString::from(s.as_bytes())
    }

    const MAPPED: u16 = FLAG_PAIRED; // paired + mapped (no UNMAPPED bit)
    const MATE_UNMAPPED: u16 = FLAG_PAIRED | FLAG_UNMAPPED;

    // ---- template_category: every branch ----

    #[test]
    fn template_category_no_primary_for_empty() {
        assert_eq!(template_category(&[]), "no_primary");
    }

    #[test]
    fn template_category_single_record_branches() {
        assert_eq!(template_category(&[FLAG_UNMAPPED]), "unmapped");
        assert_eq!(template_category(&[0]), "single_end"); // mapped, not paired
        assert_eq!(template_category(&[FLAG_PAIRED]), "orphan"); // paired+mapped, lone primary
    }

    #[test]
    fn template_category_pair_branches() {
        assert_eq!(template_category(&[MAPPED, MAPPED]), "pe_both_mapped");
        assert_eq!(template_category(&[MAPPED, MATE_UNMAPPED]), "orphan");
        assert_eq!(template_category(&[MATE_UNMAPPED, MAPPED]), "orphan");
        assert_eq!(template_category(&[MATE_UNMAPPED, MATE_UNMAPPED]), "both_unmapped");
    }

    #[test]
    fn template_category_more_than_two_is_multi_primary() {
        assert_eq!(template_category(&[MAPPED, MAPPED, MAPPED]), "multi_primary");
    }

    // ---- group_key ----

    #[test]
    fn group_key_is_single_for_one_end() {
        assert_eq!(group_key(&[bs("3:100:F")]), bs("3:100:F"));
    }

    #[test]
    fn group_key_pair_is_sorted_and_order_independent() {
        let ab = group_key(&[bs("a"), bs("b")]);
        let ba = group_key(&[bs("b"), bs("a")]);
        assert_eq!(ab, bs("a|b"));
        assert_eq!(ab, ba);
    }

    #[test]
    fn group_key_collapses_opposite_strands_at_same_coord_deterministically() {
        // Picard force-collapses RF->FR at identical coordinates. Sorting the two
        // kf strings reproduces that *because* 'F' < 'R', so the key is identical
        // regardless of which end was read1. This is the load-bearing coincidence
        // the pair-key proxy relies on — pin it down.
        let fr = group_key(&[bs("5:100:F"), bs("5:100:R")]);
        let rf = group_key(&[bs("5:100:R"), bs("5:100:F")]);
        assert_eq!(fr, rf);
        assert_eq!(fr, bs("5:100:F|5:100:R"));
        // Same-strand pairs at the same coordinate stay distinct from the FR case
        // and from each other.
        assert_ne!(group_key(&[bs("5:100:F"), bs("5:100:F")]), fr);
        assert_ne!(group_key(&[bs("5:100:R"), bs("5:100:R")]), fr);
        assert_ne!(
            group_key(&[bs("5:100:F"), bs("5:100:F")]),
            group_key(&[bs("5:100:R"), bs("5:100:R")])
        );
    }

    #[test]
    fn group_key_distinguishes_different_positions() {
        assert_ne!(
            group_key(&[bs("5:100:F"), bs("5:200:R")]),
            group_key(&[bs("5:100:F"), bs("5:201:R")])
        );
    }

    // ---- summarize_sets ----

    #[test]
    fn summarize_concordant_discordant_and_mixed_pe_sets() {
        let mut agg = Aggregator::default();
        // Concordant pe set (size 2): both tools mark the same count (1 each).
        agg.add_template(&[MAPPED, MAPPED], &[bs("A"), bs("B")], false, false);
        agg.add_template(&[MAPPED, MAPPED], &[bs("A"), bs("B")], true, true);
        // Discordant pe set (size 2): Picard marks 1, partner marks 0.
        agg.add_template(&[MAPPED, MAPPED], &[bs("C"), bs("D")], false, false);
        agg.add_template(&[MAPPED, MAPPED], &[bs("C"), bs("D")], true, false);
        let s = agg.summarize_sets();
        assert_eq!(*s.cat_total.get("pe_both_mapped").unwrap(), 2);
        assert_eq!(*s.cat_concordant.get("pe_both_mapped").unwrap(), 1);
    }

    #[test]
    fn summarize_singletons_are_excluded_from_sets_and_counted_separately() {
        let mut agg = Aggregator::default();
        // A lone pe template (size 1) where the partner over-marked: not a "set",
        // counted as a partner extra call on a Picard singleton.
        agg.add_template(&[MAPPED, MAPPED], &[bs("X"), bs("Y")], false, true);
        // A lone pe template Picard marked (sanity bucket).
        agg.add_template(&[MAPPED, MAPPED], &[bs("P"), bs("Q")], true, false);
        let s = agg.summarize_sets();
        assert!(s.cat_total.is_empty()); // nothing reached size >= 2
        assert_eq!(s.singleton_partner_marked, 1);
        assert_eq!(s.singleton_picard_marked, 1);
    }

    #[test]
    fn summarize_marks_category_mixed_when_a_set_spans_categories() {
        let mut agg = Aggregator::default();
        // Two templates share a single-end key but differ in category (single_end
        // vs orphan) -> the set is tagged "mixed".
        agg.add_template(&[0], &[bs("K")], false, false); // single_end
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("K")], false, false); // orphan, same key
        let s = agg.summarize_sets();
        assert_eq!(*s.cat_total.get("mixed").unwrap(), 1);
    }

    // ---- triage_orphans: every bucket ----

    #[test]
    fn orphan_concordant_dup_and_nondup() {
        let mut agg = Aggregator::default();
        // Two orphans share a key: both tools mark both (concordant dup x2).
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("D")], true, true);
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("D")], true, true);
        // One orphan neither tool marks (concordant nondup).
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("N")], false, false);
        let t = agg.triage_orphans();
        assert_eq!(t.total_orphans, 3);
        assert_eq!(t.concordant_dup, 2);
        assert_eq!(t.concordant_nondup, 1);
        assert_eq!(t.discordant, 0);
    }

    #[test]
    fn orphan_cat1_cross_table() {
        let mut agg = Aggregator::default();
        agg.add_template(&[MAPPED, MAPPED], &[bs("A"), bs("B")], false, false); // pair registers A,B
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("A")], true, false); // orphan@A, picard-only
        let t = agg.triage_orphans();
        assert_eq!(t.total_orphans, 1);
        assert_eq!(t.picard_marked, 1);
        assert_eq!(t.partner_marked, 0);
        assert_eq!(t.cat1_cross_table, 1);
        assert_eq!(t.discordant, 1);
        assert_eq!(
            (t.cat2_tiebreaker, t.cat3_picard_only_other, t.cat4_partner_only_other),
            (0, 0, 0)
        );
    }

    #[test]
    fn orphan_cat2_tiebreaker() {
        let mut agg = Aggregator::default();
        // No pair at key X; two orphans, each tool marks one but a different keeper.
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("X")], true, false);
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("X")], false, true);
        let t = agg.triage_orphans();
        assert_eq!(t.total_orphans, 2);
        assert_eq!(t.cat2_tiebreaker, 2);
        assert_eq!(
            (t.cat1_cross_table, t.cat3_picard_only_other, t.cat4_partner_only_other),
            (0, 0, 0)
        );
    }

    #[test]
    fn orphan_cat3_picard_only_other() {
        let mut agg = Aggregator::default();
        // Lone orphan at a key with no pair: Picard marks, partner doesn't, and it
        // is not a tiebreaker (n == 1) -> picard-only-other.
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("Z")], true, false);
        let t = agg.triage_orphans();
        assert_eq!(t.discordant, 1);
        assert_eq!(t.cat3_picard_only_other, 1);
        assert_eq!((t.cat1_cross_table, t.cat2_tiebreaker, t.cat4_partner_only_other), (0, 0, 0));
    }

    #[test]
    fn orphan_cat4_partner_only_other() {
        let mut agg = Aggregator::default();
        // Lone orphan, partner marks but Picard doesn't, no pair, not a tiebreaker
        // -> partner-only-other (e.g. samblaster strand-drop).
        agg.add_template(&[MAPPED, MATE_UNMAPPED], &[bs("W")], false, true);
        let t = agg.triage_orphans();
        assert_eq!(t.discordant, 1);
        assert_eq!(t.cat4_partner_only_other, 1);
        assert_eq!((t.cat1_cross_table, t.cat2_tiebreaker, t.cat3_picard_only_other), (0, 0, 0));
    }

    // ---- pct ----

    #[test]
    fn pct_zero_denominator_is_literal_zero() {
        assert_eq!(pct(0, 0), "0");
        assert_eq!(pct(5, 0), "0");
    }

    #[test]
    fn pct_formats_six_decimals() {
        assert_eq!(pct(1, 2), "50.000000");
        assert_eq!(pct(1997800, 2000367), "99.871674");
    }

    // ---- append_row ----

    #[test]
    fn append_row_writes_header_only_when_file_is_empty() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let header = ["a", "b"];
        append_row(&path, &header, &["1".into(), "2".into()]).unwrap();
        append_row(&path, &header, &["3".into(), "4".into()]).unwrap();
        let contents = std::fs::read_to_string(&path).unwrap();
        assert_eq!(contents, "a\tb\n1\t2\n3\t4\n");
    }

    // ---- Inheritance::observe ----

    fn rv(flag: u16, kf: Option<&[u8]>) -> RecView {
        RecView { flag, kf: kf.map(BString::from) }
    }

    #[test]
    fn inheritance_consistent_when_all_records_share_dup_flag() {
        let mut inh = Inheritance::default();
        // Group 1: primary + supplementary, both dup -> consistent, has supp.
        inh.observe(&[
            rv(FLAG_PAIRED | FLAG_DUP, Some(b"k")),
            rv(FLAG_SUPPLEMENTARY | FLAG_DUP, None),
        ]);
        // Group 2: two primaries, both non-dup, no supp -> consistent.
        inh.observe(&[rv(FLAG_PAIRED, Some(b"k")), rv(FLAG_PAIRED, Some(b"k"))]);
        assert_eq!(
            (inh.total_qnames, inh.groups_with_supp, inh.consistent, inh.inconsistent),
            (2, 1, 2, 0)
        );
    }

    #[test]
    fn inheritance_inconsistent_when_dup_flag_differs_within_group() {
        let mut inh = Inheritance::default();
        // Primary marked dup but its supplementary is not -> inconsistent.
        inh.observe(&[rv(FLAG_PAIRED | FLAG_DUP, Some(b"k")), rv(FLAG_SUPPLEMENTARY, None)]);
        assert_eq!(
            (inh.total_qnames, inh.groups_with_supp, inh.consistent, inh.inconsistent),
            (1, 1, 0, 1)
        );
    }

    // ---- end-to-end I/O: build tiny BAMs and run GroupReader + compare ----

    use noodles_sam::Header;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value as BufValue;

    const DUP: u16 = 0x400;
    const SUPP: u16 = FLAG_PAIRED | 0x800;

    /// Build an unmapped [`RecordBuf`] with the given read name, flags, and
    /// optional kf tag. Left unmapped (no reference) so it writes cleanly
    /// against an empty header; only name/flags/kf matter to the comparison.
    fn rec(name: &[u8], flags: u16, kf: Option<&[u8]>) -> RecordBuf {
        let mut rb = RecordBuf::default();
        *rb.name_mut() = Some(name.into());
        *rb.flags_mut() = Flags::from(flags);
        if let Some(k) = kf {
            rb.data_mut().insert(Tag::new(b'k', b'f'), BufValue::String(BString::from(k)));
        }
        rb
    }

    /// Write records into a tiny BAM (empty header) at `path` via noodles.
    fn write_bam(path: &Path, recs: &[RecordBuf]) {
        let header = Header::default();
        let mut w = bam::io::Writer::new(File::create(path).unwrap());
        w.write_header(&header).unwrap();
        for r in recs {
            w.write_alignment_record(&header, r).unwrap();
        }
        w.try_finish().unwrap();
    }

    fn group_reader(path: &Path) -> GroupReader<BamSource> {
        GroupReader::new(open_bam(path).unwrap())
    }

    #[test]
    fn group_reader_yields_whole_qname_groups() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("g.bam");
        write_bam(
            &p,
            &[
                rec(b"a", FLAG_PAIRED, Some(b"k1")),
                rec(b"a", FLAG_PAIRED, Some(b"k2")),
                rec(b"a", SUPP, None),
                rec(b"b", FLAG_PAIRED, Some(b"k3")),
            ],
        );
        let mut gr = group_reader(&p);
        let (q1, r1) = gr.next_group().unwrap().unwrap();
        assert_eq!(q1.as_slice(), b"a");
        assert_eq!(r1.len(), 3);
        let (q2, r2) = gr.next_group().unwrap().unwrap();
        assert_eq!(q2.as_slice(), b"b");
        assert_eq!(r2.len(), 1);
        assert!(gr.next_group().unwrap().is_none());
        assert!(gr.next_group().unwrap().is_none()); // idempotent at EOF
    }

    #[test]
    fn compare_end_to_end_over_real_bams() {
        let dir = tempfile::tempdir().unwrap();
        let pic = dir.path().join("picard.bam");
        let par = dir.path().join("partner.bam");

        let k1: &[u8] = b"0:100:F";
        let k2: &[u8] = b"0:200:R";
        // Picard: two identical pairs (a concordant pe dup set of size 2), an
        // orphan whose mapped end's kf == a pair end (cross-table dup Picard
        // catches), and a supplementary record (excluded from the comparison).
        write_bam(
            &pic,
            &[
                rec(b"pairX", FLAG_PAIRED, Some(k1)),
                rec(b"pairX", FLAG_PAIRED, Some(k2)),
                rec(b"pairX", SUPP, None), // supplementary: counts for inheritance, not comparison
                rec(b"pairY", FLAG_PAIRED | DUP, Some(k1)),
                rec(b"pairY", FLAG_PAIRED | DUP, Some(k2)),
                rec(b"orf", FLAG_PAIRED | DUP, Some(k1)), // mapped end, cross-table dup
                rec(b"orf", FLAG_PAIRED | FLAG_UNMAPPED | DUP, None), // unmapped mate
            ],
        );
        // Partner: same primaries/order, but it does NOT mark the orphan.
        write_bam(
            &par,
            &[
                rec(b"pairX", FLAG_PAIRED, None),
                rec(b"pairX", FLAG_PAIRED, None),
                rec(b"pairX", SUPP, None),
                rec(b"pairY", FLAG_PAIRED | DUP, None),
                rec(b"pairY", FLAG_PAIRED | DUP, None),
                rec(b"orf", FLAG_PAIRED, None),
                rec(b"orf", FLAG_PAIRED | FLAG_UNMAPPED, None),
            ],
        );

        let cmp = compare(&mut group_reader(&pic), &mut group_reader(&par)).unwrap();

        // 6 primary records across 3 groups (the supplementary is not a primary).
        assert_eq!(cmp.n_primary, 6);
        assert_eq!(cmp.n_groups, 3);

        let s = cmp.agg.summarize_sets();
        assert_eq!(*s.cat_total.get("pe_both_mapped").unwrap(), 1);
        assert_eq!(*s.cat_concordant.get("pe_both_mapped").unwrap(), 1);
        assert_eq!(s.singleton_picard_marked, 1);
        assert_eq!(s.singleton_partner_marked, 0);

        let t = cmp.agg.triage_orphans();
        assert_eq!(t.total_orphans, 1);
        assert_eq!(t.picard_marked, 1);
        assert_eq!(t.partner_marked, 0);
        assert_eq!(t.cat1_cross_table, 1);
        assert_eq!(t.discordant, 1);

        // Inheritance: every group is dup-flag-consistent on both sides; only
        // pairX carries a supplementary record.
        assert_eq!(cmp.picard_inh.total_qnames, 3);
        assert_eq!(cmp.picard_inh.groups_with_supp, 1);
        assert_eq!(cmp.picard_inh.inconsistent, 0);
        assert_eq!(cmp.partner_inh.total_qnames, 3);
        assert_eq!(cmp.partner_inh.groups_with_supp, 1);
        assert_eq!(cmp.partner_inh.inconsistent, 0);
    }

    #[test]
    fn compare_tolerates_differing_sec_supp_between_tools() {
        // Picard has a supplementary record for 'a' that the partner lacks. The
        // primary records match, so group-level lockstep must not choke on the
        // per-group record-count difference (the old per-record lockstep would).
        let dir = tempfile::tempdir().unwrap();
        let pic = dir.path().join("p.bam");
        let par = dir.path().join("o.bam");
        write_bam(
            &pic,
            &[
                rec(b"a", FLAG_PAIRED, Some(b"0:1:F")),
                rec(b"a", FLAG_PAIRED, Some(b"0:2:R")),
                rec(b"a", SUPP, None), // picard-only supplementary
            ],
        );
        write_bam(&par, &[rec(b"a", FLAG_PAIRED, None), rec(b"a", FLAG_PAIRED, None)]);
        let cmp = compare(&mut group_reader(&pic), &mut group_reader(&par)).unwrap();
        assert_eq!(cmp.n_groups, 1);
        assert_eq!(cmp.n_primary, 2);
        assert_eq!(cmp.picard_inh.groups_with_supp, 1);
        assert_eq!(cmp.partner_inh.groups_with_supp, 0);
    }

    #[test]
    fn compare_errors_on_qname_mismatch() {
        let dir = tempfile::tempdir().unwrap();
        let pic = dir.path().join("p.bam");
        let par = dir.path().join("o.bam");
        write_bam(&pic, &[rec(b"a", FLAG_PAIRED, Some(b"0:1:F"))]);
        write_bam(&par, &[rec(b"b", FLAG_PAIRED, None)]);
        let err = match compare(&mut group_reader(&pic), &mut group_reader(&par)) {
            Err(e) => e,
            Ok(_) => panic!("expected an error"),
        };
        assert!(err.to_string().contains("QNAME mismatch"), "got: {err}");
    }

    #[test]
    fn compare_errors_on_group_count_mismatch() {
        let dir = tempfile::tempdir().unwrap();
        let pic = dir.path().join("p.bam");
        let par = dir.path().join("o.bam");
        // Picard has two groups (a, b); partner has only one (a).
        write_bam(
            &pic,
            &[rec(b"a", FLAG_PAIRED, Some(b"0:1:F")), rec(b"b", FLAG_PAIRED, Some(b"0:2:R"))],
        );
        write_bam(&par, &[rec(b"a", FLAG_PAIRED, None)]);
        let err = match compare(&mut group_reader(&pic), &mut group_reader(&par)) {
            Err(e) => e,
            Ok(_) => panic!("expected an error"),
        };
        assert!(err.to_string().contains("EOF mismatch"), "got: {err}");
    }
}
