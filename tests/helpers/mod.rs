//! Test helpers — build SAM inputs programmatically and decode SAM/BAM
//! output back to records (via noodles) for assertions.
//!
//! Reading is done entirely in-process with noodles (no `samtools` shell-out):
//! `read_recs_and_header` auto-detects SAM vs BAM and returns owned
//! `RecordBuf`s plus the header. `RecordBuf` normalizes aux data (de-duped,
//! last value wins), which is fine for asserting on tag *values*; the one place
//! a test must observe *physically* duplicate tags uses `count_tag_occurrences`,
//! which reads via noodles-bam's lazy record instead.
#![allow(dead_code)] // helpers used by some tests but not all

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, io::Write as _, record::data::field::Tag},
};
use noodles_util::alignment::io::reader::Builder as AlignmentReaderBuilder;
use tempfile::TempDir;

/// Path to the dupblaster binary built by cargo, located via CARGO_BIN_EXE_*.
pub fn rust_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_dupblaster"))
}

/// Builder for a SAM file in memory; writes to a temp file on demand.
pub struct SamBuilder {
    pub header: String,
    pub records: Vec<String>,
}

impl SamBuilder {
    pub fn new() -> Self {
        Self { header: String::from("@HD\tVN:1.6\tSO:unsorted\n"), records: Vec::new() }
    }

    pub fn sq(mut self, name: &str, length: usize) -> Self {
        use std::fmt::Write as _;
        let _ = writeln!(self.header, "@SQ\tSN:{name}\tLN:{length}");
        self
    }

    /// Append an `@RG` read-group line. `library` (the `LB:` field) is
    /// optional. Keeps tab/field formatting in one place so tests don't
    /// hand-splice header text.
    pub fn rg(mut self, id: &str, sample: &str, library: Option<&str>) -> Self {
        use std::fmt::Write as _;
        let _ = write!(self.header, "@RG\tID:{id}\tSM:{sample}");
        if let Some(lb) = library {
            let _ = write!(self.header, "\tLB:{lb}");
        }
        self.header.push('\n');
        self
    }

    /// Append a SAM record. Builder-style positional args mirror the SAM spec.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn record(
        mut self,
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u32,
        mapq: u8,
        cigar: &str,
        rnext: &str,
        pnext: u32,
        tlen: i32,
        seq: &str,
        qual: &str,
    ) -> Self {
        let seq = if seq.is_empty() { "*" } else { seq };
        let qual = if qual.is_empty() { "*" } else { qual };
        self.records.push(format!(
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        ));
        self
    }

    /// Convenience: a simple aligned record with default sequence and quality.
    #[allow(clippy::too_many_arguments)]
    pub fn rec_simple(
        self,
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u32,
        cigar: &str,
        mate_rname: &str,
        mate_pos: u32,
        tlen: i32,
    ) -> Self {
        let seqlen = cigar_query_len(cigar);
        let seq = "A".repeat(seqlen);
        let qual = "I".repeat(seqlen);
        self.record(qname, flag, rname, pos, 60, cigar, mate_rname, mate_pos, tlen, &seq, &qual)
    }

    /// Like [`Self::rec_simple`] but attaches an `RG:Z:<rg>` read-group tag, so
    /// library-aware tests can route a record to a particular `@RG` (and hence
    /// library). Keeps the aux-field formatting out of the test bodies.
    #[allow(clippy::too_many_arguments)]
    pub fn rec_simple_rg(
        self,
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u32,
        cigar: &str,
        mate_rname: &str,
        mate_pos: u32,
        tlen: i32,
        rg: &str,
    ) -> Self {
        let mut builder =
            self.rec_simple(qname, flag, rname, pos, cigar, mate_rname, mate_pos, tlen);
        if let Some(last) = builder.records.last_mut() {
            last.push_str(&format!("\tRG:Z:{rg}"));
        }
        builder
    }

    pub fn write_to(&self, path: &Path) {
        let mut f = fs::File::create(path).expect("create sam");
        f.write_all(self.header.as_bytes()).unwrap();
        for r in &self.records {
            f.write_all(r.as_bytes()).unwrap();
            f.write_all(b"\n").unwrap();
        }
    }
}

/// Roughly compute the query length implied by a CIGAR string for synthetic
/// SAM building. Recognizes M/=/X/I/S as query-consuming.
pub fn cigar_query_len(cigar: &str) -> usize {
    if cigar == "*" {
        return 50;
    }
    let mut len = 0usize;
    let mut num = 0usize;
    for ch in cigar.chars() {
        if let Some(d) = ch.to_digit(10) {
            num = num * 10 + d as usize;
        } else {
            if matches!(ch, 'M' | '=' | 'X' | 'I' | 'S') {
                len += num;
            }
            num = 0;
        }
    }
    len.max(1)
}

/// Set up a temp dir with conventional input/output paths.
pub struct TestEnv {
    pub _tmp: TempDir,
    pub input: PathBuf,
    pub rust_out: PathBuf,
}

impl TestEnv {
    pub fn new() -> Self {
        let tmp = TempDir::new().expect("temp dir");
        let p = tmp.path();
        Self { input: p.join("in.sam"), rust_out: p.join("rust.sam"), _tmp: tmp }
    }
}

/// Read a SAM- or BAM-formatted file and return its header and all records as
/// owned `RecordBuf`s. Format is auto-detected by noodles (BGZF magic → BAM,
/// else SAM text), so this works on dupblaster's BAM output and on
/// `SamBuilder`-written SAM alike, with no external `samtools`.
pub fn read_recs_and_header(path: &Path) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader =
        AlignmentReaderBuilder::default().build_from_path(path).expect("open alignment file");
    let header = reader.read_header().expect("read header");
    let records = reader
        .records(&header)
        .map(|r| {
            let r = r.expect("read record");
            RecordBuf::try_from_alignment_record(&header, &*r).expect("record -> RecordBuf")
        })
        .collect();
    (header, records)
}

/// Records only (header discarded) — for assertions that ignore the header.
pub fn read_records(path: &Path) -> Vec<RecordBuf> {
    read_recs_and_header(path).1
}

/// The `@PG` header lines of a SAM/BAM file, as SAM text (re-serialized from
/// the parsed header via noodles). noodles injects no `@PG` of its own, so this
/// is the natural replacement for `samtools view -H --no-PG`.
pub fn read_pg_lines(path: &Path) -> Vec<String> {
    let (header, _records) = read_recs_and_header(path);
    let mut buf = Vec::new();
    sam::io::Writer::new(&mut buf).write_header(&header).expect("serialize header");
    String::from_utf8(buf)
        .expect("header is utf-8")
        .lines()
        .filter(|l| l.starts_with("@PG"))
        .map(String::from)
        .collect()
}

/// Per-record count of how many times a two-byte aux `tag` physically appears,
/// read via noodles-bam's *lazy* record (which preserves duplicate fields,
/// unlike `RecordBuf`). Used to assert mate-tag idempotency, where the bug is a
/// second physical copy of an existing tag.
pub fn count_tag_occurrences(bam_path: &Path, tag: [u8; 2]) -> Vec<usize> {
    let want = Tag::new(tag[0], tag[1]);
    let mut reader = bam::io::Reader::new(File::open(bam_path).expect("open bam"));
    reader.read_header().expect("read header");
    reader
        .records()
        .map(|r| {
            let r = r.expect("read record");
            // Fail loudly on a malformed aux field rather than silently
            // dropping it — this helper exists to count tags precisely.
            r.data()
                .iter()
                .filter(|f| {
                    let (tag, _) = f.as_ref().expect("decode aux field");
                    *tag == want
                })
                .count()
        })
        .collect()
}

/// Value of a string-typed aux tag on a record (e.g. `MC`), or None.
pub fn tag_string(rec: &RecordBuf, tag: [u8; 2]) -> Option<String> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::new(tag[0], tag[1]))? {
        Value::String(s) => Some(s.to_string()),
        _ => None,
    }
}

/// Integer value of an integer-typed aux tag on a record (e.g. `MQ`), or None.
pub fn tag_int(rec: &RecordBuf, tag: [u8; 2]) -> Option<i64> {
    rec.data().get(&Tag::new(tag[0], tag[1]))?.as_int()
}

/// Convert a SAM file to BAM in-process (replaces `samtools view -b`), used to
/// build BAM *inputs* for the BAM-handling tests.
pub fn sam_to_bam(sam_path: &Path, bam_path: &Path) {
    let mut reader = AlignmentReaderBuilder::default().build_from_path(sam_path).expect("open sam");
    let header = reader.read_header().expect("read header");
    let mut writer = bam::io::Writer::new(File::create(bam_path).expect("create bam"));
    writer.write_header(&header).expect("write header");
    for r in reader.records(&header) {
        let r = r.expect("read record");
        writer.write_alignment_record(&header, &*r).expect("write record");
    }
    writer.try_finish().expect("finish bam");
}

// ── SAM FLAG bit constants (for test fixtures and assertions) ─────────────

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST_SEGMENT: u16 = 0x40;
pub const FLAG_LAST_SEGMENT: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Run the dupblaster binary on `sam_input` writing BAM to `bam_out`, then
/// decode back and return `(qname, flag)` tuples in output order. Extra
/// arguments to dupblaster are passed verbatim.
pub fn run_and_extract_flags(
    sam_input: &Path,
    bam_out: &Path,
    extra: &[&str],
) -> Vec<(String, u16)> {
    let out = Command::new(rust_binary())
        .args(["-i"])
        .arg(sam_input)
        .args(["-o"])
        .arg(bam_out)
        .args(extra)
        .output()
        .expect("rust dupblaster ran");
    assert!(
        out.status.success(),
        "rust dupblaster failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (_header, records) = read_recs_and_header(bam_out);
    records
        .iter()
        .map(|r| {
            let qname = r.name().map(|n| n.to_string()).unwrap_or_default();
            (qname, u16::from(r.flags()))
        })
        .collect()
}
