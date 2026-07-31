//! BAM writer built on the `bgzf` crate + noodles header serialization.
//!
//! Mirrors the read path: we own the BGZF framing (so we can pick
//! compression level 0 — "stored" — matching what `samtools view -u`
//! produces) and serialize records as raw bytes via
//! [`fgumi_raw_bam::RawRecord`] without going through htslib.
//!
//! BAM-on-disk layout:
//!
//! ```text
//! magic ("BAM\1")
//! l_text   (u32 LE)
//! text     (l_text bytes — SAM @-prefixed header lines)
//! n_ref    (u32 LE)
//! refs:    n_ref × { l_name (u32 LE), name (NUL-terminated, l_name bytes),
//!                    l_ref (u32 LE) }
//! records: each: block_size (u32 LE), <block_size bytes of record payload>
//! ```
//!
//! All of the above is fed into a [`bgzf::Writer`] which compresses + CRCs
//! it into BGZF blocks. The BGZF EOF marker is emitted automatically on
//! `finish()`. Temporary output may additionally sit inside a zstd frame; see
//! [`RawBamWriter::open_temp`].

use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use bgzf::{CompressionLevel, Writer as BgzfWriter};
use fgumi_raw_bam::RawRecord;
use noodles_sam::Header;
use noodles_sam::io::Writer as SamWriter;
use rawb_io::WriteBehind;
use tempfile::NamedTempFile;

/// BGZF magic bytes "BAM\1" — first 4 bytes of every BAM file.
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// The destinations we write BAM bytes to. No buffering of their own — BGZF
/// hands down whole blocks, and a [`WriteBehind`] sits in front.
enum Sink {
    /// A regular file opened for writing.
    File(std::fs::File),
    /// Standard output (used when `--output` is absent or is `"-"`).
    Stdout(io::Stdout),
    /// A file behind a zstd stream, for temporary output only.
    ///
    /// Nesting the codec inside the sink rather than around the BGZF writer keeps
    /// the writer's type — and so every signature that carries it — unchanged, and
    /// puts the compression on the [`WriteBehind`] thread instead of the worker.
    /// The BGZF layer above still runs at level 0, so this is stored blocks inside
    /// a zstd frame; that framing costs ~26 bytes per 64 KB and compresses away.
    Zstd(Box<zstd::stream::write::Encoder<'static, std::fs::File>>),
}

impl Sink {
    /// Close the sink, writing the zstd frame epilogue when there is one.
    ///
    /// Without the epilogue a decoder rejects the file as truncated, so this must
    /// run before anything reads a compressed temp file back.
    fn finish(self) -> Result<()> {
        match self {
            Sink::File(mut f) => f.flush().context("flushing output file")?,
            Sink::Stdout(mut s) => s.flush().context("flushing stdout")?,
            Sink::Zstd(encoder) => {
                encoder.finish().context("finishing zstd temp output")?;
            }
        }
        Ok(())
    }
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Sink::File(f) => f.write(buf),
            Sink::Stdout(s) => s.write(buf),
            Sink::Zstd(e) => e.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Sink::File(f) => f.flush(),
            Sink::Stdout(s) => s.flush(),
            Sink::Zstd(e) => e.flush(),
        }
    }
}

/// BAM writer that owns its BGZF framing.
///
/// Use [`RawBamWriter::open`] with a [`CompressionLevel`] to construct.
/// Level 0 produces "stored" (uncompressed) BGZF blocks — matches
/// `samtools view -u`. Levels 1-9 are standard zlib levels; 10-12 are
/// libdeflate's extra-strong tiers (much more CPU for marginal size
/// wins). After writing all records, call [`Self::finish`] to flush
/// and emit the BGZF EOF marker.
pub struct RawBamWriter {
    /// The underlying BGZF writer wrapping a [`WriteBehind`]. Wrapped in
    /// `Option` so `finish()` can `take()` ownership to drive `BgzfWriter::finish`.
    bgzf: Option<BgzfWriter<WriteBehind<Sink>>>,
    /// Records handed to [`Self::write_record`], so a caller that reads the file
    /// back can tell a short read from a complete one.
    records: u64,
}

impl RawBamWriter {
    /// Open a BAM writer at the given path (or stdout if `None`/`"-"`),
    /// emit the BAM header (magic + SAM text + reference list), and leave
    /// the writer positioned for record output.
    ///
    /// Output goes through a [`WriteBehind`] with a `ring_bytes` ring
    /// buffer so BGZF compression on the worker thread is decoupled from
    /// the actual `write()` syscall on the underlying file/stdout.
    pub fn open(
        path: Option<&Path>,
        header: &Header,
        ring_bytes: usize,
        level: CompressionLevel,
    ) -> Result<Self> {
        let sink = match path {
            Some(p) if p.to_string_lossy() != "-" => {
                let f = std::fs::File::create(p)
                    .with_context(|| format!("creating {}", p.display()))?;
                Sink::File(f)
            }
            _ => Sink::Stdout(io::stdout()),
        };
        Self::from_sink(sink, header, ring_bytes, level)
    }

    /// Open a BAM writer for a temporary file, optionally behind zstd.
    ///
    /// Taking a [`NamedTempFile`] rather than a path is what keeps the zstd sink
    /// away from the primary output, where it would be actively harmful:
    /// [`Self::abandon`] deliberately leaves a file without its BGZF EOF marker so
    /// a reader can see it is incomplete, and that says nothing to a reader who
    /// cannot parse the outer frame at all. A temp file is never read by anything
    /// but this process, and is unlinked either way.
    ///
    /// BGZF stays at level 0 either way; `zstd_level` adds a second layer around
    /// it rather than replacing it.
    pub fn open_temp(
        temp: &NamedTempFile,
        header: &Header,
        ring_bytes: usize,
        zstd_level: Option<i32>,
    ) -> Result<Self> {
        let path = temp.path();
        let file =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let sink = match zstd_level {
            None => Sink::File(file),
            Some(level) => {
                let encoder = zstd::stream::write::Encoder::new(file, level)
                    .with_context(|| format!("starting zstd for {}", path.display()))?;
                Sink::Zstd(Box::new(encoder))
            }
        };
        let stored = CompressionLevel::new(0).expect("compression level 0 is valid");
        Self::from_sink(sink, header, ring_bytes, stored)
    }

    /// Wrap `sink` in the write-behind and BGZF layers and emit the BAM header.
    fn from_sink(
        sink: Sink,
        header: &Header,
        ring_bytes: usize,
        level: CompressionLevel,
    ) -> Result<Self> {
        let threaded = WriteBehind::with_thread_name(sink, ring_bytes, "dupblaster");
        let mut bgzf = BgzfWriter::new(threaded, level);
        write_bam_header(&mut bgzf, header)?;
        Ok(Self { bgzf: Some(bgzf), records: 0 })
    }

    /// Write one BAM record (block_size prefix + payload).
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        let w = self.bgzf.as_mut().expect("writer already finished");
        let block_size = u32::try_from(rec.len())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "record exceeds u32"))?;
        w.write_all(&block_size.to_le_bytes())?;
        w.write_all(rec.as_ref())?;
        self.records += 1;
        Ok(())
    }

    /// How many records have been written so far.
    ///
    /// A reader cannot get this from the file: a truncated stream ends with the
    /// same `UnexpectedEof` that a complete one does, so anything reading a file
    /// back needs the count from the side that wrote it.
    pub fn records_written(&self) -> u64 {
        self.records
    }

    /// Flush BGZF, emit the EOF marker, drain the IO writer thread, join it, and
    /// close the sink. After this the writer is unusable.
    pub fn finish(mut self) -> Result<()> {
        if let Some(w) = self.bgzf.take() {
            let threaded = w.finish().context("flushing BGZF writer")?;
            let sink = threaded.finish().context("flushing IO writer thread")?;
            sink.finish()?;
        }
        Ok(())
    }

    /// Give up on the output *without* marking it complete, for an aborting run.
    ///
    /// The partial file is deliberately left on disk — a failed run that also
    /// vanished its output would be worse — but it must not look finished. The
    /// BGZF EOF marker is what `samtools quickcheck` and every other reader use
    /// to tell a complete file from a truncated one, so emitting it here would
    /// make a short file pass every integrity check it has.
    ///
    /// `bgzf::Writer`'s `Drop` writes that marker unconditionally, so simply
    /// dropping is not an option: this flushes what is already buffered (keeping
    /// as much of the partial output as possible) and then **leaks** the writer
    /// to suppress `Drop`. Leaking is sound precisely here — the process is
    /// exiting non-zero, so the OS closes the file and reaps the write-behind
    /// thread, and nothing observes the leaked buffer.
    pub fn abandon(mut self) {
        if let Some(mut w) = self.bgzf.take() {
            // Best-effort: the run is already failing, so a flush error changes
            // nothing about the outcome.
            let _ = w.flush();
            std::mem::forget(w);
        }
    }
}

/// Serialize the BAM-on-disk header into `writer`:
/// magic, l_text, SAM text, n_ref, refs.
fn write_bam_header<W: Write>(writer: &mut W, header: &Header) -> Result<()> {
    writer.write_all(BAM_MAGIC).context("writing BAM magic")?;
    // SAM text: noodles serializes the typed Header back to @-prefixed lines.
    let text = serialize_sam_text(header)?;
    let l_text =
        u32::try_from(text.len()).map_err(|_| anyhow::anyhow!("BAM header text exceeds u32"))?;
    writer.write_all(&l_text.to_le_bytes()).context("writing l_text")?;
    writer.write_all(&text).context("writing SAM text")?;

    // Binary reference list: n_ref u32, then per ref { l_name u32, name NUL,
    // l_ref u32 }. noodles' write_reference_sequences is private to the
    // crate, so we inline the loop here (it's ~6 lines).
    let refs = header.reference_sequences();
    let n_ref = u32::try_from(refs.len()).map_err(|_| anyhow::anyhow!("n_ref exceeds u32"))?;
    writer.write_all(&n_ref.to_le_bytes()).context("writing n_ref")?;
    for (name, map) in refs {
        let l_name = u32::try_from(name.len() + 1)
            .map_err(|_| anyhow::anyhow!("ref name length exceeds u32"))?;
        writer.write_all(&l_name.to_le_bytes()).context("writing l_name")?;
        writer.write_all(name).context("writing ref name")?;
        writer.write_all(&[0]).context("writing ref name NUL")?;
        let l_ref = u32::try_from(usize::from(map.length()))
            .map_err(|_| anyhow::anyhow!("ref length exceeds u32"))?;
        writer.write_all(&l_ref.to_le_bytes()).context("writing l_ref")?;
    }
    Ok(())
}

/// Serialize a noodles `Header` back into its SAM text form (the `@`-prefixed
/// lines that appear between `l_text` and `n_ref` in a BAM file).
fn serialize_sam_text(header: &Header) -> Result<Vec<u8>> {
    let mut buf = SamWriter::new(Vec::new());
    buf.write_header(header).context("serializing SAM header text")?;
    Ok(buf.into_inner())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write `count` records to a temp BAM at `level`, returning its size on disk.
    fn temp_bam_bytes(level: Option<i32>, count: usize) -> u64 {
        let temp = NamedTempFile::new().expect("temp file");
        let mut writer = RawBamWriter::open_temp(&temp, &Header::default(), 1 << 20, level)
            .expect("writer opens");
        let mut record = RawRecord::new();
        record.as_mut_vec().resize(64, 0);
        for _ in 0..count {
            writer.write_record(&record).expect("write");
        }
        writer.finish().expect("finish");
        std::fs::metadata(temp.path()).expect("stat").len()
    }

    #[test]
    fn compressing_the_temp_bam_makes_it_smaller() {
        // Without this, compression could silently stop happening on both the write
        // and the read side and every round-trip test would still pass.
        let plain = temp_bam_bytes(None, 20_000);
        let compressed = temp_bam_bytes(Some(1), 20_000);
        assert!(compressed < plain, "zstd grew the temp BAM: {compressed} vs {plain}");
    }

    #[test]
    fn a_truncated_compressed_temp_bam_never_reads_back_complete() {
        // The reason `records_written` exists. Truncation shows up one of two ways
        // depending on where the surviving data ends: mid-block it raises
        // "incomplete frame", but landing exactly on a BGZF block boundary makes the
        // reader's 18-byte header read return `UnexpectedEof`, which it treats as a
        // clean end of stream. Only the second is silent, and only comparing against
        // the write-side count catches it — so what is pinned here is the property
        // the check relies on: a truncated file never yields the full record count.
        use crate::raw_reader::RawBamReader;

        let temp = NamedTempFile::new().expect("temp file");
        let written = 20_000u64;
        let mut writer = RawBamWriter::open_temp(&temp, &Header::default(), 1 << 20, Some(1))
            .expect("writer opens");
        let mut record = RawRecord::new();
        for i in 0..written {
            // Varied content, so the frame spans many blocks and a truncation
            // near the tail still leaves most of them decodable.
            let bytes = record.as_mut_vec();
            bytes.clear();
            bytes.extend((0..64u64).map(|b| (i.wrapping_mul(31).wrapping_add(b)) as u8));
            writer.write_record(&record).expect("write");
        }
        assert_eq!(writer.records_written(), written);
        writer.finish().expect("finish");

        // Drop only the tail: the epilogue and the last block or two.
        let full = std::fs::metadata(temp.path()).expect("stat").len();
        std::fs::OpenOptions::new()
            .write(true)
            .open(temp.path())
            .expect("reopen")
            .set_len(full - full / 10)
            .expect("truncate");

        let file = std::fs::File::open(temp.path()).expect("open");
        let decoder = zstd::stream::read::Decoder::new(file).expect("decoder");
        let mut reader = RawBamReader::new(std::io::BufReader::new(decoder), false);
        let mut recovered = 0u64;
        let mut scratch = RawRecord::new();
        let read_back = reader.read_header().map_err(|_| ()).and_then(|_| {
            loop {
                match reader.read_record(&mut scratch) {
                    Ok(true) => recovered += 1,
                    Ok(false) => return Ok(()),
                    Err(_) => return Err(()),
                }
            }
        });
        assert!(
            read_back.is_err() || recovered < written,
            "a truncated temp BAM reported success with all {written} records"
        );
    }

    #[test]
    fn the_writer_counts_the_records_it_wrote() {
        let temp = NamedTempFile::new().expect("temp file");
        let mut writer = RawBamWriter::open_temp(&temp, &Header::default(), 1 << 20, None)
            .expect("writer opens");
        let mut record = RawRecord::new();
        record.as_mut_vec().resize(64, 0);
        for _ in 0..7 {
            writer.write_record(&record).expect("write");
        }
        assert_eq!(writer.records_written(), 7);
    }
}
