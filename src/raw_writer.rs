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
//! `finish()`.

use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use bgzf::{CompressionLevel, Writer as BgzfWriter};
use fgumi_raw_bam::RawRecord;
use noodles_sam::Header;
use noodles_sam::io::Writer as SamWriter;

use crate::io_threading::ThreadedWriter;

/// BGZF magic bytes "BAM\1" — first 4 bytes of every BAM file.
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// Output sinks we open. Files go through a buffered writer; stdout uses
/// the line-buffered default (BGZF writes are large enough that extra
/// buffering doesn't help much).
enum Sink {
    /// A regular file opened for writing.
    File(std::fs::File),
    /// Standard output (used when `--output` is absent or is `"-"`).
    Stdout(io::Stdout),
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Sink::File(f) => f.write(buf),
            Sink::Stdout(s) => s.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Sink::File(f) => f.flush(),
            Sink::Stdout(s) => s.flush(),
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
    /// The underlying BGZF writer wrapping a [`ThreadedWriter`]. Wrapped in
    /// `Option` so `finish()` can `take()` ownership to drive `BgzfWriter::finish`.
    bgzf: Option<BgzfWriter<ThreadedWriter>>,
}

impl RawBamWriter {
    /// Open a BAM writer at the given path (or stdout if `None`/`"-"`),
    /// emit the BAM header (magic + SAM text + reference list), and leave
    /// the writer positioned for record output.
    ///
    /// Output goes through a [`ThreadedWriter`] with a `ring_bytes` ring
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
        let threaded = ThreadedWriter::new(sink, ring_bytes);
        let mut bgzf = BgzfWriter::new(threaded, level);
        write_bam_header(&mut bgzf, header)?;
        Ok(Self { bgzf: Some(bgzf) })
    }

    /// Write one BAM record (block_size prefix + payload).
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        let w = self.bgzf.as_mut().expect("writer already finished");
        let block_size = u32::try_from(rec.len())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "record exceeds u32"))?;
        w.write_all(&block_size.to_le_bytes())?;
        w.write_all(rec.as_ref())?;
        Ok(())
    }

    /// Flush BGZF, emit the EOF marker, drain the IO writer thread, and
    /// join it. After this the writer is unusable.
    pub fn finish(mut self) -> Result<()> {
        if let Some(w) = self.bgzf.take() {
            let threaded = w.finish().context("flushing BGZF writer")?;
            threaded.finish().context("flushing IO writer thread")?;
        }
        Ok(())
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
