//! BAM reader built on a custom BGZF block loop + `fgumi-raw-bam::RawRecord`.
//!
//! Two notable departures from the htslib BAM reader:
//!
//! 1. **CRC32 verification is skipped.** For uncompressed-BAM input that
//!    check is ~5-7% of wall time.
//! 2. **Stored BGZF blocks are short-circuited.** Producers like
//!    `samtools view -u` and the uncompressed-BAM output of our own writer
//!    emit deflate blocks with `BTYPE = 00` (stored). For those, the
//!    "compressed" payload is literally the uncompressed bytes wrapped in
//!    a 5-byte deflate framing. We detect `BTYPE = 00` by peeking the
//!    first payload byte and read the payload *directly* into the sliding
//!    decompressed buffer, bypassing libdeflater's stored-block memcpy.
//!    The deflate path is unchanged for compressed input.
//!
//! The BAM header binary framing (magic, l_text, text, n_ref, refs) is
//! parsed in-house; only the SAM text portion is handed to `noodles_sam`
//! for typed manipulation (so callers can add a @PG line and re-serialize).

use std::io::{self, BufRead, ErrorKind};
use std::num::NonZeroUsize;

use anyhow::{Context, Result, anyhow, bail};
use fgumi_raw_bam::{BAM_MAGIC, RawRecord};
use libdeflater::Decompressor;
use noodles_sam::Header;
use noodles_sam::header::ReferenceSequences;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::ReferenceSequence;

/// Fixed size in bytes of the BGZF block header (gzip + BGZF-specific subfields).
const BGZF_HEADER_SIZE: usize = 18;
/// Fixed size in bytes of the BGZF block footer (CRC32 + ISIZE).
const BGZF_FOOTER_SIZE: usize = 8;

/// How many BGZF blocks to process per refill. ~64 × 64 KB = ~4 MB; keeps
/// the working set small while amortizing per-block bookkeeping.
const BLOCKS_PER_REFILL: usize = 64;

/// A BGZF block decompresses to at most 64 KiB (the BGZF spec caps the
/// uncompressed payload). A footer ISIZE above this is corrupt — reject it
/// before it drives a multi-gigabyte `buf.resize`.
const MAX_BGZF_ISIZE: usize = 0x1_0000;

/// Generous upper bound on a single BAM record's `block_size`. The spec
/// doesn't fix one, but real records are far below this; a larger value
/// is corruption and would otherwise grow `buf` unboundedly.
const MAX_BAM_RECORD_SIZE: usize = 256 * 1024 * 1024;

/// Generous upper bound on header `l_text` and per-reference `l_name`. Even
/// assemblies with millions of contigs stay well under this; a near-`u32`
/// value is corruption.
const MAX_BAM_HEADER_FIELD: usize = 1024 * 1024 * 1024;

/// Streaming BAM reader.
pub struct RawBamReader<R: BufRead> {
    /// The underlying buffered source (stdin pipe or file).
    inner: R,
    /// libdeflater decompressor, reused across blocks to avoid repeated allocation.
    decompressor: Decompressor,
    /// Sliding window of decompressed bytes not yet returned to caller.
    buf: Vec<u8>,
    /// Read position within `buf`.
    pos: usize,
    /// Reusable scratch for the *compressed* payload of deflate blocks
    /// (stored blocks bypass this and read directly into `buf`).
    compressed_scratch: Vec<u8>,
    /// Set once we hit BGZF EOF.
    eof: bool,
    /// When true, verify each block's CRC32 against the footer.
    check_crc: bool,
}

impl<R: BufRead> RawBamReader<R> {
    /// Wrap a buffered reader positioned at the start of a BAM (BGZF) stream.
    /// Pass `check_crc = true` to verify each block's CRC32 against its
    /// footer (skipped when `false` — a measurable wall-time saving).
    pub fn new(reader: R, check_crc: bool) -> Self {
        Self {
            inner: reader,
            decompressor: Decompressor::new(),
            buf: Vec::with_capacity(BLOCKS_PER_REFILL * 65_536),
            pos: 0,
            compressed_scratch: Vec::with_capacity(65_536),
            eof: false,
            check_crc,
        }
    }

    /// Read and parse the BAM header. Returns the typed `Header` so the
    /// caller can inspect/modify (e.g. add @PG) and re-serialize for output.
    pub fn read_header(&mut self) -> Result<Header> {
        self.ensure_bytes(4)?;
        if &self.buf[self.pos..self.pos + 4] != BAM_MAGIC {
            bail!(
                "Not a BAM file: expected magic {:?}, got {:?}",
                BAM_MAGIC,
                &self.buf[self.pos..self.pos + 4]
            );
        }
        self.pos += 4;

        let l_text = self.read_u32()? as usize;
        if l_text > MAX_BAM_HEADER_FIELD {
            bail!("BAM header text length {l_text} is implausibly large — input is likely corrupt");
        }
        self.ensure_bytes(l_text)?;
        let text = &self.buf[self.pos..self.pos + l_text];
        let mut header: Header = if text.is_empty() {
            Header::default()
        } else {
            std::str::from_utf8(text)
                .context("BAM header text is not valid UTF-8")?
                .parse()
                .context("parsing SAM header text")?
        };
        self.pos += l_text;

        let n_ref = self.read_u32()? as usize;
        let mut binary_refs = ReferenceSequences::with_capacity(n_ref);
        for _ in 0..n_ref {
            let l_name = self.read_u32()? as usize;
            // Reference names are NUL-terminated per BAM spec, so l_name>=1.
            if l_name == 0 {
                bail!("BAM reference name length is zero (spec requires NUL terminator)");
            }
            if l_name > MAX_BAM_HEADER_FIELD {
                bail!(
                    "BAM reference name length {l_name} is implausibly large — input is likely corrupt"
                );
            }
            self.ensure_bytes(l_name)?;
            // Skip the trailing NUL when copying out the name bytes.
            let name = self.buf[self.pos..self.pos + l_name - 1].to_vec();
            self.pos += l_name;
            let l_ref = self.read_u32()? as usize;
            let length = NonZeroUsize::new(l_ref)
                .ok_or_else(|| anyhow!("reference sequence with zero length"))?;
            binary_refs.insert(name.into(), Map::<ReferenceSequence>::new(length));
        }

        if header.reference_sequences().is_empty() {
            *header.reference_sequences_mut() = binary_refs;
        } else {
            // The BAM spec says the text @SQ list and the binary reference
            // list must match. tid values in alignment records are indices
            // into the binary list, so silently disagreeing would route
            // reads to the wrong contig name. Verify and bail on mismatch.
            let text_refs = header.reference_sequences();
            if text_refs.len() != binary_refs.len() {
                bail!(
                    "BAM header @SQ count ({}) does not match binary ref list ({})",
                    text_refs.len(),
                    binary_refs.len(),
                );
            }
            for ((t_name, t_map), (b_name, b_map)) in text_refs.iter().zip(binary_refs.iter()) {
                if t_name != b_name {
                    bail!(
                        "BAM @SQ name mismatch at index: text={:?} binary={:?}",
                        String::from_utf8_lossy(t_name),
                        String::from_utf8_lossy(b_name),
                    );
                }
                if t_map.length() != b_map.length() {
                    bail!(
                        "BAM @SQ length mismatch for {}: text={} binary={}",
                        String::from_utf8_lossy(t_name),
                        usize::from(t_map.length()),
                        usize::from(b_map.length()),
                    );
                }
            }
        }

        Ok(header)
    }

    /// Read one BAM record into `rec`. Returns `Ok(true)` on success,
    /// `Ok(false)` at EOF.
    pub fn read_record(&mut self, rec: &mut RawRecord) -> io::Result<bool> {
        if !self.has_bytes(4) {
            self.try_refill()?;
            if !self.has_bytes(4) {
                return Ok(false);
            }
        }
        let block_size = u32::from_le_bytes([
            self.buf[self.pos],
            self.buf[self.pos + 1],
            self.buf[self.pos + 2],
            self.buf[self.pos + 3],
        ]) as usize;
        self.pos += 4;

        if block_size > MAX_BAM_RECORD_SIZE {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!(
                    "BAM record block_size {block_size} is implausibly large — input is likely corrupt"
                ),
            ));
        }
        self.ensure_bytes_io(block_size)?;
        // `resize(_, 0)` would write `block_size` zeros that are immediately
        // overwritten by `extend_from_slice` below — that's ~5 GB of pointless
        // zero-stores on a 30 GB BAM run. Clear + extend lets the allocator
        // hand us the slot without zeroing.
        let vec = rec.as_mut_vec();
        vec.clear();
        vec.extend_from_slice(&self.buf[self.pos..self.pos + block_size]);
        self.pos += block_size;
        Ok(true)
    }

    /// Read 4 bytes from the decompressed buffer as a little-endian `u32`.
    fn read_u32(&mut self) -> Result<u32> {
        self.ensure_bytes(4)?;
        let n = u32::from_le_bytes([
            self.buf[self.pos],
            self.buf[self.pos + 1],
            self.buf[self.pos + 2],
            self.buf[self.pos + 3],
        ]);
        self.pos += 4;
        Ok(n)
    }

    /// Returns `true` if `n` or more decompressed bytes are available in `buf`.
    #[inline]
    fn has_bytes(&self, n: usize) -> bool {
        self.buf.len() - self.pos >= n
    }

    /// Ensure at least `n` bytes are available in the decompressed buffer,
    /// refilling from BGZF blocks as needed. Returns `Err` on EOF or IO error.
    fn ensure_bytes(&mut self, n: usize) -> Result<()> {
        while !self.has_bytes(n) {
            if !self.try_refill().context("reading BGZF blocks")? {
                return Err(anyhow!(
                    "Unexpected EOF: needed {n} more bytes (have {}/{})",
                    self.buf.len() - self.pos,
                    n
                ));
            }
        }
        Ok(())
    }

    /// Same as [`Self::ensure_bytes`] but returns `io::Result` for callers in the
    /// `Read`/`BufRead` impl chain that cannot propagate `anyhow::Error`.
    fn ensure_bytes_io(&mut self, n: usize) -> io::Result<()> {
        while !self.has_bytes(n) {
            if !self.try_refill()? {
                return Err(io::Error::new(
                    ErrorKind::UnexpectedEof,
                    format!(
                        "Truncated BAM record: needed {n} more bytes (have {}/{})",
                        self.buf.len() - self.pos,
                        n
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Refill the decompressed buffer with one batch of BGZF blocks.
    /// Returns `Ok(true)` if any bytes were appended, `Ok(false)` at EOF.
    fn try_refill(&mut self) -> io::Result<bool> {
        if self.eof {
            return Ok(false);
        }
        // Drain already-consumed bytes from the front of `buf`. `drain` is
        // O(remaining), so on a large buffer with a small tail the memmove
        // cost dominates. Use `copy_within` instead — same memmove but with
        // less bookkeeping — and only run it when there's a non-trivial
        // tail to preserve.
        if self.pos > 0 {
            let remaining = self.buf.len() - self.pos;
            if remaining > 0 {
                self.buf.copy_within(self.pos.., 0);
            }
            self.buf.truncate(remaining);
            self.pos = 0;
        }

        let starting = self.buf.len();
        for _ in 0..BLOCKS_PER_REFILL {
            if !self.read_one_block()? {
                break;
            }
        }
        Ok(self.buf.len() > starting)
    }

    /// Read and process one BGZF block. Returns `Ok(true)` if a block was
    /// consumed (or an EOF marker block was skipped), `Ok(false)` at end
    /// of stream. Stored blocks (`BTYPE = 00`) bypass libdeflater.
    fn read_one_block(&mut self) -> io::Result<bool> {
        // Read the 18-byte BGZF header. UnexpectedEof on the first byte
        // means we're at clean end of stream.
        let mut header = [0u8; BGZF_HEADER_SIZE];
        match self.inner.read_exact(&mut header) {
            Ok(()) => {}
            Err(e) if e.kind() == ErrorKind::UnexpectedEof => {
                self.eof = true;
                return Ok(false);
            }
            Err(e) => return Err(e),
        }

        if header[0] != 0x1f || header[1] != 0x8b {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid BGZF magic 0x{:02x} 0x{:02x}", header[0], header[1]),
            ));
        }
        if header[12] != b'B' || header[13] != b'C' {
            return Err(io::Error::new(ErrorKind::InvalidData, "missing BGZF BC subfield"));
        }

        let bsize = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        if bsize < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("BGZF block too small: {bsize}"),
            ));
        }
        let payload_len = bsize - BGZF_HEADER_SIZE - BGZF_FOOTER_SIZE;

        // Peek the first payload byte to detect a stored (BTYPE = 00) block.
        // The two BTYPE bits are bits 1-2 of the first deflate byte.
        let first = {
            let peek = self.inner.fill_buf()?;
            if peek.is_empty() {
                return Err(io::Error::new(
                    ErrorKind::UnexpectedEof,
                    "truncated BGZF block payload",
                ));
            }
            peek[0]
        };
        let is_stored = (first & 0b110) == 0;

        if is_stored {
            // Stored block layout:
            //   1 byte  BFINAL | BTYPE | (5 unused bits aligned to next byte boundary)
            //   2 bytes LEN  (little-endian)
            //   2 bytes NLEN (one's complement of LEN, ignored here)
            //   LEN bytes payload (= uncompressed data)
            //   8-byte BGZF footer (CRC32 + ISIZE) — skipped without verification
            let mut framing = [0u8; 5];
            self.inner.read_exact(&mut framing)?;
            let len = u16::from_le_bytes([framing[1], framing[2]]) as usize;
            // Sanity: deflate stored framing should match isize. We could
            // verify against ISIZE in the footer, but that's an extra read
            // and the BGZF spec already constrains payload_len to LEN + 5.
            if len + 5 != payload_len {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("stored block size mismatch: LEN={len} payload={payload_len}"),
                ));
            }
            // Read LEN bytes directly into the sliding output buffer — no
            // intermediate scratch, no libdeflater call.
            let start = self.buf.len();
            read_exact_into_spare(&mut self.inner, &mut self.buf, len)?;
            // Read the BGZF footer. Optionally verify CRC32 (libdeflater uses
            // hardware CRC on aarch64/x86).
            let mut footer = [0u8; BGZF_FOOTER_SIZE];
            self.inner.read_exact(&mut footer)?;
            if self.check_crc {
                verify_block_crc(&self.buf[start..start + len], &footer, "stored")?;
            }
        } else {
            // Compressed: read payload into scratch, then inflate via
            // libdeflater. We still need ISIZE from the footer to size the
            // output region — read footer first into a local, then process.
            self.compressed_scratch.clear();
            read_exact_into_spare(&mut self.inner, &mut self.compressed_scratch, payload_len)?;
            let mut footer = [0u8; BGZF_FOOTER_SIZE];
            self.inner.read_exact(&mut footer)?;
            let isize = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;
            if isize > MAX_BGZF_ISIZE {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "BGZF block ISIZE {isize} exceeds the 64 KiB maximum — input is likely corrupt"
                    ),
                ));
            }
            if isize > 0 {
                let start = self.buf.len();
                // Reserve isize bytes of *initialized* slice for libdeflater
                // to overwrite. libdeflater's deflate_decompress requires a
                // `&mut [u8]`, so we still pay one zero-fill per block here —
                // but the alternative (unsafe spare-cap transmute) is risky
                // given libdeflater's API contract is just "writes output."
                self.buf.resize(start + isize, 0);
                let n = self
                    .decompressor
                    .deflate_decompress(&self.compressed_scratch, &mut self.buf[start..])
                    .map_err(|e| {
                        io::Error::new(ErrorKind::InvalidData, format!("inflate failed: {e:?}"))
                    })?;
                if n != isize {
                    return Err(io::Error::new(
                        ErrorKind::InvalidData,
                        format!("BGZF isize mismatch: header={isize} decompressed={n}"),
                    ));
                }
                if self.check_crc {
                    verify_block_crc(&self.buf[start..start + isize], &footer, "deflate")?;
                }
            }
        }

        Ok(true)
    }
}

/// Read `n` bytes from `src` into the spare capacity of `buf`, extending
/// `buf.len()` by `n`. Avoids the `Vec::resize(_, 0)` zero-fill that we would
/// otherwise pay just to hand `read_exact` an initialized slice.
fn read_exact_into_spare<R: std::io::Read>(
    src: &mut R,
    buf: &mut Vec<u8>,
    n: usize,
) -> io::Result<()> {
    use std::mem::MaybeUninit;
    buf.reserve(n);
    let spare: &mut [MaybeUninit<u8>] = &mut buf.spare_capacity_mut()[..n];
    // SAFETY: `read_exact` writes into the slice but is not allowed to read
    // from it (the std contract). On Ok return, `n` bytes are initialized;
    // `set_len` makes them part of `buf`.
    let target: &mut [u8] =
        unsafe { std::slice::from_raw_parts_mut(spare.as_mut_ptr() as *mut u8, spare.len()) };
    src.read_exact(target)?;
    unsafe { buf.set_len(buf.len() + n) };
    Ok(())
}

/// Verify that the BGZF block's CRC32 footer matches a hardware CRC over
/// `decompressed`. `label` is `"stored"` or `"deflate"` to disambiguate the
/// error message.
fn verify_block_crc(
    decompressed: &[u8],
    footer: &[u8; BGZF_FOOTER_SIZE],
    label: &'static str,
) -> io::Result<()> {
    let expected = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
    let actual = libdeflater::crc32(decompressed);
    if actual != expected {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!(
                "BGZF {label}-block CRC32 mismatch: expected 0x{expected:08x}, got 0x{actual:08x}"
            ),
        ));
    }
    Ok(())
}
