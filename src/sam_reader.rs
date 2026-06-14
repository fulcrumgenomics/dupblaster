//! SAM text input → `RawRecord` (BAM byte layout).
//!
//! The producer-side workflow for dupblaster is typically an aligner
//! emitting SAM on stdout. We read that text and convert each line into
//! BAM-on-disk record bytes so the rest of the pipeline (`dedup.rs`,
//! `RawBamWriter`) is identical for SAM and BAM inputs.
//!
//! Efficiency notes:
//! * Newlines and tabs are located with `memchr` (SIMD via NEON/SSE/AVX
//!   under the hood).
//! * Each record is built **in place** into the caller's [`RawRecord`]
//!   Vec — its capacity persists across records, so steady-state we do
//!   zero per-record heap allocations.
//! * The BAM fixed header is written **last** (at offset 0..32) so we
//!   can lay out variable-length fields incrementally without an
//!   intermediate buffer.
//! * SEQ packing uses a precomputed 256-byte ASCII→nibble table.
//! * Integer aux tags are emitted in the smallest fitting BAM subtype
//!   (`c`/`C`/`s`/`S`/`i`/`I`) for byte-level parity with `samtools view`.

use std::collections::HashMap;
use std::io::{self, BufRead, ErrorKind};

use anyhow::{Context, Result, anyhow, bail};
use fgumi_raw_bam::RawRecord;
use memchr::{memchr, memchr_iter};
use noodles_sam::Header;
use wide::u8x16;

/// Bytes read per refill from the underlying buffered reader.
const READ_CHUNK: usize = 256 * 1024;

/// ASCII base → BAM 4-bit nibble. Default = `N` (15).
static ENCODE_BASE: [u8; 256] = build_encode_table();

/// Build the `ENCODE_BASE` lookup table at compile time. All 256 entries default
/// to `N` (15); the IUPAC bases accepted by dupblaster are then set explicitly.
const fn build_encode_table() -> [u8; 256] {
    let mut t = [15u8; 256];
    t[b'=' as usize] = 0;
    t[b'A' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'M' as usize] = 3;
    t[b'G' as usize] = 4;
    t[b'R' as usize] = 5;
    t[b'S' as usize] = 6;
    t[b'V' as usize] = 7;
    t[b'T' as usize] = 8;
    t[b'W' as usize] = 9;
    t[b'Y' as usize] = 10;
    t[b'H' as usize] = 11;
    t[b'K' as usize] = 12;
    t[b'D' as usize] = 13;
    t[b'B' as usize] = 14;
    t[b'N' as usize] = 15;
    // Accept lowercase too — IUPAC allows it.
    t[b'a' as usize] = 1;
    t[b'c' as usize] = 2;
    t[b'g' as usize] = 4;
    t[b't' as usize] = 8;
    t[b'n' as usize] = 15;
    t
}

/// Streaming SAM reader. Parses one record per `read_record()` call,
/// emitting BAM-format bytes into the caller's [`RawRecord`].
pub struct SamReader<R: BufRead> {
    /// Underlying buffered source of SAM text bytes.
    inner: R,
    /// Sliding window of unconsumed SAM text. Drained from the front on
    /// each refill so consumed bytes don't accumulate.
    buf: Vec<u8>,
    /// Read position within `buf`.
    pos: usize,
    /// `RNAME` (and `RNEXT`) → BAM tid. Populated from `@SQ` lines.
    tid_by_name: HashMap<Vec<u8>, i32>,
    /// Set once the underlying reader returns 0 bytes.
    eof: bool,
    /// Scratch for parsed CIGAR ops while building a record — keeps the
    /// allocation alive across `read_record` calls.
    cigar_scratch: Vec<u32>,
}

impl<R: BufRead> SamReader<R> {
    /// Wrap a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            inner: reader,
            buf: Vec::with_capacity(READ_CHUNK),
            pos: 0,
            tid_by_name: HashMap::new(),
            eof: false,
            cigar_scratch: Vec::with_capacity(32),
        }
    }

    /// Read and parse `@`-prefixed SAM header lines. Returns when the
    /// next line does not start with `@` (or EOF), leaving the reader
    /// positioned at the first record line.
    pub fn read_header(&mut self) -> Result<Header> {
        // We assemble the header bytes for noodles to parse.
        let mut header_text: Vec<u8> = Vec::new();
        loop {
            // Make sure we can see at least one byte.
            self.ensure_some_bytes_or_eof().context("reading SAM header")?;
            if self.pos == self.buf.len() {
                break; // EOF before any record line
            }
            if self.buf[self.pos] != b'@' {
                break; // Reached first alignment record
            }
            // Find end of this header line.
            let nl = self.find_newline_or_eof().context("reading SAM header line")?;
            header_text.extend_from_slice(&self.buf[self.pos..nl]);
            header_text.push(b'\n');
            self.pos = (nl + 1).min(self.buf.len());
        }

        let header: Header = if header_text.is_empty() {
            Header::default()
        } else {
            std::str::from_utf8(&header_text)
                .context("SAM header is not valid UTF-8")?
                .parse()
                .context("parsing SAM header")?
        };

        // Build name→tid map from @SQ entries in their order in the header.
        for (idx, (name, _)) in header.reference_sequences().iter().enumerate() {
            self.tid_by_name.insert(name.to_vec(), idx as i32);
        }
        Ok(header)
    }

    /// Read one SAM record, encoding it as BAM bytes into `rec`. Returns
    /// `Ok(true)` on success, `Ok(false)` at EOF.
    pub fn read_record(&mut self, rec: &mut RawRecord) -> io::Result<bool> {
        loop {
            // Drain any blank lines.
            while self.pos < self.buf.len() && self.buf[self.pos] == b'\n' {
                self.pos += 1;
            }
            if self.pos >= self.buf.len() && !self.try_refill()? {
                return Ok(false);
            }

            // Find the next newline.
            if let Some(off) = memchr(b'\n', &self.buf[self.pos..]) {
                let line_start = self.pos;
                let mut line_end = self.pos + off;
                self.pos = line_end + 1;
                // Strip trailing \r if present.
                if line_end > line_start && self.buf[line_end - 1] == b'\r' {
                    line_end -= 1;
                }
                if line_end == line_start {
                    continue; // empty line; skip
                }
                let SamReader { buf, tid_by_name, cigar_scratch, .. } = self;
                parse_sam_line(&buf[line_start..line_end], rec, tid_by_name, cigar_scratch)
                    .map_err(|e| {
                        io::Error::new(ErrorKind::InvalidData, format!("parsing SAM record: {e:#}"))
                    })?;
                return Ok(true);
            }

            // No newline in current buffer — refill or treat trailing
            // bytes as a final (unterminated) record.
            if !self.try_refill()? {
                if self.buf.len() > self.pos {
                    let line_start = self.pos;
                    let line_end = self.buf.len();
                    self.pos = line_end;
                    let SamReader { buf, tid_by_name, cigar_scratch, .. } = self;
                    parse_sam_line(&buf[line_start..line_end], rec, tid_by_name, cigar_scratch)
                        .map_err(|e| {
                            io::Error::new(
                                ErrorKind::InvalidData,
                                format!("parsing SAM record: {e:#}"),
                            )
                        })?;
                    return Ok(true);
                }
                return Ok(false);
            }
        }
    }

    /// Make sure `buf[pos..]` contains at least one byte, or that we've
    /// reached EOF.
    fn ensure_some_bytes_or_eof(&mut self) -> io::Result<()> {
        if self.pos >= self.buf.len() {
            let _ = self.try_refill()?;
        }
        Ok(())
    }

    /// Find the next newline. If we reach EOF without finding one, return
    /// the buffer end position.
    fn find_newline_or_eof(&mut self) -> io::Result<usize> {
        loop {
            if let Some(off) = memchr(b'\n', &self.buf[self.pos..]) {
                return Ok(self.pos + off);
            }
            if !self.try_refill()? {
                return Ok(self.buf.len());
            }
        }
    }

    /// Drain consumed bytes and read up to `READ_CHUNK` more. Returns
    /// `Ok(true)` if any bytes were added, `Ok(false)` at EOF.
    fn try_refill(&mut self) -> io::Result<bool> {
        if self.eof {
            return Ok(false);
        }
        // Compact: copy the unconsumed tail to the front of buf, then
        // truncate. `copy_within` is one memmove; `Vec::drain` does the
        // same memmove plus extra bookkeeping.
        if self.pos > 0 {
            let remaining = self.buf.len() - self.pos;
            if remaining > 0 {
                self.buf.copy_within(self.pos.., 0);
            }
            self.buf.truncate(remaining);
            self.pos = 0;
        }
        // Read directly into the spare capacity to avoid zero-filling
        // bytes we're about to overwrite — see `read_into_spare` for the
        // soundness argument.
        let n = read_into_spare(&mut self.inner, &mut self.buf, READ_CHUNK)?;
        if n == 0 {
            self.eof = true;
            return Ok(false);
        }
        Ok(true)
    }
}

/// Read up to `cap` bytes from `src` into the spare capacity of `buf`,
/// extending `buf`'s length by however many bytes were actually read. This
/// avoids the zero-fill that `Vec::resize(..., 0)` would impose just to
/// hand `Read` an initialized slice — `Read::read` is forbidden from reading
/// from the buffer it's handed, so the uninit bytes never escape.
fn read_into_spare<R: std::io::Read>(
    src: &mut R,
    buf: &mut Vec<u8>,
    cap: usize,
) -> io::Result<usize> {
    use std::mem::MaybeUninit;
    buf.reserve(cap);
    let spare: &mut [MaybeUninit<u8>] = &mut buf.spare_capacity_mut()[..cap];
    // SAFETY: `MaybeUninit<u8>` and `u8` have the same layout. `Read::read`
    // promises not to read from the destination, only write. After `read`
    // returns Ok(n), bytes [0..n] of the spare slice are initialized.
    let target: &mut [u8] =
        unsafe { std::slice::from_raw_parts_mut(spare.as_mut_ptr() as *mut u8, spare.len()) };
    let n = src.read(target)?;
    // SAFETY: `read` wrote `n` initialized bytes; extend `buf.len()` to cover them.
    unsafe { buf.set_len(buf.len() + n) };
    Ok(n)
}

/// Parse one SAM line (without trailing newline) and encode as BAM bytes
/// into `rec`. Layout is built incrementally:
///
///   [32-byte placeholder] qname+NUL | cigar u32s | packed seq | qual | aux
///
/// then the 32-byte placeholder is overwritten with the actual fixed header.
fn parse_sam_line(
    line: &[u8],
    rec: &mut RawRecord,
    tid_by_name: &HashMap<Vec<u8>, i32>,
    cigar_scratch: &mut Vec<u32>,
) -> Result<()> {
    // Tokenize the mandatory fields by tabs. SAM has 11 mandatory fields,
    // separated by 10 tabs; the 11th tab (if present) is the start of aux.
    let mut tabs: [usize; 11] = [0; 11];
    let mut n_tabs = 0usize;
    for off in memchr_iter(b'\t', line) {
        if n_tabs < 11 {
            tabs[n_tabs] = off;
            n_tabs += 1;
        } else {
            break;
        }
    }
    if n_tabs < 10 {
        bail!("SAM line has only {} fields (needs at least 11)", n_tabs + 1);
    }

    let field = |i: usize| -> &[u8] {
        let start = if i == 0 { 0 } else { tabs[i - 1] + 1 };
        let end = if i < n_tabs { tabs[i] } else { line.len() };
        &line[start..end]
    };
    // Aux portion is everything after the 11th field's tab (if any).
    let aux_start = if n_tabs >= 11 { tabs[10] + 1 } else { line.len() };

    let qname = field(0);
    let flag = u16::try_from(parse_u32_ascii(field(1))?)
        .map_err(|_| anyhow!("FLAG field exceeds u16 range"))?;
    let rname = field(2);
    let pos_1based = parse_i32_ascii(field(3))?;
    let mapq = u8::try_from(parse_u32_ascii(field(4))?)
        .map_err(|_| anyhow!("MAPQ field exceeds u8 range"))?;
    let cigar_text = field(5);
    let rnext = field(6);
    let pnext_1based = parse_i32_ascii(field(7))?;
    let tlen = parse_i32_ascii(field(8))?;
    let seq = field(9);
    let qual = field(10);

    let tid = name_to_tid(rname, tid_by_name)?;
    let mtid = if rnext == b"=" { tid } else { name_to_tid(rnext, tid_by_name)? };
    let pos_0based = if pos_1based > 0 { pos_1based - 1 } else { -1 };
    let mpos_0based = if pnext_1based > 0 { pnext_1based - 1 } else { -1 };

    cigar_scratch.clear();
    if cigar_text != b"*" && !cigar_text.is_empty() {
        parse_cigar_ops(cigar_text, cigar_scratch)?;
    }
    let n_cigar = cigar_scratch.len();
    if n_cigar > u16::MAX as usize {
        bail!("CIGAR op count {} exceeds u16::MAX", n_cigar);
    }

    let l_seq = if seq == b"*" { 0usize } else { seq.len() };
    let qual_present = qual != b"*";
    if qual_present && l_seq > 0 && qual.len() != l_seq {
        bail!("QUAL length {} != SEQ length {}", qual.len(), l_seq);
    }

    let bin = if tid < 0 {
        4680u16
    } else {
        let beg = pos_0based.max(0) as i64;
        let ra_len = ref_length_from_cigar(cigar_scratch);
        let end = beg + ra_len.max(1) as i64;
        reg2bin(beg, end) as u16
    };

    let l_qname = qname.len() + 1;
    if l_qname > 255 {
        bail!("qname length {} exceeds u8::MAX", l_qname);
    }

    // Build the variable-length section, leaving 32 zero bytes at the
    // start as a placeholder for the fixed header.
    let vec = rec.as_mut_vec();
    vec.clear();
    vec.resize(32, 0);

    // QNAME + NUL terminator.
    vec.extend_from_slice(qname);
    vec.push(0);

    // CIGAR ops (u32 LE each).
    for &op in cigar_scratch.iter() {
        vec.extend_from_slice(&op.to_le_bytes());
    }

    // Packed sequence.
    if l_seq > 0 {
        pack_seq_into(seq, vec);
    }

    // Quality scores.
    if l_seq > 0 {
        if qual_present {
            decode_qual_into(qual, vec);
        } else {
            // QUAL = "*" → all 0xFF per BAM spec.
            vec.extend(std::iter::repeat_n(0xFFu8, l_seq));
        }
    }

    // Aux tags.
    if aux_start < line.len() {
        parse_aux_tags(&line[aux_start..], vec)?;
    }

    // Now write the 32-byte fixed header in place.
    write_bam_header(
        &mut vec[..32],
        tid,
        pos_0based,
        l_qname as u8,
        mapq,
        bin,
        n_cigar as u16,
        flag,
        l_seq as u32,
        mtid,
        mpos_0based,
        tlen,
    );

    Ok(())
}

/// Write all 32 bytes of the BAM fixed-length record header into the pre-
/// allocated placeholder slice `dst`. Serializes each field as little-endian
/// per the BAM specification.
#[allow(clippy::too_many_arguments)]
fn write_bam_header(
    dst: &mut [u8],
    tid: i32,
    pos: i32,
    l_qname: u8,
    mapq: u8,
    bin: u16,
    n_cigar: u16,
    flag: u16,
    l_seq: u32,
    mtid: i32,
    mpos: i32,
    tlen: i32,
) {
    dst[0..4].copy_from_slice(&tid.to_le_bytes());
    dst[4..8].copy_from_slice(&pos.to_le_bytes());
    dst[8] = l_qname;
    dst[9] = mapq;
    dst[10..12].copy_from_slice(&bin.to_le_bytes());
    dst[12..14].copy_from_slice(&n_cigar.to_le_bytes());
    dst[14..16].copy_from_slice(&flag.to_le_bytes());
    dst[16..20].copy_from_slice(&l_seq.to_le_bytes());
    dst[20..24].copy_from_slice(&mtid.to_le_bytes());
    dst[24..28].copy_from_slice(&mpos.to_le_bytes());
    dst[28..32].copy_from_slice(&tlen.to_le_bytes());
}

/// Look up a RNAME/RNEXT byte slice in the header's tid map. `"*"` returns
/// `-1` (the BAM unmapped sentinel); any other name must be present in the map
/// or an error is returned.
fn name_to_tid(name: &[u8], tid_by_name: &HashMap<Vec<u8>, i32>) -> Result<i32> {
    if name == b"*" {
        return Ok(-1);
    }
    tid_by_name
        .get(name)
        .copied()
        .ok_or_else(|| anyhow!("RNAME/RNEXT '{}' not in header", String::from_utf8_lossy(name)))
}

/// Parse a non-negative ASCII decimal into `u32`. Returns an error on
/// overflow rather than silently wrapping.
fn parse_u32_ascii(bytes: &[u8]) -> Result<u32> {
    if bytes.is_empty() {
        bail!("empty integer field");
    }
    let mut n = 0u32;
    for &b in bytes {
        if !b.is_ascii_digit() {
            bail!("invalid digit in integer: 0x{:02x}", b);
        }
        n = n
            .checked_mul(10)
            .and_then(|x| x.checked_add((b - b'0') as u32))
            .ok_or_else(|| anyhow!("integer overflow parsing u32 from {:?}", bstr_lossy(bytes)))?;
    }
    Ok(n)
}

/// Parse an optionally-signed ASCII decimal into `i32`. Returns an error
/// on overflow rather than silently wrapping.
fn parse_i32_ascii(bytes: &[u8]) -> Result<i32> {
    if bytes.is_empty() {
        bail!("empty integer field");
    }
    let (negative, digits) = match bytes[0] {
        b'-' => (true, &bytes[1..]),
        b'+' => (false, &bytes[1..]),
        _ => (false, bytes),
    };
    if digits.is_empty() {
        bail!("integer field has only a sign");
    }
    let mut n: i32 = 0;
    for &b in digits {
        if !b.is_ascii_digit() {
            bail!("invalid digit in integer: 0x{:02x}", b);
        }
        let digit = (b - b'0') as i32;
        // Build magnitude as a negative number so we can represent i32::MIN.
        n = n
            .checked_mul(10)
            .and_then(|x| x.checked_sub(digit))
            .ok_or_else(|| anyhow!("integer overflow parsing i32 from {:?}", bstr_lossy(bytes)))?;
    }
    if negative {
        Ok(n)
    } else {
        n.checked_neg()
            .ok_or_else(|| anyhow!("integer overflow parsing i32 from {:?}", bstr_lossy(bytes)))
    }
}

/// Render a byte slice for error messages without allocating on the
/// happy path; falls back to the lossy UTF-8 string for invalid bytes.
fn bstr_lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Parse a CIGAR text string (e.g. `"5S45M5S"`) into packed BAM CIGAR `u32`
/// words appended to `out`. Each word encodes `(length << 4) | op_code`
/// per the BAM spec. Returns an error on zero-length ops, overflow, or
/// unknown op characters.
fn parse_cigar_ops(text: &[u8], out: &mut Vec<u32>) -> Result<()> {
    let mut len: u32 = 0;
    let mut any_digits = false;
    for &b in text {
        if b.is_ascii_digit() {
            len = len
                .checked_mul(10)
                .and_then(|x| x.checked_add((b - b'0') as u32))
                .ok_or_else(|| anyhow!("CIGAR op length overflows u32"))?;
            any_digits = true;
        } else {
            if !any_digits {
                bail!("CIGAR op without preceding length");
            }
            if len == 0 {
                bail!("CIGAR op with zero length is invalid per SAM spec");
            }
            // The packed BAM op format reserves the low 4 bits for the
            // op code; lengths must fit in the upper 28 bits.
            if len > (u32::MAX >> 4) {
                bail!("CIGAR op length {len} exceeds 28-bit BAM limit");
            }
            let op_code = match b {
                b'M' => 0u32,
                b'I' => 1,
                b'D' => 2,
                b'N' => 3,
                b'S' => 4,
                b'H' => 5,
                b'P' => 6,
                b'=' => 7,
                b'X' => 8,
                _ => bail!("invalid CIGAR op code: 0x{:02x}", b),
            };
            out.push((len << 4) | op_code);
            len = 0;
            any_digits = false;
        }
    }
    if any_digits {
        bail!("CIGAR has trailing length without op");
    }
    Ok(())
}

/// Pack an ASCII SEQ field into BAM 4-bit nibble pairs, appending to `out`.
/// Processes 16-base chunks via SIMD ([`simd_pack_chunk`]) with a scalar
/// fallback for non-ACGTN= bytes and the remaining tail.
fn pack_seq_into(seq: &[u8], out: &mut Vec<u8>) {
    let n = seq.len();
    out.reserve(n.div_ceil(2));

    // SIMD fast path on 16-base chunks. Falls back to scalar per chunk
    // when the chunk contains anything outside {A,C,G,T,N,=} (case-
    // insensitive) — in real Illumina data that's >99.9% of chunks.
    let mut i = 0;
    while i + 16 <= n {
        // Load 16 ASCII bases.
        let bases_arr: [u8; 16] = seq[i..i + 16].try_into().unwrap();
        let bases = u8x16::new(bases_arr);

        if let Some(packed8) = simd_pack_chunk(bases) {
            // Append 8 packed bytes (16 bases → 8 nibble pairs).
            out.extend_from_slice(&packed8);
        } else {
            // Scalar fallback for this 16-base chunk.
            for j in (0..16).step_by(2) {
                let hi = ENCODE_BASE[seq[i + j] as usize];
                let lo = ENCODE_BASE[seq[i + j + 1] as usize];
                out.push((hi << 4) | lo);
            }
        }
        i += 16;
    }

    // Scalar tail for remaining bases (< 16). Pairs, then a single
    // half-byte if `n` is odd.
    while i + 2 <= n {
        let hi = ENCODE_BASE[seq[i] as usize];
        let lo = ENCODE_BASE[seq[i + 1] as usize];
        out.push((hi << 4) | lo);
        i += 2;
    }
    if i < n {
        let hi = ENCODE_BASE[seq[i] as usize];
        out.push(hi << 4);
    }
}

/// Try to pack 16 ASCII bases into 8 BAM-nibble bytes via SIMD. Returns
/// `None` if the chunk contains anything outside {A,C,G,T,N,=} case-
/// insensitive — caller must fall back to scalar.
///
/// Encoding trick: the low 4 bits of ASCII characters in our allowed set
/// are all distinct, so a 16-byte PSHUFB-style lookup (swizzle_relaxed)
/// converts low-4-bit-of-ASCII → BAM nibble. Then nibble pairs are packed
/// into bytes via even/odd swizzles + a "shift left by 4" lookup table
/// (wide's u8x16 has no shift or multiply).
#[inline]
fn simd_pack_chunk(bases: u8x16) -> Option<[u8; 8]> {
    // Force uppercase by clearing bit 5 (works for A-Z, leaves `=` alone).
    let upper = bases & u8x16::splat(0xDF);
    let is_a = upper.simd_eq(u8x16::splat(b'A'));
    let is_c = upper.simd_eq(u8x16::splat(b'C'));
    let is_g = upper.simd_eq(u8x16::splat(b'G'));
    let is_t = upper.simd_eq(u8x16::splat(b'T'));
    let is_n = upper.simd_eq(u8x16::splat(b'N'));
    let is_eq = bases.simd_eq(u8x16::splat(b'='));
    let valid = is_a | is_c | is_g | is_t | is_n | is_eq;
    // to_bitmask packs the MSB of each lane. All 16 valid → 0xFFFF.
    if valid.to_bitmask() & 0xFFFF != 0xFFFF {
        return None;
    }

    // BAM nibble lookup indexed by ASCII low-4-bits. ACGTN= map to unique
    // low-nibble values; other slots default to N (15) but those lanes
    // are filtered out by the valid-check above.
    //
    // ASCII lows: A=1 C=3 T=4 G=7 N=14 ==13. BAM: A=1 C=2 T=8 G=4 N=15 ==0.
    const BAM_LUT: u8x16 = u8x16::new([
        15, // 0  (unused, defaults to N)
        1,  // 1  → A
        15, // 2  (unused)
        2,  // 3  → C
        8,  // 4  → T
        15, // 5
        15, // 6
        4,  // 7  → G
        15, // 8
        15, // 9
        15, // 10
        15, // 11
        15, // 12
        0,  // 13 → =
        15, // 14 → N (intentional)
        15, // 15
    ]);
    let low4 = bases & u8x16::splat(0x0F);
    let nibbles = BAM_LUT.swizzle_relaxed(low4);

    // Pack adjacent nibble pairs into bytes:
    //   out[i] = (nibbles[2i] << 4) | nibbles[2i+1]
    // wide::u8x16 lacks shift/mul, but we can substitute `<< 4` with a
    // 16-byte lookup table (since nibbles are 0..15).
    const EVENS_IDX: u8x16 = u8x16::new([0, 2, 4, 6, 8, 10, 12, 14, 0, 2, 4, 6, 8, 10, 12, 14]);
    const ODDS_IDX: u8x16 = u8x16::new([1, 3, 5, 7, 9, 11, 13, 15, 1, 3, 5, 7, 9, 11, 13, 15]);
    const SHIFT4_LUT: u8x16 =
        u8x16::new([0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240]);

    let evens = nibbles.swizzle_relaxed(EVENS_IDX);
    let odds = nibbles.swizzle_relaxed(ODDS_IDX);
    let evens_shifted = SHIFT4_LUT.swizzle_relaxed(evens);
    let packed = evens_shifted | odds;

    // First 8 lanes hold our packed output (lanes 8..15 are a duplicate
    // of the same data because of the swizzle index pattern).
    let arr = packed.to_array();
    let mut result = [0u8; 8];
    result.copy_from_slice(&arr[..8]);
    Some(result)
}

/// Decode SAM ASCII quality scores (Phred+33) to BAM raw Phred bytes
/// (subtract 33), 16 bytes at a time via SIMD. Tail handled scalar.
#[inline]
fn decode_qual_into(qual: &[u8], out: &mut Vec<u8>) {
    let n = qual.len();
    let start = out.len();
    out.resize(start + n, 0);
    let dst = &mut out[start..];

    let bias = u8x16::splat(33);
    let mut i = 0;
    while i + 16 <= n {
        let arr: [u8; 16] = qual[i..i + 16].try_into().unwrap();
        let bytes = u8x16::new(arr);
        let result = bytes - bias; // wrapping per wide semantics
        dst[i..i + 16].copy_from_slice(&result.to_array());
        i += 16;
    }
    while i < n {
        dst[i] = qual[i].wrapping_sub(33);
        i += 1;
    }
}

/// Sum the reference-consuming lengths from a slice of packed BAM CIGAR words.
/// Used to compute the `reg2bin` half-open end coordinate for `[beg, beg+ra_len)`.
fn ref_length_from_cigar(ops: &[u32]) -> i32 {
    let mut len = 0i32;
    for &op in ops {
        let l = (op >> 4) as i32;
        // Reference-consuming ops: M (0), D (2), N (3), = (7), X (8).
        match op & 0xF {
            0 | 2 | 3 | 7 | 8 => len = len.saturating_add(l),
            _ => {}
        }
    }
    len
}

/// BAM `reg2bin` — same indexing scheme htslib uses. Pulled verbatim
/// from the SAMv1 spec.
/// Compute the BAM bin number for a 0-based half-open interval `[beg, end)`,
/// per SAMv1 spec section 5.1.1 (the R-tree binning scheme). Returns 0 when
/// the interval spans the largest bin (effectively "no bin").
///
/// Bin offsets are `1 + 8 + 64 + 512 + 4096 = (8^k - 1) / 7` for the five
/// non-root levels; we just inline the constants since they're standardized.
fn reg2bin(beg: i64, end: i64) -> u32 {
    // BAM-spec bin-tree level offsets — see SAMv1.pdf §5.1.1.
    const LEVEL_OFFSETS: [u32; 5] = [4681, 585, 73, 9, 1];
    const LEVEL_SHIFTS: [u32; 5] = [14, 17, 20, 23, 26];
    let end = end - 1;
    for (offset, shift) in LEVEL_OFFSETS.iter().zip(LEVEL_SHIFTS.iter()) {
        if beg >> shift == end >> shift {
            return offset.wrapping_add((beg >> shift) as u32);
        }
    }
    0
}

/// Parse the tab-separated aux-tag portion of a SAM line and append each tag's
/// BAM binary encoding to `out`. Supported types: `A`, `i`, `f`, `Z`, `H`, `B`.
fn parse_aux_tags(text: &[u8], out: &mut Vec<u8>) -> Result<()> {
    for tag in text.split(|&b| b == b'\t') {
        if tag.is_empty() {
            continue;
        }
        if tag.len() < 5 || tag[2] != b':' || tag[4] != b':' {
            bail!("malformed aux tag: '{}'", String::from_utf8_lossy(tag));
        }
        let key = &tag[0..2];
        let ty = tag[3];
        let val = &tag[5..];
        out.extend_from_slice(key);
        match ty {
            b'A' => {
                if val.len() != 1 {
                    bail!("A-type aux value must be exactly one character");
                }
                out.push(b'A');
                out.push(val[0]);
            }
            b'i' => {
                let v = parse_i32_ascii(val)?;
                emit_int_tag_smallest(v, out);
            }
            b'f' => {
                let f: f32 = std::str::from_utf8(val)
                    .context("aux f value not UTF-8")?
                    .parse()
                    .context("parsing f aux value")?;
                out.push(b'f');
                out.extend_from_slice(&f.to_le_bytes());
            }
            b'Z' => {
                out.push(b'Z');
                out.extend_from_slice(val);
                out.push(0);
            }
            b'H' => {
                out.push(b'H');
                out.extend_from_slice(val);
                out.push(0);
            }
            b'B' => emit_b_array_tag(val, out)?,
            _ => bail!("unknown aux type code: 0x{:02x}", ty),
        }
    }
    Ok(())
}

/// Emit an integer aux value using the smallest fitting BAM subtype, so
/// our output is byte-identical to what `samtools view -b` would produce.
fn emit_int_tag_smallest(v: i32, out: &mut Vec<u8>) {
    if v >= 0 {
        if v <= u8::MAX as i32 {
            out.push(b'C');
            out.push(v as u8);
        } else if v <= u16::MAX as i32 {
            out.push(b'S');
            out.extend_from_slice(&(v as u16).to_le_bytes());
        } else {
            out.push(b'I');
            out.extend_from_slice(&(v as u32).to_le_bytes());
        }
    } else if v >= i8::MIN as i32 {
        out.push(b'c');
        out.push((v as i8) as u8);
    } else if v >= i16::MIN as i32 {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else {
        out.push(b'i');
        out.extend_from_slice(&v.to_le_bytes());
    }
}

/// Tokenize the 11 mandatory tab-separated fields, returning the tab
/// positions (length-1 since there are 10 inter-field tabs plus an
/// optional 11th tab before aux). Used by `profile_tests::profile_sam_parse`.
#[cfg(test)]
fn find_tabs_for_profile(line: &[u8], tabs: &mut [usize; 11]) -> usize {
    let mut n = 0usize;
    for off in memchr_iter(b'\t', line) {
        if n < 11 {
            tabs[n] = off;
            n += 1;
        } else {
            break;
        }
    }
    n
}

/// Encode and append a SAM `B`-type aux array (`B:<subtype>,v0,v1,...`) as
/// BAM binary. Two-pass: count elements, then emit the count prefix followed
/// by each element in its fixed-width little-endian representation.
fn emit_b_array_tag(val: &[u8], out: &mut Vec<u8>) -> Result<()> {
    // B array format: "<subtype>,<v0>,<v1>,..."
    if val.is_empty() {
        bail!("empty B-type aux value");
    }
    let subtype = val[0];
    // Validate up-front — we'd otherwise silently pass through an
    // unknown subtype byte and write a corrupt BAM.
    if !matches!(subtype, b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f') {
        bail!("invalid B-array subtype: 0x{:02x}", subtype);
    }
    let rest = if val.len() > 1 && val[1] == b',' { &val[2..] } else { &val[1..] };
    out.push(b'B');
    out.push(subtype);

    // Two passes: count, then emit. We need the count as a u32 prefix.
    let count = if rest.is_empty() {
        0u32
    } else {
        u32::try_from(rest.iter().filter(|&&b| b == b',').count() + 1)
            .map_err(|_| anyhow!("B array element count exceeds u32"))?
    };
    out.extend_from_slice(&count.to_le_bytes());

    if count == 0 {
        return Ok(());
    }

    for elem in rest.split(|&b| b == b',') {
        match subtype {
            b'c' => {
                let v = i8::try_from(parse_i32_ascii(elem)?)
                    .map_err(|_| anyhow!("B:c element {:?} out of i8 range", bstr_lossy(elem)))?;
                out.push(v as u8);
            }
            b'C' => {
                let v = u8::try_from(parse_u32_ascii(elem)?)
                    .map_err(|_| anyhow!("B:C element {:?} out of u8 range", bstr_lossy(elem)))?;
                out.push(v);
            }
            b's' => {
                let v = i16::try_from(parse_i32_ascii(elem)?)
                    .map_err(|_| anyhow!("B:s element {:?} out of i16 range", bstr_lossy(elem)))?;
                out.extend_from_slice(&v.to_le_bytes());
            }
            b'S' => {
                let v = u16::try_from(parse_u32_ascii(elem)?)
                    .map_err(|_| anyhow!("B:S element {:?} out of u16 range", bstr_lossy(elem)))?;
                out.extend_from_slice(&v.to_le_bytes());
            }
            b'i' => out.extend_from_slice(&parse_i32_ascii(elem)?.to_le_bytes()),
            b'I' => out.extend_from_slice(&parse_u32_ascii(elem)?.to_le_bytes()),
            b'f' => {
                let f: f32 = std::str::from_utf8(elem)
                    .context("B-f element not UTF-8")?
                    .parse()
                    .context("parsing B-f element")?;
                out.extend_from_slice(&f.to_le_bytes());
            }
            _ => bail!("unknown B subtype: 0x{:02x}", subtype),
        }
    }
    Ok(())
}

#[cfg(test)]
mod parser_tests {
    use super::*;

    #[test]
    fn parse_u32_accepts_zero() {
        assert_eq!(parse_u32_ascii(b"0").unwrap(), 0);
    }

    #[test]
    fn parse_u32_rejects_overflow() {
        // 2^32 = 4_294_967_296, which doesn't fit in u32.
        assert!(parse_u32_ascii(b"4294967296").is_err());
        assert!(parse_u32_ascii(b"9999999999").is_err());
    }

    #[test]
    fn parse_u32_rejects_non_digits() {
        assert!(parse_u32_ascii(b"12a3").is_err());
        assert!(parse_u32_ascii(b"").is_err());
    }

    #[test]
    fn parse_i32_handles_min_value() {
        assert_eq!(parse_i32_ascii(b"-2147483648").unwrap(), i32::MIN);
    }

    #[test]
    fn parse_i32_rejects_overflow_either_direction() {
        assert!(parse_i32_ascii(b"2147483648").is_err());
        assert!(parse_i32_ascii(b"-2147483649").is_err());
    }

    #[test]
    fn parse_cigar_rejects_zero_length_op() {
        let mut out = Vec::new();
        assert!(parse_cigar_ops(b"50M0I50M", &mut out).is_err());
    }

    #[test]
    fn parse_cigar_rejects_28bit_overflow() {
        let mut out = Vec::new();
        // 2^28 = 268435456 doesn't fit in the upper 28 bits.
        assert!(parse_cigar_ops(b"268435456M", &mut out).is_err());
    }

    #[test]
    fn parse_cigar_accepts_max_28bit_length() {
        let mut out = Vec::new();
        let max = u32::MAX >> 4;
        let s = format!("{max}M");
        assert!(parse_cigar_ops(s.as_bytes(), &mut out).is_ok());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0] >> 4, max);
        assert_eq!(out[0] & 0xf, 0); // M
    }

    #[test]
    fn b_array_tag_rejects_unknown_subtype() {
        let mut out = Vec::new();
        assert!(emit_b_array_tag(b"x,1,2,3", &mut out).is_err());
    }

    #[test]
    fn b_array_tag_rejects_out_of_range_element() {
        let mut out = Vec::new();
        // 999 doesn't fit in i8.
        assert!(emit_b_array_tag(b"c,1,2,999", &mut out).is_err());
        out.clear();
        // -1 doesn't fit in u8.
        assert!(emit_b_array_tag(b"C,1,2,-1", &mut out).is_err());
    }

    #[test]
    fn b_array_tag_accepts_valid_inputs() {
        let mut out = Vec::new();
        emit_b_array_tag(b"i,-1,2,300", &mut out).unwrap();
        // 'B','i', count=3 (u32 LE), three i32 LE values.
        assert_eq!(out[0..2], *b"Bi");
        assert_eq!(u32::from_le_bytes(out[2..6].try_into().unwrap()), 3);
    }
}

#[cfg(test)]
mod profile_tests {
    use super::*;
    use std::time::Instant;

    /// Microbench each parse-section against real records. Path and
    /// record cap are tunable via env vars:
    ///   `SAM_PROFILE_PATH` (default `/tmp/blaster-bench/big.sam`, then
    ///     `/tmp/blaster-bench/small.sam`)
    ///   `SAM_PROFILE_RECORDS` (default 5_000_000)
    ///
    /// Run with:
    ///   cargo test --release --lib profile_sam_parse -- --ignored --nocapture
    ///
    /// Ignored by default so CI doesn't depend on the bench fixture.
    #[test]
    #[ignore]
    fn profile_sam_parse() {
        let path = std::env::var("SAM_PROFILE_PATH").unwrap_or_else(|_| {
            let big = "/tmp/blaster-bench/big.sam";
            if std::path::Path::new(big).exists() {
                big.to_string()
            } else {
                "/tmp/blaster-bench/small.sam".to_string()
            }
        });
        let max_records: usize = std::env::var("SAM_PROFILE_RECORDS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(5_000_000);

        eprintln!("loading {path} (cap {max_records} records)...");
        let bytes = match std::fs::read(&path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("skipping profile_sam_parse: {e}");
                return;
            }
        };

        // Parse header via SamReader to set up the tid map.
        let mut reader = SamReader::new(std::io::Cursor::new(bytes.clone()));
        let _ = reader.read_header().expect("header");
        let tid_by_name = reader.tid_by_name.clone();

        // Find where the alignment records begin (first non-@ line).
        let mut body_start = 0usize;
        for (i, &b) in bytes.iter().enumerate() {
            let at_line_start = i == 0 || bytes[i - 1] == b'\n';
            if at_line_start && b != b'@' && b != b'\n' {
                body_start = i;
                break;
            }
        }

        // Index lines as (start, end) offsets into the single body buffer.
        // Avoids 5M small Vec<u8> allocations for a big-input run.
        let mut line_offsets: Vec<(usize, usize)> = Vec::with_capacity(max_records);
        let mut p = body_start;
        while p < bytes.len() && line_offsets.len() < max_records {
            let off = memchr(b'\n', &bytes[p..]).unwrap_or(bytes.len() - p);
            let end = p + off;
            if end > p && bytes[p] != b'@' {
                line_offsets.push((p, end));
            }
            p = end + 1;
        }
        let n = line_offsets.len();
        eprintln!(
            "=== profiling SAM parsers on {n} real records ({:.1} MB body) ===",
            (bytes.len() - body_start) as f64 / 1e6
        );
        let lines = &bytes;
        let line = |i: usize| -> &[u8] {
            let (s, e) = line_offsets[i];
            &lines[s..e]
        };

        // 1) Tab tokenization
        let mut sink = 0usize;
        let mut tabs = [0usize; 11];
        let t = Instant::now();
        for i in 0..n {
            let nt = find_tabs_for_profile(line(i), &mut tabs);
            sink = sink.wrapping_add(nt);
        }
        report("tokenize tabs (memchr_iter)", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink);

        // Pre-tokenize so subsequent benches don't pay the tab cost.
        // (tabs[..nt], nt) per record, flat arrays for cache density.
        let mut all_tabs: Vec<[usize; 11]> = Vec::with_capacity(n);
        let mut all_nt: Vec<u8> = Vec::with_capacity(n);
        for i in 0..n {
            let mut tabs = [0usize; 11];
            let nt = find_tabs_for_profile(line(i), &mut tabs);
            all_tabs.push(tabs);
            all_nt.push(nt as u8);
        }

        fn field<'a>(line: &'a [u8], tabs: &[usize; 11], nt: usize, i: usize) -> &'a [u8] {
            let start = if i == 0 { 0 } else { tabs[i - 1] + 1 };
            let end = if i < nt { tabs[i] } else { line.len() };
            &line[start..end]
        }

        // 2) Integer fields: FLAG (1), POS (3), MAPQ (4), PNEXT (7), TLEN (8)
        let mut sum: i64 = 0;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            sum += parse_u32_ascii(field(l, tabs, nt, 1)).unwrap() as i64;
            sum += parse_i32_ascii(field(l, tabs, nt, 3)).unwrap() as i64;
            sum += parse_u32_ascii(field(l, tabs, nt, 4)).unwrap() as i64;
            sum += parse_i32_ascii(field(l, tabs, nt, 7)).unwrap() as i64;
            sum += parse_i32_ascii(field(l, tabs, nt, 8)).unwrap() as i64;
        }
        report("parse 5 int fields (FLAG/POS/MAPQ/PNEXT/TLEN)", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sum);

        // 3) tid HashMap lookup (RNAME and RNEXT)
        let mut sink_i = 0i64;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            let rname = field(l, tabs, nt, 2);
            let rnext = field(l, tabs, nt, 6);
            sink_i += name_to_tid(rname, &tid_by_name).unwrap() as i64;
            if rnext != b"=" {
                sink_i += name_to_tid(rnext, &tid_by_name).unwrap() as i64;
            }
        }
        report("tid HashMap lookup (RNAME + RNEXT)", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink_i);

        // 4) CIGAR parse
        let mut cigar = Vec::with_capacity(32);
        let mut sink_u = 0usize;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            let cigar_text = field(l, tabs, nt, 5);
            cigar.clear();
            if cigar_text != b"*" && !cigar_text.is_empty() {
                parse_cigar_ops(cigar_text, &mut cigar).unwrap();
            }
            sink_u = sink_u.wrapping_add(cigar.len());
        }
        report("parse_cigar_ops", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink_u);

        // 5) SEQ packing (SIMD)
        let mut packed = Vec::with_capacity(256);
        let mut sink_u = 0usize;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            let seq = field(l, tabs, nt, 9);
            if seq != b"*" {
                packed.clear();
                pack_seq_into(seq, &mut packed);
                sink_u = sink_u.wrapping_add(packed.len());
            }
        }
        report("pack_seq_into (SIMD)", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink_u);

        // 6) QUAL decode (SIMD wrapping_sub of 33)
        let mut quals = Vec::with_capacity(256);
        let mut sink_u8 = 0u8;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            let qual = field(l, tabs, nt, 10);
            if qual != b"*" {
                quals.clear();
                decode_qual_into(qual, &mut quals);
                sink_u8 = sink_u8.wrapping_add(quals.last().copied().unwrap_or(0));
            }
        }
        report("decode_qual_into (SIMD)", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink_u8);

        // 7) Aux tag parsing
        let mut aux = Vec::with_capacity(256);
        let mut sink_u = 0usize;
        let t = Instant::now();
        for i in 0..n {
            let l = line(i);
            let tabs = &all_tabs[i];
            let nt = all_nt[i] as usize;
            let aux_start = if nt >= 11 { tabs[10] + 1 } else { l.len() };
            if aux_start < l.len() {
                aux.clear();
                parse_aux_tags(&l[aux_start..], &mut aux).unwrap();
                sink_u = sink_u.wrapping_add(aux.len());
            }
        }
        report("parse_aux_tags", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(sink_u);

        // 8) Whole parse_sam_line (for baseline)
        let mut rec = fgumi_raw_bam::RawRecord::new();
        let mut cigar = Vec::with_capacity(32);
        let t = Instant::now();
        for i in 0..n {
            parse_sam_line(line(i), &mut rec, &tid_by_name, &mut cigar).unwrap();
        }
        report("full parse_sam_line", t.elapsed().as_nanos() as u64, n);
        std::hint::black_box(rec.len());

        eprintln!("===");
    }

    fn report(label: &str, total_ns: u64, n: usize) {
        let per = total_ns as f64 / n as f64;
        let total_s = total_ns as f64 / 1e9;
        eprintln!("  {:>50}: {:>9.1} ns/rec  ({:.3}s for {} records)", label, per, total_s, n);
    }
}
