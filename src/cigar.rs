//! CIGAR analysis helpers.
//!
//! dupblaster needs four numbers from a CIGAR string: the leading soft/hard
//! clip count (`sclip`), the trailing soft/hard clip count (`eclip`), the
//! reference-consuming aligned length (`ra_len`), and the query-consuming
//! aligned length (`qa_len`).

/// Output of parsing a CIGAR for dupblaster's needs.
#[derive(Debug, Default, Clone, Copy)]
pub struct CigarInfo {
    /// Soft-or-hard clip length at the 5' end of the read (leading clip).
    pub sclip: i32,
    /// Soft-or-hard clip length at the 3' end of the read (trailing clip).
    pub eclip: i32,
    /// Reference-consuming aligned length (sum of M/=/X/D/N op lengths).
    pub ra_len: i32,
    /// Query-consuming aligned length (sum of M/=/X/I op lengths, excludes
    /// clips).
    pub qa_len: i32,
}

impl CigarInfo {
    /// Build from an iterator yielding packed CIGAR ops (`len << 4 | code`,
    /// the packed format stored in a BAM record). `RawRecord::cigar_ops_iter`
    /// yields these by value; `&[u32]` callers pass `slice.iter().copied()`.
    ///
    /// The C++ implementation toggles `first` to false on the first
    /// `M/=/X`; both `S` and `H` accumulate into `sclip` while `first` is
    /// still true, and into `eclip` afterwards. Indel-only CIGARs leave
    /// `first` true throughout, which matches the C++.
    pub fn from_cigar_ops<I: IntoIterator<Item = u32>>(ops: I) -> Self {
        let mut info = Self::default();
        let mut first = true;
        for word in ops {
            let len = (word >> 4) as i32;
            let code = word & 0xf;
            match code {
                // M, =, X
                0 | 7 | 8 => {
                    info.ra_len += len;
                    info.qa_len += len;
                    first = false;
                }
                // S, H
                4 | 5 => {
                    if first {
                        info.sclip += len;
                    } else {
                        info.eclip += len;
                    }
                }
                // D, N
                2 | 3 => info.ra_len += len,
                // I
                1 => info.qa_len += len,
                // P — consumes neither query nor reference
                _ => {}
            }
        }
        info
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Op codes — same values used in the BAM spec.
    const M: u32 = 0;
    const I: u32 = 1;
    const D: u32 = 2;
    const S: u32 = 4;
    const H: u32 = 5;

    /// Build a packed CIGAR op: `(len << 4) | code`.
    const fn op(len: u32, code: u32) -> u32 {
        (len << 4) | code
    }

    #[test]
    fn match_only() {
        let info = CigarInfo::from_cigar_ops([op(100, M)]);
        assert_eq!(info.sclip, 0);
        assert_eq!(info.eclip, 0);
        assert_eq!(info.ra_len, 100);
        assert_eq!(info.qa_len, 100);
    }

    #[test]
    fn leading_soft_clip() {
        let info = CigarInfo::from_cigar_ops([op(5, S), op(95, M)]);
        assert_eq!(info.sclip, 5);
        assert_eq!(info.eclip, 0);
        assert_eq!(info.ra_len, 95);
        assert_eq!(info.qa_len, 95);
    }

    #[test]
    fn trailing_soft_clip() {
        let info = CigarInfo::from_cigar_ops([op(95, M), op(5, S)]);
        assert_eq!(info.sclip, 0);
        assert_eq!(info.eclip, 5);
        assert_eq!(info.ra_len, 95);
        assert_eq!(info.qa_len, 95);
    }

    #[test]
    fn both_clips_with_indel() {
        let info = CigarInfo::from_cigar_ops([
            op(3, S),
            op(40, M),
            op(2, I),
            op(50, M),
            op(5, D),
            op(5, M),
            op(2, H),
        ]);
        assert_eq!(info.sclip, 3);
        assert_eq!(info.eclip, 2);
        assert_eq!(info.ra_len, 100); // M+M+M = 95 + D5 = 100
        assert_eq!(info.qa_len, 97); // M+M+M = 95 + I2 = 97
    }

    #[test]
    fn s_then_h_both_count_as_sclip() {
        let info = CigarInfo::from_cigar_ops([op(3, S), op(2, H), op(50, M)]);
        assert_eq!(info.sclip, 5);
        assert_eq!(info.eclip, 0);
    }

    #[test]
    fn accepts_a_by_value_iterator() {
        // The hot path feeds `RawRecord::cigar_ops_iter` (yields `u32` by
        // value); make sure a `slice.iter().copied()` iterator parses the same.
        let ops = [op(5, S), op(40, M), op(2, I), op(50, M), op(3, S)];
        let info = CigarInfo::from_cigar_ops(ops.iter().copied());
        assert_eq!(info.sclip, 5);
        assert_eq!(info.eclip, 3);
        assert_eq!(info.ra_len, 90); // 40 + 50
        assert_eq!(info.qa_len, 92); // 40 + 2 + 50
    }
}
