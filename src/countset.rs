//! A struct-of-arrays (SoA) SwissTable that tracks a per-key occurrence count.
//!
//! This is the storage cell of the `--duplication-spectrum` count side-table (see
//! [`crate::counts`]): a set of signatures, each with a `u16` count, used to
//! build the group-size histogram (η_k). It only ever holds signatures seen
//! **≥2×** — singletons are recovered by subtraction — so it stays small.
//!
//! Layout is SoA (`ctrl` + `keys` + `counts` in three arrays) rather than a
//! `HashMap<K, u16>` because a `HashMap`'s `(K, u16)` bucket pads to the key's
//! alignment (a `(u64, u16)` bucket is 16 B; a `(u128, u16)` bucket 32 B),
//! wasting most of the payload. SoA packs a `u64` key + `u16` count + `u8` ctrl
//! byte into ~11 B/slot. The probe is a SwissTable group scan using the **SWAR**
//! strategy (8-wide, pure integer) — portable and exactly what hashbrown itself
//! uses off x86, so no SIMD dependency is pulled in.
//!
//! Keys are hashed with dupblaster's own [`crate::sig::U64Hasher`] /
//! [`crate::sig::U32Hasher`] (the golden-ratio multiply the dedup tables use),
//! invoked here rather than reimplemented.

use std::hash::Hasher;

use crate::sig::{U32Hasher, U64Hasher};

/// Empty control byte (top bit set; a stored fingerprint has its top bit clear).
const EMPTY: u8 = 0xFF;
/// Initial capacity on first growth (power of two, ≥ `GROUP`).
const INITIAL_CAP: usize = 64;
/// SWAR group width in bytes.
const GROUP: usize = 8;
const SWAR_LO: u64 = 0x0101_0101_0101_0101;
const SWAR_HI: u64 = 0x8080_8080_8080_8080;

/// A key type storable in a [`CountSet`]: hashable via dupblaster's specialized
/// hasher and default-constructible (for the backing key array).
pub(crate) trait CountKey: Copy + Eq + Default {
    /// Hash through the same golden-ratio multiply the dedup tables use.
    fn ghash(self) -> u64;
}

impl CountKey for u64 {
    #[inline(always)]
    fn ghash(self) -> u64 {
        let mut h = U64Hasher::default();
        h.write_u64(self);
        h.finish()
    }
}

impl CountKey for u32 {
    #[inline(always)]
    fn ghash(self) -> u64 {
        let mut h = U32Hasher::default();
        h.write_u32(self);
        h.finish()
    }
}

/// Load `GROUP` control bytes at `pos` as a little-endian `u64`. The control
/// array is sized `cap + GROUP - 1` with the first `GROUP - 1` bytes mirrored at
/// the end, so a group load at any `pos < cap` stays in bounds.
#[inline(always)]
fn load_group(ctrl: &[u8], pos: usize) -> u64 {
    u64::from_le_bytes(ctrl[pos..pos + GROUP].try_into().unwrap())
}

/// SWAR "which lanes equal `byte`": returns a word with `0x80` in each matching
/// byte lane, `0x00` elsewhere (classic zero-byte detection after XOR).
#[inline(always)]
fn match_lanes(group: u64, byte: u8) -> u64 {
    let x = group ^ (SWAR_LO.wrapping_mul(byte as u64));
    x.wrapping_sub(SWAR_LO) & !x & SWAR_HI
}

/// Lowest set lane of a `match_lanes` mask (0..GROUP).
#[inline(always)]
fn lowest_lane(mask: u64) -> usize {
    (mask.trailing_zeros() >> 3) as usize
}

/// A struct-of-arrays open-addressing set of `K` keys, each carrying a `u16`
/// occurrence count. Counts start at 2 (an entry is only ever created on a key's
/// *second* observation) and saturate at [`u16::MAX`].
pub(crate) struct CountSet<K> {
    /// Control bytes (`cap + GROUP - 1`; the leading `GROUP - 1` mirrored).
    ctrl: Vec<u8>,
    /// Backing key slots (`cap`).
    keys: Vec<K>,
    /// Backing count slots (`cap`), pre-initialized to 2.
    counts: Vec<u16>,
    cap: usize,
    mask: usize,
    len: usize,
}

impl<K: CountKey> CountSet<K> {
    /// An empty set (allocates nothing until the first [`Self::bump`]).
    pub(crate) fn new() -> Self {
        Self { ctrl: Vec::new(), keys: Vec::new(), counts: Vec::new(), cap: 0, mask: 0, len: 0 }
    }

    /// The occurrence count of each held key, in arbitrary order. Feeds the η_k
    /// histogram (each yielded `c` is one signature observed exactly `c` times).
    pub(crate) fn counts(&self) -> impl Iterator<Item = u16> + '_ {
        (0..self.cap).filter(move |&i| self.ctrl[i] & 0x80 == 0).map(move |i| self.counts[i])
    }

    /// Record one *repeat* observation of `key`: create it with count 2 on its
    /// first repeat (an entry only ever appears on a key's second observation),
    /// otherwise increment (saturating at [`u16::MAX`]).
    #[inline]
    pub(crate) fn bump(&mut self, key: K) {
        if (self.len + 1) * 8 > self.cap * 7 {
            self.grow();
        }
        let h = key.ghash();
        let h2 = ((h >> 57) as u8) & 0x7F;
        let mut pos = (h as usize) & self.mask;
        let mut stride = 0usize;
        loop {
            let group = load_group(&self.ctrl, pos);
            let mut m = match_lanes(group, h2);
            while m != 0 {
                let i = (pos + lowest_lane(m)) & self.mask;
                if self.keys[i] == key {
                    self.counts[i] = self.counts[i].saturating_add(1);
                    return;
                }
                m &= m - 1;
            }
            let e = match_lanes(group, EMPTY);
            if e != 0 {
                let i = (pos + lowest_lane(e)) & self.mask;
                self.set_ctrl(i, h2);
                self.keys[i] = key;
                // counts[i] stays at its pre-initialized 2 (this is key's 2nd sighting).
                self.len += 1;
                return;
            }
            stride += GROUP;
            pos = (pos + stride) & self.mask;
        }
    }

    /// Write a control byte, mirroring it into the trailing wraparound region
    /// when it lands in the first `GROUP - 1` slots.
    #[inline(always)]
    fn set_ctrl(&mut self, i: usize, v: u8) {
        self.ctrl[i] = v;
        if i < GROUP - 1 {
            self.ctrl[i + self.cap] = v;
        }
    }

    /// Double the capacity (or allocate [`INITIAL_CAP`] from empty) and re-insert
    /// every live entry, carrying its count across.
    fn grow(&mut self) {
        let new_cap = if self.cap == 0 { INITIAL_CAP } else { self.cap * 2 };
        let new_mask = new_cap - 1;
        let mut nctrl = vec![EMPTY; new_cap + GROUP - 1];
        let mut nkeys = vec![K::default(); new_cap];
        let mut ncounts = vec![2u16; new_cap];
        for i in 0..self.cap {
            if self.ctrl[i] & 0x80 != 0 {
                continue; // empty slot
            }
            let (key, cnt) = (self.keys[i], self.counts[i]);
            let h = key.ghash();
            let h2 = ((h >> 57) as u8) & 0x7F;
            let mut pos = (h as usize) & new_mask;
            let mut stride = 0usize;
            loop {
                let e = match_lanes(load_group(&nctrl, pos), EMPTY);
                if e != 0 {
                    let j = (pos + lowest_lane(e)) & new_mask;
                    nctrl[j] = h2;
                    if j < GROUP - 1 {
                        nctrl[j + new_cap] = h2;
                    }
                    nkeys[j] = key;
                    ncounts[j] = cnt;
                    break;
                }
                stride += GROUP;
                pos = (pos + stride) & new_mask;
            }
        }
        self.ctrl = nctrl;
        self.keys = nkeys;
        self.counts = ncounts;
        self.cap = new_cap;
        self.mask = new_mask;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// splitmix64 — a cheap deterministic key generator for the tests.
    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[test]
    fn count_is_observations_not_bumps() {
        // A key observed 5× is bumped on observations 2..=5, i.e. 4 times, and
        // must report a count of 5 (the first bump creates it at 2).
        let mut cs: CountSet<u64> = CountSet::new();
        for _ in 0..4 {
            cs.bump(42);
        }
        assert_eq!(cs.counts().count(), 1);
        assert_eq!(cs.counts().next().unwrap(), 5);
    }

    #[test]
    fn bump_matches_a_shadow_map_over_many_colliding_keys() {
        let mut cs: CountSet<u64> = CountSet::new();
        let mut shadow: HashMap<u64, u32> = HashMap::new(); // bumps per key
        let mut st = 0xDEAD_BEEF_u64;
        for _ in 0..300_000 {
            let key = splitmix(&mut st) & 0x1FFF; // ~8k distinct → heavy repetition
            cs.bump(key);
            *shadow.entry(key).or_insert(0) += 1;
        }
        assert_eq!(cs.counts().count(), shadow.len(), "distinct-key count");
        let mut got: Vec<u16> = cs.counts().collect();
        // stored count = bumps + 1 (an entry starts at 2 after its first bump).
        let mut want: Vec<u16> =
            shadow.values().map(|&b| (b + 1).min(u16::MAX as u32) as u16).collect();
        got.sort_unstable();
        want.sort_unstable();
        assert_eq!(got, want, "per-key counts");
    }

    #[test]
    fn counts_saturate_at_u16_max() {
        let mut cs: CountSet<u32> = CountSet::new();
        for _ in 0..70_000 {
            cs.bump(7u32);
        }
        assert_eq!(cs.counts().count(), 1);
        assert_eq!(cs.counts().next().unwrap(), u16::MAX);
    }

    #[test]
    fn u32_keys_stay_distinct() {
        let mut cs: CountSet<u32> = CountSet::new();
        cs.bump(1);
        cs.bump(1);
        cs.bump(2);
        let mut got: Vec<u16> = cs.counts().collect();
        got.sort_unstable();
        assert_eq!(cs.counts().count(), 2);
        assert_eq!(got, vec![2, 3]); // key 2 bumped once → 2; key 1 bumped twice → 3
    }
}
