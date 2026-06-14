//! Microbenchmark for the dedup-table prefetch experiment.
//!
//! Approximates the per-block hash-insert workload at 12× steady state:
//! - ~4096 partition cells (matches min-bins=32 → ~32² partitions, 4× strands)
//! - ~40,000 entries per cell at steady state
//! - 100M (cell_off, sig) ops with random cell selection
//!
//! Compares three modes:
//!  1. baseline: per-op `insert(sig)`
//!  2. batched-warmup: chunks of K, `contains(sig)` pass then `insert(sig)` pass
//!  3. batched-prefetch: chunks of K, prefetch then `insert(sig)`
//!
//! Mode 2 issues real probe loads in parallel via the OOO engine; mode 3 uses
//! explicit `prfm`/`prefetch_t0` to bring the U64Set struct cache line in.

use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

/// Identity-style hasher for `u64` keys using a single golden-ratio multiply.
/// Matches the hasher used in `sig.rs` so the benchmark exercises the same
/// cache-miss pattern as the real dup table.
#[derive(Default)]
struct U64Hasher(u64);
impl Hasher for U64Hasher {
    #[inline(always)]
    fn write_u64(&mut self, n: u64) {
        self.0 = n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!()
    }
}
/// A `HashSet<u64>` backed by the golden-ratio [`U64Hasher`] — approximates one
/// cell of the partitioned dup table in `sig.rs`.
type U64Set = HashSet<u64, BuildHasherDefault<U64Hasher>>;

/// 64-bit xorshift PRNG — fast, non-cryptographic, sufficient for generating
/// random cell offsets and signatures in the benchmark without allocation.
#[inline(always)]
fn xorshift(s: &mut u64) -> u64 {
    *s ^= *s << 13;
    *s ^= *s >> 7;
    *s ^= *s << 17;
    *s
}

/// Issue a hardware prefetch hint for the cache line holding `*ptr`.
///
/// On AArch64 emits `prfm pldl1keep`; on x86-64 emits `_MM_HINT_T0`. No-op on
/// other targets (the `unused_variables` lint is suppressed to keep the call
/// site uniform across platforms).
#[inline(always)]
#[allow(unused_variables)]
fn prefetch_read<T>(ptr: *const T) {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) ptr,
            options(nostack, preserves_flags, readonly),
        );
    }
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch(ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
    }
}

/// Build `n_cells` [`U64Set`] cells each pre-populated with `steady_size` random
/// entries seeded from `seed_start`, replicating the steady-state dup table.
fn populate(n_cells: usize, steady_size: usize, seed_start: u64) -> Vec<U64Set> {
    let mut sets: Vec<U64Set> =
        (0..n_cells).map(|_| U64Set::with_capacity_and_hasher(64, Default::default())).collect();
    let mut seed = seed_start;
    for cell in &mut sets {
        for _ in 0..steady_size {
            cell.insert(xorshift(&mut seed));
        }
    }
    sets
}

fn main() {
    const N_CELLS: usize = 4096;
    const STEADY_SIZE: usize = 40_000;
    const N_OPS: usize = 50_000_000;
    let batches = [4usize, 8, 16, 32];

    eprintln!("Populating {} cells × {} entries...", N_CELLS, STEADY_SIZE);
    let t0 = Instant::now();
    let template = populate(N_CELLS, STEADY_SIZE, 0xDEAD_BEEF);
    eprintln!("  ({:?})", t0.elapsed());

    eprintln!("Generating {} ops...", N_OPS);
    let mut ops: Vec<(u32, u64)> = Vec::with_capacity(N_OPS);
    let mut seed = 0xABCD_0123_u64;
    for _ in 0..N_OPS {
        let off = (xorshift(&mut seed) % N_CELLS as u64) as u32;
        let sig = xorshift(&mut seed);
        ops.push((off, sig));
    }

    // Baseline.
    let mut sets = template.clone();
    let mut dup_count = 0usize;
    let t1 = Instant::now();
    for &(off, sig) in &ops {
        if !sets[off as usize].insert(sig) {
            dup_count += 1;
        }
    }
    let baseline = t1.elapsed();
    eprintln!("\nbaseline:           {:.3}s ({} dups)", baseline.as_secs_f64(), dup_count);

    // Batched contains-warmup.
    for &batch in &batches {
        let mut sets = template.clone();
        let mut dup_count2 = 0usize;
        let t = Instant::now();
        for chunk in ops.chunks(batch) {
            for &(off, sig) in chunk {
                let _ = sets[off as usize].contains(&sig);
            }
            for &(off, sig) in chunk {
                if !sets[off as usize].insert(sig) {
                    dup_count2 += 1;
                }
            }
        }
        let elapsed = t.elapsed();
        let speedup =
            100.0 * (baseline.as_secs_f64() - elapsed.as_secs_f64()) / baseline.as_secs_f64();
        eprintln!(
            "contains-warmup B={:2}: {:.3}s  ({:+.2}%)  dups={}",
            batch,
            elapsed.as_secs_f64(),
            speedup,
            dup_count2,
        );
        assert_eq!(dup_count, dup_count2);
    }

    // Batched prefetch on the U64Set struct address.
    for &batch in &batches {
        let mut sets = template.clone();
        let mut dup_count3 = 0usize;
        let t = Instant::now();
        for chunk in ops.chunks(batch) {
            for &(off, _) in chunk {
                prefetch_read(&sets[off as usize] as *const U64Set);
            }
            for &(off, sig) in chunk {
                if !sets[off as usize].insert(sig) {
                    dup_count3 += 1;
                }
            }
        }
        let elapsed = t.elapsed();
        let speedup =
            100.0 * (baseline.as_secs_f64() - elapsed.as_secs_f64()) / baseline.as_secs_f64();
        eprintln!(
            "prefetch-struct B={:2}: {:.3}s  ({:+.2}%)  dups={}",
            batch,
            elapsed.as_secs_f64(),
            speedup,
            dup_count3,
        );
        assert_eq!(dup_count, dup_count3);
    }
}
