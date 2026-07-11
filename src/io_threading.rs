//! Dedicated IO threads with a 16 MB byte ring buffer in each direction.
//!
//! Production dupblaster pipelines look like
//!     aligner | dupblaster | sorter
//! where both ends are typically bursty (the aligner has variable alignment
//! cost per read; the sorter periodically flushes a chunk to disk). With a
//! single-threaded reader+worker+writer, the small OS pipe buffer (~64 KB
//! on macOS) is the only thing decoupling stages, so any blip in one stage
//! stalls the others.
//!
//! [`ThreadedReader`] and [`ThreadedWriter`] put a dedicated thread on
//! each IO end with a 16 MB user-space ring buffer in between. The worker
//! reads/writes through the ring, never blocking on the kernel pipe.
//!
//! Design choices:
//! * **`ringbuf::HeapRb`** for the bytes — single allocation up front,
//!   recycled forever, no per-chunk heap traffic.
//! * **`thread::park` / `unpark`** for blocking. Lock-free fast path when
//!   the ring isn't full/empty; only the rare contended case parks.
//! * **`read()` straight into `vacant_slices_mut`** — one memcpy
//!   (kernel→ring) instead of two (kernel→temp + temp→ring).
//! * **Symmetric on write**: worker pushes bytes into the ring directly;
//!   IO writer thread drains via `as_slices()` + `write_all` + `skip`.

use std::io::{self, BufRead, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::thread::{self, JoinHandle};

use ringbuf::traits::{Consumer, Observer, Producer, Split};
use ringbuf::{HeapCons, HeapProd, HeapRb};

/// Acquire a Mutex even if it's been poisoned by a previous panic.
///
/// Both threads communicate IO errors through the same `Mutex<Option<io::Error>>`,
/// so a panic on one side would otherwise cascade into a panic on the other
/// when it next tries to read or store an error. Treating poison as
/// "no recorded error" lets us surface the original failure instead.
fn lock_or_recover<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

/// Reconstruct an `io::Error` equivalent to `e` without consuming it.
///
/// `io::Error` is `!Clone` (it may wrap an arbitrary boxed payload), so the
/// error slot is surfaced by *copying* rather than draining: every observer
/// gets its own error and the original stays put. This is what makes the slot a
/// write-once, read-many latch — the presence of an error is itself the sticky
/// "this end has failed" state, so no separate flag is needed. OS errors
/// round-trip losslessly (kind + errno + message); for other errors we preserve
/// the kind and the `Display` message, which is the actionable part.
fn clone_io_error(e: &io::Error) -> io::Error {
    match e.raw_os_error() {
        Some(code) => io::Error::from_raw_os_error(code),
        None => io::Error::new(e.kind(), e.to_string()),
    }
}

/// Join an IO thread to completion, ignoring its `()` result.
///
/// A panic inside the IO loop is *not* lost despite the discarded join
/// payload: the [`PanicGuard`] installed in each loop records a fallback
/// `io::Error` and wakes the counterpart before unwinding, so the failure
/// still surfaces through the normal error channel (and `finish()` returns
/// `Err`). Joining here just reaps the thread and orders its teardown.
fn join_io_thread(join: Option<JoinHandle<()>>) {
    if let Some(h) = join {
        let _ = h.join();
    }
}

/// Drop guard that guarantees the counterpart IO end is woken when an IO
/// loop exits — including via an unexpected panic, where the loop's own
/// `unpark` calls never run. Without this, a panic in the IO thread while
/// the worker is parked on a full/empty ring would hang forever.
///
/// On a panicking unwind it also records a fallback error (if none is set) so
/// the panic surfaces as a clean failure rather than silent truncation. That
/// recorded error is the sticky failure state for both ends — the writer has no
/// separate flag, and the reader additionally sets `eof` so a parked reader
/// stops waiting for bytes that will never come.
struct PanicGuard<'a> {
    /// The thread to wake when this guard fires (the counterpart IO end).
    counterpart: &'a thread::Thread,
    /// Shared error slot: records a fallback error message if we're panicking
    /// and the slot is still empty. Never drained, so its presence is what tells
    /// the counterpart this end has failed.
    error: &'a Mutex<Option<io::Error>>,
    /// EOF flag to set on panic so a parked reader stops waiting; `None` on
    /// the write side.
    eof: Option<&'a AtomicBool>,
}

impl Drop for PanicGuard<'_> {
    fn drop(&mut self) {
        if thread::panicking() {
            if let Some(eof) = self.eof {
                eof.store(true, Ordering::Release);
            }
            let mut slot = lock_or_recover(self.error);
            if slot.is_none() {
                *slot = Some(io::Error::other("dupblaster IO thread panicked"));
            }
        }
        self.counterpart.unpark();
    }
}

// ─── Read side ─────────────────────────────────────────────────────────────

/// `BufRead`-compatible reader fed by an IO thread.
pub struct ThreadedReader {
    /// Consumer side of the ring buffer; the worker reads bytes from here.
    consumer: HeapCons<u8>,
    /// Handle to the IO read thread used to call `unpark` when the ring drains.
    io_thread: thread::Thread,
    /// Set by the IO thread once the underlying source reaches EOF.
    eof: Arc<AtomicBool>,
    /// Set by `Drop` to tell the IO thread to exit even if the source isn't done.
    stop: Arc<AtomicBool>,
    /// First read error from the IO thread, if any.
    error: Arc<Mutex<Option<io::Error>>>,
    /// Join handle consumed by `Drop` to reap the IO thread.
    join: Option<JoinHandle<()>>,
}

impl ThreadedReader {
    /// Spawn an IO thread that reads from `src` into a ring buffer of
    /// `ring_bytes` capacity.
    pub fn new<R: Read + Send + 'static>(src: R, ring_bytes: usize) -> Self {
        let rb = HeapRb::<u8>::new(ring_bytes.max(64 * 1024));
        let (producer, consumer) = rb.split();
        let eof = Arc::new(AtomicBool::new(false));
        let stop = Arc::new(AtomicBool::new(false));
        let error = Arc::new(Mutex::new(None));

        let eof_io = eof.clone();
        let stop_io = stop.clone();
        let error_io = error.clone();
        let consumer_thread = thread::current();

        let join = thread::Builder::new()
            .name("dupblaster-io-read".into())
            .spawn(move || io_read_loop(src, producer, eof_io, stop_io, error_io, consumer_thread))
            .expect("spawning IO read thread");
        let io_thread = join.thread().clone();

        Self { consumer, io_thread, eof, stop, error, join: Some(join) }
    }

    /// Copy the stored IO error, if any, leaving it in the slot (see
    /// [`clone_io_error`]). Idempotent: re-reads keep surfacing the failure
    /// rather than masking it as a clean EOF.
    fn peek_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).as_ref().map(clone_io_error)
    }
}

impl Read for ThreadedReader {
    fn read(&mut self, dst: &mut [u8]) -> io::Result<usize> {
        let src = self.fill_buf()?;
        let n = src.len().min(dst.len());
        dst[..n].copy_from_slice(&src[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for ThreadedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            // Surface any IO-thread error first.
            if let Some(e) = self.peek_error() {
                return Err(e);
            }

            if self.consumer.occupied_len() > 0 {
                let (first, _second) = self.consumer.as_slices();
                return Ok(first);
            }

            if self.eof.load(Ordering::Acquire) {
                // Drain any straggler bytes the producer published just
                // before setting EOF. `Acquire` above pairs with
                // `Release` in the IO thread, so re-checking occupied_len
                // here observes any final push.
                if self.consumer.occupied_len() > 0 {
                    continue;
                }
                return Ok(&[]);
            }

            // Ring is empty and producer hasn't flagged EOF yet. Park
            // until the IO thread unparks us. Spurious wakeups are
            // harmless because we re-check the loop condition.
            thread::park();
        }
    }

    fn consume(&mut self, amt: usize) {
        self.consumer.skip(amt);
        // Wake the IO thread in case it parked on a full ring.
        self.io_thread.unpark();
    }
}

impl Drop for ThreadedReader {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Release);
        self.io_thread.unpark();
        join_io_thread(self.join.take());
    }
}

/// Body of the dedicated read IO thread. Pumps bytes from `src` into the ring
/// buffer, parking when the ring is full, and waking the consumer on each push
/// or at EOF/error. The [`PanicGuard`] ensures the consumer is always woken on
/// exit, even if the thread panics.
fn io_read_loop<R: Read>(
    mut src: R,
    mut producer: HeapProd<u8>,
    eof: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    error: Arc<Mutex<Option<io::Error>>>,
    consumer_thread: thread::Thread,
) {
    // Wake the consumer on any exit (incl. panic) so it never parks against
    // a dead producer; on panic also flag EOF + record a fallback error.
    let _guard = PanicGuard { counterpart: &consumer_thread, error: &error, eof: Some(&*eof) };
    loop {
        if stop.load(Ordering::Acquire) {
            break;
        }

        let (first, _second) = producer.vacant_slices_mut();
        if first.is_empty() {
            // Wake the consumer in case it's waiting (and we've just
            // become full because it's been slow).
            consumer_thread.unpark();
            thread::park();
            continue;
        }

        // SAFETY: We're reinterpreting `&mut [MaybeUninit<u8>]` as
        // `&mut [u8]` only to pass to `Read::read`, which writes into
        // every byte it claims to have read. After the read returns
        // `Ok(n)`, exactly `n` bytes are initialized; `advance_write_index(n)`
        // exposes only those to the consumer.
        let dst: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(first.as_mut_ptr() as *mut u8, first.len()) };
        match src.read(dst) {
            Ok(0) => {
                eof.store(true, Ordering::Release);
                consumer_thread.unpark();
                break;
            }
            Ok(n) => {
                // SAFETY: `read` initialized exactly `n` bytes of the
                // vacant slice; ringbuf requires that pre-condition.
                unsafe {
                    producer.advance_write_index(n);
                }
                consumer_thread.unpark();
            }
            Err(e) => {
                *lock_or_recover(&error) = Some(e);
                eof.store(true, Ordering::Release);
                consumer_thread.unpark();
                break;
            }
        }
    }
    // Final wake so a consumer parked on `eof=false` sees the new state.
    consumer_thread.unpark();
}

// ─── Write side ────────────────────────────────────────────────────────────

/// `Write`-compatible writer that hands bytes off to an IO thread.
pub struct ThreadedWriter {
    /// Producer side of the ring buffer; the worker pushes bytes here.
    producer: HeapProd<u8>,
    /// Handle to the IO write thread used to call `unpark` when new bytes are ready.
    io_thread: thread::Thread,
    /// Set by `finish()` or `Drop` to signal the IO thread that no more data is coming.
    finished: Arc<AtomicBool>,
    /// The IO thread's write error, if any. Written once by the IO thread and
    /// never drained, so its presence is the sticky "this writer has failed"
    /// state: every `write`/`flush`/`finish` copies it out (see [`peek_error`])
    /// and rejects the operation, and a failed writer can never accept more
    /// bytes into a ring its exited IO thread would never drain.
    ///
    /// [`peek_error`]: ThreadedWriter::peek_error
    error: Arc<Mutex<Option<io::Error>>>,
    /// Join handle consumed by `Drop` to reap the IO thread.
    join: Option<JoinHandle<()>>,
}

impl ThreadedWriter {
    /// Spawn an IO thread that writes the ring contents to `dst`. Ring
    /// holds `ring_bytes` of pending output.
    pub fn new<W: Write + Send + 'static>(dst: W, ring_bytes: usize) -> Self {
        let rb = HeapRb::<u8>::new(ring_bytes.max(64 * 1024));
        let (producer, consumer) = rb.split();
        let finished = Arc::new(AtomicBool::new(false));
        let error = Arc::new(Mutex::new(None));

        let finished_io = finished.clone();
        let error_io = error.clone();
        let producer_thread = thread::current();

        let join = thread::Builder::new()
            .name("dupblaster-io-write".into())
            .spawn(move || io_write_loop(dst, consumer, finished_io, error_io, producer_thread))
            .expect("spawning IO write thread");
        let io_thread = join.thread().clone();

        Self { producer, io_thread, finished, error, join: Some(join) }
    }

    /// Copy the stored IO error, if any, leaving it in the slot (see
    /// [`clone_io_error`]). The slot is never drained, so this keeps surfacing
    /// the failure on every call — including the re-entrant writes `Drop` makes.
    fn peek_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).as_ref().map(clone_io_error)
    }

    /// Flush remaining bytes, signal the IO thread to drain, then join.
    /// Returns the IO thread's final result. Idempotent — calling twice
    /// is a no-op the second time.
    pub fn finish(mut self) -> io::Result<()> {
        self.finished.store(true, Ordering::Release);
        self.io_thread.unpark();
        join_io_thread(self.join.take());
        if let Some(e) = self.peek_error() {
            return Err(e);
        }
        Ok(())
    }
}

impl Write for ThreadedWriter {
    fn write(&mut self, mut buf: &[u8]) -> io::Result<usize> {
        // Surface any IO-thread error before touching the ring, and only once
        // per call rather than per ring-push iteration (the previous
        // per-iteration check locked the `Mutex` ~150 M times on a 30 GB run).
        // Rejecting a write on a failed writer *up front* is what prevents the
        // deadlock: otherwise we push into a ring the exited IO thread will
        // never drain, fill it, and `park` forever with no one left to `unpark`
        // us. This is exactly the re-entry `bgzf::Writer::Drop` performs (flush
        // buffered blocks + write the EOF marker) after a broken pipe surfaced.
        if let Some(e) = self.peek_error() {
            return Err(e);
        }
        let initial_len = buf.len();
        while !buf.is_empty() {
            let pushed = self.producer.push_slice(buf);
            if pushed > 0 {
                buf = &buf[pushed..];
                self.io_thread.unpark();
            } else {
                // Ring is full — let the IO thread drain. Park; the
                // IO thread will unpark us after it writes.
                thread::park();
                // The IO thread also unparks us when it dies on a write
                // error, leaving the ring permanently full. Without this
                // re-check we'd loop forever pushing into a ring nobody
                // drains. Surfacing the error here both reports the failure
                // and breaks the deadlock.
                if let Some(e) = self.peek_error() {
                    return Err(e);
                }
            }
        }
        Ok(initial_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Nothing to flush at this layer — bytes are already in the ring
        // or written to `dst`. The IO thread does its own write_all.
        if let Some(e) = self.peek_error() {
            return Err(e);
        }
        Ok(())
    }
}

impl Drop for ThreadedWriter {
    fn drop(&mut self) {
        // If finish() wasn't called, signal anyway so the IO thread can
        // shut down cleanly. Errors are silently dropped here — explicit
        // finish() is the right path for callers who care.
        self.finished.store(true, Ordering::Release);
        self.io_thread.unpark();
        join_io_thread(self.join.take());
    }
}

/// Body of the dedicated write IO thread. Drains the ring buffer into `dst`,
/// parking when the ring is empty, and waking the producer after each drain to
/// signal available space. The [`PanicGuard`] ensures the producer is always
/// woken on exit so it doesn't block forever against a dead consumer.
fn io_write_loop<W: Write>(
    mut dst: W,
    mut consumer: HeapCons<u8>,
    finished: Arc<AtomicBool>,
    error: Arc<Mutex<Option<io::Error>>>,
    producer_thread: thread::Thread,
) {
    // Wake the producer on any exit (incl. panic) so it never parks against
    // a dead consumer; on panic also record a fallback error.
    let _guard = PanicGuard { counterpart: &producer_thread, error: &error, eof: None };
    // Record an IO error into the shared slot, then wake the producer (which may
    // be parked on a full ring) so it surfaces the failure instead of blocking
    // forever. The slot is never drained, so this single write latches the
    // failure for every future `write`/`flush`/`finish`.
    let record_error = |e: io::Error| {
        *lock_or_recover(&error) = Some(e);
        producer_thread.unpark();
    };

    loop {
        if consumer.occupied_len() > 0 {
            let (first, _second) = consumer.as_slices();
            // Copy locally because `skip` borrows consumer mutably below.
            let n = first.len();
            if let Err(e) = dst.write_all(first) {
                record_error(e);
                break;
            }
            consumer.skip(n);
            producer_thread.unpark();
            continue;
        }
        if finished.load(Ordering::Acquire) {
            // Drain any stragglers — check again under acquire ordering.
            if consumer.occupied_len() > 0 {
                continue;
            }
            // Flush the underlying writer before exit.
            if let Err(e) = dst.flush() {
                record_error(e);
            }
            break;
        }
        thread::park();
    }
    producer_thread.unpark();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Round-trip a payload through a ThreadedReader: bytes in == bytes out.
    #[test]
    fn threaded_reader_round_trip_small() {
        let payload: Vec<u8> = (0..1000u32).flat_map(|i| i.to_le_bytes()).collect();
        let mut r = ThreadedReader::new(Cursor::new(payload.clone()), 64 * 1024);
        let mut out = Vec::new();
        std::io::copy(&mut r, &mut out).unwrap();
        assert_eq!(out, payload);
    }

    /// Payload much larger than the ring buffer — exercises the wrap-around.
    #[test]
    fn threaded_reader_round_trip_larger_than_ring() {
        let ring = 4096;
        let payload: Vec<u8> = (0..(ring * 8) as u32).map(|i| i as u8).collect();
        let mut r = ThreadedReader::new(Cursor::new(payload.clone()), ring);
        let mut out = Vec::new();
        std::io::copy(&mut r, &mut out).unwrap();
        assert_eq!(out, payload);
    }

    /// Write a payload through ThreadedWriter and confirm the underlying
    /// sink received every byte after `finish()`.
    #[test]
    fn threaded_writer_round_trip_with_finish() {
        // ThreadedWriter takes ownership of `W: Write + Send + 'static`, so
        // we hand it a `Sink` that mirrors bytes into a shared buffer the
        // test can inspect.
        struct Sink(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for Sink {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let payload: Vec<u8> = (0..50_000u32).map(|i| i as u8).collect();
        let captured = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
        let mut w = ThreadedWriter::new(Sink(captured.clone()), 4096);
        w.write_all(&payload).unwrap();
        w.finish().unwrap();
        assert_eq!(*captured.lock().unwrap(), payload);
    }

    /// A sink that always fails its writes.
    struct FailingSink;
    impl Write for FailingSink {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "downstream closed"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A sink-write failure must surface as an `Err` rather than hanging the
    /// producer forever against a full, never-draining ring. Regression test
    /// for the missing error re-check after `park()` in `write`.
    #[test]
    fn threaded_writer_surfaces_sink_error_without_deadlock() {
        // Payload ≫ ring (clamped to a 64 KiB minimum) forces the producer
        // to fill the ring and park while the IO thread dies on its first
        // write.
        let ring = 4096;
        let payload = vec![0u8; 1024 * 1024];
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            let mut w = ThreadedWriter::new(FailingSink, ring);
            let result = w.write_all(&payload).and_then(|()| w.finish());
            tx.send(result.is_err()).unwrap();
        });
        let errored = rx
            .recv_timeout(std::time::Duration::from_secs(10))
            .expect("writer deadlocked on a failing sink");
        worker.join().unwrap();
        assert!(errored, "a failing sink must surface as an error");
    }

    /// Once a `write` has surfaced the rich error, `finish()` must also report
    /// failure rather than falsely returning `Ok`. The error slot is a latch
    /// that is never drained, so both calls observe it.
    #[test]
    fn threaded_writer_finish_fails_after_error_already_surfaced() {
        let ring = 4096;
        let payload = vec![7u8; 1024 * 1024];
        let mut w = ThreadedWriter::new(FailingSink, ring);
        // The oversized payload guarantees the producer parks and observes
        // the error, so the first `write_all` returns `Err`.
        assert!(w.write_all(&payload).is_err());
        assert!(w.finish().is_err(), "finish must report the latched error too");
    }

    /// Number of iterations for the re-entrant-write stress test. Kept modest so
    /// the default suite stays fast (each iteration spins up and tears down a
    /// fresh IO thread); the deadlock (pre-fix) reproduces on the very first
    /// iteration, and the fix has been verified locally across hundreds of
    /// thousands of iterations via the `REENTRANT_STRESS_ITERS` override below.
    const REENTRANT_STRESS_ITERS: usize = 2_000;

    /// Re-entrant writes performed after the error is first surfaced. bgzf's
    /// `Drop` re-enters the sink several times (one write per buffered block,
    /// then the EOF marker, then a flush). The dying IO thread leaves at most a
    /// couple of stray `unpark` tokens, each of which can rescue one parked
    /// write; this count comfortably exceeds them so the pre-fix hang would
    /// reproduce reliably rather than depending on `unpark` timing.
    const REENTRANT_WRITES: usize = 16;

    /// Run `body` on its own thread and wait up to 10 s for it to finish.
    ///
    /// Returns `Some(value)` if `body` completes, or `None` if it *deadlocks*
    /// (parks forever). A non-deadlocked body finishes in microseconds, so the
    /// 10 s bound only ever elapses on a genuine hang. On deadlock we repeatedly
    /// `unpark` the stuck thread — each nudge lets `ThreadedWriter::write` fall
    /// through to its post-park error re-check and return, and a body that
    /// re-parks (a loop of writes, like bgzf's multi-write `Drop`) needs several
    /// — then join it, so the stress test never leaks threads across iterations.
    fn run_or_detect_deadlock<T, F>(body: F) -> Option<T>
    where
        T: Send + 'static,
        F: FnOnce() -> T + Send + 'static,
    {
        use std::sync::mpsc::RecvTimeoutError;
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            let _ = tx.send(body());
        });
        let handle = worker.thread().clone();
        // Fast path: a healthy body reports back near-instantly.
        if let Ok(value) = rx.recv_timeout(std::time::Duration::from_secs(10)) {
            worker.join().unwrap();
            return Some(value);
        }
        // Deadlock path: keep nudging until the body unwinds and reports, then
        // reap it. Bounded so a hang `unpark` can't clear doesn't wedge the
        // suite — we detach and report rather than block the test binary.
        for _ in 0..1000 {
            handle.unpark();
            match rx.recv_timeout(std::time::Duration::from_millis(20)) {
                Ok(_) => {
                    worker.join().unwrap();
                    return None;
                }
                Err(RecvTimeoutError::Timeout) => continue,
                Err(RecvTimeoutError::Disconnected) => {
                    worker.join().unwrap();
                    return None;
                }
            }
        }
        std::mem::forget(worker);
        None
    }

    /// Drive a `ThreadedWriter` over a failing sink to the exact deadlock state:
    /// surface the write error, leaving the ring full and the IO thread gone,
    /// then re-enter `write` repeatedly the way `bgzf::Writer::Drop` does.
    /// Returns `true` iff every write surfaced an error. Without the fix one of
    /// the re-entrant writes parks forever, so this never returns and the
    /// watchdog reports a deadlock instead.
    fn reentrant_write_after_error_surfaced(ring: usize) -> bool {
        let mut w = ThreadedWriter::new(FailingSink, ring);
        // Payload ≥ ring: the producer fills the ring and parks while the IO
        // thread dies on its first write, so this surfaces the error and leaves
        // the ring full with the IO thread already exited.
        let mut all_errored = w.write_all(&vec![0u8; 2 * ring]).is_err();
        for _ in 0..REENTRANT_WRITES {
            all_errored &= w.write_all(b"trailing BGZF EOF marker").is_err();
        }
        all_errored
    }

    /// After a sink write fails and that error has been surfaced once, further
    /// writes must return `Err` rather than park forever.
    ///
    /// This reproduces the re-entry `bgzf::Writer::Drop` performs while
    /// unwinding: it flushes its block buffer and writes the BGZF EOF marker
    /// back through this writer *after* a broken-pipe error already propagated
    /// out. If a re-entrant write parks on a full ring whose IO thread has
    /// exited, the process deadlocks in `Drop` — which is how a released
    /// dupblaster hung for hours when its downstream sorter died on ENOSPC.
    #[test]
    fn threaded_writer_reentrant_write_after_error_does_not_deadlock() {
        let all_errored =
            run_or_detect_deadlock(|| reentrant_write_after_error_surfaced(64 * 1024)).expect(
                "re-entrant write deadlocked: parked on a full ring whose IO thread had exited",
            );
        assert!(all_errored, "every write after a surfaced error must error, not succeed");
    }

    /// Stress the re-entrant-write path across many IO-thread lifetimes to shake
    /// out the timing-dependent hang. Each iteration spins up a fresh writer +
    /// IO thread, kills it via a failing sink, surfaces the error, then
    /// re-enters. The pre-fix deadlock reproduces on the first iteration; the fix
    /// must hold across all of them.
    ///
    /// Override the iteration count with `REENTRANT_STRESS_ITERS=<n>` to grind
    /// on it harder (e.g. hundreds of thousands) when auditing the fix locally.
    #[test]
    fn threaded_writer_reentrant_write_stress_no_deadlock() {
        let iters = std::env::var("REENTRANT_STRESS_ITERS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(REENTRANT_STRESS_ITERS);
        for i in 0..iters {
            let all_errored =
                run_or_detect_deadlock(|| reentrant_write_after_error_surfaced(64 * 1024))
                    .unwrap_or_else(|| panic!("iteration {i}: re-entrant write deadlocked"));
            assert!(all_errored, "iteration {i}: every re-entrant write must surface an error");
        }
    }

    /// A panic inside the IO write thread must not deadlock the producer and
    /// must surface as a failure (not a silently-successful `finish`).
    #[test]
    fn threaded_writer_panic_surfaces_without_deadlock() {
        struct PanicSink;
        impl Write for PanicSink {
            fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
                panic!("sink panicked");
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let ring = 4096;
        let payload = vec![1u8; 1024 * 1024];
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            let mut w = ThreadedWriter::new(PanicSink, ring);
            let result = w.write_all(&payload).and_then(|()| w.finish());
            tx.send(result.is_err()).unwrap();
        });
        let errored = rx
            .recv_timeout(std::time::Duration::from_secs(10))
            .expect("writer deadlocked on a panicking IO thread");
        worker.join().unwrap();
        assert!(errored, "an IO-thread panic must surface as an error");
    }
}
