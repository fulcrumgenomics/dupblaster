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
/// On a panicking unwind it also records a fallback error (if none is set)
/// so the panic surfaces as a clean failure rather than silent truncation.
struct PanicGuard<'a> {
    /// The thread to wake when this guard fires (the counterpart IO end).
    counterpart: &'a thread::Thread,
    /// Shared error slot: records a fallback error message if we're panicking
    /// and the slot is still empty.
    error: &'a Mutex<Option<io::Error>>,
    /// Sticky error flag to set on panic (writer side); `None` on the read
    /// side, which signals failure via `eof` + the error slot instead.
    errored: Option<&'a AtomicBool>,
    /// EOF flag to set on panic so a parked reader stops waiting; `None` on
    /// the write side.
    eof: Option<&'a AtomicBool>,
}

impl Drop for PanicGuard<'_> {
    fn drop(&mut self) {
        if thread::panicking() {
            if let Some(errored) = self.errored {
                errored.store(true, Ordering::Release);
            }
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

    /// Take the stored IO error out of the shared slot, leaving it empty.
    fn take_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).take()
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
            if let Some(e) = self.take_error() {
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
    let _guard = PanicGuard {
        counterpart: &consumer_thread,
        error: &error,
        errored: None,
        eof: Some(&*eof),
    };
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
    /// First write error from the IO thread, if any.
    error: Arc<Mutex<Option<io::Error>>>,
    /// Sticky "an error occurred" flag. The rich `error` is `take`n by the
    /// first `write`/`flush` that surfaces it, but this flag stays set so
    /// `finish()` can never report success after a failed write even if the
    /// error was already drained. Set by the IO thread *before* it stores
    /// the error, so observing this flag implies the error slot is (or was)
    /// populated.
    errored: Arc<AtomicBool>,
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
        let errored = Arc::new(AtomicBool::new(false));

        let finished_io = finished.clone();
        let error_io = error.clone();
        let errored_io = errored.clone();
        let producer_thread = thread::current();

        let join = thread::Builder::new()
            .name("dupblaster-io-write".into())
            .spawn(move || {
                io_write_loop(dst, consumer, finished_io, error_io, errored_io, producer_thread)
            })
            .expect("spawning IO write thread");
        let io_thread = join.thread().clone();

        Self { producer, io_thread, finished, error, errored, join: Some(join) }
    }

    /// Take the stored IO error out of the shared slot, leaving it empty.
    fn take_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).take()
    }

    /// Flush remaining bytes, signal the IO thread to drain, then join.
    /// Returns the IO thread's final result. Idempotent — calling twice
    /// is a no-op the second time.
    pub fn finish(mut self) -> io::Result<()> {
        self.finished.store(true, Ordering::Release);
        self.io_thread.unpark();
        join_io_thread(self.join.take());
        if let Some(e) = self.take_error() {
            return Err(e);
        }
        // The rich error may have already been drained by an earlier
        // `write`/`flush`; the sticky flag guarantees we still fail here.
        if self.errored.load(Ordering::Acquire) {
            return Err(io::Error::other("IO write thread reported an error"));
        }
        Ok(())
    }
}

impl Write for ThreadedWriter {
    fn write(&mut self, mut buf: &[u8]) -> io::Result<usize> {
        // Surface any pending IO-thread error once per call, not per
        // ring-push iteration. The previous per-iteration check acquired
        // the `Mutex` ~150 M times on a 30 GB run.
        if let Some(e) = self.take_error() {
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
                if let Some(e) = self.take_error() {
                    return Err(e);
                }
                if self.errored.load(Ordering::Acquire) {
                    return Err(io::Error::other("IO write thread reported an error"));
                }
            }
        }
        Ok(initial_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Nothing to flush at this layer — bytes are already in the ring
        // or written to `dst`. The IO thread does its own write_all.
        if let Some(e) = self.take_error() {
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
    errored: Arc<AtomicBool>,
    producer_thread: thread::Thread,
) {
    // Wake the producer on any exit (incl. panic) so it never parks against
    // a dead consumer; on panic also set the sticky flag + a fallback error.
    let _guard = PanicGuard {
        counterpart: &producer_thread,
        error: &error,
        errored: Some(&*errored),
        eof: None,
    };
    // Record an IO error: set the sticky flag *before* storing the rich
    // error so any thread that observes `errored` is guaranteed the error
    // slot is populated, then wake the producer (which may be parked on a
    // full ring) so it can surface the failure instead of blocking forever.
    let record_error = |e: io::Error| {
        errored.store(true, Ordering::Release);
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

    /// Once a `write` has surfaced (and drained) the rich error, `finish()`
    /// must still report failure via the sticky flag rather than falsely
    /// returning `Ok`. Regression test for `take_error` clearing the slot.
    #[test]
    fn threaded_writer_finish_fails_after_error_already_surfaced() {
        let ring = 4096;
        let payload = vec![7u8; 1024 * 1024];
        let mut w = ThreadedWriter::new(FailingSink, ring);
        // The oversized payload guarantees the producer parks and observes
        // the error, so the first `write_all` returns `Err` and drains it.
        assert!(w.write_all(&payload).is_err());
        assert!(w.finish().is_err(), "finish must stay failed after a drained error");
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
