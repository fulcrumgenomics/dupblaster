//! dupblaster closes its standard output as soon as the last record is written,
//! before it writes the metrics files, so that a reader downstream in a pipe
//! (e.g. a sort) sees end of file without waiting for the rest of the run.
#![cfg(unix)]

mod helpers;

use std::ffi::CString;
use std::io::Read;
use std::os::unix::ffi::OsStrExt;
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use helpers::*;

/// How long any one step may block before the test fails instead of hanging.
const TIMEOUT: Duration = Duration::from_secs(30);

/// Kills the child if the test fails while it is still blocked, so a failing
/// test doesn't leave a dupblaster process behind.
struct KillOnDrop(Child);

impl Drop for KillOnDrop {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// Run `read` on its own thread and return what it returns, failing with `what`
/// if it takes longer than [`TIMEOUT`]. A thread still blocked at that point is
/// abandoned; it goes away with the test process.
fn within_timeout<T: Send + 'static>(
    what: &str,
    read: impl FnOnce() -> std::io::Result<T> + Send + 'static,
) -> T {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(read());
    });
    rx.recv_timeout(TIMEOUT).unwrap_or_else(|_| panic!("{what}")).unwrap()
}

/// The duplicate-metrics file is a FIFO, so dupblaster blocks when it opens it
/// until this test opens the other end. Reading stdout to end of file first can
/// only finish if dupblaster closed stdout before writing its metrics.
///
/// The sequencing-duplicate split runs here (it's on by default), with read names
/// it can parse. On an input this small it finishes in no time, so the test can't
/// tell whether the split ran before or after stdout closed; it checks the default
/// configuration rather than that ordering.
#[test]
fn stdout_reaches_end_of_file_before_the_metrics_are_written() {
    let env = TestEnv::new();
    let name = "M1:1:FC1:1:1101:1000:2000";
    SamBuilder::new()
        .sq("chr1", 1_000_000)
        .rec_simple(name, 99, "chr1", 100, "50M", "=", 200, 150)
        .rec_simple(name, 147, "chr1", 200, "50M", "=", 100, -150)
        .write_to(&env.input);
    let prefix = metrics_prefix_in(env._tmp.path());
    let fifo = prefix.with_extension("duplicate-metrics.tsv");
    let fifo_c = CString::new(fifo.as_os_str().as_bytes()).unwrap();
    // SAFETY: `fifo_c` is a valid NUL-terminated path for the duration of the call.
    assert_eq!(unsafe { libc::mkfifo(fifo_c.as_ptr(), 0o600) }, 0, "mkfifo {}", fifo.display());

    let mut child = KillOnDrop(
        Command::new(rust_binary())
            .args(["-i"])
            .arg(&env.input)
            .args(["--metrics-prefix"])
            .arg(&prefix)
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .unwrap(),
    );

    let mut stdout = child.0.stdout.take().unwrap();
    let bytes = within_timeout(
        "stdout stayed open while dupblaster waited to write its metrics",
        move || {
            let mut bytes = Vec::new();
            stdout.read_to_end(&mut bytes).map(|_| bytes)
        },
    );
    // A complete BGZF stream ends with the 28-byte empty block that marks end of file.
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    assert!(bytes.ends_with(&BGZF_EOF), "stdout ended before the BGZF end-of-file marker");

    // Let dupblaster write its metrics and finish.
    let metrics = within_timeout("dupblaster never wrote its duplicate metrics", move || {
        let mut metrics = String::new();
        std::fs::File::open(&fifo)?.read_to_string(&mut metrics).map(|_| metrics)
    });
    assert!(metrics.lines().count() > 1, "no metrics rows: {metrics:?}");
    let deadline = Instant::now() + TIMEOUT;
    let status = loop {
        if let Some(status) = child.0.try_wait().unwrap() {
            break status;
        }
        assert!(Instant::now() < deadline, "dupblaster didn't exit after writing its metrics");
        std::thread::sleep(Duration::from_millis(20));
    };
    assert!(status.success());
}
