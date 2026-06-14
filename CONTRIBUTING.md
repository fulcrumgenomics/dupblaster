# Contributing to dupblaster

Thanks for your interest in dupblaster. This document covers the dev
loop, code-style expectations, the release flow, and conventions for
contributors and maintainers.

dupblaster operates under Fulcrum Genomics' organisation-level Code of
Conduct, which applies to all interactions in this repository.

## Getting Started

**Prerequisites:**
- Rust stable, minimum version from `rust-toolchain.toml` /
  `Cargo.toml`'s `rust-version` field (currently **1.89**).
- [`cargo-nextest`][nextest] for the test runner used in CI.
- [`cargo-deny`][cargo-deny] for the supply-chain check (optional
  locally; CI runs it on every PR).

The test suite needs no external tools — BAM input is built and read
back in-process via `noodles` (no `samtools` shell-out).

```sh
cargo build              # debug build
cargo build --release    # release build (portable, runs anywhere)
```

`cargo build --release` deliberately does **not** set
`target-cpu=native`, so binaries built locally and binaries built by
`cargo install dupblaster` from crates.io produce the same artifact.
For locally-tuned benchmarking, opt in explicitly:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

[nextest]: https://nexte.st/
[cargo-deny]: https://embarkstudios.github.io/cargo-deny/
[samtools]: http://www.htslib.org/

## Verification Checklist

Run all four before sending a PR. CI runs the same gates.

```sh
cargo ci-fmt    # rustfmt --check
cargo ci-lint   # clippy --all-targets -D warnings
cargo ci-test   # nextest, --locked
cargo deny check  # licenses, advisories, bans, sources
```

The `ci-*` aliases live in `.cargo/config.toml`. If `cargo ci-fmt`
fails, run `cargo fmt` and re-stage.

## Code Style

dupblaster follows the [Rust API Guidelines][rust-api] and a few
project-local rules:

- **Idiomatic Rust.** Don't transliterate from C or Python; write Rust.
- **Names matter.** Prefer meaningful names even if longer
  (`signature_table_cells`, not `cells`). Short names are fine in
  closures and tight loops.
- **Small, focused functions.** Extract helpers when a function's
  responsibilities start to fan out. Aim for code that makes sense
  when you come back to it in six months.
- **Doc comments on every public item.** Private items get doc
  comments when behavior is non-obvious. Comments should explain
  *why*, not *what* — let the code show the what.
- **No premature abstraction.** Solve the problem in front of you.
  Three similar lines are fine; refactor only when a real third use
  arrives.

[rust-api]: https://rust-lang.github.io/api-guidelines/

## Testing

dupblaster has three test surfaces:

1. **In-module unit tests** (`#[cfg(test)] mod tests { … }`) for pure
   functions: parsers, the dup signature, metric calculations, CIGAR
   logic, etc.
2. **Integration tests** in `tests/` for end-to-end runs of the
   binary against built-in or programmatically generated SAM/BAM
   inputs.
3. **Doctests** on documented public items where examples clarify
   usage. Keep these short; they're slower than unit tests.

### Test data is generated in code, never committed

A `SamBuilder` helper in `tests/helpers/mod.rs` builds SAM inputs
programmatically. Use it; don't commit `.sam` or `.bam` fixture files.
Reviewers can see the exact input next to the assertion:

```rust
SamBuilder::new()
    .sq("chr1", 1_000_000)
    .rec_simple("r1", 99, "chr1", 100, "50M", "=", 200, 150)
    .rec_simple("r1", 147, "chr1", 200, "50M", "=", 100, -150)
    .write_to(&env.input);
```

The same applies to BAM input: build a SAM with `SamBuilder` and
convert it to BAM in-process with the `helpers::sam_to_bam` noodles
helper (see `tests/test_bam.rs`) — no external `samtools` needed.

### Test naming and structure

- Each integration-test file should focus on one feature surface
  (e.g. `test_stats.rs`, `test_compression.rs`,
  `test_query_grouped_check.rs`).
- Name tests after the behavior they're asserting:
  `level_6_produces_smaller_output_than_default` is better than
  `test_compression_2`.
- Prefer many small tests over parameterized / table-driven ones —
  small tests are easier to debug when one fails.
- Cover the happy path, error cases, and edge cases.

## Performance

- **Correct first, fast second.** Every optimization needs a baseline
  test it doesn't regress.
- **Profile before tuning.** Use [`samply`][samply] or
  cargo-flamegraph to find real hot paths instead of guessing.
- **dupblaster's hot path** is the worker thread that pulls records
  from the read ring buffer, computes signatures, and pushes results
  into the write ring buffer. Most allocations are pooled; new
  allocations in the hot path are a regression.
- **Benchmark dataset and harness** live in `benchmark-pipeline/`.
  Run with `cd benchmark-pipeline && ./install.sh && ./run.sh`.
  Reproduces wall time, RSS, set-equivalence concordance vs Picard,
  orphan-discordance triage, and supplementary-flag inheritance across
  the dup-marking tools in the suite (dupblaster's modes, samblaster,
  Picard MarkDuplicates, samtools markdup, dupsifter) on a subsampled
  NYGC 1000G HG03953 CRAM; see `benchmark-pipeline/README.md` for the
  authoritative tool list. See
  `benchmark-pipeline/README.md` (including "Adding a new tool" if
  you want to plug another marker into the comparison); high-level
  numbers are in the top-level `README.md` Benchmarks section.

[samply]: https://github.com/mstange/samply

## Adding or upgrading dependencies

- Prefer dependencies that are already in the tree. If we already
  pull in a crate (`anyhow`, `clap`, `noodles-sam`, `bgzf`, …),
  reach for it before adding a new one.
- New direct deps need a clear justification in the PR.
- Pin to the latest stable major/minor at time of add. Dependabot
  will keep things current after that.
- After any dep change, run `cargo deny check`. New licenses get
  added to `deny.toml`'s allow-list only after deliberate review
  (no copyleft).

## Reporting Bugs

Open a GitHub issue with:

- The dupblaster version (`dupblaster --version`).
- The input format (SAM or BAM), the rough size, and any non-default
  flags used.
- The full stderr from the failing run.
- A **minimal** reproducer if at all possible — a small `samtools
  view` selection that still triggers the bug is ideal. Don't share
  proprietary data; we don't need it to debug.

## Pull Requests

- Keep PRs focused. 250–1000 LOC is a comfortable review size; bigger
  changes should be split into a stack of small commits or staged PRs.
- Commit messages explain *why*. The "what" is in the diff.
- Each PR should include tests for the behavior it adds or changes.
- Update CHANGELOG.md's `[Unreleased]` section with a one-line entry
  in the appropriate subsection (Added / Changed / Fixed / Removed).
- All four CI gates must be green before merge.

## Releasing

Releases are cut with [cargo-release]. Install it once with
`cargo install cargo-release`, and run `cargo login` once with a
crates.io API token so the publish step works. Configuration lives in
`release.toml` at the repo root.

```sh
# Dry run — review what would change. CHANGELOG.md's [Unreleased]
# section gets renamed to the new version + date, a fresh empty
# [Unreleased] section is inserted, and link refs are rewritten.
cargo release 0.1.0

# Real release — bumps Cargo.toml, updates Cargo.lock, commits with
# message "Bump version to 0.1.0", tags v0.1.0, pushes to origin,
# and publishes to crates.io.
cargo release 0.1.0 --execute
```

After the push, create the GitHub release object pointing at the new
tag (e.g. `gh release create v0.1.0 --generate-notes`), then update
the [bioconda dupblaster recipe][bioconda-recipe] with the new
version + tarball SHA256:

```sh
curl -sL https://github.com/fulcrumgenomics/dupblaster/archive/refs/tags/v0.1.0.tar.gz \
    | shasum -a 256
```

[cargo-release]: https://github.com/crate-ci/cargo-release
[bioconda-recipe]: https://github.com/bioconda/bioconda-recipes/tree/master/recipes/dupblaster
