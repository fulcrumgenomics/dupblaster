# CLAUDE.md — dupblaster

Guidance for Claude (or any coding agent) working in this repo. Human-facing
contributor docs live in `CONTRIBUTING.md`; this file captures the things an
agent must not get wrong.

## What this is

`dupblaster` — fast, streaming duplicate marking for **query-grouped** SAM/BAM,
inspired by samblaster and Picard MarkDuplicates. Short-read only (not
ONT/PacBio).

## Before calling any change done

Run all four gates and make them pass:

```sh
cargo ci-fmt     # rustfmt --check
cargo ci-lint    # clippy --all-targets -D warnings
cargo ci-test    # nextest, --locked
cargo ci-doc     # rustdoc, -D warnings (catches broken intra-doc links)
```

The `ci-*` aliases live in `.cargo/config.toml`. If `ci-fmt` fails, run
`cargo fmt`. CI runs these four plus `cargo deny check`, for five gates in
total.

## Workflow — this repo is PUBLIC

- Work on a **branch + PR onto `main`**. Never commit, push, amend, or
  force-push `main` directly. The sole exception is the release tooling below,
  which the **human** runs.
- Keep PRs focused; each PR includes tests for the behavior it changes.

## CHANGELOG — update it IN THE SAME PR

For any **user-visible** change (behavior, CLI flags, output, metrics), add a
one-line entry to the `[Unreleased]` section of `CHANGELOG.md`, in the right
subsection (Added / Changed / Fixed / Removed), **as part of that same PR** —
not later during release prep. Internal-only changes (refactors, tests, CI,
release tooling) need no entry.

## Releases — NEVER publish

Releases run via `cargo-release` (`release.toml`). The publish is the **human's**
job, always:

- **Never** run `cargo publish` or `cargo release --execute`. Dry-run only
  (`cargo release X.Y.Z` with no `--execute`).
- Prep you *can* do: confirm the CHANGELOG `[Unreleased]` entry is in place, run
  the dry-run to verify the version bump / CHANGELOG rewrite / tag, then hand the
  human the `--execute` command.
- After the human publishes, create the matching **GitHub release from the tag**
  (`gh release create vX.Y.Z --notes-file ...`) using that version's CHANGELOG
  section as the notes.
- `release.toml` uses `{{version}}` (the version being released) in the
  pre-release commit message, and that only interpolates because
  `consolidate-commits = false` is set. In a workspace's *consolidated* commit
  context cargo-release binds only `{{date}}` — `{{version}}`,
  `{{next_version}}`, and `{{prev_version}}` all render as literals (how v0.1.1
  got the message "Bump version to {{next_version}}"). If the dry run prints
  `Unrendered {{...}} present in template`, the commit message is broken —
  don't hand the human the `--execute` command until it's clean.

## Testing conventions

- **Generate test data in code; never commit data files** (BAMs, etc.). The
  suite reads/writes SAM/BAM in-process via `noodles` — no `samtools` shell-out.
- Name tests after the behavior asserted; prefer many small tests over
  table-driven ones.
- See `tests/helpers/mod.rs` (`SamBuilder`, `run_and_capture`) for the
  integration-test pattern, and the in-module `#[cfg(test)]` blocks for unit
  tests of pure functions.

## More

See `CONTRIBUTING.md` for dependency policy, performance guidance, code style,
and bug-reporting detail.
