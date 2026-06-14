#!/usr/bin/env bash
# One-shot setup for the dupblaster benchmark pipeline.
#
# Installs pixi if it isn't already on PATH, materializes the pixi env, and
# builds dupblaster in release mode from the parent crate. Idempotent.
#
# Usage:
#   ./install.sh                # default
#   ./install.sh --skip-build   # set up env only; skip the cargo build
set -euo pipefail

cd "$(dirname "$0")"
PIPELINE_DIR="$PWD"
REPO_ROOT="$(cd .. && pwd)"

DO_BUILD=1
for arg in "$@"; do
  case "$arg" in
    --skip-build) DO_BUILD=0 ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

# ---- pixi -----------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
  log "Installing pixi to ~/.pixi"
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
else
  log "pixi found: $(pixi --version)"
fi

log "Materializing pixi env (snakemake, samtools, samblaster, GNU time, ...)"
pixi install

# ---- dupblaster (release) -------------------------------------------------
if [ "$DO_BUILD" -eq 1 ]; then
  if ! command -v cargo >/dev/null 2>&1; then
    echo "ERROR: cargo not on PATH. Install rustup (https://rustup.rs) and re-run." >&2
    exit 1
  fi
  # Build the whole workspace: the pipeline needs both `dupblaster` and the
  # `bench-compare` comparison tool (the `compare` rule invokes the latter).
  log "Building dupblaster + bench-compare (release) from $REPO_ROOT"
  ( cd "$REPO_ROOT" && cargo build --release --workspace )
  log "dupblaster:    $REPO_ROOT/target/release/dupblaster"
  log "bench-compare: $REPO_ROOT/target/release/bench-compare"
fi

cat <<EOF

Setup complete.

Next:
  ./run.sh --dry-run        # preview the job graph
  ./run.sh                  # run the pipeline (downloads ~17 GB on first run)

Override knobs at the CLI:
  ./run.sh --config replicates=3    # default is 1; bump to measure variance
  ./run.sh --config samples=HG03953
  ./run.sh --config dupblaster_bin=/path/to/dupblaster
EOF
