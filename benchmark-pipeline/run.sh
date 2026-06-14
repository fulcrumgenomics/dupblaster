#!/usr/bin/env bash
# Wrapper around `pixi run snakemake` with the bench-resource pool sized so
# only one timed rule executes at a time even when --cores is large.
#
# Usage:
#   ./run.sh                       # full run, uses all cores for prep work
#   ./run.sh --dry-run             # preview the DAG
#   ./run.sh --cores 8             # cap cores
#   ./run.sh --config replicates=3 # bump from default 1 (or set samples=, dupblaster_bin=, picard_heap=)
# Anything after `--` is forwarded verbatim to snakemake.
set -euo pipefail
cd "$(dirname "$0")"

DRY=0
CORES="all"
EXTRA=()
while [ $# -gt 0 ]; do
  case "$1" in
    -n|--dry-run) DRY=1; shift ;;
    --cores)      CORES="$2"; shift 2 ;;
    --cores=*)    CORES="${1#--cores=}"; shift ;;
    --) shift; EXTRA+=("$@"); break ;;
    -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done

SNAKE_ARGS=(--cores "$CORES" --resources bench=100 --rerun-incomplete)
[ "$DRY" -eq 1 ] && SNAKE_ARGS+=(-n -p)
# Bash 3.2 (macOS default) errors under `set -u` on empty-array expansion;
# the `${name[@]+...}` form expands only if the array is set.
SNAKE_ARGS+=(${EXTRA[@]+"${EXTRA[@]}"})

echo "==> pixi run snakemake ${SNAKE_ARGS[*]}"
exec pixi run snakemake "${SNAKE_ARGS[@]}"
