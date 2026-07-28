"""Walk results/runs/{sample}/{tool}/t{nslots}/rep{N}/time.devnull.txt, parse
each, emit one wide TSV row per (sample, tool, sweep point, rep).

Two cost metrics are reported per run:

  cpu_s           actual CPU consumed (user + sys) — the work the tool really
                  did, independent of how it spread that work over threads.
  reserved_cpu_s  slots x wall — the CPU you had to *claim* to run it, which
                  is what a cloud instance or an SGE/Slurm reservation
                  actually costs. A tool whose threads hide IO latency rather
                  than adding throughput pays here without gaining: dupblaster
                  is charged for 3 slots while drawing ~1.1 cores of work.

`--spec` carries the class and reserved-slot count for each (tool, sweep
point) from the Snakefile's TOOLS registry, so those facts live in one place.
"""

import argparse
import sys
from pathlib import Path

# Snakemake puts the scripts dir on sys.path via -m, but we may also be run
# as a plain script — add our own dir explicitly to be safe.
sys.path.insert(0, str(Path(__file__).parent))
from parse_gnu_time import parse  # noqa: E402


COLUMNS = [
    "sample", "tool", "tool_class", "nslots", "reserved_slots", "rep",
    "wall_s", "cpu_s", "reserved_cpu_s",
    "user_s", "sys_s", "cpu_percent", "max_rss_kb", "exit_status",
]


def parse_spec(entries: list[str]) -> dict[tuple[str, str], tuple[str, int]]:
    """`tool:nslots:class:reserved_slots` -> {(tool, nslots): (class, slots)}."""
    spec = {}
    for entry in entries:
        tool, nslots, cls, slots = entry.split(":")
        spec[(tool, nslots)] = (cls, int(slots))
    return spec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="Root of per-run outputs, e.g. results/runs")
    ap.add_argument("--output", required=True)
    ap.add_argument("--spec", nargs="*", default=[],
                    help="tool:nslots:class:reserved_slots entries, from the "
                         "Snakefile's TOOLS registry")
    args = ap.parse_args()

    spec = parse_spec(args.spec)
    root = Path(args.root)
    rows: list[dict] = []
    for time_txt in sorted(root.glob("*/*/t*/rep*/time.devnull.txt")):
        # path = root/sample/tool/t{nslots}/repN/time.devnull.txt
        rep_dir, slots_dir, tool_dir, sample_dir = time_txt.parents[0:4]
        tool = tool_dir.name
        nslots = slots_dir.name.removeprefix("t")
        cls, reserved = spec.get((tool, nslots), ("", 0))
        m = parse(time_txt)
        wall, user, sys_s = (m.get(k, "") for k in ("wall_s", "user_s", "sys_s"))
        cpu = round(user + sys_s, 2) if user != "" and sys_s != "" else ""
        reserved_cpu = round(reserved * wall, 2) if wall != "" and reserved else ""
        rows.append({
            "sample": sample_dir.name,
            "tool": tool,
            "tool_class": cls,
            "nslots": nslots,
            "reserved_slots": reserved,
            "rep": int(rep_dir.name.removeprefix("rep")),
            "wall_s":         wall,
            "cpu_s":          cpu,
            "reserved_cpu_s": reserved_cpu,
            "user_s":      user,
            "sys_s":       sys_s,
            "cpu_percent": m.get("cpu_percent", ""),
            "max_rss_kb":  m.get("max_rss_kb", ""),
            "exit_status": m.get("exit_status", ""),
        })

    with open(args.output, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


if __name__ == "__main__":
    main()
