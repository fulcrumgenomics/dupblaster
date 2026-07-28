"""Walk results/runs/{sample}/{tool}/t{nslots}/rep{N}/time.devnull.txt, parse
each, emit one wide TSV row per (sample, tool, sweep point, rep).

Every run reports wall time and actual CPU (`cpu_s` = user + sys). Swept runs
additionally report `reserved_cpu_s` (= sweep point x wall): the CPU you had
to *claim* to run it, which is what a cloud instance or an SGE/Slurm
reservation actually costs. It is left empty for stock runs, which are
measured on a single-CPU assumption — there the claim would be 1 slot, so
reserved cost is just `wall_s`.

`--spec` carries each timed (tool, sweep point) and its class from the
Snakefile's TOOLS registry, so those facts live in one place.
"""

import argparse
import sys
from pathlib import Path

# Snakemake puts the scripts dir on sys.path via -m, but we may also be run
# as a plain script — add our own dir explicitly to be safe.
sys.path.insert(0, str(Path(__file__).parent))
from parse_gnu_time import parse  # noqa: E402


COLUMNS = [
    "sample", "tool", "tool_class", "nslots", "rep",
    "wall_s", "cpu_s", "reserved_cpu_s",
    "user_s", "sys_s", "cpu_percent", "max_rss_kb", "exit_status",
]


def parse_spec(entries: list[str]) -> dict[tuple[str, str], str]:
    """`tool:nslots:class` -> {(tool, nslots): class}."""
    spec = {}
    for entry in entries:
        tool, nslots, cls = entry.split(":")
        spec[(tool, nslots)] = cls
    return spec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="Root of per-run outputs, e.g. results/runs")
    ap.add_argument("--output", required=True)
    ap.add_argument("--spec", nargs="*", default=[],
                    help="tool:nslots:class entries, from the Snakefile's "
                         "TOOLS registry")
    args = ap.parse_args()

    spec = parse_spec(args.spec)
    root = Path(args.root)
    rows: list[dict] = []
    for time_txt in sorted(root.glob("*/*/t*/rep*/time.devnull.txt")):
        # path = root/sample/tool/t{nslots}/repN/time.devnull.txt
        rep_dir, slots_dir, tool_dir, sample_dir = time_txt.parents[0:4]
        tool = tool_dir.name
        nslots = slots_dir.name.removeprefix("t")
        m = parse(time_txt)
        wall, user, sys_s = (m.get(k, "") for k in ("wall_s", "user_s", "sys_s"))
        cpu = round(user + sys_s, 2) if user != "" and sys_s != "" else ""
        # Reserved cost only means something for a swept run, where the slots
        # claimed are an explicit input; stock runs assume a single CPU.
        reserved_cpu = round(int(nslots) * wall, 2) if nslots.isdigit() and wall != "" else ""
        rows.append({
            "sample": sample_dir.name,
            "tool": tool,
            "tool_class": spec.get((tool, nslots), ""),
            "nslots": nslots,
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
