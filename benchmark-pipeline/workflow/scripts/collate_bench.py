"""Walk results/runs/{sample}/{tool}/rep{N}/time.txt, parse each, emit one
wide TSV row per (sample, tool, rep)."""

import argparse
import sys
from pathlib import Path

# Snakemake puts the scripts dir on sys.path via -m, but we may also be run
# as a plain script — add our own dir explicitly to be safe.
sys.path.insert(0, str(Path(__file__).parent))
from parse_gnu_time import parse  # noqa: E402


COLUMNS = [
    "sample", "tool", "rep",
    "wall_s", "user_s", "sys_s", "cpu_percent", "max_rss_kb", "exit_status",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="Root of per-run outputs, e.g. results/runs")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    root = Path(args.root)
    rows: list[dict] = []
    for time_txt in sorted(root.glob("*/*/rep*/time.devnull.txt")):
        # path = root/sample/tool/repN/time.devnull.txt
        rep_dir, tool_dir, sample_dir = time_txt.parents[0:3]
        sample = sample_dir.name
        tool = tool_dir.name
        rep = int(rep_dir.name.removeprefix("rep"))
        m = parse(time_txt)
        rows.append({
            "sample": sample, "tool": tool, "rep": rep,
            "wall_s":      m.get("wall_s", ""),
            "user_s":      m.get("user_s", ""),
            "sys_s":       m.get("sys_s", ""),
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
