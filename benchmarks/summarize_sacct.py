#!/usr/bin/env python3
"""Read existing sacct -P data without submitting jobs or touching WGS files."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


FIELDS = ["JobID", "JobName", "State", "ElapsedRaw", "AllocCPUS", "TotalCPU", "MaxRSS"]


def parse_duration(s: str) -> float:
    """Parse Slurm CPU durations as seconds (D-HH:MM:SS or HH:MM:SS)."""
    s = s.strip()
    if not s or s in {"Unknown", "N/A"}:
        return 0.0
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    parts = s.split(":")
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h, (m, sec) = "0", parts
    elif len(parts) == 1:
        h, m, sec = "0", "0", parts[0]
    else:
        raise ValueError(f"Invalid CPU duration: {s!r}")
    return days * 86400 + int(h) * 3600 + int(m) * 60 + float(sec)


def parse_memory_gib(value: str) -> float | None:
    """Interpret Slurm K/M/G/T suffixes as IEC units; ignore unknown/missing."""
    v = value.strip()
    if not v or v in {"0", "Unknown", "N/A"}:
        return None
    suffix = v[-1].upper()
    if suffix in "KMGTP":
        x = float(v[:-1])
        return x * (1024 ** ("KMGTP".index(suffix) - 2))
    return float(v) / (1024 ** 3)


def _read_rows(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="|")
        first = next(reader, None)
        if first is None:
            return rows
        columns = first if first[:2] == FIELDS[:2] else FIELDS
        if first != columns:
            rows.append(dict(zip(columns, first)))
        for parts in reader:
            if parts and any(field.strip() for field in parts):
                rows.append(dict(zip(columns, parts)))
    return rows


def summarize(path: Path) -> list[dict[str, object]]:
    rows = _read_rows(path)
    steps_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    mains: list[dict[str, str]] = []
    for row in rows:
        jobid = row.get("JobID", "").strip()
        if "." in jobid:
            steps_by_id[jobid.split(".", 1)[0]].append(row)
        elif jobid:
            mains.append(row)
    by_name: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in mains:
        name = row.get("JobName", "") or "(unnamed)"
        try:
            wall = float(row.get("ElapsedRaw", "0") or 0)
            cpus = int(float(row.get("AllocCPUS", "0") or 0))
            used_raw = row.get("TotalCPU", "").strip()
            used = parse_duration(used_raw) if used_raw and used_raw not in {"Unknown", "N/A"} else None
        except ValueError:
            continue
        steps = steps_by_id.get(row.get("JobID", ""), [])
        rss = [parse_memory_gib(r.get("MaxRSS", "")) for r in [row, *steps]]
        by_name[name].append({
            "wall_s": wall,
            "cpu_s": used,
            "cpu_efficiency": used / (wall * cpus) if used is not None and wall > 0 and cpus else None,
            "rss_gib": max((x for x in rss if x is not None), default=None),
            "state": row.get("State", ""),
        })
    output = []
    for name, jobs in by_name.items():
        efficiencies = [j["cpu_efficiency"] for j in jobs if j["cpu_efficiency"] is not None]
        rss_values = [j["rss_gib"] for j in jobs if j["rss_gib"] is not None]
        output.append({
            "job_name": name,
            "jobs": len(jobs),
            "total_wall_h": sum(j["wall_s"] for j in jobs) / 3600,
            "total_used_cpu_h": (sum(j["cpu_s"] for j in jobs) / 3600
                                 if all(j["cpu_s"] is not None for j in jobs) else None),
            "median_cpu_eff_pct": 100 * statistics.median(efficiencies) if efficiencies else None,
            "max_reported_rss_gib": max(rss_values, default=None),
            "non_completed": sum(not str(j["state"]).startswith("COMPLETED") for j in jobs),
        })
    return sorted(output, key=lambda x: -x["total_wall_h"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sacct_tsv", type=Path, help="sacct -P [-n] with seven documented fields")
    args = parser.parse_args(argv)
    rows = summarize(args.sacct_tsv)
    print("job_name\tjobs\tΣwall_h\tΣCPU_used_h\tmedian_cpu_eff_%\tmax_reported_RSS_GiB\tnon_completed")
    for r in rows:
        fmt = lambda x: "NA" if x is None else f"{x:.2f}"
        print("\t".join([
            str(r["job_name"]), str(r["jobs"]), fmt(r["total_wall_h"]),
            fmt(r["total_used_cpu_h"]), fmt(r["median_cpu_eff_pct"]),
            fmt(r["max_reported_rss_gib"]), str(r["non_completed"]),
        ]))
    print("# NOTE: Σwall_h is summed job durations, not pipeline walltime; "
          "RSS may be absent or step-specific; CPU efficiency is indicative.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
