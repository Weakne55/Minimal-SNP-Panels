#!/usr/bin/env python3
"""Convenient dry-run-first launcher for workflow/Snakefile.

Only orchestrator dependencies (Snakemake, PyYAML) are required on the host.
All bioinformatics tools are executed via the configured Apptainer/Singularity image.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from minsnp.core import identifier


def command(args: argparse.Namespace) -> list[str]:
    cmd = ["snakemake", "--snakefile", "workflow/Snakefile", "--configfile", "config/snakemake.yaml",
           "--profile", f"profiles/{args.backend}"]
    overrides: list[str] = []
    if args.run_id:
        overrides.append("run_id=" + identifier(args.run_id))
    if args.fastq:
        manifests = [str(Path(item).expanduser().resolve()) for item in args.fastq]
        overrides.append("manifests=" + json.dumps(manifests, separators=(",", ":")))
    if args.image:
        overrides.append("container_image=" + str(Path(args.image).expanduser().resolve())
                         if not args.image.startswith(("docker://", "oras://", "library://"))
                         else "container_image=" + args.image)
    if args.trim:
        overrides.append("trim_fastq=true")
    if args.qc:
        overrides.append("variant_qc=true")
    if overrides:
        cmd += ["--config", *overrides]
    if args.jobs is not None:
        cmd += ["--jobs", str(args.jobs)]
    if args.cores is not None:
        cmd += ["--cores", str(args.cores)]
    if args.bind:
        binds = ",".join(f"{Path(d).expanduser().resolve()}:{Path(d).expanduser().resolve()}" for d in args.bind)
        cmd += ["--apptainer-args", f"--bind {binds}"]
    if not args.submit:
        cmd += ["--dry-run"]
    cmd += ["--show-failed-logs"]
    return cmd


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Snakemake + Apptainer WGS workflow (dry-run by default)")
    parser.add_argument("--fastq", nargs="+", metavar="TSV", help="FASTQ manifests from one or both species")
    parser.add_argument("--run-id", help="Output run ID, default from config/snakemake.yaml")
    parser.add_argument("--backend", choices=["local", "slurm"], default="slurm")
    parser.add_argument("--image", help="Absolute pre-built SIF path or docker:// OCI reference")
    parser.add_argument("--bind", action="append", metavar="DIR", help="External host root, e.g. /mnt/tank (can repeat)")
    parser.add_argument("--trim", action="store_true", help="Enable fastp trimming")
    parser.add_argument("--qc", action="store_true", help="Run technical variant QC after joint genotyping")
    parser.add_argument("--jobs", type=int, help="Max concurrent SLURM jobs")
    parser.add_argument("--cores", type=int, help="Local core limit")
    parser.add_argument("--submit", action="store_true", help="Actually execute (default only dry-run)")
    args = parser.parse_args(argv)
    cmd = command(args)
    print("$ " + shlex.join(cmd), flush=True)
    try:
        return subprocess.run(cmd, cwd=ROOT, check=False).returncode
    except FileNotFoundError:
        print("Snakemake not found. Install envs/snakemake.yaml first.", file=sys.stderr)
        return 127


if __name__ == "__main__":
    sys.exit(main())
