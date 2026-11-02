#!/usr/bin/env bash
# Create a small, reproducible handoff archive from the LIVE server repository.
# Uses only Bash + Python 3 standard library; does not change/delete data or jobs.
set -euo pipefail

ROOT="$(cd "${1:-.}" && pwd -P)"
if [[ ! -f "$ROOT/run_pipeline.sh" || ! -f "$ROOT/wgs_test/wgs_pilot.py" ]]; then
  echo "ERROR: provide the live Minimal-SNP-Panels repository root (run_pipeline.sh and wgs_test/wgs_pilot.py required)." >&2
  exit 1
fi
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${2:-$(dirname "$ROOT")/Minimal-SNP-Panels-handoff-${STAMP}.zip}"
if [[ "$OUT" != /* ]]; then
  OUT="$(pwd -P)/$OUT"
fi
if [[ -e "$OUT" ]]; then
  echo "ERROR: archive already exists: $OUT" >&2
  exit 1
fi
mkdir -p "$(dirname "$OUT")"

python3 - "$ROOT" "$OUT" <<'PY'
from __future__ import annotations

import datetime as dt
import hashlib
import io
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

root, target = Path(sys.argv[1]), Path(sys.argv[2])
prefix = "Minimal-SNP-Panels/"

# Include the working tree (tracked AND untracked code), but never traverse FASTQ,
# BAM, GVCF, GenomicsDB, reference indices or the huge original results tree.
source_dirs = [
    "config", "data/metadata", "docs", "envs", "inputs",
    "literature", "src", "stages", "tests", "wgs_test",
    "results/pilot_orchestration/pilot_20261030",
    "results/beluga/pilot_20261030_SAMN51337163",
    "results/beluga/pilot_20261030_SAMN51337170",
    "results/polar_bear/pilot_20261030_SAMN16454148",
    "results/polar_bear/pilot_20261030_SAMN16454149",
]
ignored_dirs = {
    ".git", ".venv", "venv", "__pycache__", ".pytest_cache", ".mypy_cache",
    ".ruff_cache", ".tox", ".nox", ".ipynb_checkpoints", "node_modules",
    ".snakemake", ".next", "dist", "build", "genomicsdb", "tmp", "temp",
}
large_exts = (
    ".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam", ".bai", ".cram", ".crai",
    ".sam", ".vcf", ".vcf.gz", ".g.vcf.gz", ".bcf", ".csi", ".tbi",
    ".fna", ".fasta", ".fa", ".bwt", ".bwt.2bit.64", ".pac", ".ann", ".amb",
    ".sa", ".sra", ".sralite", ".zip", ".tar", ".tar.gz", ".tgz", ".7z",
    ".sqlite", ".db", ".part", ".incomplete.bam", ".jar", ".class", ".pyc",
)
secrets = re.compile(r"(^|/)(\.env($|\.)|\.ssh|\.aws|\.netrc|id_rsa|id_ed25519|credentials?($|\.)|secrets?($|\.))", re.I)
max_file_bytes = 10 * 1024 * 1024
max_total_bytes = 140 * 1024 * 1024
log_limit_bytes = 3 * 1024 * 1024
selected: dict[str, Path] = {}
skipped: list[str] = []

def eligible(p: Path) -> bool:
    if not p.is_file() or p.is_symlink():
        return False
    rel = p.relative_to(root).as_posix()
    rel_lower = rel.lower()
    if secrets.search(rel) or any(rel_lower.endswith(x) for x in large_exts):
        return False
    if p.stat().st_size > max_file_bytes:
        skipped.append(f"{rel}\tfile > 10 MiB")
        return False
    # Result directories: include only metadata and diagnostic text/logs, not heavy outputs.
    if rel.startswith("results/"):
        if rel.startswith("results/pilot_orchestration/pilot_20261030/"):
            return True
        if rel.endswith(("/run.json", "/versions.json", "/wgs/gvcf_manifest.tsv")):
            return True
        if "/logs/" in rel and rel.endswith((".log", ".out", ".err", ".txt")):
            return True
        return False
    return True

def collect_tree(d: Path) -> None:
    if not d.is_dir():
        return
    for current, dirs, filenames in os.walk(d, followlinks=False):
        dirs[:] = [name for name in sorted(dirs)
                   if name not in ignored_dirs
                   and not name.endswith(".egg-info")
                   and not (Path(current) / name).is_symlink()]
        for name in sorted(filenames):
            p = Path(current) / name
            if eligible(p):
                selected[p.relative_to(root).as_posix()] = p

for dirname in source_dirs:
    collect_tree(root / dirname)

for p in root.iterdir():
    if p.is_file() and eligible(p):
        selected[p.name] = p

if "run_pipeline.sh" not in selected or "wgs_test/wgs_pilot.py" not in selected:
    raise SystemExit("Key scripts not collected; refusing to create incomplete archive")


def command(args: list[str], *, limit: int = 8_000_000) -> bytes:
    try:
        completed = subprocess.run(
            args, cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=45, check=False,
        )
        content = completed.stdout
        return content[:limit] + (b"\n[TRUNCATED]\n" if len(content) > limit else b"")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"Not available: {exc}\n".encode()

meta = {
    "_handoff/README.txt": (
        "ARCHIVE PURPOSE: code/config/docs/tests + small diagnostic logs only.\n"
        "Made from the live project directory, including untracked code.\n"
        "FASTQ, BAM, gVCF, GenomicsDB, reference files, and old ZIPs deliberately excluded.\n"
        "_handoff/git-diff.patch may contain uncommitted local edits.\n"
        "Do not confuse this snapshot with the live server or start another pilot automatically.\n"
        f"Root on source server: {root}\n"
        f"Created: {dt.datetime.now().astimezone().isoformat()}\n"
    ).encode(),
    "_handoff/git-branch-and-status.txt": command(["git", "status", "--short", "--branch", "--untracked-files=no"]),
    "_handoff/git-last-commit.txt": command(["git", "log", "-1", "--format=fuller"]),
    "_handoff/git-diff.patch": command(["git", "diff", "--binary"], limit=12_000_000),
    "_handoff/git-staged-diff.patch": command(["git", "diff", "--cached", "--binary"], limit=12_000_000),
    "_handoff/slurm-active.txt": command(["squeue", "-u", os.environ.get("USER", "bkarpovich"),
                                          "-o", "%.20i %.25j %.12T %.15M %.25R"]),
    "_handoff/tool-versions.txt": b"",  # Filled below
}
versions = []
for tool, cmd in [
    ("python", ["python3", "--version"]),
    ("bwa-mem3", ["bwa-mem3", "version"]),
    ("samtools", ["samtools", "--version"]),
    ("bcftools", ["bcftools", "--version"]),
    ("gatk", ["gatk", "--version"]),
    ("slurm", ["squeue", "--version"]),
]:
    result = command(cmd, limit=2000).decode(errors="replace")
    versions.append(f"{tool}: {result[:400].strip()}\n")
meta["_handoff/tool-versions.txt"] = "\n".join(versions).encode()

rows = ["relative_path\toriginal_bytes\tsha256\tpacked_as"]
with zipfile.ZipFile(target, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6,
                     allowZip64=True) as archive:
    total_bytes = 0
    for rel, path in sorted(selected.items()):
        size = path.stat().st_size
        if total_bytes + min(size, log_limit_bytes) > max_total_bytes:
            skipped.append(f"{rel}\tarchive budget exhausted")
            continue
        is_log = rel.startswith("results/") and rel.endswith((".log", ".out", ".err"))
        if is_log and size > log_limit_bytes:
            # Preserve the diagnostic tail instead of an enormous progress-bar log.
            with path.open("rb") as handle:
                handle.seek(-min(size, 128 * 1024), 2)
                tail = handle.read()
            entry = rel + ".last128KiB.txt"
            archive.writestr(prefix + entry,
                             f"TRUNCATED LOG: {rel}; original bytes: {size}; last 128KiB follows\n".encode() + tail)
            rows.append(f"{rel}\t{size}\t{hashlib.sha256(tail).hexdigest()}\t{entry} (tail only)")
            total_bytes += len(tail)
        else:
            archive.write(path, prefix + rel)
            # Everything included is small (<=10MiB), safe to hash in chunks.
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            rows.append(f"{rel}\t{size}\t{digest.hexdigest()}\t{rel}")
            total_bytes += size

    archive.writestr(prefix + "_handoff/files_manifest.tsv", "\n".join(rows) + "\n")
    archive.writestr(prefix + "_handoff/skipped_files.tsv", "relative_path\treason\n" + "\n".join(skipped) + "\n")
    for name, data in meta.items():
        archive.writestr(prefix + name, data)

with zipfile.ZipFile(target) as check:
    bad = check.testzip()
    if bad is not None:
        raise SystemExit(f"ERROR: corrupt archive member: {bad}")

print(f"Created: {target}")
print(f"Files included: {len(rows) - 1}")
print(f"Uncompressed selected bytes: {total_bytes:,}")
print(f"Archive bytes: {target.stat().st_size:,}")
print("ZIP integrity: OK")
print("Note: inspect _handoff/git-diff.patch and configuration files for secrets before sharing.")
PY

sha256sum "$OUT"
echo "HANDOFF_ARCHIVE=$OUT"
