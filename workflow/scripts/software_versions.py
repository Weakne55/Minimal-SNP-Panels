#!/usr/bin/env python3
"""Record actual executable versions from inside the WGS container."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


def main(dest: str):
    commands = {"bwa-mem3": ["bwa-mem3", "version"], "samtools": ["samtools", "--version"],
                "bcftools": ["bcftools", "--version"], "gatk": ["gatk", "--version"],
                "fastp": ["fastp", "--version"], "mosdepth": ["mosdepth", "--version"],
                "plink2": ["plink2", "--version"]}
    data = {"python": sys.version.split()[0]}
    for name, cmd in commands.items():
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            lines = (result.stdout + "\n" + result.stderr).strip().splitlines()
            data[name] = {"exit_code": result.returncode, "version": lines[0][:250] if lines else "<empty>"}
        except (OSError, subprocess.TimeoutExpired) as exc:
            data[name] = {"error": str(exc)}
    target = Path(dest)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: software_versions.py OUTPUT.json")
    main(sys.argv[1])
