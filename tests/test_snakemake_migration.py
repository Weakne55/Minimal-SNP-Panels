"""Snakemake migration contracts, independent of a cluster and large FASTQ data."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def launcher():
    spec = importlib.util.spec_from_file_location("minsnp_snakemake_launcher", ROOT / "scripts/run_snakemake.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def args(**updates):
    values = dict(run_id="test_20261102", fastq=["inputs/fastq.template.tsv"],
                  backend="slurm", image=None, bind=["/mnt/tank"], trim=False, qc=False,
                  jobs=3, cores=None, submit=False)
    values.update(updates)
    return argparse.Namespace(**values)


def test_launcher_dry_run_and_manifests():
    cmd = launcher().command(args())
    assert cmd[0] == "snakemake"
    assert "--dry-run" in cmd
    assert "--profile" in cmd
    assert "profiles/slurm" in cmd
    assert "--apptainer-args" in cmd
    assert "--bind /mnt/tank:/mnt/tank" in cmd
    assert "--jobs" in cmd
    value = next(item for item in cmd if item.startswith("manifests="))
    parsed = json.loads(value.partition("=")[2])
    assert parsed == [str((Path.cwd() / "inputs/fastq.template.tsv").resolve())]


def test_launcher_submit_and_optional_settings():
    cmd = launcher().command(args(submit=True, trim=True, qc=True, backend="local", bind=[]))
    assert "--dry-run" not in cmd
    assert "profiles/local" in cmd
    assert "trim_fastq=true" in cmd
    assert "variant_qc=true" in cmd


def test_configs_and_rules_present():
    workflow = (ROOT / "workflow/Snakefile").read_text()
    cfg = yaml.safe_load((ROOT / "config/snakemake.yaml").read_text())
    for rule in ("preflight", "trim", "align", "merge_bam", "markdup", "bam_qc",
                 "haplotypecaller", "sample_map", "genomicsdb", "genotype", "concat",
                 "export_vcf", "validate_cohort", "software_versions", "variant_qc", "variant_qc_stats", "all"):
        assert f"rule {rule}:" in workflow
    assert "container: IMAGE" in workflow
    assert "cohort.lock.json" in workflow
    assert cfg["variant_qc"] is False
    assert cfg["container_image"].endswith(".sif")
    for stage, resources in cfg["stage_resources"].items():
        assert min(resources["threads"], resources["mem_mb"], resources["runtime"]) >= 1, stage
    profile = yaml.safe_load((ROOT / "profiles/slurm/config.yaml").read_text())
    assert profile["executor"] == "slurm"
    assert profile["use-apptainer"] is True
    assert "--nodelist" in profile["default-resources"]["slurm_extra"]


def test_container_has_explicit_path_for_singularity_exec():
    dockerfile = (ROOT / "containers/Dockerfile.snakemake").read_text()
    assert "ENV PATH=/opt/conda/envs/wgs-gatk-mem3/bin:$PATH" in dockerfile
    assert "envs/wgs-gatk-mem3.yaml" in dockerfile
    assert "ENTRYPOINT" not in "\n".join(line for line in dockerfile.splitlines() if not line.startswith("#"))
