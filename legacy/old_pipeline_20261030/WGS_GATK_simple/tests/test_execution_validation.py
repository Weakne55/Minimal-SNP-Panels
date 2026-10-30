"""Safe state/outputs regression tests; no WGS subprocesses are launched."""

import json
from pathlib import Path

import pytest

import pipeline


def task_for(path, stage="haplotypecaller"):
    return {
        "id": "test", "stage": stage, "deps": [], "inputs": [],
        "outputs": [str(path)], "signature": "signature", "body": "unused",
    }


def test_vcf_index_is_required_for_completion(tmp_path):
    workdir = tmp_path / "work"
    marker = workdir / "checkpoints" / "test.json"
    marker.parent.mkdir(parents=True)
    vcf = tmp_path / "calls.vcf.gz"
    vcf.write_bytes(b"nonempty")
    task = task_for(vcf)
    marker.write_text(json.dumps({"signature": "signature"}))

    assert pipeline.missing_output(task) == str(vcf) + " (index)"
    assert not pipeline.is_completed(task, workdir)

    csi = Path(str(vcf) + ".csi")
    csi.write_bytes(b"index")
    assert pipeline.missing_output(task) is None
    assert pipeline.is_completed(task, workdir)

    vcf.write_bytes(b"")
    assert pipeline.missing_output(task) == str(vcf)
    assert not pipeline.is_completed(task, workdir)


def test_non_gatk_stage_does_not_need_vcf_index(tmp_path):
    output = tmp_path / "report.html"
    task = task_for(output, stage="multiqc")
    assert pipeline.missing_output(task) == str(output)
    output.write_text("ready")
    assert pipeline.missing_output(task) is None


def test_failed_task_has_no_checkpoint(tmp_path):
    workdir = tmp_path / "work"
    vcf = workdir / "call.vcf.gz"
    task = task_for(vcf)
    task["body"] = "set -euo pipefail\nprintf 'data' > " + str(vcf) + "\n"
    plan = pipeline.write_plan([task], {}, tmp_path / "config.yaml", workdir)
    with pytest.raises(RuntimeError, match="missing VCF index"):
        pipeline.run_one_task(plan, "test")
    assert not (workdir / "checkpoints" / "test.json").exists()
