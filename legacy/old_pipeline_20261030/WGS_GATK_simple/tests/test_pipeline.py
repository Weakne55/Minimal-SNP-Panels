"""Unit tests that do NOT require BWA, GATK or a SLURM cluster.

Run from the repository root: python3 -m pytest -q
"""
import json
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import model
import pipeline
import workflow


def make_example(tmp_path, mode="fastq", trim=False):
    """Build a tiny fake reference and sample table for planning tests."""
    ref = tmp_path / "species.fa"
    ref.write_text(">chr1\nACGT\n>chr2\nACGT\n")
    Path(str(ref) + ".fai").write_text("chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n")
    (tmp_path / "species.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n@SQ\tSN:chr2\tLN:4\n"
    )
    for ext in ("amb", "ann", "bwt.2bit.64", "pac"):
        Path(str(ref) + "." + ext).write_text("index")

    entries = []
    for name in ("S1", "S2"):
        data = dict(sample=name, species="Test_species", input_type=mode,
                    library="L1", lane="lane1", r1="", r2="", bam="", cram="")
        if mode == "fastq":
            for mate in ("r1", "r2"):
                f = tmp_path / f"{name}_{mate}.fq.gz"
                f.write_text("test")
                data[mate] = str(f)
        else:
            file = tmp_path / f"{name}.{mode}"
            file.write_text("test")
            Path(str(file) + (".bai" if mode == "bam" else ".crai")).write_text("index")
            data[mode] = str(file)
        entries.append(data)

    tsv = tmp_path / "samples.tsv"
    with tsv.open("w") as out:
        out.write("\t".join(model.COLUMNS) + "\n")
        for entry in entries:
            out.write("\t".join(entry[key] for key in model.COLUMNS) + "\n")

    config = tmp_path / "config.yaml"
    cfg = dict(
        samples=str(tsv), workdir=str(tmp_path / "output"),
        references={"Test_species": {"fasta": str(ref), "contigs": ["chr1", "chr2"]}},
        preprocessing={"trim_fastq": trim},
        duplicate_marking={"optical_distance": 0},
        gatk=dict(hc_threads=4, hc_heap_gb=16, db_heap_gb=24,
                  genotype_heap_gb=12, db_batch_size=50, db_readers=1),
        slurm={"max_array_parallel": 2},
    )
    config.write_text(yaml.safe_dump(cfg))
    return config, cfg


@pytest.mark.parametrize("mode", ["fastq", "bam", "cram"])
def test_build_expected_stages(tmp_path, mode):
    path, _ = make_example(tmp_path, mode=mode)
    _, cfg, workdir = workflow.config_load(str(path))
    units, samples, groups = workflow.validate(cfg)
    tasks, maps = workflow.make_tasks(cfg, workdir, units, samples, groups)
    stages = [t["stage"] for t in tasks]
    assert stages.count("haplotypecaller") == 4  # 2 samples x 2 contigs
    assert stages.count("genomicsdb") == 2
    assert stages.count("genotype") == 2
    assert stages.count("concat") == 1
    assert stages.count("multiqc") == 1
    assert stages.count("align") == (2 if mode == "fastq" else 0)
    assert stages.count("markdup") == (2 if mode == "fastq" else 0)
    assert stages.count("import") == (0 if mode == "fastq" else 2)
    assert len(maps) == 2
    assert not workdir.exists(), "Planning must not create outputs"

    by_id = {t["id"]: t for t in tasks}
    assert all(all(dep in by_id for dep in t["deps"]) for t in tasks)
    assert len(by_id) == len(tasks)
    assert [t["stage"] for t in tasks] == sorted(stages, key=workflow.STAGES.index)
    final = by_id["concat.Test_species"]
    assert len(final["deps"]) == 2
    assert final["outputs"][0].endswith("cohort.raw.vcf.gz")


def test_trim_and_lane_merging(tmp_path):
    path, _ = make_example(tmp_path, trim=True)
    cfg = yaml.safe_load(path.read_text())
    sample_table = Path(cfg["samples"])
    first = sample_table.read_text().splitlines()[1].split("\t")
    first[4] = "lane2"  # change lane field
    new_r1 = tmp_path / "second_R1.fq.gz"
    new_r2 = tmp_path / "second_R2.fq.gz"
    new_r1.write_text("reads")
    new_r2.write_text("reads")
    first[5], first[6] = str(new_r1), str(new_r2)
    sample_table.write_text(sample_table.read_text() + "\t".join(first) + "\n")
    _, settings, workdir = workflow.config_load(str(path))
    units, samples, groups = workflow.validate(settings)
    tasks, _ = workflow.make_tasks(settings, workdir, units, samples, groups)
    assert sum(t["stage"] == "trim" for t in tasks) == 3
    assert sum(t["stage"] == "align" for t in tasks) == 3
    assert len(next(t for t in tasks if t["id"] == "markdup.S1")["deps"]) == 2


def test_missing_data_validation(tmp_path):
    path, _ = make_example(tmp_path)
    Path(yaml.safe_load(path.read_text())["samples"]).unlink()
    _, cfg, _ = workflow.config_load(str(path))
    with pytest.raises(ValueError, match="Missing sample table"):
        workflow.validate(cfg)


def test_plan_resume_and_signatures(tmp_path):
    config, _ = make_example(tmp_path)
    path, cfg, workdir = workflow.config_load(str(config))
    units, samples, groups = workflow.validate(cfg)
    tasks, maps = workflow.make_tasks(cfg, workdir, units, samples, groups)
    assert all(pipeline.pending_tasks(tasks, workdir).values())
    plan = pipeline.write_plan(tasks, maps, path, workdir)
    saved = json.loads(plan.read_text())
    assert len(saved["tasks"]) == len(tasks)
    assert (workdir / "jobs" / "sample_maps" / "Test_species" / "chr1.tsv").exists()
    assert len(list((workdir / "jobs" / "scripts").glob("*.sh"))) == len(tasks)
    assert all(pipeline.pending_tasks(tasks, workdir).values())  # no stage was executed

    # A completion record must correspond to nonempty outputs and correct signature.
    first = tasks[0]
    for out in first["outputs"]:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text("ok")
    marker = workdir / "checkpoints" / f"{first['id']}.json"
    marker.write_text(json.dumps({"signature": first["signature"]}))
    assert pipeline.is_completed(first, workdir)
    marker.write_text(json.dumps({"signature": "obsolete-code"}))
    assert not pipeline.is_completed(first, workdir)


def test_invalid_resource_settings(tmp_path):
    path, cfg = make_example(tmp_path)
    cfg["slurm"]["resources"] = {"align": {"cpus": 2}}
    path.write_text(yaml.safe_dump(cfg))
    _, loaded, workdir = workflow.config_load(str(path))
    units, samples, groups = workflow.validate(loaded)
    with pytest.raises(ValueError, match="at least 8"):
        workflow.make_tasks(loaded, workdir, units, samples, groups)


def test_slurm_arrays_and_dependencies(tmp_path, monkeypatch):
    """Mock sbatch: verify stage arrays and afterok dependencies, no real jobs."""
    config, _ = make_example(tmp_path, mode="fastq", trim=True)
    path, cfg, workdir = workflow.config_load(str(config))
    units, samples, groups = workflow.validate(cfg)
    tasks, sample_maps = workflow.make_tasks(cfg, workdir, units, samples, groups)
    plan = pipeline.write_plan(tasks, sample_maps, path, workdir)

    calls = []

    class Result:
        returncode = 0
        stderr = ""

        def __init__(self, job_id):
            self.stdout = f"{job_id}\n"

    def fake_run(args, **kwargs):
        calls.append(args)
        return Result(1000 + len(calls))

    monkeypatch.setattr(pipeline.shutil, "which", lambda _: "/usr/bin/sbatch")
    monkeypatch.setattr(pipeline.subprocess, "run", fake_run)
    pending = pipeline.pending_tasks(tasks, workdir)
    job_ids = pipeline.submit_to_slurm(tasks, pending, workdir, plan, cfg)
    assert list(job_ids) == [stage for stage in workflow.STAGES
                             if any(t["stage"] == stage for t in tasks)]
    assert len(calls) == len(job_ids)
    assert all("--parsable" in cmd for cmd in calls)
    assert any("--dependency=afterok:1001" in arg for arg in calls[1])
    assert len((workdir / "jobs" / "array_align.txt").read_text().splitlines()) == 2
    assert len((workdir / "jobs" / "array_haplotypecaller.txt").read_text().splitlines()) == 4
    assert json.loads((workdir / "jobs" / "submission.json").read_text())["job_ids"] == job_ids


def test_local_task_execution_and_resume(tmp_path):
    """Execute a fake, harmless task to test checkpoints without scientific tools."""
    workdir = tmp_path / "results"
    result = workdir / "dummy.txt"
    task = {
        "id": "dummy", "stage": "multiqc", "deps": [], "inputs": [],
        "outputs": [str(result)],
        "signature": "demo-signature",
        "body": f"set -euo pipefail\nprintf 'ok\\n' > {result}\n",
    }
    plan = pipeline.write_plan([task], {}, tmp_path / "config.yaml", workdir)
    assert pipeline.pending_tasks([task], workdir)["dummy"]
    pipeline.run_one_task(plan, "dummy")
    assert result.read_text() == "ok\n"
    assert not pipeline.pending_tasks([task], workdir)["dummy"]
    result.unlink()
    assert pipeline.pending_tasks([task], workdir)["dummy"]
