"""Offline BWA-MEM3 migration checks (no alignment or SLURM required)."""
import json
import os
from pathlib import Path
import subprocess

import pytest

from test_pipeline import make_example
import workflow
import pipeline


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "stages" / "align_lane.sh"


def _plan(tmp_path, mode="fastq"):
    config, _ = make_example(tmp_path, mode=mode)
    _, cfg, workdir = workflow.config_load(str(config))
    units, samples, groups = workflow.validate(cfg)
    tasks, _ = workflow.make_tasks(cfg, workdir, units, samples, groups)
    return tasks, workdir


def test_bwa_v1_index_alone_is_not_accepted(tmp_path):
    config, cfg = make_example(tmp_path)
    ref = Path(cfg["references"]["Test_species"]["fasta"])
    (tmp_path / "species.fa.bwt.2bit.64").unlink()
    Path(str(ref) + ".bwt").write_text("BWA v1 index")
    Path(str(ref) + ".sa").write_text("BWA v1 index")
    _, loaded, _ = workflow.config_load(str(config))
    with pytest.raises(ValueError, match=r"bwt\.2bit\.64"):
        workflow.validate(loaded)


def test_ready_bam_does_not_require_aligner_index(tmp_path):
    config, cfg = make_example(tmp_path, mode="bam")
    ref = Path(cfg["references"]["Test_species"]["fasta"])
    for suffix in (".amb", ".ann", ".pac", ".bwt.2bit.64"):
        Path(str(ref) + suffix).unlink()
    _, loaded, workdir = workflow.config_load(str(config))
    units, samples, groups = workflow.validate(loaded)
    tasks, _ = workflow.make_tasks(loaded, workdir, units, samples, groups)
    assert not any(t["stage"] == "align" for t in tasks)


def test_bwa_mem3_index_is_an_alignment_input(tmp_path):
    tasks, workdir = _plan(tmp_path)
    align = [t for t in tasks if t["stage"] == "align"]
    assert align
    for task in align:
        assert any(str(p).endswith(".bwt.2bit.64") for p in task["inputs"])
        assert not any(str(p).endswith(".bwt") or str(p).endswith(".sa")
                       for p in task["inputs"])
    # Future checkpoint validation must observe this index as an external input.
    for task in tasks:
        for output in task["outputs"]:
            path = Path(output)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("ok")
        if task["stage"] in ("haplotypecaller", "genotype"):
            Path(task["outputs"][0] + ".tbi").write_text("index")
        cp = workdir / "checkpoints" / (task["id"] + ".json")
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps({"signature": task["signature"]}))
    assert not any(pipeline.pending_tasks(tasks, workdir).values())
    index = tmp_path / "species.fa.bwt.2bit.64"
    later = index.stat().st_mtime + 600
    os.utime(index, (later, later))
    pending = pipeline.pending_tasks(tasks, workdir)
    assert all(pending[t["id"]] for t in align)
    assert pending["concat.Test_species"]
    assert pending["multiqc"]  # QC depends on markdup, and MultiQC depends on QC


def test_alignment_invokes_mem3_compat_and_preserves_samtools_chain(tmp_path):
    """Mock executables, no scientific processes or FASTQ reads."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "bwa_args.json"
    bwa = bindir / "bwa-mem3"
    bwa.write_text("""#!/usr/bin/env python3
import json, os, sys
with open(os.environ['MOCK_BWA_ARGS'], 'w') as f:
    json.dump(sys.argv[1:], f)
sys.stdout.write('@HD\\tVN:1.6\\n')
""")
    bwa.chmod(0o755)
    samtools = bindir / "samtools"
    samtools.write_text("""#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
if args[0] == 'sort':
    pathlib.Path(args[args.index('-o') + 1]).write_bytes(sys.stdin.buffer.read())
else:
    sys.stdout.buffer.write(sys.stdin.buffer.read())
""")
    samtools.chmod(0o755)
    out = tmp_path / "lane.bam"
    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env["PATH"]
    env["MOCK_BWA_ARGS"] = str(log)
    proc = subprocess.run([
        "bash", str(SCRIPT), "/mock/ref.fa", "/mock/R1.fq.gz",
        "/mock/R2.fq.gz", str(out), "S1", "L1", "lane1", "8",
    ], capture_output=True, text=True, env=env, check=True)
    assert out.read_bytes().startswith(b"@HD")
    args = json.loads(log.read_text())
    assert args[0:4] == ["mem", "--compat=bwa-mem", "-t", "4"]
    assert args[4] == "-R"
    assert r"\tSM:S1" in args[5]
    assert args[6:] == ["/mock/ref.fa", "/mock/R1.fq.gz", "/mock/R2.fq.gz"]
