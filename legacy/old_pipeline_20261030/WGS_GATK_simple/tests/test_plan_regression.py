"""Golden DAG snapshots after the intentional BWA-MEM3 aligner migration.

No bioinformatics commands run. Fixture paths are normalized for portability;
the raw task signatures are tested separately because they include PROJECT paths.
"""
import hashlib
import json
from pathlib import Path

import pytest
import yaml

from test_pipeline import make_example
import workflow

CASES = (
    "fastq_untrimmed", "fastq_trimmed", "fastq_multilane",
    "bam", "cram", "beluga_bear_autosomes",
)
GOLDEN = Path(__file__).with_name("golden_plan_digests.json")


def make_case(tmp_path, case):
    mode = "bam" if case == "bam" else "cram" if case == "cram" else "fastq"
    trim = case in ("fastq_trimmed", "fastq_multilane")
    config, settings = make_example(tmp_path, mode=mode, trim=trim)
    if case == "fastq_multilane":
        sample_table = Path(settings["samples"])
        first = sample_table.read_text().splitlines()[1].split("\t")
        first[4] = "lane2"
        for mate, col in (("R1", 5), ("R2", 6)):
            path = tmp_path / f"S1_lane2_{mate}.fq.gz"
            path.write_text("test")
            first[col] = str(path)
        with sample_table.open("a") as out:
            out.write("\t".join(first) + "\n")
    if case == "beluga_bear_autosomes":
        # Synthetic reference headers: actual numbers of autosomes for both species.
        references = {}
        for species, count in (("Beluga_species", 21), ("Bear_species", 36)):
            fasta = tmp_path / f"{species}.fa"
            contigs = [f"ctg{i:02}" for i in range(1, count + 1)]
            fasta.write_text("".join(f">{name}\nA\n" for name in contigs))
            Path(str(fasta) + ".fai").write_text("".join(
                f"{name}\t1\t0\t1\t2\n" for name in contigs
            ))
            fasta.with_suffix(".dict").write_text(
                "@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{name}\tLN:1\n" for name in contigs)
            )
            for suffix in ("amb", "ann", "bwt.2bit.64", "pac"):
                Path(str(fasta) + "." + suffix).write_text("index")
            references[species] = {"fasta": str(fasta), "contigs": contigs}
        settings["references"] = references
        sample_table = Path(settings["samples"])
        sample_table.write_text(sample_table.read_text().replace("Test_species", "Beluga_species", 1).replace("Test_species", "Bear_species", 1))
        config.write_text(yaml.safe_dump(settings))
    return config


def build_plan(tmp_path, case):
    config = make_case(tmp_path, case)
    _, cfg, workdir = workflow.config_load(str(config))
    units, samples, groups = workflow.validate(cfg)
    return workflow.make_tasks(cfg, workdir, units, samples, groups)


def normalized_snapshot(tmp_path, tasks, maps):
    # Signatures depend on absolute code path; bodies and all other fields are
    # compared across machines, including dependencies, ordering and map text.
    serializable = {"tasks": [{**task, "signature": "<path-dependent SHA256>"}
                              for task in tasks], "sample_maps": maps}
    raw = json.dumps(serializable, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return raw.replace(str(tmp_path), "<FIXTURE>").replace(str(workflow.PROJECT), "<PROJECT>")


@pytest.mark.parametrize("case", CASES)
def test_golden_dag(tmp_path, case):
    tasks, maps = build_plan(tmp_path, case)
    snapshot = normalized_snapshot(tmp_path, tasks, maps)
    digest = hashlib.sha256(snapshot.encode()).hexdigest()
    golden = json.loads(GOLDEN.read_text())[case]
    assert len(tasks) == golden["task_count"]
    assert len(maps) == golden["map_count"]
    assert digest == golden["sha256"], f"Workflow DAG changed in {case}"


def test_signatures_match_original_calculation(tmp_path):
    tasks, maps = build_plan(tmp_path, "fastq_trimmed")
    for task in tasks:
        stage_script = workflow.PROJECT / "scripts" / "stages" / workflow.SCRIPT[task["stage"]]
        if task["stage"] == "genomicsdb":
            # Original signature incorporates the full contents of its sample map.
            species, contig = task["id"].split(".", 2)[1:]
            map_file = str(tmp_path / "output" / "jobs" / "sample_maps" / species / f"{contig}.tsv")
            extra = maps[map_file]
        else:
            extra = ""
        digest = hashlib.sha256(task["body"].encode() + stage_script.read_bytes() + extra.encode())
        if task["stage"] == "import":
            digest.update((workflow.PROJECT / "scripts" / "validate_alignment.py").read_bytes())
        assert task["signature"] == digest.hexdigest(), task["id"]


def test_resume_propagates_modified_input(tmp_path):
    """A stale FASTQ invalidates all affected descendants, not other samples."""
    import json
    import os
    import time
    import pipeline

    tasks, _ = build_plan(tmp_path, "fastq_trimmed")
    workdir = tmp_path / "output"
    for task in tasks:
        for filename in task["outputs"]:
            output = Path(filename)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text("ok")
        if task["stage"] in ("haplotypecaller", "genotype"):
            Path(task["outputs"][0] + ".tbi").write_text("index")
        marker = workdir / "checkpoints" / (task["id"] + ".json")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(json.dumps({"signature": task["signature"]}))

    assert not any(pipeline.pending_tasks(tasks, workdir).values())
    target = tmp_path / "S1_r1.fq.gz"
    later = time.time() + 30
    os.utime(target, (later, later))
    pending = pipeline.pending_tasks(tasks, workdir)
    assert pending["trim.S1.L1.lane1"]
    assert pending["align.S1.L1.lane1"]
    assert pending["markdup.S1"]
    assert pending["hc.Test_species.S1.chr1"]
    assert pending["db.Test_species.chr1"]
    assert pending["genotype.Test_species.chr1"]
    assert pending["concat.Test_species"]
    assert pending["multiqc"]
    assert not pending["trim.S2.L1.lane1"]
    assert not pending["align.S2.L1.lane1"]
    assert not pending["markdup.S2"]


def test_legacy_reference_non_aligning_snapshots_unchanged():
    """BAM/CRAM plans must remain identical to legacy BWA-based pipeline."""
    legacy = json.loads(Path(__file__).with_name("golden_plan_digests_bwa_legacy.json").read_text())
    new = json.loads(GOLDEN.read_text())
    for case in ("bam", "cram"):
        assert new[case] == legacy[case]
    for case in ("fastq_untrimmed", "fastq_trimmed", "fastq_multilane", "beluga_bear_autosomes"):
        assert new[case]["sha256"] != legacy[case]["sha256"]
