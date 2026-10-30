"""Offline tests for stage tool names, packaging, and trim outputs."""
import os
from pathlib import Path
import subprocess

import yaml

ROOT = Path(__file__).resolve().parents[1]
STAGES = ROOT / "scripts" / "stages"


def test_updated_environment_versions():
    cfg = yaml.safe_load((ROOT / "envs" / "wgs.yaml").read_text())
    deps = set(cfg["dependencies"])
    assert {"bwa-mem3=0.15.1", "gatk4-main=4.7.0.0", "fastp=1.4.0",
            "mosdepth=0.3.14", "multiqc=1.35", "samtools=1.24",
            "bcftools=1.24"} <= deps
    assert not any(x.startswith("gatk4=") for x in deps)


def test_all_ten_stage_scripts_parse():
    scripts = sorted(STAGES.glob("*.sh"))
    assert len(scripts) == 10
    for script in scripts:
        subprocess.run(["bash", "-n", str(script)], check=True)


def test_trim_creates_all_output_directories_without_real_fastp(tmp_path):
    mock_bin = tmp_path / "bin"
    mock_bin.mkdir()
    fastp = mock_bin / "fastp"
    fastp.write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
for flag in ('-o', '-O', '-j', '-h'):
    Path(args[args.index(flag) + 1]).write_text(flag)
""")
    fastp.chmod(0o755)
    r1 = tmp_path / "input_R1.fq.gz"
    r2 = tmp_path / "input_R2.fq.gz"
    r1.write_text('')
    r2.write_text('')
    outputs = [tmp_path / "r1" / "1.fq.gz", tmp_path / "r2" / "2.fq.gz",
               tmp_path / "qc" / "trim.json", tmp_path / "html" / "trim.html"]
    env = os.environ.copy()
    env["PATH"] = str(mock_bin) + os.pathsep + env["PATH"]
    subprocess.run(["bash", str(STAGES / "trim.sh"), str(r1), str(r2),
                    *map(str, outputs)], check=True, env=env)
    assert all(x.exists() for x in outputs)


def test_stage_binary_names_intact():
    expected = {
        "trim.sh": "fastp -i",
        "align_lane.sh": "bwa-mem3 mem --compat=bwa-mem",
        "markdup.sh": "samtools markdup",
        "bam_qc.sh": "mosdepth --threads",
        "import_aligned.sh": "samtools quickcheck",
        "haplotypecaller.sh": "HaplotypeCaller",
        "genomicsdb.sh": "GenomicsDBImport",
        "genotype.sh": "GenotypeGVCFs",
        "concat.sh": "bcftools concat",
        "multiqc.sh": "multiqc --force",
    }
    for filename, command in expected.items():
        assert command in (STAGES / filename).read_text()
