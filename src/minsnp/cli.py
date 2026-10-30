"""Единая точка входа: локальный последовательный режим и простой SLURM job."""
from __future__ import annotations
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from . import joint, qc, wgs
from .resources import configure_resources
from .core import ROOT, capture_versions, config, file_exists, prepare_run, require_commands, samples_in_vcf, source_fingerprint


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Minimal SNP Panels: FASTQ или multi-sample VCF/BCF")
    p.add_argument("command", choices=["run", "wgs", "wgs-vcf", "joint", "qc", "structure", "panels", "check-tools"])
    p.add_argument("--species", choices=["beluga", "polar_bear"])
    source = p.add_mutually_exclusive_group()
    source.add_argument("--fastq", help="TSV: sample_id species run_id library fastq_r1 fastq_r2")
    source.add_argument("--vcf", help="Индексированный multi-sample VCF.gz или BCF")
    source.add_argument("--gvcf-manifest", help="Внутренний вход для independent joint")
    p.add_argument("--labels", help="TSV: sample_id population [group]")
    p.add_argument("--mode", choices=["structure", "panels", "both"], default="both")
    p.add_argument("--run-id", help="Имя версии запуска")
    p.add_argument("--outdir", default="results")
    p.add_argument("--references", default=str(ROOT / "config/references.yaml"))
    p.add_argument("--settings", default=str(ROOT / "config/pipeline.yaml"))
    p.add_argument("--backend", choices=["local", "slurm"], default="local")
    p.add_argument("--threads", type=int, help="CPU >= 1, переопределяет threads из YAML")
    p.add_argument("--memory", help="RAM для job: например 32G, 64G, 128G")
    p.add_argument("--resume", action="store_true")
    return p


def submit_slurm(args, cfg: dict, run: Path) -> None:
    """Один SLURM job на WGS→VCF, с явно выделенными CPU и RAM."""
    require_commands("sbatch")
    slurm_cfg = cfg["slurm"]
    logs = run / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    script = run / "slurm_job.sh"
    previous = run / "slurm_job_id.txt"
    if previous.is_file():
        prior_id = previous.read_text().strip().split(";")[0]
        active = subprocess.run(["squeue", "-h", "-j", prior_id], capture_output=True, text=True)
        if active.returncode == 0 and active.stdout.strip():
            raise ValueError(f"SLURM job {prior_id} уже выполняется для run-id={run.name}")

    # Внутри job запускаем тот же CLI в local mode. Путь и среда остаются прежними.
    arguments = sys.argv[1:]
    arguments = ["local" if x == "slurm" else x.replace("--backend=slurm", "--backend=local") for x in arguments]
    if "--run-id" not in arguments and not any(x.startswith("--run-id=") for x in arguments):
        arguments.extend(["--run-id", run.name])
    if "--resume" not in arguments:
        arguments.append("--resume")
    if args.threads is None:
        arguments.extend(["--threads", str(cfg["_requested_threads"])])
    if args.memory is None:
        arguments.extend(["--memory", cfg["_requested_memory"]])
    script.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + f"export PYTHONPATH={shlex.quote(str(ROOT / 'src'))}\n"
        + shlex.join([sys.executable, "-m", "minsnp.cli", *arguments]) + "\n"
    )
    script.chmod(0o750)
    command = [
        "sbatch", "--parsable", "--nodes=1", "--ntasks=1",
        "--job-name", f"minsnp_{args.species}",
        "--cpus-per-task", str(cfg["_requested_threads"]),
        "--mem", cfg["_requested_memory"],
        "--time", str(slurm_cfg["time"]),
        "--chdir", str(ROOT),
        "--output", str(logs / "slurm_%j.out"),
        "--error", str(logs / "slurm_%j.err"),
    ]
    if slurm_cfg.get("partition"):
        command.extend(["--partition", str(slurm_cfg["partition"])])
    command.append(str(script))
    job_id = subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip()
    previous.write_text(job_id + "\n")
    print("SLURM job:", job_id)

def main() -> None:
    args = parser().parse_args()
    if args.command == "check-tools":
        print(json.dumps(capture_versions(["python", "bwa-mem3", "samtools", "bcftools", "gatk", "fastp", "mosdepth", "plink2", "sbatch"]), indent=2))
        return
    if not args.species:
        raise ValueError("Укажите --species")
    cfg, spec = config(Path(args.references), Path(args.settings), args.species)
    cfg = configure_resources(args, cfg)
    if (args.command in {"wgs", "wgs-vcf"} and not args.fastq) or (args.command == "joint" and not (args.gvcf_manifest or args.fastq)):
        raise ValueError("Для wgs нужен --fastq; для joint нужен --gvcf-manifest или --fastq")
    if args.command in {"qc", "structure", "panels"} and not args.vcf:
        raise ValueError("Для этого модуля нужен --vcf")
    if args.command == "run" and not (args.fastq or args.vcf):
        raise ValueError("Для run используйте --fastq или --vcf")
    if args.command == "panels" and not args.labels:
        raise ValueError("Для panels нужен --labels")
    inputs = {"command": args.command, "fastq": str(Path(args.fastq).resolve()) if args.fastq else None,
              "vcf": str(Path(args.vcf).resolve()) if args.vcf else None,
              "gvcf": str(Path(args.gvcf_manifest).resolve()) if args.gvcf_manifest else None,
              "labels": str(Path(args.labels).resolve()) if args.labels else None,
              "mode": args.mode, "ref": str(Path(args.references).resolve()), "settings": str(Path(args.settings).resolve()),
              "checksums": {k: source_fingerprint(v) for k, v in {"fastq": args.fastq, "vcf": args.vcf, "gvcf": args.gvcf_manifest, "labels": args.labels, "references": args.references, "settings": args.settings}.items()}}
    inputs["resources"] = {"threads": cfg["_requested_threads"], "memory": cfg["_requested_memory"]}
    run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    run = prepare_run(Path(args.outdir).resolve(), args.species, run_id, inputs, args.resume)
    if args.backend == "slurm":
        submit_slurm(args, cfg, run)
        return
    logs = run / "logs"
    logs.mkdir(exist_ok=True)
    (run / "versions.json").write_text(json.dumps(capture_versions(["bwa-mem3", "samtools", "bcftools", "gatk", "fastp", "mosdepth", "plink2"]), indent=2))
    if args.fastq:
        require_commands("bwa-mem3", "samtools", "gatk", "bcftools", "mosdepth")
        manifest = wgs.run(args.fastq, args.species, spec, cfg, run / "wgs", logs, args.resume)
        if args.command == "wgs":
            print(manifest)
            return
        vcf = joint.run(manifest, spec, run / "joint", logs, args.resume)
    elif args.gvcf_manifest:
        require_commands("bcftools", "gatk")
        vcf = joint.run(Path(args.gvcf_manifest), spec, run / "joint", logs, args.resume)
    else:
        vcf = file_exists(args.vcf)
        if len(samples_in_vcf(vcf)) < 2:
            raise ValueError("Ожидается multi-sample VCF/BCF с минимум двумя образцами")
    if args.command == "joint":
        print(vcf)
        return
    if args.command == "wgs-vcf":
        from .core import check_index, stage
        exported = run / "joint" / "cohort.raw.vcf.gz"
        if not args.resume or not exported.is_file():
            stage("export_vcf.sh", [vcf, exported], logs, "export_vcf")
        check_index(exported)
        print(exported)
        return
    require_commands("bcftools")
    filtered = qc.run(vcf, spec["fasta"], spec["contigs"], cfg, run / "qc", logs, args.resume)
    if args.command == "qc":
        print(filtered)
        return
    labels = Path(args.labels) if args.labels else None
    if args.command == "structure" or (args.command == "run" and labels is None and args.mode in {"structure", "both"}):
        require_commands("plink2")
        from . import structure
        labels = structure.run(filtered, spec["contigs"], cfg, run / "structure", logs, args.resume)
    if args.command == "structure" or (args.command == "run" and args.mode == "structure"):
        print(labels)
        return
    if labels is None:
        raise ValueError("Для поиска панелей нужны --labels или режим --mode both")
    from . import panels
    result = panels.run(filtered, labels, cfg, run / "panels", logs, spec["assembly"])
    print(result)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ОШИБКА: {exc}", file=sys.stderr)
        sys.exit(1)
