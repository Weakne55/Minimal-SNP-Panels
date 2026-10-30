#!/usr/bin/env python3

"""
Управление WGS-пайплайном: проверка, планирование и запуск задач.

Команды для пользователя:
    ./run.sh validate  — проверить входные данные и конфигурацию;
    ./run.sh plan      — показать количество запланированных задач;
    ./run.sh local     — выполнить задачи последовательно на текущем узле;
    ./run.sh submit    — отправить задачи в SLURM;
    ./run.sh status    — показать активные задания SLURM.

Команда run-task используется внутри пайплайна для запуска одной задачи.
Биоинформатические команды определены в scripts/stages/*.sh.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from workflow import PROJECT, STAGES, config_load, validate, make_tasks

RESOURCE_DEFAULTS = {
    "trim": (4, 8000, 4),
    "align": (8, 20000, 16),
    "markdup": (8, 20000, 12),
    "import": (4, 12000, 8),
    "qc": (4, 8000, 8),
    "haplotypecaller": (4, 24000, 24),
    "genomicsdb": (2, 32000, 18),
    "genotype": (2, 16000, 12),
    "concat": (2, 8000, 4),
    "multiqc": (2, 8000, 2),
}

def missing_output(task):
    """Return the first absent/empty output, including an index for GATK VCFs.

    A VCF index is accepted in any of the three formats used by the
    original pipeline; all other stages only require their declared outputs.
    """
    for filename in task["outputs"]:
        path = Path(filename)
        if not path.is_file() or path.stat().st_size == 0:
            return filename

    if task["stage"] in ("haplotypecaller", "genotype"):
        vcf = task["outputs"][0]
        if not any(
            Path(vcf + suffix).is_file() and Path(vcf + suffix).stat().st_size > 0
            for suffix in (".tbi", ".idx", ".csi")
        ):
            return f"{vcf} (index)"

    return None


def is_completed(task, workdir):
    """
    Проверить, можно ли использовать результат ранее выполненной задачи.

    Задача считается завершённой, если:
      1. существует checkpoint с подходящей сигнатурой;
      2. все ожидаемые выходные файлы существуют и не пусты;
      3. для GVCF/VCF после GATK есть файл индекса;
      4. входные файлы не новее checkpoint.

    Сигнатура формируется в workflow.py и позволяет обнаруживать изменения
    параметров или команд. При любой ошибке чтения проверка возвращает False.
    """

    marker = workdir / "checkpoints" / f"{task['id']}.json"

    try:
        record = json.loads(marker.read_text(encoding="utf8"))

        if record["signature"] != task["signature"]:
            return False

        finished_at = marker.stat().st_mtime

        if missing_output(task) is not None:
            return False

        if any(
            not Path(filename).is_file()
            or Path(filename).stat().st_mtime > finished_at
            for filename in task["inputs"]
        ):
            return False

        return True

    except (OSError, ValueError, KeyError):
        return False

def pending_tasks(tasks, workdir, force=False):
    """
    Определить, какие задачи необходимо выполнить (или повторить).

    Задача попадает в pending, если она не завершена либо зависит от задачи,
    которую нужно пересчитать. Список tasks должен идти в порядке зависимостей.
    При force=True все задачи помечаются для повторного выполнения.
    """

    pending = {}

    for task in tasks:
        pending[task["id"]] = (
            force
            or not is_completed(task, workdir)
            or any(pending[dependency] for dependency in task["deps"])
        )

    return pending

def write_plan(tasks, sample_maps, config_file, workdir):
    """
    Сохранить скрипты задач, карты образцов GVCF и общий JSON-план.

    Вызывается только перед local/submit, а не при validate/plan.
    Возвращает путь к pipeline_plan.json.
    """

    (workdir / "jobs" / "scripts").mkdir(parents=True, exist_ok=True)
    (workdir / "jobs" / "logs").mkdir(parents=True, exist_ok=True)
    (workdir / "checkpoints").mkdir(parents=True, exist_ok=True)

    for task in tasks:
        for output in task["outputs"]:
            Path(output).parent.mkdir(parents=True, exist_ok=True)

        filename = workdir / "jobs" / "scripts" / f"{task['id']}.sh"
        filename.write_text(task["body"], encoding="utf8")
        filename.chmod(0o750)

    for filename, contents in sample_maps.items():
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf8")

    plan = workdir / "pipeline_plan.json"
    temporary = workdir / "pipeline_plan.json.tmp"

    temporary.write_text(
        json.dumps(
            {
                "root": str(PROJECT),
                "workdir": str(workdir),
                "config": str(config_file),
                "tasks": tasks,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf8",
    )
    temporary.replace(plan)

    return plan

def run_one_task(plan_file, task_id):
    """
    Выполнить одну задачу по её идентификатору из сохранённого плана.

    Используется как при последовательном запуске, так и на узле SLURM.
    Отметка об успешном завершении создаётся только после проверки выходов.
    """

    plan = json.loads(Path(plan_file).read_text(encoding="utf8"))

    if Path(plan["root"]).resolve() != PROJECT:
        raise ValueError("Plan belongs to another project location; regenerate it")

    workdir = Path(plan["workdir"])
    task = next((t for t in plan["tasks"] if t["id"] == task_id), None)

    if task is None:
        raise ValueError(f"Unknown task: {task_id}")

    script = workdir / "jobs" / "scripts" / f"{task_id}.sh"
    if script.read_text(encoding="utf8") != task["body"]:
        raise ValueError(f"Saved task script has changed: {script}")

    marker = workdir / "checkpoints" / f"{task_id}.json"
    marker.unlink(missing_ok=True)

    for filename in task["outputs"]:
        Path(filename).unlink(missing_ok=True)

    print(f"RUN  {task_id}", flush=True)
    subprocess.run(["bash", str(script)], cwd=PROJECT, check=True)

    missing = missing_output(task)
    if missing is not None:
        if missing.endswith(" (index)"):
            raise RuntimeError(f"{task_id}: missing VCF index")
        raise RuntimeError(f"{task_id}: missing or empty output {missing}")

    temporary = marker.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(
            {
                "signature": task["signature"],
                "finished": time.time(),
            }
        )
    )
    temporary.replace(marker)

    print(f"DONE {task_id}", flush=True)

def active_jobs(workdir):
    """
    Получить список активных SLURM Job ID из последней отправки.

    Список Job ID берётся из jobs/submission.json, затем проверяется через
    squeue. Пустой результат означает отсутствие видимых активных заданий.
    """

    record = workdir / "jobs" / "submission.json"

    if not record.exists() or not shutil.which("squeue"):
        return []

    ids = list(json.loads(record.read_text()).get("job_ids", {}).values())
    if not ids:
        return []

    result = subprocess.run(
        ["squeue", "-h", "-j", ",".join(map(str, ids)), "-o", "%A"],
        capture_output=True,
        text=True,
        check=True,
    )

    return sorted(set(result.stdout.split()))

def submit_to_slurm(tasks, pending, workdir, plan, cfg):
    """
    Отправить невыполненные задачи в SLURM, сгруппировав их по этапам.

    Для каждого этапа создаётся один SLURM job array. Зависимость afterok
    разрешает следующий этап только после успешного завершения предыдущих
    массивов, от которых зависят его задачи.

    Такая схема проста, но синхронизирует задачи на уровне целого массива,
    а не отдельного образца или хромосомы.
    """

    if not shutil.which("sbatch"):
        raise RuntimeError("sbatch not found: run on a SLURM login node")

    slurm = cfg.get("slurm", {})
    parallel = int(slurm.get("max_array_parallel", 4))

    if parallel < 1:
        raise ValueError("slurm.max_array_parallel must be >= 1")

    overrides = slurm.get("resources", {})
    jobs = {}  # Название этапа -> полученный SLURM Job ID.
    by_id = {task["id"]: task for task in tasks}

    # Group pending tasks once, keeping the order used by job arrays.
    pending_by_stage = {stage: [] for stage in STAGES}
    for task in tasks:
        if pending[task["id"]]:
            pending_by_stage[task["stage"]].append(task)

    for stage in STAGES:
        selected = pending_by_stage[stage]
        if not selected:
            continue

        previous_stages = {
            by_id[dependency]["stage"]
            for task in selected
            for dependency in task["deps"]
            if pending[dependency] and by_id[dependency]["stage"] != stage
        }

        dependencies = [
            jobs[name] for name in STAGES if name in previous_stages
        ]

        task_list = workdir / "jobs" / f"array_{stage}.txt"
        task_list.write_text("\n".join(t["id"] for t in selected) + "\n")

        default_cpus, default_mem, default_hours = RESOURCE_DEFAULTS[stage]
        settings = overrides.get(stage, {})

        cpus = int(settings.get("cpus", default_cpus))
        mem = int(settings.get("mem_mb", default_mem))
        hours = int(settings.get("hours", default_hours))

        if min(cpus, mem, hours) < 1:
            raise ValueError(f"{stage}: SLURM resources must be positive")

        command = [
            "sbatch",
            "--parsable",
            f"--job-name=WGS_{stage}",
            f"--cpus-per-task={cpus}",
            f"--mem={mem}M",
            f"--time={hours}:00:00",
            f"--array=0-{len(selected) - 1}%{parallel}",
            f"--output={workdir}/jobs/logs/{stage}_%A_%a.out",
            f"--error={workdir}/jobs/logs/{stage}_%A_%a.err",
            f"--export=ALL,PIPELINE_PYTHON={sys.executable}",
        ]

        for key in ("partition", "account", "qos"):
            if slurm.get(key):
                command.append(f"--{key}={slurm[key]}")

        if dependencies:
            command.append("--dependency=afterok:" + ":".join(dependencies))

        command.extend(
            [
                str(PROJECT / "scripts" / "run_array.sh"),
                str(plan),
                str(task_list),
            ]
        )

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode:
            raise RuntimeError(
                f"SLURM rejected stage {stage}: {result.stderr.strip()}"
            )

        job_id = result.stdout.strip().split(";")[0]
        if not re.fullmatch(r"\d+", job_id):
            raise RuntimeError(f"Unexpected sbatch response: {result.stdout!r}")

        jobs[stage] = job_id

        (workdir / "jobs" / "submission.json").write_text(
            json.dumps(
                {
                    "job_ids": jobs,
                    "plan": str(plan),
                },
                indent=2,
            )
        )

        print(f"{stage:16} {len(selected):5} tasks -> SLURM {job_id}")

    return jobs

def main():
    """Обработать команду пользователя и выбрать режим работы пайплайна."""

    parser = argparse.ArgumentParser(prog="./run.sh", description=__doc__)
    subcommands = parser.add_subparsers(dest="action", required=True)

    descriptions = {
        "validate": "проверить входные данные и настройки (ничего не запускает)",
        "plan": "показать план и число задач (ничего не запускает)",
        "local": "выполнить задачи последовательно на текущем узле",
        "submit": "запустить или продолжить задачи через SLURM",
        "status": "показать активные задания SLURM",
    }

    for action, description in descriptions.items():
        sub = subcommands.add_parser(action, help=description)
        sub.add_argument("--config", default="config/config.yaml")

        if action in ("plan", "local", "submit"):
            sub.add_argument(
                "--force",
                action="store_true",
                help="Rerun every task (never during active jobs)",
            )

        if action == "plan":
            sub.add_argument(
                "--allow-missing-inputs",
                action="store_true",
                help="Preview before downloading FASTQs and references",
            )

    worker = subcommands.add_parser("run-task", help="служебная команда для SLURM")
    worker.add_argument("--plan", required=True)
    worker.add_argument("--task-id", required=True)

    args = parser.parse_args()

    if args.action == "run-task":
        run_one_task(args.plan, args.task_id)
        return

    config_file, cfg, workdir = config_load(args.config)

    if args.action == "status":
        print("Active SLURM jobs:", active_jobs(workdir))
        return

    check_files = not getattr(args, "allow_missing_inputs", False)
    units, samples, groups = validate(cfg, check_files=check_files)
    tasks, sample_maps = make_tasks(cfg, workdir, units, samples, groups)

    pending = pending_tasks(tasks, workdir, force=getattr(args, "force", False))

    print(f"Samples: {len(samples)} | Species: {', '.join(groups)}")

    for stage in STAGES:
        stage_tasks = [t for t in tasks if t["stage"] == stage]
        remaining = sum(pending[t["id"]] for t in stage_tasks)
        print(f"{stage:16} {remaining:5} pending / {len(stage_tasks):5} total")

    if args.action in ("validate", "plan"):
        print("No work was started and no output files were written.")
        return

    running = active_jobs(workdir)
    if running:
        raise RuntimeError(f"SLURM jobs are still running: {running}")

    if not any(pending.values()):
        print("All tasks have completed; nothing to run.")
        return

    plan = write_plan(tasks, sample_maps, config_file, workdir)

    if args.action == "local":
        for task in tasks:
            if pending[task["id"]]:
                run_one_task(plan, task["id"])

    elif args.action == "submit":
        submit_to_slurm(tasks, pending, workdir, plan, cfg)

if __name__ == "__main__":
    try:
        main()
    except (
        ValueError,
        KeyError,
        RuntimeError,
        OSError,
        subprocess.CalledProcessError,
    ) as exc:
        sys.exit(f"ERROR: {exc}")
