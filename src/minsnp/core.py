"""Контракты, запуск программ и воспроизводимость."""
from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import yaml

ROOT = Path(__file__).resolve().parents[2]
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def identifier(s: str) -> str:
    if not SAFE_ID.fullmatch(s) or s in {".", ".."}:
        raise ValueError(f"Некорректный идентификатор: {s!r}")
    return s


def read_tsv(path: str | Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or not set(columns).issubset(reader.fieldnames):
            raise ValueError(f"{path}: ожидаются столбцы {', '.join(columns)}")
        rows = list(reader)
    if not rows or any(not all(row.get(key, "").strip() for key in columns) for row in rows):
        raise ValueError(f"{path}: пустая таблица или обязательные поля")
    return rows


def write_tsv(path: Path, columns: list[str], rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def file_exists(path: str | Path) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise FileNotFoundError(f"Файл не найден или пуст: {p}")
    return p


def config(refs: Path, settings: Path, species: str) -> tuple[dict, dict]:
    cfg = yaml.safe_load(file_exists(settings).read_text())
    refs_yaml = yaml.safe_load(file_exists(refs).read_text())
    try:
        spec = refs_yaml["species"][species]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Вид {species!r} отсутствует в references.yaml") from exc
    if not isinstance(spec.get("contigs"), list) or not spec["contigs"]:
        raise ValueError(f"В references.yaml не задан проверенный список contigs вида {species}")
    if len(spec["contigs"]) != len(set(spec["contigs"])):
        raise ValueError("Повторяющиеся contigs в references.yaml")
    fasta = file_exists(spec["fasta"])
    fai = file_exists(str(fasta) + ".fai")
    actual = {line.split("\t")[0] for line in fai.read_text().splitlines()}
    missing = set(spec["contigs"]) - actual
    if missing:
        raise ValueError(f"Contigs отсутствуют в FASTA: {sorted(missing)}")
    spec = dict(spec, fasta=str(fasta))
    return cfg, spec


def capture_versions(names: list[str]) -> dict[str, str]:
    versions = {}
    for name in names:
        exe = shutil.which(name)
        if exe is None:
            versions[name] = "НЕ УСТАНОВЛЕН"
            continue
        try:
            cmd = {"gatk": [name, "--version"], "bwa-mem3": [name, "version"]}.get(name, [name, "--version"])
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            versions[name] = (result.stdout + result.stderr).strip().splitlines()[0][:240]
        except (OSError, subprocess.TimeoutExpired, IndexError):
            versions[name] = str(exe)
    return versions


def call(args: list[str], logdir: Path, name: str) -> None:
    logdir.mkdir(parents=True, exist_ok=True)
    logfile = logdir / f"{name}.log"
    with logfile.open("a", encoding="utf-8") as log:
        log.write("\n$ " + " ".join(__import__("shlex").quote(x) for x in args) + "\n")
        log.flush()
        proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=False)
    if proc.returncode:
        raise RuntimeError(f"Ошибка этапа {name}, код {proc.returncode}; лог: {logfile}")


def stage(script: str, args: list[str], logdir: Path, name: str) -> None:
    call(["bash", str(ROOT / "stages" / script), *map(str, args)], logdir, name)


def source_fingerprint(path: str | Path | None) -> str | None:
    if path is None:
        return None
    p = file_exists(path)
    # Для крупных VCF/FASTQ не считаем полный SHA256: для них фиксируем размер и mtime.
    if p.stat().st_size > 10_000_000:
        stat = p.stat()
        return f"size={stat.st_size};mtime_ns={stat.st_mtime_ns}"
    return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare_run(base: Path, species: str, run_id: str, inputs: dict, resume: bool) -> Path:
    run = base / identifier(species) / identifier(run_id)
    run.mkdir(parents=True, exist_ok=True)
    state = run / "run.json"
    obj = {"species": species, "inputs": inputs}
    if state.exists():
        prev = json.loads(state.read_text())
        if prev != obj or not resume:
            raise ValueError(f"Каталог уже существует; для того же запуска используйте --resume: {run}")
    else:
        state.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    return run


def require_commands(*commands: str) -> None:
    missing = [c for c in commands if shutil.which(c) is None]
    if missing:
        raise RuntimeError("Не найдены инструменты в PATH: " + ", ".join(missing))


def check_index(vcf: str | Path) -> None:
    p = Path(vcf)
    if not (Path(str(p) + ".csi").is_file() or Path(str(p) + ".tbi").is_file()):
        raise FileNotFoundError(f"Отсутствует индекс VCF/BCF: {p}")


def samples_in_vcf(vcf: str | Path) -> list[str]:
    require_commands("bcftools")
    result = subprocess.run(["bcftools", "query", "-l", str(vcf)], capture_output=True, text=True, check=True)
    return result.stdout.splitlines()
