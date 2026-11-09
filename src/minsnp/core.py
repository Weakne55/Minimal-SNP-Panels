"""Общие проверки, файловые операции и запуск внешних программ."""

from __future__ import annotations

import csv
import logging
import os
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def identifier(s: str) -> str:
    """Проверить идентификатор на безопасность использования в путях."""
    if not SAFE_ID.fullmatch(s) or s in {".", ".."}:
        raise ValueError(f"Некорректный идентификатор: {s!r}")
    return s


def read_tsv(path: str | Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    """Прочитать TSV с обязательными столбцами и непустыми значениями."""
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
    """Записать TSV, создав родительский каталог при необходимости."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def file_exists(path: str | Path) -> Path:
    """Проверить существование непустого файла и вернуть абсолютный путь."""
    p = Path(path).expanduser().resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise FileNotFoundError(f"Файл не найден или пуст: {p}")
    return p


def config(refs: Path, settings: Path, species: str) -> tuple[dict, dict]:
    """Загрузить YAML и проверить референс и список аутосом вида."""
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

    # В первом столбце .fai перечислены contig исходного FASTA.
    actual = {line.split("\t")[0] for line in fai.read_text().splitlines()}
    missing = set(spec["contigs"]) - actual
    if missing:
        raise ValueError(f"Contigs отсутствуют в FASTA: {sorted(missing)}")
    spec = dict(spec, fasta=str(fasta))
    return cfg, spec


def capture_versions(names: list[str]) -> dict[str, str]:
    """Собрать версии внешних инструментов для воспроизводимости запуска."""
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


def call(args: list[str], logdir: Path, name: str, env: dict[str, str] | None = None) -> None:
    """Выполнить команду, записать stdout/stderr в лог и проверить код возврата."""
    logdir.mkdir(parents=True, exist_ok=True)
    logfile = logdir / f"{name}.log"
    logger = logging.getLogger("minsnp")
    logger.info("START command stage=%s log=%s", name, logfile)
    started = time.monotonic()
    with logfile.open("a", encoding="utf-8") as log:
        log.write("\n$ " + shlex.join(map(str, args)) + "\n")
        log.flush()
        proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=False,
                              **({"env": env} if env is not None else {}))
    logger.log(logging.INFO if not proc.returncode else logging.ERROR,
               "END command stage=%s exit=%d elapsed_s=%.1f",
               name, proc.returncode, time.monotonic() - started)
    if proc.returncode:
        raise RuntimeError(f"Ошибка этапа {name}, код {proc.returncode}; лог: {logfile}")


def stage(script: str, args: list[str], logdir: Path, name: str, cfg: dict | None = None) -> None:
    """Выполнить Bash-стадию с лимитами ресурсов из конфигурации."""
    env = None
    if cfg is not None:
        from .resources import stage_effective
        cpus, memory = stage_effective(cfg, script)
        env = dict(os.environ, MIN_SNP_THREADS=str(cpus),
                   MIN_SNP_MEMORY_MB=str(memory), OMP_NUM_THREADS=str(cpus))
    call(["bash", str(ROOT / "stages" / script), *map(str, args)], logdir, name, env=env)


def require_commands(*commands: str) -> None:
    """Проверить наличие необходимых программ в PATH."""
    missing = [c for c in commands if shutil.which(c) is None]
    if missing:
        raise RuntimeError("Не найдены инструменты в PATH: " + ", ".join(missing))


def check_index(vcf: str | Path) -> None:
    """Проверить наличие .csi или .tbi (без проверки содержимого индекса)."""
    p = Path(vcf)
    if not (Path(str(p) + ".csi").is_file() or Path(str(p) + ".tbi").is_file()):
        raise FileNotFoundError(f"Отсутствует индекс VCF/BCF: {p}")


def samples_in_vcf(vcf: str | Path) -> list[str]:
    """Получить список образцов из VCF/BCF с помощью bcftools."""
    require_commands("bcftools")
    result = subprocess.run(["bcftools", "query", "-l", str(vcf)],
                            capture_output=True, text=True, check=True)
    return result.stdout.splitlines()
