#!/usr/bin/env python3
"""
Построение плана обработки полногеномных данных (WGS).

Этот модуль определяет:
    - какие этапы необходимо выполнить;
    - сколько отдельных задач нужно создать;
    - от каких предыдущих задач зависит каждая новая задача;
    - какие команды и выходные файлы соответствуют задачам.

Модуль сам не выполняет вычисления и не отправляет задания в SLURM.
Непосредственные команды BWA, samtools, GATK и остальных программ
находятся в scripts/stages/*.sh.

Запуском и возобновлением вычислений управляет scripts/pipeline.py.
"""
from __future__ import annotations

import hashlib
import re
import shlex
from pathlib import Path
from typing import TypedDict

import yaml

# Функции для проверки таблицы образцов и подготовленных референсов.
from model import abspath, check_references, load_manifest


# ============================================================
# 1. Константы и расположение проекта
# ============================================================

# Корневая директория проекта — родитель директории scripts/.
# Это позволяет запускать скрипт из любой текущей директории.
PROJECT = Path(__file__).resolve().parents[1]

# Порядок вычислительных этапов. Он используется при построении
# плана и при отправке групп задач через SLURM.
#
# ВАЖНО: это порядок стадий, а не самостоятельный список зависимостей.
# Зависимости конкретной задачи задаются отдельно через поле "deps".
STAGES = (
    "trim", "align", "markdup", "import", "qc",
    "haplotypecaller", "genomicsdb", "genotype", "concat", "multiqc",
)

# Соответствие: имя этапа -> Bash-скрипт с реальной командой.
# Изменять научные команды нужно в scripts/stages/, а не здесь.
SCRIPT = {
    "trim": "trim.sh",
    "align": "align_lane.sh",
    "markdup": "markdup.sh",
    "import": "import_aligned.sh",
    "qc": "bam_qc.sh",
    "haplotypecaller": "haplotypecaller.sh",
    "genomicsdb": "genomicsdb.sh",
    "genotype": "genotype.sh",
    "concat": "concat.sh",
    "multiqc": "multiqc.sh",
}

# Допустимые символы для идентификаторов задач.
# Идентификаторы используются также в именах файлов и checkpoint.
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class Task(TypedDict):
    """Serialized task schema used by workflow, resume, and SLURM.

    Keep a plain dictionary: existing checkpoints and JSON plans depend
    on these keys and must remain compatible with the old implementation.
    """

    id: str
    stage: str
    deps: list[str]
    inputs: list[str]
    outputs: list[str]
    body: str
    signature: str


# ============================================================
# 2. Загрузка настроек из config.yaml
# ============================================================

def config_load(config_file):
    """
    Загрузить YAML-конфигурацию и привести файловые пути к абсолютным.

    Аргументы:
        config_file: путь к config.yaml (абсолютный или относительно проекта).

    Возвращает:
        path:    абсолютный путь к config.yaml;
        cfg:     словарь настроек;
        workdir: абсолютный путь к директории результатов.

    Функция проверяет базовую структуру настроек, но не создаёт
    рабочую директорию и не запускает вычисления.
    """

    # Относительный путь к конфигурации считается от корня проекта.
    path = Path(config_file).expanduser()
    if not path.is_absolute():
        path = PROJECT / path
    path = path.resolve()
    # PyYAML читает файл в обычные Python-словари и списки.
    cfg = yaml.safe_load(path.read_text(encoding="utf8"))

    # Проверяем обязательные разделы конфигурации.
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a YAML dictionary")
    if not isinstance(cfg.get("workdir"), str) or not cfg["workdir"].strip():
        raise ValueError("Set workdir in config/config.yaml")

    # Результаты должны располагаться ВНЕ репозитория.
    # Это защищает код от смешения с большими BAM/GVCF/VCF файлами.
    workdir = Path(abspath(cfg["workdir"], PROJECT))
    if workdir == Path("/") or workdir == PROJECT or PROJECT in workdir.parents:
        raise ValueError("workdir must be outside the project directory and cannot be /")
    if not isinstance(cfg.get("references"), dict):
        raise ValueError("Set references in config/config.yaml")

    # Приводим пути к таблице образцов и FASTA к абсолютным.
    # Проверка референсов добавит в каждый раздел также путь к .dict.
    cfg["samples"] = abspath(cfg["samples"], PROJECT)
    for reference in cfg["references"].values():
        reference["fasta"] = abspath(reference["fasta"], PROJECT)
    # На этом этапе проверяем лишь формат настроек и имён contigs,
    # а не существование файлов: полная проверка выполняется в validate().
    check_references(cfg["references"], PROJECT, check_files=False)

    return path, cfg, workdir


# ============================================================
# 3. Проверка входных данных и параметров
# ============================================================

def validate(cfg, check_files=True):
    """
    Проверить образцы, референсы и настройки до запуска вычислений.

    При check_files=True проверяются подготовленные индексы референсов,
    наличие входных FASTQ/BAM/CRAM и индексов готовых BAM/CRAM.

    Возвращает:
        units:   список отдельных FASTQ-дорожек или файлов BAM/CRAM;
        samples: словарь уникальных образцов;
        groups:  словарь вида {species: [sample_id, ...]}.
    """

    # Проверки файлов референсов и формата samples.tsv реализованы
    # в model.py; здесь используется их результат.
    check_references(cfg["references"], PROJECT, check_files=check_files)
    units, samples, groups = load_manifest(
        cfg["samples"], cfg["references"], PROJECT, check_files=check_files
    )

    # Опциональная обработка FASTQ должна включаться строго через
    # логическое true/false в YAML (а не через строку "true").
    preprocessing = cfg.get("preprocessing", {})
    if not isinstance(preprocessing, dict) or not isinstance(
        preprocessing.get("trim_fastq", False), bool
    ):
        raise ValueError("preprocessing.trim_fastq must be true or false")

    # Оптическая дистанция для маркировки дубликатов неотрицательна.
    duplicate_marking = cfg.get("duplicate_marking", {})
    if not isinstance(duplicate_marking, dict) or int(
        duplicate_marking.get("optical_distance", 0)
    ) < 0:
        raise ValueError("Invalid duplicate_marking config")

    # Числовые параметры GATK (память, потоки, batch size и др.)
    # должны быть положительными целыми числами.
    gatk = cfg.get("gatk", {})
    for key in (
        "hc_heap_gb", "db_heap_gb", "genotype_heap_gb",
        "db_batch_size", "db_readers", "hc_threads",
    ):
        value = gatk.get(key, 1)
        if not isinstance(value, int) or value < 1:
            raise ValueError(f"gatk.{key} must be a positive integer")

    # Для готовых BAM/CRAM требуется индекс, иначе переходить
    # к дальнейшей обработке такого файла нельзя.
    # Для FASTQ эта проверка не нужна.
    if check_files:
        for unit in units:
            if unit.input_type not in ("bam", "cram"):
                continue
            file = Path(unit.source)
            if unit.input_type == "bam":
                indexes = [Path(str(file) + ".bai"), file.with_suffix(".bai"),
                           Path(str(file) + ".csi")]
            else:
                indexes = [Path(str(file) + ".crai"), file.with_suffix(".crai")]
            if not any(p.is_file() and p.stat().st_size > 0 for p in indexes):
                raise ValueError(f"Missing index for input {file}")

    return units, samples, groups


# ============================================================
# 4. Описание одной вычислительной задачи
# ============================================================

def make_task(stage, name, arguments, outputs, deps=(), inputs=(), extra="") -> Task:
    """
    Создать описание одной задачи для последующего выполнения.

    Аргументы:
        stage:     название стадии из STAGES;
        name:      уникальный идентификатор задачи;
        arguments: позиционные аргументы Bash-скрипта стадии;
        outputs:   ожидаемые выходные файлы;
        deps:      ID задач, от которых зависит текущая;
        inputs:    внешние входные файлы для проверки актуальности;
        extra:     дополнительные данные для контрольной подписи.

    Возвращает словарь с командой, зависимостями, путями и SHA-256.

    Здесь формируется ТОЛЬКО текст будущей команды — ничего
    не исполняется. Фактический запуск выполняет pipeline.py.
    """

    # Не допускаем специальные символы в ID: они попадут в пути
    # к рабочим файлам и в идентификаторы SLURM-задач.
    if not SAFE_ID.fullmatch(name):
        raise ValueError(f"Unsafe task ID: {name}")

    # shlex.join безопасно экранирует пробелы и спецсимволы
    # в аргументах команды Bash.
    command = shlex.join([
        "bash", str(PROJECT / "scripts" / "stages" / SCRIPT[stage]),
        *(str(arg) for arg in arguments),
    ])
    # Команда записывается как отдельный Bash-скрипт с
    # остановкой при ошибке любой команды.
    body = "set -euo pipefail\n" + command + "\n"

    # Сигнатура зависит от текста команды, содержимого скрипта
    # соответствующей стадии и дополнительного параметра extra.
    # Если код стадии/аргументы изменятся, checkpoint станет неактуален.
    digest = hashlib.sha256(
        body.encode()
        + (PROJECT / "scripts" / "stages" / SCRIPT[stage]).read_bytes()
        + str(extra).encode()
    )
    # Для импорта BAM/CRAM учитываем также код проверки заголовка.
    # Изменение validate_alignment.py потребует повторного импорта.
    if stage == "import":
        digest.update((PROJECT / "scripts" / "validate_alignment.py").read_bytes())
    # deps очищаем от повторов, сохраняя порядок.
    return {
        "id": name,
        "stage": stage,
        "deps": list(dict.fromkeys(deps)),
        "inputs": [str(x) for x in inputs],
        "outputs": [str(x) for x in outputs],
        "body": body,
        "signature": digest.hexdigest(),
    }


# ============================================================
# 5. Проверка запрошенных вычислительных ресурсов
# ============================================================

def check_resource_limits(cfg):
    """
    Проверить, что ресурсы SLURM достаточны для настроенных команд.

    Проверяются минимальное число CPU для alignment/markdup,
    совместимость потоков HaplotypeCaller с выделенными CPU
    и запас оперативной памяти сверх Java heap.

    Возвращает число CPU, выделяемых этапу align.
    """

    # Настройки можно переопределить в разделе slurm.resources.
    limits = cfg.get("slurm", {}).get("resources", {})
    gatk = cfg.get("gatk", {})

    # Значение ресурса для стадии либо пользовательское, либо default.
    def resource(stage, key, default):
        return int(limits.get(stage, {}).get(key, default))

    # Для текущих Bash-скриптов выравнивания и маркировки
    # дубликатов предусмотрено не менее 8 процессорных ядер.
    align_cpus = resource("align", "cpus", 8)
    if align_cpus < 8:
        raise ValueError("Alignment requires at least 8 SLURM CPUs")
    if resource("markdup", "cpus", 8) < 8:
        raise ValueError("Markdup requires at least 8 SLURM CPUs")
    if gatk.get("hc_threads", 4) > resource("haplotypecaller", "cpus", 4):
        raise ValueError("gatk.hc_threads exceeds SLURM haplotypecaller CPUs")

    # Системная RAM должна быть БОЛЬШЕ Java heap: Java и прочие
    # процессы расходуют память сверх размера кучи (-Xmx).
    for stage, setting, default_heap, default_mem in (
        ("haplotypecaller", "hc_heap_gb", 16, 24000),
        ("genomicsdb", "db_heap_gb", 24, 32000),
        ("genotype", "genotype_heap_gb", 12, 16000),
    ):
        heap = gatk.get(setting, default_heap)
        if resource(stage, "mem_mb", default_mem) <= heap * 1024:
            raise ValueError(f"{stage}: SLURM RAM must exceed Java heap")
    return align_cpus


# ============================================================
# 6. Построение полного графа задач (workflow)
# ============================================================

def _add_fastq_tasks(cfg, workdir, units, add, align_cpus):
    """Plan optional fastp and alignment of every FASTQ lane.

    Returns {sample_id: [(lane_bam, alignment_task_id), ...]}.
    """
    refs = cfg["references"]
    trim = cfg.get("preprocessing", {}).get("trim_fastq", False)
    # lane_bams хранит для каждого образца пары:
    #     (путь_к_BAM_дорожки, ID_задачи_alignment).
    # Эти пары понадобятся для объединения дорожек и markdup.
    lane_bams = {}
    for unit in units:
        # Готовые BAM/CRAM будут обработаны во втором блоке.
        if unit.input_type != "fastq":
            continue
        # Уникальный ID дорожки: образец + библиотека + lane.
        lane_id = f"{unit.sample}.{unit.library}.{unit.lane}"
        ref = refs[unit.species]["fasta"]
        r1, r2 = unit.r1, unit.r2
        dependencies = []
        inputs = [r1, r2]

        # Необязательный fastp: если включён trim_fastq,
        # сначала создаём очищенные FASTQ и QC-отчёты fastp.
        if trim:
            trimmed = workdir / "trimmed" / unit.sample
            r1 = str(trimmed / f"{unit.library}.{unit.lane}_R1.fastq.gz")
            r2 = str(trimmed / f"{unit.library}.{unit.lane}_R2.fastq.gz")
            qc_dir = workdir / "qc" / "fastp" / unit.sample
            json_report = qc_dir / f"{unit.library}.{unit.lane}.json"
            html_report = qc_dir / f"{unit.library}.{unit.lane}.html"
            trim_job = add(
                "trim", f"trim.{lane_id}",
                [unit.r1, unit.r2, r1, r2, json_report, html_report],
                [r1, r2, json_report, html_report],
                inputs=inputs,
            )
            # Alignment будет зависеть от успешного fastp.
            # Исходные FASTQ уже зарегистрированы как вход trim-задачи.
            dependencies = [trim_job]
            inputs = []  # trim job already watches source FASTQs

        # Создаём задачу картирования на референс и сортировки BAM.
        # Сигнатура и проверка актуальности учитывают FASTA и индексы BWA.
        lane_bam = workdir / "lanes" / unit.sample / f"{unit.library}.{unit.lane}.bam"
        align_job = add(
            "align", f"align.{lane_id}",
            [ref, r1, r2, lane_bam, unit.sample, unit.library, unit.lane, align_cpus],
            [lane_bam],
            deps=dependencies,
            inputs=[*inputs, ref, ref + ".fai", *(
                ref + "." + ext for ext in ("amb", "ann", "bwt", "pac", "sa")
            )],
        )
        # Запоминаем готовящийся BAM и ID задачи для объединения дорожек.
        lane_bams.setdefault(unit.sample, []).append((lane_bam, align_job))

    return lane_bams


def _add_sample_tasks(cfg, workdir, samples, lane_bams, add):
    """Plan BAM preparation/import and QC; return sample BAM and task IDs."""
    refs = cfg["references"]
    optical = int(cfg.get("duplicate_marking", {}).get("optical_distance", 0))
    # При входных FASTQ объединяем BAM дорожек и отмечаем дубликаты.
    # При входных BAM используем исходный файл без модификации.
    # При входных CRAM создаём BAM для дальнейшей обработки.
    # Значение: sample -> (путь_к_подготовленному_BAM, ID_задачи).
    prepared_samples = {}

    for sample, unit in samples.items():
        ref = refs[unit.species]["fasta"]
        if unit.input_type == "fastq":
            # У всех дорожек этого образца должна завершиться align,
            # только после этого начнётся объединение/markdup.
            bam = workdir / "bam" / f"{sample}.bam"
            metrics = workdir / "qc" / "duplicates" / f"{sample}.json"
            lane_files, lane_jobs = zip(*lane_bams[sample])
            prepare_job = add(
                "markdup", f"markdup.{sample}",
                [bam, metrics, optical, 8, *lane_files],
                [bam, str(bam) + ".bai", metrics],
                deps=lane_jobs,
            )
        elif unit.input_type == "bam":
            # Готовый BAM не перезаписываем; отдельная import-задача
            # проверяет его заголовок и создаёт отчёт idxstats.
            bam = unit.source  # validated external BAM is not modified
            report = workdir / "qc" / "input_validation" / f"{sample}.idxstats.tsv"
            prepare_job = add(
                "import", f"import.{sample}",
                ["bam", unit.source, report, sample, ref],
                [report],
                inputs=[unit.source, ref, ref + ".fai", refs[unit.species]["dict"]],
            )
        else:  # pre-existing CRAM -> BAM
            # Конвертируем CRAM в BAM с использованием выбранного
            # референса; исходный CRAM остаётся нетронутым.
            bam = workdir / "bam" / f"{sample}.bam"
            prepare_job = add(
                "import", f"import.{sample}",
                ["cram", unit.source, bam, sample, ref],
                [bam, str(bam) + ".bai"],
                inputs=[unit.source, ref, ref + ".fai", refs[unit.species]["dict"]],
            )

        # Дополнительно создаём QC-задачу для BAM.
        # Она зависит от подготовки BAM, но не от GATK.
        prepared_samples[sample] = (str(bam), prepare_job)
        qc_prefix = workdir / "qc" / "bam" / sample
        qc_outputs = [
            str(qc_prefix) + ext for ext in (
                ".flagstat.txt", ".stats.txt", ".mosdepth.summary.txt",
                ".regions.bed.gz", ".thresholds.bed.gz",
            )
        ]
        add("qc", f"qc.{sample}", [bam, qc_prefix, 4], qc_outputs,
            deps=[prepare_job])

    return prepared_samples


def _add_joint_calling_tasks(cfg, workdir, groups, prepared_samples, add):
    """Plan per-contig GVCFs, joint genotyping, and species VCFs.

    Sample-map contents are returned for the execution layer to write later.
    """
    refs = cfg["references"]
    gatk = cfg.get("gatk", {})
    hc_threads = gatk.get("hc_threads", 4)
    hc_heap = gatk.get("hc_heap_gb", 16)
    db_heap = gatk.get("db_heap_gb", 24)
    gt_heap = gatk.get("genotype_heap_gb", 12)
    db_batch = gatk.get("db_batch_size", 50)
    readers = gatk.get("db_readers", 1)
    # Виды обрабатываются независимо друг от друга.
    # Для каждого contig каждого образца создаётся отдельная GVCF.
    # После этого все GVCF данного contig объединяются в GenomicsDB
    # и генотипируются совместно.
    sample_maps = {}
    for species, sample_ids in sorted(groups.items()):
        ref = refs[species]["fasta"]

        for contig in refs[species]["contigs"]:
            # Соответствие sample -> GVCF и ID задач HaplotypeCaller.
            gvcfs = {}
            hc_jobs = []

            for sample in sample_ids:
                bam, prepare_job = prepared_samples[sample]
                gvcf = workdir / "gvcf" / species / sample / f"{contig}.g.vcf.gz"
                temp = workdir / "tmp" / "hc" / species / sample / contig
                # HaplotypeCaller не начнётся раньше подготовки BAM.
                hc_job = add(
                    "haplotypecaller", f"hc.{species}.{sample}.{contig}",
                    [ref, bam, contig, gvcf, hc_heap, hc_threads, temp],
                    [gvcf], deps=[prepare_job],
                    inputs=[ref, ref + ".fai", refs[species]["dict"]],
                )
                gvcfs[sample] = str(gvcf)
                hc_jobs.append(hc_job)

            # GATK GenomicsDBImport принимает sample-name map:
            # каждая строка — ID образца и путь к его GVCF.
            # Файл будет записан позднее при сохранении плана.
            map_content = "".join(f"{s}\t{gvcfs[s]}\n" for s in sorted(sample_ids))
            map_file = workdir / "jobs" / "sample_maps" / species / f"{contig}.tsv"
            sample_maps[str(map_file)] = map_content
            # GenomicsDBImport ждёт успешного завершения всех
            # HaplotypeCaller-задач для данного contig.
            db = workdir / "genomicsdb" / species / contig
            db_job = add(
                "genomicsdb", f"db.{species}.{contig}",
                [ref, contig, map_file, db, db_heap, db_batch, readers,
                 workdir / "tmp" / "genomicsdb" / species / contig],
                [db / ".completed"], deps=hc_jobs, extra=map_content,
            )
            # GenotypeGVCFs преобразует GenomicsDB в общий VCF.
            # Запуск зависит от готовности соответствующей GenomicsDB.
            vcf = workdir / "joint" / species / f"{contig}.vcf.gz"
            add(
                "genotype", f"genotype.{species}.{contig}",
                [ref, contig, db, vcf, gt_heap,
                 workdir / "tmp" / "genotype" / species / contig],
                [vcf], deps=[db_job],
            )

        # Собираем пути к VCF всех выбранных contigs данного вида.
        # Важен порядок contigs из config.yaml / индекса FASTA.
        chromosome_vcfs = [
            workdir / "joint" / species / f"{contig}.vcf.gz"
            for contig in refs[species]["contigs"]
        ]
        # Одна concat-задача объединяет VCF contigs в итоговый
        # cohort.raw.vcf.gz для данного вида.
        result = workdir / "vcf" / species / "cohort.raw.vcf.gz"
        add(
            "concat", f"concat.{species}", [result, *chromosome_vcfs],
            [result, str(result) + ".csi"],
            deps=[f"genotype.{species}.{contig}" for contig in refs[species]["contigs"]],
        )

    return sample_maps


def _add_multiqc_task(tasks, workdir, add):
    """Report BAM/fastp QC without introducing a dependency on variant calling."""
    # --------------------------------------------------------
    # 6.4. Сводный отчёт MultiQC и возврат плана
    # --------------------------------------------------------
    # MultiQC зависит от задач контроля качества и (если включён)
    # fastp, но НЕ от завершения совместного генотипирования.
    report = workdir / "qc" / "multiqc" / "multiqc_report.html"
    qc_jobs = [t["id"] for t in tasks if t["stage"] in ("trim", "qc")]
    add("multiqc", "multiqc", [workdir / "qc", report.parent], [report],
        deps=qc_jobs)



def make_tasks(cfg, workdir, units, samples, groups) -> tuple[list[Task], dict[str, str]]:
    """Build a complete, side-effect-free WGS task graph.

    Returns (tasks, sample_maps), using the existing task dictionaries and
    stable IDs, signatures, commands and stage order. Only pipeline.py writes
    the plan and submits tasks.
    """
    align_cpus = check_resource_limits(cfg)
    tasks = []
    task_ids = set()

    def add(stage, name, arguments, outputs, deps=(), inputs=(), extra=""):
        if name in task_ids:
            raise ValueError(f"Duplicate task ID: {name}")
        task_ids.add(name)
        tasks.append(make_task(stage, name, arguments, outputs, deps, inputs, extra))
        return name

    lane_bams = _add_fastq_tasks(cfg, workdir, units, add, align_cpus)
    prepared_samples = _add_sample_tasks(cfg, workdir, samples, lane_bams, add)
    sample_maps = _add_joint_calling_tasks(cfg, workdir, groups, prepared_samples, add)
    _add_multiqc_task(tasks, workdir, add)

    # Topological stage order is required by resume and SLURM submission.
    tasks.sort(key=lambda task: STAGES.index(task["stage"]))
    return tasks, sample_maps
