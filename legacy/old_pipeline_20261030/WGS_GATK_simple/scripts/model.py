"""
Чтение и проверка таблицы образцов и референсных геномов.

Каждая строка samples.tsv описывает:
    - одну дорожку секвенирования FASTQ (R1 и R2);
    - либо один готовый файл BAM/CRAM.

Если образец представлен несколькими дорожками FASTQ,
они объединяются после картирования.

Модуль также проверяет:
    - корректность названий образцов и видов;
    - наличие входных файлов;
    - отсутствие дублирующихся записей;
    - наличие референсных геномов и их индексов;
    - соответствие последовательностей в .fai и .dict.

Модуль не создаёт и не изменяет индексы референсных геномов.
"""

from __future__ import annotations

import csv
import re

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")

COLUMNS = (
    "sample",
    "species",
    "input_type",
    "library",
    "lane",
    "r1",
    "r2",
    "bam",
    "cram",
)

@dataclass(frozen=True)
class Unit:
    """
    Одна единица входных данных.

    Для FASTQ — одна дорожка секвенирования с парой R1/R2.
    Для BAM/CRAM — один готовый файл выровненных прочтений.

    Атрибуты:
        sample:     идентификатор образца;
        species:    название вида;
        input_type: тип данных (fastq, bam или cram);
        library:    идентификатор библиотеки;
        lane:       идентификатор дорожки;
        r1:         путь к FASTQ с первыми прочтениями;
        r2:         путь к FASTQ со вторыми прочтениями;
        bam:        путь к готовому BAM;
        cram:       путь к готовому CRAM.

    Неиспользуемые пути остаются пустыми строками.

    frozen=True запрещает изменение атрибутов после
    создания объекта.
    """

    sample: str
    species: str
    input_type: str
    library: str
    lane: str

    r1: str = ""
    r2: str = ""
    bam: str = ""
    cram: str = ""

    @property
    def key(self):
        """
        Уникальный ключ единицы входных данных.

        Позволяет различать несколько дорожек FASTQ
        одного образца.
        """
        return self.sample, self.library, self.lane

    @property
    def source(self):
        """
        Путь к готовому BAM или CRAM.

        Для FASTQ возвращает пустую строку.
        """
        return self.bam or self.cram

def abspath(filename: str, root: Path) -> str:
    """
    Преобразовать путь к абсолютному.

    Абсолютные пути остаются абсолютными.

    Относительные пути интерпретируются относительно
    корневой директории проекта (root), а не относительно
    текущей рабочей директории терминала.

    Символ ~ заменяется домашней директорией пользователя.

    Возвращает:
        Абсолютный путь в виде строки.
    """

    path = Path(filename).expanduser()

    if path.is_absolute():
        return str(path.resolve())

    return str((root / path).resolve())

def load_manifest(tsv_path, references, root, check_files=True):
    """
    Прочитать и проверить таблицу образцов.

    Аргументы:
        tsv_path:    путь к samples.tsv;
        references:  словарь референсных геномов из config.yaml;
        root:        корневая директория проекта;
        check_files: проверять существование входных файлов.

    Возвращает три объекта:
        1. Список Unit — все дорожки FASTQ или BAM/CRAM.
        2. Словарь unique_samples — по одному Unit на образец.
        3. Словарь by_species — списки образцов по видам.

    Проверяет:
        - наличие таблицы и обязательных столбцов;
        - допустимость идентификаторов;
        - соответствие видов конфигурации;
        - правильность указания FASTQ/BAM/CRAM;
        - существование и непустоту входных файлов;
        - отсутствие дублирующихся дорожек;
        - согласованность данных одного образца;
        - отсутствие нескольких библиотек FASTQ у образца.
    """

    source = Path(tsv_path)

    if not source.is_file():
        raise ValueError(f"Missing sample table: {source}")

    with source.open(newline="", encoding="utf8") as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if (
            not reader.fieldnames
            or not set(COLUMNS).issubset(reader.fieldnames)
        ):
            raise ValueError(
                "samples.tsv requires columns: " + ", ".join(COLUMNS)
            )

        if "coverage_type" in reader.fieldnames:
            raise ValueError(
                "Legacy coverage_type is not supported (WGS only)"
            )

        rows = list(reader)

    if not rows:
        raise ValueError(
            "samples.tsv contains only a header; add real samples"
        )

    units = {}

    unique_samples = {}

    libraries = defaultdict(set)

    for line_number, row in enumerate(rows, start=2):

        values = {
            name: (row.get(name) or "").strip()
            for name in COLUMNS
        }

        for field in ("sample", "species", "library", "lane"):

            if not SAFE_NAME.fullmatch(values[field]):
                raise ValueError(
                    f"Line {line_number}: "
                    f"invalid {field}: {values[field]!r}"
                )

        species = values["species"]

        if species not in references:
            raise ValueError(
                f"Line {line_number}: unknown species {species}"
            )

        kind = values["input_type"]

        required = {
            "fastq": {"r1", "r2"},
            "bam": {"bam"},
            "cram": {"cram"},
        }.get(kind)

        if required is None:
            raise ValueError(
                f"Line {line_number}: "
                "input_type must be fastq/bam/cram"
            )

        supplied = {
            field
            for field in ("r1", "r2", "bam", "cram")
            if values[field]
        }

        if supplied != required:
            raise ValueError(
                f"Line {line_number}: {kind} requires "
                f"{sorted(required)}, but found {sorted(supplied)}"
            )

        for field in required:

            filename = abspath(values[field], root)
            values[field] = filename

            if check_files:

                if (
                    not Path(filename).is_file()
                    or Path(filename).stat().st_size == 0
                ):
                    raise ValueError(
                        f"Line {line_number}: "
                        f"missing/empty input {filename}"
                    )

        unit = Unit(**values)

        if unit.key in units:
            raise ValueError(
                f"Line {line_number}: "
                f"duplicate sample/library/lane {unit.key}"
            )

        if unit.sample in unique_samples:

            previous = unique_samples[unit.sample]

            if (
                previous.species,
                previous.input_type
            ) != (
                unit.species,
                unit.input_type
            ):
                raise ValueError(
                    f"Line {line_number}: "
                    f"mixed species/input types for {unit.sample}"
                )

            if kind != "fastq":
                raise ValueError(
                    f"Line {line_number}: "
                    "only one BAM/CRAM per sample"
                )

        units[unit.key] = unit

        unique_samples[unit.sample] = unit

        if kind == "fastq":
            libraries[unit.sample].add(unit.library)

    for sample, names in libraries.items():

        if len(names) > 1:
            raise ValueError(
                f"{sample}: multiple libraries are not supported "
                "(duplicate marking must be library-aware)"
            )

    by_species = defaultdict(list)

    for sample, unit in unique_samples.items():
        by_species[unit.species].append(sample)

    return (
        list(units.values()),
        unique_samples,
        {
            name: sorted(sample_ids)
            for name, sample_ids in by_species.items()
        },
    )

def reference_dict_path(fasta: Path) -> Path:
    """
    Определить путь к словарю последовательностей .dict.

    GATK/Picard обычно используют словарь с тем же базовым
    именем, что и FASTA, но с расширением .dict.

    Пример:
        genome.fna -> genome.dict

    Функция только формирует путь.
    Файл не создаётся.
    """

    for extension in (
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".fasta",
        ".fna",
        ".fa",
    ):

        if fasta.name.endswith(extension):

            return fasta.with_name(
                fasta.name[:-len(extension)] + ".dict"
            )

    return Path(str(fasta) + ".dict")

def check_references(references, root, check_files=True, alignment_species=None):
    """
    Проверить референсные геномы и подготовленные индексы.

    Аргументы:
        references:  словарь референсов из config.yaml;
        root:        корневая директория проекта;
        check_files: проверять существование и содержимое файлов.
        alignment_species: множество видов, для которых нужен индекс
                           BWA-MEM3. None означает все виды.

    Проверяет:
        - наличие хотя бы одного референса;
        - допустимость названий видов и contigs;
        - отсутствие повторяющихся contigs;
        - наличие FASTA, FAI, DICT и при необходимости индексов BWA-MEM3;
        - наличие выбранных contigs в .fai;
        - правильный порядок contigs;
        - совпадение имён и длин последовательностей
          между .fai и .dict.

    Функция добавляет в конфигурацию абсолютные пути
    entry["fasta"] и entry["dict"].

    Не создаёт индексы и не изменяет файлы референса.
    """

    if not isinstance(references, dict) or not references:
        raise ValueError(
            "references must list at least one species"
        )

    for species, entry in references.items():

        if (
            not SAFE_NAME.fullmatch(species)
            or not isinstance(entry, dict)
        ):
            raise ValueError(
                f"Invalid reference specification for {species}"
            )

        contigs = entry.get("contigs")

        if (
            not isinstance(contigs, list)
            or not contigs
            or any(
                not isinstance(name, str)
                or not SAFE_NAME.fullmatch(name)
                for name in contigs
            )
        ):
            raise ValueError(
                f"{species}: contigs must be a list of valid names"
            )

        if len(contigs) != len(set(contigs)):
            raise ValueError(
                f"{species}: contigs repeated"
            )

        fasta = Path(
            abspath(str(entry.get("fasta", "")), root)
        )

        entry["fasta"] = str(fasta)

        dictionary = reference_dict_path(fasta)

        entry["dict"] = str(dictionary)

        if not check_files:
            continue

        files = [
            fasta,
            Path(str(fasta) + ".fai"),
            dictionary,
        ]

        if alignment_species is None or species in alignment_species:
            files += [
                Path(str(fasta) + "." + ext)
                for ext in ("amb", "ann", "bwt.2bit.64", "pac")
            ]

        for filename in files:

            if (
                not filename.is_file()
                or filename.stat().st_size == 0
            ):
                raise ValueError(
                    f"Missing reference/index {filename}"
                )

        fai_rows = [
            line.split("\t")
            for line in Path(
                str(fasta) + ".fai"
            ).read_text().splitlines()
            if line.strip()
        ]

        sequences = [
            (cols[0], int(cols[1]))
            for cols in fai_rows
        ]

        positions = {
            name: i
            for i, (name, length) in enumerate(sequences)
        }

        if len(positions) != len(sequences):
            raise ValueError(
                f"{fasta}: repeated contig names"
            )

        if any(name not in positions for name in contigs):
            raise ValueError(
                f"{species}: selected contigs missing from {fasta}.fai"
            )

        order = [
            positions[name]
            for name in contigs
        ]

        if order != sorted(order):
            raise ValueError(
                f"{species}: contig order must match FASTA index"
            )

        dict_lengths = {}

        for line in dictionary.read_text().splitlines():

            if not line.startswith("@SQ\t"):
                continue

            tags = dict(
                field.split(":", 1)
                for field in line.split("\t")[1:]
                if ":" in field
            )

            dict_lengths[tags.get("SN", "")] = int(
                tags.get("LN", "0")
            )

        if dict_lengths != dict(sequences):
            raise ValueError(
                f"{species}: .dict and .fai disagree "
                "on reference contigs"
            )
