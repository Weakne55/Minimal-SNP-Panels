# Minimal SNP Panels — Beluga & Polar Bear

Воспроизводимый конвейер для получения совместных VCF/BCF из WGS и отбора
минимальных наборов SNP для классификации белухи (*Delphinapterus leucas*) и
белого медведя (*Ursus maritimus*) по заданным популяционным меткам.

## Режимы и схема

| Режим | Команда | Результат |
|---|---|---|
| FASTQ → VCF | `--fastq MANIFEST.tsv` | BAM/gVCF для образцов, общий BCF и VCF на вид |
| VCF → structure | `--vcf COHORT.bcf` | Заглушка |
| VCF + metadata → SNP panel | `--vcf COHORT.bcf --metadata LABELS.tsv --species beluga` | Заглушка |

```text
FASTQ R1/R2
  → [trim] → align                         по sequencing unit (run/library)
  → merge_bam → markdup → bam_qc           по образцу
  → HaplotypeCaller                       один gVCF на образец
  → GenomicsDBImport → GenotypeGVCFs       один workspace/callset на вид
  → export_vcf                            один VCF.gz на вид

VCF/BCF → population structure            заглушка, без фиктивных результатов
VCF/BCF + metadata → QC → SNP ranking → LD selection → CV → minimal panel
```

## Быстрый старт на HPC

### 1. Получить код и установить окружение

```bash
git clone --branch develop --single-branch https://github.com/Weakne55/Minimal-SNP-Panels.git
cd Minimal-SNP-Panels
conda env create -f envs/wgs-gatk-mem3.yaml
conda activate wgs-gatk-mem3
python pipeline.py --help
python pipeline.py --check-tools
python -m pytest -q tests/
```

Потребуются Linux, Python 3.11, GATK4, BWA-MEM3 (AVX2), samtools, bcftools,
fastp, mosdepth и, для кластерного запуска, SLURM (`sbatch`, `squeue`).
Окружение задано в `envs/wgs-gatk-mem3.yaml`. Опция `--check-tools` печатает
версии, но не останавливается при отсутствии части программ: перед запуском
проверьте их доступность в `PATH` на вычислительных узлах.

### 2. Подготовить референсы

Используются следующие сборки NCBI:

| `species` | Genome accession | Число аутосом |
|---|---|---:|
| `beluga` | `GCA_059928025.1` | 21 |
| `polar_bear` | `GCA_062004589.1` | 36 |

Установка NCBI Datasets CLI и скачивание FASTA:

```bash
conda create -y -n ncbi-datasets -c conda-forge ncbi-datasets-cli

mkdir -p data/references/beluga data/references/polar_bear

conda run -n ncbi-datasets datasets download genome accession GCA_059928025.1 \
  --include genome,seq-report --filename data/references/beluga/ncbi_dataset.zip

conda run -n ncbi-datasets datasets download genome accession GCA_062004589.1 \
  --include genome,seq-report --filename data/references/polar_bear/ncbi_dataset.zip

unzip -n data/references/beluga/ncbi_dataset.zip -d data/references/beluga
unzip -n data/references/polar_bear/ncbi_dataset.zip -d data/references/polar_bear

find data/references -name '*_genomic.fna'
```

Для каждого полученного FASTA выполнить **на compute node**, а не на login node: 

```bash
conda activate wgs-gatk-mem3

samtools faidx /absolute/path/to/reference_genomic.fna

gatk CreateSequenceDictionary -R /absolute/path/to/reference_genomic.fna \
  -O /absolute/path/to/reference_genomic.dict

bwa-mem3 index /absolute/path/to/reference_genomic.fna
```

Повторите для обоих видов. Задайте реальные **абсолютные пути** `species.*.fasta`
в `config/references.yaml`. Существующие списки `contigs` в этом YAML не
меняйте: они определяют нужные аутосомы. Скрипт
`wgs_test/genome_reference_setup/setup_genomes.sh` использует старый BWA index
и **не заменяет** команду `bwa-mem3 index`.

Если индексы уже существуют и соответствуют этим сборкам, не создавайте их
повторно.

### 3. Подготовить FASTQ-манифест

Скачайте paired-end FASTQ из SRA/ENA или используйте локальные файлы. Для
публичных данных проверяйте соответствие run/BioSample, пары R1/R2 и MD5.
Полная процедура скачивания с ENA и проверкой MD5 приведена в
[руководстве по эксплуатации](docs/USER_GUIDE.md#4-скачать-fastq-и-создать-входной-tsv).

Файл `inputs/my_cohort.tsv` (разделитель — табуляция):

```tsv
sample_id	species	run_id	library	fastq_r1	fastq_r2
B01	beluga	RUN01	LIB01	/absolute/reads/B01_R1.fastq.gz	/absolute/reads/B01_R2.fastq.gz
B01	beluga	RUN02	LIB01	/absolute/reads/B01_R1_2.fastq.gz	/absolute/reads/B01_R2_2.fastq.gz
```

Одна строка — одна sequencing unit, несколько строк на образец допустимы.
Тройка `(sample_id, run_id, library)` должна быть уникальной, каждый FASTQ
используется один раз. Пример для редактирования: `inputs/fastq.template.tsv`.
Относительные FASTQ-пути разрешаются относительно расположения TSV.

### 4. Настроить SLURM

Измените в `config/pipeline.yaml` значения под свой кластер:

```yaml
trim_fastq: false
logging:
  level: INFO
slurm:
  partition: main
  nodelist: meduza-1  # удалите, если узел не нужно фиксировать
  max_concurrent_jobs: 3
  stages:
    align:           {cpus: 16, memory: 64G, time: '2-00:00:00'}
    haplotypecaller: {cpus: 8,  memory: 32G, time: '2-00:00:00'}
    genomicsdb:      {cpus: 8,  memory: 64G, time: '2-00:00:00'}
    genotype:        {cpus: 8,  memory: 48G, time: '2-00:00:00'}
```

Это **фрагмент**, не полный YAML: сохраните настройки остальных стадий.
Время, RAM и CPU задаются отдельно для каждого этапа; вся последовательность
берётся из `slurm.stages`. На стадиях joint genotyping одна задача обрабатывает
**все** указанные аутосомы вида — оцените её ресурсы заранее.

### 5. Проверить и запустить

```bash
# Без --submit: проверить манифест, референсы и распределение ресурсов.
python pipeline.py --fastq inputs/my_cohort.tsv \
  --backend slurm --run-id pilot_01

# Отправить задачи в SLURM.
python pipeline.py --fastq inputs/my_cohort.tsv \
  --backend slurm --run-id pilot_01 --submit

# Посмотреть прогресс по завершённым задачам.
python pipeline.py status --run-id pilot_01
```

**Важно:** `--backend local` — не dry run: без `--submit` он **сразу выполняет**
все стадии на текущей машине. Для возобновления локального запуска используйте
`--resume`. В SLURM повторная отправка arrays для существующего run-id
заблокирована; в случае частичной отправки сначала изучите `slurm_jobs.partial.json`
и `sacct`, затем примите решение об исправлении, отмене заданий или новом run-id.

### 6. Мониторинг и логи

```bash
python pipeline.py status --run-id pilot_01

squeue -u "$USER"

cat results/pilot_01/slurm_jobs.json

tail -f results/pilot_01/logs/pipeline.log

find results/pilot_01/logs -type f | sort
```

- `logs/pipeline.log` — CLI, планирование и отправка заданий.
- `logs/slurm/<stage>_<jobid>_<task>.out/.err` — SLURM stdout/stderr.
- `logs/workers/<stage>/<index>.log` — действия конкретного Python worker.
- `samples/<species>/<sample>/logs/*.log` и `joint/<species>/logs/*.log`
  — команды и объединённый stdout/stderr инструментов.
- `completed/<stage>/<index>.json` — маркер успешно проверенной стадии.

**Отсутствие задания в `squeue` не означает успех:** проверяйте `sacct`, маркер
`completed/` и целостность итоговых VCF/BCF.

## Результаты

| Путь внутри `results/<run-id>/` | Что содержит |
|---|---|
| `slurm_plan.json`, `versions.json` | План и версии доступных инструментов |
| `samples/<species>/<sample>/wgs/lanes/` | BAM по sequencing units |
| `samples/<species>/<sample>/wgs/bam/` | Объединённый BAM и BAM с отмеченными дубликатами |
| `samples/<species>/<sample>/wgs/qc/` | `flagstat`, `stats`, `mosdepth`, метрики дубликатов |
| `samples/<species>/<sample>/wgs/gvcf/` | Один индексированный gVCF образца по всем аутосомам |
| `joint/<species>/genomicsdb/` | Один GenomicsDB workspace вида |
| `joint/<species>/cohort.raw.bcf` | Совместные генотипы, бинарный VCF |
| `joint/<species>/cohort.raw.vcf.gz` | Совместные генотипы, сжатый VCF (с индексом) |

Значения `*.raw.*` — **ещё не отфильтрованные по качеству варианты**.
Для оценки качества проверьте sample IDs, покрытие, число вариантов и
совпадение `REF` с нужной сборкой.

## Панели и популяционная структура

```bash
# Пока возвращает NotImplementedError; файл labels не создаётся.
python pipeline.py --vcf results/pilot_01/joint/beluga/cohort.raw.bcf

# Реально запускает локальный QC и отбор панелей (НЕ SLURM).
python pipeline.py \
  --vcf results/pilot_01/joint/beluga/cohort.raw.bcf \
  --metadata inputs/my_population_labels.tsv \
  --species beluga --run-id beluga_panel_01
```

Для второго режима замените шаблон метаданных на фактический TSV с колонками
`sample_id`, `population` и опционально `group`. Все образцы VCF должны иметь
ровно по одной строке. Отбор SNP и LD выполняются **внутри training folds**.
Пути результатов: `qc/filtered.bcf`, `qc/qc_summary.tsv`,
`panels/minimal_panel.tsv`, `panels/cv_metrics.tsv`,
`panels/nested_cv_metrics.tsv` (если возможна), `panels/model.joblib`,
графики качества и матрицы ошибок. `model.joblib` загружайте только из
доверенного источника.

Полная интерпретация QC, CV, ограничений и восстановление после ошибок:
[руководство по эксплуатации](docs/USER_GUIDE.md).

## Для разработчиков

- [Архитектура и контракты](docs/PROJECT_GUIDE.md)
- [Подробная инструкция по установке и запуску](docs/USER_GUIDE.md)
- [Описание Python-модулей](src/minsnp/README.md)
- [Журнал рефакторинга](docs/REFACTORING.md)
- `python -m pytest -q tests/` — контрактные тесты без реального GATK/SLURM
- `python -m compileall -q pipeline.py src/minsnp` — синтаксис Python
- `for f in stages/*.sh; do bash -n "$f" || exit; done` — синтаксис Bash

Разработка и сравнение SNP-панелей — исследовательская задача; качество на
внешних данных и достаточность выборки необходимо проверять отдельно.
