# Minimal SNP Panels — Beluga & Polar Bear

**Reproducible WGS processing and minimal SNP-panel design for population assignment** in beluga whales (*Delphinapterus leucas*) and polar bears (*Ursus maritimus*).

Проект разрабатывает два независимых пути:

1. **Подготовка WGS (необязательная для пользователей готового VCF):** paired-end FASTQ → BAM → per-sample/per-autosome gVCF → joint genotyping → multi-sample BCF/VCF.gz.
2. **Основной анализ:** multi-sample VCF/BCF + `population_labels.tsv` → variant QC → SNP informativeness / LD filtering → candidate panels → classification / cross-validation → минимальная панель среди проверенных размеров.


## Устройство репозитория

```text
run_pipeline.sh             CLI (Python package src/minsnp/)
src/minsnp/                 wgs, joint, qc, structure, panels, reporting
stages/                     Bash: BWA-MEM3, samtools, GATK, bcftools
config/                     сборки, ресурсы и параметры (YAML)
inputs/                     TSV-шаблоны и примеры пилота
data/metadata/              инвентаризация WGS и её описание
data/wgs/                   локальные FASTQ, не публикуются в Git
results/                    локальные BAM, gVCF, BCF, логи, не публикуются
wgs_test/                   вспомогательные инструменты WGS и референсов
literature/                 обзор литературы и методов
AGENTS.md                   правила безопасной разработки
docs/PROJECT_GUIDE.md       техническая документация и история пилота
```

В уже запущенном пилоте SLURM-скрипты содержат путь **`tools/wgs_pilot.py`**. На сервере необходимо сохранить эту локальную совместимую копию `tools/` до завершения массива и joint jobs. Для будущих запусков из актуального репозитория предназначен `wgs_test/`. Каталог `legacy/` сохранён только локально на сервере.

## Установка

Требуется Linux/HPC, Python >=3.11, Conda, BWA-MEM3, samtools, bcftools, GATK4, fastp (если включён trimming), mosdepth; для структуры — PLINK2. Для запуска на SLURM нужен `sbatch`.

```bash
# На sphinx: текущее окружение, не обновлять во время вычислений
conda activate wgs-gatk-mem3
./run_pipeline.sh check-tools
./run_pipeline.sh --help
bash -n run_pipeline.sh stages/*.sh
```

Локальные тесты (если каталог `tests/` присутствует на сервере): `PYTHONPATH=src python -m pytest -q`. Тесты исключены из опубликованного Git-дерева по текущему решению о структуре репозитория.

### Референсы

`config/references.yaml` задаёт пути к FASTA и аутосомы каждого вида:

| Вид | Сборка | Аутосомы |
|---|---|---:|
| Белуха | `GCA_059928025.1` (BTC_Dleu_Yukka2) | 21 |
| Белый медведь | `GCA_062004589.1` (BTC_Umar_Aika3) | 36 |

Для BWA-MEM3 требуется, в частности, индекс `<FASTA>.bwt.2bit.64`; старый индекс BWA-MEM не является его заменой. BWA-MEM3 требует AVX2: на `sphinx` для индексации и картирования пилота использовался `meduza-1`, а не `kraken`.

## Интерфейс

Выполняйте команды из корня проекта. **Примеры ниже документируют CLI; не отправляйте повторную обработку пилота `pilot_20261030`.**

```bash
# FASTQ отдельной когорты → BAM + per-autosome gVCF
./run_pipeline.sh wgs --fastq inputs/beluga_pilot_10.fastq.tsv \
  --species beluga --run-id my_new_wgs --backend local \
  --threads 8 --memory 64G

# Уже полученные gVCF → многосэмпловый BCF
./run_pipeline.sh joint --gvcf-manifest /path/to/gvcf_manifest.tsv \
  --species beluga --run-id my_joint --backend local \
  --threads 8 --memory 64G

# Готовый multi-sample BCF → variant QC
./run_pipeline.sh qc --vcf /path/to/cohort.raw.bcf \
  --species beluga --run-id my_qc --backend local

# Основной анализ по VCF/BCF и заранее подготовленным меткам
./run_pipeline.sh run --vcf /path/to/cohort.raw.bcf \
  --labels /path/to/population_labels.tsv --mode panels \
  --species beluga --run-id my_panels --backend local
```

Команды с FASTQ приведены **только как примеры синтаксиса**, не как предложение запускать их на тех же данных. Для нового WGS потребуется отдельный `run-id`, подготовленные FASTQ и референсы; `--threads` для WGS — **не менее 2**. Для отдельного исследования структуры: `./run_pipeline.sh structure --vcf ... --species ... --run-id ...`, но **сначала необходимо исправить `structure.py` и проверить тестами**.

**Ограничения на TSV:** FASTQ — `sample_id`, `species`, `run_id`, `library`, `fastq_r1`, `fastq_r2`; labels — `sample_id`, `population` и необязательный `group`. Пути FASTQ относительно расположения TSV; шаблоны без строк не запускаются. Входной VCF/BCF должен быть индексирован, относиться к одному виду/сборке и содержать минимум два образца.