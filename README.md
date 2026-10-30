# Minimal SNP Panels — Beluga & Polar Bear

Пять модулей: `wgs` → `joint` → `qc` → `structure` → `panels`. Основные пользовательские входы: FASTQ с TSV либо готовый индексированный multi-sample VCF.gz/BCF. gVCF используется только между WGS и joint.

**Статус:** код и локальные тесты актуализированы; реальный end-to-end FASTQ → VCF на `sphinx` пока не подтверждён.

## Установка и проверка на sphinx

Основной Git-репозиторий: `Minimal-SNP-Panels` (ветка разработки `develop`).
Текущее окружение — `wgs-gatk-mem3` с BWA-MEM3, samtools, bcftools,
GATK4, fastp, mosdepth, PLINK2 и Python 3.11. Не обновляйте инструменты
посреди активного WGS-запуска.

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
conda activate wgs-gatk-mem3
python -m pip install --no-deps -e .
./run_pipeline.sh check-tools
python -m pytest -q
bash -n run_pipeline.sh stages/*.sh
```

Исходные FASTQ находятся в `data/wgs/`; исходная инвентаризация —
`data/metadata/WGS_final.xlsx`. Результаты прежних запусков сохраняются
в `results/`, старый код — в `legacy/`. Отсутствующие скрипты подготовки
данных и референсов искать в `tools/`.

### Референсы

Список аутосом и пути FASTA задаёт `config/references.yaml`:

- Белуха: `GCA_059928025.1` (21 аутосома).
- Белый медведь: `GCA_062004589.1` (36 аутосом).

До запуска FASTQ → VCF нужны FASTA, `.fai`, GATK sequence dictionary и
полный индекс BWA-MEM3, включая `<FASTA>.bwt.2bit.64`. **Индекс BWA
не заменяет индекс BWA-MEM3.** Предыдущая проверка на `sphinx` обнаружила
отсутствие `.bwt.2bit.64` для обоих видов; повторно проверяйте его перед
запуском. Подготовка индексов выполняется отдельным заданием с достаточной RAM.

Пустой `inputs/fastq.template.tsv` — **только шаблон**, а не готовый manifest.
Один sample может иметь несколько строк для разных sequencing runs.
Подробности формата — в [`docs/CONTRACTS.md`](docs/CONTRACTS.md).

## Команды запуска

Все примеры выполняются из корня репозитория. Файл `inputs/fastq.tsv`
нужно **сначала создать** из фактических FASTQ; его нет среди шаблонов.
Вместо `beluga` допустимо `polar_bear`.

```bash
# 1. Только FASTQ → BAM → per-contig gVCF и gvcf_manifest.tsv
./run_pipeline.sh wgs --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_gvcf_001 --backend slurm --threads 8 --memory 64G

# 2. FASTQ → BAM → gVCF → joint BCF → joint VCF.gz (без variant QC и панелей)
./run_pipeline.sh wgs-vcf --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_vcf_001 --backend slurm --threads 16 --memory 96G

# 3. Только совместное генотипирование из существующего gVCF manifest
./run_pipeline.sh joint --gvcf-manifest /path/to/gvcf_manifest.tsv \
  --species beluga --run-id beluga_joint_001 --backend slurm \
  --threads 8 --memory 64G

# 4. Только QC индексированного multi-sample VCF/BCF
./run_pipeline.sh qc --vcf /path/to/cohort.raw.bcf --species beluga \
  --run-id beluga_qc_001 --backend slurm --threads 4 --memory 32G

# 5. Полный анализ по FASTQ (WGS, joint, QC, structure, panels)
./run_pipeline.sh run --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_full_001 --backend slurm --threads 16 --memory 96G

# 6. Готовый VCF + внешние population labels → QC и панели
./run_pipeline.sh run --vcf /path/to/cohort.bcf \
  --labels /path/to/population_labels.tsv --species beluga \
  --run-id beluga_panels_001 --backend slurm --threads 8 --memory 64G

# 7. Локальная отладка с одним CPU (без отправки в SLURM)
./run_pipeline.sh wgs-vcf --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_local_001 --backend local --threads 1 --memory 48G
```

**Параметры:** `--threads` — целое число не меньше 1; `--memory` — бюджет RAM,
например `32G`, `65536M`. Они независимы. По умолчанию используются
`threads` и `slurm.memory` из `config/pipeline.yaml`. Внутри задания
учитывается фактическое выделение SLURM. Дополнительные CPU не обязаны
ускорять все инструменты. Выделение 1 CPU не отменяет требования BWA-MEM3
к достаточному объёму RAM. Подробнее —
[`docs/WGS_SLURM_RESOURCES.md`](docs/WGS_SLURM_RESOURCES.md).

**SLURM:** один `sbatch` на весь вызов (без job arrays); число CPU/RAM
задаётся `--threads`/`--memory`, partition и time — в `config/pipeline.yaml`.

```bash
squeue -u "$USER"
# Подставить полученный job ID
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

Результаты: `results/<species>/<run-id>/`; SLURM job ID:
`results/<species>/<run-id>/slurm_job_id.txt`; логи в `logs/`.
Для продолжения **того же запуска с теми же входами и ресурсами** применяется
`--resume`; для новых CPU/RAM, manifest или референса используйте новый `run-id`.
Не запускайте WGS по пустому TSV-шаблону.

### Основные выходы `wgs-vcf`

```text
results/beluga/beluga_vcf_001/
├── wgs/
│   ├── lanes/                # выравнивания отдельных runs
│   ├── bam/                  # объединённые BAM и BAM после markdup
│   ├── qc/                   # flagstat/stats/mosdepth
│   ├── gvcf/                 # per-sample × per-contig gVCF
│   └── gvcf_manifest.tsv
├── joint/
│   ├── genomicsdb/
│   ├── per_contig/
│   ├── cohort.raw.bcf
│   └── cohort.raw.vcf.gz
└── logs/
```

`wgs-vcf` **не** запускает SNP QC, population structure или поиск панелей.

## Научные ограничения

* `structure`: PCA по LD-pruned SNP + KMeans и silhouette; это **предварительные** генетические группы. Нужны контроль родства, batch-effects, ancestry/альтернативные методы и ручная проверка перед интерпретацией меток.
* `panels`: используется Nei-подобная мера дифференциации `(Ht-Hs)/Ht` (не Weir–Cockerham FST), плюс LD filtering по `r²`. Отбор SNP и LD рассчитываются отдельно внутри каждого training fold.
* При достаточном N сохраняется `nested_cv_metrics.tsv` с подбором размера на внутренних folds. При малых классах файл заменяется `nested_cv_unavailable.txt`.
* Выбор размера финальной панели по тем же CV без независимого теста не доказывает её внешнюю точность. Если labels получены кластеризацией того же WGS, метрики отражают воспроизведение кластеров, а не географическую валидацию.
* Нет импутации сложных low-pass genotype likelihoods, sex chromosome-specific ploidy, опубликованных SNP-панелей и независимого inference на новых образцах. Они требуют отдельной проверки перед финальной сдачей.

## Тестирование

```bash
python -m pytest -q
bash -n run_pipeline.sh stages/*.sh
python -m compileall -q src
```

В проекте **нет синтетических генотипов и искусственных биологических наборов данных**. Для реального smoke-test используйте сохранённые gVCF/VCF двух белух и реальные WGS, на отдельном `run_id`.

## Реальный smoke-test

Команды проверки подготовленного VCF белух, независимого joint и SLURM см. `docs/SPHINX_SMOKE_TEST.md`.
