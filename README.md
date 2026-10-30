# Minimal SNP Panels — новая реализация v0.1

Пять модулей: `wgs` → `joint` → `qc` → `structure` → `panels`. Основные пользовательские входы: FASTQ с TSV либо готовый индексированный multi-sample VCF.gz/BCF. gVCF используется только между WGS и joint.

**Статус:** проектный код с офлайн-тестами. Никакого реального end-to-end теста на `sphinx` не было. Не заменяйте им работающую версию и не удаляйте существующие результаты.

## Установка на sphinx

```bash
conda activate wgs-gatk-mem3
# В исходном conda list присутствуют bwa-mem3 0.15.1, samtools/bcftools 1.24, GATK 4.7.0.0, fastp, mosdepth.
# Для аналитических модулей дополнительно требуются PLINK 2, scikit-learn, matplotlib.
conda install -n wgs-gatk-mem3 -c conda-forge -c bioconda plink2 scikit-learn matplotlib
python -m pip install --no-deps -e .
./run_pipeline.sh check-tools
```

Не устанавливайте зависимости во время обработки WGS; перед запуском проверьте, что окружение разрешает все зависимости. Изменение существующего Conda-окружения выполняется вручную. Можно создать его копию для проверки совместимости.

### Референсы

`inputs/` содержит пустые TSV-шаблоны **только с заголовками**, без выдуманных данных.

`config/references.yaml` содержит FASTA и списки аутосом белухи и белого медведя, перенесённые из исходного WGS-архива. Проверьте пути, `.fai`, `.dict`, индексы BWA-MEM3 до запуска:

```bash
samtools faidx /path/reference.fna
gatk CreateSequenceDictionary -R /path/reference.fna
bwa-mem3 index /path/reference.fna
```

**Осторожно:** `bwa-mem3 index` выполните только если индексы ещё не созданы. Не смешивайте индексы разных сборок.

## Команды

```bash
# 1. FASTQ → joint BCF → QC → PCA/clustering → SNP panels
./run_pipeline.sh run --fastq /path/fastq.tsv --species beluga --run-id beluga_001

# 2. Внешний VCF/BCF → QC → PCA/clustering → SNP panels
./run_pipeline.sh run --vcf /path/cohort.vcf.gz --species beluga --run-id beluga_vcf_001

# 3. Внешний VCF/BCF + метки → QC → панели, без PCA
./run_pipeline.sh run --vcf /path/cohort.bcf --labels /path/labels.tsv --species beluga --run-id beluga_panel_001

# 4. Только анализ структуры (без автоматического отбора SNP)
./run_pipeline.sh run --vcf /path/cohort.bcf --mode structure --species beluga --run-id beluga_structure_001

# 5. Только WGS (до per-contig gVCF)
./run_pipeline.sh wgs --fastq /path/fastq.tsv --species beluga --run-id beluga_wgs_001

# 6. Отдельный joint по уже созданному gvcf_manifest.tsv
./run_pipeline.sh joint --gvcf-manifest /path/gvcf_manifest.tsv --species beluga --run-id beluga_joint_001
```

Результаты сохраняются в `results/<species>/<run_id>/`, журналы в `<run_id>/logs/`. При повторном запуске неизменного конфигурационного набора используйте `--resume`. Повторный вызов на том же `run_id` без `--resume` завершится ошибкой.

### SLURM

Пока один `sbatch` на **весь вызов CLI** (последовательный анализ внутри job), без массивов задач или разбиения по образцам. Это базовая совместимость со SLURM, не масштабируемый HPC-оркестратор.

```bash
./run_pipeline.sh run --fastq /path/fastq.tsv --species beluga \
  --run-id beluga_slurm_001 --backend slurm
squeue -u "$USER"
```

ID задачи: `results/beluga/beluga_slurm_001/slurm_job_id.txt`. stdout: `logs/slurm_<job_id>.out`. Настраиваем `config/pipeline.yaml -> slurm` (CPU, RAM, partition, walltime). Job запускает тот же Python interpreter из активированного окружения, поэтому SLURM должен экспортировать PATH и переменные окружения.

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
