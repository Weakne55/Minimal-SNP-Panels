# Загрузка **полных** WGS для Minimal SNP Panels

Этот пакет предназначен для первой **полноразмерной** когорты WGS: две белухи и два белых медведя. Он не нарезает FASTQ, не использует тестовые интервалы и не изменяет исходный `WGS_Bash_GATK`.

Референсы: белуха `GCA_059928025.1` (Yukka2) и белый медведь `GCA_062004589.1` (Aika3). Готовые FASTA, `.fai`, `.dict`, BWA-индексы переиспользуются без изменений.

Данные: белухи — BioProject `PRJNA1328577` (полноразмерные Illumina WGS, обычно ≥20x); медведи — `SRR15170741` (46,3 Gb оснований) и `SRR15170742` (71,9 Gb оснований). Скачиваются **все чтения** выбранных ENA run accession; если у выбранного BioSample несколько ENA runs одной библиотеки, скачиваются **все эти runs** как отдельные ланы. Отбор белух выполняется по ENA metadata (больше 50 млрд прочтённых оснований суммарно, разные BioSample, одна библиотека; из подходящих выбираются меньшие по размеру сжатых FASTQ). Выбор не обеспечивает представительность по географическим популяциям — это технический первый набор WGS.

Перед большой загрузкой скрипт создаёт `ena_selection.json` с фиксированными URL, размерами и MD5. `download` не меняет выбор. Он использует `curl -C -` для возобновления незавершённых `.part`, сверяет полный MD5 **обоих FASTQ** и не принимает повреждённые файлы. Успешные загрузки не перезаписываются. Можно запускать снова после прерывания.

## 1. Распаковка, переменные

Разместите этот архив в `/mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels`, затем:

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
unzip full_wgs_setup.zip
conda activate wgs-gatk

BASE=/mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
PROJ="$BASE/WGS_Bash_GATK"
HELP="$BASE/full_wgs_setup"
REFROOT=/mnt/tank/scratch/plem/data/genomes
DATA="$BASE/data/full_wgs_20261022"
CFG="$BASE/config/full_wgs_20261022"
WORK="$BASE/results/full_wgs_20261022"
PLAN="$CFG/ena_selection.json"
```

## 2. Проверка полного набора аутосом по NCBI (до большой загрузки)

```bash
python "$HELP/prepare_full_wgs.py" references \
  --plan "$PLAN" --reference-root "$REFROOT"
```

Должно быть **21** и **36** аутосом. Если есть ошибка, остановитесь и сначала исправьте соответствие NCBI sequence report / FASTA. Скрипт ничего не изменяет в каталогах референсов.

## 3. Запросить актуальные ENA accession и оценить объём **без скачивания**

```bash
python "$HELP/prepare_full_wgs.py" discover --plan "$PLAN"
python "$HELP/prepare_full_wgs.py" report --plan "$PLAN"
df -h "$BASE" "$REFROOT"
```

Список accession — не результат гадания. Для белух конкретные SRR будут получены из ENA BioProject `PRJNA1328577`. Проверяем названия вида, Illumina, `library_strategy=WGS`, `library_layout=PAIRED`, уникальность BioSample, наличие двух FASTQ, byte size и MD5. Для медведей фиксированы два ранее подтверждённых SRR. **Проверьте выведенный список и сумму гигабайт** до запуска загрузки.

Скрипт автоматически не продолжит работу, если ENA не ответит, данные недоступны или не пройдут критерии. Для повторного выбора accession сохранённый PLAN следует явно переименовать, чтобы не потерять provenance.

## 4. Скачать **все** FASTQ

Длинную загрузку лучше запускать в допустимой для кластера долгоживущей сессии (например, `tmux` на download/login узле, если правилами это разрешено) или отдельным SLURM-заданием **только если вычислительные узлы имеют доступ в интернет**.

```bash
python "$HELP/prepare_full_wgs.py" download \
  --plan "$PLAN" \
  --project "$PROJ" \
  --reference-root "$REFROOT" \
  --reads "$DATA" \
  --output "$CFG" \
  --workdir "$WORK" \
  2>&1 | tee "$CFG/download.log"
```

Замечание: для длительной загрузки **не следует полагаться на один pipeline с `tee` без `pipefail`**. Безопаснее:

```bash
set -o pipefail
python "$HELP/prepare_full_wgs.py" download \
  --plan "$PLAN" --project "$PROJ" --reference-root "$REFROOT" \
  --reads "$DATA" --output "$CFG" --workdir "$WORK" \
  2>&1 | tee "$CFG/download.log"
```

Если `curl` оборвался — повторить `download` с теми же путями. Не удалять проверенные готовые FASTQ. При MD5 mismatch остановиться и выяснить причину; повреждённый `.part` можно удалить вручную после проверки.

Требуется свободное пространство не только под сжатые FASTQ, но и под BAM, gVCF, GenomicsDB и промежуточные файлы. Грубая оценка — несколько сотен гигабайт сверх размера FASTQ; проверить свободный объём до загрузки и до вычислений. Не использовать `/home` для больших WGS.

## 5. Полный список аутосом

`prepare_full_wgs.py` прочитает `sequence_report.jsonl`, сохранённый при скачивании NCBI Datasets (`--include genome,seq-report`). Учитываются только последовательности с `role=assembled-molecule`, `assignedMoleculeLocationType=Chromosome` и **числовым** `chrName`. Проверяются соответствие accessions/длин `.fai` и число аутосом: **21 для белухи, 36 для медведя**. Порядок соответствует `.fai`. X, Y, митохондрия и unplaced scaffolds в анализ не включаются, так как здесь предусмотрена стандартная диплоидная модель GATK.

Если метаданные chromosome names не соответствуют ожиданиям, скрипт **остановится**, а не выберет первые N хромосом по размеру. Тогда нужно изучить исходный `sequence_report.jsonl` и уточнить карту хромосом. Файлы референса не изменяются.

После всех загрузок появятся:

```text
config/full_wgs_20261022/
    ena_selection.json        # BioSample + ENA run + MD5 + размеры + источники
    full_samples.tsv          # Все FASTQ и строки лан
    full_wgs.yaml             # Полные аутосомы, FASTA, workdir
    autosome_selection.tsv   # Проверяемая таблица выбора аутосом

data/full_wgs_20261022/<species>/<run>/<run>_{1,2}.fastq.gz
```

## 6. Запуск настоящего WGS pipeline

```bash
cd "$PROJ"
python -m unittest discover -s tests -v
bash run_pipeline.sh validate --config "$CFG/full_wgs.yaml"
bash run_pipeline.sh plan --config "$CFG/full_wgs.yaml"
bash run_pipeline.sh submit --config "$CFG/full_wgs.yaml"
squeue -u "$USER"
```

Никаких `test_intervals`, subsampling, low-pass режима, фильтрации SNP или QC-gate здесь нет. Этапы: весь FASTQ → BWA-MEM → BAM, samtools markdup → **GATK HaplotypeCaller per (sample, autosome)** → gVCF → GenomicsDBImport → GenotypeGVCFs → VCF по каждому виду. BAM QC/MultiQC выполняется как отчёт, но не фильтрует.

Ожидаемые результаты:

```text
results/full_wgs_20261022/
  bam/<sample>.bam, .bam.bai
  gvcf/<species>/<sample>/<autosome>.g.vcf.gz
  genomicsdb/<species>/<autosome>/
  joint/<species>/<autosome>.vcf.gz
  vcf/Delphinapterus_leucas/cohort.raw.vcf.gz
  vcf/Ursus_maritimus/cohort.raw.vcf.gz
  qc/bam/, qc/multiqc/, jobs/logs/, checkpoints/
```

Проверить samples итоговых VCF:

```bash
bcftools query -l "$WORK/vcf/Delphinapterus_leucas/cohort.raw.vcf.gz"
bcftools query -l "$WORK/vcf/Ursus_maritimus/cohort.raw.vcf.gz"
bcftools index -n "$WORK/vcf/Delphinapterus_leucas/cohort.raw.vcf.gz"
bcftools index -n "$WORK/vcf/Ursus_maritimus/cohort.raw.vcf.gz"
```

В конце должно быть по 2 уникальных sample ID на вид. Если job завершился с ошибкой, сперва посмотрите `jobs/logs/*err` и `jobs/submission.json`; не запускайте второй submit пока первый массив ещё работает. После остановки неуспешных заданий можно использовать resume — обычный `submit` без `--force`.

## Условия и ограничения

- Это **настоящие полноценные WGS**, а не репрезентативная population panel. Первые четыре образца нужны для проверки работы инфраструктуры и полной цепочки GATK; потом выберем остальные WGS с учётом популяций и качества метаданных.
- У медведя `SRR15170741` ~46,3 Гб последовательностей (~19× номинального покрытия относительно сборки 2,415 Гб); `SRR15170742` ~71,9 Гб (~30×). Это отношение до фильтрации чтений/дубликатов, не оценка фактического depth.
- Для белух критерий ≥50 млрд оснований означает ориентировочно ≥17,6× относительно 2,847 Гб до фильтрации, но отдельные BioSample могут иметь выше.
- Полный запуск на вашем HPC здесь не проверялся: технические и биологические показатели оценим после завершения первых BAM/gVCF.
- Для полноценных научных панелей понадобится отдельный QC этап и выбор состава когорт. Это сознательно отложено по договорённости.
- Не копируйте референсные FASTA и не стройте индексы повторно.

Источники: NCBI BioProject PRJNA1328577; NCBI SRA SRR15170741/SRR15170742; ENA filereport API; NCBI Datasets sequence_report.
