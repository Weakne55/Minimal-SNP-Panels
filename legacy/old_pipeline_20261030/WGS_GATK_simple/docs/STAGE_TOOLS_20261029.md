# WGS stages: проверка и обновление инструментов (29.10.2026)

Входные аргументы и научные команды `scripts/stages/*.sh` сохранены, кроме
создания всех каталогов вывода в `trim.sh`. Версии обновляются в `envs/wgs.yaml`,
а не вставляются в каждую Bash-команду. Пользовательская точка запуска и
SLURM-планировщик не изменены. Выполнение WGS не запускалось.

| Стадия | Команды | Решение | Основание и ограничение |
|---|---|---|---|
| `trim.sh` | fastp | 0.24.0 → 1.4.0 | Актуальный Bioconda; сохранить `--detect_adapter_for_pe`, выходы FASTQ/JSON/HTML; возможны изменения trimming |
| `align_lane.sh` | bwa-mem3, samtools collate/fixmate/sort | BWA-MEM3 0.15.1; samtools 1.24 | У BWA-MEM3 в Bioconda зависимость `htslib>=1.24,<1.25`; сохраняем текущий порядок tag/coordinate sort |
| `markdup.sh` | samtools merge/markdup/index | samtools 1.24; не менять | FastDup перспективен, но может менять определение дубликатов и JSON-отчёт; не drop-in replacement |
| `bam_qc.sh` | samtools flagstat/stats, mosdepth | mosdepth 0.3.11 → 0.3.14 | Отчёты и критерии `--mapq`, `--thresholds` без изменений |
| `import_aligned.sh` | samtools quickcheck/idxstats/view/index | samtools 1.24 | BAM анализируется без копирования; CRAM конвертируется в BAM без повторного markdup |
| `haplotypecaller.sh` | GATK HaplotypeCaller GVCF/IndexFeatureFile | GATK 4.6.2.0 → 4.7.0.0 | Обновление библиотеки I/O; возможны изменения результатов генотипирования |
| `genomicsdb.sh` | GATK GenomicsDBImport | GATK 4.7.0.0 | Сохраняем параметры; `--bypass-feature-reader` перспективен только после проверки валидности всех gVCF |
| `genotype.sh` | GATK GenotypeGVCFs/IndexFeatureFile | GATK 4.7.0.0 | Сохраняем GVCF-совместимый joint genotyping |
| `concat.sh` | bcftools view/concat/index | bcftools 1.24 | Версия 1.25 вышла 27.10.2026, но ABI и HTSlib конфликтуют с BWA-MEM3 0.15.1 |
| `multiqc.sh` | multiqc | 1.30 → 1.35 | Отчёт HTML может визуально измениться |

## Блокирующая зависимость HTSlib

В актуальной сборке Bioconda BWA-MEM3 0.15.1 опубликована зависимость
`htslib>=1.24,<1.25`; добавлять samtools/bcftools 1.25 в один YAML без
проверки solver нельзя. Объявление SAMtools/BCFtools/HTSlib 1.25
(27.10.2026) содержит исправления безопасности при обработке повреждённых
VCF/BCF/индексов. Ветка с HTSlib 1.24 в Bioconda требует дополнительного
решения: совместимая пересборка BWA-MEM3, разделение окружений по этапам
либо доступная исправленная сборка 1.24.1. **Не обрабатывать недоверенные
файлы старыми библиотеками, пока конфликт не разрешён.**

## Почему не заменяем научные инструменты сразу

- `FastDup`: производительность в опубликованной статье сравнивается в
  частности с Picard и Sambamba, а не обязательно с текущей командой
  `samtools markdup`; вход/выход и воспроизводимость надо проверять отдельно.
- `DeepVariant+GLnexus`, `Strelka2`, `bcftools mpileup/call`: не являются
  прозрачной заменой GATK HaplotypeCaller → GenomicsDB → GenotypeGVCFs.
  Меняют модель variant calling, формат промежуточных данных или процедуру
  joint genotyping; для белухи и медведя нужна самостоятельная валидация.
- `samtools sort`, `markdup`, mosdepth: уже потоковые инструменты с низкими
  накладными расходами; модификации числа потоков, памяти и I/O требуют
  измерений на кластере.

## Применение в отдельной среде

```bash
conda env create -f WGS_GATK_simple/envs/wgs.yaml --dry-run
conda env create -f WGS_GATK_simple/envs/wgs.yaml
conda activate wgs-gatk-mem3
bash WGS_GATK_simple/scripts/check_stack.sh
python -m pytest -q WGS_GATK_simple/tests tests
```

Команды создания окружения зависят от доступности пакетов и кэша Conda:
**разрешение зависимостей не проверено** в текущем автономном контейнере.
**Внимание к resume:** существующие сигнатуры задач хешируют команды и Bash-файлы,
но **не версии исполняемых программ, установленных в Conda**. Обновление только
`envs/wgs.yaml` не гарантирует пересчёт ранее завершённых HaplotypeCaller,
GenomicsDBImport или GenotypeGVCFs. Запускать обновлённый WGS необходимо
с новым `workdir`; никогда не подменять пакеты внутри старого активного
окружения/расчёта. Отдельно запланировано добавление toolchain manifest в
механизм подписей при рефакторинге управляющего кода.

Все версии должны быть записаны в отчёте каждого реального запуска.

## Источники

- https://github.com/fg-labs/bwa-mem3
- https://bioconda.github.io/recipes/bwa-mem3/README.html
- https://github.com/broadinstitute/gatk/releases
- https://bioconda.github.io/recipes/gatk4/README.html
- https://github.com/OpenGene/fastp
- https://anaconda.org/bioconda/fastp/files
- https://github.com/brentp/mosdepth/releases
- https://github.com/MultiQC/MultiQC/releases
- https://sourceforge.net/p/samtools/mailman/message/59422401/
- https://gatk.broadinstitute.org/hc/en-us/articles/4413056242203-GenomicsDBImport
