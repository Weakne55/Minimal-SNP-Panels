# BWA-MEM → BWA-MEM3 (29.10.2026)

## Что именно поменялось

- В `align_lane.sh`: `bwa mem` заменён на `bwa-mem3 mem --compat=bwa-mem`.
- Samtools-конвейер `collate | fixmate | sort` и параметры `-t`, `-R` не изменены.
- В проверке референса и входах задач `align`: индекс BWA v1 (`.bwt`, `.sa`) заменён на BWA-MEM3 (`.bwt.2bit.64`); `.amb`, `.ann`, `.pac` сохраняются.
- Входы BAM/CRAM не требуют индекса aligner (но по-прежнему требуют FASTA/FAI/DICT).
- В `envs/wgs.yaml` BWA заменён на BWA-MEM3 0.15.1, окружение названо `wgs-gatk-mem3`.
- **Из-за HTSlib-совместимости сборки Bioconda** samtools и bcftools в новом YAML изменены с 1.21 на 1.24. В `envs/wgs_bwa_legacy.yaml` сохранена прежняя конфигурация для контрольного сравнения. Это означает, что сравнение целого WGS-пайплайна — не строго одна переменная; для изоляции эффекта aligner сравнивать именно выравнивание при одинаковой версии samtools.

Публичные источники: https://github.com/fg-labs/bwa-mem3 ; https://bwa-mem3.readthedocs.io/en/latest/getting-started/migrating.html ; https://bioconda.github.io/recipes/bwa-mem3/README.html

## 1. Порядок развёртывания на sphinx

Работать в отдельной ветке и проверить `git status` до распаковки. Не перезаписывать работающий `pipeline_plan.json` и не удалять checkpoints. Остановить/дождаться завершения **своих** задач выравнивания перед переиндексацией общих референсов, если они используют те же индексы. Не запускать `./run.sh submit` до тестов.

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
# В своей ветке для миграции:
# tar -xzf /PATH/bwa_mem3_migration_patch_20261029.tar.gz

# Не меняет существующее окружение wgs-gatk
conda env create -f WGS_GATK_simple/envs/wgs.yaml
conda activate wgs-gatk-mem3
bwa-mem3 version
samtools --version | head -1
bcftools --version | head -1

python -m pytest -q WGS_GATK_simple/tests tests
bash WGS_GATK_simple/scripts/check_stack.sh
```

При невозможности разрешить зависимости нового окружения не модифицируйте старое `wgs-gatk`. Сначала проверьте конфликт пакетов и дайте лог solver; не снимайте автоматически версии всех инструментов.

## 2. Индексирование референсов

Индексы BWA v1 (`.bwt`, `.sa`) **не совместимы** с BWA-MEM3. BWA-MEM3 создаёт рядом с той же FASTA `.amb`, `.ann`, `.pac`, `.bwt.2bit.64`. `FASTA`, `.fai`, `.dict` остаются прежними; для нового индекса `samtools faidx` и `GATK CreateSequenceDictionary` повторно не нужны.

Для референсов белухи `GCA_059928025.1` и белого медведя `GCA_062004589.1` используйте пути `references.*.fasta` **из рабочей** конфигурации, а не имена из примеров. Например, в выделенной вычислительной сессии с достаточной памятью:

```bash
# Команду выполнять на вычислительном узле с разрешёнными ресурсами;
# не строить индексы на login-узле sphinx.
conda activate wgs-gatk-mem3
bwa-mem3 index -t 8 /absolute/path/to/reference.fa
ls -lh /absolute/path/to/reference.fa.{amb,ann,pac,bwt.2bit.64}
```

Внимание: `bwa-mem3 index` перезапишет общие `.amb/.ann/.pac`, используемые также старым BWA. Не выполнять его одновременно со старым выравниванием на том же референсе. До миграции сохраните отдельную резервную копию старых индексов или индексируйте отдельную копию FASTA. Заложите десятки гигабайт диска для каждого WGS-референса и достаточную память/ресурсы вычислительного узла.

## 3. Проверка без запуска реальных задач

```bash
cd WGS_GATK_simple
./run.sh validate --config /ABS/PATH/config.yaml
./run.sh plan --config /ABS/PATH/config.yaml
```

Проверить, что `validate` не ругается на отсутствующий `.bwt.2bit.64`; это **не** доказательство научной корректности выравнивания.

## 4. Пилотное сравнение

На одинаковой небольшой паре FASTQ сравнить:

- BWA-MEM 0.7.18 с его индексом (исходный commit, старое окружение) и BWA-MEM3 0.15.1 `--compat=bwa-mem` с новым индексом.
- Один и тот же `-t`, `-R`, reference FASTA и по возможности одинаковая версия samtools. Параметры `--compat=bwa-mem` документированы для BWA 0.7.19, поэтому **не ожидайте гарантированной bitwise identity с 0.7.18**. `@PG` будет отличаться.
- Метрики `samtools flagstat`, `samtools stats`, `samtools idxstats`, доли mapped/paired/duplicate, различия позиционирования/MAPQ/CIGAR, при необходимости GVCF/VCF.
- `sacct`: wall clock, CPU time, peak RSS, стоимость CPU·h. Один быстрый запуск не доказывает общую производительность на полном WGS.

Выполнять сравнение в **новых отдельных рабочих директориях** и не менять существующие BAM/GVCF/VCF. После выбора нового aligner зафиксировать версию бинарника, версии остальных инструментов и пути к референсам.

## 5. Возобновление задач

Изменение `align_lane.sh` меняет signature задач `align`, а новые индексы попадают в `inputs`. На старом `workdir` это пометит соответствующие BAM и downstream-производные как устаревшие и может вызвать полный перерасчёт. **По умолчанию для пилота используйте новый `workdir`.** Старые checkpoint и готовые результаты не удалять.
