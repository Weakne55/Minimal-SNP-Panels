# Ветка `snakemake`: перенос WGS на Snakemake + Singularity/Apptainer

## Что реализовано

- `FASTQ TSV` (оба вида, несколько runs на animal) → **fastp (опционально)** → **BWA-MEM3 + samtools** → run BAM → merge → markdup BAM/BAI → `samtools`/`mosdepth` QC → **GATK HaplotypeCaller** для каждой аутосомы → отдельные gVCF/CSI.
- По каждому `species × contig`: sample map → **GATK GenomicsDBImport** → **GenotypeGVCFs** → индексированный VCF.gz.
- По каждому виду: **bcftools concat** → `cohort.raw.bcf`/CSI → `cohort.raw.vcf.gz`/CSI.
- Опционально: технический variant QC через существующий `stages/qc.sh` (`variant_qc: true`), выход `qc/<species>/filtered.bcf`/CSI и `bcftools.stats.txt`.

Исходные Bash-стадии **не переписаны по алгоритму**. Перенесены оркестрация, параллельность, зависимости, повторный запуск и контейнерный runtime. Каждый `sample × contig` становится независимым заданием Snakemake. В отличие от старого `slurm_wgs.py`, это **не один job-array с последовательными хромосомами**.

Структура по существу соответствует исходному WGS-пайплайну: разработка генетических меток и минимальных SNP-панелей не подключена к CLI исходного snapshot; `structure.py` содержит недостижимый код. Не следует трактовать этот перенос как научно валидированную реализацию PCA, классификации или nested CV. Эти блоки требуют отдельной ревизии.

## 1. Создать ветку на сервере

```bash
git fetch origin
git switch develop
git pull --ff-only origin develop
git switch -c snakemake
# Применить предоставленный snakemake_migration.patch:
git apply /path/to/snakemake_migration.patch
# Контроль:
git status --short
```

Если `snakemake` существует на сервере, можно использовать `git switch snakemake` и `git pull` вместо создания новой ветки. Сохраните незакоммиченные изменения перед переключением. После проверок `git add`, `git commit`, `git push -u origin snakemake`. Исходная `develop` не изменяется.

## 2. Установить *только* Snakemake-оркестратор на login-узле

```bash
conda env create -f envs/snakemake.yaml
conda activate minsnp-snakemake
snakemake --version
which singularity || which apptainer
```

Вся тяжёлая биоинформатика исполняется **внутри контейнера**; Snakemake в установленном на login-узле окружении только составляет DAG и отправляет задания через SLURM. Не нужен `fakeroot`, не нужно собирать контейнер на вычислительном узле.

## 3. Сформировать контейнер без fakeroot на HPC

В репозитории `containers/Dockerfile.snakemake` (на базе существующего `envs/wgs-gatk-mem3.yaml`). Важно: в отличие от прежнего `Dockerfile`, он добавляет все бинарные файлы окружения в системный `PATH`: `apptainer exec` **не** запускает Docker entrypoint.

**Предпочтительно — через GitHub Actions.** Ветка `snakemake` содержит `.github/workflows/build-wgs-container.yml`. Workflow собирает и публикует OCI-образ:

```text
ghcr.io/weakne55/minimal-snp-panels-wgs:snakemake
```

После первого успешного workflow и при доступности образа из кластера:

```bash
bash scripts/pull_sif.sh \
  docker://ghcr.io/weakne55/minimal-snp-panels-wgs:snakemake \
  containers/wgs-gatk-mem3.sif
```

**Альтернативно**, если имеется Docker/Podman на отдельном компьютере с доступом к интернету и OCI registry:

```bash
docker build -f containers/Dockerfile.snakemake -t ghcr.io/weakne55/minimal-snp-panels-wgs:snakemake .
docker push ghcr.io/weakne55/minimal-snp-panels-wgs:snakemake
# На HPC выполнить scripts/pull_sif.sh из примера выше.
```

Если GHCR-образ приватный, перед `pull` потребуется авторизация (`apptainer registry login` / `singularity remote login` согласно конфигурации). При закрытом доступе HPC к сети SIF можно собрать на другой машине и передать на сервер по `scp`. Не запускайте `singularity build` из `.def` с `--fakeroot` на кластере, если нет настроенных `subuid/subgid`.

**Проверить SIF:**

```bash
singularity exec containers/wgs-gatk-mem3.sif sh -c \
  'command -v bwa-mem3 && command -v samtools && command -v gatk && command -v bcftools && command -v mosdepth'
```

Для reproducibility предпочтительно указывать **тег с Git SHA или immutable digest**, а не `:snakemake`. Можно передать собственный SIF через `--image /path/to/<sha>.sif`. Существование образа и успешную загрузку нельзя считать подтверждёнными до реальной проверки на кластере.

## 4. Параметры данных и ресурсов

`config/snakemake.yaml` управляет `run_id`, TSV-манифестами, входными конфигами, SIF-путём, включением `fastp`/QC и ресурсами отдельных правил. `config/references.yaml` и `config/pipeline.yaml` сохраняют прежнюю семантику сборок, путей и базовых ресурсов.

Для полного референса ожидаются: FASTA, `.fai`, `.dict` (с заменённым суффиксом FASTA), BWA-MEM3 индексы `.amb`, `.ann`, `.pac`, `.bwt.2bit.64`. Они уже могут быть подготовлены старым pipeline. Автоматическая перестройка индексов на HPC не включена намеренно, чтобы не запустить дорогостоящую задачу без согласования.

**Внешние FASTQ и референсы находятся вне корня репозитория.** При выполнении через Apptainer нужно явно сделать их доступными из контейнера. Пример для путей `/mnt/tank`:

```bash
--bind /mnt/tank
```

В wrapper это преобразуется в `--apptainer-args '--bind /mnt/tank:/mnt/tank'`. Для путей на других файловых системах добавьте `--bind` повторно. На кластере должен быть общий доступ к `.sif`, референсам и outputs со всех SLURM-узлов.

### Первый dry-run (не отправляет задания)

```bash
python scripts/run_snakemake.py \
  --fastq inputs/beluga_two.fastq.tsv inputs/polar_bear_two.fastq.tsv \
  --run-id snakemake_pilot_001 \
  --backend slurm --bind /mnt/tank
```

Если хотите протестировать на одном образце — создайте небольшой FASTQ TSV со строками только выбранного животного (манифесты `*_two.fastq.tsv` содержат **по два животных**, не один).

### Реальный пилот на SLURM

```bash
python scripts/run_snakemake.py \
  --fastq inputs/beluga_two.fastq.tsv inputs/polar_bear_two.fastq.tsv \
  --run-id snakemake_pilot_001 \
  --backend slurm --jobs 3 --bind /mnt/tank --submit
```

Локальный запуск:

```bash
python scripts/run_snakemake.py --fastq inputs/beluga_two.fastq.tsv \
  --run-id local_trial --backend local --cores 8 --bind /mnt/tank --submit
```

При необходимости вызвать Snakemake напрямую:

```bash
snakemake -s workflow/Snakefile --configfile config/snakemake.yaml \
  --profile profiles/slurm --dry-run
```

Параметры `--submit` и `--dry-run` wrapper — **не** команда `sbatch`: wrapper передаёт управление `snakemake`. Повторное выполнение той же команды `--submit` после ошибки возобновляет незавершённые правила. `--rerun-incomplete` включён в профили. Запускайте только одну Snakemake-инстанцию на `run_id` одновременно.

### SLURM и ограничения параллельности

`profiles/slurm/config.yaml`: `executor: slurm`, `jobs: 6`, `slurm_partition: main`, `slurm_extra: --nodelist=meduza-1` (запланированный AVX2-узел для BWA-MEM3). Проверьте `sinfo`, `scontrol show node meduza-1`; если нода недоступна или узкое место по RAM — поправьте профиль. В разных правилах из `config/snakemake.yaml` можно изменять `threads`, `mem_mb` и `runtime` (минуты). **В отличие от старого массива `max_parallel_samples: 3`, ограничение задаёт `--jobs`, а не число животных**. Например, `--jobs 3` допускает до трёх одновременных заданий любой стадии.

Здесь BWA-MEM3 и `samtools` идут в одном конвейере и уже ограничивают потоки и сортировку по выделенной RAM. GATK `-Xmx` ограничен через `MIN_SNP_MEMORY_MB`. Каждая аутосома имеет отдельное время и лимит, из-за чего можно эффективно перераспределять ресурсы между видами.

## 5. Результаты

```text
results/<run-id>/
├── cohort.lock.json                  # фиксирует cohort + конфигурацию
├── software_versions.json            # версии инструментов из контейнера
├── manifests/<species>/gvcf_manifest.tsv
├── samples/<species>/<sample>/wgs/
│   ├── lanes/<run>/<library>/fixed.bam
│   ├── bam/<sample>.merged.bam, <sample>.bam(.bai)
│   ├── qc/<sample>.flagstat.txt, .stats.txt, .mosdepth.summary.txt
│   └── gvcf/<sample>/<contig>.g.vcf.gz(.csi)
├── joint/<species>/
│   ├── sample_maps/<contig>.tsv
│   ├── genomicsdb/<contig>/
│   ├── per_contig/<contig>.vcf.gz(.csi)
│   ├── cohort.raw.bcf(.csi)
│   └── cohort.raw.vcf.gz(.csi)
├── qc/<species>/filtered.bcf(.csi)   # если variant_qc: true
└── logs/
```

Задачи внутри каждого вида разделены по аутосомам. Порядок contigs при `bcftools concat` точно такой, как в `config/references.yaml`. При изменении cohorts, references.yaml или pipeline.yaml для существующего `run-id` preflight выдаёт ошибку вместо тихой замены cohort. Входные FASTQ контролируются по размеру и `mtime`, не по полной хеш-сумме.

## 6. Проверки и известные ограничения

```bash
PYTHONPATH=src python -m pytest -q
for f in stages/*.sh scripts/pull_sif.sh; do bash -n "$f"; done
python scripts/run_snakemake.py --help
# После получения SIF:
python scripts/run_snakemake.py --run-id dryrun_001 --backend slurm --bind /mnt/tank
```

- **Dry-run и модульные тесты не равнозначны сквозному WGS-run.** Нужен пилот FASTQ → BAM → gVCF → GenomicsDB → multi-sample BCF/VCF с проверкой заголовков и индексов.
- Стадия `stages/qc.sh` сохранена из архивной версии и требует научной проверки порогов (особенно для небольшой когорты и low-coverage WGS); `variant_qc: false` по умолчанию.
- Существующие `stages/align.sh`, `merge_bam.sh`, `markdup.sh`, `concat.sh`, `export_vcf.sh` запрещают перезаписывать корректно существующие outputs. При изменении кода, входов или параметров уже завершённого запуска используйте **новый `run-id`**. Для продолжения **прерванных** задач без изменения конфигурации используйте тот же `run-id`; Snakemake заново построит только недостающие outputs.
- Результаты VCF/BCF сами по себе **не** означают надёжного определения подпопуляций или валидной SNP-панели.
- Подгрузка OCI с `docker://` требует доступного реестра; SIF предпочтительнее для разделённой HPC-среды. Состав референсных индексов должен совпадать с FASTA.
- В контейнере GATK GenomicsDB может создавать native-memory/temporary нагрузку сверх `Xmx`. Проверяйте `sacct`/`seff` и меняйте `mem_mb` после пилота.
