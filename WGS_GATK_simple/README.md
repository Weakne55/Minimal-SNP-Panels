# WGS → cohort VCF (BWA-MEM + GATK4)

Простой Bash/Python/SLURM-конвейер для WGS белух (`Delphinapterus_leucas`) и белых медведей (`Ursus_maritimus`). Вход — парные FASTQ **или** готовые координатно отсортированные BAM/CRAM; выход — **один сырой совместно генотипированный VCF по виду**.

**Важно:** конвейер делает *variant calling / joint genotyping*, но **не выполняет финальную фильтрацию SNP, PCA или отбор минимальной панели**. Эти задачи идут после получения VCF.

## Быстрый старт

На сервере с установленным SLURM и Conda:

```bash
cd WGS_GATK_simple
conda env create -f envs/wgs.yaml   # один раз
conda activate wgs-gatk

cp config/config.example.yaml config/config.yaml
# 1. Изменить config/config.yaml (путь к данным, папка результатов,
#    референсы и реальный список хромосом каждого вида).
# 2. Заполнить config/samples.tsv (примеры в config/samples.example.tsv).

bash scripts/check_stack.sh            # проверить установленные инструменты
./run.sh validate                    # проверить входы и референсы
./run.sh plan                        # увидеть число задач, ничего не запускать
./run.sh submit                      # отправить задачи в SLURM
./run.sh status                      # увидеть активные SLURM job ID
```

Если задача упала: посмотрите логи, устраните причину и повторите `./run.sh submit`. Готовые этапы будут пропущены, а зависимые от незавершённых этапы повторены.

Для тестового запуска без SLURM: `./run.sh local` (обрабатывает все задачи последовательно; полногеномные данные так запускать обычно не следует).

## Что нужно настроить

**1. `config/config.yaml`.** Обычно меняются только четыре вещи:

- `samples` — путь к таблице образцов, по умолчанию `config/samples.tsv`;
- `workdir` — **отдельная** папка результатов, вне папки с кодом;
- `references.*.fasta` — путь к собранному FASTA соответствующего вида;
- `references.*.contigs` — названия целевых аутосом **в порядке `.fai`**.

Для получения списка настоящих имён хромосом:

```bash
python3 scripts/list_reference_contigs.py /data/references/beluga.fa.fai
```

Нельзя оставлять `[CONTIG1, CONTIG2]` в рабочей конфигурации: это заглушки. Для текущего исследования ожидаются 21 аутосома белухи и 36 белого медведя, но **названия** берутся именно из выбранного FASTA `.fai`, а не придумываются вручную. Наличие всех остальных контигов в FASTA допустимо: в анализ попадут только выбранные `contigs`.

Для каждого FASTA заранее должны существовать `ref.fa.fai`, `ref.dict`, а также индексы BWA (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`). Конвейер **не** перестраивает эти индексы автоматически. Если они ещё не созданы:

```bash
bwa index /data/references/beluga.fa
samtools faidx /data/references/beluga.fa
gatk CreateSequenceDictionary -R /data/references/beluga.fa
```

**2. `config/samples.tsv`.** Таблица с TAB-разделителями и заголовком:

```text
sample  species  input_type  library  lane  r1  r2  bam  cram
```

Здесь столбцы разделены табуляциями (не пробелами). Полный пример — `config/samples.example.tsv`.

- FASTQ: одна строка = одна парная lane; заполнены `r1` и `r2`, остальные поля файлов пусты. Несколько lane одного образца задаются несколькими строками с одинаковыми `sample`, `species`, `library` и разными `lane`.
- BAM: одна строка на образец, заполнен только `bam` (из путей). Нужен индекс `.bai`/`.csi`.
- CRAM: одна строка на образец, заполнен только `cram`. Нужен индекс `.crai`. Файл будет конвертирован в BAM.

Готовые BAM/CRAM должны быть coordinate-sorted, иметь корректный `@RG SM`, соответствовать **точно тому же референсу**, что указан в конфигурации, и быть готовы к variant calling. Для BAM/CRAM пайплайн **не** повторяет markdup и не исправляет уже существующие ошибки подготовки.

Имена образцов, видов, библиотек, дорожек: буквы латиницы, цифры, `_`, `-`, `.`. В одном образце FASTQ поддерживается только одна библиотека; несколько lane этой библиотеки поддерживаются.

## Куда складываются результаты

Внутри `workdir`:

```text
workdir/
├── bam/SAMPLE.bam                    # BAM после markdup (либо конвертации CRAM)
├── qc/                               # flagstat, samtools stats, mosdepth, MultiQC
│   └── multiqc/multiqc_report.html
├── gvcf/SPECIES/SAMPLE/CONTIG.g.vcf.gz
├── genomicsdb/SPECIES/CONTIG/        # рабочая база GATK
├── joint/SPECIES/CONTIG.vcf.gz       # совместный VCF по хромосоме
├── vcf/SPECIES/cohort.raw.vcf.gz     # ОСНОВНОЙ РЕЗУЛЬТАТ
├── vcf/SPECIES/cohort.raw.vcf.gz.csi # индекс
├── jobs/logs/                       # логи stdout/stderr каждого задания SLURM
├── jobs/scripts/                    # фактические команды каждого задания
├── jobs/submission.json             # номера отправленных SLURM jobs
├── checkpoints/                     # отметки успешно завершённых заданий
└── pipeline_plan.json               # сохранённый план выполнения
```

**Внимание:** `cohort.raw.vcf.gz` — *raw* VCF, ещё не окончательно QC-filtered. Для downstream SNP-анализа нужна отдельная стадия QC и фильтрации.

## Как работает конвейер

```text
samples.tsv + config.yaml
          |
          +--> FASTQ --> [опциональный fastp] --> BWA-MEM
          |                                       |
          |                         samtools sort + markdup
          |                                       |
          +--> готовый BAM/CRAM ------------------+
                                                  |
                                                BAM QC
                                                  |
                                  HaplotypeCaller: GVCF/sample/contig
                                                  |
                                     GenomicsDBImport: per contig
                                                  |
                                     GenotypeGVCFs: per contig
                                                  |
                                      bcftools concat: per species
                                                  |
                                           cohort.raw.vcf.gz
```

- `./run.sh` — единственная команда, которую запускает пользователь.
- `scripts/model.py` — проверка TSV и референсов.
- `scripts/workflow.py` — описание порядка этапов и формирование плана.
- `scripts/pipeline.py` — проверка готовности задач, отправка в SLURM, повторный запуск.
- `scripts/run_array.sh` — вызов одной задачи внутри SLURM array.
- `scripts/stages/*.sh` — **сами биоинформатические команды**. Чтобы понять, что делает конкретная стадия, откройте один соответствующий `.sh`.

Дополнительно см. [объяснение всех стадий](docs/PIPELINE.md) и [инструкцию для сопровождающего](docs/MAINTENANCE.md).

## Проверка и передача другим людям

```bash
python3 -m pytest -q              # структурные тесты без реальных WGS и без SLURM
bash -n run.sh scripts/run_array.sh scripts/stages/*.sh
./run.sh validate                 # на вашей настоящей конфигурации
./run.sh plan
```

Тесты используют `pytest` (поставляется как отдельная тестовая зависимость). **Структурные тесты не заменяют реальный end-to-end запуск с BWA/GATK и отдельную проверку на SLURM.** После изменений рекомендуется сначала выполнить запуск на маленьких тестовых WGS и двух коротких контигах.

При передаче проекта не включайте сырые FASTQ, большие BAM, рабочие GenomicsDB и приватные данные в репозиторий. Архив кода и отдельно настроенные пути/данные достаточно для воспроизводимого развёртывания.

## Архитектура после рефакторинга 29.10.2026

`model.py` читает и проверяет `samples.tsv` и референсы. `workflow.py`
только формирует задачи, а `pipeline.py` запускает, отслеживает и возобновляет
их через локальный режим либо SLURM. В `workflow.py` функция `make_tasks()`
явно вызывает четыре независимых этапа планирования: FASTQ, подготовка BAM,
совместное генотипирование и MultiQC. Объект `Task` — типизированный словарь,
совместимый с существующими JSON-планами. Команды биоинформатики по-прежнему
находятся в `scripts/stages/*.sh` и при рефакторинге не изменялись.

Отдельно от WGS находится `../population_structure/prepare.py` — только
заглушка формирования плана будущего анализа популяционной структуры.
`../snp_panel/` содержит контракт будущего независимого конвейера, который
принимает VCF/BCF и `population_labels.tsv` без необходимости запускать WGS.

Проверки без вычислений: `python3 -m pytest -q` в `WGS_GATK_simple`.
