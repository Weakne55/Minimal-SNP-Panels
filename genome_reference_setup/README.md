# Подготовка референсов белухи и белого медведя (NCBI, октябрь 2026)

Это **отдельная одноразовая операция**. `WGS_Bash_GATK` только читает готовые референсы и индексы, никогда не копирует и не пересоздаёт их.

## Сборки

| Вид | NCBI GenBank Assembly | Название | Источник |
|---|---|---|---|
| Белуха — *Delphinapterus leucas* | `GCA_059928025.1` | `BTC_Dleu_Yukka2` | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_059928025.1/ |
| Белый медведь — *Ursus maritimus* | `GCA_062004589.1` | `BTC_Umar_Aika3` | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_062004589.1/ |

Сборки выбраны в соответствии с рекомендацией участников проекта. Ранее предложенная сборка медведя из Dryad **больше не используется**. В каталогах оставляем исходные FASTA с accession-содержащими именами и metadata NCBI, не переименовывая хромосомы.

## Подготовка окружения

На сервере `sphinx`, в `wgs-gatk`:

```bash
conda activate wgs-gatk
command -v datasets || conda install -c conda-forge ncbi-datasets-cli
command -v bwa samtools gatk unzip python3
```

Оба FASTA скачиваются через `datasets download genome accession ... --include genome,seq-report`, без генных аннотаций и без других сборок.

## 1. Скачивание

```bash
REFROOT=/mnt/tank/scratch/plem/data/genomes
bash setup_genomes.sh download --root "$REFROOT"
```

Команда сохраняет ZIP и распакованную FASTA в отдельных каталогах. Существующие каталоги КРС не затрагиваются. Если FASTA уже есть, скачивание пропускается. Перезапись готовых файлов не выполняется.

## 2. Индексирование — один раз, через SLURM

Для работы WGS-пайплайна требуются:
- `<FASTA>.fai` от `samtools faidx`;
- `<FASTA-без-расширения>.dict` от `gatk CreateSequenceDictionary`;
- `<FASTA>.amb`, `.ann`, `.bwt`, `.pac`, `.sa` от `bwa index -a bwtsw`.

```bash
sbatch --job-name=genomes_index --cpus-per-task=2 --mem=24G \
  --time=12:00:00 --output=genomes_index_%j.out \
  --wrap="bash $(pwd)/setup_genomes.sh index --root $REFROOT"
```

Укажите `--partition` или `--account`, только если они требуются вашим кластером. Индексирование может быть ресурсоёмким — его следует выполнять на вычислительном узле. Скрипт сохраняет существующие непустые индексы и строит только отсутствующие.

## 3. Проверка

```bash
bash setup_genomes.sh check --root "$REFROOT"
bash setup_genomes.sh inspect --root "$REFROOT" | tee reference_contigs_report.txt
```

`inspect` печатает реальные FASTA, порядок и длины последовательностей и предлагает **только одну длинную последовательность для smoke test**, без утверждения, что она обязательно аутосома. Для полного анализа состав аутосом следует установить по `sequence_report.jsonl`, который включён в пакет NCBI. Не выбирайте автоматически первые N записей `.fai`.

## Ожидаемые каталоги

```text
/mnt/tank/scratch/plem/data/genomes/
├── Delphinapterus_leucas_GCA_059928025.1/
│   ├── ncbi_dataset.zip
│   └── ncbi_dataset/data/GCA_059928025.1/
│       ├── GCA_059928025.1_*_genomic.fna
│       ├── GCA_059928025.1_*_genomic.fna.fai
│       ├── GCA_059928025.1_*_genomic.dict
│       └── GCA_059928025.1_*_genomic.fna.{amb,ann,bwt,pac,sa}
└── Ursus_maritimus_GCA_062004589.1/
    ├── ncbi_dataset.zip
    └── ncbi_dataset/data/GCA_062004589.1/
        ├── GCA_062004589.1_*_genomic.fna
        ├── GCA_062004589.1_*_genomic.fna.fai
        ├── GCA_062004589.1_*_genomic.dict
        └── GCA_062004589.1_*_genomic.fna.{amb,ann,bwt,pac,sa}
```

Фактическое имя FASTA зависит от названия сборки в архиве NCBI. Для обнаружения всех референсов:

```bash
find "$REFROOT/Delphinapterus_leucas_GCA_059928025.1" "$REFROOT/Ursus_maritimus_GCA_062004589.1" \
  -name '*_genomic.fna' -type f -print
```

В итоговом `WGS_Bash_GATK/config/config.yaml` используйте **точные абсолютные пути к этим файлам**, а `contigs` впишите по результатам анализа sequence report. Смешивать координаты VCF разных сборок нельзя без отдельного remapping/liftover.

## Ограничения

- Эти индексы создаются вне WGS-пайплайна; он не модифицируется.
- Генная аннотация (GFF3/GTF) сейчас не нужна. Если авторы предоставят её, убедитесь в соответствии сборке и идентификаторам последовательностей.
- Скрипт проверяет наличие индексов и соответствие `.fai`/`.dict` (имена, длины, порядок). Это ещё не реальный тест картирования.
- Сборка с лучшими N50/BUSCO не обязательно автоматически даёт лучший результат популяционной генетики: при тестировании полезно измерить mapping rate, coverage, uniformity и возможное смещение относительно происхождения референсной особи.
- ZIP сохраняется после распаковки; после проверки целостности его можно удалить для экономии места.
