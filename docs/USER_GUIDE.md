# Руководство по эксплуатации — Minimal SNP Panels

Подробные инструкции к Bash/Python-конвейеру для получения совместного VCF/BCF из полногеномных FASTQ и разработки минимальных SNP-панелей, разделяющих генетические популяции белухи (*Delphinapterus leucas*) и белого медведя (*Ursus maritimus*).

**Версия документации:** ветка `develop`, архитектура поэтапных SLURM arrays (ноябрь 2026). Все команды ниже выполняются **из корня репозитория**, если не указано иное. Замените пути, раздел SLURM и ресурсы на значения своего кластера.

> **Статус:** FASTQ → VCF успешно проверен на одном WGS; оба аналитических режима — заглушки.

## Содержание

1. [Режимы работы и архитектура](#1-режимы-работы-и-архитектура)
2. [Скачать репозиторий и подготовить окружение](#2-скачать-репозиторий-и-подготовить-окружение)
3. [Скачать и подготовить референсные геномы](#3-скачать-и-подготовить-референсные-геномы)
4. [Скачать FASTQ и создать входной TSV](#4-скачать-fastq-и-создать-входной-tsv)
5. [Настроить SLURM, ресурсы и фильтры](#5-настроить-slurm-ресурсы-и-фильтры)
6. [Проверить и запустить FASTQ → VCF](#6-проверить-и-запустить-fastq--vcf)
7. [Отслеживать задания и читать логи](#7-отслеживать-задания-и-читать-логи)
8. [Разобраться в результатах FASTQ → VCF](#8-разобраться-в-результатах-fastq--vcf)
9. [Режим VCF → population structure](#9-режим-vcf--population-structure)
10. [Режим VCF + metadata → SNP panel](#10-режим-vcf--metadata--snp-panel)
11. [Ошибки, возобновление и воспроизводимость](#11-ошибки-возобновление-и-воспроизводимость)
12. [Источники и технические ограничения](#12-источники-и-технические-ограничения)

## 1. Режимы работы и архитектура

| Команда | Вход | Результат | Текущее состояние |
|---|---|---|---|
| `--fastq` | TSV со ссылками на paired-end FASTQ | BAM, QC, gVCF по образцу, совместные VCF.gz и BCF по виду | Проверен на одном WGS; научная валидация не завершена |
| `--vcf` | Индексированный multi-sample VCF.gz/BCF | В будущем: PCA/кластеризация и `population_labels.tsv` | Заглушка: код возврата `2` |
| `--vcf --metadata --species` | Индексированный VCF/BCF + таблица популяций | QC, оценка SNP-панелей, CV, модель | Не реализован (заглушка) |

Вариант `--labels` — синоним `--metadata`. Виды не объединяются в одном joint VCF. Стадии WGS:

```text
FASTQ R1/R2
  ├─ trim                [опционально; одна задача на sequencing run/library]
  ├─ align               [одна задача на sequencing run/library]
  ├─ merge_bam           [одна задача на sample]
  ├─ markdup             [одна задача на sample]
  ├─ bam_qc              [одна задача на sample]
  ├─ haplotypecaller     [ОДИН вызов на sample: ВСЕ указанные аутосомы → один gVCF]
  ├─ genomicsdb          [ОДНА задача на species: все аутосомы → общий workspace]
  ├─ genotype            [ОДНА задача на species: все аутосомы → общий BCF]
  └─ export_vcf          [одна задача на species: bgzip VCF]
```

Для `genomicsdb` GATK получает все аутосомы вида через несколько `-L` и запускается с `--reader-threads 1` (ограничение GATK при нескольких интервалах). `genotype` также получает все аутосомы в одном вызове. Вызовов `concat` нет. **Учитывайте повышенные требования к RAM, диску и времени** на этих стадиях.

**SLURM:** каждая активная стадия — отдельный array, следующий array зависит от завершения предыдущего через `afterok`. `slurm.max_concurrent_jobs` ограничивает число одновременно исполняющихся array tasks **в пределах одного запуска (`run-id`)**. Он не ограничивает несколько независимо запущенных `run-id`.

## 2. Скачать репозиторий и подготовить окружение

### 2.1. Требования

- Linux и файловая система, доступная вычислительным узлам; достаточно места для FASTQ, BAM, gVCF, GenomicsDB и совместных VCF. Нужны права записи в `results/`.
- Для SLURM: команды `sbatch`, `squeue`, желательно `sacct`, и настроенный доступ к compute nodes. Для BWA-MEM3 процессор узла должен поддерживать **AVX2**.
- `git`, `conda` (Miniforge/Miniconda), `curl`, `unzip`; для работы `python`, `bwa-mem3`, `samtools`, `bcftools`, `gatk`, `fastp`, `mosdepth`.
- В момент установки и скачивания референсов/FASTQ нужен доступ к NCBI/ENA. Уже скачанные данные можно использовать без повторной загрузки.

### 2.2. Загрузка кода

```bash
git clone --branch develop --single-branch \
  https://github.com/Weakne55/Minimal-SNP-Panels.git
cd Minimal-SNP-Panels
pwd
```

Если получен ZIP ветки `develop`, распакуйте архив и перейдите в полученную папку с `pipeline.py`. Повторно клонировать репозиторий не требуется. **Пока изменения из обновлённого ZIP не закоммичены в ветку GitHub `develop`, `git clone` может скачать более старый код:** для воспроизведения именно описанной здесь архитектуры используйте обновлённый архив или опубликуйте его изменения в репозитории. Проверка версии: в `config/pipeline.yaml` должно быть поле `slurm.max_concurrent_jobs`.

### 2.3. Conda-окружение WGS/анализа

```bash
conda env create -f envs/wgs-gatk-mem3.yaml
conda activate wgs-gatk-mem3

python --version
python pipeline.py --help
python pipeline.py --check-tools
```

Окружение включает Python 3.11, BWA-MEM3 0.15.1, GATK 4.7.0.0, samtools 1.24, bcftools 1.24, fastp 1.4.0, mosdepth 0.3.14, plink2, NumPy, scikit-learn, matplotlib, PyYAML и pytest. Отдельная установка Python-пакета не требуется: `pipeline.py` подключает `src/` самостоятельно.

`--check-tools` **печатает версии, но не завершается с ошибкой из-за отсутствующего инструмента**. Поэтому явно проверьте наличие программ:

```bash
for cmd in python bwa-mem3 samtools bcftools gatk fastp mosdepth; do
  command -v "$cmd" || echo "ОТСУТСТВУЕТ: $cmd"
done
command -v sbatch || echo 'SLURM sbatch отсутствует (для local backend допустимо)'
```

При запуске SLURM-заданий окружение должно быть активировано **перед** `--submit`: отправляемые jobs наследуют переменные окружения (`--export=ALL`). Убедитесь, что conda-пути и зависимости доступны с compute nodes. Для проверки AVX2 **на выделенном узле**: `grep -qw avx2 /proc/cpuinfo && echo AVX2_OK`.

### 2.4. Минимальная структура файлов

```text
Minimal-SNP-Panels/
├── pipeline.py                       # CLI
├── config/
│   ├── pipeline.yaml                # ресурсы, SLURM, QC, панели
│   └── references.yaml              # FASTA и список аутосом для каждого вида
├── envs/wgs-gatk-mem3.yaml
├── inputs/                           # шаблоны и примеры TSV
├── stages/                           # Bash-стадии
├── src/minsnp/                       # Python-модули
├── data/                             # FASTQ/референсы (не коммитить в Git)
└── results/                          # результаты по run-id (не коммитить в Git)
```

## 3. Скачать и подготовить референсные геномы

Используются **две строго зафиксированные сборки**:

| Вид / ключ в конфиге | Сборка NCBI | Аутосомы |
|---|---|---:|
| Белуха / `beluga` | [`GCA_059928025.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_059928025.1/) | 21 |
| Белый медведь / `polar_bear` | [`GCA_062004589.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_062004589.1/) | 36 |

**Необходимые внешние файлы:** `*_genomic.fna` каждой сборки и локально создаваемые индексы FASTA (`.fai`), GATK (`.dict`), BWA-MEM3 (`.amb`, `.ann`, `.pac`, `.bwt.2bit.64`). NCBI `sequence_report.jsonl` полезен для независимой проверки соответствия аутосом. **Дополнительные annotation GFF/GTF, dbSNP, known-sites VCF и наборы гаплотипов текущим кодом не требуются.**

Если референсы и индексы уже находятся в общем хранилище, **не скачивайте и не индексируйте их повторно**: перейдите к разделу 3.3 и проверьте `config/references.yaml`.

### 3.1. Установка NCBI Datasets CLI (отдельное окружение)

```bash
conda create -y -n ncbi-datasets -c conda-forge ncbi-datasets-cli
conda run -n ncbi-datasets datasets --version
```

Скачивание геномов в пределах репозитория (на сервере можно выбрать любое другое общее хранилище):

```bash
REF_ROOT="$PWD/data/references"
mkdir -p "$REF_ROOT/beluga" "$REF_ROOT/polar_bear"

conda run -n ncbi-datasets datasets download genome accession GCA_059928025.1 \
  --include genome,seq-report --filename "$REF_ROOT/beluga/ncbi_dataset.zip"
conda run -n ncbi-datasets datasets download genome accession GCA_062004589.1 \
  --include genome,seq-report --filename "$REF_ROOT/polar_bear/ncbi_dataset.zip"

unzip -n "$REF_ROOT/beluga/ncbi_dataset.zip" -d "$REF_ROOT/beluga"
unzip -n "$REF_ROOT/polar_bear/ncbi_dataset.zip" -d "$REF_ROOT/polar_bear"

find "$REF_ROOT" -type f \( -name '*_genomic.fna' -o -name 'sequence_report.jsonl' \)
```

Официальная документация: [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).

### 3.2. Индексация (один раз на сборку)

Сначала активируйте вычислительное окружение WGS:

```bash
conda activate wgs-gatk-mem3
REF_ROOT="${REF_ROOT:-$PWD/data/references}"

BELUGA_FA=$(find "$REF_ROOT/beluga/ncbi_dataset/data/GCA_059928025.1" \
  -maxdepth 1 -type f -name '*_genomic.fna' -print -quit)
BEAR_FA=$(find "$REF_ROOT/polar_bear/ncbi_dataset/data/GCA_062004589.1" \
  -maxdepth 1 -type f -name '*_genomic.fna' -print -quit)

test -s "$BELUGA_FA" && test -s "$BEAR_FA" || { echo 'FASTA не найден'; exit 1; }
printf 'Белуха: %s\nМедведь: %s\n' "$BELUGA_FA" "$BEAR_FA"
```

**Индексация BWA-MEM3 ресурсоёмкая; запускайте команды на выделенном вычислительном узле с достаточной RAM, а не на login-node.** Например, если на вашем кластере разрешены интерактивные задания, сначала запросите ресурсы командой `srun --partition=main --cpus-per-task=8 --mem=128G --time=12:00:00 --pty bash` (параметры примерные и должны быть разрешены администратором). В полученной shell-сессии заново выполните `conda activate wgs-gatk-mem3` и команды из начала раздела 3.2 для определения `REF_ROOT`, `BELUGA_FA`, `BEAR_FA`; переменные из прежней shell не гарантированно перенесутся. Затем:

```bash
for FA in "$BELUGA_FA" "$BEAR_FA"; do
  samtools faidx "$FA"
  gatk CreateSequenceDictionary -R "$FA" -O "${FA%.*}.dict"
  bwa-mem3 index "$FA"
done
```

Для уже существующих `.fai`, `.dict`, BWA-MEM3 индексов **пропустите** соответствующую команду; не перезаписывайте проиндексированные референсы, пока ими пользуются запущенные задачи. Проверьте комплект файлов:

```bash
for FA in "$BELUGA_FA" "$BEAR_FA"; do
  for ext in .fai .amb .ann .pac .bwt.2bit.64; do
    test -s "${FA}${ext}" || echo "ОТСУТСТВУЕТ: ${FA}${ext}"
  done
  test -s "${FA%.*}.dict" || echo "ОТСУТСТВУЕТ: ${FA%.*}.dict"
done
```

> **Важно:** исторический `wgs_test/genome_reference_setup/setup_genomes.sh` использует **`bwa index`**, а не `bwa-mem3 index`. Его нельзя считать готовой процедурой для нынешнего BWA-MEM3-конвейера: старые `.bwt`/`.sa` **не заменяют** `.bwt.2bit.64`. Пользуйтесь приведёнными выше командами или уже проверенными индексами.

### 3.3. Прописать FASTA в `config/references.yaml`

Для каждого вида измените только `fasta` на **абсолютный путь** к скачанному `.fna`:

```yaml
species:
  beluga:
    assembly: GCA_059928025.1
    fasta: /absolute/path/to/beluga/ncbi_dataset/data/GCA_059928025.1/ASSEMBLY_genomic.fna
    contigs: # сохраните полный проверенный список 21 accession из репозитория
      - CM187741.1
      # ... остальные 20 из исходного config/references.yaml
  polar_bear:
    assembly: GCA_062004589.1
    fasta: /absolute/path/to/polar_bear/ncbi_dataset/data/GCA_062004589.1/ASSEMBLY_genomic.fna
    contigs: # сохраните полный проверенный список 36 accession из репозитория
      - CM197969.1
      # ... остальные 35 из исходного config/references.yaml
```

**Это фрагмент для объяснения структуры, не замена полного YAML**: не удаляйте существующие `contigs`. Списки аутосом уже записаны в репозитории. Не заменяйте их всеми contigs FASTA и не переносите индексы от другой версии сборки.

Проверка списка accession из конфига против FASTA `.fai`:

```bash
python - <<'PY'
from pathlib import Path
import yaml
refs = yaml.safe_load(Path('config/references.yaml').read_text())['species']
for name, spec in refs.items():
    fasta = Path(spec['fasta'])
    fai = Path(str(fasta) + '.fai')
    assert fasta.is_file() and fai.is_file(), (name, fasta)
    available = {line.split('\t')[0] for line in fai.read_text().splitlines()}
    missing = set(spec['contigs']) - available
    assert not missing, (name, sorted(missing))
    print(f"{name}: assembly={spec['assembly']}, autosomes={len(spec['contigs'])}, OK")
PY
```

## 4. Скачать FASTQ и создать входной TSV

### 4.1. Откуда брать FASTQ

Вход WGS — **полные парные Illumina FASTQ** (R1 и R2) для каждого `run_id`. Репозиторий содержит **только примеры манифестов**, не сами reads. Источник примеров — публичные run accessions NCBI SRA / ENA; сведения об образцах дополнительно собраны в `data/metadata/`.

Проверяйте принадлежность run к BioSample и виду, платформу, парность, наличие R1/R2 и доступность файлов по [ENA file report API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html). Пример для белухи `SRR35396908`:

```bash
RUN=SRR35396908
SPECIES=Delphinapterus_leucas
FASTQ_DIR="$PWD/data/wgs/$SPECIES/$RUN"
mkdir -p "$FASTQ_DIR"

curl -fsSLG 'https://www.ebi.ac.uk/ena/portal/api/filereport' \
  --data-urlencode "accession=$RUN" \
  --data-urlencode 'result=read_run' \
  --data-urlencode 'fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes' \
  -o "$FASTQ_DIR/ena.tsv"
cat "$FASTQ_DIR/ena.tsv"
```

Следующая команда скачивает **ровно два** FASTQ, проверяет размер и MD5 по ENA. Для другого run повторите блок с новым `RUN`/`SPECIES`; не запускайте несколько загрузчиков в одну директорию:

```bash
python - "$FASTQ_DIR/ena.tsv" "$FASTQ_DIR" <<'PY'
import csv, hashlib, os, subprocess, sys
from pathlib import Path
from urllib.parse import urlparse

metadata, output = Path(sys.argv[1]), Path(sys.argv[2])
with metadata.open(newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
if len(rows) != 1:
    raise SystemExit(f'Ожидался ровно один ENA run, получено {len(rows)}')
row = rows[0]
urls = row['fastq_ftp'].split(';') if row['fastq_ftp'] else []
md5s = row['fastq_md5'].split(';') if row['fastq_md5'] else []
sizes = row['fastq_bytes'].split(';') if row['fastq_bytes'] else []
if not (len(urls) == len(md5s) == len(sizes) == 2):
    raise SystemExit('Нет ровно двух FASTQ с MD5/size: проверьте ENA и тип библиотеки')
for link, md5, size_text in zip(urls, md5s, sizes):
    url = link if '://' in link else 'https://' + link
    path = output / Path(urlparse(url).path).name
    temp = path.with_name(path.name + '.part')
    expected_bytes = int(size_text)
    def valid(p):
        if not p.is_file() or p.stat().st_size != expected_bytes:
            return False
        digest = hashlib.md5()
        with p.open('rb') as fh:
            for chunk in iter(lambda: fh.read(16 * 1024 * 1024), b''):
                digest.update(chunk)
        return digest.hexdigest().lower() == md5.lower()
    if not valid(path):
        if path.exists():
            raise SystemExit(f'Существующий FASTQ повреждён: {path}')
        subprocess.run(['curl', '-fL', '--retry', '5', '-C', '-', '-o', str(temp), url], check=True)
        if not valid(temp):
            raise SystemExit(f'Ошибка размера/MD5 для {temp}; не используйте файл')
        os.replace(temp, path)
    print('VERIFIED', path, expected_bytes)
PY
```

Образцы в `inputs/beluga_two.fastq.tsv` и `inputs/polar_bear_two.fastq.tsv` — примеры **реальных имён файлов и ID**; их пути могут не существовать на другом сервере. Для внутренних данных белых медведей используйте свои FASTQ и корректные метаданные, не публикуя закрытые данные в Git.

### 4.2. FASTQ-манифест для `--fastq`

Обязательные столбцы **именно в таком смысле**:

| Столбец | Что указать |
|---|---|
| `sample_id` | Уникальный идентификатор животного, совпадающий с BAM `SM` и VCF sample name |
| `species` | **`beluga`** или **`polar_bear`** — точные ключи из `references.yaml` |
| `run_id` | ID sequencing run (например, `SRR35396908`) |
| `library` | ID библиотеки, например `SRX30479613` или внутренний `LIB01` |
| `fastq_r1` | Существующий R1 FASTQ.gz |
| `fastq_r2` | Существующий R2 FASTQ.gz |

Файл — **TSV с табуляциями**, не CSV. Пути к FASTQ могут быть абсолютными или относительными **от папки, содержащей TSV**. Каждая комбинация `(sample_id, run_id, library)` уникальна. Несколько строк с одинаковым `sample_id` означают несколько sequencing runs/libraries **одного животного**; они будут объединены в один BAM. Один FASTQ нельзя повторно использовать в других строках.

Пример `inputs/my_pilot.fastq.tsv` (приведены две реально использовавшиеся публичные пары; до запуска скачайте их по указанным путям):

```tsv
sample_id	species	run_id	library	fastq_r1	fastq_r2
SAMN51337163	beluga	SRR35396908	SRX30479613	../data/wgs/Delphinapterus_leucas/SRR35396908/SRR35396908_1.fastq.gz	../data/wgs/Delphinapterus_leucas/SRR35396908/SRR35396908_2.fastq.gz
SAMN16454148	polar_bear	SRR15170742	SRX11478086	../data/wgs/Ursus_maritimus/SRR15170742/SRR15170742_1.fastq.gz	../data/wgs/Ursus_maritimus/SRR15170742/SRR15170742_2.fastq.gz
```

Создайте `inputs/my_pilot.fastq.tsv` в текстовом редакторе (без Excel-замены табов пробелами), либо скопируйте и скорректируйте `inputs/beluga_two.fastq.tsv` / `inputs/polar_bear_two.fastq.tsv`. Можно передать **несколько TSV одновременно**, включая оба вида:

```bash
python pipeline.py --fastq inputs/beluga_two.fastq.tsv inputs/polar_bear_two.fastq.tsv \
  --backend slurm --run-id pilot_01
```

### 4.3. Файлы, которые **не** нужно скачивать для FASTQ → VCF

`population_labels.tsv`, опубликованные SNP-панели, SNP-аннотации, VCF известных вариантов и другие метки не требуются для WGS variant calling. Они понадобятся только на следующих исследовательских этапах. Готовый индексированный VCF/BCF можно подать сразу в аналитический режим без повторного FASTQ → VCF.

## 5. Настроить SLURM, ресурсы и фильтры

Основной конфиг — `config/pipeline.yaml`. До запуска посмотрите доступные разделы/узлы и ограничения:

```bash
sinfo -o '%P %a %l %D'
sinfo -N -o '%N %P %f'
```

Пример **схемы** (в репозитории уже содержатся все стадии; не заменяйте полный файл этим сокращённым фрагментом):

```yaml
trim_fastq: false        # true → дополнительная SLURM стадия fastp
logging:
  level: INFO            # DEBUG / INFO / WARNING / ERROR
slurm:
  partition: main        # имя раздела на вашем кластере
  nodelist: null         # либо конкретный узел с AVX2, например meduza-1
  max_concurrent_jobs: 3
  stages:
    trim:            {cpus: 4,  memory: 8G,  time: '12:00:00'}
    align:           {cpus: 16, memory: 64G, time: '2-00:00:00'}
    merge_bam:       {cpus: 4,  memory: 16G, time: '12:00:00'}
    markdup:         {cpus: 8,  memory: 32G, time: '1-00:00:00'}
    bam_qc:          {cpus: 4,  memory: 16G, time: '12:00:00'}
    haplotypecaller: {cpus: 8,  memory: 32G, time: '2-00:00:00'}
    genomicsdb:      {cpus: 8,  memory: 64G, time: '2-00:00:00'}
    genotype:        {cpus: 8,  memory: 48G, time: '2-00:00:00'}
    export_vcf:      {cpus: 4,  memory: 16G, time: '12:00:00'}
    qc:              {cpus: 4,  memory: 16G, time: '12:00:00'}
```

- `partition`: должен совпадать с доступным разделом SLURM. `nodelist` можно удалить или задать `null`, если хотите разрешить планировщику выбирать узел; в таком случае обеспечьте AVX2 на выделяемых для `align` узлах.
- `max_concurrent_jobs: N`: каждый array отправляется с `%N`. Например, в массиве `align` из 20 элементов при `N=3` одновременно запускаются максимум 3. **Это не лимит на весь кластер или несколько `run-id`**.
- `stages.<name>.cpus/memory/time`: **ресурсы на одну задачу**. В примере `align` при `N=3` теоретически запрашивает до 48 CPU и 192 ГБ памяти одновременно. Это верхняя граница запроса, не измеренное потребление.
- `trim_fastq: false`: исходные reads картируются без fastp. При `true` отдельная стадия пишет trimmed FASTQ и JSON/HTML отчёты.
- `logging.level`: уровень Python-логирования. Вывод samtools/BWA/GATK/fastp записывается в отдельные файлы независимо от уровня.

Для анализа VCF с метаданными дополнительно используются разделы:

```yaml
qc:
  min_dp: 5              # маскирование генотипов с низким DP, если FORMAT/DP существует
  min_gq: 20             # маскирование с низким GQ, если FORMAT/GQ существует
  max_missing: 0.10      # макс. доля отсутствующих генотипов на SNP
  min_maf: 0.01          # мин. minor allele frequency
panels:
  sizes: [10, 20, 50, 100, 200]
  folds: 5
  seed: 42
  ld_r2: 0.2
  selection_pool: 10000
```

Это **исходные настройки**, а не автоматически оптимальные параметры для каждого вида. Отбор панели опирается на частоты генотипов и labels; при малом числе особей в популяции число fold автоматически сокращается до размера наименьшей группы (но не ниже 2).

## 6. Проверить и запустить FASTQ → VCF

### 6.1. Проверить код, окружение и исходные данные

```bash
conda activate wgs-gatk-mem3
python -m compileall -q pipeline.py src/minsnp
for file in stages/*.sh; do bash -n "$file" || exit 1; done
PYTHONPATH=src python -m pytest -q
python pipeline.py --check-tools
```

**Известное замечание по тестам:** `test_structure_is_explicit_stub` в текущем архиве иногда падает из-за сравнения списка `Path.iterdir()` в фиксированном порядке (порядок файлов в каталоге не гарантирован). Это дефект самой проверки порядка, не реализация population structure; его следует исправить отдельно. Все остальные ошибки тестов требуют разбора.

Перед производственным запуском убедитесь, что все пути FASTQ существуют, `.fai` и BWA-MEM3 индексы соответствуют выбранным FASTA, а `partition` и `nodelist` существуют на кластере.

### 6.2. Dry run: **не отправляет задания**

```bash
python pipeline.py \
  --fastq inputs/beluga_two.fastq.tsv inputs/polar_bear_two.fastq.tsv \
  --backend slurm \
  --run-id pilot_01
```

Вывод содержит количество образцов, sequencing runs, видов, лимит одновременных SLURM tasks и оценку **запрашиваемых** CPU/RAM для стадий. В режиме SLURM без `--submit` **не создаются** план, результаты и задания. **Это не проверка работоспособности GATK/SLURM и не полный dry-run каждой внешней команды.**

### 6.3. Отправить поэтапные arrays в SLURM

```bash
python pipeline.py \
  --fastq inputs/beluga_two.fastq.tsv inputs/polar_bear_two.fastq.tsv \
  --backend slurm \
  --run-id pilot_01 \
  --submit
```

При успешной отправке создаются `results/pilot_01/slurm_plan.json` и `slurm_jobs.json`. Следующие стадии автоматически отправляются **с зависимостью** `afterok`, то есть начнутся только после успешного завершения **всех** задач предыдущего array.

Если нужно другое расположение результатов, **на обеих командах** используйте одинаковые `--outdir /path/to/shared/results`. Для альтернативных конфигов есть `--settings /path/to/pipeline.yaml` и `--references /path/to/references.yaml`.

Важные правила:

1. Для нового запуска выбирайте **новый** `--run-id`; не переиспользуйте ранее отправленный SLURM run.
2. До завершения запуска не изменяйте `pipeline.py`, `src/minsnp/*.py`, `stages/*.sh`, соответствующие TSV и YAML: worker проверяет их контрольные суммы и характеристики FASTQ.
3. Первым запуском лучше обработать **небольшую когорту** из 2–4 образцов; затем проверить BAM, gVCF и совместный VCF/BCF, прежде чем запускать полную выборку.

### 6.4. Локальный режим (без SLURM)

```bash
python pipeline.py \
  --fastq inputs/beluga_two.fastq.tsv \
  --backend local --run-id local_test_01
```

Все стадии выполняются последовательно в текущей shell-сессии. Это **не** лёгкий тест: полные WGS требуют больших CPU/RAM/дискового пространства. Для остановленного **локального** run допускается повторить команду с `--resume` (см. раздел 11). `--submit` с local backend недоступен.

## 7. Отслеживать задания и читать логи

### 7.1. Статус пайплайна

```bash
python pipeline.py status --run-id pilot_01
# Если использовали --outdir, обязательно повторите его:
python pipeline.py status --run-id pilot_01 --outdir /path/to/shared/results
```

Команда показывает по каждой стадии `completed/total` и известные состояния SLURM. Файл `completed/<stage>/<task_index>.json` появляется **только при успешном завершении worker**. Исчезновение Job ID из очереди ещё не означает успешный результат.

### 7.2. Состояние в SLURM

```bash
squeue -u "$USER" -r -o '%.18i %.24j %.10T %.12M %.20R'

# Все Job IDs конкретного запуска:
cat results/pilot_01/slurm_jobs.json

# Пример: информация по одному Job ID, взятому из slurm_jobs.json:
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
scontrol show job JOB_ID
```

Состояния `PENDING`, `RUNNING`, `COMPLETED`, `FAILED`, `OUT_OF_MEMORY`, `TIMEOUT`, `CANCELLED` относятся к SLURM; `DependencyNeverSatisfied` — **причина ожидания** задания, а не отдельное состояние. При отказе **любой** задачи стадии следующая array с `afterok` может остаться неисполненной. Для завершившихся заданий полезнее `sacct`, чем `squeue`.

### 7.3. Какие логи читать

| Файл/каталог | Содержимое |
|---|---|
| `results/<run-id>/logs/pipeline.log` | Запуск CLI, отправка arrays и основные сообщения оркестратора |
| `results/<run-id>/logs/workers/<stage>/<index>.log` | `START`, `DONE`, `FAILED`, задача/входные данные, Python traceback |
| `results/<run-id>/logs/slurm/<stage>_<job>_<task>.out` / `.err` | stdout/stderr планировщика для конкретного array element |
| `results/<run-id>/samples/<species>/<sample>/logs/*.log` | Выполнявшиеся команды и полный вывод BWA/samtools/GATK/fastp/mosdepth |
| `results/<run-id>/joint/<species>/logs/*.log` | GenomicsDBImport, GenotypeGVCFs, export |
| `results/<run-id>/logs/qc/*.log` / `logs/panels/` | При локальном аналитическом запуске: командные логи QC и журналы анализа |

Примеры просмотра **во время** вычислений:

```bash
tail -n 100 results/pilot_01/logs/pipeline.log
tail -F results/pilot_01/logs/workers/align/0.log
# Другой терминал, после появления файла:
tail -F results/pilot_01/logs/slurm/align_JOBID_0.err

# Найти ошибки во всех уже существующих worker-логах:
grep -RniE 'FAILED|ERROR|Traceback|Exception' results/pilot_01/logs/workers/ || true
```

В конкретной ошибке сначала найдите `stage` и `task index` в `status`, `sacct` или SLURM `.err`, затем откройте `logs/workers/<stage>/<index>.log` и **путь внешнего командного лога**, напечатанный worker. Для соответствия индекса животному/виду см. `results/<run-id>/slurm_plan.json` → `tasks[stage][index]` и `samples[sample_index]`.

## 8. Разобраться в результатах FASTQ → VCF

Структура для запуска `pilot_01` (каталоги создаются **по мере прохождения этапов**):

```text
results/pilot_01/
├── slurm_plan.json                    # фиксированный список стадий, задач и SHA256
├── slurm_jobs.json                    # stage → SLURM Job ID
├── versions.json                      # сведения об установленных инструментах
├── manifests/<species>/<sample>.tsv   # нормализованный input каждого sample
├── completed/<stage>/<index>.json     # успешно завершённые tasks
├── logs/
│   ├── pipeline.log
│   ├── workers/<stage>/<index>.log
│   └── slurm/<stage>_<job>_<task>.out/.err
├── samples/<species>/<sample>/
│   ├── wgs/lanes/<sample>/<run>/<library>/fixed.bam
│   ├── wgs/bam/<sample>.merged.bam
│   ├── wgs/bam/<sample>.bam           # итоговый BAM с markdup
│   ├── wgs/bam/<sample>.bam.bai
│   ├── wgs/qc/<sample>.dup.txt
│   ├── wgs/qc/<sample>.flagstat.txt
│   ├── wgs/qc/<sample>.stats.txt
│   ├── wgs/qc/<sample>.mosdepth.summary.txt
│   ├── wgs/gvcf/<sample>.g.vcf.gz    # ОДИН sample-wide gVCF (+ индекс)
│   ├── wgs/gvcf_manifest.tsv
│   └── logs/*.log
└── joint/<species>/
    ├── sample_map.tsv                # имена sample → путь к полному gVCF
    ├── genomicsdb/                  # одна база для всех аутосом
    ├── cohort.raw.bcf                # общий BCF (+ .csi)
    ├── cohort.raw.vcf.gz             # общий bgzip VCF (+ .csi)
    └── logs/*.log
```

### Что означает каждый биологический результат

| Артефакт | Интерпретация |
|---|---|
| `fixed.bam` | Картированные paired-end reads одной библиотеки/run; подготовлены для merge/markdup |
| `<sample>.bam` | Итоговое, координатно отсортированное картирование животного с помеченными (не удалёнными) дубликатами |
| `*.dup.txt` | Метрики дубликатов от `samtools markdup` |
| `flagstat.txt`, `stats.txt` | Числа reads, доли mapped/paired, различные статистики BAM |
| `mosdepth.summary.txt` | Среднее покрытие по хромосомам/регионам отчёта mosdepth; QC глубины |
| `<sample>.g.vcf.gz` | Результат GATK HaplotypeCaller в режиме `-ERC GVCF` для **одного** животного по перечисленным аутосомам; включает reference blocks, **не** готовую SNP-панель |
| `cohort.raw.bcf`, `cohort.raw.vcf.gz` | **Основные входы последующего анализа**: совместно генотипированные варианты всех животных одного вида; разные форматы одних и тех же данных, **ещё без аналитического QC** |

Проверка итогового VCF (замените вид на `polar_bear`, если необходимо):

```bash
VCF=results/pilot_01/joint/beluga/cohort.raw.vcf.gz
BCF=results/pilot_01/joint/beluga/cohort.raw.bcf

test -s "$VCF" && test -s "${VCF}.csi"
bcftools index -n "$VCF"               # число записей вариантов
bcftools query -l "$VCF"               # sample IDs из заголовка
bcftools index -s "$VCF"               # contig и число записей по contigs
bcftools stats -s - "$VCF" > results/pilot_01/beluga.bcftools.stats.txt
bcftools index -n "$BCF"               # проверка доступности BCF через индекс
```

**Проверьте отдельно**: ожидаемое число и имена животных; отсутствие смешения видов; сборку/contig accession; непустой VCF и разумное покрытие. `cohort.raw` — исходный joint callset, **не** автоматически очищенная панель или доказательство качества вызова SNP.

## 9. Режим VCF → population structure

```bash
python pipeline.py --vcf /absolute/path/to/cohort.raw.bcf
# или индексированный .vcf.gz
```

Пока `src/minsnp/structure.py` — **явная заглушка**: проверяет наличие входа и индекса (`.csi`/`.tbi`), возвращает **код выхода 2** с сообщением `НЕ РЕАЛИЗОВАНО`, **не создаёт** PCA, clusters или `population_labels.tsv`. Опции `--submit`, `--resume` и SLURM backend для неё не поддерживаются. **Не ожидайте результатов этого режима до реализации и тестирования алгоритма.**

Для следующего этапа исследования метки должны быть отдельно получены/подтверждены по полногеномному набору SNP (QC, LD pruning, PCA, кластеризация/ancestry и проверка устойчивости), а затем переданы в третий режим как metadata.

## 10. Режим VCF + metadata → SNP panel

Этот режим пока не реализован.

```bash
python pipeline.py \
  --vcf /absolute/path/to/cohort.raw.bcf \
  --metadata /absolute/path/to/population_labels.tsv \
  --species polar_bear
```

Команда завершается `NotImplementedError`.

SNP-level QC, ранжирование SNP, LD filtering,
classification и CV пока не выполняются.

Выходные панели, модели и графики не создаются.

## 11. Ошибки, возобновление и воспроизводимость

### 11.1. Частые причины сбоев

| Симптом | Проверка / действие |
|---|---|
| `Не найден FASTQ` | Проверить путь относительно **папки TSV**; убедиться, что файл не пуст |
| `Референс не проиндексирован для BWA-MEM3` | Нужны `.amb/.ann/.pac/.bwt.2bit.64` именно от `bwa-mem3 index`; BWA-MEM индексов недостаточно |
| `Exec format error` / instruction set / AVX2 | Убедиться в поддержке AVX2 на **compute node** |
| `OUT_OF_MEMORY` или `TIMEOUT` | Посмотреть `sacct`/`.err`, скорректировать **соответствующую** `slurm.stages.<name>` в **новом** запуске |
| `DependencyNeverSatisfied` / следующая стадия не стартует | Найти FAILED/OOM/TIMEOUT у **предыдущего** stage array; `afterok` не разрешит продолжение |
| `Изменился код/TSV/конфиг` | Сверить `slurm_plan.json`, не менять файлы активного запуска; начать новый run при правках |
| `Каталог запуска уже существует` | Взять новый `--run-id`; не удалять прежние результаты без проверки |
| `Нет индекса VCF/BCF` | Построить `bcftools index -f` и проверить соответствие сборке |
| `sample_id ... не совпадают` | Сверить `bcftools query -l` и таблицу metadata; один ID = одно животное |
| `НЕ РЕАЛИЗОВАНО: ... structure` | Ожидаемое поведение заглушки; PCA/labels пока отсутствуют |

### 11.2. Resume и повторный запуск

- **SLURM:** после отправки сохраняется `slurm_jobs.json`; повторное `--submit` с тем же `run-id` **заблокировано**, даже если стадия упала. Если отправка оборвалась, `slurm_jobs.partial.json` показывает уже созданные Job IDs; **не отправляйте вручную те же задачи**, пока не выяснили их состояние. Для нового запуска используйте другой `run-id`.
- **Local WGS:** можно повторить ту же команду с `--resume`, если не менялись исходные TSV/FASTQ, код и конфиги. Уже валидированные результаты пропускаются, неполные могут перемещаться в файлы `.invalid_N` для диагностики.
- **Local panel:** `--resume` проверяет неизменность входов/кода/конфига; это **не** восстановление отдельного прерванного fold: аналитические вычисления могут повториться.
- **Старые запуски:** план с `schema_version: 1` и gVCF по `sample × contig` несовместим с новой схемой `schema_version: 3`, где один gVCF на образец и один совместный GenomicsDB/BCF на вид. Используйте новый `run-id`. Не пытайтесь продолжать планы версии 2: изменился набор стадий, выходные пути и состав задач.

Для воспроизводимости сохраняйте версию кода (Git commit), `slurm_plan.json`/`analysis_plan.json`, файлы конфигурации, данные по сборкам, `versions.json`, FASTQ/metadata checksums, используемые population labels и отчёты QC. **Не добавляйте крупные FASTQ/BAM/VCF/результаты в Git**.

## 12. Источники и технические ограничения

- [NCBI Datasets — установка CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) и [скачивание геномов](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/).
- [ENA — программный доступ к file reports](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html): URL/размер/MD5 публичных FASTQ.
- [GATK GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/27007914857627-GenomicsDBImport): импорт индивидуальных gVCF по sample-name-map.
- [`docs/PROJECT_GUIDE.md`](docs/PROJECT_GUIDE.md): архитектура, контракты данных и порядок задач.
- [`src/minsnp/README.md`](src/minsnp/README.md): состав Python-модулей.
- [`wgs_test/genome_reference_setup/README.md`](wgs_test/genome_reference_setup/README.md): происхождение старых reference-helper и предупреждение о BWA-MEM3.

**Не реализовано:** автоматическое определение популяций, независимая научная валидация SNP-панелей, автоматический restart упавших SLURM stages, общий лимит заданий между разными `run-id`. Полноценный FASTQ → VCF и VCF + metadata → панель необходимо протестировать на реальном пилоте перед использованием итогов в отчёте.
