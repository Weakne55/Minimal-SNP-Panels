# Minimal SNP Panels — технический справочник


## 1. Назначение и границы анализа

**Виды:** белуха (*Delphinapterus leucas*) и белый медведь (*Ursus maritimus*). **Дедлайн:** 25.12.2026.

- **WGS pipeline:** paired-end FASTQ → optional fastp → BWA-MEM3 → samtools (fixmate/sort/merge/markdup) → BAM QC → GATK HaplotypeCaller `-ERC GVCF` по аутосомам → gVCF manifest.
- **Joint:** GATK GenomicsDBImport → GenotypeGVCFs по аутосомам → `bcftools concat` → multi-sample `cohort.raw.bcf`, при `wgs-vcf` дополнительно индексированный `cohort.raw.vcf.gz`.
- **Основной pipeline:** индексированный multi-sample VCF.gz/BCF + `population_labels.tsv` → variant/sample QC → SNP ranking и LD filtering → панели разных размеров → classifier + CV и отчётные таблицы.
- **Исследовательская ветка:** LD-pruned SNP → PCA / clustering → предварительные population labels. Её результаты требуют независимой проверки, а вход готовых labels не должен зависеть от реализации PCA.

`wgs.py` обрабатывает животных отдельно; `joint.py` совместно генотипирует когорту. Нельзя заменить joint genotyping простым `bcftools merge` отдельных окончательных VCF. `bcftools concat` в этом проекте соединяет **непересекающиеся участки генома**, а не особей.

## 2. Актуальное состояние на момент снимка

| Показатель | Подтверждённое / ограничение |
|---|---|
| Ветка | `develop`  |
| Инвентаризация WGS | 109 белух (121 run) + 144 медведя (287 runs) = **253 WGS-записи / 408 runs**; возможны дубли животных |
| Пилот `pilot_20261030` | **20 животных** (10+10), **36 runs**, около 585 GiB FASTQ.gz |
| Загрузка | **36/36 runs подтверждены** размером и MD5 по предыдущим логам; физические FASTQ в Git отсутствуют |
| Первые gVCF manifests | `SAMN51337163` (21 аутосома), `SAMN16454148` (36 аутосом) |
| Joint pilot 10+10 | **Не подтверждён**; совместные BCF на момент снимка не готовы |
| Старый результат | `results/beluga_wgs_20261023/.../cohort.raw.vcf.gz`: 2 белухи, отдельный запуск, не текущий пилот |

**12 собственных карских медведей** предусмотрены ТЗ, но их фактическое включение в пилот не подтверждено. Не публиковать закрытые образцы или их конфиденциальные метаданные.

### Референсы и окружение

| Вид | Assembly | Аутосом | Референсная FASTA (серверный путь) |
|---|---|---:|---|
| Белуха | `GCA_059928025.1` (BTC_Dleu_Yukka2) | 21 | `/mnt/tank/scratch/plem/data/genomes/Delphinapterus_leucas_GCA_059928025.1/ncbi_dataset/data/GCA_059928025.1/GCA_059928025.1_BTC_Dleu_Yukka2_genomic.fna` |
| Медведь | `GCA_062004589.1` (BTC_Umar_Aika3) | 36 | `/mnt/tank/scratch/plem/data/genomes/Ursus_maritimus_GCA_062004589.1/ncbi_dataset/data/GCA_062004589.1/GCA_062004589.1_BTC_Umar_Aika3_genomic.fna` |

Источник истины для contigs — `config/references.yaml`; нельзя предполагать непрерывность номеров `CM*`. Используются **только аутосомы**. Для BWA-MEM3 обязательна индексная часть `.bwt.2bit.64` наряду с `.amb`, `.ann`, `.pac`. Индексы были успешно подготовлены заданиями `943426`, `943427`, переиндексация без причины не требуется.

На `sphinx` окружение `wgs-gatk-mem3` (Python 3.11.18, BWA-MEM3 0.15.1, samtools/bcftools 1.24, GATK 4.7.0.0, SLURM 21.08.5). BWA-MEM3 требует AVX2, который отсутствует на части узлов (`kraken`); для актуальных WGS задач выбран `meduza-1`. `sacct` недоступен из-за `slurmdbd`.

## 3. Структура файлов и режимы исполнения

```text
run_pipeline.sh           пользовательский CLI
src/minsnp/               cli.py, core.py, wgs.py, joint.py, qc.py,
                          structure.py, panels.py, reporting.py, resources.py
stages/                   Bash-стадии BWA/samtools/GATK/bcftools
config/                   pipeline.yaml и references.yaml
inputs/                   FASTQ, gVCF, labels TSV-шаблоны/примеры
data/metadata/            WGS_final.xlsx + README_metadata.md
data/wgs/                 серверные реальные FASTQ, исключены из Git
results/                  серверные BAM/gVCF/GenomicsDB/VCF/логи, исключены
wgs_test/                 инструменты подготовки WGS/референсов и pilot launcher
literature/               literature_methods.tsv — отдельная рабочая таблица
legacy/                   локальная архивная реализация, не в публичном Git
tests/                    локальные тесты, по принятому решению исключены из Git
docs/PROJECT_GUIDE.md     единая техническая документация
```

**Важная совместимость:** уже отправленные SLURM scripts вызывают `.../tools/wgs_pilot.py`, поэтому серверная локальная копия `tools/` должна существовать, даже когда в Git публикуется `wgs_test/`. Переименование каталогов Git само по себе не изменяет эти сгенерированные скрипты. Не запускать `update_docs_inplace.py` или `configure_wgs_slurm.py` в режиме перезаписи документации: их шаблоны содержат старые пути и могут восстановить удалённые MD.

`./run_pipeline.sh` вызывает `minsnp.cli`, принимает команды `run`, `wgs`, `wgs-vcf`, `joint`, `qc`, `structure`, `panels`, `check-tools`. Режим `--backend slurm` запускает **один SLURM job на вызов CLI**; pilot launcher — отдельный внешний DAG с массивами задач, внутри которых CLI использует `--backend local`. Нельзя смешивать эти две схемы и считать, что каждый CLI вызов создаёт массивы.

## 4. Входы и выходы

### FASTQ: `--fastq`

Обязательный TSV-заголовок (tab-separated):

```tsv
sample_id	species	run_id	library	fastq_r1	fastq_r2
```

Одна строка соответствует одной R1/R2 паре. У одного sample допускаются несколько runs, в том числе нескольких библиотек. `species`: `beluga` либо `polar_bear`. Пути FASTQ могут быть абсолютными или **относительными к каталогу самого TSV**. Скрипт проверяет существование обоих FASTQ и запрещает повторную пару `sample_id/run_id/library` и повторное использование одного файла. Шаблон из одного заголовка не является валидным манифестом.

Для текущих данных подготовлены примеры `inputs/beluga_pilot_10.fastq.tsv` (10 samples/20 runs) и `inputs/polar_bear_pilot_10.fastq.tsv` (10 samples/16 runs); они предполагают серверную директорию `data/wgs/`. **Не запускать повторный WGS по этим манифестам, пока активен пилот.**

### GVCF: `--gvcf-manifest`

```tsv
sample_id	contig	gvcf_path
```

Ожидается полная матрица sample × autosome, уникальные sample IDs, один sample в каждом gVCF, корректные индексы и соответствие contigs/длин референсу. Путь к gVCF должен существовать. Файл `wgs/gvcf_manifest.tsv` сам по себе ещё не гарантирует, что все gVCF физически корректны: после завершения задачи проверять индексы, логи и валидаторы.

### VCF/BCF: `--vcf`

Индексированный `VCF.gz` либо `BCF` с минимум двумя образцами, один вид/сборка. Перед population assignment нужен отфильтрованный набор SNP и независимый контроль sample IDs. Нельзя смешивать BCF разных сборок без явного liftover/повторного calling.

### Population labels: `--labels`

```tsv
sample_id	population	group
```

`group` опционален; годится для групповых разбиений, если применимо. Указываемые классы должны быть обоснованы; географические и опубликованные метки в Excel — **не** автоматически подтверждённые генетические группы. Обязательно документировать случаи неоднозначной принадлежности и повторных BioSample одного животного.

### Выходы

| Путь внутри `results/<species>/<run-id>/` | Содержимое |
|---|---|
| `wgs/lanes/` | `fixed.bam` по run/library |
| `wgs/bam/`, `wgs/qc/` | объединённые и markdup BAM, индексы, BAM QC |
| `wgs/gvcf/`, `wgs/gvcf_manifest.tsv` | gVCF по sample × autosome и manifest |
| `joint/genomicsdb/`, `joint/per_contig/` | промежуточное совместное генотипирование |
| `joint/cohort.raw.bcf` + `.csi` | итоговый многосэмпловый BCF |
| `joint/cohort.raw.vcf.gz` + индекс | экспорт `wgs-vcf` или оркестратора пилота |
| `qc/`, `structure/`, `panels/` | downstream QC, исследовательская структура и SNP panels (не считать завершёнными до проверки) |
| `logs/`, `versions.json` | логи и зафиксированные версии инструментов |

`--resume` использовать лишь для того же запуска с идентичными исходными файлами и параметрами. Замена референса/manifest/ресурсов требует нового `run-id` и отдельной валидации. Для уже отправленного пилота запрещено повторно запускать `submit --confirm`.

## 5. Запуск и безопасный мониторинг

Проверка окружения без запуска анализа:

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
conda activate wgs-gatk-mem3
./run_pipeline.sh check-tools
bash -n run_pipeline.sh stages/*.sh
python -m compileall -q src
# tests/ локальные и могут отсутствовать в новом clone
PYTHONPATH=src python -m pytest -q
```

**Диагностика уже работающего пилота, только чтение:**

```bash
PILOT=results/pilot_orchestration/pilot_20261030
squeue -r -u "$USER" -o '%.20i %.25j %.12T %.15M %.25R'
python tools/wgs_pilot.py status --run-id pilot_20261030
find "$PILOT/verified" -maxdepth 1 -name '*.json' | wc -l
find results/beluga results/polar_bear \
  -path '*/pilot_20261030_*/wgs/gvcf_manifest.tsv' -type f -print
scontrol show job 943429
scontrol show job 943431
df -h /mnt/tank/scratch/plem/data
```

Используйте `scontrol show job JOBID` для реального `StdOut`/`StdErr`, а не придумывайте имя лога. `squeue` не показывает завершённые задачи: отсутствие job из очереди само по себе не доказывает `COMPLETED`. Доступные признаки — завершение stdout, manifest, количество gVCF/индексов, текущая последующая задача и собственные marker files. На этом кластере `sacct` часто выдаёт `Unable to connect to database`.

Первый WGS beluga: `results/beluga/pilot_20261030_SAMN51337163/wgs/gvcf_manifest.tsv` (21); первый WGS bear: `results/polar_bear/pilot_20261030_SAMN16454148/wgs/gvcf_manifest.tsv` (36). Проверки BAM можно делать на готовых файлах через `samtools quickcheck -v`, `samtools flagstat`, `samtools idxstats`; `quickcheck` **не** равен биологическому QC. Полный QC по BAM в момент нагрузки на NAS выполнять умеренно.

До запуска joint по 10 образцам должны завершиться все 10 sample WGS tasks соответствующего массива. Подтверждать итоговый VCF/BCF командами `bcftools query -l`, `bcftools index -n`, `bcftools stats` и контролем числа samples (=10 на вид), SNP/indel counts, состава аутосом и ненулевого размера. Это **проверки после завершения**, не команды повторного calling.

### Пилот: DAG и история восстановления

```text
36 ENA/SRA runs ──> verification (size+MD5) ──> FASTQ sample manifests
                 ├── BWA-MEM3 index beluga  (943426) ──> WGS array 943428 (10, %1) ──> joint 943429
                 └── BWA-MEM3 index bear    (943427) ──> WGS array 943430 (10, %1) ──> joint 943431
```

- Изначальные index jobs `943413`/`943414` упали из-за отсутствия AVX2; переиндексация на `meduza-1` успешно завершилась.
- Первоначальные загрузки имели прерывания `curl` (`unexpected EOF`, эпизодический HTTP 403). В `tools/wgs_pilot.py` использованы resume `.part` и явные повторные попытки с обязательной проверкой size/MD5. Восстановление завершилось проверкой **36/36 runs**.
- Из-за частичных сбоев загрузки зависимость оригинального manifest `943415` не выполнилась; заменён на `943506`, после чего WGS arrays стартовали. Исходные задания **не** нужно повторно отправлять.
- На снимке `943428_1` и `943430_1` RUNNING, `943429`/`943431` PENDING (dependency); это исторические сведения, а не текущая очередь.
- Не запускать `recover_pilot_avx2.sh` и `wgs_pilot.py submit --confirm` повторно. При сбое исправлять конкретный sample/contig, не повторять весь DAG.

### CPU, RAM, storage

В CLI `--threads` и `--memory` независимы; WGS в текущем `wgs.py` требует **threads >= 2**. Память BWA-MEM3 зависит от размера индекса референса, не только от CPU. В pilot используется лимит 32 CPU / 96 GiB на WGS task; индексация потребовала около 125 GiB peak RSS. `align.sh` распределяет CPU между BWA и samtools, `fastp` ограничен 16 потоками, HaplotypeCaller/GenomicsDB ограничивают Java heap, поэтому увеличение CPU не гарантирует ускорения всех стадий.

Первый BAM белухи около **21.71 GB**, медведя **25.55 GB**, при этом lane BAM тоже занимает десятки GB; gVCF дополняют расход. FASTQ пилота ~585 GiB, предварительный резерв для результатов ~1.5–3 TiB (оценка, не измерение). **Ничего не удалять из активно обрабатываемого каталога**. Для полной когорты ~6.70 TiB сжатых FASTQ; требования NAS уточнять по завершённому пилоту.

## 6. Variant QC и анализ структуры — требования

Параметры по умолчанию в `config/pipeline.yaml`: `min_gq=20`, `min_dp=5`, `max_missing=0.2`, `min_maf=0.01`. Это **исходные настройки для проверки, не утверждённые QC-пороги**. Особенно важно исследовать медведей с оценённой глубиной ниже 2×: фильтрация GQ/DP может превратить их генотипы в missing и исказить структуру. Отдельно оценить mapping/depth, missingness, heterozygosity, allele balance, batch/BioProject, родственность, дубликаты и размер класса.

Для PCA использовать LD-pruned SNP. `structure.py` концептуально содержит PLINK2 PCA и кластеризацию KMeans/silhouette, но на снимке 02.11 есть блокирующая ошибка первичного вызова: нельзя писать, что модуль прошёл real-data анализ, до исправления и тестирования. Кластеры должны быть проверены независимым методом, устойчивостью числа групп и внешними данными; SNP panel pipeline должен принимать уже зафиксированные метки независимо от PCA.

## 7. SNP-panel selection и валидация

В `panels.py` начально реализованы Nei-подобная величина `(Ht-Hs)/Ht`, LD selection по `r²`, логистическая регрессия, CV и в доступных условиях nested CV. **Это не Weir–Cockerham FST**. Размеры панели в конфигурации: **10, 20, 50, 100, 200 SNP**; они являются кандидатами, не заранее доказанным минимумом. Предпочтительные выходы: ranked SNP, allele frequencies / population contrasts, candidate panels, per-class recall/precision, balanced accuracy, macro-F1, confusion matrices, зависимость метрики от размера панели, оценка устойчивости SNP.

Обязательные условия:

1. SNP selection, train-dependent QC/статистики, LD filtering, fitted transformations, classifier tuning и выбор размера панели должны быть внутри training folds. Для выбора размера использовать внутренние folds, качество оценивать на внешних.
2. Избегать разделения разных runs/BioSample одного животного по разным folds; где возможно применять grouped splits по животному/родству.
3. Для маленьких классов nested CV может быть невозможна; нельзя выдавать успешный `pytest` за достаточность N.
4. Кластеры, определённые по тем же WGS, дают метрики воспроизведения этих групп, **не** независимую точность географической принадлежности.
5. Отсутствуют подтверждённая независимая валидация, реализация опубликованных SNP-панелей и полная поддержка likelihood-based low-pass genotype calling.

`literature/literature_methods.tsv` поддерживается отдельно и **не меняется** при реорганизации Markdown. Сравнивать воспроизводимые фиксированные SNP-панели по координатам, аллелям и сборке либо документировать невозможность сравнения.


## 8. GitHub, воспроизводимость и безопасность

Публичный GitHub должен содержать исходный код, конфигурации и только разрешённую инвентаризацию (`data/metadata/`), а также примеры входных TSV. **Не публиковать** FASTQ/BAM/gVCF/GenomicsDB/results, локальный `legacy/`, закрытые данные, токены и ключи. Ветка разработки — `develop`; `main` защищать от случайных force-push и удаления. Контрольная дата релиза: 25.12.2026. Не добавлять лицензию автоматически до решения прав на код и данных.

Для каждого результата фиксировать commit, `run-id`, исходные accession/manifest, сборку/версии, параметры QC, число животных и SNP, локальные пути результатов и limitations. Рабочий файл `literature_methods.tsv` использовать как источник для будущих Introduction/Methods/Discussion.

**Документация была обновлена без изменения кода и без отправки задач на HPC.** Сведения о выполнении в любой последующий день следует перепроверять на `sphinx`.
