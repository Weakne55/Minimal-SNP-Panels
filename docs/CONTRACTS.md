# Контракты входных и выходных данных

## Пользовательские входы

**FASTQ (`--fastq`)** — TSV с обязательным заголовком:

```text
sample_id	species	run_id	library	fastq_r1	fastq_r2
```

Каждая строка — одна paired-end FASTQ-пара; у образца допустимы несколько
`run_id`/`library`. `species` — `beluga` либо `polar_bear`.
Пути к FASTQ могут быть абсолютными или **относительными к каталогу TSV**.
Повторное использование одного FASTQ и дубликаты комбинаций
`sample_id`/`run_id`/`library` недопустимы. Manifest должен содержать
**реальные** файлы, а не только заголовок шаблона.

**VCF (`--vcf`)** — индексированный multi-sample `VCF.gz` или `BCF`
из той же референсной сборки, что и `--species`.

**Метки (`--labels`)** — TSV с `sample_id`, `population`; дополнительный
`group` используют для группового разделения обучающей и тестовой выборок.
Географические поля из `data/metadata/WGS_final.xlsx` сами по себе не
являются подтверждёнными генетическими population labels.

## Внутренние контракты WGS → joint

`gvcf_manifest.tsv`: `sample_id`, `contig`, `gvcf_path`.
Ожидается полная матрица sample × contig, однообразцовые индексированные
хромосомные gVCF с одинаковой сборкой и корректными sample IDs.

## Выходы

- `wgs/lanes/` — BAM по sequencing runs;
- `wgs/bam/` — merged и markdup BAM, индексы;
- `wgs/qc/` — BAM QC и duplicate metrics;
- `wgs/gvcf/` и `wgs/gvcf_manifest.tsv` — gVCF для joint;
- `joint/genomicsdb/`, `joint/per_contig/` — промежуточный joint calling;
- `joint/cohort.raw.bcf` и `.csi` — индексированный multi-sample BCF;
- `joint/cohort.raw.vcf.gz` и `.csi` — экспорт VCF в режиме `wgs-vcf`;
- `qc/` — нормализованные и фильтрованные SNP (отдельная команда `qc` или `run`);
- `structure/`, `panels/` — последующие исследовательские стадии.

Выходы сохраняются в `results/<species>/<run-id>/` (по умолчанию).
Не менять или удалять старые `results/beluga_wgs_20261023/` при тестировании.
