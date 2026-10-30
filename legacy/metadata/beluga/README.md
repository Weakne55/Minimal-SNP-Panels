# Метаданные WGS белух: поля и GitHub

В каталоге собраны метаданные 106 BioSample (*Delphinapterus leucas*) из публичных проектов **PRJNA1182041** (78) и **PRJNA1328577** (28). Все 106 соответствуют 106 SRA Run (SRR). Здесь **нет** FASTQ, BAM, VCF и низкопокрытийной выборки Montana et al. (2024).

## Структура

```text
data/metadata/beluga/
├── README.md
├── metadata_column_classification.tsv
├── raw/
│   ├── beluga_WGS_metadata_final_106.tsv
│   └── ncbi_sra_runs_all_106.tsv
└── curated/
    ├── beluga_samples_working_106.tsv
    └── beluga_runs_working_106.tsv
```

Как в `polar_bear/`, `raw/` хранит полные таблицы с техническими/аудитными полями, а `curated/` — сокращённые рабочие таблицы образцов и запусков. `metadata_column_classification.tsv` описывает поля двух исходных TSV и их роль. Изменения относительно полного предыдущего комплекта: не включены отдельные сводки и выборки по порогам покрытия, поскольку они выводятся фильтрами из `beluga_samples_working_106.tsv`.

В отличие от медведей, для белух выгрузки `ncbi_sra_cohort_WGS_106.tsv` и `ncbi_sra_runs_all_106.tsv` **побайтово совпадают**: в полученных двух SRA-экспортах нет дополнительных запусков, исключённых из когорты. Поэтому оставлена только одна таблица SRA, без дублирования данных.

## Состав и правила интерпретации

| BioProject | BioSample | SRA Run | Assay Type в SRA | Возраст |
|---|---:|---:|---|---|
| PRJNA1182041 | 78 | 78 | `OTHER` | 78 взрослых |
| PRJNA1328577 | 28 | 28 | `WGS` | 20 взрослых, 8 новорождённых |
| **Итого** | **106** | **106** | | **98 взрослых, 8 новорождённых** |

- В `PRJNA1182041` SRA указывает `Assay Type = OTHER`, несмотря на опубликованное полногеномное секвенирование, и `LibrarySource = GENOMIC`. Не исключайте все 78 по одному полю `Assay Type`.
- `estimated_raw_depth_x` — оценка по `Bases / 2 350 000 000`, а **не фактическое среднее покрытие выровненных BAM**. `published_mapped_depth_x` пуст, если индивидуального опубликованного измерения нет.
- Сырым покрытием ≥10× обладают 104 образца, ≥15× — 98, ≥20× — 60. Среди взрослых ≥15× — 90, ≥20× — 52. Для получения подвыборок фильтруйте `curated/beluga_samples_working_106.tsv` по `development_stage`, `raw_depth_ge15` или `raw_depth_ge20`.
- `collection_locality_sra`, `ecotype_sra`, `geographic_group_from_metadata` — метаданные о происхождении, **не** подтверждённые генетические популяции. `genetic_population_label` оставлено пустым до проверки структуры. Для 10 записей локализация не указана (`Canada: NA`).
- Новорождённых учитывайте отдельно от взрослых при формировании обучающих выборок. Порог по `estimated_raw_depth_x` — предварительный; окончательный отбор после mapping/QC.

## Основные поля рабочих таблиц

- **Образцы:** `biosample`, `sample_id`, `run_accessions`, `n_runs`, `bioproject`, `primary_paper`, `development_stage`, `sex`, `estimated_raw_depth_x`, `raw_depth_ge15`, `raw_depth_ge20`, `collection_locality_sra`, `ecotype_sra`, `genetic_population_label`, `training_preliminary_status`.
- **Запуски:** `run_accession`, `biosample`, `bioproject`, `experiment`, `assay_type`, `library_layout`, `library_source`, `instrument`, `bases`, `bytes`.

## Источники

- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1182041 — de Greef et al. (2026), https://doi.org/10.1086/738889
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1328577 — Orton et al. (2026), https://doi.org/10.1111/mec.70521
- Таблицы собраны из предоставленных экспортов NCBI SRA Run Selector и географических метаданных.

`population_labels.tsv` не создаётся из этих метаданных автоматически: для него требуется отдельная проверка популяционной структуры.
