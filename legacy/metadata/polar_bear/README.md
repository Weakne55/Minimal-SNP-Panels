# Метаданные WGS белых медведей: поля и GitHub

## Рекомендованная структура репозитория

```text
data/metadata/polar_bear/
├── raw/
│   ├── polar_bear_WGS_metadata_final_144.tsv
│   ├── ncbi_sra_cohort_WGS_286.tsv
│   └── ncbi_sra_runs_all_390.tsv
├── curated/
│   ├── polar_bear_samples_working_144.tsv
│   └── polar_bear_runs_working_286.tsv
└── README.md
```


**Важно:** `published_group` и `geographic_candidate_group` — опубликованные/географические аннотации, а **не** `population_labels.tsv`. Генетические метки задаются отдельно после PCA/структуры.

Служебные поля нельзя удалять из *исходных* файлов без потери происхождения или контроля качества. Сокращённая копия не подменяет исходный каталог.

`ENA_raw_depth_estimate_x` — предварительная оценка из количества прочтений, не эквивалентная глубине после картирования. Сравнения глубины между исследованиями требуют единого pipeline.

В таблицах SRA нет прямых `fastq_ftp`/`fastq_md5`; для скачивания используйте SRR через SRA Toolkit или ранее подготовленный отдельный ENA FASTQ manifest.

## Классы полей в основной таблице (103)

### WORKING — В работе с образцами, группами, покрытием и допуском к анализу (30)

`biosample`, `animal_id_in_paper`, `high_coverage_id`, `run_accessions`, `n_runs`, `bioproject`, `primary_paper`, `published_mapped_depth_x`, `later_reanalysis_depth_x`, `ENA_raw_depth_estimate_x`, `depth_evidence_class`, `geographic_region`, `district`, `collection_locality`, `published_group`, `published_group_code`, `published_group_provenance`, `sex`, `genetic_affinity_exception`, `is_confirmed_cross_biosample_duplicate`, `duplicate_partner_biosample`, `possible_alias_duplicate_group`, `QC_flags`, `geographic_candidate_group`, `geographic_candidate_confidence`, `metadata_review_priority`, `metadata_review_reason`, `sample_tissue_published`, `ncbi_sra_instruments`, `metadata_integrity_review`

### TECHNICAL_QC — Технические показатели из публикаций и баз; нужны для контроля/аудита, но не для финальных панелей (31)

`cahill_S2_locality_mapped_depth_x`, `raw_depth_vs_mapped_depth_ratio`, `raw_vs_mapped_depth_review`, `miller_published_read_pairs_or_fragments`, `miller_published_bases_Gb`, `miller_published_read_layout`, `liu_raw_data_Gb`, `liu_mapped_data_Gb`, `liu_mapped_reads_percent`, `liu_genome_coverage_percent`, `laidre_heterozygosity`, `laidre_ROH_gt2Mb_value_reported`, `ena_fastq_files_count`, `ena_fastq_compressed_bytes`, `ena_fastq_compressed_GiB`, `ncbi_sra_experiments`, `ncbi_sra_library_selection`, `ncbi_sra_library_names`, `ncbi_sra_sequencing_centers`, `ncbi_sra_total_bases`, `ncbi_sra_total_download_bytes`, `ncbi_sra_tissue`, `ncbi_sra_development_stage`, `ncbi_sra_collection_date`, `ncbi_sra_additional_runs_not_in_wgs_manifest`, `ncbi_sra_additional_run_types`, `lan2022_mapped_depth_x`, `lan2022_aligned_reads_percent`, `lan2022_deduplicated_aligned_reads`, `lan2022_pcr_duplicate_percent`, `lan2022_reference_covered_fraction`

### AUDIT_SOURCE — Ссылки, источники, альтернативные идентификаторы, детали сопоставлений (31)

`doi`, `published_depth_source`, `source_paper_url`, `biosample_url`, `notes`, `sex_evidence`, `cahill_S2_candidate_depth_note`, `geographic_candidate_basis`, `geographic_candidate_source_url`, `extra_source_urls`, `source_verified_sex`, `sex_source_ref`, `liu_collection_id`, `liu_sample_institution`, `liu_S2_row`, `laidre_age_class`, `laidre_S7_sample_id`, `cahill_candidate_collection_accession`, `cahill_candidate_accession_status`, `source_supplementary_new`, `notes_evidence_verified`, `ncbi_sra_sample_name`, `ncbi_sra_geographic_location`, `ncbi_sra_sex`, `lan2022_sample_id`, `lan2022_collection_id`, `lan2022_locality`, `lan2022_sex`, `lan2022_population_broad_code`, `lan2022_mitogenome_accession`, `lan2022_evidence`

### REDUNDANT_OR_EMPTY — Повторяющиеся, пустые и константные поля, не требуются в сокращённой таблице (11)

`normalized_sex`, `geographic_label_usable_for_training`, `liu_latitude`, `liu_longitude`, `ena_manifest_complete`, `public_metadata_status`, `ncbi_sra_assay_type`, `ncbi_sra_library_layout`, `ncbi_sra_platform`, `ncbi_sra_metadata_evidence`, `wgs_cohort_run_status`

## Классы полей run manifest

**Рабочие (14):** `run_accession`, `biosample`, `bioproject`, `experiment`, `assay_type`, `organism`, `library_layout`, `library_selection`, `library_source`, `library_name`, `instrument`, `platform`, `bases`, `bytes`

**Аудит (7):** `sample_name`, `geo_loc_name`, `sex`, `tissue`, `dev_stage`, `collection_date`, `ncbi_source_csv`

**Контроль отбора (3):** `original_wgs_manifest`, `in_144_biosample`, `exclusion_reason`

## Что не включать в классификацию подпопуляций

SRA технические поля (`bases`, `bytes`, `instrument`, `library_*`), географические кандидаты, источник публикации и само имя BioSample не должны становиться признаками классификатора. Модель обучается на SNP-матрице. Био- и технические метаданные нужны для QC, разбиения, проверки возможного batch effect и интерпретации.
