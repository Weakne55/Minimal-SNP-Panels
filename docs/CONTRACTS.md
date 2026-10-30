# Контракты входов и выходов

## Пользовательские входы

* `--fastq fastq.tsv`: обязательные столбцы `sample_id,species,run_id,library,fastq_r1,fastq_r2` (разделитель **TAB**). Одна строка — одна paired-end FASTQ-пара/секвенирующий запуск. Допустимы несколько `run_id` и `library` для одного `sample_id`. Все пути — абсолютные; идентификаторы должны совпадать со значениями SM в VCF.
* `--vcf cohort.vcf.gz` либо `--vcf cohort.bcf`: **индексированный** multi-sample VCF/BCF. Исходный VCF должен относиться к сборке, выбранной по `--species`; совпадение только названий хромосом недостаточно, нужна референсная совместимость.
* `--labels labels.tsv`: два обязательных столбца `sample_id,population`; необязательный `group` (известные родственники/повторения). Строго один label на sample; набор sample IDs обязан совпадать с VCF.

## Внутренняя граница WGS → joint

`gvcf_manifest.tsv`: `sample_id,contig,gvcf_path`. Для каждого sample есть **все** contigs из references.yaml; gVCF однообразцовый и индексированный; sample ID проверяется по заголовку.

## Выходы

* `wgs/gvcf/<sample>/<contig>.g.vcf.gz` + индексы, `wgs/bam/`, `wgs/qc/`, `wgs/gvcf_manifest.tsv`.
* `joint/cohort.raw.bcf` и `.csi`.
* `qc/normalized.bcf`, `qc/filtered.bcf`, `qc/bcftools.stats.txt`, `qc/qc_summary.tsv`.
* `structure/pca.tsv`, `structure/pca.png`, `structure/population_labels.tsv`, `structure/cluster_diagnostics.tsv`.
* `panels/cv_metrics.tsv`, `panels/nested_cv_metrics.tsv` (при достаточном N), `panels/cv_predictions.tsv`, `panels/fold_selections.tsv`, `panels/minimal_panel.tsv`, `panels/model.joblib`, рисунки.

Промежуточные результаты никогда не записываются в исходный репозиторий; каждый `run_id` создаёт отдельную папку.
