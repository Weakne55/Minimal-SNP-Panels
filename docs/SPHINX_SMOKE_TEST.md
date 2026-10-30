# Проверка на sphinx без изменения исходного WGS/VCF

## 1. Окружение

```bash
conda activate wgs-gatk-mem3
cd /path/to/Minimal-SNP-Panels-greenfield
./run_pipeline.sh check-tools
pytest -q
```

Обратите внимание, что в переданном `conda list` отсутствуют PLINK2, scikit-learn и matplotlib; используйте инструкцию из README, не обновляйте GATK/bwa/samtools без причины. Проверить доступность индексов FASTA и BWA-MEM3, убедиться, что `bcftools +setGT -h` и `bcftools +fill-tags -h` работают.

## 2. Первый реальный тест: QC уже полученных двух белух

**Это не тест population structure**: два sample слишком мало для статистически осмысленного кластерного анализа и cross-validation.

```bash
VCF=/mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels/results/beluga_wgs_20261023/vcf/Delphinapterus_leucas/cohort.raw.vcf.gz
bcftools query -l "$VCF"
bcftools index -n "$VCF"
./run_pipeline.sh qc --vcf "$VCF" --species beluga \
    --run-id smoke_beluga_qc_v1 --outdir /mnt/tank/scratch/plem/data/Bogdan/minsnp_new_runs
```

Проверить: `qc/normalized.bcf`, `qc/filtered.bcf`, `.csi`, `bcftools.stats.txt`, `sample_qc.tsv`, `qc_summary.tsv`; число sample должно остаться 2; хромосомы только белухи; ошибки REF/ALT отсутствуют.

```bash
OUT=/mnt/tank/scratch/plem/data/Bogdan/minsnp_new_runs/beluga/smoke_beluga_qc_v1/qc/filtered.bcf
bcftools query -l "$OUT"
bcftools index -n "$OUT"
bcftools stats -s - "$OUT" | head -50
```

## 3. Тест отдельного joint по двум белухам

Используйте существующие `gVCF` без копирования или перезаписи, но создайте `gvcf_manifest.tsv` из реальных путей для 2 samples × 21 contigs. Проверить sample ID каждого gVCF по `bcftools query -l` и наличие индексов. Запустить `./run_pipeline.sh joint --gvcf-manifest ... --species beluga --run-id smoke_joint_v1 --outdir ...`, сравнить с baseline через `bcftools stats`, количество sample/contigs и генотипы на пересекающихся позициях. Пока сравнение не проведено, не заменять baseline.

## 4. Анализ Structure и Panels

Выбрать реальную когорту с достаточным числом независимых образцов и проверенными метками; не использовать двух белух как тест классификатора. Для nested CV требуется несколько образцов в каждой популяции и, при наличии родства, достаточно независимых групп. Сначала выполнить `structure` и проверить PCA, кластерную интерпретацию и batch effects, затем `panels`.

## 5. Проверка SLURM

Сначала успешно выполнить локальный QC из п.2; затем повторить на **новом run_id** с `--backend slurm`. Проверить `slurm_job_id.txt`, `squeue`, stdout, выходной индекс, сравнить число SNP с локальным результатом. Не подавать полномасштабную FASTQ когорту до этого.
