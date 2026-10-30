# Проверка FASTQ → VCF на sphinx

Документ описывает **проверки**, а не подтверждённые результаты реального
end-to-end запуска новой реализации. Не перезаписывать старые результаты.

## 1. Окружение и код

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
conda activate wgs-gatk-mem3
python -m pip install --no-deps -e .
./run_pipeline.sh --help
./run_pipeline.sh check-tools
python -m pytest -q
bash -n run_pipeline.sh stages/*.sh
```

## 2. Референсы

Проверить пути из `config/references.yaml`. Для BWA-MEM3 обязательно
`<FASTA>.bwt.2bit.64` наряду с `.amb`, `.ann`, `.pac`. По предыдущему
контрольному запуску индекса `.bwt.2bit.64` не было у обоих видов.
Индексировать на вычислительном узле отдельным заданием, не в login-сессии.

## 3. FASTQ manifest

Реальные FASTQ уже находятся в `data/wgs/Delphinapterus_leucas/` и
`data/wgs/Ursus_maritimus/`. Сформировать TSV в формате
`inputs/fastq.template.tsv`, проверить sample/run/library и пути.
**Сам шаблон использовать нельзя:** он содержит только заголовок.

## 4. Первое реальное испытание

Начать с реального ограниченного объёма чтений и отдельного `run-id`.
Возможность `wgs-vcf` и новый ресурсный SLURM-код пока прошли лишь
офлайн-тесты — реальный pipeline необходимо проверить.

```bash
./run_pipeline.sh wgs-vcf --fastq inputs/fastq.tsv --species beluga \
  --run-id smoke_beluga_new_001 --backend slurm --threads 4 --memory 64G
```

Эта команда предполагает, что `inputs/fastq.tsv` **создан** и индексы готовы.
Контролировать задание:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

После завершения проверить sample IDs и индексы:

```bash
RESULT=results/beluga/smoke_beluga_new_001/joint/cohort.raw.vcf.gz
bcftools query -l "$RESULT"
bcftools index -n "$RESULT"
bcftools stats "$RESULT" | head -50
```

## 5. Сравнение со старым результатом

Старый VCF двух белух находится в
`results/beluga_wgs_20261023/vcf/Delphinapterus_leucas/cohort.raw.vcf.gz`.
Не путать его с новыми путями `results/beluga/<run-id>/...`.
Сравнить sample IDs, сборку, число SNP и согласие генотипов на
совпадающих позициях. Старый VCF не перезаписывать.

## 6. Последующие исследования

Два образца не подходят для корректной population assignment CV.
Сначала проверить FASTQ → VCF и QC, затем подбирать достаточно большую,
независимую выборку с верифицированными population labels.
