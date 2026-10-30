# CPU и память для FASTQ → VCF

Параметры `--threads` и `--memory` независимы. RAM нельзя определять только по числу CPU: размер BWA-MEM3 индекса зависит от референса.

## Usage

```bash
# Один CPU, достаточно памяти для индекса генома
./run_pipeline.sh wgs-vcf --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_cpu1 --backend slurm --threads 1 --memory 48G

# 32 CPU, отдельный RAM-бюджет
./run_pipeline.sh wgs-vcf --fastq inputs/fastq.tsv --species beluga \
  --run-id beluga_cpu32 --backend slurm --threads 32 --memory 96G
```

Без `--threads` и `--memory` используются `threads` и `slurm.memory` из `config/pipeline.yaml`. Внутри `salloc`/`sbatch` при отсутствии `--threads` используется `SLURM_CPUS_PER_TASK`.

## Распределение ресурсов

- `align.sh`: одновременно запускает BWA-MEM3, samtools collate, fixmate и sort. CPU распределяются между BWA и sort; при 1 CPU однопоточные процессы делят выделенный процессор. `samtools sort -m` ограничивается по RAM на поток.
- `merge_bam.sh`, `markdup.sh`, `bam_qc.sh`: samtools получает максимум `CPU - 1` дополнительных рабочих потоков, но не больше 32 (`-@ 0` при CPU=1).
- `trim.sh`: fastp не использует более 16 потоков, даже если выделено больше CPU.
- `haplotypecaller.sh`: Java heap не больше 16 GiB, PairHMM ограничен 8 потоками.
- `genomicsdb.sh`: heap до 24 GiB, 8 reader threads максимум, batch-size до 50; остаётся RAM под native libraries.
- `genotype.sh`: heap до 16 GiB; JVM учитывает выделенное число CPU.
- `concat.sh`, `export_vcf.sh`: bcftools получает до 16 дополнительных потоков для сжатия.

Ограничения 8/16 потоков являются ограничениями конкретных стадий, **не** ограничением допустимого `--threads`. При больших CPU-бюджетах отдельные стадии не обязаны загружать все выделенные процессоры.

## Проверка Slurm

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

Запрос 1 CPU не означает, что весь FASTQ → VCF запустится с любой памятью. Индекс референса и объём данных задают минимальный фактический RAM-бюджет. При недостатке RAM Slurm может завершить job с `OUT_OF_MEMORY`.

Создаётся один Slurm job на весь запуск: оптимизация job arrays / отдельных стадий в эту правку не входит. Для нового набора ресурсов используйте новый `--run-id`.

## Граница управления ресурсами

`--threads N` допускает `N >= 1`; `--memory` задаётся независимо от
CPU. Внутри одного SLURM job стадии выполняются **последовательно**:
выделение 64 CPU не означает, что каждая стадия действительно использует 64.
BWA-MEM3 и samtools используют доли CPU-бюджета; Java heap ограничен
отдельно от native memory. Конфигурация времени и partition — в
`config/pipeline.yaml`, блок `slurm`.

Если `--threads` или `--memory` изменились, используйте новый `run-id`,
а не `--resume` старого запуска. Прежде чем запускать полную выборку,
сверьте `ReqMem`, `AllocCPUS`, `MaxRSS` через `sacct`.
