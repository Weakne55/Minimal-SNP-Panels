# Модули `src/minsnp`

Точка входа — `pipeline.py` в корне репозитория.

| Модуль | Назначение |
| --- | --- |
| `core.py` | Ввод/вывод, референсы, запуск команд и логирование |
| `resources.py` | Порядок стадий, CPU/RAM, лимиты SLURM |
| `slurm_wgs.py` | SLURM arrays, зависимости, логирование, status |
| `wgs.py` | FASTQ → BAM → один gVCF на образец по всем аутосомам |
| `joint.py` | GenomicsDB → GenotypeGVCFs по всем аутосомам вида → BCF/VCF |
| `structure.py` | Заглушка VCF → population structure |
| `panels.py` | Заглушка VCF + metadata → SNP panel |

Аналитический QC вариантов и отбор SNP-панелей пока не реализованы.

`bam_qc.sh` сохранён: он необходим для технического QC BAM.
