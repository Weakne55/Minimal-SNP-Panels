# AGENTS — Minimal SNP Panels

## Контракты и приоритеты

Проект имеет **три** режима: `--fastq` (FASTQ→VCF), `--vcf` (structure, заглушка), `--vcf --metadata/--labels --species` (QC→SNP panel). Не выдавать фиктивные population labels и результаты анализа; полный end-to-end тест на сервере обязателен перед заявлением production readiness.

`pipeline.py` — одна CLI-точка входа. При `--backend slurm` на каждый этап отдельный array с `afterok` предыдущего; `slurm.max_concurrent_jobs` ограничивает все SLURM workers одного run-id. Ресурсы — в `slurm.stages.<name>`. HaplotypeCaller **один раз на sample для всех аутосом**; GenomicsDB/GenotypeGVCFs один вызов на вид по всем выбранным аутосомам; виды не объединять в общий joint.

## Безопасность

- Plan schema v3, контрольные суммы inputs/config/code; не обновлять активный код при запущенных задачах.
- Прерванные BAM/VCF нельзя считать готовыми только из-за наличия файла: проверять индекс и быстрые проверки целостности, сохранять подозрительные выходы как `.invalid_N`.
- Избегать общих writable outputs между параллельными workers.
- Не допускать повторной отправки SLURM jobs того же run-id, в том числе при частичной отправке.
- `--submit` обязателен для реальной отправки; без него для SLURM только dry-run.
- Полный ход исполнения сохранять в `logs/`, stdout/stderr GATK/samtools/bcftools — в operation log каждого stage.
- Отбор SNP и LD selection при cross-validation исключительно внутри training folds. Metadata — только фиксированные метки; `structure.py` заглушка.

## Процедура разработки

Менять адресно; дополнять тестами `tests/`; после исправления выполнять `compileall`, проверять Bash через `bash -n`, запускать `pytest`. Не коммитить референсы, FASTQ, BAM, gVCF, BCF и результаты анализа. Описание поведения обновлять в README и docs/PROJECT_GUIDE.md. Не считать unit tests заменой реального запуска.
