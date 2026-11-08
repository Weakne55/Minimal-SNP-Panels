# Рефакторинг — 2026-11-08

Основа: архив `Minimal-SNP-Panels-all-autosomes.zip` (schema v3).
Цель: сократить дублирование и привести документацию к фактическому коду,
**не меняя** три режима CLI и схему вызова GATK по всем аутосомам.

## Журнал изменений и причин

| Действие | Причина | Проверка |
|---|---|---|
| Единые `wgs.FASTQ_COLUMNS` и `load_cohort()` для полного и sample TSV | Устранить два независимых валидатора и несовпадение их правил | `test_fastq_validation_shared_between_cohort_and_sample` |
| Разделение WGS на `_run_lane()`, `_run_sample()` и `_run_checked()` | Упростить ветвление, унифицировать обработку некорректных результатов | `test_hc_worker_runs_once_and_writes_manifest`, `test_prepare_plan_and_resume_with_multiple_lanes` |
| `resources.LANE_STAGES`, `SAMPLE_STAGES`, `JOINT_STAGES` → `STAGE_ORDER` | Одна точка определения порядка и уровня гранулярности SLURM-стадий | `test_build_tasks_has_one_hc_per_sample_and_joint_per_species` |
| Удалены неиспользуемые `refs/settings` из `make_tasks()` | Интерфейс должен содержать только необходимые аргументы | Тесты плана и scheduler |
| `.submission.lock` вокруг `sbatch` | Не допустить конкурентной отправки arrays одного run-id | `test_concurrent_submission_lock_blocks_duplicate_sbatch` |
| `panels.py`: раздельные функции для fixed CV, nested CV и финальной модели | Уменьшить смешение расчётов и сохранения файлов без изменения алгоритмов | `test_panel_fold_selection_and_outputs` |
| Консистентные docstrings и комментарии; упрощён `core.call()` | Описывать назначение, контракты и причины, а не повторять код | `compileall` и тесты |
| Dry run FASTQ дополнительно проверяет доступность референсов/аутосом | Не выдавать успешный dry run при неправильном FASTA | Тесты CLI / проверка на примере |
| Удалён неиспользуемый `reporting.plot_structure()` | Анализ структуры — заглушка, избегаем мёртвого кода и ложного впечатления о PCA | Поиск usages и тест заглушки |
| `README.md` → краткий старт, детали → `docs/USER_GUIDE.md` | Упростить навигацию без потери ранее написанной инструкции | Проверка CLI-флагов, путей, YAML и ссылок |
| Удалён `/tests/` из `.gitignore`, обновлён `AGENTS.md` | Тесты должны попадать в репозиторий, архитектура по всем аутосомам | Проверка .gitignore и поиск устаревших терминов |

## Контракты, которые сохранены

- `pipeline.py`: `--fastq`, `--vcf` (только заглушка), `--vcf --metadata/--labels --species`.
- **Один HaplotypeCaller на образец** для всех указанных аутосом.
- **Один GenomicsDBImport и один GenotypeGVCFs на вид** для всех аутосом,
  без `concat` и без per-chromosome SLURM tasks.
- Один SLURM array на стадию, последовательность `afterok`, глобальный для
  данного запуска лимит `slurm.max_concurrent_jobs`.
- Имена основных выходов `cohort.raw.bcf`, `cohort.raw.vcf.gz`,
  `minimal_panel.tsv`, `cv_metrics.tsv`, `nested_cv_metrics.tsv`.
- Контроль входных хешей и plan schema v3.

## Нерешённые ограничения

- Реальный end-to-end запуск на GATK/SLURM не проводился в этой среде.
- Популяционная структура намеренно остаётся заглушкой.
- Режим SNP-панелей остаётся локальным; научное качество и переносимость
  классификации требуют отдельной внешней валидации.
- Модификации Python-кода делают предыдущие run-id несовместимыми с
  `code_sha256`: требуется новый run-id.
