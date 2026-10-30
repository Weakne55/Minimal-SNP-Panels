# Организация репозитория

```text
Minimal-SNP-Panels/
├── run_pipeline.sh              # CLI
├── src/minsnp/                 # Python: CLI, WGS, joint, QC, structure, panels
├── stages/                     # Bash-стадии FASTQ → VCF и QC
├── config/                     # Референсы, ресурсы, параметры анализа
├── envs/                       # Conda
├── inputs/                     # Шаблоны TSV, не реальные манифесты
├── tools/
│   ├── full_wgs_setup/         # Подготовка WGS-выборок и загрузок
│   └── genome_reference_setup/ # Подготовка референсов/индексов
├── data/
│   ├── metadata/               # WGS_final.xlsx + README_metadata.md
│   └── wgs/                    # Реальные FASTQ, не включать в Git
├── results/                    # Реальные результаты, не включать в Git
├── tests/                      # Тесты текущего pipeline
├── docs/                       # Актуальные инструкции и решения
├── literature/                 # Обзор методов
├── benchmarks/                 # Производительность
└── legacy/                     # Старая реализация и исторические файлы
```

Выполняемый код — `run_pipeline.sh`, `src/minsnp/`, `stages/`.
В `tools/` живут утилиты подготовки данных и референсов, которые могут
иметь собственные тесты. `legacy/` не является активным pipeline.

`data/wgs/` содержит десятки гигабайт FASTQ. `results/` содержит
результаты старого эксперимента и новых запусков. Эти каталоги **не
перемещать, не удалять, не коммитить в Git**.

Пути FASTA и список аутосом находятся в `config/references.yaml`.
`inputs/fastq.template.tsv` — шаблон без образцов. Рабочий manifest
создаётся отдельно перед запуском.
