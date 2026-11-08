# FASTQ-манифесты

Рабочие TSV с абсолютными путями не публикуются в Git.

Используйте `fastq.template.tsv`.

Обязательные столбцы:

- sample_id
- species
- run_id
- library
- fastq_r1
- fastq_r2

Одна строка соответствует sequencing unit.

Несколько sequencing units могут принадлежать
одному образцу.

Файлы FASTQ должны быть доступны с compute node.
