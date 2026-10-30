# Stage tool upgrade — 29.10.2026

Пакет рассчитан на `develop` после коммита `b598526` (BWA-MEM3).
Он также совместим с последующим рефакторингом `wgs_refactor_next_20261029`:
файлы `pipeline.py`, `run_array.sh` и `workflow.py` не изменяются.

Изменённые файлы:
- `WGS_GATK_simple/envs/wgs.yaml` — обновлены fastp, GATK, mosdepth, MultiQC;
- `WGS_GATK_simple/scripts/stages/trim.sh` — создаются каталоги для всех выходов;
- `WGS_GATK_simple/tests/test_stage_tools_upgrade.py` — 4 offline-теста;
- `WGS_GATK_simple/docs/STAGE_TOOLS_20261029.md` — аудит всех 10 стадий.

## Применение

```bash
cd /mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
# убедитесь, что работаем в develop и текущие изменения сохранены
# git status --short
tar -tzf /PATH/TO/wgs_stages_tools_20261029.tar.gz
tar -xzf /PATH/TO/wgs_stages_tools_20261029.tar.gz
python -m pytest -q WGS_GATK_simple/tests tests
for stage in WGS_GATK_simple/scripts/stages/*.sh; do bash -n "$stage" || exit 1; done

git diff --check
# НЕ обновлять старое окружение в месте и НЕ повторно использовать старый workdir
conda env create --dry-run -f WGS_GATK_simple/envs/wgs.yaml
conda env create -f WGS_GATK_simple/envs/wgs.yaml
conda activate wgs-gatk-mem3
bash WGS_GATK_simple/scripts/check_stack.sh
```

**Важное ограничение:** несовместимость BWA-MEM3/HTSlib 1.25 в текущей сборке.
Производительность и разрешение зависимостей на `sphinx` ещё не проверены.
