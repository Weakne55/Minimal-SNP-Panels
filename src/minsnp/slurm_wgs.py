"""Поэтапные SLURM arrays для смешанной WGS-когорты.

Каждая active stage имеет один array, зависимый от успешного окончания
всего предыдущего array (afterok). Так ограничение %max_concurrent_jobs
действует на все задания данного run-id, включая совместное генотипирование по виду.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import shlex
import subprocess

from . import joint, wgs
from .core import ROOT, config, file_exists, write_tsv
from .resources import (
    JOINT_STAGES, LANE_STAGES, SAMPLE_STAGES, STAGE_ORDER,
    get_parallelism, job_allocation, validate_resources,
)

LOG = logging.getLogger('minsnp')
COLUMNS = wgs.FASTQ_COLUMNS


def sha256(path: str | Path) -> str:
    return hashlib.sha256(file_exists(path).read_bytes()).hexdigest()


def code_hash() -> str:
    files = [ROOT / 'pipeline.py', *sorted((ROOT / 'src/minsnp').glob('*.py')),
             *sorted((ROOT / 'stages').glob('*.sh')),
             ROOT / 'stages/resource_lib.bash']
    digest = hashlib.sha256()
    for path in sorted(set(files)):
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(file_exists(path).read_bytes())
    return digest.hexdigest()


def task_at(items: list, index: int):
    if not 0 <= index < len(items):
        raise ValueError(f'Недопустимый индекс SLURM array: {index}')
    return items[index]


def _write_json_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(path)


def _check_partial_plan_files(run: Path, manifests: list[tuple[Path, list]]) -> None:
    allowed_files = {run / 'slurm_plan.json.tmp'}
    allowed_dirs = {run / 'manifests'}
    for manifest, _ in manifests:
        allowed_files.update((manifest, manifest.with_name(manifest.name + '.tmp')))
        allowed_dirs.add(manifest.parent)
    for entry in run.rglob('*'):
        if (entry.is_symlink() or
            (entry.is_dir() and entry not in allowed_dirs) or
            (entry.is_file() and entry not in allowed_files)):
            raise ValueError(f'Каталог запуска без плана содержит посторонний файл: {entry}')


def make_tasks(samples: list[dict], trim: bool) -> dict[str, list[dict]]:
    """Зафиксировать порядок элементов массива для каждой стадии."""
    tasks = {name: [] for name in STAGE_ORDER}
    for ix, sample in enumerate(samples):
        for run_id, library in sample['lanes']:
            task = {'sample_index': ix, 'run_id': run_id, 'library': library}
            if trim:
                tasks['trim'].append(task)
            tasks['align'].append(task)
        for name in SAMPLE_STAGES:
            tasks[name].append({'sample_index': ix})
    species_list = sorted({s['species'] for s in samples})
    for species in species_list:
        # Все задачи joint genotyping относятся ко всему набору аутосом вида.
        for name in JOINT_STAGES:
            tasks[name].append({'species': species})
    return tasks


def prepare_plan(run: Path, sheets: list[str], refs: str, settings: str,
                 cfg: dict, resume: bool = False) -> Path:
    validate_resources(cfg)
    cohort = wgs.load_cohort(sheets)
    run = run.expanduser().resolve()
    samples, manifests = [], []
    for species in sorted(cohort):
        config(Path(refs), Path(settings), species)
        for sample in sorted(cohort[species]):
            records = cohort[species][sample]
            manifest = run / 'manifests' / species / f'{sample}.tsv'
            values = [tuple(row[key] for key in COLUMNS) for row in records]
            lanes = sorted([row['run_id'], row['library']] for row in records)
            samples.append({'species': species, 'sample': sample,
                            'manifest': str(manifest), 'lanes': lanes, 'records': values})
            manifests.append((manifest, values))
    inputs = [{'path': str(Path(sheet).expanduser().resolve()), 'sha256': sha256(sheet)}
              for sheet in sheets]
    read_files = sorted({row[4] for item in samples for row in item['records']} |
                        {row[5] for item in samples for row in item['records']})
    fastq_stat = []
    for filename in read_files:
        st = file_exists(filename).stat()
        fastq_stat.append([filename, st.st_size, st.st_mtime_ns])
    tasks = make_tasks(samples, bool(cfg.get('trim_fastq', False)))
    plan = {
        'schema_version': 3, 'run': str(run),
        'references': str(Path(refs).resolve()), 'reference_sha256': sha256(refs),
        'settings': str(Path(settings).resolve()), 'settings_sha256': sha256(settings),
        'code_sha256': code_hash(), 'inputs': inputs, 'fastq_stat': fastq_stat,
        'samples': [{k: s[k] for k in ('species', 'sample', 'manifest', 'lanes')} for s in samples],
        'tasks': tasks,
        'stages': [s for s in STAGE_ORDER if tasks[s]],
    }
    path = run / 'slurm_plan.json'
    if path.exists():
        if not resume:
            raise ValueError(f'Каталог запуска уже существует: {run}; нужен --resume или новый --run-id')
        if json.loads(path.read_text()) != plan:
            raise ValueError(f'Изменился план существующего запуска: {path}')
    else:
        if run.exists() and any(run.iterdir()):
            if not resume:
                raise ValueError(f'Каталог запуска не пуст, но нет плана: {run}; нужен --resume')
            _check_partial_plan_files(run, manifests)
        run.mkdir(parents=True, exist_ok=True)
        for manifest, rows in manifests:
            tmp = manifest.with_name(manifest.name + '.tmp')
            write_tsv(tmp, list(COLUMNS), rows)
            tmp.replace(manifest)
        _write_json_atomic(path, plan)
    return path


def verify_plan(plan: dict) -> None:
    if plan.get('schema_version') != 3:
        raise ValueError('Устаревший SLURM-план: используйте новый --run-id')
    for name, key in (('references', 'reference_sha256'), ('settings', 'settings_sha256')):
        if sha256(plan[name]) != plan[key]:
            raise ValueError(f'Изменён конфиг активного анализа: {plan[name]}')
    if code_hash() != plan['code_sha256']:
        raise ValueError('Изменился код с момента создания SLURM-плана')
    for item in plan['inputs']:
        if sha256(item['path']) != item['sha256']:
            raise ValueError(f'Изменён FASTQ TSV: {item["path"]}')
    for filename, size, mtime in plan['fastq_stat']:
        st = file_exists(filename).stat()
        if (st.st_size, st.st_mtime_ns) != (size, mtime):
            raise ValueError(f'Изменён исходный FASTQ: {filename}')


def _request(run: Path, cfg: dict, stage_name: str, python: str,
             dependency: str | None, count: int) -> list[str]:
    """array на конкретный этап, не смешиваем CPU budgets разных этапов."""
    if count < 1:
        raise ValueError(f'Пустой этап {stage_name}')
    cpus, memory = job_allocation(cfg, stage_name)
    slurm = cfg['slurm']
    command = [python, str(ROOT / 'pipeline.py'), '_worker', '--plan', str(run / 'slurm_plan.json'),
               '--stage', stage_name]
    wrapped = shlex.join(command) + ' --index "${SLURM_ARRAY_TASK_ID:?}"'
    req = [
        'sbatch', '--parsable', '--nodes=1', '--ntasks=1', '--export=ALL',
        f'--job-name=minsnp_{stage_name}', f'--cpus-per-task={cpus}', f'--mem={memory}',
        f'--time={slurm["stages"][stage_name]["time"]}', f'--chdir={ROOT}',
        f'--output={run / "logs" / "slurm" / (stage_name + "_%A_%a.out")}',
        f'--error={run / "logs" / "slurm" / (stage_name + "_%A_%a.err")}',
        f'--array=0-{count-1}%{get_parallelism(cfg)}',
    ]
    if slurm.get('partition'):
        req.append(f'--partition={slurm["partition"]}')
    if slurm.get('nodelist'):
        req.append(f'--nodelist={slurm["nodelist"]}')
    if dependency:
        req.append(f'--dependency=afterok:{dependency}')
    req.append('--wrap=' + wrapped)
    return req


def submit_jobs(run: Path, cfg: dict, python: str) -> dict[str, str]:
    plan = json.loads(file_exists(run / 'slurm_plan.json').read_text())
    verify_plan(plan)
    validate_resources(cfg)
    jobs_file = run / 'slurm_jobs.json'
    partial_file = run / 'slurm_jobs.partial.json'
    lock = run / '.submission.lock'
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise ValueError(f'Отправка заданий для {run} уже выполняется или была прервана') from exc
    try:
        if jobs_file.exists() or partial_file.exists():
            raise ValueError(f'Задания уже отправлялись для {run}; повторная отправка заблокирована')
        return _submit_locked(run, cfg, python, plan, jobs_file, partial_file)
    finally:
        lock.rmdir()


def _submit_locked(run: Path, cfg: dict, python: str, plan: dict,
                   jobs_file: Path, partial_file: Path) -> dict[str, str]:
    """Отправить набор arrays после эксклюзивной блокировки run-id."""
    (run / 'logs' / 'slurm').mkdir(parents=True, exist_ok=True)
    jobs: dict[str, str] = {}
    dependency = None
    try:
        for stage_name in plan['stages']:
            count = len(plan['tasks'][stage_name])
            request = _request(run, cfg, stage_name, python, dependency, count)
            result = subprocess.run(request, capture_output=True, text=True, check=True)
            job = result.stdout.strip().split(';')[0]
            if not job.isdecimal():
                raise RuntimeError(f'sbatch не вернул Job ID: {result.stdout!r}')
            jobs[stage_name] = job
            dependency = job
            _write_json_atomic(partial_file, jobs)
            LOG.info('Submitted stage=%s job=%s tasks=%s max_active=%s',
                     stage_name, job, count, get_parallelism(cfg))
    except Exception:
        if not partial_file.exists():
            _write_json_atomic(partial_file, jobs)
        raise
    _write_json_atomic(jobs_file, jobs)
    partial_file.unlink(missing_ok=True)
    return jobs


def setup_logging(level: str, filename: Path | None = None) -> None:
    """Console + UTF-8 file; один worker пишет только в собственный лог."""
    logger = logging.getLogger('minsnp')
    logger.setLevel(getattr(logging, level))
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    formatter = logging.Formatter('%(asctime)s %(levelname)s [%(process)d] %(message)s',
                                  datefmt='%Y-%m-%dT%H:%M:%S')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)
    if filename is not None:
        filename.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename, encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def worker(plan_path: str, stage_name: str, index: int = 0) -> Path:
    plan = json.loads(file_exists(plan_path).read_text())
    verify_plan(plan)
    if stage_name not in plan['stages']:
        raise ValueError(f'Неизвестная или отключённая стадия: {stage_name}')
    task = task_at(plan['tasks'][stage_name], index)
    run = Path(plan['run'])
    # Конфигурация одинаковая для всех видов, поэтому можно читать из YAML.
    import yaml
    cfg = yaml.safe_load(Path(plan['settings']).read_text())
    setup_logging(cfg.get('logging', {}).get('level', 'INFO'),
                  run / 'logs' / 'workers' / stage_name / f'{index}.log')
    LOG.info('START stage=%s task=%s details=%s', stage_name, index, task)
    done = run / 'completed' / stage_name / f'{index}.json'
    try:
        if stage_name in (*LANE_STAGES, *SAMPLE_STAGES):
            item = task_at(plan['samples'], task['sample_index'])
            species = item['species']
            _, spec = config(Path(plan['references']), Path(plan['settings']), species)
            out = run / 'samples' / species / item['sample']
            lane = (task['run_id'], task['library']) if stage_name in LANE_STAGES else None
            result = wgs.run_stage(stage_name, item['manifest'], species, spec, cfg,
                                   out / 'wgs', out / 'logs', lane)
        elif stage_name in JOINT_STAGES:
            species = task['species']
            _, spec = config(Path(plan['references']), Path(plan['settings']), species)
            out = run / 'joint' / species
            result = joint.run_stage(stage_name, plan, species, spec,
                                     run, cfg, out, out / 'logs')
        else:
            raise ValueError(f'Недопустимый этап: {stage_name}')
        _write_json_atomic(done, {'stage': stage_name, 'index': index,
                                  'task': task, 'output': str(result), 'status': 'completed'})
        LOG.info('DONE stage=%s task=%s output=%s', stage_name, index, result)
        return result
    except Exception:
        done.unlink(missing_ok=True)
        LOG.exception('FAILED stage=%s task=%s', stage_name, index)
        raise


def status(run: Path) -> None:
    plan = json.loads(file_exists(run / 'slurm_plan.json').read_text())
    jobs_path = run / 'slurm_jobs.json'
    if not jobs_path.exists():
        jobs_path = run / 'slurm_jobs.partial.json'
    jobs = json.loads(jobs_path.read_text()) if jobs_path.exists() else {}
    print(f'Run: {run}; образцов: {len(plan["samples"])}')
    for stage_name in plan['stages']:
        total = len(plan['tasks'][stage_name])
        done = sum((run / 'completed' / stage_name / f'{i}.json').is_file()
                   for i in range(total))
        job = jobs.get(stage_name)
        summary = f'{stage_name}: {done}/{total} завершено'
        if job:
            try:
                result = subprocess.run(['squeue', '-h', '-r', '-j', job, '-o', '%T'],
                                        capture_output=True, text=True, check=False)
                states = result.stdout.strip().splitlines()
                from collections import Counter
                distribution = dict(Counter(states)) if result.returncode == 0 else {}
                summary += f' | SLURM {job}: {distribution or "нет в squeue"}'
            except FileNotFoundError:
                summary += f' | SLURM {job}: squeue недоступен'
        print(summary)
    for species in sorted({s['species'] for s in plan['samples']}):
        vcf = run / 'joint' / species / 'cohort.raw.vcf.gz'
        print(f'{species}: {vcf} — {"файл и индекс есть" if joint.valid_variant(vcf) else "не подтверждён"}')
