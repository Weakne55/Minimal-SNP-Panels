"""Оркестрация одной смешанной WGS-когорты и joint отдельно по видам.

В SLURM отправляется один array по всем животным и одно joint-задание
на каждый вид. Ни одно задание по отдельной хромосоме не создаётся.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

from . import joint, wgs
from .core import ROOT, check_index, config, file_exists, read_tsv, samples_in_vcf, stage, write_tsv
from .resources import get_parallelism, job_allocation, validate_resources

COLUMNS = ('sample_id', 'species', 'run_id', 'library', 'fastq_r1', 'fastq_r2')


def sha256(path: str | Path) -> str:
    """Хешировать небольшие манифесты/конфиги/код, но не гигабайтные FASTQ."""
    return hashlib.sha256(file_exists(path).read_bytes()).hexdigest()


def code_hash() -> str:
    """Фиксировать исполняемый код, чтобы следующие workers не сменили версию."""
    modules = ('analysis_stub.py', 'core.py', 'joint.py', 'resources.py',
               'slurm_wgs.py', 'wgs.py')
    files = [ROOT / 'pipeline.py',
             *(ROOT / 'src/minsnp' / name for name in modules),
             *sorted((ROOT / 'stages').glob('*.sh')),
             ROOT / 'stages/resource_lib.bash']
    digest = hashlib.sha256()
    for path in files:
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(file_exists(path).read_bytes())
    return digest.hexdigest()


def task_at(items: list, index: int):
    if not 0 <= index < len(items):
        raise ValueError(f'Недопустимый индекс SLURM array: {index}')
    return items[index]


def _write_json_atomic(path: Path, value: dict) -> None:
    """Не оставлять частично записанные служебные JSON-файлы."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def _check_partial_plan_files(run: Path, manifests: list[tuple[Path, list]]) -> None:
    """Разрешить resume только после незавершённой записи ожидаемых manifests.

    Без опубликованного slurm_plan.json происхождение чужих файлов неизвестно:
    ничего постороннего не удаляем и не перезаписываем.
    """
    allowed_files = {run / 'slurm_plan.json.tmp'}
    allowed_dirs = {run / 'manifests'}
    for manifest, _ in manifests:
        allowed_files.add(manifest)
        allowed_files.add(manifest.with_name(manifest.name + '.tmp'))
        allowed_dirs.add(manifest.parent)
    for entry in run.rglob('*'):
        if entry.is_symlink() or (entry.is_dir() and entry not in allowed_dirs) or (entry.is_file() and entry not in allowed_files):
            raise ValueError(f'Каталог запуска без плана содержит посторонний файл: {entry}')


def prepare_plan(run: Path, sheets: list[str], refs: str, settings: str,
                 cfg: dict, resume: bool = False) -> Path:
    """Проверить весь состав когорты и сохранить неизменяемый план.

    Применяется как при local, так и при SLURM, чтобы оба backend вызывали
    один и тот же код обработки животного и joint genotyping.
    """
    validate_resources(cfg)
    cohort = wgs.load_cohort(sheets)
    run = run.expanduser().resolve()
    samples = []
    manifests = []
    for species in sorted(cohort):
        # Референс должен быть задан для каждого вида, который есть в TSV.
        config(Path(refs), Path(settings), species)
        for sample in sorted(cohort[species]):
            rows = cohort[species][sample]
            manifest = run / 'manifests' / species / f'{sample}.tsv'
            values = [tuple(row[key] for key in COLUMNS) for row in rows]
            samples.append({'species': species, 'sample': sample,
                            'manifest': str(manifest), 'records': values})
            manifests.append((manifest, values))

    # Контрольные суммы исходных TSV и размер/mtime FASTQ для обнаружения подмены.
    inputs = [{'path': str(Path(sheet).expanduser().resolve()), 'sha256': sha256(sheet)}
              for sheet in sheets]
    read_files = sorted({row[4] for sample in samples for row in sample['records']} |
                        {row[5] for sample in samples for row in sample['records']})
    fastq_stat = []
    for filename in read_files:
        stat = file_exists(filename).stat()
        fastq_stat.append([filename, stat.st_size, stat.st_mtime_ns])
    plan = {'schema_version': 1, 'run': str(run), 'references': str(Path(refs).resolve()),
            'reference_sha256': sha256(refs), 'settings': str(Path(settings).resolve()),
            'settings_sha256': sha256(settings), 'code_sha256': code_hash(),
            'inputs': inputs, 'fastq_stat': fastq_stat,
            'samples': [{'species': s['species'], 'sample': s['sample'],
                         'manifest': s['manifest']} for s in samples]}
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
            # Сначала пишем отдельный временный TSV, затем публикуем атомарно.
            # При обрыве --resume перезапишет только ожидаемые временные файлы.
            temporary = manifest.with_name(manifest.name + '.tmp')
            write_tsv(temporary, list(COLUMNS), rows)
            temporary.replace(manifest)
        _write_json_atomic(path, plan)
    return path


def verify_plan(plan: dict) -> None:
    """Не допускать обработки разными версиями кода или метаданных."""
    for name, key in (('references', 'reference_sha256'), ('settings', 'settings_sha256')):
        if sha256(plan[name]) != plan[key]:
            raise ValueError(f'Изменён файл активного анализа: {plan[name]}')
    if code_hash() != plan['code_sha256']:
        raise ValueError('Изменился исходный код с момента создания SLURM-плана')
    for item in plan['inputs']:
        if sha256(item['path']) != item['sha256']:
            raise ValueError(f'Изменён FASTQ TSV: {item["path"]}')
    for filename, size, mtime in plan['fastq_stat']:
        stat = file_exists(filename).stat()
        if (stat.st_size, stat.st_mtime_ns) != (size, mtime):
            raise ValueError(f'Изменён исходный FASTQ: {filename}')


def _request(run: Path, cfg: dict, group: str, species: str | None,
             python: str, dependency: str | None = None, count: int = 0) -> list[str]:
    """Построить sbatch без shell-скрипта-обёртки."""
    slurm = cfg['slurm']
    resource_group = 'sample' if group == 'sample' else 'joint'
    cpus, memory = job_allocation(cfg, resource_group)
    name = 'wgs' if group == 'sample' else f'joint_{species}'
    command = [python, str(ROOT / 'pipeline.py'), '_worker', '--plan', str(run / 'slurm_plan.json'),
               '--group', group]
    if species:
        command += ['--worker-species', species]
    wrapped = shlex.join(command)
    if group == 'sample':
        wrapped += ' --index "${SLURM_ARRAY_TASK_ID:?}"'
    request = ['sbatch', '--parsable', '--nodes=1', '--ntasks=1', '--export=ALL',
               f'--job-name=minsnp_{name}', f'--cpus-per-task={cpus}', f'--mem={memory}',
               f'--time={slurm[resource_group]["time"]}', f'--chdir={ROOT}',
               f'--output={run / "logs" / (name + "_%A_%a.out")}',
               f'--error={run / "logs" / (name + "_%A_%a.err")}']
    if slurm.get('partition'):
        request.append(f'--partition={slurm["partition"]}')
    if slurm.get('nodelist'):
        request.append(f'--nodelist={slurm["nodelist"]}')
    if group == 'sample':
        request.append(f'--array=0-{count - 1}%{get_parallelism(cfg)}')
    if dependency:
        request.append(f'--dependency=afterok:{dependency}')
    request.append('--wrap=' + wrapped)
    return request


def submit_jobs(run: Path, cfg: dict, python: str) -> dict[str, str]:
    """Один общий WGS array и одно независимое joint-задание на каждый вид."""
    plan = json.loads(file_exists(run / 'slurm_plan.json').read_text())
    validate_resources(cfg)
    if not plan['samples']:
        raise ValueError('Пустая когорта')
    jobs_file = run / 'slurm_jobs.json'
    partial_file = run / 'slurm_jobs.partial.json'
    if jobs_file.exists() or partial_file.exists():
        raise ValueError(f'Задания уже отправлялись для {run}; повторная отправка заблокирована')
    (run / 'logs').mkdir(parents=True, exist_ok=True)
    jobs: dict[str, str] = {}
    species_list = sorted({s['species'] for s in plan['samples']})
    try:
        for key, group, species in [('samples', 'sample', None)] + [
            (f'joint_{species}', 'joint', species) for species in species_list]:
            request = _request(run, cfg, group, species, python,
                               dependency=jobs.get('samples') if group == 'joint' else None,
                               count=len(plan['samples']))
            result = subprocess.run(request, capture_output=True, text=True, check=True)
            job = result.stdout.strip().split(';')[0]
            if not job.isdecimal():
                raise RuntimeError(f'sbatch не вернул Job ID: {result.stdout!r}')
            jobs[key] = job
            # Фиксируем каждый успешно принятый job до следующего sbatch.
            _write_json_atomic(partial_file, jobs)
    except Exception:
        if not partial_file.exists():
            _write_json_atomic(partial_file, jobs)
        raise
    _write_json_atomic(jobs_file, jobs)
    partial_file.unlink(missing_ok=True)
    return jobs


def _collect_joint_manifest(plan: dict, species: str, spec: dict, run: Path) -> Path:
    rows = []
    selected = [s for s in plan['samples'] if s['species'] == species]
    for item in selected:
        manifest = run / 'samples' / species / item['sample'] / 'wgs/gvcf_manifest.tsv'
        entries = read_tsv(manifest, ('sample_id', 'contig', 'gvcf_path'))
        if len(entries) != len(spec['contigs']) or {e['contig'] for e in entries} != set(spec['contigs']):
            raise ValueError(f'Неполный gVCF manifest: {manifest}')
        if any(e['sample_id'] != item['sample'] for e in entries):
            raise ValueError(f'Чужие sample_id в gVCF manifest: {manifest}')
        rows.extend((e['sample_id'], e['contig'], e['gvcf_path']) for e in entries)
    if not rows:
        raise ValueError(f'Нет образцов для joint: {species}')
    target = run / 'manifests' / species / 'cohort_gvcf.tsv'
    if target.is_file():
        prior = read_tsv(target, ('sample_id', 'contig', 'gvcf_path'))
        if [tuple(entry[k] for k in ('sample_id', 'contig', 'gvcf_path')) for entry in prior] != rows:
            raise ValueError(f'Изменился состав joint: {target}')
    else:
        write_tsv(target, ['sample_id', 'contig', 'gvcf_path'], rows)
    return target


def worker(plan_path: str, group: str, index: int = 0, species: str | None = None) -> Path:
    """Выполнить ровно одно животное либо joint одного вида."""
    plan = json.loads(file_exists(plan_path).read_text())
    verify_plan(plan)
    run = Path(plan['run'])
    if group == 'sample':
        item = task_at(plan['samples'], index)
        cfg, spec = config(Path(plan['references']), Path(plan['settings']), item['species'])
        out = run / 'samples' / item['species'] / item['sample']
        manifest = wgs.run(item['manifest'], item['species'], spec, cfg,
                           out / 'wgs', out / 'logs', True)
        print(f'Завершён WGS {item["species"]}/{item["sample"]}: {manifest}', flush=True)
        return manifest
    if group != 'joint' or species not in {s['species'] for s in plan['samples']}:
        raise ValueError(f'Недопустимая группа или вид: {group}/{species}')
    cfg, spec = config(Path(plan['references']), Path(plan['settings']), species)
    sheet = _collect_joint_manifest(plan, species, spec, run)
    out = run / 'joint' / species
    logs = run / 'logs' / species
    bcf = joint.run(sheet, spec, out, logs, True, cfg=cfg)
    vcf = out / 'cohort.raw.vcf.gz'
    if not vcf.is_file():
        stage('export_vcf.sh', [bcf, vcf], logs, 'export_vcf', cfg=cfg)
    check_index(vcf)
    expected = {s['sample'] for s in plan['samples'] if s['species'] == species}
    found = samples_in_vcf(vcf)
    if len(found) != len(expected) or set(found) != expected:
        raise ValueError(f'Неверный состав образцов в {vcf}: {found}')
    print(f'Готовый совместный VCF для {species}: {vcf}', flush=True)
    return vcf


def status(run: Path) -> None:
    """Показать статус без предположения, что исчезнувший из squeue job успешен."""
    plan = json.loads(file_exists(run / 'slurm_plan.json').read_text())
    jobs_path = run / 'slurm_jobs.json'
    jobs = json.loads(jobs_path.read_text()) if jobs_path.exists() else {}
    for group, job in jobs.items():
        try:
            result = subprocess.run(['squeue', '-h', '-r', '-j', job, '-o', '%i|%T|%R'],
                                    capture_output=True, text=True, check=False)
        except FileNotFoundError:
            print(f'{group} ({job}): squeue недоступен')
            continue
        if result.returncode:
            print(f'{group} ({job}): ошибка squeue: {result.stderr.strip()}')
        else:
            print(f'{group} ({job}): {result.stdout.strip() or "нет в очереди; результат не подтверждён"}')
    total = len(plan['samples'])
    done = sum((run / 'samples' / s['species'] / s['sample'] / 'wgs/gvcf_manifest.tsv').is_file()
               for s in plan['samples'])
    print(f'Найдено WGS-манифестов: {done}/{total} (не полная проверка gVCF)')
    for species in sorted({s['species'] for s in plan['samples']}):
        vcf = run / 'joint' / species / 'cohort.raw.vcf.gz'
        print(f'{species}: {vcf} — ' + ('файл и индекс присутствуют' if vcf.is_file() and
              (Path(str(vcf) + '.csi').is_file() or Path(str(vcf) + '.tbi').is_file())
              else 'ещё не подтверждён'))
