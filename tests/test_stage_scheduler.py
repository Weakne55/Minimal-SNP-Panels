"""Контракты нового scheduler, проверяются без SLURM/GATK/FASTQ обработки."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from minsnp import resources, slurm_wgs


@pytest.fixture
def settings():
    root = Path(__file__).resolve().parents[1]
    return yaml.safe_load((root / 'config/pipeline.yaml').read_text())


def test_resources_per_stage_and_limit(settings):
    resources.validate_resources(settings)
    assert resources.job_allocation(settings, 'haplotypecaller') == (8, '32G')
    assert resources.get_parallelism(settings) == 3
    wrong = copy.deepcopy(settings)
    wrong['slurm']['max_concurrent_jobs'] = 0
    with pytest.raises(ValueError, match='max_concurrent_jobs'):
        resources.validate_resources(wrong)
    wrong = copy.deepcopy(settings)
    wrong['slurm']['stages']['align']['memory'] = '3banana'
    with pytest.raises(ValueError, match='памяти'):
        resources.validate_resources(wrong)


def test_build_tasks_has_one_hc_per_sample_and_joint_per_species(tmp_path, monkeypatch):
    samples = [
        {'species': 'beluga', 'sample': 'B1', 'lanes': [['R1', 'L1'], ['R2', 'L1']]},
        {'species': 'polar_bear', 'sample': 'P1', 'lanes': [['R3', 'L2']]},
    ]
    monkeypatch.setattr(slurm_wgs, 'config',
                        lambda refs, settings, species: ({}, {'contigs': ['c1', 'c2'] if species == 'beluga' else ['c3']}))
    tasks = slurm_wgs.make_tasks(samples, trim=True)
    assert len(tasks['trim']) == len(tasks['align']) == 3
    assert len(tasks['haplotypecaller']) == 2
    assert len(tasks['markdup']) == len(tasks['merge_bam']) == 2
    assert len(tasks['genomicsdb']) == len(tasks['genotype']) == 2
    assert len(tasks['export_vcf']) == 2
    assert 'concat' not in tasks
    assert tasks['genomicsdb'][0] == {'species': 'beluga'}
    assert tasks['genotype'] == tasks['genomicsdb']
    assert slurm_wgs.make_tasks(samples, trim=False)['trim'] == []


def test_sbatch_stage_order_afterok_and_global_throttle(settings, monkeypatch, tmp_path):
    run = tmp_path / 'run'
    run.mkdir()
    stages = ['align', 'merge_bam', 'markdup', 'bam_qc', 'haplotypecaller',
              'genomicsdb', 'genotype', 'export_vcf']
    tasks = {stage: [{'i': i} for i in range(2 if stage in ('align', 'genomicsdb') else 1)]
             for stage in stages}
    plan = {'stages': stages, 'tasks': tasks}
    (run / 'slurm_plan.json').write_text(json.dumps(plan))
    monkeypatch.setattr(slurm_wgs, 'verify_plan', lambda plan: None)
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, f'{900 + len(calls)}\n', '')
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', fake_run)
    jobs = slurm_wgs.submit_jobs(run, settings, '/path/to/python')
    assert list(jobs) == stages
    assert len(calls) == 8
    for i, call in enumerate(calls):
        assert any(arg.startswith('--array=') and arg.endswith('%3') for arg in call)
        assert f'--cpus-per-task={settings["slurm"]["stages"][stages[i]]["cpus"]}' in call
        dependencies = [arg for arg in call if arg.startswith('--dependency=')]
        assert dependencies == ([] if i == 0 else [f'--dependency=afterok:{900+i}'])
    assert json.loads((run / 'slurm_jobs.json').read_text()) == jobs
    with pytest.raises(ValueError, match='уже отправлялись'):
        slurm_wgs.submit_jobs(run, settings, '/path/to/python')


def test_sbatch_partial_submission_preserves_ids(settings, monkeypatch, tmp_path):
    run = tmp_path / 'run'
    run.mkdir()
    plan = {'stages': ['align', 'markdup'], 'tasks': {'align': [{}], 'markdup': [{}]}}
    (run / 'slurm_plan.json').write_text(json.dumps(plan))
    monkeypatch.setattr(slurm_wgs, 'verify_plan', lambda plan: None)
    calls = []
    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        if len(calls) == 2:
            raise subprocess.CalledProcessError(1, cmd, stderr='quota')
        return subprocess.CompletedProcess(cmd, 0, '510\n', '')
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', fake_run)
    with pytest.raises(subprocess.CalledProcessError):
        slurm_wgs.submit_jobs(run, settings, sys.executable)
    assert json.loads((run / 'slurm_jobs.partial.json').read_text()) == {'align': '510'}


def test_cli_routes():
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    import pipeline
    parser = pipeline.build_parser()
    assert pipeline.resolve_mode(parser.parse_args(['--fastq', 'x.tsv'])) == 'fastq_to_vcf'
    assert pipeline.resolve_mode(parser.parse_args(['--vcf', 'x.bcf'])) == 'population_structure'
    assert pipeline.resolve_mode(parser.parse_args(['--vcf', 'x.bcf', '--metadata', 'meta.tsv'])) == 'snp_panel'
    assert pipeline.resolve_mode(parser.parse_args(['--vcf', 'x.bcf', '--labels', 'meta.tsv'])) == 'snp_panel'


def test_haplotypecaller_has_no_per_contig_loop():
    root = Path(__file__).resolve().parents[1]
    shell = (root / 'stages/haplotypecaller.sh').read_text()
    assert 'for contig in "$@"; do intervals+=(-L "$contig"); done' in shell
    assert shell.count('    HaplotypeCaller ') == 1


def test_prepare_plan_and_resume_with_multiple_lanes(tmp_path):
    refs = tmp_path / 'ref.fa'
    refs.write_text('>c1\nACTG\n>c2\nACTG\n')
    (tmp_path / 'ref.fa.fai').write_text('c1\t4\t4\t4\t5\nc2\t4\t13\t4\t5\n')
    refs_yaml = tmp_path / 'refs.yaml'
    refs_yaml.write_text(yaml.safe_dump({'species': {'beluga': {
        'assembly': 'test', 'fasta': str(refs), 'contigs': ['c1', 'c2']}}}))
    config_yaml = tmp_path / 'pipeline.yaml'
    base = Path(__file__).resolve().parents[1] / 'config/pipeline.yaml'
    config_yaml.write_text(base.read_text())
    inputs = tmp_path / 'fastq.tsv'
    fastqs = []
    for name in ('R1_1', 'R1_2', 'R2_1', 'R2_2'):
        fastq = tmp_path / f'{name}.fq.gz'
        fastq.write_bytes(b'not real fastq: test manifests only\n')
        fastqs.append(fastq)
    inputs.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n' +
                      f'B1\tbeluga\tR1\tL1\t{fastqs[0]}\t{fastqs[1]}\n' +
                      f'B1\tbeluga\tR2\tL1\t{fastqs[2]}\t{fastqs[3]}\n')
    run = tmp_path / 'results' / 'run1'
    cfg = yaml.safe_load(config_yaml.read_text())
    plan_file = slurm_wgs.prepare_plan(run, [str(inputs)], str(refs_yaml),
                                        str(config_yaml), cfg)
    plan = json.loads(plan_file.read_text())
    assert plan['schema_version'] == 3
    assert plan['samples'][0]['sample'] == 'B1'
    assert len(plan['tasks']['align']) == 2
    assert len(plan['tasks']['haplotypecaller']) == 1
    assert len(plan['tasks']['genomicsdb']) == 1
    assert plan['tasks']['genomicsdb'] == [{'species': 'beluga'}]
    assert 'concat' not in plan['stages']
    assert plan['stages'][0] == 'align'
    slurm_wgs.verify_plan(plan)
    assert slurm_wgs.prepare_plan(run, [str(inputs)], str(refs_yaml), str(config_yaml), cfg,
                                  resume=True) == plan_file
    with pytest.raises(ValueError, match='уже существует'):
        slurm_wgs.prepare_plan(run, [str(inputs)], str(refs_yaml), str(config_yaml), cfg)
    inputs.write_text(inputs.read_text() + '\\n')
    with pytest.raises(ValueError, match='FASTQ TSV'):
        slurm_wgs.verify_plan(plan)


def test_hc_worker_runs_once_and_writes_manifest(tmp_path, monkeypatch, settings):
    from minsnp import wgs
    sample = 'B1'
    table = tmp_path / 'manifest.tsv'
    r1, r2 = tmp_path / 'r1.gz', tmp_path / 'r2.gz'
    r1.write_text('dummy')
    r2.write_text('dummy')
    table.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n'
                     f'{sample}\tbeluga\tR1\tLIB1\t{r1}\t{r2}\n')
    out = tmp_path / 'wgs'
    bam = out / 'bam' / f'{sample}.bam'
    bam.parent.mkdir(parents=True)
    bam.write_text('fake-bam')
    Path(str(bam) + '.bai').write_text('fake-index')
    monkeypatch.setattr(wgs, '_valid_bam', lambda path, index=False: path.is_file())
    monkeypatch.setattr(wgs, '_valid_gvcf', lambda path: path.is_file() and Path(str(path) + '.tbi').is_file())
    commands = []
    def fake_stage(script, args, logdir, name, cfg):
        commands.append((script, args))
        Path(args[2]).write_text('dummy-gvcf')
        Path(str(args[2]) + '.tbi').write_text('dummy-index')
    monkeypatch.setattr(wgs, 'stage', fake_stage)
    result = wgs.run_stage('haplotypecaller', str(table), 'beluga',
                            {'fasta': '/ref.fa', 'contigs': ['c1', 'c2']},
                            settings, out, tmp_path / 'logs')
    assert len(commands) == 1
    assert commands[0][0] == 'haplotypecaller.sh'
    assert commands[0][1][-2:] == ['c1', 'c2']
    assert result == out / 'gvcf' / 'B1.g.vcf.gz'
    assert (out / 'gvcf_manifest.tsv').read_text().count('\n') == 2
    wgs.run_stage('haplotypecaller', str(table), 'beluga',
                  {'fasta': '/ref.fa', 'contigs': ['c1', 'c2']},
                  settings, out, tmp_path / 'logs')
    assert len(commands) == 1  # resume: повторный запуск не нужен


def test_structure_is_explicit_stub(tmp_path):
    from minsnp import structure
    vcf = tmp_path / 'cohort.bcf'
    vcf.write_bytes(b'dummy')
    (tmp_path / 'cohort.bcf.csi').write_bytes(b'dummy-index')
    with pytest.raises(NotImplementedError, match='заглушка'):
        structure.run(vcf)
    assert set(tmp_path.iterdir()) == {vcf, tmp_path / 'cohort.bcf.csi'}


def test_worker_logs_and_completion_marker(monkeypatch, settings, tmp_path):
    from minsnp import wgs
    run = tmp_path / 'run'
    run.mkdir()
    plan = {
        'run': str(run), 'settings': str(tmp_path / 'config.yaml'),
        'references': str(tmp_path / 'references.yaml'),
        'stages': ['haplotypecaller'],
        'samples': [{'species': 'beluga', 'sample': 'B1', 'manifest': 'any.tsv'}],
        'tasks': {'haplotypecaller': [{'sample_index': 0}]},
    }
    plan_file = run / 'slurm_plan.json'
    plan_file.write_text(json.dumps(plan))
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    monkeypatch.setattr(slurm_wgs, 'verify_plan', lambda x: None)
    monkeypatch.setattr(slurm_wgs, 'config', lambda *args: ({}, {'fasta': '/ref.fa', 'contigs': ['c1']}))
    output = run / 'some.g.vcf.gz'
    output.write_text('test')
    monkeypatch.setattr(wgs, 'run_stage', lambda *args: output)
    assert slurm_wgs.worker(str(plan_file), 'haplotypecaller', 0) == output
    marker = run / 'completed' / 'haplotypecaller' / '0.json'
    log = run / 'logs' / 'workers' / 'haplotypecaller' / '0.log'
    assert json.loads(marker.read_text())['status'] == 'completed'
    assert 'START stage=haplotypecaller' in log.read_text()
    assert 'DONE stage=haplotypecaller' in log.read_text()
    def fail(*args):
        raise RuntimeError('simulated failure')
    monkeypatch.setattr(wgs, 'run_stage', fail)
    with pytest.raises(RuntimeError, match='simulated'):
        slurm_wgs.worker(str(plan_file), 'haplotypecaller', 0)
    assert not marker.exists()
    assert 'FAILED stage=haplotypecaller' in log.read_text()


def test_joint_runs_one_import_and_genotype_for_all_autosomes(monkeypatch, settings, tmp_path):
    from minsnp import joint
    spec = {'fasta': '/ref.fa', 'contigs': ['c1', 'c2', 'c3']}
    plan = {'samples': [
        {'species': 'beluga', 'sample': 'B1'},
        {'species': 'beluga', 'sample': 'B2'},
        {'species': 'polar_bear', 'sample': 'P1'},
    ]}
    out = tmp_path / 'joint' / 'beluga'
    logs = out / 'logs'
    seen = []

    monkeypatch.setattr(joint, 'sample_map_from_plan', lambda *args: {
        'B1': '/gvcf/B1.g.vcf.gz', 'B2': '/gvcf/B2.g.vcf.gz'})
    monkeypatch.setattr(joint, 'valid_variant', lambda path: path.is_file() and Path(str(path)+'.csi').is_file())
    monkeypatch.setattr(joint, 'samples_in_vcf', lambda path: ['B1', 'B2'])

    def simulate_stage(script, args, logdir, name, cfg):
        seen.append((script, list(map(str, args)), name))
        if script == 'genomicsdb.sh':
            db = Path(args[2]); db.mkdir(parents=True); (db / '.completed').write_text('complete')
        else:
            dest = Path(args[2] if script == 'genotype.sh' else args[1])
            dest.write_text('dummy'); Path(str(dest) + '.csi').write_text('dummy index')

    monkeypatch.setattr(joint, 'stage', simulate_stage)
    db = joint.run_stage('genomicsdb', plan, 'beluga', spec, tmp_path, settings, out, logs)
    bcf = joint.run_stage('genotype', plan, 'beluga', spec, tmp_path, settings, out, logs)
    vcf = joint.run_stage('export_vcf', plan, 'beluga', spec, tmp_path, settings, out, logs)
    assert db == out / 'genomicsdb'
    assert bcf == out / 'cohort.raw.bcf'
    assert vcf == out / 'cohort.raw.vcf.gz'
    assert [x[0] for x in seen] == ['genomicsdb.sh', 'genotype.sh', 'export_vcf.sh']
    assert seen[0][1] == ['/ref.fa', str(out/'sample_map.tsv'), str(db), 'c1', 'c2', 'c3']
    assert seen[1][1] == ['/ref.fa', str(db), str(bcf), 'c1', 'c2', 'c3']
    assert (out / 'sample_map.tsv').read_text() == (
        'B1\t/gvcf/B1.g.vcf.gz\nB2\t/gvcf/B2.g.vcf.gz\n')
    assert not (out / 'per_contig').exists()
    assert not (out / 'genomicsdb' / 'c1').exists()
    for stage in ('genomicsdb', 'genotype', 'export_vcf'):
        result = joint.run_stage(stage, plan, 'beluga', spec, tmp_path, settings, out, logs)
        assert result.exists()
    assert len(seen) == 3  # корректные результаты пропускаются при повторе


def test_joint_shell_scripts_use_every_autosome_in_single_gatk_call():
    root = Path(__file__).resolve().parents[1]
    import os
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        bin_dir = tmp / 'bin'; bin_dir.mkdir()
        calls = tmp / 'calls.txt'
        gatk = bin_dir / 'gatk'
        gatk.write_text('''#!/usr/bin/env python3
import os,sys,pathlib
args=sys.argv[1:]
with open(os.environ['MOCK_CALLS'],'a') as f: f.write('gatk ' + ' '.join(args) + '\\n')
if 'GenomicsDBImport' in args:
    pathlib.Path(args[args.index('--genomicsdb-workspace-path')+1]).mkdir()
else:
    pathlib.Path(args[args.index('-O')+1]).write_text('VCF')
''')
        gatk.chmod(0o755)
        bcftools = bin_dir / 'bcftools'
        bcftools.write_text('''#!/usr/bin/env python3
import os,sys,pathlib
args=sys.argv[1:]
with open(os.environ['MOCK_CALLS'],'a') as f: f.write('bcftools ' + ' '.join(args) + '\\n')
if args[0]=='view': pathlib.Path(args[args.index('-o')+1]).write_text('BCF')
if args[0]=='index' and '-c' in args: pathlib.Path(args[-1]+'.csi').write_text('IDX')
''')
        bcftools.chmod(0o755)
        sample_map=tmp/'samples.tsv'; sample_map.write_text('B1\t/B1.g.vcf.gz\n')
        env=dict(os.environ, PATH=f'{bin_dir}:{os.environ["PATH"]}', MOCK_CALLS=str(calls),
                 MIN_SNP_THREADS='2', MIN_SNP_MEMORY_MB='16384')
        genomicsdb = root/'stages/genomicsdb.sh'
        genotype = root/'stages/genotype.sh'
        subprocess.run(['bash',str(genomicsdb),str(tmp/'ref.fa'),str(sample_map),str(tmp/'db'),
                       'c1','c2','c3'],check=True,env=env)
        subprocess.run(['bash',str(genotype),str(tmp/'ref.fa'),str(tmp/'db'),str(tmp/'cohort.raw.bcf'),
                       'c1','c2','c3'],check=True,env=env)
        lines=calls.read_text().splitlines()
        gatk_calls=[line for line in lines if line.startswith('gatk ')]
        assert len(gatk_calls)==2
        for line in gatk_calls:
            assert '-L c1 -L c2 -L c3' in line
        assert '--reader-threads 1' in gatk_calls[0]
        assert (tmp/'db/.completed').is_file()
        assert (tmp/'cohort.raw.bcf').is_file()
        assert (tmp/'cohort.raw.bcf.csi').is_file()
