#!/usr/bin/env python3
"""SLURM orchestration for the 20-animal WGS pilot.

Usage:
  python wgs_test/wgs_pilot.py plan --selection inputs/wgs_pilot_20.tsv
  python wgs_test/wgs_pilot.py submit --selection inputs/wgs_pilot_20.tsv --run-id pilot_20261030 --auto-resources
  python wgs_test/wgs_pilot.py status --run-id pilot_20261030

Download -> sample manifest -> per-sample WGS -> per-species joint calling.
No synthetic read/genotype data are ever created.  All fastq_ftp paths and MD5s
are resolved from ENA's live filereport at execution time.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "wgs"
RESULTS = REPO / "results"
VALID_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
ACCESSION = re.compile(r"^SRR[0-9]+$")
SPECIES_DIR = {"beluga": "Delphinapterus_leucas", "polar_bear": "Ursus_maritimus"}
REQUIRED = ("species", "sample_id", "bioproject", "coverage_x", "fastq_size_gib", "sra_runs")


def validate_id(value: str) -> str:
    if not VALID_ID.fullmatch(value) or value in (".", ".."):
        raise ValueError(f"Invalid identifier: {value!r}")
    return value


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(REQUIRED).issubset(reader.fieldnames):
            raise ValueError(f"Missing columns in {path}: {REQUIRED}")
        records = list(reader)
    if len(records) != 20:
        raise ValueError(f"The pilot must contain 20 samples, found {len(records)}")
    if {k: sum(x["species"] == k for x in records) for k in SPECIES_DIR} != {"beluga": 10, "polar_bear": 10}:
        raise ValueError("Exactly 10 belugas and 10 polar bears are required")
    samples, run_to_sample = set(), {}
    for row in records:
        sp, sample = row["species"], validate_id(row["sample_id"])
        if sp not in SPECIES_DIR or (sp, sample) in samples:
            raise ValueError(f"Unknown species / duplicated sample: {sp} {sample}")
        samples.add((sp, sample))
        if float(row["fastq_size_gib"]) <= 0 or float(row["coverage_x"]) <= 0:
            raise ValueError(f"Missing/invalid coverage or size: {sample}")
        for run in runs_for(row):
            if not ACCESSION.fullmatch(run) or run in run_to_sample:
                raise ValueError(f"Invalid/repeated SRA run: {run}")
            run_to_sample[run] = (sp, sample)
    return records


def runs_for(row: dict[str, str]) -> list[str]:
    return [x.strip() for x in row["sra_runs"].split(";") if x.strip()]


def tasks(records: list[dict[str, str]]) -> list[tuple[str, str, str]]:
    return [(row["species"], row["sample_id"], run) for row in records for run in runs_for(row)]


def workdir(run_id: str) -> Path:
    return RESULTS / "pilot_orchestration" / validate_id(run_id)


def atomic_text(path: Path, content: str) -> None:
    """Write only after full content is available; never silently rewrite an input manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        if path.read_text() != content:
            raise ValueError(f"Existing file differs, refusing to overwrite: {path}")
        return
    temp = path.with_name(path.name + f".{os.getpid()}.tmp")
    temp.write_text(content)
    os.replace(temp, path)


def eno_run(run: str) -> dict:
    fields = "run_accession,experiment_accession,library_layout,fastq_ftp,fastq_md5,fastq_bytes"
    params = urlencode({"accession": run, "result": "read_run", "fields": fields, "format": "tsv"})
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + params
    request = Request(url, headers={"User-Agent": "MinimalSNP-Panels-pilot/1.0"})
    # Only run metadata is downloaded here; FASTQ content goes through curl.
    with urlopen(request, timeout=90) as response:
        text = response.read().decode("utf-8-sig")
    records = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(records) != 1 or records[0].get("run_accession") != run:
        raise ValueError(f"ENA returned no unique record for {run}: {text[:180]}")
    r = records[0]
    if r.get("library_layout", "").upper() != "PAIRED":
        raise ValueError(f"ENA {run}: library_layout is not PAIRED")
    paths, md5s, sizes = [r.get(k, "").split(";") for k in ("fastq_ftp", "fastq_md5", "fastq_bytes")]
    if not paths or len(paths) != len(md5s) or len(paths) != len(sizes):
        raise ValueError(f"ENA {run}: incomplete FASTQ URLs/MD5/size")
    files = {}
    for ftp_path, md5, size in zip(paths, md5s, sizes):
        name = ftp_path.rsplit("/", 1)[-1]
        match = re.fullmatch(r"(.+)_([12])\.fastq\.gz", name)
        if not match:
            raise ValueError(f"ENA {run}: unpaired or unexpected FASTQ filename: {name}")
        key, mate = match.group(1), int(match.group(2))
        if (key, mate) in files or not re.fullmatch(r"[a-fA-F0-9]{32}", md5):
            raise ValueError(f"ENA {run}: duplicated mate or invalid MD5")
        if not size.isdecimal() or int(size) <= 0:
            raise ValueError(f"ENA {run}: invalid FASTQ byte count")
        host_path = ftp_path.replace("ftp://", "").replace("https://", "", 1)
        parsed = urlparse("https://" + host_path)
        if parsed.hostname not in ("ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk"):
            raise ValueError(f"Unexpected ENA FTP host: {parsed.hostname}")
        files[(key, mate)] = {"url": parsed.geturl(), "md5": md5.lower(), "bytes": int(size), "name": name}
    keys = sorted({key for key, mate in files})
    if any((key, 1) not in files or (key, 2) not in files for key in keys):
        raise ValueError(f"ENA {run}: unmatched R1/R2")
    return {"run_accession": run, "library": validate_id(r.get("experiment_accession") or run),
            "pairs": [{"r1": files[(key, 1)], "r2": files[(key, 2)]} for key in keys]}


def checksum(path: Path) -> str:
    md5 = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            md5.update(chunk)
    return md5.hexdigest()


def download_file(remote: dict, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    size, digest = remote["bytes"], remote["md5"]
    if target.is_file():
        if target.stat().st_size != size or checksum(target) != digest:
            raise ValueError(f"Existing FASTQ disagrees with ENA MD5/size: {target}; do not overwrite")
        print("Existing FASTQ verified:", target, flush=True)
        return
    partial = target.with_name(target.name + ".part")
    if partial.exists() and partial.stat().st_size > size:
        raise ValueError(f"Oversized partial download: {partial}")
    # Каждый повторный запуск curl продолжает существующий .part.
    # Не используем встроенный --retry: он может удалять скачанные байты.
    print("Downloading", remote["url"], "->", partial, flush=True)

    for attempt in range(1, 9):
        if partial.exists() and partial.stat().st_size == size:
            break

        cmd = [
            "curl", "--fail", "--location",
            "--connect-timeout", "60",
            "--speed-limit", "1024", "--speed-time", "300",
            "--continue-at", "-",
            "--output", str(partial),
            remote["url"],
        ]

        result = subprocess.run(cmd, check=False)

        if result.returncode == 0:
            break

        if result.returncode == 33:
            raise RuntimeError(
                f"Сервер не поддерживает возобновление загрузки: {remote['url']}"
            )

        if attempt == 8:
            raise RuntimeError(
                f"Не удалось скачать FASTQ за 8 попыток: {target}"
            )

        print(
            f"curl exit {result.returncode}; повтор {attempt}/8, "
            f"сохранено {partial.stat().st_size if partial.exists() else 0} байт",
            flush=True,
        )
        time.sleep(min(60, 10 * attempt))
    if partial.stat().st_size != size or checksum(partial) != digest:
        raise ValueError(f"Bad FASTQ size/MD5 for {partial}; retaining .part for inspection")
    os.replace(partial, target)
    print("Verified", target, flush=True)


def download_task(selected: list[dict], run_id: str, index: int) -> None:
    items = tasks(selected)
    sp, sample, run = items[index]
    dest = DATA / SPECIES_DIR[sp] / run
    state = workdir(run_id) / "ena" / f"{run}.json"
    if state.exists():
        meta = json.loads(state.read_text())
    else:
        meta = eno_run(run)
        atomic_text(state, json.dumps(meta, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    pairs = []
    for part, pair in enumerate(meta["pairs"], 1):
        # For ordinary single-pair runs keep the established SRR..._1/2.fastq.gz layout.
        prefix = run if len(meta["pairs"]) == 1 else f"{run}_part{part:02d}"
        p1 = dest / f"{prefix}_1.fastq.gz"
        p2 = dest / f"{prefix}_2.fastq.gz"
        download_file(pair["r1"], p1)
        download_file(pair["r2"], p2)
        pairs.append({"sample_id": sample, "species": sp,
                      "run_id": validate_id(prefix), "library": meta["library"],
                      "fastq_r1": str(p1.resolve()), "fastq_r2": str(p2.resolve())})
    verified = workdir(run_id) / "verified" / f"{run}.json"
    atomic_text(verified, json.dumps(pairs, ensure_ascii=False, sort_keys=True, indent=2) + "\n")


def tsv(header: list[str], contents: list[list[str]]) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(contents)
    return buf.getvalue()


def make_manifests(selected: list[dict], run_id: str) -> None:
    dest = workdir(run_id) / "manifests"
    header = ["sample_id", "species", "run_id", "library", "fastq_r1", "fastq_r2"]
    for row in selected:
        all_pairs = []
        for run in runs_for(row):
            validated = workdir(run_id) / "verified" / f"{run}.json"
            if not validated.is_file():
                raise FileNotFoundError(f"Run not fully verified: {validated}")
            pairs = json.loads(validated.read_text())
            if any(r["species"] != row["species"] or r["sample_id"] != row["sample_id"] for r in pairs):
                raise ValueError(f"Wrong species or sample in {validated}")
            for r in pairs:
                for mate in ("fastq_r1", "fastq_r2"):
                    if not Path(r[mate]).is_file():
                        raise FileNotFoundError(r[mate])
            all_pairs.extend(pairs)
        content = tsv(header, [[r[k] for k in header] for r in all_pairs])
        atomic_text(dest / row["species"] / f"{row['sample_id']}.tsv", content)
        print("Sample manifest:", row["species"], row["sample_id"], len(all_pairs), "read pairs")
    atomic_text(dest / "complete.ok", "all manifests created\n")


def species_samples(selected: list[dict], species: str) -> list[dict]:
    return [r for r in selected if r["species"] == species]


def run_cmd(args: list[str]) -> None:
    print("$", shlex.join(args), flush=True)
    subprocess.run(args, cwd=REPO, check=True)


def index_reference(species: str) -> None:
    import yaml
    paths = yaml.safe_load((REPO / "config" / "references.yaml").read_text())
    fasta = Path(paths["species"][species]["fasta"])
    if not fasta.is_file():
        raise FileNotFoundError(fasta)
    required = [Path(str(fasta) + suffix) for suffix in (".bwt.2bit.64", ".pac", ".ann", ".amb")]
    if all(p.is_file() and p.stat().st_size for p in required):
        print("BWA-MEM3 index already exists:", fasta)
        return
    # Fail on an existing BWT before touching reference; never overwrite live indexes.
    if required[0].exists():
        raise ValueError(f"Incomplete existing BWA-MEM3 BWT: {required[0]}")
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK") or "8")
    memory_mb = int(os.environ.get("SLURM_MEM_PER_NODE") or "65536")
    max_memory_gb = max(4, (memory_mb - 4096) // 1024)
    # Build the index in its final location; do not run concurrently for the same species.
    run_cmd(["bwa-mem3", "index", "-t", str(threads), "--max-memory", f"{max_memory_gb}G", str(fasta)])
    if not all(p.is_file() and p.stat().st_size for p in required):
        raise RuntimeError(f"Index incomplete: {fasta}")


def sample_task(selected: list[dict], run_id: str, species: str, index: int) -> None:
    sample = species_samples(selected, species)[index]["sample_id"]
    manifest = workdir(run_id) / "manifests" / species / f"{sample}.tsv"
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "1")
    mem_mb = int(os.environ.get("SLURM_MEM_PER_NODE", "65536"))
    memory = f"{mem_mb}M"
    run_cmd([str(REPO / "run_pipeline.sh"), "wgs", "--fastq", str(manifest),
             "--species", species, "--run-id", f"{run_id}_{sample}",
             "--threads", cpus, "--memory", memory, "--backend", "local", "--resume"])


def joint_task(selected: list[dict], run_id: str, species: str) -> None:
    # One row per sample x autosome; all samples are included in a single joint cohort.
    source_header = ["sample_id", "contig", "gvcf_path"]
    entries = []
    for row in species_samples(selected, species):
        sample = row["sample_id"]
        path = RESULTS / species / f"{run_id}_{sample}" / "wgs" / "gvcf_manifest.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames != source_header:
                raise ValueError(f"Unexpected gVCF manifest: {path}")
            for r in reader:
                if r["sample_id"] != sample or not Path(r["gvcf_path"]).is_file():
                    raise ValueError(f"Wrong sample or missing gVCF: {path}")
                entries.append([r[x] for x in source_header])
    if not entries or len({e[0] for e in entries}) != 10:
        raise ValueError(f"Incomplete cohort for {species}")
    manifest = workdir(run_id) / "manifests" / f"{species}_combined_gvcf.tsv"
    atomic_text(manifest, tsv(source_header, entries))
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "1")
    mem_mb = int(os.environ.get("SLURM_MEM_PER_NODE", "65536"))
    joint_run = f"{run_id}_joint"
    run_cmd([str(REPO / "run_pipeline.sh"), "joint", "--gvcf-manifest", str(manifest),
             "--species", species, "--run-id", joint_run, "--threads", cpus,
             "--memory", f"{mem_mb}M", "--backend", "local", "--resume"])
    output = RESULTS / species / joint_run / "joint"
    bcf, vcf = output / "cohort.raw.bcf", output / "cohort.raw.vcf.gz"
    if not bcf.is_file() or not Path(str(bcf) + ".csi").is_file():
        raise RuntimeError(f"Missing joint BCF and index: {bcf}")
    if not vcf.exists():
        # export_vcf.sh accepts exactly two positional arguments; resource settings
        # are passed via the same MIN_SNP_* environment as in the CLI stages.
        os.environ["MIN_SNP_THREADS"] = cpus
        os.environ["MIN_SNP_MEMORY_MB"] = str(mem_mb)
        run_cmd(["bash", str(REPO / "stages" / "export_vcf.sh"), str(bcf), str(vcf)])
    if not Path(str(vcf) + ".csi").is_file():
        raise RuntimeError(f"VCF missing index: {vcf}")
    for file in (bcf, vcf):
        actual = subprocess.run(["bcftools", "query", "-l", str(file)],
                                check=True, capture_output=True, text=True).stdout.splitlines()
        expected = {r["sample_id"] for r in species_samples(selected, species)}
        if len(actual) != 10 or set(actual) != expected:
            raise ValueError(f"Wrong joint sample set in {file}: {actual}")
        result = subprocess.run(["bcftools", "index", "-n", str(file)],
                                check=True, capture_output=True, text=True)
        print("Verified:", file, "samples=10", "indexed_variants=", result.stdout.strip())
    atomic_text(workdir(run_id) / f"{species}.complete.ok", str(vcf) + "\n")


def capacity_snapshot() -> None:
    print("\n=== Cluster state (informational; SLURM is the allocation authority) ===")
    for title, cmd in (("Nodes", ["sinfo", "-N", "-h", "-o", "%N|%C|%m|%T|%P"]),
                       ("Queue", ["squeue", "-h", "-o", "%i|%u|%C|%m|%T|%P"]),
                       ("Disk", ["df", "-h", str(DATA)])):
        if shutil.which(cmd[0]):
            proc = subprocess.run(cmd, text=True, capture_output=True)
            print(title + ":\n" + (proc.stdout[:3500] if proc.returncode == 0 else proc.stderr[:300]))


def auto_resources(max_cpus: int, max_mem_gib: int) -> tuple[int, int]:
    """Suggest a conservative single-job request from current idle CPU and AllocMem.

    This is only a snapshot, not a reservation. sbatch is still responsible for
    enforcing actual availability and limits on the worker node.
    """
    command = ["sinfo", "-N", "-h", "-o", "%N|%C|%m|%T|%P"]
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    candidates = []
    for line in output.splitlines():
        parts = line.split("|")
        if len(parts) < 4 or not any(x in parts[3].lower() for x in ("idle", "mix")):
            continue
        try:
            node = parts[0]
            idle = int(parts[1].split("/")[1])
            details = subprocess.run(["scontrol", "show", "node", node], check=True,
                                     capture_output=True, text=True).stdout
            real = int(re.search(r"\bRealMemory=(\d+)", details).group(1))
            alloc = int(re.search(r"\bAllocMem=(\d+)", details).group(1))
            free = max(0, real - alloc)
        except (ValueError, IndexError, AttributeError, subprocess.CalledProcessError):
            continue
        cpus = min(max_cpus, max(1, idle // 2))
        memory_gib = min(max_mem_gib, max(0, (free // 2) // 1024))
        if cpus >= 2 and memory_gib >= 32:
            candidates.append((memory_gib * cpus, cpus, memory_gib, node))
    if not candidates:
        raise RuntimeError("No node with >=2 idle CPUs and >=64GiB unallocated memory. Specify --cpus/--memory explicitly, or submit later")
    _, cpus, memory, node = max(candidates)
    print(f"Auto suggestion from {node}: --cpus {cpus} --memory {memory}G (snapshot only)")
    return cpus, memory


def job_script(run_id: str, stage: str, selection: Path, species: str = "", array: bool = False) -> Path:
    path = workdir(run_id) / "jobs" / f"{stage}{'_' + species if species else ''}.sh"
    argv = [str(Path(sys.executable).resolve()), str(Path(__file__).resolve()), "worker",
            "--selection", str(selection), "--run-id", run_id, "--stage", stage]
    if species:
        argv.extend(["--species", species])
    command = shlex.join(argv)
    if array:
        command += ' --index "${SLURM_ARRAY_TASK_ID:?}"'
    body = "#!/usr/bin/env bash\nset -euo pipefail\n" + f"cd {shlex.quote(str(REPO))}\n" + \
           'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"\n' + command + "\n"
    atomic_text(path, body)
    path.chmod(0o750)
    return path


def sbatch(run_id: str, jobname: str, script: Path, cpus: int, mem: str, time_limit: str,
           depend: str = "", array: str = "", partition: str = "") -> str:
    logs = workdir(run_id) / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    args = ["sbatch", "--parsable", "--nodes=1", "--ntasks=1", "--export=ALL",
            "--chdir", str(REPO), "--job-name", jobname,
            "--cpus-per-task", str(cpus), "--mem", mem, "--time", time_limit,
            "--output", str(logs / (jobname + "_%A_%a.out")),
            "--error", str(logs / (jobname + "_%A_%a.err"))]
    if depend:
        args.append("--dependency=afterok:" + depend)
    if array:
        args.append("--array=" + array)
    if partition:
        args.extend(["--partition", partition])
    args.append(str(script))
    print("SUBMIT:", shlex.join(args), flush=True)
    result = subprocess.run(args, check=True, capture_output=True, text=True)
    return result.stdout.strip().split(";")[0]


def disk_report(selected: list[dict], required_factor: float = 2.5) -> bool:
    total = sum(float(x["fastq_size_gib"]) for x in selected)
    existing = sum(float(x["fastq_size_gib"]) for x in selected if x.get("existing_fastq_runs"))
    missing = max(0, total - existing)
    DATA.mkdir(parents=True, exist_ok=True)
    free_gib = shutil.disk_usage(DATA).free / (1024 ** 3)
    minimum = missing * required_factor
    print(f"Estimated FASTQ total={total:.1f} GiB, nominal existing={existing:.1f} GiB")
    print(f"Estimated NEW FASTQ={missing:.1f} GiB; free on data filesystem={free_gib:.1f} GiB")
    print(f"Planning threshold={minimum:.1f} GiB (new FASTQ x {required_factor}); approximate, not a guarantee")
    return free_gib >= minimum


def submit(selected: list[dict], args: argparse.Namespace) -> None:
    if not shutil.which("sbatch") or not shutil.which("curl"):
        raise RuntimeError("Both sbatch and curl must be in PATH")
    wd = workdir(args.run_id)
    wd.mkdir(parents=True, exist_ok=True)
    selected_file = args.selection.expanduser().resolve()
    snapshot = {"selection": str(selected_file), "sha256": hashlib.sha256(selected_file.read_bytes()).hexdigest()}
    atomic_text(wd / "selection.lock.json", json.dumps(snapshot, indent=2) + "\n")
    if (wd / "jobs.json").exists() or (wd / "jobs.partial.json").exists():
        raise ValueError(f"Jobs were already submitted (possibly partially) for {args.run_id}; inspect {wd}")
    if not disk_report(selected) and not args.allow_low_disk:
        raise RuntimeError("Not enough disk space by the pilot safety threshold; use --allow-low-disk only after manual review")
    if args.auto_resources:
        cpus, mem_gib = auto_resources(args.max_cpus, args.max_memory_gib)
        memory = f"{mem_gib}G"
    else:
        cpus = args.cpus
        memory = args.memory
    if cpus < 1 or not re.fullmatch(r"[1-9]\d*[GM]", memory.upper()):
        raise ValueError("Invalid CPUs or memory")
    tasks_count = len(tasks(selected))
    print(f"{tasks_count} runs, 20 samples; WGS per sample CPU={cpus} RAM={memory}")
    if not args.confirm:
        print("DRY RUN: no jobs submitted. Add --confirm to submit all dependencies.")
        return
    selection = args.selection.expanduser().resolve()
    jobs = {}
    try:
        jobs["download"] = sbatch(args.run_id, "pilot_dl", job_script(args.run_id,"download",selection,array=True),
                                   2,"4G",args.download_time,
                                   array=f"0-{tasks_count - 1}%{args.download_parallel}",partition=args.partition)
        for sp in SPECIES_DIR:
            jobs[f"index_{sp}"] = sbatch(args.run_id, f"pilot_idx_{sp}",
                                         job_script(args.run_id,"index",selection,species=sp),
                                         args.index_cpus,args.index_memory,args.index_time,partition=args.partition)
        jobs["manifests"] = sbatch(args.run_id,"pilot_manifest",job_script(args.run_id,"manifest",selection),
                                   1,"4G","01:00:00",depend=jobs["download"],partition=args.partition)
        for sp in SPECIES_DIR:
            jobs[f"wgs_{sp}"] = sbatch(args.run_id,f"pilot_wgs_{sp}",
                                      job_script(args.run_id,"sample",selection,sp,array=True),
                                      cpus,memory,args.wgs_time,
                                      depend=f"{jobs['manifests']}:{jobs['index_' + sp]}",
                                      array=f"0-9%{args.wgs_parallel}",partition=args.partition)
            jobs[f"joint_{sp}"] = sbatch(args.run_id,f"pilot_joint_{sp}",
                                        job_script(args.run_id,"joint",selection,sp),
                                        args.joint_cpus,args.joint_memory,args.joint_time,
                                        depend=jobs[f"wgs_{sp}"],partition=args.partition)
    except Exception as exc:
        atomic_text(wd / "jobs.partial.json", json.dumps(jobs, indent=2) + "\n")
        raise RuntimeError(f"Partial SLURM submission; inspect {wd}/jobs.partial.json before resubmitting: {exc}") from exc
    atomic_text(wd / "jobs.json", json.dumps(jobs, indent=2) + "\n")
    print("Submitted pipeline:", wd / "jobs.json")
    print(json.dumps(jobs, indent=2))


def status(run_id: str) -> None:
    """Статус по squeue и контрольным файлам, без sacct/slurmdbd."""
    wd = workdir(run_id)
    jobs_file = wd / "jobs.json"
    if not jobs_file.is_file():
        raise FileNotFoundError(jobs_file)

    print("=== Задания SLURM ===")
    for name, jobid in json.loads(jobs_file.read_text()).items():
        result = subprocess.run(
            ["squeue", "-h", "-j", str(jobid), "-o", "%i|%T|%R"],
            capture_output=True, text=True,
        )
        print(f"\n{name} (job {jobid})")
        if result.returncode != 0:
            print("Ошибка squeue:", result.stderr.strip())
            continue

        lines = result.stdout.splitlines()
        if not lines:
            print("Не находится в очереди; исход проверяйте по логам.")
        else:
            for line in lines[:5]:
                print(" ", line.replace("|", "  "))
            if len(lines) > 5:
                print(f"  ...ещё {len(lines) - 5} строк")

    print("\n=== Результаты ===")
    verified_dir = wd / "verified"
    count = len(list(verified_dir.glob("*.json"))) if verified_dir.is_dir() else 0
    print(f"Проверено SRA runs: {count}/36")

    for species in SPECIES_DIR:
        outcome = wd / f"{species}.complete.ok"
        if outcome.is_file():
            print(species, "VCF готов:", outcome.read_text().strip())
        else:
            print(species, "финальный VCF ещё не подтверждён")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    actions = p.add_subparsers(dest="action", required=True)
    for command in ("plan", "submit", "worker"):
        q = actions.add_parser(command)
        q.add_argument("--selection", required=True, type=Path)
        if command != "plan":
            q.add_argument("--run-id", required=True)
        if command == "plan":
            q.add_argument("--check-ena", action="store_true", help="Verify live ENA FASTQ availability for all 36 runs")
        if command == "worker":
            q.add_argument("--stage", required=True, choices=["download","index","manifest","sample","joint"])
            q.add_argument("--species", choices=list(SPECIES_DIR))
            q.add_argument("--index", type=int, default=0)
        if command == "submit":
            q.add_argument("--confirm", action="store_true", help="Actually call sbatch; otherwise dry run")
            q.add_argument("--auto-resources", action="store_true", help="Suggest WGS allocation from current SLURM idle resources")
            q.add_argument("--max-cpus", type=int, default=32)
            q.add_argument("--max-memory-gib", type=int, default=96)
            q.add_argument("--cpus", type=int, default=16)
            q.add_argument("--memory", default="64G")
            q.add_argument("--index-cpus", type=int, default=8)
            q.add_argument("--index-memory", default="64G")
            q.add_argument("--joint-cpus", type=int, default=8)
            q.add_argument("--joint-memory", default="64G")
            q.add_argument("--download-parallel", type=int, default=3)
            q.add_argument("--wgs-parallel", type=int, default=1)
            q.add_argument("--partition", default="")
            q.add_argument("--download-time", default="2-00:00:00")
            q.add_argument("--index-time", default="1-00:00:00")
            q.add_argument("--wgs-time", default="7-00:00:00")
            q.add_argument("--joint-time", default="7-00:00:00")
            q.add_argument("--allow-low-disk", action="store_true")
    s = actions.add_parser("status")
    s.add_argument("--run-id", required=True)
    return p


def main() -> None:
    args = parser().parse_args()
    if args.action == "status":
        status(args.run_id)
        return
    selected = table(args.selection)
    if args.action == "plan":
        print("Selected animals", len(selected), "SRA runs", len(tasks(selected)))
        for sp in SPECIES_DIR:
            s = species_samples(selected, sp)
            print(sp, len(s), "samples", sum(len(runs_for(r)) for r in s), "runs",
                  "studies", len({r["bioproject"] for r in s}), "bioprojects")
        disk_report(selected)
        capacity_snapshot()
        if args.check_ena:
            total_bytes = 0
            for sp, sample, run in tasks(selected):
                meta = eno_run(run)
                byte_count = sum(pair["r1"]["bytes"] + pair["r2"]["bytes"] for pair in meta["pairs"])
                total_bytes += byte_count
                print(f"ENA OK {sp} {sample} {run}: {len(meta['pairs'])} pair(s), {byte_count / 1024**3:.2f} GiB")
            print(f"ENA actual total compressed FASTQ estimate: {total_bytes / 1024**3:.1f} GiB")
        return
    if args.action == "submit":
        for name in ("download_parallel", "wgs_parallel", "index_cpus", "joint_cpus"):
            if getattr(args, name) < 1:
                raise ValueError(f"{name} must be >= 1")
        submit(selected, args)
        return
    # Never allow queued workers to consume a modified cohort selection.
    lock = workdir(args.run_id) / "selection.lock.json"
    if not lock.is_file():
        raise FileNotFoundError(f"Selection lock missing: {lock}")
    recorded = json.loads(lock.read_text())
    if hashlib.sha256(args.selection.read_bytes()).hexdigest() != recorded["sha256"]:
        raise ValueError("Cohort selection changed after jobs were submitted")
    if args.stage == "download":
        download_task(selected, args.run_id, args.index)
    elif args.stage == "index":
        index_reference(args.species)
    elif args.stage == "manifest":
        make_manifests(selected, args.run_id)
    elif args.stage == "sample":
        sample_task(selected, args.run_id, args.species, args.index)
    elif args.stage == "joint":
        joint_task(selected, args.run_id, args.species)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        sys.exit(1)
