"""Unit tests for orchestration contracts. No actual ENA/SLURM/WGS execution."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from argparse import Namespace
from unittest.mock import patch

script = Path(__file__).resolve().parent.parent / "tools" / "wgs_pilot.py"
if not script.exists():
    script = Path(__file__).with_name("wgs_pilot.py")
spec = importlib.util.spec_from_file_location("pilot", script)
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)
selection = Path(__file__).resolve().parent.parent / "inputs" / "wgs_pilot_20.tsv"
if not selection.exists():
    selection = Path(__file__).with_name("wgs_pilot_20.tsv")


class PilotTest(unittest.TestCase):
    def test_selected_cohort_contract(self):
        rows = pilot.table(selection)
        self.assertEqual(len(rows), 20)
        self.assertEqual(len(pilot.tasks(rows)), 36)
        self.assertEqual(len({r["bioproject"] for r in rows}), 9)
        self.assertEqual([len(pilot.runs_for(r)) for r in rows if r["sample_id"] == "SAMN12238539"], [8])

    def test_ena_pair_field_parsing(self):
        import urllib.request
        md5 = hashlib.md5(b"mock data bytes").hexdigest()
        line = ("run_accession\texperiment_accession\tlibrary_layout\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
                "SRR123\tSRX321\tPAIRED\t"
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123_1.fastq.gz;"
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123_2.fastq.gz\t"
                f"{md5};{md5}\t15;15\n")
        response = io.BytesIO(line.encode())
        class Response:
            def __enter__(self): return response
            def __exit__(self, *a): return False
        with patch.object(pilot, "urlopen", return_value=Response()):
            output = pilot.eno_run("SRR123")
        self.assertEqual(output["library"], "SRX321")
        self.assertEqual(len(output["pairs"]), 1)
        self.assertTrue(output["pairs"][0]["r1"]["url"].startswith("https://ftp.sra.ebi.ac.uk/"))

    def test_existing_fastq_verifies_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "SRR123_1.fastq.gz"
            target.write_bytes(b"some local bytes")
            meta = {"bytes": target.stat().st_size, "md5": hashlib.md5(target.read_bytes()).hexdigest()}
            pilot.download_file(meta, target)
            target.write_bytes(b"different content")
            with self.assertRaisesRegex(ValueError, "disagrees with ENA"):
                pilot.download_file(meta, target)

    def test_manifest_preserves_multi_runs_and_libraries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch.object(pilot, "RESULTS", root), patch.object(pilot, "DATA", root / "data"):
                rows = pilot.table(selection)
                run_id = "pilot_test"
                wd = pilot.workdir(run_id)
                for sp, sample, run in pilot.tasks(rows):
                    p1 = root / "data" / f"{run}_1.fastq.gz"
                    p2 = root / "data" / f"{run}_2.fastq.gz"
                    p1.parent.mkdir(parents=True, exist_ok=True)
                    p1.touch(); p2.touch()
                    meta = [{"sample_id": sample, "species": sp, "run_id": run,
                             "library": "SRX" + run[3:], "fastq_r1": str(p1), "fastq_r2": str(p2)}]
                    pilot.atomic_text(wd / "verified" / f"{run}.json", json.dumps(meta) + "\n")
                pilot.make_manifests(rows, run_id)
                path = wd / "manifests" / "beluga" / "SAMN12238539.tsv"
                content = path.read_text().splitlines()
                self.assertEqual(len(content), 9)
                self.assertTrue((wd / "manifests" / "complete.ok").is_file())
                # same selected sample -> same stable manifest on repeated invocation
                pilot.make_manifests(rows, run_id)

    def test_submission_builds_array_dependencies_without_contacting_slurm(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(pilot, "RESULTS", Path(tmp)), \
                 patch.object(pilot, "disk_report", return_value=True), \
                 patch.object(pilot.shutil, "which", return_value="/usr/bin/dummy"), \
                 patch.object(pilot, "sbatch", side_effect=[str(x) for x in range(100, 108)]) as sb:
                args = Namespace(run_id="test_batch", selection=selection, auto_resources=False,
                                 cpus=16, memory="64G", index_cpus=8, index_memory="64G",
                                 joint_cpus=8, joint_memory="64G", download_time="2-00:00:00",
                                 index_time="1-00:00:00", wgs_time="7-00:00:00",
                                 joint_time="7-00:00:00", download_parallel=3, wgs_parallel=1,
                                 partition="", confirm=True, allow_low_disk=False)
                pilot.submit(pilot.table(selection), args)
                jobs = json.loads((Path(tmp)/"pilot_orchestration"/"test_batch"/"jobs.json").read_text())
                self.assertEqual(len(jobs), 8)
                self.assertEqual(jobs["download"], "100")
                self.assertEqual(jobs["joint_polar_bear"], "107")
                # Manifest is after the complete download array; each WGS array
                # waits for the matching index; joint waits for whole WGS array.
                self.assertEqual(sb.call_args_list[3].kwargs["depend"], "100")
                self.assertEqual(sb.call_args_list[4].kwargs["depend"], "103:101")
                self.assertEqual(sb.call_args_list[5].kwargs["depend"], "104")
                self.assertEqual(sb.call_args_list[6].kwargs["depend"], "103:102")
                self.assertEqual(sb.call_args_list[7].kwargs["depend"], "106")
                self.assertEqual(sb.call_args_list[0].kwargs["array"], "0-35%3")

    def test_array_job_script_and_repeated_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(pilot, "RESULTS", Path(tmp)):
                a = pilot.job_script("pilot_1", "sample", selection, "beluga", array=True)
                content = a.read_text()
                self.assertIn('"${SLURM_ARRAY_TASK_ID:?}"', content)
                self.assertIn("--stage sample", content)
                pilot.job_script("pilot_1", "sample", selection, "beluga", array=True)


if __name__ == "__main__":
    unittest.main()
