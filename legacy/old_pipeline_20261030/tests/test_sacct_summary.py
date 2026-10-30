import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("sacct_summary", ROOT / "benchmarks" / "summarize_sacct.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class TestSacctSummary(unittest.TestCase):
    def test_cpu_parser(self):
        self.assertEqual(module.parse_duration("01:30:00"), 5400)
        self.assertEqual(module.parse_duration("1-00:00:01"), 86401)
        self.assertEqual(module.parse_duration("00:30"), 30)

    def test_memory_parser(self):
        self.assertAlmostEqual(module.parse_memory_gib("1024M"), 1.0)
        self.assertAlmostEqual(module.parse_memory_gib("1G"), 1.0)

    def test_aggregate_job_and_batch_memory(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "sacct.tsv"
            path.write_text(
                "23|align_lane|COMPLETED|3600|4|02:00:00|\n"
                "23.batch|batch|COMPLETED|3600|4|02:00:00|8192M\n"
                "24|genomicsdb|FAILED|600|2|00:02:00|512M\n")
            rows = module.summarize(path)
            align = next(r for r in rows if r["job_name"] == "align_lane")
            self.assertEqual(align["jobs"], 1)
            self.assertEqual(align["max_reported_rss_gib"], 8)
            self.assertAlmostEqual(align["median_cpu_eff_pct"], 50.0)
            gen = next(r for r in rows if r["job_name"] == "genomicsdb")
            self.assertEqual(gen["non_completed"], 1)
    def test_missing_cpu_is_unknown_not_zero(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "sacct.tsv"
            path.write_text("JobID|JobName|State|ElapsedRaw|AllocCPUS|ReqMem\n"
                            "1|WGS_align|COMPLETED|3600|16|32000M\n")
            row = module.summarize(path)[0]
            self.assertIsNone(row["total_used_cpu_h"])
            self.assertIsNone(row["median_cpu_eff_pct"])
            self.assertIsNone(row["max_reported_rss_gib"])


if __name__ == "__main__":
    unittest.main()
