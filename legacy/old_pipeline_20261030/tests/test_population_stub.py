import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("population_prepare", ROOT / "population_structure" / "prepare.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class TestPopulationStub(unittest.TestCase):
    def test_plans_only_no_inferred_labels(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            vcf = root / "variants.vcf.gz"
            vcf.write_bytes(b"test placeholder file, not a real vcf")
            out = root / "plan.json"
            self.assertEqual(module.main([
                "--vcf", str(vcf), "--species", "Ursus_maritimus",
                "--assembly", "GCA_062004589.1", "--out", str(out)]), 0)
            p = json.loads(out.read_text())
            self.assertEqual(p["status"], "not_implemented")
            self.assertIsNone(p["population_labels_tsv"])
            self.assertEqual(len(p["vcfs"]), 1)
            self.assertFalse((root / "population_labels.tsv").exists())

    def test_duplicate_variants_error(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.bcf"
            p.write_bytes(b"\x00")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                module.make_plan([str(p), str(p)], "species", "assembly")

    def test_missing_and_unknown_input(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.csv"
            p.write_text("none")
            with self.assertRaisesRegex(ValueError, "Unsupported"):
                module.make_plan([str(p)], "species", "assembly")
            with self.assertRaisesRegex(ValueError, "not found"):
                module.make_plan([str(p) + ".vcf"], "species", "assembly")


if __name__ == "__main__":
    unittest.main()
