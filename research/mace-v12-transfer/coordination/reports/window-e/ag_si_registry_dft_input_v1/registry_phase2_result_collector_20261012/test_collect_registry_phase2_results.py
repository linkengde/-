import importlib.util
import hashlib
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("collect_registry_phase2_results.py")
SPEC = importlib.util.spec_from_file_location("phase2_collector", SCRIPT)
collector = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(collector)

RUNNER_SPEC = importlib.util.spec_from_file_location("phase2_runner", SCRIPT.with_name("phase2_runner.py"))
runner = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
RUNNER_SPEC.loader.exec_module(runner)


class CollectorParsingTests(unittest.TestCase):
    def test_phase2_json_collects_all_label_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "R2"
            out = case / "results" / "medium"
            out.mkdir(parents=True)
            (case / "input.py").write_text("# frozen input\n")
            (case / "structure.xyz").write_text("0\nempty test fixture\n")
            (case / "C_EXECUTION_RELEASE_RECORD.json").write_text(json.dumps({
                "configuration_sha256": "c" * 64,
                "kpoint_mesh": [8, 8, 1],
                "scf_converged": True,
                "scf_steps": 12,
            }))
            (out / "result.json").write_text(json.dumps({
                "total_energy_eV": -1.25,
                "force_consistent_free_energy_eV": -1.3,
                "forces_eV_A": [[0.1, 0.2, 0.3]],
                "stress_eV_A3": [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
            }))
            item = collector.collect_one("R2", root)
            self.assertEqual(item["data_status"], "RESULTS_COLLECTED")
            self.assertEqual(item["config_sha"], "c" * 64)
            self.assertEqual(item["energy"], -1.25)
            self.assertEqual(item["free_energy"], -1.3)
            self.assertEqual(item["forces_status"], "AVAILABLE")
            self.assertEqual(item["stress_status"], "AVAILABLE")
            self.assertEqual(item["scf_steps"], 12)

    def test_other_mesh_is_not_misreported_as_phase2(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "R1"
            out = case / "results" / "coarse"
            out.mkdir(parents=True)
            (case / "input.py").write_text("# test\n")
            (case / "structure.xyz").write_text("0\nfixture\n")
            (out / "result.json").write_text(json.dumps({
                "kpoint_mesh": [6, 6, 1], "total_energy_eV": -2.0,
                "force_consistent_free_energy_eV": -2.1,
            }))
            item = collector.collect_one("R1", root)
            self.assertEqual(item["data_status"], "KPOINT_MISMATCH")
            self.assertIsNone(item["energy"])
            self.assertIsNone(item["free_energy"])

    def test_gpaw_log_parser_reads_energy_forces_and_scf(self):
        text = """kpts=MonkhorstPack(size=array([8, 8, 1]))
Converged in 9 steps
Free energy: -3.500
Extrapolated: -3.250
Forces  # eV/Ang
|        x |        y |        z |
|----------|----------|----------|
|  0.1 | 0.2 | 0.3 | # Ag 0
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gpaw.log"
            path.write_text(text)
            parsed = collector.parse_log(path)
        self.assertEqual(parsed["kpoint"], [8, 8, 1])
        self.assertTrue(parsed["scf_converged"])
        self.assertEqual(parsed["scf_steps"], 9)
        self.assertEqual(parsed["energy"], -3.25)
        self.assertEqual(parsed["free_energy"], -3.5)
        self.assertEqual(parsed["forces"], [[0.1, 0.2, 0.3]])

    def test_input_hash_mismatch_blocks_result_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "R2"
            out = case / "results" / "medium"
            out.mkdir(parents=True)
            (case / "input.py").write_text("# actual input\n")
            (case / "structure.xyz").write_text("0\nfixture\n")
            (case / "C_EXECUTION_RELEASE_RECORD.json").write_text(json.dumps({
                "input_sha256": "0" * 64,
                "configuration_sha256": "c" * 64,
                "kpoint_mesh": [8, 8, 1],
            }))
            (out / "result.json").write_text(json.dumps({
                "total_energy_eV": -1.25,
                "force_consistent_free_energy_eV": -1.3,
            }))
            item = collector.collect_one("R2", root)
        self.assertEqual(item["data_status"], "PROVENANCE_MISMATCH")
        self.assertIsNone(item["energy"])
        self.assertIsNone(item["free_energy"])

    def test_auth_request_binds_runner_collector_inputs_and_structures(self):
        request = json.loads(SCRIPT.with_name("PHASE2_EXEC_AUTH_REQUEST.json").read_text())
        actual = runner.requested_hash_bindings()
        self.assertEqual(request["bindings"], actual)
        preflight = runner.preflight()
        self.assertFalse(preflight["execution_authorized"])
        self.assertTrue(preflight["source_hashes_match"])

    def test_prepared_not_issued_auth_cannot_pass_execution_gate(self):
        auth = {"authorization_status": "PREPARED_NOT_ISSUED", "execution_authorized": False}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "REGISTRY_PHASE2_EXEC_AUTH.json"
            raw = json.dumps(auth).encode()
            path.write_bytes(raw)
            digest = hashlib.sha256(raw).hexdigest()
            with self.assertRaises(SystemExit):
                runner.verify_auth(path, digest)

    def test_summary_json_has_required_phase2_fields_and_only_reviewed_ids(self):
        summary = json.loads(SCRIPT.with_name("REGISTRY_PHASE2_RESULTS.json").read_text())
        self.assertEqual([item["structure_id"] for item in summary["results"]], ["R2", "R1"])
        required = {
            "structure_id", "input_sha", "structure_sha", "config_sha", "kpoint",
            "energy", "free_energy", "forces_status", "stress_status", "scf_steps",
        }
        for item in summary["results"]:
            self.assertTrue(required.issubset(item))


if __name__ == "__main__":
    unittest.main()
