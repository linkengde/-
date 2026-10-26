#!/usr/bin/env python3
"""Preflight the separately registered Window A B1 mesh task; read-only."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from ase.io import read

ROOT = Path(__file__).resolve().parents[6]
HERE = Path(__file__).resolve().parent
BASE = ROOT / "research/mace-v12-transfer"
TASK = json.loads((HERE / "task_manifest.json").read_text())
SOURCE_PACK = ROOT / TASK["source_pack"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

assert TASK["schema"] == "WINDOW_A_B1_PHASE_KMESH_TASK_V1"
assert TASK["task_id"] == "window-a-b1-kmesh-20261026"
assert TASK["status"] == "registered_ready"
assert TASK["global_dft_authorized"] is True and TASK["global_launch_enabled"] is True
assert TASK["training_enabled"] is False and TASK["dataset_integration_enabled"] is False
assert TASK["method_decision"]["status"] == "FROZEN_FOR_THIS_MESH_TASK_ONLY"
assert sha(SOURCE_PACK) == TASK["source_pack_sha256"]
source = json.loads(SOURCE_PACK.read_text())
assert len(TASK["records"]) == 12 == source["record_count"]
assert TASK["resource_cap"] == {
    "mpi_ranks_per_endpoint": 4,
    "maximum_concurrent_endpoints": 1,
    "gpu_count": 0,
    "mandatory_endpoints": 12,
    "conditional_endpoints_maximum": 2,
    "absolute_endpoint_maximum": 14,
    "automatic_retries": 0,
    "execution_order": "Sequential; run each coarse/fine pair and compare it before moving to the next pair.",
}
ledger = BASE / "coordination/reports/window-c/b1_candidate_role_freeze_20261024/role_ledger.csv"
roles = {r["config_type"]: r for r in csv.DictReader(ledger.open())}
pairs = defaultdict(list)
for record, original in zip(TASK["records"], source["records"]):
    assert record["label"] == original["label"]
    for key in ("input", "input_sha256", "config_type", "natoms", "pair_id", "pair_member", "kpts", "method"):
        assert record[key] == original[key], (record["label"], key)
    assert record["candidate_role"] == "TRAINING_CANDIDATE_ONLY"
    assert record["launch_enabled"] is True and record["dft_authorized"] is True
    assert record["dataset_integration_enabled"] is False
    input_path = (ROOT / record["input"]).resolve()
    assert input_path.is_relative_to(ROOT)
    assert sha(input_path) == record["input_sha256"]
    atoms = read(input_path)
    assert len(atoms) == record["natoms"] and atoms.pbc.tolist() == record["pbc"]
    assert not any(k in atoms.info for k in ("REF_energy", "PW_PBE_energy_eV", "free_energy_eV", "energy"))
    assert not any(k in atoms.arrays for k in ("REF_forces", "PW_PBE_forces", "forces"))
    role = roles[record["config_type"]]
    assert role["file_sha256"] == record["input_sha256"]
    assert role["candidate_role"] == "TRAINING_CANDIDATE_ONLY"
    assert role["dev_test_role"] == "NONE"
    pairs[record["pair_id"]].append(record)

assert len(pairs) == 6
for pair_id, members in pairs.items():
    assert len(members) == 2 and {r["pair_member"] for r in members} == {"coarse", "fine"}
    assert len({r["input_sha256"] for r in members}) == 1
    for key in ("input", "input_sha256", "config_type", "natoms", "pbc"):
        assert members[0][key] == members[1][key]
    m0, m1 = (dict(r["method"]) for r in members)
    m0.pop("kpts"); m1.pop("kpts")
    assert m0 == m1

assert len(TASK["conditional_records"]) == 2
for record in TASK["conditional_records"]:
    assert record["phase"] == "liquid_Ag" and record["kpts"] == [3, 3, 3]
    assert record["launch_enabled"] is False
    assert record["compare_to_label"] in {r["label"] for r in TASK["records"]}
    assert sha((ROOT / record["input"]).resolve()) == record["input_sha256"]

print(json.dumps({"status":"PASS_B1_TASK_PREFLIGHT","mandatory_endpoints":12,"same_geometry_pairs":6,"conditional_endpoints":2,"launch_enabled":True,"training_or_integration":False}, indent=2))
