#!/usr/bin/env python3
"""Validate archived arrays, proposal arithmetic and explicit source hashes.

Uses the Python standard library. Does not import/load MACE/GPAW or run inference,
DFT, MD or training. Reads only this delivery and named public dependencies.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[5]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def table(name):
    with (OUT / name).open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), name
    assert all(None not in row and None not in row.values() for row in rows), name
    assert b"\r\n" not in (OUT / name).read_bytes(), name
    return rows


def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=2e-10, abs_tol=2e-10), (a, b)


atoms = table("interface_per_atom.csv")
metrics = table("interface_metrics.csv")
elements = table("interface_per_element.csv")
response = table("interface_response_comparison.csv")
assert (len(atoms), len(metrics), len(elements), len(response)) == (312, 6, 24, 14)
groups = defaultdict(list)
keys = set()
for atom in atoms:
    key = atom["label"], atom["model"], int(atom["atom_id"])
    assert key not in keys
    keys.add(key)
    groups[key[:2]].append(atom)
    errors = [float(atom[f"model_F{axis}_eV_A"]) - float(atom[f"DFT_F{axis}_eV_A"])
              for axis in "xyz"]
    for axis, err in zip("xyz", errors):
        close(err, atom[f"error_F{axis}_eV_A"])
    close(math.sqrt(sum(x * x for x in errors)), atom["error_norm_eV_A"])


def check_metrics(rows, recorded):
    count = len(rows)
    close(math.sqrt(sum(float(r["error_norm_eV_A"]) ** 2 for r in rows) / count),
          recorded["force_vector_RMSE_eV_A"])
    close(max(float(r["error_norm_eV_A"]) for r in rows),
          recorded["max_atom_force_vector_error_eV_A"])
    for axis in "xyz":
        close(math.sqrt(sum(float(r[f"error_F{axis}_eV_A"]) ** 2 for r in rows) / count),
              recorded[f"{axis}_RMSE_eV_A"])


by_shift = {}
for m in metrics:
    rows = groups[m["label"], m["model"]]
    assert len(rows) == int(m["n_atoms"]) == 52
    assert sorted(int(r["index_zero_based"]) for r in rows) == list(range(52))
    assert [r["element"] for r in rows].count("Ag") == 4
    check_metrics(rows, m)
    ag = [r for r in rows if r["element"] == "Ag"]
    for method, prefix in [("DFT", "DFT"), ("model", "model")]:
        gz = sum(float(r[f"{prefix}_Fz_eV_A"]) for r in ag)
        close(gz, m[f"whole_Ag_layer_{method}_Gz_eV_A"])
        close(gz / 4, m[f"mean_Ag_{method}_Fz_eV_A"])
    close(float(m["pair_model_projection_eV_A"]) - float(m["pair_DFT_projection_eV_A"]),
          m["pair_signed_error_eV_A"])
    close(float(m["pair_DFT_projection_eV_A"]) / 2, m["pair_half_coordinate_DFT_eV_A"])
    close(float(m["pair_model_projection_eV_A"]) / 2, m["pair_half_coordinate_model_eV_A"])
    by_shift[float(m["ag_shift_z_A"]), m["model"]] = m

for elem in elements:
    rows = [r for r in groups[elem["label"], elem["model"]] if r["element"] == elem["element"]]
    assert len(rows) == int(elem["n_atoms"])
    check_metrics(rows, elem)

for r in response:
    model = r["model"]
    minus, plus = by_shift[-0.2, model], by_shift[0.2, model]
    assert float(r["span_A"]) == 0.4
    if r["quantity"] == "rigid_Ag_layer_z_restoring_force_secant":
        for method, output in [("DFT", "DFT_value"), ("model", "model_value")]:
            value = -(float(plus[f"whole_Ag_layer_{method}_Gz_eV_A"])
                      - float(minus[f"whole_Ag_layer_{method}_Gz_eV_A"])) / 0.4
            close(value, r[output])
    if model == "MACE-MP-0b3-medium" and "energy" in r["quantity"]:
        assert r["model_value"] == r["status"] == "NOT_ARCHIVED"

metadata = json.loads((OUT / "inference_metadata.json").read_text())
assert metadata["model_sha256_before"] == metadata["model_sha256_after"]
assert metadata["inference_counts"] == {
    "new_MP_configurations": 2, "reused_C_MP_configurations": 1,
    "reused_A_V18_configurations": 3, "new_V18_configurations": 0,
    "new_backbone_configurations": 0,
}
for paths in [metadata["public_source_dependencies_SHA256"]]:
    for rel, expected in paths.items():
        assert digest(ROOT / rel) == expected, rel
for model, rel in metadata["model_paths"].items():
    assert digest(ROOT / rel) == metadata["model_sha256_after"][model]
independent = json.loads((OUT / "independent_metric_review.json").read_text())
assert independent["status"] == "PASS" and all(independent["checks"].values())
for name, expected in independent["source_sha256"].items():
    assert digest(OUT / name) == expected, name
plan_review = json.loads((OUT / "independent_plan_review.json").read_text())
assert plan_review["result"] == "PASS" and not plan_review["issue"]
for name, expected in plan_review["evidence_sha256"].items():
    assert digest(OUT / name) == expected, name

for dest, source in [
    ("backbone_response_reused.csv", "response_comparison.csv"),
    ("backbone_metrics_reused.csv", "per_configuration_metrics.csv"),
]:
    original = ROOT / "research/mace-v12-transfer/coordination/reports/window-c/ti_c_backbone_response_final_audit_20261011" / source
    assert digest(OUT / dest) == digest(original), dest

pure = table("pure_phase_role_audit.csv")
assert len(pure) == 26 and len({r["label"] for r in pure}) == 26
for r in pure:
    assert r["train_eligible_now"] == r["role_change_authorized_this_turn"] == "NO"
    assert r["scf_converged"] == r["label_energy_both_present"] == "True"
    assert r["archive_SHA"] == r["archive_verifier"] == "PASS"
    assert r["dataset_role_as_published"] == "numerical_pure_phase_control_not_training_or_independent_validation"
for material, count, geometries in [("solid_Ag", 15, 7), ("solid_Ti3SiC2", 11, 5)]:
    entries = [r for r in pure if r["material"] == material]
    assert len(entries) == count and len({r["source_sha256"] for r in entries}) == geometries

proposals = table("candidate_data_requirements.csv")
assert len({r["proposal_group"] for r in proposals}) == len(proposals)
minimum = [r for r in proposals if r["proposal_group"].startswith(("P-MIN-DEV-", "P-MIN-TEST-"))]
assert len(minimum) == 8
for prefix in ["P-MIN-DEV-", "P-MIN-TEST-"]:
    assert sum(int(r["proposed_DFT_labels"]) for r in minimum if r["proposal_group"].startswith(prefix)) == 11
assert sum(int(r["proposed_DFT_labels"]) for r in minimum) == 22
assert sum(int(r["proposed_parent_roots"]) for r in minimum) == 8
assert all(r["assignment"] == "NOT_ASSIGNED" for r in proposals)
assert all(r["authorization"] in ["NOT_AUTHORIZED", "ROLE_CHANGE_NOT_AUTHORIZED"] for r in proposals)
liquid = json.loads((OUT / "liquid_pilot_preservation_review.json").read_text())
assert liquid["scientific_runs_started"] == 0
assert not liquid["sealed_structures_or_labels_read_or_enumerated"]
prov = json.loads((OUT / "PROVENANCE.json").read_text())
assert prov["new_operations"]["registered_compute_tasks"] == 0
assert prov["new_operations"]["changed_data_roles"] == 0
assert not prov["sampler_and_DFT_pilot_gates_changed"]
assert not prov["sealed_structure_or_label_reads"] and not prov["sealed_structure_enumeration"]

checksum = OUT / "SHA256SUMS.txt"
if checksum.exists():
    for line in checksum.read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert Path(name).name == name, name
        assert digest(OUT / name) == expected, name

print(json.dumps({"status": "PASS", "interface_rows": [312, 24, 6, 14],
                  "pure_records": 26, "training_eligible_pure_records": 0,
                  "minimum_proposed_labels": 22, "registered_new_tasks": 0,
                  "scientific_operations_by_validator": 0}, ensure_ascii=False))
