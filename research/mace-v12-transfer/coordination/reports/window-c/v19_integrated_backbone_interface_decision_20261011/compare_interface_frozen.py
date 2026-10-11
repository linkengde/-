#!/usr/bin/env python3
"""Read public A archives, reuse verified forces, infer only two missing MP-0b3 cases.

No DFT, MD, training, network operation, or sealed-data access is performed.
Run with the project Python 3.12/MACE 0.3.16/PyTorch 2.5.1 CPU environment.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/cloud-c-integrated-interface-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/cloud-c-integrated-interface-cache")

import numpy as np
import torch
from ase.geometry import find_mic
from ase.io import read
from mace.calculators import MACECalculator

OUT = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("CLOUD_C_REPO_ROOT", OUT.parents[5])).resolve()
BASE = ROOT / "research/mace-v12-transfer"
GAP = BASE / "periodic_interface_v4/pbe_interface_v19_gap_scan"
A = GAP / "calculations/v18_gap_scan_force_screen"
C = BASE / "coordination/reports/window-c/three_model_static_comparison_20261010"
LABELS = [f"AgSi_rigid_Ag_z_{s}_k6x6x1_sigma0p10_dipoleXY_v19" for s in ("m020", "parent", "p020")]
SHIFTS = dict(zip(LABELS, [-0.2, 0.0, 0.2]))
FAMILY = "AgSi_rigid_gap_parent_A_v19"
MP = "MACE-MP-0b3-medium"
V18 = "V18"
MODEL_PATHS = {
    MP: BASE / "periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model",
    V18: BASE / "periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model",
}
MODEL_SHAS = {
    MP: "2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012",
    V18: "757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62",
}
EXPECTED_INPUTS = dict(zip(LABELS, [
    "ea1c62e43ece0339f7047ac85b09d95fe6c6f07cf3d58fd0eab444ce9a088c43",
    "73eff17b1913b653017647e9af9c9681e340f0b8067936419f3377963a3bc02c",
    "8c6ccb44585ba8b0e9db9777c657ae2e92286a280b042a5f9bc1665cb9637d01",
]))
EXPECTED_RESULTS = dict(zip(LABELS, [
    "e6cc182e99f82b1a1df52c874a4fa30d54dfb10fa1a65551615c99ba79de8cef",
    "b3528f0f3bcafa7f8e90b5a8ee64250c9b2364de9b2e215ed184b077839a786a",
    "7254a9dd11d98a07032a2b6acebb429c499e00322e997dd0288e4514319049d6",
]))


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def relative(p):
    return str(Path(p).relative_to(ROOT))


def load_csv(p):
    with Path(p).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows):
    assert rows and all(list(row) == list(rows[0]) for row in rows)
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_selected_inventory(folder, names):
    inventory = {}
    for line in (folder / "SHA256SUMS.txt").read_text().splitlines():
        if line.strip():
            expected, name = line.split(maxsplit=1)
            inventory[name.lstrip("*")] = expected
    result = {}
    for name in names:
        assert name in inventory
        p = (folder / name).resolve()
        assert p.is_relative_to(folder.resolve()) and p.is_file()
        actual = sha(p)
        assert actual == inventory[name], f"Published dependency changed: {name}"
        result[relative(p)] = actual
    result[relative(folder / "SHA256SUMS.txt")] = sha(folder / "SHA256SUMS.txt")
    return result


def errors(delta):
    norm = np.linalg.norm(delta, axis=1)
    axis = np.sqrt(np.mean(delta**2, axis=0))
    total = float(np.sum(delta**2))
    return {
        "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(norm**2))),
        "x_RMSE_eV_A": float(axis[0]), "y_RMSE_eV_A": float(axis[1]),
        "z_RMSE_eV_A": float(axis[2]),
        "max_atom_force_vector_error_eV_A": float(norm.max()),
        "z_squared_error_fraction": float(np.sum(delta[:, 2]**2) / total) if total else 0.0,
    }


def pair_geometry(atoms):
    marked = np.flatnonzero(atoms.arrays["central_pair"])
    ag = [int(i) for i in marked if atoms[i].symbol == "Ag"]
    si = [int(i) for i in marked if atoms[i].symbol == "Si"]
    assert len(marked) == 2 and len(ag) == len(si) == 1
    i, j = ag[0], si[0]
    vector, distance = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, atoms.pbc)
    assert float(distance) > 0
    return i, j, vector / float(distance), float(distance)


def main():
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    torch.set_default_dtype(torch.float64)
    assert sys.version_info[:2] == (3, 12)
    assert importlib.metadata.version("mace-torch") == "0.3.16"
    assert torch.__version__ == "2.5.1+cpu"
    assert importlib.metadata.version("ase") == "3.29.0"
    hashes_before = {model: sha(path) for model, path in MODEL_PATHS.items()}
    assert hashes_before == MODEL_SHAS
    public_dependencies = verify_selected_inventory(A, ["metrics.json", "per_atom_force_errors.csv", "analyze_gap_scan.py"])
    public_dependencies.update(verify_selected_inventory(C, ["per_atom_errors.csv", "reference_inventory.csv", "model_inventory.csv", "structure_metrics.csv", "compare_static.py"]))
    manifest = json.loads((GAP / "input_manifest.json").read_text())
    records = {r["label"]: r for r in manifest["records"] if r["label"] in LABELS}
    assert set(records) == set(LABELS)
    archives, atoms_by_label, summaries = {}, {}, {}
    for label in LABELS:
        record = records[label]
        folder = GAP / "calculations" / label
        verify = json.loads(subprocess.check_output([sys.executable, "-B", str(GAP / "verify_result.py"), label, str(folder)], text=True, env=os.environ.copy()))
        assert verify["status"] == "PASS" and len(verify["checks"]) == 18 and all(verify["checks"].values())
        result = folder / f"{label}_PW_PBE.extxyz"
        source = GAP / record["input"]
        assert sha(source) == EXPECTED_INPUTS[label] == record["input_sha256"]
        assert sha(result) == EXPECTED_RESULTS[label]
        parent = ROOT / record["parent_geometry_path"]
        assert sha(parent) == record["parent_geometry_sha256"]
        atoms, original = read(result), read(parent)
        assert original.get_chemical_symbols() == atoms.get_chemical_symbols()
        assert np.array_equal(original.arrays["lammps_id"], atoms.arrays["lammps_id"])
        assert np.array_equal(original.arrays["central_pair"], atoms.arrays["central_pair"])
        parent_positions = original.positions.copy()
        parent_positions[:, 2] += record["global_z_translation_A"]
        parent_positions[np.array(original.get_chemical_symbols()) == "Ag", 2] += SHIFTS[label]
        parent_cell = original.cell.array.copy()
        parent_cell[2, 2] += record["vacuum_expansion_A"]
        assert np.allclose(atoms.positions, parent_positions, rtol=0, atol=1e-12)
        assert np.allclose(atoms.cell.array, parent_cell, rtol=0, atol=1e-12)
        assert np.array_equal(atoms.arrays["lammps_id"], record["ordered_ids"])
        assert np.array_equal(atoms.arrays["central_pair"], record["pair_markings"])
        assert atoms.get_chemical_symbols() == record["ordered_symbols"]
        assert atoms.pbc.tolist() == [True, True, False] == record["pbc"]
        assert len(atoms) == 52 and sum(atoms.numbers == 47) == 4
        summary = json.loads((folder / "summary.json").read_text())
        assert summary["dataset_role"] == record["role"]
        atoms_by_label[label], summaries[label] = atoms, summary
        archives[label] = {"result_path": relative(result), "source_path": relative(source), "source_sha256": sha(source), "result_sha256": sha(result), "summary_sha256": sha(folder / "summary.json"), "parent_geometry_path": relative(parent), "parent_geometry_sha256": sha(parent), "ordered_atom_ids": atoms.arrays["lammps_id"].tolist(), "ordered_symbols": atoms.get_chemical_symbols(), "pbc": atoms.pbc.tolist(), "cell_A": atoms.cell.array.tolist(), "coordinate_order": "unchanged source/result lammps_id order", "ag_layer_shift_z_A": SHIFTS[label], "dataset_role": record["role"], "archive_verifier": verify, "source_parent_transform_verified": True}
    method = summaries[LABELS[0]]["method"]
    assert all(s["method"] == method for s in summaries.values())
    assert all(s["energy_convention"] == manifest["energy_convention"] for s in summaries.values())
    assert method["xc"] == "PBE" and method["cutoff_eV"] == 500 and method["kpts"] == [6, 6, 1] and method["smearing_eV"] == 0.10
    ids = atoms_by_label[LABELS[1]].arrays["lammps_id"]
    assert all(np.array_equal(a.arrays["lammps_id"], ids) for a in atoms_by_label.values())

    # Full precision A JSON is authoritative; its 12-digit CSV is independently cross-checked.
    a_report = json.loads((A / "metrics.json").read_text())
    assert a_report["candidate"]["model_sha256"] == MODEL_SHAS[V18]
    assert all(a_report["candidate"]["lineage_checks"].values())
    a_by_label = {r["label"]: r for r in a_report["per_structure"]}
    a_csv = load_csv(A / "per_atom_force_errors.csv")
    predictions, energies, origins, inference_times = {}, {}, {}, {}
    reuse_verification = {}
    for label in LABELS:
        a, r = atoms_by_label[label], a_by_label[label]
        assert r["source_sha256"] == EXPECTED_INPUTS[label] and r["result_sha256"] == EXPECTED_RESULTS[label]
        assert r["pbc"] == a.pbc.tolist()
        rows = r["per_atom"]
        assert len(rows) == len(a)
        assert [x["atom_id"] for x in rows] == a.arrays["lammps_id"].tolist()
        assert [x["species"] for x in rows] == a.get_chemical_symbols()
        assert np.array_equal([x["position_A"] for x in rows], a.positions)
        assert np.array_equal([x["dft_force_eV_A"] for x in rows], a.arrays["PW_PBE_forces"])
        assert np.array_equal([x["central_pair"] for x in rows], a.arrays["central_pair"].astype(bool))
        prediction = np.asarray([x["model_force_eV_A"] for x in rows], dtype=float)
        csv_rows = [x for x in a_csv if x["label"] == label]
        assert len(csv_rows) == len(a)
        for index, x in enumerate(csv_rows):
            assert int(x["atom_id"]) == int(a.arrays["lammps_id"][index])
            assert x["species"] == a[index].symbol
            assert np.allclose([float(x[k]) for k in ("x", "y", "z")], a.positions[index], rtol=0, atol=1e-8)
            for prefix, expected in (("dft", a.arrays["PW_PBE_forces"][index]), ("mace", prediction[index])):
                assert np.allclose([float(x[prefix + "_f" + axis]) for axis in "xyz"], expected, rtol=0, atol=1e-10)
        metric = errors(prediction - a.arrays["PW_PBE_forces"])
        assert abs(metric["force_vector_RMSE_eV_A"] - r["all_atom_force_vector_RMSE_eV_A"]) < 1e-12
        assert abs(metric["max_atom_force_vector_error_eV_A"] - r["max_atom_force_vector_error_eV_A"]) < 1e-12
        predictions[(label, V18)] = prediction
        energies[(label, V18)] = float(r["model_energy_eV_cell"])
        origins[(label, V18)] = "REUSED_A_FULL_PRECISION_METRICS_JSON"
        reuse_verification[label + ":" + V18] = {"status": "PASS", "model_hash": True, "geometry_order_ids_PBC": True, "DFT_vectors": True, "CSV_rounding_crosscheck": True, "metrics_recomputed": True}

    # Minus MP-0b3 was already inferred by C: exact archive/geometry/atom-ID identity only.
    minus = LABELS[0]
    old_inventory = next(r for r in load_csv(C / "reference_inventory.csv") if r["label"] == minus)
    assert old_inventory["source_sha256"] == EXPECTED_INPUTS[minus] and old_inventory["result_sha256"] == EXPECTED_RESULTS[minus]
    assert old_inventory["summary_sha256"] == archives[minus]["summary_sha256"]
    assert json.loads(old_inventory["pbc"]) == atoms_by_label[minus].pbc.tolist()
    assert old_inventory["parent_family_id"] == FAMILY
    old_model = next(r for r in load_csv(C / "model_inventory.csv") if r["model"] == MP)
    assert old_model["sha256"] == MODEL_SHAS[MP]
    c_rows = [r for r in load_csv(C / "per_atom_errors.csv") if r["label"] == minus and r["model"] == MP]
    a = atoms_by_label[minus]
    assert len(c_rows) == len(a)
    assert [int(r["index_zero_based"]) for r in c_rows] == list(range(len(a)))
    assert [int(r["atom_id"]) for r in c_rows] == a.arrays["lammps_id"].tolist()
    assert [r["element"] for r in c_rows] == a.get_chemical_symbols()
    assert np.array_equal([[float(r[k]) for k in ("x_A", "y_A", "z_A")] for r in c_rows], a.positions)
    assert np.array_equal([[float(r[f"DFT_F{x}_eV_A"]) for x in "xyz"] for r in c_rows], a.arrays["PW_PBE_forces"])
    predictions[(minus, MP)] = np.array([[float(r[f"model_F{x}_eV_A"]) for x in "xyz"] for r in c_rows])
    old_metric = next(r for r in load_csv(C / "structure_metrics.csv") if r["label"] == minus and r["model"] == MP)
    for key, value in errors(predictions[(minus, MP)] - a.arrays["PW_PBE_forces"]).items():
        assert abs(value - float(old_metric[key])) < 1e-12
    origins[(minus, MP)] = "REUSED_C_THREE_MODEL_FULL_PRECISION_FORCE_CSV"
    energies[(minus, MP)] = None
    reuse_verification[minus + ":" + MP] = {"status": "PASS", "model_hash": True, "source_result_summary_hashes": True, "geometry_order_ids_PBC": True, "DFT_vectors": True, "metrics_recomputed": True, "energy": "NOT_ARCHIVED; no duplicate inference"}

    # Only parent and plus require fresh inference. Frozen CPU double model; no optimizers.
    start_load = time.perf_counter()
    calculator = MACECalculator(model_paths=str(MODEL_PATHS[MP]), device="cpu", default_dtype="float64")
    load_seconds = time.perf_counter() - start_load
    model = calculator.models[0]
    model.eval()
    assert not model.training and all(not p.requires_grad for p in model.parameters())
    assert all(p.dtype == torch.float64 for p in model.parameters())
    assert float(model.r_max) == 6.0
    supported = model.atomic_numbers.detach().cpu().tolist()
    assert {6, 14, 22, 47}.issubset(supported)
    for label in LABELS[1:]:
        a = atoms_by_label[label].copy()
        a.calc = calculator
        started = time.perf_counter()
        predictions[(label, MP)] = np.asarray(a.get_forces(), dtype=float).copy()
        energies[(label, MP)] = float(a.get_potential_energy())  # Same cached calculator evaluation.
        inference_times[label] = time.perf_counter() - started
        assert predictions[(label, MP)].shape == (52, 3) and np.isfinite(predictions[(label, MP)]).all()
        assert np.isfinite(energies[(label, MP)])
        origins[(label, MP)] = "NEW_FROZEN_CPU_FLOAT64_INFERENCE_2_THREADS"
    assert all(not p.requires_grad for p in model.parameters()) and not model.training
    hashes_after = {m: sha(p) for m, p in MODEL_PATHS.items()}
    assert hashes_after == hashes_before

    atom_rows, element_rows, metric_rows, force_quantities = [], [], [], {}
    pair_definitions = {}
    parent = LABELS[1]
    for label in LABELS:
        a = atoms_by_label[label]
        symbols = np.array(a.get_chemical_symbols())
        ids = a.arrays["lammps_id"]
        ag_mask = symbols == "Ag"
        ag_indices = np.flatnonzero(ag_mask)
        ag_mean_z = float(a.positions[ag_mask, 2].mean())
        distances = np.array([float(a.get_distances(i, ag_indices, mic=True).min()) for i in range(len(a))])
        i, j, normal, pair_distance = pair_geometry(a)
        dft = np.asarray(a.arrays["PW_PBE_forces"], dtype=float)
        dft_gz = float(dft[ag_mask, 2].sum())
        dft_pair = float((dft[j] - dft[i]) @ normal)
        force_quantities[(label, "DFT")] = {"Gz": dft_gz, "pair": dft_pair}
        pair_definitions[label] = {"Ag_index_zero_based": i, "Si_index_zero_based": j, "Ag_atom_id": int(ids[i]), "Si_atom_id": int(ids[j]), "normal_Ag_to_Si_MIC": normal.tolist(), "distance_A": pair_distance}
        for name in (MP, V18):
            prediction = predictions[(label, name)]
            delta = prediction - dft
            metrics = errors(delta)
            model_gz = float(prediction[ag_mask, 2].sum())
            model_pair = float((prediction[j] - prediction[i]) @ normal)
            force_quantities[(label, name)] = {"Gz": model_gz, "pair": model_pair}
            parent_model_energy = energies[(parent, name)]
            energy = energies[(label, name)]
            relative_energy = energy - parent_model_energy if energy is not None else "NOT_ARCHIVED"
            native_relative = summaries[label]["energy_eV_cell"] - summaries[parent]["energy_eV_cell"]
            free_relative = summaries[label]["free_energy_eV_cell"] - summaries[parent]["free_energy_eV_cell"]
            metric_rows.append({"label": label, "ag_shift_z_A": SHIFTS[label], "parent_family_id": FAMILY, "dataset_role": records[label]["role"], "model": name, "model_sha256": MODEL_SHAS[name], "prediction_origin": origins[(label, name)], "source_sha256": EXPECTED_INPUTS[label], "result_sha256": EXPECTED_RESULTS[label], "n_atoms": len(a), **metrics, "marked_Ag_atom_id": int(ids[i]), "marked_Si_atom_id": int(ids[j]), "marked_pair_distance_A": pair_distance, "pair_DFT_projection_eV_A": dft_pair, "pair_model_projection_eV_A": model_pair, "pair_signed_error_eV_A": model_pair - dft_pair, "pair_half_coordinate_DFT_eV_A": dft_pair / 2, "pair_half_coordinate_model_eV_A": model_pair / 2, "whole_Ag_layer_DFT_Gz_eV_A": dft_gz, "whole_Ag_layer_model_Gz_eV_A": model_gz, "whole_Ag_layer_Gz_signed_error_eV_A": model_gz - dft_gz, "mean_Ag_DFT_Fz_eV_A": dft_gz / 4, "mean_Ag_model_Fz_eV_A": model_gz / 4, "DFT_native_delta_from_parent_eV_cell": native_relative, "DFT_free_energy_delta_from_parent_eV_cell": free_relative, "model_energy_delta_from_own_parent_eV_cell": relative_energy, "relative_energy_error_vs_native_eV_cell": relative_energy - native_relative if energy is not None else "NOT_ARCHIVED", "relative_energy_error_vs_free_eV_cell": relative_energy - free_relative if energy is not None else "NOT_ARCHIVED", "Ti_C_far_ge7A_squared_error_fraction": float(np.sum(delta[(distances >= 7) & np.isin(symbols, ["Ti", "C"])]**2) / np.sum(delta**2))})
            for element in sorted(set(symbols)):
                mask = symbols == element
                element_rows.append({"label": label, "ag_shift_z_A": SHIFTS[label], "model": name, "element": element, "n_atoms": int(mask.sum()), **errors(delta[mask])})
            for index in range(len(a)):
                atom_rows.append({"label": label, "ag_shift_z_A": SHIFTS[label], "parent_family_id": FAMILY, "dataset_role": records[label]["role"], "model": name, "prediction_origin": origins[(label, name)], "index_zero_based": index, "atom_id": int(ids[index]), "element": str(symbols[index]), "x_A": float(a.positions[index, 0]), "y_A": float(a.positions[index, 1]), "z_A": float(a.positions[index, 2]), "central_pair": bool(a.arrays["central_pair"][index]), "z_from_Ag_plane_A": float(a.positions[index, 2] - ag_mean_z), "nearest_Ag_distance_A": float(distances[index]), **{f"DFT_F{x}_eV_A": float(dft[index, k]) for k, x in enumerate("xyz")}, **{f"model_F{x}_eV_A": float(prediction[index, k]) for k, x in enumerate("xyz")}, **{f"error_F{x}_eV_A": float(delta[index, k]) for k, x in enumerate("xyz")}, "error_norm_eV_A": float(np.linalg.norm(delta[index]))})
            if name == V18:
                arow = a_by_label[label]
                assert abs(model_pair - arow["marked_pair"]["legacy_force_difference_definition_eV_A"]["model"]) < 1e-12
                assert abs(model_gz - arow["whole_Ag_layer_normal_force_eV_A"]["model"]) < 1e-12

    response_rows = []
    def add_response(name, quantity, unit, dft, predicted, definition, caveat, status="AVAILABLE"):
        error = predicted - dft if status == "AVAILABLE" else "NOT_ARCHIVED"
        response_rows.append({"model": name, "quantity": quantity, "minus_shift_A": -0.2, "plus_shift_A": 0.2, "span_A": 0.4, "unit": unit, "DFT_value": dft, "model_value": predicted if status == "AVAILABLE" else "NOT_ARCHIVED", "signed_model_minus_DFT": error, "relative_signed_error_percent": 100 * error / dft if status == "AVAILABLE" and dft != 0 else "NOT_ARCHIVED", "status": status, "definition": definition, "interpretation_limit": caveat})
    for name in (MP, V18):
        dft_layer = -(force_quantities[(LABELS[2], "DFT")]["Gz"] - force_quantities[(LABELS[0], "DFT")]["Gz"]) / 0.4
        pred_layer = -(force_quantities[(LABELS[2], name)]["Gz"] - force_quantities[(LABELS[0], name)]["Gz"]) / 0.4
        add_response(name, "rigid_Ag_layer_z_restoring_force_secant", "eV/A^2", dft_layer, pred_layer, "-[Gz(+0.20)-Gz(-0.20)]/0.40; Gz=sum of four Ag Fz", "Finite-interval layer-coordinate response; not bulk elasticity or an infinitesimal derivative")
        add_response(name, "mean_Ag_Fz_restoring_secant", "eV/A^2", dft_layer / 4, pred_layer / 4, "layer-coordinate force secant divided by four Ag atoms", "Per-Ag average of rigid-layer response")
        dft_pair = -(force_quantities[(LABELS[2], "DFT")]["pair"] - force_quantities[(LABELS[0], "DFT")]["pair"]) / 0.4
        pred_pair = -(force_quantities[(LABELS[2], name)]["pair"] - force_quantities[(LABELS[0], name)]["pair"]) / 0.4
        add_response(name, "marked_pair_projection_negative_secant_per_rigid_Ag_z", "eV/A^2", dft_pair, pred_pair, "-[P(+0.20)-P(-0.20)]/0.40; P=(F_Si-F_Ag).unit_MIC(Ag->Si)", "Normals vary by geometry; this local projection is neither layer force nor pair-distance stiffness")
        available = all(energies[(label, name)] is not None for label in LABELS)
        minus_e, parent_e, plus_e = [energies[(label, name)] for label in LABELS]
        pred_gz = -(plus_e - minus_e) / 0.4 if available else None
        pred_curvature = (plus_e + minus_e - 2 * parent_e) / 0.04 if available else None
        for key, summary_key in (("native", "energy_eV_cell"), ("free_energy", "free_energy_eV_cell")):
            dm, dc, dp = [summaries[label][summary_key] for label in LABELS]
            add_response(name, "central_energy_slope_Gz_vs_" + key, "eV/A", -(dp - dm) / 0.4, pred_gz, "-[E(+0.20)-E(-0.20)]/0.40 within each energy zero", "Free energy is the force-consistent DFT key; native energy retained separately; finite-interval secant differs from direct center force", "AVAILABLE" if available else "NOT_ARCHIVED")
            add_response(name, "three_point_energy_curvature_vs_" + key, "eV/A^2", (dp + dm - 2 * dc) / 0.04, pred_curvature, "[E(+0.20)+E(-0.20)-2E(0)]/0.20^2 within each energy zero", "Finite-interval energy curvature, not the endpoint force secant; no cross-model absolute energies compared", "AVAILABLE" if available else "NOT_ARCHIVED")

    assert len(atom_rows) == 312 and len(element_rows) == 24 and len(metric_rows) == 6 and len(response_rows) == 14
    write_csv("interface_per_atom.csv", atom_rows)
    write_csv("interface_per_element.csv", element_rows)
    write_csv("interface_metrics.csv", metric_rows)
    write_csv("interface_response_comparison.csv", response_rows)
    metadata = {"scope": "Only three public A AgSi rigid-layer archived geometries; frozen-model development/numerical diagnostics", "parent_family_id": FAMILY, "original_source_CIF_SHA256_recorded_by_A": manifest["source_CIF_sha256"], "original_CIF_bytes_not_reopened": True, "model_paths": {m: relative(p) for m, p in MODEL_PATHS.items()}, "model_sha256_before": hashes_before, "model_sha256_after": hashes_after, "model_weights_unchanged": True, "inference_counts": {"new_MP_configurations": 2, "reused_C_MP_configurations": 1, "reused_A_V18_configurations": 3, "new_V18_configurations": 0, "new_backbone_configurations": 0}, "new_inference_runtime": {"python": sys.version.split()[0], "torch": torch.__version__, **{p: importlib.metadata.version(p) for p in ("mace-torch", "ase", "numpy", "scipy")}, "device": "CPU", "dtype": "float64", "torch_threads": torch.get_num_threads(), "torch_interop_threads": torch.get_num_interop_threads(), "requires_grad_parameters": False, "network_eval_mode": True, "r_max_A": float(model.r_max), "atomic_numbers_order": supported, "model_load_seconds": load_seconds, "inference_seconds": inference_times}, "reused_V18_runtime": a_report["candidate"]["runtime"], "reused_C_MP_runtime": {"definition": old_model["inference_runtime"], "threads": 4}, "fresh_and_reused_force_dtype": "float64; full-precision archived predictions", "scientific_DFT_method": method, "energy_convention": manifest["energy_convention"], "public_source_dependencies_SHA256": public_dependencies, "pair_definition": "P=(F_Si-F_Ag).unit_MIC(Ag->Si); P/2 is generalized force for equal/opposite half-displacements", "layer_definition": "Gz=sum of four Ag Fz; mean Ag Fz=Gz/4; local marked-pair projection is not layer force", "pair_geometry_by_label": pair_definitions, "energy_comparison": "Only differences from the same model's parent; both DFT native and free_energy differences retained; reused MP minus energy NOT_ARCHIVED, with no duplicate inference", "metric_definition": {"vector_RMSE": "sqrt(mean_i(sum_xyz((F_model-F_DFT)^2)))", "axis_RMSE": "sqrt(mean_i((F_model,axis-F_DFT,axis)^2))", "max_atom_error": "max_i(norm_2(F_model-F_DFT))"}, "scientific_limitations": ["Same scored source parent family; no independent-test claim", "No new DFT, MD, training, or sealed structure/label reads", "Original numerical diagnostic roles preserved; none added to training", "No melt or solidification capability established by static forces", "Backbone evidence reused from prior published four-point audit; no new backbone inference"], "archives": archives}
    verification = {"status": "PASS", "checks": {"three_DFT_archive_verifiers_18_each": True, "input_result_hashes_pinned": True, "parent_geometry_hash_and_rigid_transform": True, "same_order_IDs_PBC_and_pair": True, "same_DFT_scientific_method_and_energy_convention": True, "four_archived_model_configuration_predictions_reused": True, "two_missing_MP_predictions_only": True, "V18_never_loaded_or_inferred": True, "model_hashes_unchanged": True, "full_precision_V18_JSON_matches_published_rounded_CSV": True, "six_vector_and_element_metrics_recomputed": True, "no_energy_zero_mixing": True, "original_data_roles_unchanged": True}, "reuse_validation": reuse_verification, "output_rows": {"interface_per_atom.csv": len(atom_rows), "interface_per_element.csv": len(element_rows), "interface_metrics.csv": len(metric_rows), "interface_response_comparison.csv": len(response_rows)}, "missing_quantities": ["MP minus energy, MP central energy slope, MP three-point energy curvature: NOT_ARCHIVED; no duplicate inference"]}
    (OUT / "inference_metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    (OUT / "verification.json").write_text(json.dumps(verification, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": "PASS", "fresh_inference_configurations": 2, "reused_prediction_configurations": 4, "new_inference_seconds": inference_times, "force_metrics": [{k: r[k] for k in ("ag_shift_z_A", "model", "force_vector_RMSE_eV_A", "max_atom_force_vector_error_eV_A", "pair_model_projection_eV_A", "whole_Ag_layer_model_Gz_eV_A")} for r in metric_rows]}, indent=2))


if __name__ == "__main__":
    main()
