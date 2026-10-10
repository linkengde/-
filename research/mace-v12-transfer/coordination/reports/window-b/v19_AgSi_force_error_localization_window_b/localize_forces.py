#!/usr/bin/env python3
"""Reproduce V18 force errors and localize them on two public AgSi CIF references."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mace-cpu-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/mace-cpu-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import torch
from ase.geometry import find_mic
from ase.io import read, write
from mace.calculators import MACECalculator


MODEL = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model")
CHECKPOINT = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/checkpoints/MACE_periodic_v18_clean_core_run-45_epoch-79.pt")
TRAIN = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/data/train.extxyz")
DATASET_MANIFEST = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/data/dataset_manifest.json")
SELECTION_RECORD = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/selection_record.json")
ARTIFACT_MANIFEST = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/artifact_manifest.json")
PREVIOUS = Path("research/mace-v12-transfer/coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/metrics.json")
REPORT = Path(__file__).resolve().parent
REFS = {
    "AgSi_COD9009647_pilot_k6x6": Path("research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot/calculations/AgSi_COD9009647_pilot_k6x6"),
    "AgSi_COD9009647_vacuum26_k6x6_sigma0p10": Path("research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_vacuum_pilot/calculations/AgSi_COD9009647_vacuum26_k6x6_sigma0p10"),
}
GATES = {"force_vector_RMSE_eV_A": 0.05, "separating_force_abs_error_eV_A": 0.10}
NEIGHBOR_CUTOFF_A = 4.5


def repo_root() -> Path:
    for p in (REPORT, *REPORT.parents):
        if (p / ".git").exists():
            return p
    raise RuntimeError("Repository root not found")


ROOT = repo_root()


def path(relative: Path) -> Path:
    p = (ROOT / relative).resolve()
    if ROOT.resolve() not in p.parents or not p.is_file():
        raise RuntimeError(f"Missing or out-of-repository input: {relative}")
    return p


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ids(atoms) -> np.ndarray:
    return np.asarray(atoms.arrays.get("lammps_id", np.arange(len(atoms))), dtype=np.int64)


def marked_pair(atoms) -> dict:
    mask = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    idx = np.flatnonzero(mask)
    ag = [int(i) for i in idx if atoms[int(i)].symbol == "Ag"]
    x = [int(i) for i in idx if atoms[int(i)].symbol != "Ag"]
    if len(idx) != 2 or len(ag) != 1 or len(x) != 1:
        raise RuntimeError("Public reference must contain one marked Ag-X pair")
    i, j = ag[0], x[0]
    vec, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    return {"Ag_index": i, "X_index": j, "Ag_id": int(ids(atoms)[i]), "X_id": int(ids(atoms)[j]), "X_species": atoms[j].symbol, "distance_A": float(np.linalg.norm(vec))}


def archive(label: str, folder_rel: Path, prior: dict) -> dict:
    folder = (ROOT / folder_rel).resolve()
    sums = folder / "SHA256SUMS.txt"
    checked = []
    for line in sums.read_text().splitlines():
        if not line.strip():
            continue
        expected, name = line.strip().split(maxsplit=1)
        name = name.lstrip("*")
        file = (folder / name).resolve()
        if folder not in file.parents or not file.is_file():
            raise RuntimeError(f"Invalid archive member {label}/{name}")
        actual = sha256(file)
        checked.append({"path": name, "sha256": actual, "pass": expected == actual})
    if not checked or not all(row["pass"] for row in checked):
        raise RuntimeError(f"Archive checksum failure: {label}")
    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Archive is not a converged match: {label}")
    checks = verification.get("checks", {})
    if verification.get("status") != "PASS" or not checks or not all(checks.values()):
        raise RuntimeError(f"Stored archive verification failed: {label}")
    result = folder / f"{label}_PW_PBE.extxyz"
    prior_entry = next(row for row in prior["per_reference_results"] if row["label"] == label)
    if sha256(result) != prior_entry["reference_archive_result_sha256"]:
        raise RuntimeError(f"Reference result differs from completed inference record: {label}")
    atoms = read(result, index=0)
    if len(atoms) != 52 or atoms.get_chemical_formula() != "C16Ag4Si8Ti24" or atoms.pbc.tolist() != [True, True, True]:
        raise RuntimeError(f"Reference identity/PBC check failed: {label}")
    forces = atoms.arrays.get("PW_PBE_forces")
    if forces is None or np.asarray(forces).shape != (len(atoms), 3) or not np.isfinite(forces).all():
        raise RuntimeError(f"Public DFT force array missing or invalid: {label}")
    atom_ids = ids(atoms)
    if len(np.unique(atom_ids)) != len(atoms):
        raise RuntimeError(f"Nonunique atom IDs: {label}")
    source_text = str(summary.get("source", "")).replace("\\", "/")
    marker = "research/mace-v12-transfer/"
    offset = source_text.find(marker)
    if offset < 0:
        raise RuntimeError(f"Reference source path cannot be resolved: {label}")
    source_rel = Path(source_text[offset:])
    source_path = path(source_rel)
    if sha256(source_path) != summary.get("source_sha256"):
        raise RuntimeError(f"Reference source hash mismatch: {label}")
    source_atoms = read(source_path, index=0)
    geometry_matches = (
        np.array_equal(atoms.numbers, source_atoms.numbers)
        and np.array_equal(ids(atoms), ids(source_atoms))
        and np.array_equal(atoms.positions, source_atoms.positions)
        and np.array_equal(atoms.cell.array, source_atoms.cell.array)
        and np.array_equal(atoms.pbc, source_atoms.pbc)
    )
    if not geometry_matches:
        raise RuntimeError(f"Archived result differs from exact source atom order/geometry: {label}")
    return {"label": label, "folder": folder_rel, "summary": summary, "verification": verification, "checksum_checks": checked, "result_path": result, "result_sha256": sha256(result), "atoms": atoms, "source_sha256": summary.get("source_sha256"), "source_path": source_rel}


def force_projection(forces: np.ndarray, atoms, pair: dict) -> float:
    i, j = pair["Ag_index"], pair["X_index"]
    vec, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    return float(np.dot(forces[j] - forces[i], vec / np.linalg.norm(vec)))


def xyz_csv_row(label, atoms, i, dft, model, delta, rank):
    pair_mask = np.asarray(atoms.arrays["central_pair"], dtype=bool)
    return {
        "reference": label,
        "atom_index_zero_based": i,
        "atom_id": int(ids(atoms)[i]),
        "species": atoms[i].symbol,
        "x_A": atoms.positions[i, 0], "y_A": atoms.positions[i, 1], "z_A": atoms.positions[i, 2],
        "dft_Fx_eV_A": dft[i, 0], "dft_Fy_eV_A": dft[i, 1], "dft_Fz_eV_A": dft[i, 2],
        "model_Fx_eV_A": model[i, 0], "model_Fy_eV_A": model[i, 1], "model_Fz_eV_A": model[i, 2],
        "error_model_minus_dft_Fx_eV_A": delta[i, 0], "error_model_minus_dft_Fy_eV_A": delta[i, 1], "error_model_minus_dft_Fz_eV_A": delta[i, 2],
        "force_vector_error_norm_eV_A": float(np.linalg.norm(delta[i])),
        "squared_vector_error_eV2_A2": float(np.dot(delta[i], delta[i])),
        "central_pair_atom": bool(pair_mask[i]),
        "rank_by_error_within_reference": rank,
    }


def neighbor_records(label: str, atoms, error_norms: np.ndarray, top_indices: list[int]) -> tuple[list[dict], list[dict]]:
    atom_ids = ids(atoms)
    pair_mask = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    rows, pairs = [], []
    symbols = atoms.get_chemical_symbols()
    for rank, center in enumerate(top_indices, 1):
        local = {center}
        neighbor_data = []
        all_distances = []
        for j in range(len(atoms)):
            if j == center:
                continue
            vec, _ = find_mic(atoms.positions[j] - atoms.positions[center], atoms.cell, pbc=atoms.pbc)
            distance = float(np.linalg.norm(vec))
            all_distances.append((j, distance, vec))
            if distance <= NEIGHBOR_CUTOFF_A:
                local.add(j)
                neighbor_data.append((j, distance, vec))
        counts = {s: sum(1 for j, _, _ in neighbor_data if symbols[j] == s) for s in sorted(set(symbols))}
        nearest = {}
        for species in sorted(set(symbols)):
            candidates = [(d, int(atom_ids[j])) for j, d, _ in all_distances if symbols[j] == species]
            if candidates:
                distance, atom_id = min(candidates)
                nearest[species] = {"distance_A": distance, "atom_id": atom_id}
        marked_indices = np.flatnonzero(pair_mask)
        marked_distances = {}
        for j in marked_indices:
            vec, _ = find_mic(atoms.positions[int(j)] - atoms.positions[center], atoms.cell, pbc=atoms.pbc)
            marked_distances["marked_Ag" if symbols[int(j)] == "Ag" else "marked_X"] = float(np.linalg.norm(vec))

        agx = []
        local_sorted = sorted(local)
        for pos, i in enumerate(local_sorted):
            for j in local_sorted[pos + 1:]:
                if (symbols[i] == "Ag") == (symbols[j] == "Ag"):
                    continue
                vec, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
                distance = float(np.linalg.norm(vec))
                if distance <= NEIGHBOR_CUTOFF_A:
                    agx.append({"Ag_id": int(atom_ids[i] if symbols[i] == "Ag" else atom_ids[j]), "X_id": int(atom_ids[j] if symbols[i] == "Ag" else atom_ids[i]), "X_species": symbols[j] if symbols[i] == "Ag" else symbols[i], "distance_A": distance, "marked_pair": bool(pair_mask[i] and pair_mask[j])})
        agx.sort(key=lambda r: r["distance_A"])
        rows.append({
            "reference": label, "center_rank": rank, "center_atom_index_zero_based": center,
            "center_atom_id": int(atom_ids[center]), "center_species": symbols[center],
            "center_force_vector_error_eV_A": float(error_norms[center]),
            "neighbor_cutoff_A": NEIGHBOR_CUTOFF_A,
            "neighbor_counts_by_species_json": json.dumps(counts, sort_keys=True),
            "nearest_neighbor_by_species_json": json.dumps(nearest, sort_keys=True),
            "local_Ag_X_pair_count_within_cutoff": len(agx),
            "closest_local_Ag_X_pairs_json": json.dumps(agx[:5], sort_keys=True),
            "distance_to_marked_Ag_A": marked_distances.get("marked_Ag"),
            "distance_to_marked_X_A": marked_distances.get("marked_X"),
            "marked_Ag_X_distance_A": marked_pair(atoms)["distance_A"],
            "center_is_marked_pair": bool(pair_mask[center]),
        })
        for j, distance, vec in sorted(neighbor_data, key=lambda item: item[1]):
            pairs.append({
                "reference": label, "center_rank": rank,
                "center_atom_id": int(atom_ids[center]), "center_species": symbols[center],
                "neighbor_atom_id": int(atom_ids[j]), "neighbor_species": symbols[j],
                "distance_A": distance, "dx_min_image_A": float(vec[0]), "dy_min_image_A": float(vec[1]), "dz_min_image_A": float(vec[2]),
                "center_marked": bool(pair_mask[center]), "neighbor_marked": bool(pair_mask[j]),
            })
    return rows, pairs


def write_csv(name: str, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {name}")
    with (REPORT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    prior_path = path(PREVIOUS)
    prior = json.loads(prior_path.read_text())
    expected_hashes = prior["candidate"]["model_and_data_lineage_hashes"]
    actual = {
        "model_sha256": sha256(path(MODEL)),
        "epoch79_checkpoint_sha256": sha256(path(CHECKPOINT)),
        "training_input_sha256": sha256(path(TRAIN)),
        "dataset_manifest_sha256": sha256(path(DATASET_MANIFEST)),
        "selection_record_sha256": sha256(path(SELECTION_RECORD)),
        "artifact_manifest_sha256": sha256(path(ARTIFACT_MANIFEST)),
    }
    for key, value in actual.items():
        if expected_hashes[key] != value:
            raise RuntimeError(f"Pinned V18 input hash changed: {key}")
    selection = json.loads(path(SELECTION_RECORD).read_text())
    dataset_manifest = json.loads(path(DATASET_MANIFEST).read_text())
    if selection.get("selected_epoch") != 79 or selection.get("blind_labels_used_for_selection") is not False or selection.get("test_file_passed_to_training") is not False:
        raise RuntimeError("V18 selection record no longer meets the published safeguards")
    if dataset_manifest.get("withheld_labels_opened") is not False or dataset_manifest.get("train_input_sha256") != actual["training_input_sha256"]:
        raise RuntimeError("V18 non-sealed training-manifest checks failed")
    if torch.cuda.is_available():
        raise RuntimeError("CPU-only task; unexpected CUDA device")
    torch.set_num_threads(4)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    refs = {label: archive(label, rel, prior) for label, rel in REFS.items()}
    train_atoms = read(path(TRAIN), index=":")
    if len(train_atoms) != 30:
        raise RuntimeError(f"Expected the published 30-frame V18 training set, got {len(train_atoms)}")
    train_pair_distances = []
    for frame_index, atoms in enumerate(train_atoms):
        marks = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
        pair_idx = np.flatnonzero(marks)
        if len(pair_idx) == 2:
            ag_idx = [int(i) for i in pair_idx if atoms[int(i)].symbol == "Ag"]
            x_idx = [int(i) for i in pair_idx if atoms[int(i)].symbol != "Ag"]
            if len(ag_idx) == len(x_idx) == 1 and atoms[x_idx[0]].symbol == "Si":
                vec, _ = find_mic(atoms.positions[x_idx[0]] - atoms.positions[ag_idx[0]], atoms.cell, pbc=atoms.pbc)
                train_pair_distances.append({"frame_index_zero_based": frame_index, "config_type": atoms.info.get("config_type"), "Ag_id": int(ids(atoms)[ag_idx[0]]), "Si_id": int(ids(atoms)[x_idx[0]]), "distance_A": float(np.linalg.norm(vec))})
    if not train_pair_distances:
        raise RuntimeError("No AgSi marked pairs found in allowed V18 train.extxyz")

    calculator = MACECalculator(model_paths=str(path(MODEL)), device="cpu", default_dtype="float64")
    all_atom_rows, species_rows, neighborhoods, neighbor_pairs = [], [], [], []
    summaries = []
    top_error_maps = {}
    original_atoms = {}
    for label, record in refs.items():
        atoms = record["atoms"]
        original_atoms[label] = atoms.copy()
        dft = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float).copy()
        atoms.calc = calculator
        model = np.asarray(atoms.get_forces(), dtype=float)
        atoms.calc = None
        if model.shape != dft.shape or not np.isfinite(model).all() or not np.isfinite(dft).all():
            raise RuntimeError(f"Invalid model/DFT forces: {label}")
        delta = model - dft
        norms = np.linalg.norm(delta, axis=1)
        ranking = np.argsort(norms)[::-1]
        rank_by_index = {int(i): rank + 1 for rank, i in enumerate(ranking)}
        for i in range(len(atoms)):
            all_atom_rows.append(xyz_csv_row(label, atoms, i, dft, model, delta, rank_by_index[i]))
        top = [int(i) for i in ranking[:10]]
        top_error_maps[label] = {int(ids(atoms)[i]): delta[i].copy() for i in top[:8]}
        nrows, prows = neighbor_records(label, atoms, norms, top)
        neighborhoods.extend(nrows)
        neighbor_pairs.extend(prows)

        total_sse = float(np.sum(norms**2))
        species_metric = {}
        for species in sorted(set(atoms.get_chemical_symbols())):
            mask = np.asarray(atoms.get_chemical_symbols()) == species
            species_norm = norms[mask]
            sse = float(np.sum(species_norm**2))
            comp = np.sqrt(np.mean(delta[mask] ** 2, axis=0))
            metric = {
                "atom_count": int(np.sum(mask)),
                "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(species_norm**2))),
                "max_atom_vector_error_eV_A": float(np.max(species_norm)),
                "mean_atom_vector_error_eV_A": float(np.mean(species_norm)),
                "component_RMSE_xyz_eV_A": comp.tolist(),
                "share_of_total_squared_error": sse / total_sse,
            }
            species_metric[species] = metric
            species_rows.append({"reference": label, "species": species, **metric, "component_RMSE_x_eV_A": float(comp[0]), "component_RMSE_y_eV_A": float(comp[1]), "component_RMSE_z_eV_A": float(comp[2])})
        pair = marked_pair(atoms)
        fsep_dft = force_projection(dft, atoms, pair)
        fsep_model = force_projection(model, atoms, pair)
        ag_mask = np.asarray(atoms.get_chemical_symbols()) == "Ag"
        agz_dft = float(np.sum(dft[ag_mask, 2]))
        agz_model = float(np.sum(model[ag_mask, 2]))
        force_rmse = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        prior_row = next(r for r in prior["per_reference_results"] if r["label"] == label)
        if abs(force_rmse - prior_row["force_metrics"]["force_vector_RMSE_eV_A"]) > 1e-6:
            raise RuntimeError(f"Reproduced force-vector RMSE differs from prior inference: {label}")
        if abs(abs(fsep_model - fsep_dft) - prior_row["force_metrics"]["marked_pair_separating_force_eV_A"]["absolute_error"]) > 1e-6:
            raise RuntimeError(f"Reproduced marked-pair force error differs from prior inference: {label}")
        summaries.append({
            "reference": label,
            "result_sha256": record["result_sha256"],
            "source_sha256": record["source_sha256"],
            "model_sha256": actual["model_sha256"],
            "checkpoint_sha256": actual["epoch79_checkpoint_sha256"],
            "training_input_sha256": actual["training_input_sha256"],
            "atom_count": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "all_atom_force_vector_RMSE_eV_A": force_rmse,
            "maximum_atom_vector_error_eV_A": float(np.max(norms)),
            "marked_pair": pair,
            "marked_pair_separation_force": {"definition": "dot(F_X-F_Ag, minimum_image_unit_vector_Ag_to_X)", "dft_signed_eV_A": fsep_dft, "model_signed_eV_A": fsep_model, "signed_error_model_minus_dft_eV_A": fsep_model-fsep_dft, "absolute_error_eV_A": abs(fsep_model-fsep_dft)},
            "Ag_layer_normal_force": {"definition": "sum of F_z over all Ag atoms; +z positive", "atom_count": int(np.sum(ag_mask)), "dft_net_eV_A": agz_dft, "model_net_eV_A": agz_model, "signed_error_model_minus_dft_eV_A": agz_model-agz_dft, "per_species": species_metric["Ag"]},
            "per_species": species_metric,
            "force_gate_pass": force_rmse <= GATES["force_vector_RMSE_eV_A"],
            "separation_gate_pass": abs(fsep_model-fsep_dft) <= GATES["separating_force_abs_error_eV_A"],
            "largest_error_atom_ids": [{"rank": rank+1, "atom_id": int(ids(atoms)[i]), "species": atoms[i].symbol, "error_eV_A": float(norms[i]), "residual_xyz_eV_A": delta[i].tolist()} for rank, i in enumerate(top)],
        })

    # Construct small signed geometry-only probes only for top-5 error sites that recur in both references.
    common_ids = set.intersection(*(set(m) for m in top_error_maps.values()))
    direction_support = []
    for atom_id in common_ids:
        vectors = [mapping[atom_id] for mapping in top_error_maps.values()]
        units = [v / np.linalg.norm(v) for v in vectors if np.linalg.norm(v) > 0]
        if len(units) == 2 and float(np.dot(units[0], units[1])) >= 0.90:
            direction_support.append((min(float(np.linalg.norm(v)) for v in vectors), atom_id, vectors))
    direction_support.sort(key=lambda row: (-row[0], row[1]))
    representative_directions = []
    for candidate in direction_support:
        vectors = candidate[2]
        units = [v / np.linalg.norm(v) for v in vectors if np.linalg.norm(v) > 0]
        direction = np.mean(units, axis=0)
        direction /= np.linalg.norm(direction)
        if all(float(np.dot(direction, previous)) < 0.95 for previous in representative_directions):
            representative_directions.append(direction)
            if len(representative_directions) == 2:
                break
    proposals = []
    proposal_atoms = []
    used_proposal_directions = []
    source_label = "AgSi_COD9009647_pilot_k6x6"
    base = original_atoms[source_label]
    for residual_norm, atom_id, vectors in direction_support:
        units = [v / np.linalg.norm(v) for v in vectors if np.linalg.norm(v) > 0]
        candidate_direction = np.mean(units, axis=0)
        candidate_direction /= np.linalg.norm(candidate_direction)
        if not any(float(np.dot(candidate_direction, direction)) >= 0.95 for direction in representative_directions):
            continue
        # Keep one representative per direction cluster; equivalent sites do not need duplicate probes.
        if any(float(np.dot(candidate_direction, used)) >= 0.95 for used in used_proposal_directions):
            continue
        used_proposal_directions.append(candidate_direction)
        matches = np.flatnonzero(ids(base) == atom_id)
        if len(matches) != 1:
            continue
        idx = int(matches[0])
        direction = np.mean([v / np.linalg.norm(v) for v in vectors], axis=0)
        direction /= np.linalg.norm(direction)
        for sign in (-1, 1):
            delta_pos = sign * 0.03 * direction
            candidate = base.copy()
            candidate.positions[idx] += delta_pos
            candidate.arrays.pop("PW_PBE_forces", None)
            for key in ("PW_PBE_energy_eV", "energy", "free_energy"):
                candidate.info.pop(key, None)
            candidate.calc = None
            proposal_id = f"disabled_probe_atom{atom_id}_{'m' if sign < 0 else 'p'}030"
            candidate.info["config_type"] = proposal_id
            candidate.info["proposal_id"] = proposal_id
            candidate.info["label_status"] = "UNLABELED_GEOMETRY_ONLY"
            candidate.info["launch_enabled"] = False
            candidate.info["training_eligible"] = False
            candidate.info["source_label"] = source_label
            candidate.info["changed_atom_id"] = int(atom_id)
            candidate.info["delta_xyz_A"] = delta_pos.tolist()
            proposal_atoms.append(candidate)
            proposals.append({
                "proposal_id": proposal_id,
                "source_label": source_label,
                "source_reference_sha256": refs[source_label]["result_sha256"],
                "changed_atom_id": int(atom_id),
                "changed_species": base[idx].symbol,
                "delta_xyz_A": delta_pos.tolist(),
                "step_length_A": float(np.linalg.norm(delta_pos)),
                "shared_top8_residual_norm_eV_A_minimum": residual_norm,
                "residual_direction_cosine_between_references": float(np.dot(vectors[0]/np.linalg.norm(vectors[0]), vectors[1]/np.linalg.norm(vectors[1]))),
                "launch_enabled": False,
                "owner": None,
                "training_eligible": False,
                "interpretation": "Unowned geometry-only local response probe; requires review/registration before any label acquisition.",
            })

    for summary in summaries:
        reference_distance = summary["marked_pair"]["distance_A"]
        train_pair_distance = min(train_pair_distances, key=lambda row: abs(row["distance_A"] - reference_distance))
        summary["training_AgSi_contact_coverage"] = {
            "training_frames_read": len(train_atoms),
            "local_pair_frame_count": len(train_pair_distances),
            "nearest_clean30_AgSi_pair_to_reference_distance": train_pair_distance,
            "reference_marked_pair_distance_A": reference_distance,
            "local_distance_overlap_within_0p05A": abs(train_pair_distance["distance_A"] - reference_distance) <= 0.05,
            "development_or_test_files_read": False,
            "sealed_labels_opened": False,
        }

    write_csv("per_atom_force_errors.csv", all_atom_rows)
    write_csv("species_force_breakdown.csv", species_rows)
    write_csv("top_error_neighborhoods.csv", neighborhoods)
    write_csv("top_error_neighbor_pairs.csv", neighbor_pairs)
    result = {
        "status": "LOCALIZATION_COMPLETE",
        "scope": "MACE CPU inference/error localization on only the two prior public AgSi archives; train.extxyz read for local contact coverage only; no DFT/MPI, training, sealed labels, dev/test structures, or model edits.",
        "runtime": {"python": sys.version.split()[0], "torch": torch.__version__, "mace_torch": "0.3.16", "device": "cpu", "torch_num_threads": torch.get_num_threads()},
        "input_hashes": {**actual, "previous_inference_metrics_sha256": sha256(prior_path)},
        "gate_definitions": GATES,
        "neighbor_cutoff_A": NEIGHBOR_CUTOFF_A,
        "reference_summaries": summaries,
        "training_contact_coverage": {"training_path": str(TRAIN), "training_sha256": actual["training_input_sha256"], "frames_read": len(train_atoms), "marked_AgSi_pairs": train_pair_distances, "withheld_labels_opened": False, "dev_or_test_structure_files_read": False},
        "geometry_suggestions": {"status": "UNOWNED_LAUNCH_DISABLED" if proposals else "NO_SUPPORTED_LOCAL_PROBE", "owner": None, "launch_enabled": False, "training_eligible": False, "proposal_count": len(proposals), "proposals": proposals, "extxyz_path": "disabled_geometry_proposals.extxyz" if proposals else None, "note": "Suggestions are report artifacts only; they do not alter any registered DFT or gap-scan entry."},
    }
    (REPORT / "localization.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    (REPORT / "disabled_geometry_proposals.json").write_text(json.dumps(result["geometry_suggestions"], indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    if proposal_atoms:
        write(REPORT / "disabled_geometry_proposals.extxyz", proposal_atoms, format="extxyz")
    print(json.dumps({"status": result["status"], "references": [{"label": r["reference"], "force_RMSE": r["all_atom_force_vector_RMSE_eV_A"], "max_error": r["maximum_atom_vector_error_eV_A"], "top_atoms": r["largest_error_atom_ids"][:5], "species": r["per_species"]} for r in summaries], "geometry_suggestions": result["geometry_suggestions"]}, indent=2))


if __name__ == "__main__":
    main()
