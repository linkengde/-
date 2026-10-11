#!/usr/bin/env python3
"""Score the pinned V18 model on the archive-verified AgSi parent and +/-0.20 A slab scan."""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mace-v18-gap-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/mace-v18-gap-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import torch
from ase.geometry import find_mic
from ase.io import read
from mace.calculators import MACECalculator


ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[3]
LABELS = [
    "AgSi_rigid_Ag_z_m020_k6x6x1_sigma0p10_dipoleXY_v19",
    "AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19",
    "AgSi_rigid_Ag_z_p020_k6x6x1_sigma0p10_dipoleXY_v19",
]
MODEL = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
LINEAGE_SCRIPT = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/infer_and_evaluate.py"
VERIFY_SCRIPT = ROOT / "verify_result.py"
OUT = Path(__file__).resolve().parent
GATES = {"all_atom_force_vector_RMSE_eV_A": 0.05, "marked_pair_separation_force_abs_error_eV_A": 0.10}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_lineage_module():
    spec = importlib.util.spec_from_file_location("v18_reference_screen", LINEAGE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def score(label: str, calculator):
    folder = ROOT / "calculations" / label
    verification = json.loads(subprocess.check_output(
        [sys.executable, str(VERIFY_SCRIPT), label, str(folder)], text=True
    ))
    if verification.get("status") != "PASS" or not all(verification["checks"].values()):
        raise RuntimeError(f"DFT archive verification failed for {label}")
    result_path = folder / f"{label}_PW_PBE.extxyz"
    atoms = read(result_path)
    dft = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float).copy()
    atoms.calc = calculator
    t0 = time.perf_counter()
    model_energy = float(atoms.get_potential_energy())
    model = np.asarray(atoms.get_forces(), dtype=float)
    inference_s = time.perf_counter() - t0
    if model.shape != dft.shape or not np.isfinite(model).all() or not np.isfinite(model_energy):
        raise RuntimeError(f"Non-finite or malformed MACE prediction for {label}")

    delta = model - dft
    atom_error = np.linalg.norm(delta, axis=1)
    symbols = np.asarray(atoms.get_chemical_symbols())
    by_species = {}
    for element in sorted(set(symbols)):
        mask = symbols == element
        by_species[element] = {
            "n": int(mask.sum()),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta[mask] ** 2, axis=1)))),
            "max_atom_vector_error_eV_A": float(atom_error[mask].max()),
            "mean_atom_vector_error_eV_A": float(atom_error[mask].mean()),
        }

    pair_indices = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
    if len(pair_indices) != 2:
        raise RuntimeError(f"Expected two central_pair atoms for {label}")
    ag = [int(i) for i in pair_indices if atoms[int(i)].symbol == "Ag"]
    contact = [int(i) for i in pair_indices if atoms[int(i)].symbol != "Ag"]
    if len(ag) != 1 or len(contact) != 1:
        raise RuntimeError(f"Expected one marked Ag and one marked contact atom for {label}")
    i, j = ag[0], contact[0]
    vector, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    unit = vector / np.linalg.norm(vector)
    pair_model = float(np.dot(model[j] - model[i], unit))
    pair_dft = float(np.dot(dft[j] - dft[i], unit))
    pair_error = pair_model - pair_dft
    ag_mask = symbols == "Ag"
    ag_fz_model = float(model[ag_mask, 2].sum())
    ag_fz_dft = float(dft[ag_mask, 2].sum())
    summary = json.loads((folder / "summary.json").read_text())
    ag_indices = np.flatnonzero(ag_mask)
    ag_plane_z = float(np.mean(atoms.positions[ag_indices, 2]))
    nearest_ag_distance = np.array([
        float(np.min(atoms.get_distances(k, ag_indices, mic=True))) for k in range(len(atoms))
    ])
    z_from_ag_plane = atoms.positions[:, 2] - ag_plane_z
    distance_bins = [(0.0, 3.0, "0-3"), (3.0, 5.0, "3-5"), (5.0, 7.0, "5-7"), (7.0, float("inf"), ">7")]
    localization = {}
    for low, high, name in distance_bins:
        mask = (nearest_ag_distance >= low) & (nearest_ag_distance < high)
        if not mask.any():
            continue
        localization[name] = {
            "atom_count": int(mask.sum()),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta[mask] ** 2, axis=1)))),
            "z_component_RMSE_eV_A": float(np.sqrt(np.mean(delta[mask, 2] ** 2))),
            "species": {},
        }
        for element in sorted(set(symbols[mask])):
            submask = mask & (symbols == element)
            localization[name]["species"][str(element)] = {
                "n": int(submask.sum()),
                "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta[submask] ** 2, axis=1)))),
                "z_component_RMSE_eV_A": float(np.sqrt(np.mean(delta[submask, 2] ** 2))),
            }
    row = {
        "label": label,
        "dft_archive_verification": "PASS",
        "result_sha256": sha(result_path),
        "source_sha256": summary["source_sha256"],
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "pbc": atoms.pbc.tolist(),
        "scf_iterations": int(summary["scf_iterations"]),
        "dft_native_energy_eV_cell": float(summary["energy_eV_cell"]),
        "dft_free_energy_eV_cell": float(summary["free_energy_eV_cell"]),
        "model_energy_eV_cell": model_energy,
        "energy_model_minus_dft_eV_cell": {
            "native": model_energy - float(summary["energy_eV_cell"]),
            "free_energy": model_energy - float(summary["free_energy_eV_cell"]),
        },
        "inference_seconds": inference_s,
        "all_atom_force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1)))),
        "max_atom_force_vector_error_eV_A": float(atom_error.max()),
        "per_species_force_errors": by_species,
        "force_error_localization_by_nearest_Ag_distance": {
            "distance_definition": "minimum-image 3D distance to any Ag atom; z is non-periodic",
            "bins_A": "[0,3), [3,5), [5,7), [7,infinity)",
            "Ag_mean_z_A": ag_plane_z,
            "bins": localization,
        },
        "marked_pair": {
            "Ag_atom_id": int(atoms.arrays["lammps_id"][i]),
            "contact_atom_id": int(atoms.arrays["lammps_id"][j]),
            "contact_species": atoms[j].symbol,
            "distance_A": float(np.linalg.norm(vector)),
            "legacy_force_difference_definition_eV_A": {
                "definition": "(F_contact - F_Ag) dot unit(Ag to contact)",
                "model": pair_model,
                "dft": pair_dft,
                "signed_model_minus_dft": pair_error,
                "absolute_error": abs(pair_error),
                "gate_eV_A": GATES["marked_pair_separation_force_abs_error_eV_A"],
                "pass": abs(pair_error) <= GATES["marked_pair_separation_force_abs_error_eV_A"],
            },
            "symmetric_relative_coordinate_generalized_force_eV_A": {
                "definition": "0.5 * legacy pair force difference; equal/opposite half-displacement convention",
                "model": pair_model / 2,
                "dft": pair_dft / 2,
                "signed_model_minus_dft": pair_error / 2,
            },
        },
        "whole_Ag_layer_normal_force_eV_A": {
            "definition": "sum of F_z over all four Ag atoms; generalized force for rigid +z Ag-layer translation",
            "model": ag_fz_model,
            "dft": ag_fz_dft,
            "signed_model_minus_dft": ag_fz_model - ag_fz_dft,
        },
    }
    row["all_atom_gate"] = {
        "limit_eV_A": GATES["all_atom_force_vector_RMSE_eV_A"],
        "pass": row["all_atom_force_vector_RMSE_eV_A"] <= GATES["all_atom_force_vector_RMSE_eV_A"],
    }
    row["per_atom"] = [
        {
            "atom_id": int(atom_id), "species": str(symbol),
            "dft_force_eV_A": dft[k].tolist(), "model_force_eV_A": model[k].tolist(),
            "error_vector_eV_A": delta[k].tolist(), "error_norm_eV_A": float(atom_error[k]),
            "position_A": atoms.positions[k].tolist(),
            "z_from_Ag_plane_A": float(z_from_ag_plane[k]),
            "nearest_Ag_distance_A": float(nearest_ag_distance[k]),
            "central_pair": bool(atoms.arrays["central_pair"][k]),
        }
        for k, (atom_id, symbol) in enumerate(zip(atoms.arrays["lammps_id"], symbols))
    ]
    atoms.calc = None
    return row


def main():
    torch.set_num_threads(4)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass
    lineage = import_lineage_module().verify_model_lineage()[0]
    calculator = MACECalculator(model_paths=str(MODEL), device="cpu", default_dtype="float64")
    rows = [score(label, calculator) for label in LABELS]
    by_label = {row["label"]: row for row in rows}
    minus = by_label[LABELS[0]]
    parent = by_label[LABELS[1]]
    plus = by_label[LABELS[2]]
    delta_A = 0.4
    finite_difference = {}
    for key, field in (("dft_native", "dft_native_energy_eV_cell"), ("dft_free_energy", "dft_free_energy_eV_cell"), ("mace", "model_energy_eV_cell")):
        force = -(plus[field] - minus[field]) / delta_A
        finite_difference[key + "_Ag_layer_generalized_force_eV_A"] = force
    finite_difference["dft_force_consistent_quantity"] = "free_energy; native extrapolated energy reported separately"
    finite_difference["mace_minus_dft_free_energy_force_eV_A"] = (
        finite_difference["mace_Ag_layer_generalized_force_eV_A"]
        - finite_difference["dft_free_energy_Ag_layer_generalized_force_eV_A"]
    )
    results = {
        "title": "V18 MACE force screen on the fixed-geometry AgSi parent and rigid Ag-layer +/-0.20 A gap scan",
        "scope": "Same uploaded-CIF family diagnostic only; not independent validation, not a production qualification, and no training or sealed labels read.",
        "candidate": {
            "name": "V18 clean-core epoch 79",
            "model_path": str(MODEL.relative_to(REPO)),
            "model_sha256": sha(MODEL),
            "lineage_checks": lineage["checks"],
            "lineage_hashes": lineage["hashes"],
            "runtime": {
                "python": sys.version.split()[0],
                "torch": importlib.metadata.version("torch"),
                "mace-torch": importlib.metadata.version("mace-torch"),
                "ase": importlib.metadata.version("ase"),
                "numpy": importlib.metadata.version("numpy"),
                "device": "CPU",
                "torch_num_threads": torch.get_num_threads(),
            },
            "inference_device": "CPU",
        },
        "methods": {
            "dft": "PBE plane-wave 500 eV, 6x6x1, sigma 0.10 eV, pbc TTF, xy dipole, GPAW 26.7.0, four MPI ranks; fixed ion geometries",
            "all_atom_force_vector_RMSE": "sqrt(mean_i(sum_xyz((F_MACE-F_DFT)^2)))",
            "legacy_pair_separation_force": "(F_contact-F_Ag) dot unit(Ag to contact); existing provisional gate 0.10 eV/A",
            "whole_layer_generalized_force": "sum(F_Ag,z); direct derivative of rigid Ag-layer +z coordinate",
            "finite_difference": "minus central energy slope from +/-0.20 A endpoints; free_energy is force-consistent, native energy retained separately",
            "gates": GATES,
        },
        "energy_scan_deltas_from_parent_eV_cell": {
            key: {
                "minus_0p20": minus[key] - parent[key],
                "plus_0p20": plus[key] - parent[key],
                "plus_minus": plus[key] - minus[key],
            }
            for key in ("dft_native_energy_eV_cell", "dft_free_energy_eV_cell", "model_energy_eV_cell")
        },
        "central_difference_generalized_force_eV_A": finite_difference,
        "per_structure": rows,
        "gate_summary": {
            row["label"]: {
                "all_atom_force_vector_gate_pass": row["all_atom_gate"]["pass"],
                "marked_pair_separation_gate_pass": row["marked_pair"]["legacy_force_difference_definition_eV_A"]["pass"],
            }
            for row in rows
        },
        "interpretation": "Failures quantify this related-CIF force diagnostic. Do not call it an independent holdout; the common geometry family and V18's known local Ag-Si training overlap limit generalization claims.",
    }
    (OUT / "metrics.json").write_text(json.dumps(results, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    with (OUT / "per_atom_force_errors.csv").open("w") as f:
        f.write("label,atom_id,species,x,y,z,dz_from_Ag_plane,nearest_Ag_distance,dft_fx,dft_fy,dft_fz,mace_fx,mace_fy,mace_fz,error_fx,error_fy,error_fz,error_norm,central_pair\n")
        for row in rows:
            for atom in row["per_atom"]:
                vals = [*atom["position_A"], atom["z_from_Ag_plane_A"], atom["nearest_Ag_distance_A"], *atom["dft_force_eV_A"], *atom["model_force_eV_A"], *atom["error_vector_eV_A"]]
                f.write(f"{row['label']},{atom['atom_id']},{atom['species']}," + ",".join(f"{x:.12g}" for x in vals) + f",{atom['error_norm_eV_A']:.12g},{int(atom['central_pair'])}\n")
    print(json.dumps({
        "status": "complete",
        "rows": [{"label": r["label"], "force_RMSE_eV_A": r["all_atom_force_vector_RMSE_eV_A"], "force_pass": r["all_atom_gate"]["pass"], "pair_abs_error_eV_A": r["marked_pair"]["legacy_force_difference_definition_eV_A"]["absolute_error"], "pair_pass": r["marked_pair"]["legacy_force_difference_definition_eV_A"]["pass"]} for r in rows],
        "central_difference": finite_difference,
        "output": str(OUT),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
