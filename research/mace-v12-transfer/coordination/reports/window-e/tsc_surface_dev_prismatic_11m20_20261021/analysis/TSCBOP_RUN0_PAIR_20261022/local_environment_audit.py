#!/usr/bin/env python3
"""Offline local-environment analysis of the largest frozen-force errors."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CASES = [
    ("REF", "TSC120_DEV_PRISM_11m20_COD9009647_10L_REF", "attempt-0003"),
    ("RESP_TIZ020", "TSC120_DEV_PRISM_11m20_COD9009647_10L_RESP_TIZ020", "attempt-0001"),
]
TARGET_IDS = [29, 30, 46]
CUTOFFS = [2.8, 3.0, 3.3, 3.5]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_force_table(path):
    with open(path, newline="") as f:
        records = list(csv.DictReader(f))
    result = {}
    for row in records:
        atom_id = int(row["source_atom_id"])
        result[atom_id] = {
            "element": row["element"],
            "dft": np.array([float(row[f"dft_f{a}_eV_A"]) for a in "xyz"]),
            "tscbop": np.array([float(row[f"tscbop_f{a}_eV_A"]) for a in "xyz"]),
            "error": np.array([float(row[f"error_f{a}_eV_A"]) for a in "xyz"]),
            "error_norm": float(row["force_error_norm_eV_A"]),
        }
    return result


def main():
    high_rows = []
    atom_rows = []
    case_summaries = {}
    for label, structure_id, attempt in CASES:
        input_path = ROOT / "inputs" / f"{structure_id}.extxyz"
        result_dir = ROOT / "results" / structure_id / attempt
        result = json.loads((result_dir / "RESULT.json").read_text())
        status = json.loads((result_dir / "RUN_STATUS.json").read_text())
        atoms = read(input_path)
        labelled = read(result_dir / "labelled.extxyz")
        comp = load_force_table(OUT / f"{label}_FORCE_COMPARISON.csv")
        ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
        syms = np.asarray(atoms.get_chemical_symbols())
        if status.get("status") != "SUCCESS" or result.get("scf_converged") is False:
            raise ValueError(f"{label}: source DFT did not converge")
        if not np.array_equal(atoms.numbers, labelled.numbers) or not np.array_equal(atoms.positions, labelled.positions):
            raise ValueError(f"{label}: labeled structure geometry differs from input")
        if set(comp) != set(ids.tolist()):
            raise ValueError(f"{label}: force table source IDs differ from input")
        for i, atom_id in enumerate(ids):
            if comp[int(atom_id)]["element"] != syms[i]:
                raise ValueError(f"{label}: element mismatch for source ID {atom_id}")
            if not np.allclose(comp[int(atom_id)]["dft"], result["forces_eV_per_angstrom"][i], atol=1e-12, rtol=0):
                raise ValueError(f"{label}: CSV force differs from RESULT.json at ID {atom_id}")

        z = atoms.positions[:, 2]
        zmin, zmax = float(z.min()), float(z.max())
        surface_distance = np.minimum(z - zmin, zmax - z)
        nearest = np.full(len(atoms), np.inf)
        coordination = {cutoff: np.zeros(len(atoms), dtype=int) for cutoff in CUTOFFS}
        counts = {cutoff: {el: np.zeros(len(atoms), dtype=int) for el in ("Ti", "Si", "C")} for cutoff in CUTOFFS}
        neighbor_lists = {}
        for i, atom_id in enumerate(ids):
            d = atoms.get_distances(i, np.arange(len(atoms)), mic=True)
            d[i] = np.inf
            nearest[i] = float(d.min())
            js = np.flatnonzero(d <= max(CUTOFFS))
            neighbor_lists[int(atom_id)] = sorted(
                [(float(d[j]), int(ids[j]), str(syms[j]), atoms.get_distance(i, j, mic=True, vector=True).tolist()) for j in js],
                key=lambda x: (x[0], x[1]),
            )
            for cutoff in CUTOFFS:
                js_cut = np.flatnonzero(d <= cutoff)
                coordination[cutoff][i] = len(js_cut)
                for j in js_cut:
                    counts[cutoff][str(syms[j])][i] += 1

        err = np.array([comp[int(i)]["error"] for i in ids])
        err_norm = np.linalg.norm(err, axis=1)
        surface_mask = surface_distance <= 3.3
        interior_mask = ~surface_mask
        case_summaries[label] = {
            "structure_id": structure_id,
            "input_sha256": sha(input_path),
            "result_json_sha256": sha(result_dir / "RESULT.json"),
            "atom_count": len(atoms),
            "z_atomic_extent_A": [zmin, zmax],
            "surface_band_definition": "distance in Cartesian z to the minimum or maximum atomic z <= 3.3 A; geometric heuristic, not a Tersoff cutoff",
            "surface_atom_count": int(surface_mask.sum()),
            "interior_atom_count": int(interior_mask.sum()),
            "force_error_vector_rmse_surface_eV_A": float(np.sqrt(np.mean(np.sum(err[surface_mask]**2, axis=1)))) if surface_mask.any() else None,
            "force_error_vector_rmse_interior_eV_A": float(np.sqrt(np.mean(np.sum(err[interior_mask]**2, axis=1)))) if interior_mask.any() else None,
            "pearson_error_norm_vs_coordination3p3": float(np.corrcoef(err_norm, coordination[3.3])[0, 1]) if np.std(coordination[3.3]) > 0 else None,
            "pearson_error_norm_vs_surface_distance": float(np.corrcoef(err_norm, surface_distance)[0, 1]) if np.std(surface_distance) > 0 else None,
            "mean_force_error_components_eV_A": {axis: float(err[:, k].mean()) for k, axis in enumerate("xyz")},
            "mean_force_error_components_by_element_eV_A": {
                element: {axis: float(err[syms == element, k].mean()) for k, axis in enumerate("xyz")}
                for element in ("Ti", "Si", "C")
            },
            "force_error_vector_rmse_by_coordination3p3": {
                str(int(cn)): float(np.sqrt(np.mean(np.sum(err[coordination[3.3] == cn] ** 2, axis=1))))
                for cn in sorted(set(coordination[3.3].tolist()))
                if int(np.sum(coordination[3.3] == cn)) >= 2
            },
            "target_atoms": {},
        }
        for i, atom_id in enumerate(ids):
            entry = {
                "geometry": label,
                "source_atom_id": int(atom_id),
                "element": str(syms[i]),
                "x_A": float(atoms.positions[i, 0]),
                "y_A": float(atoms.positions[i, 1]),
                "z_A": float(atoms.positions[i, 2]),
                "distance_to_nearest_atomic_z_surface_A": float(surface_distance[i]),
                "surface_band_le_3p3A": bool(surface_mask[i]),
                "nearest_neighbor_distance_A": float(nearest[i]),
                "coordination_le_3p3A": int(coordination[3.3][i]),
                "neighbors_Ti_le_3p3A": int(counts[3.3]["Ti"][i]),
                "neighbors_Si_le_3p3A": int(counts[3.3]["Si"][i]),
                "neighbors_C_le_3p3A": int(counts[3.3]["C"][i]),
                "dft_force_norm_eV_A": float(np.linalg.norm(comp[int(atom_id)]["dft"])),
                "tscbop_force_norm_eV_A": float(np.linalg.norm(comp[int(atom_id)]["tscbop"])),
                "force_error_norm_eV_A": float(err_norm[i]),
            }
            for cutoff in CUTOFFS:
                entry[f"coordination_le_{cutoff:.1f}A"] = int(coordination[cutoff][i])
            atom_rows.append(entry)
        for atom_id in TARGET_IDS:
            i = int(np.flatnonzero(ids == atom_id)[0])
            rec = comp[atom_id]
            coord = {f"cn_{cutoff:.1f}A": int(coordination[cutoff][i]) for cutoff in CUTOFFS}
            elem_coord = {f"{el}_neighbors_le_3p3A": int(counts[3.3][el][i]) for el in ("Ti", "Si", "C")}
            case_summaries[label]["target_atoms"][str(atom_id)] = {
                "element": str(syms[i]),
                "position_A": atoms.positions[i].tolist(),
                "distance_to_nearest_atomic_z_surface_A": float(surface_distance[i]),
                "nearest_neighbor_distance_A": float(nearest[i]),
                "coordination_counts": coord,
                "neighbor_species_counts_3p3A": elem_coord,
                "dft_force_eV_A": rec["dft"].tolist(),
                "tscbop_force_eV_A": rec["tscbop"].tolist(),
                "force_error_eV_A": rec["error"].tolist(),
                "force_error_norm_eV_A": rec["error_norm"],
                "neighbors_within_3p5A": [
                    {"source_atom_id": n_id, "element": n_el, "distance_A": dist, "minimum_image_vector_A": vec}
                    for dist, n_id, n_el, vec in neighbor_lists[atom_id]
                ],
            }
            for dist, n_id, n_el, vec in neighbor_lists[atom_id]:
                high_rows.append({
                    "geometry": label, "target_source_atom_id": atom_id, "target_element": str(syms[i]),
                    "target_surface_distance_A": float(surface_distance[i]), "target_force_error_norm_eV_A": rec["error_norm"],
                    "neighbor_source_atom_id": n_id, "neighbor_element": n_el, "distance_A": dist,
                    "minimum_image_dx_A": vec[0], "minimum_image_dy_A": vec[1], "minimum_image_dz_A": vec[2],
                    "within_2p8A": dist <= 2.8, "within_3p0A": dist <= 3.0, "within_3p3A": dist <= 3.3,
                })
    with open(OUT / "ATOM_ENVIRONMENT_SUMMARY.csv", "w", newline="") as f:
        cols = list(atom_rows[0])
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader(); w.writerows(atom_rows)
    with open(OUT / "HIGH_ERROR_ATOM_NEIGHBORS.csv", "w", newline="") as f:
        cols = list(high_rows[0])
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader(); w.writerows(high_rows)
    summary = {
        "scope": "Offline geometry and force-label analysis only; no new DFT, LAMMPS, MD, or fitting.",
        "cutoffs_A": CUTOFFS,
        "cutoff_note": "Coordination radii are geometric diagnostics only; they are not asserted to equal Tersoff switching cutoffs.",
        "target_source_atom_ids": TARGET_IDS,
        "potential_sha256": "729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7",
        "cases": case_summaries,
    }
    (OUT / "LOCAL_ENVIRONMENT_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
