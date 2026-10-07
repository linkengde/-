#!/usr/bin/env python3
"""Localize frozen-V16 force residuals on the three authorized V17 dev labels."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read


LABELS = (
    "AgC_v17_dev_validation_01",
    "AgSi_v17_dev_validation_01",
    "AgTi_v17_dev_validation_01",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", type=Path, required=True)
    args = ap.parse_args()
    repo = args.repo_root.resolve()
    out = Path(__file__).resolve().parent
    base = repo / "research/mace-v12-transfer/periodic_interface_v4"
    archive_root = base / "pbe_interface_v17_main_targeted_acquisition/development_validation_archives"
    atom_csv = base / "pbe_interface_v17_main_targeted_acquisition/development_validation_evaluation_v16/per_atom.csv"
    evaluation_path = base / "pbe_interface_v17_main_targeted_acquisition/development_validation_evaluation_v16/evaluation.json"
    evaluation = json.loads(evaluation_path.read_text())
    if evaluation["model_sha256"] != "918aad24341c6892d65ec0d31b19baf610ac0efde1cfd9d2e1c0567320e0cf98":
        raise SystemExit("Unexpected frozen V16 model hash")
    rows_by_label = {}
    with atom_csv.open(newline="") as f:
        for row in csv.DictReader(f):
            rows_by_label.setdefault(row["label"], []).append(row)

    results = []
    for label in LABELS:
        folder = archive_root / label
        structure = folder / f"{label}_PW_PBE.extxyz"
        atoms = read(structure)
        rows = rows_by_label[label]
        require_ids = set(int(x) for x in atoms.arrays["lammps_id"])
        if len(require_ids) != len(atoms) or {int(r["atom_id"]) for r in rows} != require_ids:
            raise SystemExit(f"Atom ID mapping mismatch: {label}")
        idx_by_id = {int(atom_id): i for i, atom_id in enumerate(atoms.arrays["lammps_id"])}
        central = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
        if len(central) != 2:
            raise SystemExit(f"Expected exactly two marked atoms: {label}")
        ag = [int(i) for i in central if atoms[int(i)].symbol == "Ag"]
        contact = [int(i) for i in central if atoms[int(i)].symbol != "Ag"]
        if len(ag) != 1 or len(contact) != 1:
            raise SystemExit(f"Marked pair must be Ag-X: {label}")
        ia, ix = ag[0], contact[0]
        axis = atoms.positions[ix] - atoms.positions[ia]
        axis, _ = find_mic(axis, atoms.cell, atoms.pbc)
        distance = float(np.linalg.norm(axis))
        axis /= distance

        atom_results = []
        local_sq, outer_sq = [], []
        for row in rows:
            atom_id = int(row["atom_id"])
            i = idx_by_id[atom_id]
            error = np.array([float(row[k]) for k in ("error_fx_eV_A", "error_fy_eV_A", "error_fz_eV_A")], dtype=float)
            error_norm = float(np.linalg.norm(error))
            parallel = float(np.dot(error, axis))
            transverse = float(np.linalg.norm(error - parallel * axis))
            da, _ = find_mic(atoms.positions[i] - atoms.positions[ia], atoms.cell, atoms.pbc)
            dx, _ = find_mic(atoms.positions[i] - atoms.positions[ix], atoms.cell, atoms.pbc)
            shell_distance = float(min(np.linalg.norm(da), np.linalg.norm(dx)))
            (local_sq if shell_distance <= 4.5 else outer_sq).append(error_norm ** 2)
            atom_results.append({
                "atom_id": atom_id,
                "symbol": atoms[i].symbol,
                "error_vector_norm_eV_A": error_norm,
                "error_parallel_to_marked_pair_eV_A": parallel,
                "error_transverse_to_marked_pair_eV_A": transverse,
                "distance_to_marked_pair_shell_A": shell_distance,
            })
        atom_results.sort(key=lambda x: x["error_vector_norm_eV_A"], reverse=True)
        max_atom = atom_results[0]
        max_i = idx_by_id[max_atom["atom_id"]]
        delta = atoms.positions - atoms.positions[max_i]
        delta, _ = find_mic(delta, atoms.cell, atoms.pbc)
        neighbor_distances = np.linalg.norm(delta, axis=1)
        neighbors = sorted(
            ({"atom_id": int(atoms.arrays["lammps_id"][j]),
              "symbol": atoms[j].symbol,
              "distance_A": float(neighbor_distances[j])}
             for j in range(len(atoms)) if j != max_i),
            key=lambda x: x["distance_A"],
        )[:6]
        result = {
            "label": label,
            "role": "development_validation; retained outside V17 training for now",
            "archive_structure_sha256": sha(structure),
            "marked_pair": {
                "Ag_id": int(atoms.arrays["lammps_id"][ia]),
                "X_symbol": atoms[ix].symbol,
                "X_id": int(atoms.arrays["lammps_id"][ix]),
                "distance_A": distance,
            },
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean([x["error_vector_norm_eV_A"] ** 2 for x in atom_results]))),
            "separating_force_abs_error_eV_A": evaluation["by_interface"][label.split("_", 1)[0]]["metrics"]["separating_force_abs_error_eV_A"],
            "shell_definition": "within 4.5 A MIC of either marked-pair atom",
            "shell_force_error_RMSE_eV_A": float(np.sqrt(np.mean(local_sq))),
            "outside_shell_force_error_RMSE_eV_A": float(np.sqrt(np.mean(outer_sq))),
            "largest_error_atom": max_atom,
            "largest_error_atom_nearest_neighbors": neighbors,
            "top_five_atoms": atom_results[:5],
        }
        results.append(result)

    payload = {
        "scope": "Frozen V16 predictions on only the three already authorized and scored V17 development-validation structures.",
        "model_sha256": evaluation["model_sha256"],
        "shell_metric": "RMSE of per-atom force-error vector norms; local shell is within 4.5 A MIC of either marked Ag-X atom.",
        "marked_pair_axis": "direct Ag-to-X axis, matching the frozen V16 published evaluator.",
        "results": results,
        "interpretation_limit": "Three development geometries only. Localization supports choosing future DFT environments; it does not prove causal mechanism or transferability. The three labels remain outside V17 training in this report.",
    }
    json_path = out / "localization.json"
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    lines = [
        "# V16 force-error localization on V17 development structures",
        "",
        "This is a diagnostic of the frozen V16 model on the three authorized V17 development structures. The labels remain validation data for V17; this report does not move them into training.",
        "",
        "The local shell is atoms within 4.5 Å MIC of either marked Ag–X atom. Shell RMSE is the root mean square of per-atom force-error vector norms. Directional components are projected onto the marked Ag-to-X axis. These diagnostics identify where new data may help; they do not establish the cause by themselves.",
        "",
        "| Interface | Pair distance (Å) | Vector RMSE | Separation error | 4.5 Å shell RMSE | Outside-shell RMSE | Largest residual |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        a = r["largest_error_atom"]
        lines.append(
            f"| {r['label'].split('_')[0]} | {r['marked_pair']['distance_A']:.3f} | "
            f"{r['force_vector_RMSE_eV_A']:.4f} | {r['separating_force_abs_error_eV_A']:.4f} | "
            f"{r['shell_force_error_RMSE_eV_A']:.4f} | {r['outside_shell_force_error_RMSE_eV_A']:.4f} | "
            f"{a['symbol']}#{a['atom_id']} ({a['error_vector_norm_eV_A']:.4f}) |"
        )
    lines += [
        "",
        "## What the residuals point to",
        "",
        "- **Ag–C:** the largest error is Ti#604600 (0.4590 eV/Å), with a 0.4542 eV/Å component transverse to the marked Ag–C axis. Its local neighbors include Ag#751743 at 2.154 Å and several framework C atoms at about 2.10–2.22 Å. The marked C#604804 also has 0.3024 eV/Å error. The local shell RMSE is much larger than outside the shell.",
        "- **Ag–Si:** the largest error is Si#366936 (0.2035 eV/Å); marked Si#366934 and nearby Ag#685172 also have predominantly transverse residuals. Local shell error is only moderately above the outside shell, so this is not confined to the marked pair alone.",
        "- **Ag–Ti:** marked Ti#479027 has 0.3154 eV/Å error and nearby framework C#480805 has 0.2829 eV/Å. The Ti–C distance between them is 1.869 Å. The local shell RMSE is roughly twice the outside-shell value.",
        "",
        "## Data decision",
        "",
        "Prioritize the already verified eight V17 positive/negative training diagnostics: their Ag–C, Ag–Si transverse/framework, and Ag–Ti framework moves target the same residual families. Keep these three development structures out of training for the first V17 comparison. If V17 still misses the gates, add paired DFT perturbations of the implicated framework-neighbor environments and transverse registry, then reserve a new untouched family for validation. Do not tune force weight before checking this data-coverage change with the V16 recipe fixed.",
        "",
        "The six withheld V17 structures and their labels were not accessed.",
    ]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    inputs = [atom_csv]
    inputs.extend(archive_root / name / f"{name}_PW_PBE.extxyz" for name in LABELS)
    manifest_lines = [f"{sha(p)}  {p.relative_to(repo)}" for p in [Path(__file__).resolve(), *inputs, json_path, out / "report.md"]]
    (out / "SHA256SUMS.txt").write_text("\n".join(manifest_lines) + "\n")
    print(json.dumps({r["label"]: {k: r[k] for k in ("force_vector_RMSE_eV_A", "separating_force_abs_error_eV_A", "shell_force_error_RMSE_eV_A", "outside_shell_force_error_RMSE_eV_A")} for r in results}, indent=2))


if __name__ == "__main__":
    main()
