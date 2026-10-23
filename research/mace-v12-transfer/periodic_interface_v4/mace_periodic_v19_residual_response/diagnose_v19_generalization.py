#!/usr/bin/env python3
"""Compare the foundation and fixed epoch-79 model on approved train36/dev3 only."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
RUN = ROOT / "runs/v19_train36_dev3_20261023_A"
FOUNDATION = ROOT.parent / "mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
DATASET_MANIFEST = ROOT / "data/dataset_manifest.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def main() -> None:
    models = {
        "foundation": FOUNDATION,
        "v19_epoch79": RUN / "models/MACE_periodic_v19_residual_response.model",
    }
    torch.set_default_dtype(torch.float64)
    calculators = {
        name: MACECalculator(model_paths=str(path), device="cpu", default_dtype="float64")
        for name, path in models.items()
    }
    manifest = json.loads(DATASET_MANIFEST.read_text())
    result: dict[str, object] = {
        "status": "DIAGNOSTIC_ONLY",
        "scope": "approved train36 and dev3 only; no test or sealed data opened",
        "models": {
            name: {"path": relative(path), "sha256": sha256(path)}
            for name, path in models.items()
        },
        "split_method_context": {
            "training_source_method": "GPAW PW-PBE 500 eV, Gamma single point",
            "training_pbc": "T T T",
            "development_method": "GPAW PW-PBE 500 eV, 6x6x1, sigma 0.10 eV",
            "development_pbc": "T T F",
            "development_energy_convention": manifest["energy_convention"],
            "development_source_group_count": manifest["B_audit"]["shared_source_group_count"],
        },
        "splits": {},
    }

    for split, filename in (("train", "train.extxyz"), ("dev", "valid.extxyz")):
        frames = read(ROOT / "data" / filename, index=":")
        element_order = ("Ag", "Ti", "Si", "C")
        composition = np.array(
            [[frame.get_chemical_symbols().count(element) for element in element_order] for frame in frames],
            dtype=float,
        )
        singular_values = np.linalg.svd(composition, compute_uv=False)
        rows: list[dict[str, object]] = []
        split_sum_sq = {name: {"force_component": 0.0, "force_vector": 0.0} for name in models}
        split_n_atoms = 0
        split_energy_sq = {name: 0.0 for name in models}
        split_energy_n = 0
        species_sum_sq = {name: Counter() for name in models}
        species_n = Counter()

        for atoms in frames:
            ref_e = float(atoms.info["REF_energy"])
            ref_f = np.asarray(atoms.arrays["REF_forces"], dtype=float)
            symbols = np.asarray(atoms.get_chemical_symbols())
            counts = {element: int(np.sum(symbols == element)) for element in sorted(set(symbols))}
            row: dict[str, object] = {
                "config_type": atoms.info.get("config_type", ""),
                "formula": atoms.get_chemical_formula(),
                "counts": counts,
                "n_atoms": len(atoms),
                "pbc": [bool(value) for value in atoms.pbc],
                "cell_A": np.asarray(atoms.cell).tolist(),
            }
            split_n_atoms += len(atoms)
            split_energy_n += 1
            for element, count in counts.items():
                species_n[element] += count
            for name, calculator in calculators.items():
                sample = atoms.copy()
                sample.calc = calculator
                pred_e = float(sample.get_potential_energy())
                pred_f = np.asarray(sample.get_forces(), dtype=float)
                delta = pred_f - ref_f
                vector_sq = np.sum(delta**2, axis=1)
                vector_norm = np.sqrt(vector_sq)
                energy_error = pred_e - ref_e
                split_sum_sq[name]["force_component"] += float(np.sum(delta**2))
                split_sum_sq[name]["force_vector"] += float(np.sum(vector_sq))
                split_energy_sq[name] += (energy_error / len(atoms)) ** 2
                per_species: dict[str, float] = {}
                for element, count in counts.items():
                    mask = symbols == element
                    species_sum_sq[name][element] += float(np.sum(vector_sq[mask]))
                    per_species[element] = float(np.sqrt(np.mean(vector_sq[mask])))
                row[name] = {
                    "energy_error_eV_cell": energy_error,
                    "energy_error_meV_atom": energy_error * 1000.0 / len(atoms),
                    "force_component_RMSE_eV_A": float(np.sqrt(np.mean(delta**2))),
                    "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(vector_sq))),
                    "force_vector_max_eV_A": float(np.max(vector_norm)),
                    "per_species_force_vector_RMSE_eV_A": per_species,
                }
            rows.append(row)

        aggregates: dict[str, object] = {}
        for name in models:
            aggregates[name] = {
                "energy_RMSE_of_structure_meV_atom": float(np.sqrt(split_energy_sq[name] / split_energy_n) * 1000.0),
                "force_component_RMSE_eV_A": float(np.sqrt(split_sum_sq[name]["force_component"] / (3 * split_n_atoms))),
                "force_vector_RMSE_eV_A": float(np.sqrt(split_sum_sq[name]["force_vector"] / split_n_atoms)),
                "per_species_force_vector_RMSE_eV_A": {
                    element: float(np.sqrt(species_sum_sq[name][element] / species_n[element]))
                    for element in sorted(species_n)
                },
            }
        result["splits"][split] = {
            "n_frames": len(frames),
            "formula_counts": dict(Counter(frame.get_chemical_formula() for frame in frames)),
            "source_method_counts": dict(Counter(str(frame.info.get("source_method", "<missing>")) for frame in frames)),
            "pbc_counts": dict(Counter(str(tuple(bool(value) for value in frame.pbc)) for frame in frames)),
            "composition_elements": list(element_order),
            "composition_matrix_rank": int(np.linalg.matrix_rank(composition)),
            "composition_matrix_singular_values": [float(value) for value in singular_values],
            "aggregates": aggregates,
            "frames": rows,
        }

    output = RUN / "results/v19_train_dev_diagnostic.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value["aggregates"] for key, value in result["splits"].items()}, indent=2))
    print(f"Wrote {relative(output)}")


if __name__ == "__main__":
    main()
