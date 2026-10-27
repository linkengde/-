#!/usr/bin/env python3
"""Geometry-only reproducibility checks for the two paper-derived CIF variants."""

from pathlib import Path
import json

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
CASES = {
    "abstract": {
        "path": "Ti3SiC2_Jeitschko1967_abstract_reconstruction.cif",
        "minimum_distance_A": 2.1354093859584813,
    },
    "fulltext_fourier": {
        "path": "Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif",
        "minimum_distance_A": 2.11781059527386,
    },
}


def inspect(path: Path) -> dict:
    atoms = read(path)
    assert len(atoms) == 12, (path.name, len(atoms))
    assert atoms.get_chemical_formula() == "C4Si2Ti6"
    sg = atoms.info.get("spacegroup")
    assert sg is not None and sg.no == 194, (path.name, sg)
    cellpar = atoms.cell.cellpar()
    assert np.allclose(cellpar, [3.068, 3.068, 17.669, 90, 90, 120], atol=1e-8)
    distances = atoms.get_all_distances(mic=True)
    distances[distances < 1e-8] = np.inf
    return {
        "file": path.name,
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "space_group_number": sg.no,
        "cell_parameters": cellpar.tolist(),
        "volume_A3": float(atoms.get_volume()),
        "minimum_periodic_distance_A": float(distances.min()),
    }


def main() -> None:
    results = {}
    for name, spec in CASES.items():
        result = inspect(ROOT / spec["path"])
        assert abs(result["minimum_periodic_distance_A"] - spec["minimum_distance_A"]) < 1e-9
        results[name] = result
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
