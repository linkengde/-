#!/usr/bin/env python3
"""Search fixed-orientation Ag/W integer repeat matches under frozen v3 limits."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.build import bulk, make_supercell


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
LAT = {"Ag": 4.0900219316, "W": 3.1599967891}
MAX_REPEAT = 40
MAX_AREA_A2 = 1800.0
MAX_ABS_STRAIN = 0.01

ORIENTATIONS = {
    "Ag100_W100": {
        "Ag": [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
        "W": [[0, 1, 1], [0, -1, 1], [1, 0, 0]],
        "relationship": "[010]Ag || [011]W",
    },
    "Ag111_W110": {
        "Ag": [[-2, 1, 1], [0, -1, 1], [1, 1, 1]],
        "W": [[1, -1, 0], [0, 0, -1], [1, 1, 0]],
        "relationship": "[-2,1,1]Ag || [1,-1,0]W",
    },
}


def orient_cell_axes(atoms):
    rotation = (atoms.cell.array / atoms.cell.lengths()[:, None]).T
    atoms.positions = atoms.positions @ rotation
    atoms.set_cell(atoms.cell.array @ rotation)
    return atoms


def oriented_vectors(element: str, matrix: list[list[int]]) -> np.ndarray:
    structure = "fcc" if element == "Ag" else "bcc"
    base = orient_cell_axes(
        make_supercell(bulk(element, structure, a=LAT[element], cubic=True), matrix)
    )
    return base.cell.lengths()[:2]


def axis_pairs(ag_len: float, w_len: float) -> list[dict]:
    candidates = []
    for m_ag in range(1, MAX_REPEAT + 1):
        for m_w in range(1, MAX_REPEAT + 1):
            natural_ag = m_ag * ag_len
            natural_w = m_w * w_len
            common = float(np.sqrt(natural_ag * natural_w))
            strain_ag = common / natural_ag - 1.0
            strain_w = common / natural_w - 1.0
            candidates.append(
                {
                    "Ag_repeat": m_ag,
                    "W_repeat": m_w,
                    "Ag_natural_length_A": natural_ag,
                    "W_natural_length_A": natural_w,
                    "common_length_A": common,
                    "Ag_strain": strain_ag,
                    "W_strain": strain_w,
                    "max_abs_strain": max(abs(strain_ag), abs(strain_w)),
                }
            )
    return candidates


def main() -> None:
    report = {
        "title": "Ag-W v3 low-strain coincidence-cell search",
        "equilibrium_lattice_A": LAT,
        "search": {
            "maximum_repeat_per_direction": MAX_REPEAT,
            "maximum_area_A2": MAX_AREA_A2,
            "maximum_absolute_phase_strain": MAX_ABS_STRAIN,
            "common_length_rule": "geometric mean of natural repeat lengths (equal log-strain)",
            "search_scope": "independent integer repeats along the two fixed orthogonal in-plane direction pairs",
        },
        "orientations": {},
    }
    for name, spec in ORIENTATIONS.items():
        ag_vec = oriented_vectors("Ag", spec["Ag"])
        w_vec = oriented_vectors("W", spec["W"])
        x_pairs = axis_pairs(float(ag_vec[0]), float(w_vec[0]))
        y_pairs = axis_pairs(float(ag_vec[1]), float(w_vec[1]))
        cells = []
        for x in x_pairs:
            for y in y_pairs:
                area = x["common_length_A"] * y["common_length_A"]
                max_strain = max(x["max_abs_strain"], y["max_abs_strain"])
                if area <= MAX_AREA_A2 and max_strain <= MAX_ABS_STRAIN:
                    cells.append(
                        {
                            "Ag_repeats_xy": [x["Ag_repeat"], y["Ag_repeat"]],
                            "W_repeats_xy": [x["W_repeat"], y["W_repeat"]],
                            "natural_lengths_A": {
                                "Ag": [x["Ag_natural_length_A"], y["Ag_natural_length_A"]],
                                "W": [x["W_natural_length_A"], y["W_natural_length_A"]],
                            },
                            "common_lengths_A": [x["common_length_A"], y["common_length_A"]],
                            "strains": {
                                "Ag": [x["Ag_strain"], y["Ag_strain"]],
                                "W": [x["W_strain"], y["W_strain"]],
                            },
                            "max_abs_strain": max_strain,
                            "area_A2": area,
                        }
                    )
        cells.sort(key=lambda item: (item["area_A2"], item["max_abs_strain"]))
        report["orientations"][name] = {
            "orientation_relationship": spec["relationship"],
            "natural_inplane_vectors_A": {"Ag": ag_vec.tolist(), "W": w_vec.tolist()},
            "feasible_cell_count": len(cells),
            "selected_cell": cells[0] if cells else None,
            "next_lowest_area_candidates": cells[1:6],
        }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "cell_match_report.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
