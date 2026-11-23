#!/usr/bin/env python3
"""Summarize the full-size exploratory QSGS EFS static-relaxation output."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "research/ag_w_reference_pack/scripts"))
from project_full_qsgs_v3 import read_atoms, distance_and_pairs, topology_summary  # noqa: E402


def read_dump(path: Path):
    with path.open() as stream:
        lines = stream.readlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    columns = lines[header].split()[2:]
    index = {name: columns.index(name) for name in ("id", "type", "x", "y", "z", "fx", "fy", "fz")}
    records = np.asarray([[float(value) for value in line.split()] for line in lines[header + 1 :]], dtype=float)
    order = np.argsort(records[:, index["id"]])
    records = records[order]
    pos = records[:, [index[k] for k in ("x", "y", "z")]]
    types = records[:, index["type"]].astype(np.int8)
    forces = records[:, [index[k] for k in ("fx", "fy", "fz")]]
    return records[:, index["id"]].astype(int), types, pos, forces


def main() -> None:
    source = HERE / "ag_w_qsgs_v3_exploratory_geometry.data.gz"
    ids, types, pos, forces = read_dump(HERE / "minimized.dump")
    initial = read_atoms(source)
    initial_by_id = initial[np.argsort(initial[:, 0])]
    if not np.array_equal(ids, initial_by_id[:, 0].astype(int)):
        raise ValueError("atom IDs do not match between input and minimized dump")
    if not np.array_equal(types, initial_by_id[:, 1].astype(int)):
        raise ValueError("atom types changed during minimization")

    anchor = initial_by_id[:, 4] < 20.0
    mobile = ~anchor
    force_norm = np.linalg.norm(forces, axis=1)
    w, ag, cross, same_w, same_ag, min_cross, min_w, min_ag = distance_and_pairs(pos, types)
    topo_pos = pos.copy()
    topo_pos[:, :2] %= 200.0
    topo_z_offset = max(0.0, 1.0e-6 - float(pos[:, 2].min()))
    topo_pos[:, 2] += topo_z_offset  # keep z nonperiodic while satisfying cKDTree's positive box domain
    topology = topology_summary(topo_pos, types)
    for phase in topology.values():
        phase["largest_component_z_extent_A"] = [z - topo_z_offset for z in phase["largest_component_z_extent_A"]]
    topology_pass_95 = all(
        phase["xy_winding_rank"] == 2
        and phase["largest_component_spans_z"]
        and phase["largest_component_share"] >= 0.95
        for phase in topology.values()
    )

    displacement = pos - initial_by_id[:, 2:5]
    displacement[:, :2] = (displacement[:, :2] + 100.0) % 200.0 - 100.0
    displacement_norm = np.linalg.norm(displacement, axis=1)
    mass_ratio = np.count_nonzero(types == 2) * 107.86819997225152 / (
        np.count_nonzero(types == 1) * 183.83999995270821
    )
    summary = {
        "title": "Ag-W full-size QSGS exploratory EFS static relaxation audit",
        "status": "NOT_READY_FOR_THERMAL_DYNAMICS",
        "input_geometry": source.name,
        "minimizer": "LAMMPS conjugate-gradient; 10000 iterations; stopped at iteration cap",
        "potential": "Dai 2006 EFS converted to LAMMPS eam/fs",
        "boundary": "p p f",
        "fixed_region_A": "z=0..20",
        "counts": {"W": int(len(w)), "Ag": int(len(ag)), "total": int(len(ids))},
        "mass_ratio_Ag_over_W": float(mass_ratio),
        "atoms_outside_material_z_0_200_A": int(np.count_nonzero((pos[:, 2] < -1e-8) | (pos[:, 2] > 200.0 + 1e-8))),
        "position_bounds_A": {axis: [float(pos[:, i].min()), float(pos[:, i].max())] for i, axis in enumerate(("x", "y", "z"))},
        "minimum_distances_A": {"Ag_W": min_cross, "W_W": min_w, "Ag_Ag": min_ag},
        "pairs_below_local_scan_extrema": {
            "Ag_W_below_2_54_A": int(len(cross[0])),
            "W_W_below_2_60_A": int(len(same_w)),
            "Ag_Ag_below_2_48_A": int(len(same_ag)),
        },
        "forces_eV_A": {
            "mobile_max_vector": float(force_norm[mobile].max()),
            "mobile_rms_vector": float(np.sqrt(np.mean(force_norm[mobile] ** 2))),
            "mobile_99th_percentile_vector": float(np.percentile(force_norm[mobile], 99)),
            "anchored_max_vector": float(force_norm[anchor].max()),
            "criterion_max_eV_A": 0.1,
            "criterion_rms_eV_A": 0.02,
            "pass": bool(force_norm[mobile].max() <= 0.1 and np.sqrt(np.mean(force_norm[mobile] ** 2)) <= 0.02),
        },
        "energy_eV": -2421889.9132037,
        "energy_eV_atom": -2421889.9132037 / len(ids),
        "maximum_displacement_A": float(displacement_norm.max()),
        "displacement_percentiles_A": {str(q): float(np.percentile(displacement_norm, q)) for q in (50, 90, 95, 99)},
        "atoms_displaced_over_1_A": int(np.count_nonzero(displacement_norm > 1.0)),
        "topology": topology,
        "bicontinuous_practical_95_percent_screen": bool(topology_pass_95),
        "static_force_screen_pass": bool(force_norm[mobile].max() <= 0.1 and np.sqrt(np.mean(force_norm[mobile] ** 2)) <= 0.02),
        "thermal_dynamics_started": False,
        "interpretation": "The minimizer reached its 10000-step cap with residual forces above the stated criterion. It also produced distances below the local-collapse screen and atoms above the material surface; do not start TTM or melt/quench dynamics from this structure.",
    }
    (HERE / "static_efs_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
