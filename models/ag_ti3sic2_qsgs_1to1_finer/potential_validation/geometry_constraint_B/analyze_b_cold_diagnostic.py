#!/usr/bin/env python3
"""Compare the start/end snapshots from the bounded Model B cold diagnostic."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
sys.path.insert(0, str(MODEL))

from melt_metrics import largest_connected_cluster


GEOMETRY = json.loads(
    (MODEL / "potential_validation/relaxed_backbone_candidate/geometry_report.json").read_text()
)
LX, LY, MATERIAL_LZ = GEOMETRY["material_dimensions_A"]


def read_dump(path: Path) -> tuple[float, np.ndarray, dict[str, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        if stream.readline().strip() != "ITEM: TIMESTEP":
            raise ValueError(f"Missing timestep header in {path}")
        timestep = float(stream.readline())
        if stream.readline().strip() != "ITEM: NUMBER OF ATOMS":
            raise ValueError(f"Missing atom-count header in {path}")
        count = int(stream.readline())
        bounds_header = stream.readline().split()
        if bounds_header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
            raise ValueError(f"Missing box bounds in {path}")
        bounds = np.asarray(
            [[float(v) for v in stream.readline().split()[:2]] for _ in range(3)]
        )
        box = bounds[:, 1] - bounds[:, 0]
        atom_header = stream.readline().split()[2:]
        columns = {name: index for index, name in enumerate(atom_header)}
        needed = {"id", "type", "x", "y", "z", "c_cna", "c_peatom"}
        if not needed.issubset(columns):
            raise ValueError(f"Missing columns in {path}: {needed - columns.keys()}")
        data = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
    order = np.argsort(data[:, columns["id"]])
    data = data[order]
    return float(timestep), data, {**columns, "_box_x": float(box[0]), "_box_y": float(box[1]), "_box_z": float(box[2])}


def summarize_state(data: np.ndarray, columns: dict[str, int], hot: np.ndarray) -> dict:
    cna = data[:, columns["c_cna"]]
    pe = data[:, columns["c_peatom"]]
    xyz = data[:, [columns[axis] for axis in ("x", "y", "z")]]
    cluster = largest_connected_cluster(
        xyz,
        np.ones(len(xyz), dtype=bool),
        (LX, LY, MATERIAL_LZ),
        cutoff_A=3.98,
        periodic_xy=True,
    )
    return {
        "atoms": int(len(data)),
        "cna_unknown_fraction_all_mobile_ag": float(np.mean(cna == 5)),
        "cna_unknown_fraction_hot_core": float(np.mean(cna[hot] == 5)),
        "pe_sum_mobile_ag_eV": float(np.sum(pe)),
        "pe_sum_hot_core_eV": float(np.sum(pe[hot])),
        "ag_connectivity_cutoff_A": 3.98,
        "ag_largest_component_atoms": cluster.largest_cluster_atoms,
        "ag_atoms_outside_largest_component": int(len(data) - cluster.largest_cluster_atoms),
    }


def displacement_stats(displacement: np.ndarray) -> dict[str, float]:
    p50, p90, p99, p999 = np.percentile(displacement, [50, 90, 99, 99.9])
    return {
        "rms_A": float(np.sqrt(np.mean(displacement**2))),
        "p50_A": float(p50),
        "p90_A": float(p90),
        "p99_A": float(p99),
        "p99_9_A": float(p999),
        "max_A": float(np.max(displacement)),
    }


def main() -> None:
    step0, start, columns0 = read_dump(HERE / "cold_diagnostic_start.dump.gz")
    step1, final, columns1 = read_dump(HERE / "cold_diagnostic_final.dump.gz")
    ids = start[:, columns0["id"]]
    if not np.array_equal(ids, final[:, columns1["id"]]):
        raise ValueError("Atom IDs differ between start and final snapshots")
    if not np.array_equal(start[:, columns0["type"]], final[:, columns1["type"]]):
        raise ValueError("Atom types differ between start and final snapshots")

    xyz0 = start[:, [columns0[axis] for axis in ("x", "y", "z")]]
    xyz1 = final[:, [columns1[axis] for axis in ("x", "y", "z")]]
    delta = xyz1 - xyz0
    for axis, length in ((0, LX), (1, LY)):
        delta[:, axis] -= length * np.rint(delta[:, axis] / length)
    displacement = np.linalg.norm(delta, axis=1)

    center_x, center_y = LX / 2.0, LY / 2.0
    radial2 = (xyz0[:, 0] - center_x) ** 2 + (xyz0[:, 1] - center_y) ** 2
    depth = MATERIAL_LZ - xyz0[:, 2]
    hot_core = (radial2 <= 30.0**2) & (depth >= 0.0) & (depth <= 60.0)

    start_state = summarize_state(start, columns0, hot_core)
    final_state = summarize_state(final, columns1, hot_core)
    report = {
        "status": "BOUNDED_NO_SOURCE_COLD_DIAGNOSTIC; NOT AN EQUILIBRATION",
        "duration_ps": float((step1 - step0) * 0.0005),
        "start_step": int(step0),
        "final_step": int(step1),
        "hot_core_definition": "initial Ag within r_xy<=30 A and surface depth 0-60 A",
        "hot_core_atom_count": int(hot_core.sum()),
        "mobile_ag_displacement_A": displacement_stats(displacement),
        "hot_core_ag_displacement_A": displacement_stats(displacement[hot_core]),
        "start": start_state,
        "final": final_state,
        "changes": {
            "cna_unknown_fraction_all_mobile_ag": final_state["cna_unknown_fraction_all_mobile_ag"]
            - start_state["cna_unknown_fraction_all_mobile_ag"],
            "cna_unknown_fraction_hot_core": final_state["cna_unknown_fraction_hot_core"]
            - start_state["cna_unknown_fraction_hot_core"],
            "pe_sum_mobile_ag_eV": final_state["pe_sum_mobile_ag_eV"]
            - start_state["pe_sum_mobile_ag_eV"],
            "pe_sum_hot_core_eV": final_state["pe_sum_hot_core_eV"]
            - start_state["pe_sum_hot_core_eV"],
        },
        "interpretation_limits": [
            "CNA disorder over 0.25 ps is not a liquid classification; the provisional liquid flag requires a 1 ps displacement window.",
            "The 3.98 A Ag connectivity cutoff is a geometric reference from pure-Ag calibration, not proof of a liquid pool.",
            "The geometry-only WCA obstacle and rigid ceramic are exploratory assumptions.",
        ],
    }
    output = HERE / "cold_diagnostic_metrics.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
