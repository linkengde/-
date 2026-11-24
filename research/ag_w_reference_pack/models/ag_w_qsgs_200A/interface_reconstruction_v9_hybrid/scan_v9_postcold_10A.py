#!/usr/bin/env python3
"""Rescan 10 A CIC phase on the final V9 cold NVE restart with shared flow."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
TTM = MODEL / "ttm"
sys.path.insert(0, str(TTM))
sys.path.insert(0, str(HERE))
from analyze_v9_cold_baseline import read_dump_loose  # noqa: E402
from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3,
    W_ATOMIC_VOLUME_A3,
    _raw_cic_volume_power,
    gaussian_ag_source_power,
)
from diagnose_agw_ttm_mapping_v7 import local_capacity_flow, read_dump, read_lammps_data_atoms  # noqa: E402

REPO = Path(__file__).resolve().parents[5]
DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
RESTART = HERE / "v9_nve2_end.dump.gz"
OUT = HERE / "v9_postcold_10A_phase_scan.json"
BIN_EDGES = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 150.0])


def source_summary(power: np.ndarray, dx: float, dy: float) -> dict:
    nx, ny, nz = power.shape
    xx = (dx + (np.arange(nx) + 0.5) * 10.0) % 200.0
    yy = (dy + (np.arange(ny) + 0.5) * 10.0) % 200.0
    x, y = np.meshgrid(xx, yy, indexing="ij")
    x = (x - 100.0 + 100.0) % 200.0 - 100.0
    y = (y - 100.0 + 100.0) % 200.0 - 100.0
    radius = np.sqrt(x * x + y * y)[:, :, None]
    r = np.broadcast_to(radius, power.shape)
    total = float(power.sum())
    bins = np.digitize(r.ravel(), BIN_EDGES[1:-1])
    profile = np.bincount(bins, weights=power.ravel(), minlength=len(BIN_EDGES) - 1)
    return {
        "power_eV_ps": total,
        "radial_bins_A": BIN_EDGES.tolist(),
        "radial_power_fraction": (profile / total).tolist() if total else None,
        "radial_RMS_A": float(np.sqrt(np.sum(power * r * r) / total)) if total else None,
    }


def main() -> None:
    ids0, types0, xyz0 = read_lammps_data_atoms(DATA)
    ids, types, xyz = read_dump_loose(RESTART)
    if not (np.array_equal(ids0, ids) and np.array_equal(types0, types)):
        raise RuntimeError("V9 final restart changed IDs/types")
    static_newton = (xyz0[:, 2] >= 40.0 - 1e-9) & (xyz0[:, 2] <= 300.0 + 1e-9)
    mapped = static_newton & (xyz[:, 2] >= 40.0 - 1e-9) & (xyz[:, 2] <= 200.0 + 1e-9)
    i_ag = np.flatnonzero(mapped & (types == 2))
    i_w = np.flatnonzero(mapped & (types == 1))
    q_ag = gaussian_ag_source_power(xyz[i_ag])
    exact_power = float(np.dot(np.full(len(i_ag), AG_ATOMIC_VOLUME_A3), q_ag))
    rows = []
    offsets = np.arange(0.0, 10.0, 2.5)
    for dx in offsets:
        for dy in offsets:
            origin = np.array([dx, dy, 40.0])
            _, va, source = _raw_cic_volume_power(
                xyz[i_ag] - origin, q_ag, (20, 20, 16), (10.0,)*3, AG_ATOMIC_VOLUME_A3
            )
            _, vw, _ = _raw_cic_volume_power(
                xyz[i_w] - origin, np.zeros(len(i_w)), (20, 20, 16), (10.0,)*3, W_ATOMIC_VOLUME_A3
            )
            total = va + vw
            flow = local_capacity_flow(total, (10.0,)*3, 20.0)
            rows.append({
                "offset_xy_A": [float(dx), float(dy)],
                "raw_fill_max_Ag": float(va.max() / 1000.0),
                "raw_fill_max_W": float(vw.max() / 1000.0),
                "raw_fill_max_combined": float(total.max() / 1000.0),
                "raw_overflow_A3": float(np.maximum(total - 1000.0, 0.0).sum()),
                "Ag_volume_A3": float(va.sum()),
                "W_volume_A3": float(vw.sum()),
                "source": source_summary(source, float(dx), float(dy)),
                "source_quadrature_error_eV_ps": float(source.sum() - exact_power),
                "shared_direct_max_flow": flow,
            })
    best = min(rows, key=lambda r: r["shared_direct_max_flow"]["locally_unassignable_overflow_A3"])
    result = {
        "title": "V9 post-cold-continuation 10 A grid-phase scan",
        "input": {
            "restart": RESTART.name,
            "restart_sha256": hashlib.sha256(RESTART.read_bytes()).hexdigest(),
            "data_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
            "atoms_total": int(len(ids)),
            "mapped_Ag": int(len(i_ag)),
            "mapped_W": int(len(i_w)),
            "outside_static_Newton_region": int(np.count_nonzero(~static_newton)),
            "static_Newton_atoms_outside_40_to_200_A": int(np.count_nonzero(static_newton & ~((xyz[:, 2] >= 40.0 - 1e-9) & (xyz[:, 2] <= 200.0 + 1e-9)))),
            "box_A": [200.0, 200.0, 300.0],
            "grid_A": 10.0,
            "phase_step_A": 2.5,
            "max_transfer_A": 20.0,
            "atomic_volumes_unchanged": True,
            "atom_count_or_coordinates_modified": False,
            "source_reference_power_eV_ps": exact_power,
        },
        "phase_count": len(rows),
        "phase_results": rows,
        "best_by_shared_capacity_flow": {
            "offset_xy_A": best["offset_xy_A"],
            "unassigned_overflow_A3": best["shared_direct_max_flow"]["locally_unassignable_overflow_A3"],
            "overflow_demand_A3": best["shared_direct_max_flow"]["overflow_demand_A3"],
            "source": best["source"],
        },
        "decision": "PASS only if any phase receives a complete direct assignment; no phase uses cascaded transfers.",
        "ttm_or_laser_run": False,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "best_offset": best["offset_xy_A"],
        "best_unassigned_A3": best["shared_direct_max_flow"]["locally_unassignable_overflow_A3"],
        "all_pass": any(r["shared_direct_max_flow"]["status"] == "PASS" for r in rows),
        "mapped_Ag_W": [len(i_ag), len(i_w)],
        "source_eV_ps": exact_power,
    }, indent=2))


if __name__ == "__main__":
    main()
