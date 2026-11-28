#!/usr/bin/env python3
"""独立运行一个RVE有效热导率算例并保存残差、热流分解和可选场快照。"""

from __future__ import annotations

import argparse
import json
import resource
from pathlib import Path
import sys
import time

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
sys.path.insert(0, str(HERE))
from steady_periodic_fvm import PeriodicConductivityFVM  # noqa: E402

CONFIG = json.loads((HERE / "VALIDATION_CONFIG.json").read_text(encoding="utf-8"))

MATERIALS = {
    "source_file": "geometry_generation/source_extracts/single-temperature-study/solver/thermophysical_candidates.csv",
    "Ag": {"k_300K_W_m_K": CONFIG["materials_at_temperature_K"]["Ag"]["thermal_conductivity_W_m_K"]},
    "Ti3SiC2": {"k_300K_W_m_K": CONFIG["materials_at_temperature_K"]["Ti3SiC2"]["thermal_conductivity_W_m_K"]},
}
SIZES_M = {name: item["domain_edge_m"] for name, item in CONFIG["rves"].items()}


def load_geometry(grid_n: int, seed: int, kind: str) -> tuple[np.ndarray, int]:
    if kind == "fixed_nested_300":
        source = PROJECT / "geometry_generation" / f"qsgs_seed_{seed}.npz"
        with np.load(source, allow_pickle=False) as archive:
            mask = np.asarray(archive["skeleton"], dtype=bool)
        if mask.shape != (150, 150, 150) or grid_n != 300:
            raise ValueError("fixed_nested_300 requires a committed 150³ base and target N=300")
        mask = np.repeat(np.repeat(np.repeat(mask, 2, axis=0), 2, axis=1), 2, axis=2)
        return mask, 150

    if kind == "performance" and grid_n in (50, 100):
        base = PROJECT / "geometry_generation" / f"qsgs_seed_{seed}.npz"
        with np.load(base, allow_pickle=False) as archive:
            source_mask = np.asarray(archive["skeleton"], dtype=bool)
        # Nearest source-cell-centre sampling is for runtime/RSS timing only.
        indices = np.floor((np.arange(grid_n) + 0.5) * 150.0 / grid_n).astype(np.intp)
        mask = source_mask[np.ix_(indices, indices, indices)]
        return mask, 150

    if grid_n == 150:
        source = PROJECT / "geometry_generation" / f"qsgs_seed_{seed}.npz"
    else:
        source = PROJECT / "resolution_study" / "geometries" / f"n{grid_n}" / f"qsgs_seed_{seed}.npz"
    with np.load(source, allow_pickle=False) as archive:
        mask = np.asarray(archive["skeleton"], dtype=bool)
    return mask, grid_n


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-n", type=int, required=True, choices=(50, 100, 150, 200, 300))
    parser.add_argument("--seed", type=int, required=True, choices=(1105, 1136, 1182))
    parser.add_argument("--size-name", choices=tuple(SIZES_M), required=True)
    parser.add_argument("--interface-resistance", type=float, default=0.0, help="m² K/W")
    parser.add_argument("--kind", choices=("performance", "qsgs_morphology", "fixed_nested_300"), required=True)
    parser.add_argument("--snapshot", action="store_true", help="save x/y/z imposed-gradient temperature and flux fields")
    parser.add_argument("--relative-tolerance", type=float,
                        default=CONFIG["boundary_and_solver"]["relative_tolerance"])
    args = parser.parse_args()

    wall_start = time.perf_counter()
    geometry_start = time.perf_counter()
    mask, source_grid_n = load_geometry(args.grid_n, args.seed, args.kind)
    geometry_seconds = time.perf_counter() - geometry_start
    size_m = SIZES_M[args.size_name]
    spacing = np.full(3, size_m / args.grid_n, dtype=np.float64)
    solver = PeriodicConductivityFVM(
        mask,
        spacing,
        MATERIALS["Ag"]["k_300K_W_m_K"],
        MATERIALS["Ti3SiC2"]["k_300K_W_m_K"],
        args.interface_resistance,
        fft_workers=int(CONFIG["boundary_and_solver"]["FFT_workers"]),
    )
    result = solver.solve_effective_tensor(
        relative_tolerance=args.relative_tolerance,
        maximum_iterations=int(CONFIG["boundary_and_solver"]["maximum_CG_iterations"]),
        snapshot_axis="all" if args.snapshot else None,
        snapshot_delta_temperature_K=float(CONFIG["boundary_and_solver"]["imposed_macro_temperature_difference_K"]),
        cold_temperature_K=float(CONFIG["boundary_and_solver"]["reference_temperature_K"]),
    )
    tensor = np.asarray(result["effective_conductivity_tensor_W_m_K"])
    diagonal = np.diag(tensor)
    offdiagonal = tensor - np.diag(diagonal)
    symmetry_abs = float(np.max(np.abs(tensor - tensor.T)))
    record: dict[str, object] = {
        "status": "PASS" if max(s["relative_residual"] for s in result["solver"]) <= args.relative_tolerance else "FAIL",
        "case_id": f"{args.kind}_n{args.grid_n}_seed{args.seed}_{args.size_name}_R{args.interface_resistance:.0e}",
        "geometry_kind": args.kind,
        "geometry_file_grid_n": source_grid_n,
        "grid_n": args.grid_n,
        "seed": args.seed,
        "shape_xyz": list(mask.shape),
        "control_volumes": int(mask.size),
        "size_name": args.size_name,
        "domain_size_m": [size_m] * 3,
        "spacing_m": spacing.tolist(),
        "skeleton_volume_fraction": float(mask.mean()),
        "ag_volume_fraction": float((~mask).mean()),
        "interface_resistance_m2_K_W": args.interface_resistance,
        "boundary_conditions": CONFIG["boundary_and_solver"]["boundary_condition"] + " in x/y/z",
        "temperature_for_material_k_K": CONFIG["materials_at_temperature_K"]["temperature_K"],
        "imposed_temperature_delta_K_for_snapshots": CONFIG["boundary_and_solver"]["imposed_macro_temperature_difference_K"],
        "thermal_conductivity_inputs_W_m_K": MATERIALS,
        "effective_conductivity_tensor_W_m_K": tensor.tolist(),
        "diagonal_W_m_K": diagonal.tolist(),
        "maximum_absolute_offdiagonal_W_m_K": float(np.max(np.abs(offdiagonal))),
        "maximum_tensor_symmetry_error_W_m_K": symmetry_abs,
        "phase_flux_contributions": result["phase_flux_contributions"],
        "interface_heat_flux_statistics": result["interface_heat_flux_statistics"],
        "solver_by_imposed_axis": result["solver"],
        "geometry_load_seconds": geometry_seconds,
        "total_wall_seconds": time.perf_counter() - wall_start,
        "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "snapshot_saved": bool(args.snapshot),
        "performance_geometry_mapping": (
            "nearest source-cell-centre sample; timing only, not used for convergence"
            if args.kind == "performance" and args.grid_n in (50, 100) else None
        ),
    }
    output_dir = HERE / "results" / "cases"
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{record['case_id']}.json"
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    if args.snapshot:
        for axis, snapshot in result["snapshots_by_imposed_axis"].items():
            arrays = {k: v for k, v in snapshot.items() if isinstance(v, np.ndarray)}
            np.savez_compressed(output_dir / f"{record['case_id']}_axis{axis}_fields.npz", **arrays)
            (output_dir / f"{record['case_id']}_axis{axis}_interface_flux.json").write_text(
                json.dumps(snapshot["interface_stats"], indent=2) + "\n", encoding="utf-8"
            )
    print(json.dumps({
        "status": record["status"], "case_id": record["case_id"],
        "tensor_W_m_K": tensor.tolist(),
        "iterations": [s["iterations"] for s in result["solver"]],
        "residuals": [s["relative_residual"] for s in result["solver"]],
        "wall_seconds": record["total_wall_seconds"],
        "peak_rss_mib": record["process_peak_rss_mib"],
        "output": str(output),
    }), flush=True)
    if record["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
