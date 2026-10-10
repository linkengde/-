#!/usr/bin/env python3
"""Compare frozen-geometry pilot result.json files without scientific execution."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition: raise ValueError(message)


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def load(path):
    r = json.loads(path.read_text())
    required = ("task_id", "input_sha256", "geometry_sha256", "method", "method_without_kpoints_sha256", "atom_count",
                "atom_ids", "numbers", "positions_A", "cell_A", "pbc", "native_energy_eV",
                "free_energy_eV", "forces_eV_A", "scf_converged", "scf_iterations", "kpts",
                "energy_key", "force_consistent_energy_key")
    require(all(key in r for key in required), "Incomplete result schema")
    require(r["scf_converged"] is True and 0 < r["scf_iterations"] <= 160, "SCF did not converge within approved limit")
    require(r["atom_count"] == 64 and len(r["forces_eV_A"]) == 64, "Not a 64-atom liquid Ag pilot")
    require(r["numbers"] == [47] * 64 and r["pbc"] == [True] * 3, "Wrong composition/PBC")
    require(r["kpts"] in ([1, 1, 1], [2, 2, 2]), "Unapproved mesh in current pilot")
    require(r["energy_key"] == "GPAW energy (force_consistent=False)"
            and r["force_consistent_energy_key"] == "GPAW free_energy (force_consistent=True)",
            "Incorrect native/free-energy convention")
    import numpy as np
    require(np.asarray(r["positions_A"]).shape == (64, 3)
            and np.asarray(r["cell_A"]).shape == (3, 3)
            and np.isfinite(r["positions_A"]).all() and np.isfinite(r["cell_A"]).all(),
            "Incomplete/nonfinite frozen geometry")
    require(len(r["atom_ids"]) == 64 and len(set(r["atom_ids"])) == 64
            and all(type(i) is int for i in r["atom_ids"]), "Wrong immutable atom IDs")
    geometry = {"numbers": r["numbers"], "positions": np.round(r["positions_A"], 12).tolist(),
                "cell": np.round(r["cell_A"], 12).tolist(), "pbc": r["pbc"], "atom_ids": r["atom_ids"]}
    require(canonical_sha(geometry) == r["geometry_sha256"], "Frozen geometry identity hash does not match result")
    require(canonical_sha(r["method"]) == r["method_without_kpoints_sha256"], "Method identity hash mismatch")
    for energy in ("native_energy_eV", "free_energy_eV"):
        require(isinstance(r[energy], (int, float)) and math.isfinite(r[energy]), f"Missing/nonfinite {energy}")
    require(all(len(f) == 3 and all(isinstance(v, (int, float)) and math.isfinite(v) for v in f)
                for f in r["forces_eV_A"]), "Incomplete/nonfinite full force vectors")
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("low", type=Path)
    ap.add_argument("high", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    low, high = load(args.low), load(args.high)
    for key in ("task_id", "input_sha256", "geometry_sha256", "method_without_kpoints_sha256", "atom_count", "atom_ids",
                "numbers", "positions_A", "cell_A", "pbc", "energy_key", "force_consistent_energy_key"):
        require(low[key] == high[key], f"Cannot compare different geometry/order/method/convention: {key}")
    require(low["kpts"] == [1, 1, 1] and high["kpts"] == [2, 2, 2], "Compare Gamma low versus standard MP 2x2x2 high")
    squared = [sum((b - a)**2 for a, b in zip(flo, fhi))
               for flo, fhi in zip(low["forces_eV_A"], high["forces_eV_A"])]
    vector_rmse = math.sqrt(sum(squared) / len(squared))
    max_atom_delta = math.sqrt(max(squared))
    de_native = abs(high["native_energy_eV"] - low["native_energy_eV"]) * 1000 / 64
    de_free = abs(high["free_energy_eV"] - low["free_energy_eV"]) * 1000 / 64
    gates = {"native_energy_meV_atom": de_native <= 2,
             "force_consistent_free_energy_meV_atom": de_free <= 2,
             "all_atom_force_vector_RMSE_eV_A": vector_rmse <= 0.01,
             "max_atom_force_delta_eV_A": max_atom_delta <= 0.03}
    result = {"schema_version": 1, "task_id": low["task_id"], "geometry_sha256": low["geometry_sha256"],
              "low_result_sha256": file_sha(args.low), "high_result_sha256": file_sha(args.high),
              "native_energy_delta_meV_atom": de_native, "free_energy_delta_meV_atom": de_free,
              "force_vector_RMSE_eV_A": vector_rmse, "max_atom_force_delta_eV_A": max_atom_delta,
              "force_RMSE_definition": "sqrt(mean_i(sum_xyz((Fhigh-Flow)^2))); vector, not component mean",
              "gates": gates, "NUMERICAL_COMPARISON_PASS": all(gates.values()),
              "C_KPOINT_RELEASE": False, "PILOT_ACCEPTED": False,
              "note": "Same-geometry discretization control; C review required. No independent test or phase-change claim."}
    data = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as f: f.write(data)
    print(data, end="")
    return 0 if all(gates.values()) else 3


if __name__ == "__main__":
    try: raise SystemExit(main())
    except Exception as exc: raise SystemExit(f"COMPARISON_BLOCKED: {type(exc).__name__}: {exc}")
