#!/usr/bin/env python3
"""Summarize the short exploratory Ag-W NVE pilot."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
BOX = np.array([200.0, 200.0, 1000.0])  # x/y periodic; z effectively nonperiodic at these distances
ELEMENT = {1: "W", 2: "Ag"}
FLOOR = {"W_W": 2.60, "Ag_Ag": 2.48, "Ag_W": 2.54}


def read_dump(path: Path) -> np.ndarray:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        a = np.loadtxt(f, skiprows=9)
    return a


def thermo_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    log = (HERE / "log.nve_pilot_0p025ps").read_text(errors="replace")
    names = ("step", "time_ps", "atoms", "PE_eV", "KE_eV", "Etot_eV", "T_mobile_K", "T_W_K", "T_Ag_K", "Fmax_mobile_eV_A", "Frms_mobile_eV_A", "zmin_A", "zmax_A")
    for line in log.splitlines():
        fields = line.split()
        if len(fields) != len(names) or not fields[0].isdigit():
            continue
        try:
            vals = [float(x) for x in fields]
        except ValueError:
            continue
        rows.append(dict(zip(names, vals)))
    if len(rows) < 2:
        raise RuntimeError(f"Expected initial/final thermo rows, found {len(rows)}")
    return rows


def dist_min(points: np.ndarray, same: bool = False) -> float:
    tree = cKDTree(points, boxsize=BOX)
    if same:
        return float(tree.query(points, k=2, workers=-1)[0][:, 1].min())
    raise ValueError("Use cross_min for unlike-species distances")


def cross_min(a: np.ndarray, b: np.ndarray) -> float:
    tree = cKDTree(a, boxsize=BOX)
    return float(tree.query(b, workers=-1)[0].min())


def main() -> None:
    start = read_dump(HERE / "nve_start.dump.gz")
    end = read_dump(HERE / "nve_final.dump.gz")
    if start.shape != end.shape or not np.array_equal(start[:, :2], end[:, :2]):
        raise RuntimeError("Atom IDs/types changed between pilot snapshots")
    types = start[:, 1].astype(int)
    pos0, pos1 = start[:, 2:5].copy(), end[:, 2:5].copy()
    pos0[:, :2] %= 200.0
    pos1[:, :2] %= 200.0
    delta = pos1 - pos0
    for ax in (0, 1):
        delta[:, ax] -= 200.0 * np.rint(delta[:, ax] / 200.0)
    displacement = np.linalg.norm(delta, axis=1)
    force = np.linalg.norm(end[:, 8:11], axis=1)
    # The LAMMPS block region includes its boundary faces; z=20.0 A atoms are fixed.
    mobile = pos0[:, 2] > 20.0
    species: dict[str, dict[str, object]] = {}
    for tid, name in ELEMENT.items():
        mask = (types == tid) & mobile
        species[name] = {
            "mobile_atoms": int(mask.sum()),
            "displacement_p50_p95_p99_max_A": [float(x) for x in np.percentile(displacement[mask], [50, 95, 99, 100])],
            "fraction_moved_over_1A": float(np.mean(displacement[mask] > 1.0)),
            "end_force_p50_p95_p99_max_eV_A": [float(x) for x in np.percentile(force[mask], [50, 95, 99, 100])],
        }
    rows = thermo_rows()
    initial, final = rows[0], rows[-1]
    min_dists = {}
    for name, p0, p1, same in (
        ("W_W", pos0[types == 1], pos1[types == 1], True),
        ("Ag_Ag", pos0[types == 2], pos1[types == 2], True),
    ):
        min_dists[name] = {"start_A": dist_min(p0, same), "end_A": dist_min(p1, same), "screen_A_context_only": FLOOR[name]}
    min_dists["Ag_W"] = {
        "start_A": cross_min(pos0[types == 1], pos0[types == 2]),
        "end_A": cross_min(pos1[types == 1], pos1[types == 2]),
        "screen_A_context_only": FLOOR["Ag_W"],
    }
    report = {
        "title": "Ag-W exploratory Dai-EFS 300K no-source NVE pilot",
        "status": "RUN_COMPLETED_EXPLORATORY_BASELINE_UNSTABLE",
        "decision": "LEVEL_1_PROGRAM_CHECK_PASSED; COLD_BASELINE_NOT_ESTABLISHED; DO_NOT_START_TTM_OR_LASER_FROM_THIS_ENDPOINT",
        "source": "v5 signed-normal/same-phase precondition candidate followed by finite spring-free Dai EFS minimization with only bottom z<20A fixed",
        "boundary": "p p f",
        "fixed_region": "atoms initially at z<20 A",
        "timestep_fs": 0.25,
        "steps": int(final["step"] - initial["step"]),
        "duration_ps": float(final["time_ps"] - initial["time_ps"]),
        "atoms_start_end": [int(start.shape[0]), int(end.shape[0])],
        "mobile_atom_count": int(mobile.sum()),
        "NVE_thermo_start": initial,
        "NVE_thermo_end": final,
        "energy_drift_eV": float(final["Etot_eV"] - initial["Etot_eV"]),
        "energy_drift_fraction_abs_start": float(abs(final["Etot_eV"] - initial["Etot_eV"]) / abs(initial["Etot_eV"])),
        "thermo_energy_display_resolution_eV": 0.1,
        "PE_change_eV": float(final["PE_eV"] - initial["PE_eV"]),
        "KE_change_eV": float(final["KE_eV"] - initial["KE_eV"]),
        "mobile_temperature_change_K": float(final["T_mobile_K"] - initial["T_mobile_K"]),
        "maximum_mobile_force_change_eV_A": float(final["Fmax_mobile_eV_A"] - initial["Fmax_mobile_eV_A"]),
        "minimum_pair_distances": min_dists,
        "mobile_species_motion_and_force": species,
        "z_bounds_A": {"start": [float(pos0[:, 2].min()), float(pos0[:, 2].max())], "end": [float(pos1[:, 2].min()), float(pos1[:, 2].max())]},
        "limitations": [
            "No external heating, thermostat, TTM, pulse, or heat-flux analysis was run.",
            "The original QSGS inputs remain unchanged; the test candidate uses the v5 volume-preserving mask regularization and is not a byte-for-byte reproduction of the original voxel mask.",
            "NVE energy conservation verifies integration bookkeeping, not physical validity of the EFS alloy potential.",
            "The strong PE-to-KE exchange and temperature drop indicate ongoing structural relaxation, so this is not an equilibrated 300 K baseline.",
            "Distance screens are reported for context and are not treated as pass/fail gates for this exploratory test.",
        ],
    }
    (HERE / "nve_pilot_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
