#!/usr/bin/env python3
"""Summarize the bounded 0.05 ps raw-FCC, no-source Model B screen."""

from __future__ import annotations

import json
import gzip
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read_dump(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        if stream.readline().strip() != "ITEM: TIMESTEP":
            raise ValueError(f"Bad dump header: {path}")
        step = int(stream.readline())
        stream.readline()  # NUMBER OF ATOMS
        count = int(stream.readline())
        box_header = stream.readline().split()
        if box_header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
            raise ValueError(f"Bad box header: {path}")
        bounds = np.asarray([[float(x) for x in stream.readline().split()[:2]] for _ in range(3)])
        header = stream.readline().split()[2:]
        columns = {name: i for i, name in enumerate(header)}
        raw = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
    order = np.argsort(raw[:, columns["id"]])
    raw = raw[order]
    return step, bounds, columns, raw


def thermo(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 7:
            continue
        try:
            values = [float(x) for x in fields]
        except ValueError:
            continue
        if values[0].is_integer() and values[0] in (0, 100):
            rows.append(dict(zip(("step", "time_ps", "T_Ag_K", "PE_eV", "KE_eV", "Etot_eV", "max_force_eV_A"), values)))
    unique = {int(row["step"]): row for row in rows}
    if not {0, 100}.issubset(unique):
        raise ValueError(f"Could not parse step 0 and 100 thermo rows from {path}")
    return [unique[0], unique[100]]


def main():
    step0, bounds0, col0, start = read_dump(HERE / "raw_fcc_cold_start.dump.gz")
    step1, bounds1, col1, final = read_dump(HERE / "raw_fcc_cold_final.dump.gz")
    if step0 != 0 or step1 != 100 or not np.array_equal(start[:, col0["id"]], final[:, col1["id"]]):
        raise ValueError("Unexpected dump steps or atom ID mismatch")

    p0 = start[:, [col0[k] for k in ("x", "y", "z")]]
    p1 = final[:, [col1[k] for k in ("x", "y", "z")]]
    delta = p1 - p0
    lengths = bounds0[:, 1] - bounds0[:, 0]
    for axis in (0, 1):
        delta[:, axis] -= np.rint(delta[:, axis] / lengths[axis]) * lengths[axis]
    disp = np.linalg.norm(delta, axis=1)
    f0 = np.linalg.norm(start[:, [col0[k] for k in ("fx", "fy", "fz")]], axis=1)
    f1 = np.linalg.norm(final[:, [col1[k] for k in ("fx", "fy", "fz")]], axis=1)
    t = thermo(HERE / "log.raw_fcc_cold_screen_0p05ps")
    metrics = {
        "status": "RAW_GENERATED_FCC_AG_NO_SOURCE_NVE_0P05PS_SCREEN",
        "steps": 100,
        "dt_fs": 0.5,
        "duration_ps": 0.05,
        "mobile_Ag_atoms": int(len(disp)),
        "thermo_step_0": t[0],
        "thermo_step_100": t[1],
        "total_energy_drift_eV": float(t[1]["Etot_eV"] - t[0]["Etot_eV"]),
        "start_max_force_eV_A": float(f0.max()),
        "end_max_force_eV_A": float(f1.max()),
        "displacement_A": {
            "rms": float(np.sqrt(np.mean(disp**2))),
            "p50": float(np.percentile(disp, 50)),
            "p90": float(np.percentile(disp, 90)),
            "p99": float(np.percentile(disp, 99)),
            "max": float(disp.max()),
        },
        "interpretation": [
            "No external heat source or TTM was active.",
            "NVE energy conservation tests integration, not cold-state stability.",
            "This is the raw generated Ag FCC coordinate set, before the prior Ag-only FIRE-100 preparation.",
        ],
    }
    (HERE / "raw_fcc_cold_screen_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
