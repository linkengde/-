#!/usr/bin/env python3
"""Analyze local-ROI cold NVE displacements and force changes."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LX = 301.9555789947458
LY = 304.2005152739243
LZ = 159.39156144562426
CX, CY = LX / 2.0, LY / 2.0


def read_dump(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        stream.readline()
        step = int(stream.readline())
        stream.readline()
        count = int(stream.readline())
        stream.readline()
        bounds = np.asarray([[float(v) for v in stream.readline().split()[:2]] for _ in range(3)])
        names = stream.readline().split()[2:]
        cols = {name: i for i, name in enumerate(names)}
        values = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
    order = np.argsort(values[:, cols["id"]])
    return step, bounds, cols, values[order]


def displacement_stats(disp):
    return {
        "atoms": int(len(disp)),
        "rms_A": float(np.sqrt(np.mean(disp**2))),
        "p90_A": float(np.percentile(disp, 90)),
        "p99_A": float(np.percentile(disp, 99)),
        "max_A": float(np.max(disp)),
    }


def main():
    s0, bounds, c0, start = read_dump(HERE / "local_roi_cold_start.dump.gz")
    s1, _, c1, final = read_dump(HERE / "local_roi_cold_final.dump.gz")
    if s0 != 100 or s1 != 600 or not np.array_equal(start[:, c0["id"]], final[:, c1["id"]]):
        raise ValueError("Unexpected steps or atom IDs in local ROI dumps")
    xyz0 = start[:, [c0[a] for a in ("x", "y", "z")]]
    xyz1 = final[:, [c1[a] for a in ("x", "y", "z")]]
    delta = xyz1 - xyz0
    lengths = bounds[:, 1] - bounds[:, 0]
    for axis in (0, 1):
        delta[:, axis] -= np.rint(delta[:, axis] / lengths[axis]) * lengths[axis]
    disp = np.linalg.norm(delta, axis=1)
    r = np.sqrt((xyz0[:, 0] - CX) ** 2 + (xyz0[:, 1] - CY) ** 2)
    depth = LZ - xyz0[:, 2]
    source_core = (r <= 30.0) & (depth >= 0.0) & (depth <= 60.0)
    f0 = np.linalg.norm(start[:, [c0[a] for a in ("fx", "fy", "fz")]], axis=1)
    f1 = np.linalg.norm(final[:, [c1[a] for a in ("fx", "fy", "fz")]], axis=1)

    thermo_rows = []
    for line in (HERE / "log.local_roi_cold_nve_0p25ps").read_text().splitlines():
        fields = line.split()
        if len(fields) == 7:
            try:
                row = [float(x) for x in fields]
            except ValueError:
                continue
            if row[0] in (100.0, 600.0):
                thermo_rows.append(row)
    rows = {int(row[0]): row for row in thermo_rows}
    if not {100, 600}.issubset(rows):
        raise ValueError("Could not parse start/end thermo from local ROI log")
    keys = ("step", "time_ps", "T_local_K", "PE_eV", "KE_eV", "Etot_eV", "max_force_eV_A")
    t0, t1 = (dict(zip(keys, rows[step])) for step in (100, 600))
    metrics = {
        "status": "GEOMETRY_CONSTRAINED_LOCAL_ROI_NO_SOURCE_NVE_0P25PS",
        "dt_fs": 0.5,
        "duration_ps": 0.25,
        "roi_definition": "Ag within cylinder r<=75 A around x/y center and 20<=z<=159.3916 A; all other Ag and all Ti/Si/C fixed",
        "source_core_definition": "Ag within r<=30 A and surface depth<=60 A",
        "source_core_atoms": int(source_core.sum()),
        "thermo_start": t0,
        "thermo_end": t1,
        "total_energy_drift_eV": float(t1["Etot_eV"] - t0["Etot_eV"]),
        "max_force_start_mobile_Ag_eV_A": float(f0.max()),
        "max_force_end_mobile_Ag_eV_A": float(f1.max()),
        "roi_displacement": displacement_stats(disp),
        "source_core_displacement": displacement_stats(disp[source_core]),
        "source_core_force_start_p95_max_eV_A": [float(np.percentile(f0[source_core], 95)), float(f0[source_core].max())],
        "source_core_force_end_p95_max_eV_A": [float(np.percentile(f1[source_core], 95)), float(f1[source_core].max())],
        "interpretation": [
            "No external heat source or TTM was active.",
            "A large potential-to-kinetic exchange with good total-energy conservation indicates structural relaxation rather than integration error.",
            "Rigid ceramic and fixed Ag outside the ROI are exploratory boundary assumptions.",
        ],
    }
    (HERE / "local_roi_cold_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
