#!/usr/bin/env python3
"""Analyze the bounded 0.1 ps source-core-only cold NVE screen."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


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
        columns = {name: i for i, name in enumerate(names)}
        data = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
    order = np.argsort(data[:, columns["id"]])
    return step, bounds, columns, data[order]


def thermo_rows(path: Path):
    found = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 7:
            continue
        try:
            row = [float(x) for x in fields]
        except ValueError:
            continue
        if row[0] in (100.0, 300.0):
            found[int(row[0])] = row
    if set(found) != {100, 300}:
        raise ValueError("Missing step 100/300 thermodynamic rows")
    keys = ("step", "time_ps", "T_hot_K", "PE_eV", "KE_eV", "Etot_eV", "max_force_eV_A")
    return [dict(zip(keys, found[k])) for k in (100, 300)]


def stats(x):
    return {
        "atoms": int(len(x)),
        "rms_A": float(np.sqrt(np.mean(x**2))),
        "p90_A": float(np.percentile(x, 90)),
        "p99_A": float(np.percentile(x, 99)),
        "max_A": float(np.max(x)),
    }


def main():
    s0, bounds, c0, start = read_dump(HERE / "hot_core_cold_start.dump.gz")
    s1, _, c1, end = read_dump(HERE / "hot_core_cold_final.dump.gz")
    if s0 != 100 or s1 != 300 or not np.array_equal(start[:, c0["id"]], end[:, c1["id"]]):
        raise ValueError("Unexpected hot-core dump steps or IDs")
    xyz0 = start[:, [c0[k] for k in ("x", "y", "z")]]
    xyz1 = end[:, [c1[k] for k in ("x", "y", "z")]]
    d = xyz1 - xyz0
    lengths = bounds[:, 1] - bounds[:, 0]
    for axis in (0, 1):
        d[:, axis] -= np.rint(d[:, axis] / lengths[axis]) * lengths[axis]
    disp = np.linalg.norm(d, axis=1)
    f0 = np.linalg.norm(start[:, [c0[k] for k in ("fx", "fy", "fz")]], axis=1)
    f1 = np.linalg.norm(end[:, [c1[k] for k in ("fx", "fy", "fz")]], axis=1)
    t0, t1 = thermo_rows(HERE / "log.hot_core_cold_nve_0p1ps")
    result = {
        "status": "SOURCE_CORE_ONLY_NO_SOURCE_NVE_0P1PS",
        "mobile_Ag_atoms": int(len(disp)),
        "dt_fs": 0.5,
        "duration_ps": 0.1,
        "thermo_start": t0,
        "thermo_end": t1,
        "total_energy_drift_eV": float(t1["Etot_eV"] - t0["Etot_eV"]),
        "potential_energy_change_eV": float(t1["PE_eV"] - t0["PE_eV"]),
        "kinetic_energy_change_eV": float(t1["KE_eV"] - t0["KE_eV"]),
        "force_max_start_eV_A": float(f0.max()),
        "force_max_end_eV_A": float(f1.max()),
        "force_p95_start_end_eV_A": [float(np.percentile(f0, 95)), float(np.percentile(f1, 95))],
        "displacement": stats(disp),
        "model_limits": [
            "Only the nominal source-core Ag is mobile; all other Ag and the ceramic are fixed.",
            "No heat source or TTM was active.",
            "The fixed surrounding atoms are an exploratory local-boundary condition.",
        ],
    }
    (HERE / "hot_core_cold_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
