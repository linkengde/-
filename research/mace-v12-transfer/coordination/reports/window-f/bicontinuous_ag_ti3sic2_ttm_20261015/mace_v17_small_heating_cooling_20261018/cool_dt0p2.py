#!/usr/bin/env python3
"""Turn off the hot Ag reservoir and continue the saved state through the TSC sink."""
from pathlib import Path
import json, time
import numpy as np
from ase.io import read
import run_stages as rs

root = Path(__file__).resolve().parent
atoms = read(root / "HEATING_DT0P2FS_CHECKPOINT.extxyz", format="extxyz")
rs.attach_model(atoms)
heat_start = read(root / "HEAT_COOL.extxyz", index=0, format="extxyz")
heat_layer = np.asarray(heat_start.arrays["Ag111_layer"], dtype=int)
heat_phase = np.asarray(heat_start.arrays["phase"])
heat_radial = rs.lateral_distance_to_center(heat_start)
hot_start_ids = np.asarray(heat_start.arrays["id"])[(heat_phase == "Ag") & (heat_layer >= 3) & (heat_radial <= 8.0)]
hot = np.isin(np.asarray(atoms.arrays["id"]), hot_start_ids)
phase = np.asarray(atoms.arrays["phase"])
sink = (phase == "Ti3SiC2") & (atoms.positions[:, 2] > 0.05) & (atoms.positions[:, 2] <= 4.50)
sink &= ~np.asarray(atoms.arrays["support_fixed"], dtype=bool)
if not atoms.has("momenta"):
    raise RuntimeError("Saved heating checkpoint has no velocities")

offset = float(atoms.info.get("time_ps", 0.11))
start_wall = time.perf_counter()
result, _ = rs.base_run(
    atoms, nsteps=100, dt_fs=0.2, stage="cooling_dt0p2fs",
    trajpath=root / "COOLING_DT0P2FS.extxyz",
    metrics_path=root / "COOLING_DT0P2FS_METRICS.csv",
    temperature_K=None, friction_per_fs=0.02, hot_mask=hot, sink_mask=sink,
    hot_temperature_K=1900.0, sample_every=10, seed=20261022,
    time_offset_ps=offset,
)
result.update({
    "branch_from": "HEATING_DT0P2FS_CHECKPOINT.extxyz",
    "heat_source_off": True,
    "only_heat_removal": "300 K Langevin sink on mobile Ti3SiC2 support-side region; other atoms NVE",
    "wall_seconds_including_final_report": float(time.perf_counter() - start_wall),
    "interpretation": "thermal cooldown trajectory only; not evidence of resolidification unless prior melting is established",
})
(root / "COOLING_DT0P2FS_REPORT.json").write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
print(json.dumps(result, allow_nan=True), flush=True)
