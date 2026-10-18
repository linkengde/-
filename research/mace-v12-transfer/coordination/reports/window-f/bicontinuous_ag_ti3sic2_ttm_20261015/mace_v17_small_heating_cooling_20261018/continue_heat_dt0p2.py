#!/usr/bin/env python3
"""Continue the saved local-heating state on an independent 0.2 fs branch."""
from pathlib import Path
import json, time
import numpy as np
from ase.io import read
import run_stages as rs

root = Path(__file__).resolve().parent
atoms = read(root / "HEAT_COOL.extxyz", index=-1, format="extxyz")
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
    raise RuntimeError("Saved frame has no velocities; refusing to fabricate restart velocities")

start_wall = time.perf_counter()
result, _ = rs.base_run(
    atoms, nsteps=100, dt_fs=0.2, stage="heating_dt0p2fs",
    trajpath=root / "HEATING_DT0P2FS.extxyz",
    metrics_path=root / "HEATING_DT0P2FS_METRICS.csv",
    temperature_K=None, friction_per_fs=0.02, hot_mask=hot, sink_mask=sink,
    hot_temperature_K=1900.0, sample_every=10, seed=20261021,
    time_offset_ps=float(atoms.info.get("time_ps", 0.09)),
)
result.update({
    "branch_from": "HEAT_COOL.extxyz frame at 0.09 ps",
    "start_state_preserved_velocities": True,
    "dt_fs": 0.2,
    "hot_source_atom_count": int(hot.sum()),
    "sink_atom_count": int(sink.sum()),
    "wall_seconds_including_final_report": float(time.perf_counter() - start_wall),
    "purpose": "short timestep sensitivity / continued local-heating diagnostic; not evidence of melting",
})
(root / "HEATING_DT0P2FS_REPORT.json").write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
print(json.dumps(result, allow_nan=True), flush=True)
