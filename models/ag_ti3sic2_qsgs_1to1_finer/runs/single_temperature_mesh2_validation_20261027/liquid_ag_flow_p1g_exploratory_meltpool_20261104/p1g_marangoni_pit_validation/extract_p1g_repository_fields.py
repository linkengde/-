#!/usr/bin/env python3
"""Write a GitHub-sized HDF5 subset with the P1G thermal/flow handoff fields."""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py


KEEP_FIELDS = [
    "vof", "level_set_m", "ag_marker", "skeleton_marker", "liquid_fraction",
    "temperature_K", "velocity_x_m_s", "velocity_y_m_s", "velocity_z_m_s",
    "pressure_Pa", "surface_tension_candidate_N_m", "vof_weighted_surface_tension_N_m",
    "marangoni_acceleration_candidate_x_m_s2", "marangoni_acceleration_candidate_y_m_s2",
    "marangoni_acceleration_candidate_z_m_s2",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(a.source, "r") as src, h5py.File(a.output, "w") as dst:
        for key, value in src.attrs.items():
            dst.attrs[key] = value
        dst.attrs["repository_subset"] = True
        dst.attrs["repository_subset_note"] = "Selected thermal, phase, geometry, pressure, velocity, and candidate Marangoni fields. Full HDF5 remains in the cloud run directory."
        src.copy("time_ps", dst)
        src.copy("mesh", dst)
        fields = dst.create_group("fields")
        for name in KEEP_FIELDS:
            ds = src["fields/" + name]
            out = fields.create_dataset(
                name, shape=ds.shape, dtype=ds.dtype,
                chunks=(1,) + ds.shape[1:], compression="gzip", compression_opts=4, shuffle=True,
            )
            for k, v in ds.attrs.items():
                out.attrs[k] = v
            for i in range(ds.shape[0]):
                out[i] = ds[i]
        for name in ["free_surface_height_m", "free_surface_temperature_K"]:
            src.copy(name, dst)
    with h5py.File(a.output, "r") as check:
        frames = len(check["time_ps"])
    print(f"output={a.output} bytes={a.output.stat().st_size} frames={frames}")


if __name__ == "__main__":
    main()
