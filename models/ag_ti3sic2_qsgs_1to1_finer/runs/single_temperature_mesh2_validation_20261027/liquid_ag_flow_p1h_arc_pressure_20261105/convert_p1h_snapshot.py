#!/usr/bin/env python3
"""Combine one P1H thermal/VOF frame and its pressure sidecar into HDF5."""
from __future__ import annotations

import argparse
import gzip
import struct
from pathlib import Path

import h5py
import numpy as np

HEADER = struct.Struct("<8s6I5d")
MAIN_NAMES = [
    "vof", "level_set_m", "ag_marker", "skeleton_marker", "liquid_fraction",
    "temperature_K", "enthalpy_density_J_m3", "velocity_x_m_s", "velocity_y_m_s",
    "velocity_z_m_s", "pressure_Pa", "curvature_m-1", "surface_tension_candidate_N_m",
    "vof_weighted_surface_tension_N_m", "marangoni_acceleration_candidate_x_m_s2",
    "marangoni_acceleration_candidate_y_m_s2", "marangoni_acceleration_candidate_z_m_s2",
    "source_power_density_W_m3", "skeleton_and_mushy_drag_s-1", "material_id",
]
PRESSURE_NAMES = [
    "arc_pressure_Pa", "arc_pressure_mask", "interface_normal_pressure_Pa",
    "total_normal_stress_potential_Pa",
]


def read_frame(path: Path, magic: bytes, version: int, nfields: int):
    raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    if len(raw) < HEADER.size:
        raise ValueError(f"short frame: {path}")
    fields_magic, got_version, nx, ny, nz, got_fields, case_flag, time_ps, dx, x0, y0, z0 = HEADER.unpack_from(raw)
    if (fields_magic, got_version, got_fields) != (magic, version, nfields):
        raise ValueError(f"unexpected P1H header in {path}")
    expected = HEADER.size + 4*nfields*nx*ny*nz
    if len(raw) != expected:
        raise ValueError(f"partial frame {path}: expected {expected} bytes, got {len(raw)}")
    fields = np.frombuffer(raw, dtype="<f4", offset=HEADER.size).reshape(nfields, nx, ny, nz)
    if not np.isfinite(fields).all():
        raise ValueError(f"non-finite values in {path}")
    meta = dict(shape=(nx, ny, nz), time_ps=time_ps, dx_m=dx,
                origin_m=(x0, y0, z0), case_flag=case_flag)
    return fields, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", type=Path, required=True)
    ap.add_argument("--pressure", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    main, mm = read_frame(args.main, b"P1HSNAP1", 1, len(MAIN_NAMES))
    pressure, pm = read_frame(args.pressure, b"P1HPRS01", 1, len(PRESSURE_NAMES))
    for key in ("shape", "time_ps", "dx_m", "origin_m"):
        if mm[key] != pm[key]:
            raise ValueError(f"main and pressure frame metadata differ at {key}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["time_ps"] = mm["time_ps"]
        h5.attrs["grid_shape_xyz"] = mm["shape"]
        h5.attrs["grid_spacing_m"] = mm["dx_m"]
        h5.attrs["origin_xyz_m"] = mm["origin_m"]
        h5.attrs["marangoni_enabled"] = bool(mm["case_flag"])
        h5.attrs["pressure_frame_enabled"] = bool(pm["case_flag"])
        fields = h5.create_group("fields")
        for name, data in zip(MAIN_NAMES, main):
            ds = fields.create_dataset(name, data=data, compression="gzip", compression_opts=4, shuffle=True)
            ds.attrs["axis_order"] = "x,y,z"
        pressure_group = h5.create_group("arc_pressure")
        units = ("Pa", "1", "Pa", "Pa")
        for name, unit, data in zip(PRESSURE_NAMES, units, pressure):
            ds = pressure_group.create_dataset(name, data=data, compression="gzip", compression_opts=4, shuffle=True)
            ds.attrs["unit"] = unit
            ds.attrs["axis_order"] = "x,y,z"
    print(f"wrote {args.output} ({args.output.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
