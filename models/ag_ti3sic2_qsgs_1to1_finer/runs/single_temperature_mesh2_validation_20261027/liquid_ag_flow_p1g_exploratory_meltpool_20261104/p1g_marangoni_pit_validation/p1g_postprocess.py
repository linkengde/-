#!/usr/bin/env python3
"""Convert atomic P1G fixed-grid snapshots to HDF5 and derive pit audits."""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import struct
from pathlib import Path

import h5py
import numpy as np

HEADER = struct.Struct("<8s6I5d")
FIELD_NAMES = [
    "vof", "level_set_m", "ag_marker", "skeleton_marker", "liquid_fraction",
    "temperature_K", "enthalpy_density_J_m3", "velocity_x_m_s", "velocity_y_m_s",
    "velocity_z_m_s", "pressure_Pa", "curvature_m-1", "surface_tension_candidate_N_m",
    "vof_weighted_surface_tension_N_m", "marangoni_acceleration_candidate_x_m_s2",
    "marangoni_acceleration_candidate_y_m_s2", "marangoni_acceleration_candidate_z_m_s2",
    "source_power_density_W_m3", "skeleton_and_mushy_drag_s-1", "material_id",
]
FIELD_UNITS = {
    "vof": "1",
    "level_set_m": "m",
    "ag_marker": "1",
    "skeleton_marker": "1",
    "liquid_fraction": "1",
    "temperature_K": "K",
    "enthalpy_density_J_m3": "J/m^3",
    "velocity_x_m_s": "m/s",
    "velocity_y_m_s": "m/s",
    "velocity_z_m_s": "m/s",
    "pressure_Pa": "Pa",
    "curvature_m-1": "1/m",
    "surface_tension_candidate_N_m": "N/m",
    "vof_weighted_surface_tension_N_m": "N/m",
    "marangoni_acceleration_candidate_x_m_s2": "m/s^2",
    "marangoni_acceleration_candidate_y_m_s2": "m/s^2",
    "marangoni_acceleration_candidate_z_m_s2": "m/s^2",
    "source_power_density_W_m3": "W/m^3",
    "skeleton_and_mushy_drag_s-1": "1/s",
    "material_id": "1",
}


def read_snapshot(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            raw = f.read()
    else:
        raw = path.read_bytes()
    if len(raw) < HEADER.size:
        raise ValueError(f"short snapshot: {path}")
    magic, version, nx, ny, nz, nfields, marangoni, time_ps, dx, x0, y0, z0 = HEADER.unpack_from(raw)
    if magic != b"P1GSNAP2" or version != 2 or nfields != len(FIELD_NAMES):
        raise ValueError(f"bad header/version/field count in {path}")
    expected = HEADER.size + nfields * nx * ny * nz * 4
    if len(raw) != expected:
        raise ValueError(f"partial snapshot {path}: {len(raw)} bytes, expected {expected}")
    fields = np.frombuffer(raw, dtype="<f4", offset=HEADER.size).reshape((nfields, nx, ny, nz))
    if not np.isfinite(fields).all():
        raise ValueError(f"non-finite field value in {path}")
    return dict(time_ps=time_ps, dx_m=dx, origin_m=(x0, y0, z0),
                shape=(nx, ny, nz), marangoni_enabled=bool(marangoni)), fields


def grad_periodic(a: np.ndarray, axis: int, dx: float) -> np.ndarray:
    return (np.roll(a, -1, axis=axis) - np.roll(a, 1, axis=axis))/(2.*dx)


def derive(fields: np.ndarray, meta):
    dx = meta["dx_m"]
    levelset = fields[1].astype(np.float64)
    grad = [grad_periodic(levelset, 0, dx), grad_periodic(levelset, 1, dx),
            np.gradient(levelset, dx, axis=2, edge_order=1)]
    # d = z_surface - z, so -grad(d) points from liquid to gas.
    normal = [-g for g in grad]
    nmag = np.sqrt(sum(g*g for g in normal))
    normal = [np.divide(g, nmag, out=np.zeros_like(g), where=nmag > 1.e-20) for g in normal]

    ux, uy, uz = [fields[q].astype(np.float64) for q in (7, 8, 9)]
    dvx = [grad_periodic(ux, 0, dx), grad_periodic(ux, 1, dx), np.gradient(ux, dx, axis=2)]
    dvy = [grad_periodic(uy, 0, dx), grad_periodic(uy, 1, dx), np.gradient(uy, dx, axis=2)]
    dvz = [grad_periodic(uz, 0, dx), grad_periodic(uz, 1, dx), np.gradient(uz, dx, axis=2)]
    vort = [dvz[1] - dvy[2], dvx[2] - dvz[0], dvy[0] - dvx[1]]

    vof = fields[0].astype(np.float64)
    nz = vof.shape[2]
    occupied = vof > 0.5
    # Match the solver's surface-height diagnostic: use the highest occupied
    # z cell in each x-y column.  Summing VOF through z gives a volume-equivalent
    # column height, which is not a free-surface elevation for porous/composite
    # columns containing internal voids.
    kk = np.max(np.where(occupied, np.arange(nz)[None, None, :], -1), axis=2)
    has_surface = kk >= 0
    z0 = meta["origin_m"][2]
    height = z0 + (kk + 0.5)*dx
    height = height.astype(np.float64)
    height[~has_surface] = np.nan
    safe_kk = np.clip(kk, 0, nz-1)
    temp_surface = np.take_along_axis(fields[5], safe_kk[..., None], axis=2)[..., 0].copy()
    temp_surface[~has_surface] = np.nan
    return normal, vort, height, temp_surface


def csv_rows(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def process(snapshot_dir: Path, output: Path):
    files = sorted(list(snapshot_dir.glob("p1g_*_t*_ps.bin")) +
                   list(snapshot_dir.glob("p1g_*_t*_ps.bin.gz")))
    if not files:
        raise SystemExit(f"no P1GSNAP2 snapshots in {snapshot_dir}")
    frames = []
    for p in files:
        meta, fields = read_snapshot(p)
        frames.append((p, meta, fields))
    frames.sort(key=lambda x: (x[1]["time_ps"], x[0].name))
    first = frames[0][1]
    if any(fr[1]["shape"] != first["shape"] or fr[1]["dx_m"] != first["dx_m"] for fr in frames):
        raise ValueError("snapshot mesh changed; this postprocessor expects frozen-grid output")
    nx, ny, nz = first["shape"]
    output.parent.mkdir(parents=True, exist_ok=True)
    pit_rows = []
    height0 = None
    with h5py.File(output, "w") as h5:
        h5.attrs["solver"] = "P1G audited exploratory CLSVOF/enthalpy prototype"
        h5.attrs["grid_shape_xyz"] = first["shape"]
        h5.attrs["cell_size_m"] = first["dx_m"]
        h5.attrs["axis_order"] = "time,x,y,z"
        h5.attrs["pit_height_method"] = "highest z cell with VOF>0.5 per x-y column; upper-envelope height; columns without VOF are NaN"
        mesh = h5.create_group("mesh")
        dx = first["dx_m"]
        x0, y0, z0 = first["origin_m"]
        mesh.create_dataset("x_center_m", data=x0+(np.arange(nx)+0.5)*dx)
        mesh.create_dataset("y_center_m", data=y0+(np.arange(ny)+0.5)*dx)
        mesh.create_dataset("z_center_m", data=z0+(np.arange(nz)+0.5)*dx)
        mesh.attrs["periodic_axes"] = "x,y"
        fields_group = h5.create_group("fields")
        h5.create_dataset("time_ps", shape=(len(frames),), dtype="f8")
        for name in FIELD_NAMES:
            ds = fields_group.create_dataset(name, shape=(len(frames),nx,ny,nz), dtype="f4",
                chunks=(1,nx,ny,nz), compression="gzip", compression_opts=4, shuffle=True)
            ds.attrs["units"] = FIELD_UNITS[name]
        normal_ds = fields_group.create_dataset("interface_normal_xyz", shape=(len(frames),3,nx,ny,nz),
            dtype="f4", chunks=(1,3,nx,ny,nz), compression="gzip", compression_opts=4, shuffle=True)
        normal_ds.attrs["units"] = "1"
        vort_ds = fields_group.create_dataset("vorticity_xyz_s-1", shape=(len(frames),3,nx,ny,nz),
            dtype="f4", chunks=(1,3,nx,ny,nz), compression="gzip", compression_opts=4, shuffle=True)
        vort_ds.attrs["units"] = "1/s"
        surf_ds = h5.create_dataset("free_surface_height_m", shape=(len(frames),nx,ny), dtype="f4",
            chunks=(1,nx,ny), compression="gzip", compression_opts=4, shuffle=True)
        surf_ds.attrs["units"] = "m"
        tsurf_ds = h5.create_dataset("free_surface_temperature_K", shape=(len(frames),nx,ny), dtype="f4",
            chunks=(1,nx,ny), compression="gzip", compression_opts=4, shuffle=True)
        tsurf_ds.attrs["units"] = "K"
        for idx, (path, meta, data) in enumerate(frames):
            h5["time_ps"][idx] = meta["time_ps"]
            for fi, name in enumerate(FIELD_NAMES):
                fields_group[name][idx] = data[fi]
            normal, vort, height, tsurf = derive(data, meta)
            normal_ds[idx] = np.stack(normal).astype("f4")
            vort_ds[idx] = np.stack(vort).astype("f4")
            surf_ds[idx] = height.astype("f4")
            tsurf_ds[idx] = tsurf.astype("f4")
            if height0 is None:
                height0 = height.copy()
            dep = height0 - height
            threshold = 0.5*dx
            valid = np.isfinite(dep)
            pit = valid & (dep > threshold)
            rim = valid & (dep < -threshold)
            pit_rows.append({
                "time_ps": meta["time_ps"],
                "marangoni_enabled": int(meta["marangoni_enabled"]),
                "pit_depth_max_nm": float(max(0.,np.max(dep[pit]) if np.any(pit) else 0.)*1.e9),
                "pit_area_nm2": float(np.count_nonzero(pit)*dx*dx*1.e18),
                "pit_volume_nm3": float(np.where(pit,dep,0.).sum()*dx*dx*1.e27),
                "rim_volume_nm3": float(np.where(rim,-dep,0.).sum()*dx*dx*1.e27),
                "max_surface_rise_nm": float(max(0.,np.max(-dep[rim]) if np.any(rim) else 0.)*1.e9),
                "surface_height_method": "highest z with VOF>0.5 per x-y column; upper envelope; pit/rim threshold=0.5dx",
                "snapshot": path.name,
            })
        h5.attrs["marangoni_enabled"] = first["marangoni_enabled"]
        h5.attrs["frame_count"] = len(frames)
    write_csv(snapshot_dir/"P1G_PIT_GEOMETRY.csv",list(pit_rows[0]),pit_rows)

    hist = csv_rows(snapshot_dir/"p1g_history.csv")
    if hist:
        write_csv(snapshot_dir/"P1G_CASE_SUMMARY.csv",list(hist[0]),hist)
        energy_cols=["time_ps","source_power_W","enthalpy_J","input_energy_J","bottom_sink_energy_J","energy_residual_J"]
        write_csv(snapshot_dir/"P1G_ENERGY_BALANCE.csv",energy_cols,[{k:r[k] for k in energy_cols} for r in hist])
        mass_cols=["time_ps","Ag_volume_nm3","Ag_volume_drift_pct","liquid_Ag_volume_nm3"]
        write_csv(snapshot_dir/"P1G_AG_MASS_BALANCE.csv",mass_cols,[{k:r[k] for k in mass_cols} for r in hist])
    print(f"frames={len(frames)} grid={nx}x{ny}x{nz} hdf5={output} bytes={output.stat().st_size}")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("snapshot_dir",type=Path)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    output=args.output or args.snapshot_dir/"thermal_flow_fields.h5"
    process(args.snapshot_dir,output)


if __name__=="__main__":
    main()
