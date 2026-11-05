#!/usr/bin/env python3
"""Validate compressed P1H rectangular final snapshots using only stdlib."""
import array
import csv
import gzip
import math
import pathlib
import struct
import sys


HEADER = struct.Struct("<8s6I5d")


def load_fields(path, expected_magic):
    with gzip.open(path, "rb") as f:
        header_bytes = f.read(HEADER.size)
        if len(header_bytes) != HEADER.size:
            raise ValueError(f"short header: {path}")
        magic, version, nx, ny, nz, nfields, flags, time_ps, dx_m, x0, y0, z0 = HEADER.unpack(header_bytes)
        if magic != expected_magic or version != 1:
            raise ValueError(f"unexpected format {magic!r} v{version}: {path}")
        ncells = nx * ny * nz
        payload = f.read()
    if len(payload) != nfields * ncells * 4:
        raise ValueError(f"payload size mismatch: {path}: {len(payload)} bytes")
    values = array.array("f")
    values.frombytes(payload)
    if values.itemsize != 4:
        raise ValueError("host float size is not 32 bits")
    if sys.byteorder != "little":
        values.byteswap()
    if any(not math.isfinite(v) for v in values):
        raise ValueError(f"NaN/Inf in {path}")
    return (nx, ny, nz, nfields, flags, time_ps, dx_m, values, ncells)


def field_range(values, field, ncells):
    start = field * ncells
    segment = values[start:start + ncells]
    return min(segment), max(segment)


def main(root):
    rows = []
    for nz in (64, 76, 86):
        for boundary in ("closed", "open"):
            folder = root / f"nz{nz}_{boundary}"
            snap = folder / "p1h_marangoni_off_max_time_t000.160000_ps.bin.gz"
            pressure = folder / "p1h_pressure_marangoni_off_max_time_t000.160000_ps.bin.gz"
            nx, ny, nzz, nfields, flags, t, dx, values, ncells = load_fields(snap, b"P1HSNAP1")
            if (nx, ny, nzz, nfields, t) != (64, 64, nz, 20, 0.16):
                raise ValueError(f"unexpected dimensions/fields/time in {snap}")
            pnx, pny, pnz, pnfields, pflags, pt, pdx, pvalues, pncells = load_fields(pressure, b"P1HPRS01")
            if (pnx, pny, pnz, pnfields, pt, pncells) != (64, 64, nz, 4, 0.16, ncells):
                raise ValueError(f"unexpected dimensions/fields/time in {pressure}")
            fmin, fmax = field_range(values, 0, ncells)
            lmin, lmax = field_range(values, 4, ncells)
            Tmin, Tmax = field_range(values, 5, ncells)
            matmin, matmax = field_range(values, 19, ncells)
            pmin, pmax = field_range(pvalues, 2, ncells)
            if lmin < -1.e-6 or lmax > 1.000001:
                raise ValueError(f"liquid fraction out of range: {snap}: {lmin}..{lmax}")
            if matmin < -1.e-6 or matmax > 2.000001:
                raise ValueError(f"material code out of range: {snap}: {matmin}..{matmax}")
            rows.append({
                "grid": f"{nx}x{ny}x{nzz}", "top_boundary": boundary,
                "time_ps": t, "dx_nm": dx * 1.e9, "field_count": nfields,
                "finite_values": len(values), "vof_min": fmin, "vof_max": fmax,
                "liquid_fraction_min": lmin, "liquid_fraction_max": lmax,
                "temperature_min_K": Tmin, "temperature_max_K": Tmax,
                "material_code_min": matmin, "material_code_max": matmax,
                "effective_pressure_min_Pa": pmin, "effective_pressure_max_Pa": pmax,
                "field_snapshot_sha256": __import__("hashlib").sha256(snap.read_bytes()).hexdigest(),
                "pressure_snapshot_sha256": __import__("hashlib").sha256(pressure.read_bytes()).hexdigest(),
                "validation": "PASS",
            })
    out = root / "snapshot_validation.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"validated {len(rows)} final field/pressure snapshot pairs; wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_p1h_rectangular_snapshots.py OUTPUT_ROOT")
    main(pathlib.Path(sys.argv[1]))
