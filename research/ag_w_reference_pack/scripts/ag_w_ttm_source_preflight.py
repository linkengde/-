#!/usr/bin/env python3
"""Map the Ag-only candidate TTM source onto the actual Ag-W atom geometry.

This is a static source/mask preflight. It does not run LAMMPS, TTM, or MD.
The script is self-contained and does not import the legacy Ag-Ti3SiC2 project.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "models/ag_w_qsgs_200A/ag_w_qsgs.data.gz"
DEFAULT_OUT = ROOT / "models/ag_w_qsgs_200A/ttm_source_preflight"

BOX_A = np.array([200.0, 200.0, 300.0])
MATERIAL_Z_A = (0.0, 200.0)
NEWTON_Z_A = (40.0, 200.0)
SOURCE_CENTER_XY_A = np.array([100.0, 100.0])
SOURCE_SURFACE_Z_A = 200.0
I0_EV_PS_A2 = 5.0
R_A = 30.0
SKIN_DEPTH_A = 20.0
AG_LATTICE_A = 4.09
W_LATTICE_A = 3.16
AG_ATOMIC_VOLUME_A3 = AG_LATTICE_A**3 / 4.0
W_ATOMIC_VOLUME_A3 = W_LATTICE_A**3 / 2.0
GRID_SPACINGS_A = (5.0, 10.0, 20.0)
OUTPUT_GRID_A = 10.0


def read_lammps_atomic_data(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read type/x/y/z and box lengths from a LAMMPS atom_style atomic file."""
    opener = gzip.open if path.suffix == ".gz" else open
    bounds: dict[str, tuple[float, float]] = {}
    atom_types: list[int] = []
    positions: list[tuple[float, float, float]] = []
    in_atoms = False
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if len(fields) >= 4 and fields[2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
                axis = fields[2][0]
                bounds[axis] = (float(fields[0]), float(fields[1]))
                continue
            if fields[0] == "Atoms":
                in_atoms = True
                continue
            if not in_atoms:
                continue
            try:
                atom_type = int(fields[1])
                xyz = tuple(float(value) for value in fields[2:5])
            except (ValueError, IndexError):
                # Atom section data are numeric; section headers are not.
                if atom_types:
                    break
                continue
            if atom_type not in (1, 2):
                raise ValueError(f"Expected type 1=W or 2=Ag, got type {atom_type}")
            atom_types.append(atom_type)
            positions.append(xyz)
    if len(bounds) != 3 or not atom_types:
        raise ValueError(f"Could not parse box bounds and atoms from {path}")
    origin = np.array([bounds[axis][0] for axis in "xyz"], dtype=float)
    lengths = np.array([bounds[axis][1] - bounds[axis][0] for axis in "xyz"], dtype=float)
    return np.asarray(atom_types, dtype=np.int8), np.asarray(positions, dtype=float) - origin, lengths


def source_density_eV_ps_A3(xyz_A: np.ndarray) -> np.ndarray:
    """S(t)=1 spatial source density using periodic x/y minimum-image distance."""
    dxy = (xyz_A[:, :2] - SOURCE_CENTER_XY_A + BOX_A[:2] / 2.0) % BOX_A[:2] - BOX_A[:2] / 2.0
    depth_A = SOURCE_SURFACE_Z_A - xyz_A[:, 2]
    if np.any(depth_A < -1.0e-9):
        raise ValueError("Atoms above z_surface found; define a moving-surface source before heating")
    return (I0_EV_PS_A2 / SKIN_DEPTH_A) * np.exp(
        -np.sum(dxy * dxy, axis=1) / (R_A * R_A) - np.maximum(depth_A, 0.0) / SKIN_DEPTH_A
    )


def cloud_in_cell_map(
    xyz_A: np.ndarray, shape: tuple[int, int, int], spacing_A: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Trilinear weights with periodic x/y and open/clamped z faces."""
    n = np.asarray(shape, dtype=int)
    u = xyz_A / spacing_A - 0.5
    lower = np.floor(u).astype(int)
    frac = u - lower
    cell_ids: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for bx, by, bz in itertools.product((0, 1), repeat=3):
        ix = (lower[:, 0] + bx) % n[0]
        iy = (lower[:, 1] + by) % n[1]
        iz = np.clip(lower[:, 2] + bz, 0, n[2] - 1)
        wx = frac[:, 0] if bx else 1.0 - frac[:, 0]
        wy = frac[:, 1] if by else 1.0 - frac[:, 1]
        wz = frac[:, 2] if bz else 1.0 - frac[:, 2]
        cell_ids.append(np.ravel_multi_index((ix, iy, iz), shape))
        weights.append(wx * wy * wz)
    ids = np.column_stack(cell_ids)
    wts = np.column_stack(weights)
    if not np.allclose(wts.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
        raise RuntimeError("CIC weights do not sum to one")
    return ids, wts


def _redistribute_overfill(
    volume: np.ndarray,
    power: np.ndarray,
    shape: tuple[int, int, int],
    spacing_A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int, float, float]:
    """Bound the Ag fill fraction and preserve volume and source power locally."""
    cell_volume = float(np.prod(spacing_A))
    raw_fill = volume / cell_volume
    raw_overfilled = int(np.count_nonzero(raw_fill > 1.0 + 1.0e-12))
    raw_max_fill = float(raw_fill.max(initial=0.0))
    flat_v = volume.ravel().copy()
    flat_p = power.ravel().copy()
    nx, ny, nz = shape
    moved = 0.0
    tolerance = max(cell_volume * 1.0e-12, 1.0e-12)
    for donor in np.flatnonzero(flat_v > cell_volume + tolerance):
        excess = flat_v[donor] - cell_volume
        source_per_volume = flat_p[donor] / flat_v[donor]
        flat_v[donor] = cell_volume
        flat_p[donor] = cell_volume * source_per_volume
        ix, iy, iz = donor // (ny * nz), (donor // nz) % ny, donor % nz
        visited = {int(donor)}
        for radius in range(1, max(shape) // 2 + 1):
            offsets = [
                (di, dj, dk)
                for di, dj, dk in itertools.product(range(-radius, radius + 1), repeat=3)
                if max(abs(di), abs(dj), abs(dk)) == radius
            ]
            offsets.sort(
                key=lambda d: (d[0] * spacing_A[0]) ** 2
                + (d[1] * spacing_A[1]) ** 2
                + (d[2] * spacing_A[2]) ** 2
            )
            for di, dj, dk in offsets:
                z = iz + dk
                if z < 0 or z >= nz:
                    continue
                receiver = int(np.ravel_multi_index(((ix + di) % nx, (iy + dj) % ny, z), shape))
                if receiver in visited or flat_v[receiver] <= 0.0:
                    continue
                visited.add(receiver)
                capacity = cell_volume - flat_v[receiver]
                if capacity <= tolerance:
                    continue
                amount = min(excess, capacity)
                flat_v[receiver] += amount
                flat_p[receiver] += amount * source_per_volume
                excess -= amount
                moved += amount
                if excess <= tolerance:
                    break
            if excess <= tolerance:
                break
        if excess > tolerance:
            raise RuntimeError(f"Could not redistribute Ag overfill ({excess:.6g} A^3) at cell {donor}")
    return flat_v.reshape(shape), flat_p.reshape(shape), raw_overfilled, raw_max_fill, moved


def map_ag_fields(
    xyz_A: np.ndarray,
    q_density: np.ndarray,
    shape: tuple[int, int, int],
    spacing_A: np.ndarray,
) -> dict[str, np.ndarray | int | float]:
    ids, weights = cloud_in_cell_map(xyz_A, shape, spacing_A)
    ncell = int(np.prod(shape))
    cell_volume = float(np.prod(spacing_A))
    atom_source = AG_ATOMIC_VOLUME_A3 * q_density
    raw_volume = np.bincount(
        ids.ravel(),
        weights=np.repeat(np.full(len(xyz_A), AG_ATOMIC_VOLUME_A3), ids.shape[1]) * weights.ravel(),
        minlength=ncell,
    )
    raw_power = np.bincount(
        ids.ravel(),
        weights=np.repeat(atom_source, ids.shape[1]) * weights.ravel(),
        minlength=ncell,
    )
    volume, power, overfilled, max_fill, moved = _redistribute_overfill(
        raw_volume, raw_power, shape, spacing_A
    )
    expected_volume = len(xyz_A) * AG_ATOMIC_VOLUME_A3
    expected_power = float(atom_source.sum())
    if not np.isclose(volume.sum(), expected_volume, rtol=1.0e-12, atol=1.0e-8):
        raise RuntimeError("Ag volume mapping is not conservative")
    if not np.isclose(power.sum(), expected_power, rtol=1.0e-12, atol=1.0e-10):
        raise RuntimeError("Ag source mapping is not conservative")
    fraction = volume / cell_volume
    if fraction.max(initial=0.0) > 1.0 + 1.0e-12:
        raise RuntimeError("Ag fill fraction exceeds one after redistribution")
    return {
        "fraction_ag": fraction,
        "source_power_eV_ps": power,
        "raw_overfilled_cells": overfilled,
        "raw_max_fill_fraction": max_fill,
        "redistributed_ag_volume_A3": moved,
        "atom_quadrature_power_eV_ps": expected_power,
    }


def grid_for_spacing(spacing_nominal_A: float) -> tuple[tuple[int, int, int], np.ndarray]:
    lengths = np.array([BOX_A[0], BOX_A[1], SOURCE_SURFACE_Z_A], dtype=float)
    shape_arr = np.ceil(lengths / spacing_nominal_A).astype(int)
    shape = tuple(int(value) for value in shape_arr)
    spacing = lengths / shape_arr
    return shape, spacing


def continuum_power_eV_ps() -> float:
    xhalf = BOX_A[0] / 2.0
    yhalf = BOX_A[1] / 2.0
    lateral_x = R_A * math.sqrt(math.pi) * math.erf(xhalf / R_A)
    lateral_y = R_A * math.sqrt(math.pi) * math.erf(yhalf / R_A)
    depth_integral = SKIN_DEPTH_A * (1.0 - math.exp(-SOURCE_SURFACE_Z_A / SKIN_DEPTH_A))
    return I0_EV_PS_A2 / SKIN_DEPTH_A * lateral_x * lateral_y * depth_integral


def integrate_envelope(path: Path) -> tuple[float, float, int]:
    times: list[float] = []
    values: list[float] = []
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not {"time_ps", "S"}.issubset(reader.fieldnames or []):
            raise ValueError("Envelope CSV must have columns time_ps,S")
        for row in reader:
            times.append(float(row["time_ps"]))
            values.append(float(row["S"]))
    t = np.asarray(times, dtype=float)
    s = np.asarray(values, dtype=float)
    if len(t) < 2 or not np.all(np.isfinite(t)) or not np.all(np.isfinite(s)):
        raise ValueError("Envelope must contain at least two finite rows")
    if t[0] < 0.0 or np.any(np.diff(t) <= 0.0) or np.any(s < 0.0):
        raise ValueError("time_ps must increase from a nonnegative time and S must be nonnegative")
    return float(np.trapezoid(s, t)), float(t[-1] - t[0]), len(t)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--envelope-csv", type=Path, default=None, help="optional CSV columns time_ps,S")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    atom_types, xyz, lengths = read_lammps_atomic_data(args.data)
    if not np.allclose(lengths, BOX_A, rtol=0.0, atol=1.0e-8):
        raise ValueError(f"Expected box {BOX_A.tolist()} A, found {lengths.tolist()} A")
    ag_mask = atom_types == 2
    w_mask = atom_types == 1
    in_material = (xyz[:, 2] >= MATERIAL_Z_A[0]) & (xyz[:, 2] <= MATERIAL_Z_A[1])
    newton = (xyz[:, 2] >= NEWTON_Z_A[0]) & (xyz[:, 2] <= NEWTON_Z_A[1])
    ag_all = ag_mask & in_material
    ag_newton = ag_all & newton
    if not np.any(ag_all) or not np.any(ag_newton) or np.any(ag_mask & ~in_material):
        raise ValueError("Ag atoms are absent from the material/Newton region or lie outside the defined material")

    q_all = source_density_eV_ps_A3(xyz[ag_all])
    q_newton = source_density_eV_ps_A3(xyz[ag_newton])
    p_ag_all_atom = float(AG_ATOMIC_VOLUME_A3 * q_all.sum())
    p_ag_newton_atom = float(AG_ATOMIC_VOLUME_A3 * q_newton.sum())
    p_continuum = continuum_power_eV_ps()
    grid_reports: list[dict[str, object]] = []
    saved_fields: dict[str, np.ndarray] | None = None
    all_mask_stats: dict[str, float | int] = {}
    for nominal in GRID_SPACINGS_A:
        shape, spacing = grid_for_spacing(nominal)
        mapped_newton = map_ag_fields(xyz[ag_newton], q_newton, shape, spacing)
        p_grid = float(np.asarray(mapped_newton["source_power_eV_ps"]).sum())
        relative_error = p_grid / p_ag_newton_atom - 1.0
        grid_reports.append(
            {
                "nominal_spacing_A": nominal,
                "shape_xyz": list(shape),
                "actual_spacing_xyz_A": spacing.tolist(),
                "nonzero_Ag_cells": int(np.count_nonzero(mapped_newton["fraction_ag"])),
                "raw_overfilled_Ag_cells": int(mapped_newton["raw_overfilled_cells"]),
                "raw_max_Ag_fill_fraction": float(mapped_newton["raw_max_fill_fraction"]),
                "redistributed_Ag_volume_A3": float(mapped_newton["redistributed_ag_volume_A3"]),
                "Newton_Ag_atom_quadrature_power_eV_ps": p_ag_newton_atom,
                "Newton_Ag_grid_power_eV_ps": p_grid,
                "grid_vs_atom_power_relative_error": relative_error,
            }
        )
        if not np.isclose(p_grid, p_ag_newton_atom, rtol=1.0e-12, atol=1.0e-9):
            raise RuntimeError(f"{nominal:g} A grid does not conserve mapped Ag source power")
        if nominal == OUTPUT_GRID_A:
            mapped_all = map_ag_fields(xyz[ag_all], q_all, shape, spacing)
            saved_fields = {
                "fraction_ag_all_material": np.asarray(mapped_all["fraction_ag"], dtype=np.float32),
                "fraction_ag_newton_only": np.asarray(mapped_newton["fraction_ag"], dtype=np.float32),
                "source_power_newton_ag_eV_ps": np.asarray(mapped_newton["source_power_eV_ps"], dtype=np.float64),
                "cell_spacing_A": spacing.astype(np.float64),
                "grid_shape_xyz": np.asarray(shape, dtype=np.int32),
            }
            all_grid_power = float(np.asarray(mapped_all["source_power_eV_ps"]).sum())
            all_grid_volume = float(np.asarray(mapped_all["fraction_ag"]).sum() * np.prod(spacing))
            all_mask_stats = {
                "raw_overfilled_cells": int(mapped_all["raw_overfilled_cells"]),
                "raw_max_fill_fraction": float(mapped_all["raw_max_fill_fraction"]),
                "redistributed_Ag_volume_A3": float(mapped_all["redistributed_ag_volume_A3"]),
            }

    if saved_fields is None:
        raise RuntimeError("Requested output grid not found in grid candidates")
    np.savez_compressed(args.out / "ag_mask_source_10A.npz", **saved_fields)

    envelope = None
    energy_on_source_eV = None
    if args.envelope_csv is not None:
        integral_s_ps, envelope_span_ps, sample_count = integrate_envelope(args.envelope_csv)
        energy_on_source_eV = p_ag_newton_atom * integral_s_ps
        envelope = {
            "file": str(args.envelope_csv),
            "time_span_ps": envelope_span_ps,
            "sample_count": sample_count,
            "integral_S_dt_ps": integral_s_ps,
        }
    try:
        input_data_ref = str(args.data.resolve().relative_to(ROOT))
    except ValueError:
        input_data_ref = str(args.data.resolve())

    report = {
        "status": "STATIC_ACTUAL_GEOMETRY_SOURCE_PREFLIGHT_ONLY",
        "input_data": input_data_ref,
        "source_reference": "Migrated spatial source form/parameters from models/ag_ti3sic2_qsgs_1to1_finer/TTM_SOURCE_STATIC_PREFLIGHT.md; legacy project files remain unchanged.",
        "equation": "q_e = f_Ag(x,y,z) * (I0/l_skin) * exp(-r_xy^2/R^2) * exp(-(z_surface-z)/l_skin) * S(t)",
        "source_parameters": {
            "I0_eV_ps_A-2": I0_EV_PS_A2,
            "R_1e_A": R_A,
            "l_skin_A": SKIN_DEPTH_A,
            "center_xy_A": SOURCE_CENTER_XY_A.tolist(),
            "surface_z_A": SOURCE_SURFACE_Z_A,
            "q0_at_center_surface_eV_ps_A-3_when_S1": I0_EV_PS_A2 / SKIN_DEPTH_A,
            "S_t": "unspecified; all reported on-state powers use S(t)=1",
        },
        "geometry": {
            "box_A": lengths.tolist(),
            "boundary": "p p f",
            "material_z_A": list(MATERIAL_Z_A),
            "Newton_source_exchange_z_A": list(NEWTON_Z_A),
            "atom_counts": {"total": int(len(atom_types)), "Ag": int(ag_mask.sum()), "W": int(w_mask.sum()), "Newton_Ag": int(ag_newton.sum())},
            "mask_definition": "CIC mapping of actual type-2 Ag atom coordinates with FCC atomic-volume proxy; local cell fill is bounded to [0,1]. It is not the global Ag volume fraction.",
            "Ag_atomic_volume_proxy_A3": AG_ATOMIC_VOLUME_A3,
            "W_atomic_volume_proxy_A3": W_ATOMIC_VOLUME_A3,
            "global_Ag_volume_fraction_target_not_used_as_mask": 0.6473,
        },
        "on_state_power_S1": {
            "unmasked_continuum_over_material_box_eV_ps": p_continuum,
            "actual_all_material_Ag_atom_quadrature_eV_ps": p_ag_all_atom,
            "actual_Newton_layer_Ag_atom_quadrature_eV_ps": p_ag_newton_atom,
            "actual_Newton_layer_Ag_grid_10A_eV_ps": grid_reports[1]["Newton_Ag_grid_power_eV_ps"],
            "all_material_Ag_grid_10A_eV_ps": all_grid_power,
            "W_direct_source_eV_ps": 0.0,
            "all_material_Ag_mapped_volume_A3": all_grid_volume,
            "all_material_Ag_fraction_of_material_volume": all_grid_volume / (BOX_A[0] * BOX_A[1] * SOURCE_SURFACE_Z_A),
            "all_material_Ag_mask_10A_stats": all_mask_stats,
            "Ag_share_of_unmasked_continuum_percent": 100.0 * p_ag_all_atom / p_continuum,
            "Newton_Ag_share_of_unmasked_continuum_percent": 100.0 * p_ag_newton_atom / p_continuum,
        },
        "cumulative_input_energy": {
            "time_envelope": envelope,
            "value_eV": energy_on_source_eV,
            "formula": "E_Ag(t)=P_Newton_Ag,on * integral_0^t S(t') dt'",
            "energy_per_1ps_unit_S_integral_eV": p_ag_newton_atom,
            "interpretation": "No numerical pulse energy is assigned without a specified S(t); with this source the coefficient per ps of unit S-integral equals P_Newton_Ag,on.",
        },
        "grid_sensitivity": grid_reports,
        "limitations": [
            "This maps only the spatial source and Ag mask. It does not run or validate TTM/MD.",
            "Only Newton-layer Ag is given a direct source, matching the source-only convention in the Ag-Ti3SiC2 reference; W receives zero direct electronic source in this exploratory assumption.",
            "The QSGS raw geometry has failed the EFS minimization preflight and is not an MD-ready structure; this report is geometry-only.",
            "The fixed I0, R and skin depth remain exploratory candidates, not Ag-W optical/absorption measurements. No pulse width or temporal envelope was specified.",
        ],
    }
    (args.out / "source_preflight_report.json").write_text(json.dumps(report, indent=2) + "\n")
    md = [
        "# Ag–W TTM source preflight",
        "",
        "This is a static mapping of the migrated Ag-only candidate source onto the actual Ag–W QSGS atom coordinates. It does not run TTM or MD.",
        "",
        "## Source and mask",
        "",
        "```text",
        report["equation"],
        "```",
        "",
        f"Candidate values retained from the Ag–Ti₃SiC₂ source note: `I0={I0_EV_PS_A2:g} eV/(ps·Å²)`, `R={R_A:g} Å`, `l_skin={SKIN_DEPTH_A:g} Å`, center=({SOURCE_CENTER_XY_A[0]:g},{SOURCE_CENTER_XY_A[1]:g}) Å, z_surface={SOURCE_SURFACE_Z_A:g} Å.",
        "",
        "Here `r_xy` is the x/y periodic minimum-image distance from the beam center and `d=z_surface−z` is depth into the material. The numerical powers below use `S(t)=1` and the actual Ag atom locations.",
        "",
        "The Ag mask is mapped from actual LAMMPS atom type 2 coordinates, with an FCC atomic-volume proxy and cloud-in-cell deposition onto a 10 Å grid. Cell occupancy is locally bounded to [0,1]; it is not replaced by the global 0.6473 Ag volume fraction. The source/exchange convention follows the legacy reference: only Ag atoms in the Newton layer (40–200 Å) receive direct source. W has zero direct source in this candidate.",
        "",
        "## On-state power, S(t)=1",
        "",
        "| Quantity | Power (eV/ps) |",
        "|---|---:|",
        f"| Unmasked continuum over the material box | {p_continuum:.6f} |",
        f"| Actual Ag atom quadrature, full material | {p_ag_all_atom:.6f} |",
        f"| Actual Ag atom quadrature, Newton layer | {p_ag_newton_atom:.6f} |",
        f"| Newton Ag 10 Å grid sum | {grid_reports[1]['Newton_Ag_grid_power_eV_ps']:.6f} |",
        "| W direct source | 0 |",
        "",
        f"The Newton-layer Ag mask captures {100.0 * p_ag_newton_atom / p_continuum:.4f}% of the unmasked spatial source power for this actual morphology. The 10 Å grid sum agrees with atom quadrature by {abs(float(grid_reports[1]['grid_vs_atom_power_relative_error'])):.3e} relative error.",
        "",
        "## Grid sensitivity",
        "",
        "| Nominal spacing (Å) | Grid | Raw overfilled Ag cells | Raw max Ag fill | Volume redistributed (Å³) | Newton Ag power (eV/ps) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for item in grid_reports:
        dims = "×".join(str(x) for x in item["shape_xyz"])
        md.append(
            f"| {item['nominal_spacing_A']:g} | {dims} | {item['raw_overfilled_Ag_cells']} | {item['raw_max_Ag_fill_fraction']:.6f} | {item['redistributed_Ag_volume_A3']:.3f} | {item['Newton_Ag_grid_power_eV_ps']:.6f} |"
        )
    if envelope is None:
        energy_summary = (
            f"Cumulative input energy is `E_Ag = {p_ag_newton_atom:.6f} eV/ps × integral S(t) dt`. "
            f"No pulse width or temporal envelope is specified, so the actual cumulative energy is **undetermined**. "
            f"The coefficient per 1 ps of unit `S` integral is {p_ag_newton_atom:.6f} eV; this is a normalization reference, not a selected pulse."
        )
    else:
        energy_summary = (
            f"For the supplied envelope `{args.envelope_csv}`, `integral S(t) dt={envelope['integral_S_dt_ps']:.6f} ps` "
            f"over {envelope['time_span_ps']:.6f} ps; mapped cumulative Newton-layer Ag source energy is "
            f"{energy_on_source_eV:.6f} eV. This follows the supplied envelope and is not an independently validated laser pulse."
        )
    md.extend(
        [
            "",
            "## Time integral and interpretation",
            "",
            energy_summary,
            "",
            "This assumption heats Ag electrons directly; W can gain lattice energy only through atomistic Ag–W interactions in a later MD model. Electronic transport and electron–lattice coupling in W or across Ag/W are not represented. The Ag/W property closure remains unverified, and the current QSGS geometry failed its static EFS initialization preflight. Do not interpret this report as an Ag–W laser simulation or as a production heating prescription.",
            "",
            f"Machine-readable values and the actual 10 Å mask/source arrays are in `source_preflight_report.json` and `ag_mask_source_10A.npz`. Source reference: `{report['source_reference']}`",
            "",
            "Reproduce from the repository root with `python research/ag_w_reference_pack/scripts/ag_w_ttm_source_preflight.py`. An optional `--envelope-csv` accepts `time_ps,S` columns to integrate a specified temporal envelope.",
            "",
        ]
    )
    (args.out / "TTM_SOURCE_PREFLIGHT.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps({"output": str(args.out), "status": report["status"], "power": report["on_state_power_S1"], "energy": report["cumulative_input_energy"], "grid_sensitivity": grid_reports}, indent=2))


if __name__ == "__main__":
    main()
