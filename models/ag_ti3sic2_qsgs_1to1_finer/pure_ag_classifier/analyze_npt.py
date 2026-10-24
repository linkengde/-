#!/usr/bin/env python3
"""Analyze pure-Ag NPT trajectory: thermo, density, MSD, and liquid RDF."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from analyze_dumps import read_dump

AMU_A3_TO_G_CM3 = 1.66053906660


def parse_thermo(log_path: Path):
    rows = []
    active = False
    with log_path.open() as stream:
        for line in stream:
            fields = line.split()
            if fields[:2] == ["Step", "Time"]:
                active = True
                continue
            if not active:
                continue
            if len(fields) == 7:
                try:
                    rows.append([float(value) for value in fields])
                    continue
                except ValueError:
                    pass
            if rows and not fields:
                active = False
    if not rows:
        raise ValueError(f"Could not find thermo table in {log_path}")
    return np.asarray(rows)


def analyze_dump(dump_path: Path, first_rdf_time_ps: float, mass_amu: float,
                 r_max_A: float, bin_width_A: float):
    frames = list(read_dump(dump_path))
    series = []
    rdf_hist = None
    rdf_norm = None
    rdf_frame_count = 0
    n_atoms = None
    edges = np.arange(0.0, r_max_A + bin_width_A, bin_width_A)
    centers = 0.5 * (edges[:-1] + edges[1:])

    unwrapped_displacement_A = None
    previous_fractional = None
    previous_box = None
    for step, origin, box, data, col in frames:
        n_atoms = len(data)
        time_ps = step * 0.0005
        volume = float(np.prod(box))
        density = n_atoms * mass_amu / volume * AMU_A3_TO_G_CM3
        cna_unknown = float(np.mean(data[:, col["c_cna"]] == 5))
        # Integrate fractional displacements with minimum images and the mean
        # box over each interval. This removes homogeneous barostat strain.
        xyz = data[:, [col[axis] for axis in ("x", "y", "z")]]
        fractional = ((xyz - origin) / box) % 1.0
        if previous_fractional is None:
            unwrapped_displacement_A = np.zeros_like(fractional)
        else:
            delta_fractional = fractional - previous_fractional
            delta_fractional -= np.rint(delta_fractional)
            unwrapped_displacement_A += delta_fractional * (0.5 * (box + previous_box))
        msd = float(np.mean(np.sum(unwrapped_displacement_A ** 2, axis=1)))
        series.append((step, time_ps, volume, density, cna_unknown, msd))
        previous_fractional = fractional
        previous_box = box

        if time_ps + 1.0e-12 < first_rdf_time_ps:
            continue
        xyz = data[:, [col[axis] for axis in ("x", "y", "z")]] - origin
        xyz %= box
        pairs = cKDTree(xyz, boxsize=box).query_pairs(r_max_A, output_type="ndarray")
        if pairs.size:
            delta = xyz[pairs[:, 0]] - xyz[pairs[:, 1]]
            delta -= np.rint(delta / box) * box
            distance = np.linalg.norm(delta, axis=1)
            counts, _ = np.histogram(distance, bins=edges)
        else:
            counts = np.zeros(len(centers), dtype=np.int64)
        shell = (4.0 * np.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)
        number_density = n_atoms / volume
        expected = 0.5 * n_atoms * number_density * shell
        g_r = counts / expected
        rdf_hist = g_r if rdf_hist is None else rdf_hist + g_r
        rdf_frame_count += 1

    series = np.asarray(series, dtype=float)
    rdf = rdf_hist / rdf_frame_count if rdf_frame_count else None
    return series, centers, rdf, rdf_frame_count, n_atoms


def write_csv(path: Path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dump", type=Path)
    parser.add_argument("log", type=Path)
    parser.add_argument("--prefix", type=Path, default=Path("analysis_1600K_npt"))
    parser.add_argument("--rdf-start-ps", type=float, default=5.0)
    parser.add_argument("--mass-amu", type=float, default=107.8682)
    parser.add_argument("--rdf-max-A", type=float, default=10.0)
    parser.add_argument("--rdf-bin-A", type=float, default=0.05)
    args = parser.parse_args()

    series, radii, rdf, rdf_frames, n_atoms = analyze_dump(
        args.dump, args.rdf_start_ps, args.mass_amu, args.rdf_max_A, args.rdf_bin_A
    )
    thermo = parse_thermo(args.log)
    thermo_by_step = {int(row[0]): row[1:] for row in thermo}
    merged = []
    for step, time_ps, volume, density, unknown, msd in series:
        rest = thermo_by_step.get(int(step), [np.nan] * 6)
        # The logged fields after Step are Time, Temp, PotEng, KinEng,
        # TotEng, and Press.
        _, temp, pe, ke, etotal, press = rest
        merged.append((int(step), time_ps, temp, pe, ke, etotal, press,
                       volume, density, unknown, msd))
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    series_path = args.prefix.with_name(args.prefix.name + "_series.csv")
    rdf_path = args.prefix.with_name(args.prefix.name + "_rdf.csv")
    write_csv(
        series_path,
        ["step", "time_ps", "temperature_K", "potential_energy_eV", "kinetic_energy_eV",
         "total_energy_eV", "pressure_bar", "volume_A3", "density_g_cm3",
         "CNA_unknown_fraction", "MSD_A2"], merged,
    )
    if rdf is not None:
        write_csv(rdf_path, ["r_A", "g_r"], zip(radii, rdf))

    t = series[:, 1]
    fit_mask = (t >= 2.5) & (t <= t[-1])
    slope, intercept = np.polyfit(t[fit_mask], series[fit_mask, 5], 1)
    tail = (t >= max(2.5, t[-1] - 5.0))
    tail_temp = np.array([row[2] for row in merged])[tail]
    tail_density = series[tail, 3]
    tail_volume = series[tail, 2]
    tail_pressure = np.array([row[6] for row in merged])[tail]
    tail_pe = np.array([row[3] for row in merged])[tail]
    tail_etotal = np.array([row[5] for row in merged])[tail]
    fit_msd = series[fit_mask, 5]
    msd_fit = slope * t[fit_mask] + intercept
    r_squared = 1.0 - np.sum((fit_msd - msd_fit) ** 2) / np.sum((fit_msd - fit_msd.mean()) ** 2)
    print(f"atoms={n_atoms}; frames={len(series)}; RDF frames={rdf_frames}")
    print(f"MSD fit 2.5–{t[-1]:g} ps: slope={slope:.6g} A^2/ps, D=slope/6={slope/6:.6g} A^2/ps, R2={r_squared:.5f}")
    print(f"last-5-ps temperature mean/std={tail_temp.mean():.3f}/{tail_temp.std(ddof=1):.3f} K")
    print(f"last-5-ps density mean/std={tail_density.mean():.6f}/{tail_density.std(ddof=1):.6f} g/cm^3")
    print(f"last-5-ps volume mean/std={tail_volume.mean():.3f}/{tail_volume.std(ddof=1):.3f} A^3")
    print(f"last-5-ps pressure mean/std={tail_pressure.mean():.1f}/{tail_pressure.std(ddof=1):.1f} bar")
    print(f"last-5-ps PE mean/std={tail_pe.mean():.3f}/{tail_pe.std(ddof=1):.3f} eV; Etot mean/std={tail_etotal.mean():.3f}/{tail_etotal.std(ddof=1):.3f} eV")
    print(f"last frame: T={merged[-1][2]:.3f} K, CNA-unknown={merged[-1][9]:.6f}, MSD={merged[-1][10]:.6f} A^2")
    if rdf is not None:
        peak = (radii >= 2.5) & (radii <= 3.5)
        peak_i = np.flatnonzero(peak)[np.argmax(rdf[peak])]
        minimum = (radii > radii[peak_i] + 0.15) & (radii <= 5.2)
        min_i = np.flatnonzero(minimum)[np.argmin(rdf[minimum])]
        print(f"RDF({args.rdf_start_ps:g}–{t[-1]:g} ps): first peak r={radii[peak_i]:.3f} A, g={rdf[peak_i]:.4f}; first minimum r={radii[min_i]:.3f} A, g={rdf[min_i]:.4f}")
    print(f"wrote {series_path} and {rdf_path if rdf is not None else '(no RDF)'}")


if __name__ == "__main__":
    main()
