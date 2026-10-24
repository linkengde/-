#!/usr/bin/env python3
"""Report one-ps CNA/displacement classifier statistics from LAMMPS dumps."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read_dump(path: Path):
    with path.open() as stream:
        while True:
            line = stream.readline()
            if not line:
                return
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Malformed dump before timestep in {path}")
            step = int(stream.readline())
            if stream.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"Missing atom count in {path} at step {step}")
            count = int(stream.readline())
            bounds_header = stream.readline().split()
            if bounds_header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
                raise ValueError(f"Missing box bounds in {path} at step {step}")
            bounds = np.array([[float(v) for v in stream.readline().split()[:2]] for _ in range(3)])
            origin = bounds[:, 0]
            box = bounds[:, 1] - origin
            atom_header = stream.readline().split()[2:]
            columns = {name: index for index, name in enumerate(atom_header)}
            needed = {"id", "x", "y", "z", "c_cna"}
            if not needed.issubset(columns):
                raise ValueError(f"Missing columns {needed - columns.keys()} in {path}")
            data = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
            data = data[np.argsort(data[:, columns["id"]])]
            yield step, origin, box, data, columns


def summarize(path: Path, timestep_ps: float, frame_interval_steps: int, window_ps: float):
    frames = list(read_dump(path))
    lag_frames = int(round(window_ps * 1000.0 / (frame_interval_steps * timestep_ps * 1000.0)))
    if lag_frames < 1:
        raise ValueError("Window must be at least one dump interval")
    print(f"\n{path}: {len(frames)} frames; displacement window={window_ps:g} ps")
    print("time_ps unknown_fraction disp_p50_A disp_p90_A disp_p99_A liquid_flag_fraction")
    for index in range(lag_frames, len(frames)):
        step, origin, box, data, columns = frames[index]
        old_step, old_origin, old_box, old, old_columns = frames[index - lag_frames]
        if columns["id"] != old_columns["id"] or not np.array_equal(
            data[:, columns["id"]], old[:, old_columns["id"]]
        ):
            raise ValueError("Atom ID ordering changed")
        xyz = data[:, [columns[axis] for axis in ("x", "y", "z")]]
        old_xyz = old[:, [old_columns[axis] for axis in ("x", "y", "z")]]
        frac = (xyz - origin) / box
        old_frac = (old_xyz - old_origin) / old_box
        delta_frac = frac - old_frac
        delta_frac -= np.rint(delta_frac)
        delta = delta_frac * (0.5 * (box + old_box))
        displacement = np.linalg.norm(delta, axis=1)
        cna = data[:, columns["c_cna"]]
        liquid = (cna == 5) & (displacement >= 1.0)
        time_ps = step * timestep_ps
        p50, p90, p99 = np.percentile(displacement, [50, 90, 99])
        print(
            f"{time_ps:7.3f} {np.mean(cna == 5):.6f} {p50:.6f} {p90:.6f} "
            f"{p99:.6f} {np.mean(liquid):.6f} (lag from step {old_step})"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dumps", nargs="+", type=Path)
    parser.add_argument("--timestep-ps", type=float, default=0.0005)
    parser.add_argument("--frame-interval-steps", type=int, default=1000)
    parser.add_argument("--window-ps", type=float, default=1.0)
    args = parser.parse_args()
    for dump in args.dumps:
        summarize(dump, args.timestep_ps, args.frame_interval_steps, args.window_ps)


if __name__ == "__main__":
    main()
