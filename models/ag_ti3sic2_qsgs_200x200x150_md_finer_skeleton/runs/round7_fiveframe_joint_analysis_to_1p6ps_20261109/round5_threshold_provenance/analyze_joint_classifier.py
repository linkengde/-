#!/usr/bin/env python3
"""Calibrate Q6/PTM/1-ps MSD jointly from same-MTP pure-Ag references."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier, SteinhardtParametersModifier

HERE = Path(__file__).resolve().parent
FRAME_INTERVAL_PS = 0.05
WINDOW_PS = 1.0
Q6_CUTOFF_A = 4.0
FCC_CODE = 1


def read_dump(path: Path):
    frames = []
    with path.open() as stream:
        while True:
            line = stream.readline()
            if not line:
                break
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Malformed dump: {path}")
            step = int(stream.readline())
            if stream.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"Missing atom count at step {step}")
            count = int(stream.readline())
            bounds_header = stream.readline().split()
            if bounds_header[:3] != ["ITEM:", "BOX", "BOUNDS"] or "xy" in bounds_header:
                raise ValueError(f"Expected orthogonal box in {path}")
            bounds = np.array([[float(v) for v in stream.readline().split()[:2]] for _ in range(3)])
            origin = bounds[:, 0]
            box = bounds[:, 1] - origin
            columns = stream.readline().split()[2:]
            col = {name: i for i, name in enumerate(columns)}
            data = np.loadtxt((stream.readline() for _ in range(count)), ndmin=2)
            data = data[np.argsort(data[:, col["id"]])]
            ids = data[:, col["id"]].astype(np.int64)
            xyz = data[:, [col[axis] for axis in ("x", "y", "z")]]
            frames.append({"step": step, "time_ps": step * 0.0001, "origin": origin, "box": box, "ids": ids, "xyz": xyz})
    return frames


def per_frame_q6_ptm(path: Path):
    pipeline = import_file(str(path))
    pipeline.modifiers.append(SelectTypeModifier(types={4}))
    q6 = SteinhardtParametersModifier()
    q6.max_l = 6
    q6.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    q6.cutoff = Q6_CUTOFF_A
    q6.only_selected = True
    pipeline.modifiers.append(q6)
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    result = []
    for frame_index in range(pipeline.source.num_frames):
        data = pipeline.compute(frame_index)
        ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
        order = np.argsort(ids)
        type_id = np.asarray(data.particles["Particle Type"], dtype=np.int32)[order]
        order_parameters = np.asarray(data.particles["Steinhardt"], dtype=float)[order]
        if order_parameters.ndim != 2 or order_parameters.shape[1] < 5:
            raise RuntimeError(f"Unexpected Steinhardt output shape {order_parameters.shape}")
        q6_values = order_parameters[:, 4]
        structure = np.asarray(data.particles["Structure Type"], dtype=np.int32)[order]
        result.append({"ids": ids[order], "types": type_id, "q6": q6_values, "structure": structure})
    return result


def drift_corrected_window_displacements(frames, window_ps=WINDOW_PS):
    lag = int(round(window_ps / FRAME_INTERVAL_PS))
    if len(frames) <= lag:
        raise ValueError(f"Need more than {lag} frames for a {window_ps:g} ps MSD window; got {len(frames)}")
    all_displacements = []
    for end in range(lag, len(frames)):
        start = end - lag
        increments = []
        for k in range(start + 1, end + 1):
            prev, curr = frames[k - 1], frames[k]
            if not np.array_equal(prev["ids"], curr["ids"]):
                raise ValueError("Atom IDs changed across reference frames")
            frac_prev = (prev["xyz"] - prev["origin"]) / prev["box"]
            frac_curr = (curr["xyz"] - curr["origin"]) / curr["box"]
            delta_frac = frac_curr - frac_prev
            delta_frac -= np.rint(delta_frac)
            delta = delta_frac * (0.5 * (prev["box"] + curr["box"]))
            increments.append(delta)
        net = np.sum(increments, axis=0)
        net -= np.mean(net, axis=0)
        all_displacements.append({"end_frame": end, "time_ps": frames[end]["time_ps"], "magnitude_A": np.linalg.norm(net, axis=1)})
    return all_displacements


def summarize_reference(name: str, dump_path: Path):
    frames = read_dump(dump_path)
    structural = per_frame_q6_ptm(dump_path)
    if len(frames) != len(structural):
        raise RuntimeError("Dump parser and OVITO frame counts disagree")
    if any(not np.array_equal(f["ids"], s["ids"]) for f, s in zip(frames, structural)):
        raise RuntimeError("Dump and OVITO atom ordering mismatch")
    displacement = drift_corrected_window_displacements(frames)
    final_start = max(0, len(frames) - 21)
    q6_window = np.mean(np.stack([x["q6"] for x in structural[final_start:]]), axis=0)
    nonfcc_window = np.mean(np.stack([x["structure"] != FCC_CODE for x in structural[final_start:]]), axis=0)
    final_disp = displacement[-1]["magnitude_A"]
    return {
        "name": name,
        "dump": str(dump_path),
        "frames": len(frames),
        "time_first_ps": frames[0]["time_ps"],
        "time_last_ps": frames[-1]["time_ps"],
        "window_start_ps": displacement[-1]["time_ps"] - WINDOW_PS,
        "window_end_ps": displacement[-1]["time_ps"],
        "atom_count": len(frames[-1]["ids"]),
        "mean_q6_over_last_1ps": float(np.mean(q6_window)),
        "q6_quantiles": np.percentile(q6_window, [5, 25, 50, 75, 95]).tolist(),
        "nonfcc_fraction_over_last_1ps_mean": float(np.mean(nonfcc_window)),
        "nonfcc_fraction_quantiles": np.percentile(nonfcc_window, [5, 25, 50, 75, 95]).tolist(),
        "fcc_fraction_final": float(np.mean(structural[-1]["structure"] == FCC_CODE)),
        "drift_corrected_msd_window_displacement_A_quantiles": np.percentile(final_disp, [5, 50, 90, 95, 99]).tolist(),
        "msd_window_time_ps": displacement[-1]["time_ps"],
        "box_final_A": frames[-1]["box"].tolist(),
        "q6_cutoff_A": Q6_CUTOFF_A,
    }


def main():
    solid = summarize_reference("round5 MTP Ag FCC/NVT 300 K", HERE / "solid_300K.dump")
    liquid = summarize_reference("round5 MTP Ag candidate/NPT 2600 K", HERE / "liquid_2600K.dump")
    dense_liquid = summarize_reference("round5 MTP Ag candidate/NVT 2600 K, compressed reference density", HERE / "liquid_NVT_density_matched.dump")
    result = {
        "status": "REFERENCE_METRICS_COMPUTED; CLASSIFIER_SEPARATION_NOT_YET_DECIDED",
        "definitions": {
            "q6": "Ag-only Steinhardt q6, 4.0 A cutoff, averaged per atom over final 1 ps",
            "ptm": "FCC code 1; non-FCC fraction per atom over final 1 ps",
            "msd": "1 ps unwrapped displacement from periodic minimum-image increments; affine NPT expansion and whole-cell drift removed",
        },
        "solid": solid,
        "liquid_candidate": liquid,
        "dense_liquid_candidate": dense_liquid,
    }
    (HERE / "reference_metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
