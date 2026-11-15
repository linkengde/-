#!/usr/bin/env python3
"""Generate a reproducible 150^3 random four-parameter QSGS geometry preview.

This creates geometry only; it does not run the thermal or flow solver.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "preview"


def shifted(mask: np.ndarray, axis: int, direction: int) -> np.ndarray:
    """Shift without wrapping across a physical boundary."""
    out = np.roll(mask, direction, axis=axis)
    edge = [slice(None)] * mask.ndim
    edge[axis] = 0 if direction == 1 else -1
    out[tuple(edge)] = False
    return out


def qsgs_four_parameter(
    shape: tuple[int, int, int],
    p0: float,
    growth: tuple[float, float, float],
    target_fraction: float,
    seed: int,
) -> tuple[np.ndarray, int]:
    """Grow a connected ceramic phase from random nuclei using P0/Px/Py/Pz."""
    if any(n <= 0 for n in shape):
        raise ValueError("All grid dimensions must be positive")
    if not 0.0 < p0 < 1.0 or any(not 0.0 < p < 1.0 for p in growth):
        raise ValueError("P0 and Px/Py/Pz must lie strictly between zero and one")
    if not 0.0 < target_fraction < 1.0:
        raise ValueError("target skeleton fraction must lie strictly between zero and one")

    rng = np.random.default_rng(seed)
    occupied = rng.random(shape) < p0
    target_count = int(round(target_fraction * occupied.size))
    if not occupied.any():
        occupied.flat[int(rng.integers(occupied.size))] = True

    generations = 0
    while int(occupied.sum()) < target_count:
        frontier_accept = np.zeros(shape, dtype=bool)
        for axis, probability in enumerate(growth):
            for direction in (-1, 1):
                neighbor = shifted(occupied, axis, direction)
                frontier_accept |= neighbor & (rng.random(shape) < probability)

        frontier = np.flatnonzero(frontier_accept & ~occupied)
        needed = target_count - int(occupied.sum())
        if frontier.size == 0:
            raise RuntimeError("QSGS growth frontier became empty before reaching target fraction")
        if frontier.size > needed:
            frontier = rng.choice(frontier, size=needed, replace=False)
        occupied.flat[frontier] = True
        generations += 1

    return occupied, generations


def component_metrics(mask: np.ndarray) -> dict[str, object]:
    structure = ndimage.generate_binary_structure(3, 1)
    labels, count = ndimage.label(mask, structure=structure)
    if count == 0:
        return {"component_count": 0, "largest_component_fraction": 0.0, "spans_xyz": [False] * 3}

    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    largest = int(sizes.max())
    spans: list[bool] = []
    for axis in range(3):
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis] = 0
        high[axis] = -1
        low_ids = np.unique(labels[tuple(low)])
        high_ids = np.unique(labels[tuple(high)])
        spans.append(bool(np.intersect1d(low_ids[low_ids != 0], high_ids[high_ids != 0]).size))
    return {
        "component_count": int(count),
        "largest_component_fraction": float(largest / int(mask.sum())),
        "spans_xyz": spans,
    }


def render_slices(skeleton: np.ndarray, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    center = tuple(n // 2 for n in skeleton.shape)
    sections = [
        ("XY midplane", skeleton[:, :, center[2]].T),
        ("XZ midplane", skeleton[:, center[1], :].T),
        ("YZ midplane", skeleton[center[0], :, :].T),
    ]
    cmap = ListedColormap(["#c7c7c7", "#315f8f"])
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, (title, section) in zip(axes, sections):
        ax.imshow(section, origin="lower", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xlabel("voxel index")
        ax.set_ylabel("voxel index")
    fig.suptitle("150³ QSGS preview: Ag (gray), Ti₃SiC₂ skeleton (blue)")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", default="150,150,150", help="three comma-separated grid dimensions")
    parser.add_argument("--p0", type=float, default=0.012, help="random nucleus probability")
    parser.add_argument("--px", type=float, default=0.11, help="growth probability along x")
    parser.add_argument("--py", type=float, default=0.11, help="growth probability along y")
    parser.add_argument("--pz", type=float, default=0.11, help="growth probability along z")
    parser.add_argument("--skeleton-fraction", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=1106)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    shape = tuple(int(value) for value in args.shape.split(","))
    if len(shape) != 3:
        parser.error("--shape must contain exactly three dimensions")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    skeleton, generations = qsgs_four_parameter(
        shape,
        args.p0,
        (args.px, args.py, args.pz),
        args.skeleton_fraction,
        args.seed,
    )
    ag = ~skeleton
    sc_metrics = component_metrics(skeleton)
    ag_metrics = component_metrics(ag)
    thickness_field = ndimage.distance_transform_edt(skeleton)
    medial_voxels = skeleton & (
        thickness_field == ndimage.maximum_filter(thickness_field, size=3, mode="constant")
    )
    thickness_diameters = 2.0 * thickness_field[medial_voxels]
    npz_path = args.output_dir / "qsgs_150x150x150.npz"
    np.savez_compressed(
        npz_path,
        skeleton=skeleton,
        ag=ag,
        shape_xyz=np.asarray(shape, dtype=np.int32),
        qsgs_parameters=np.asarray([args.p0, args.px, args.py, args.pz], dtype=np.float64),
        random_seed=np.asarray(args.seed, dtype=np.int64),
    )
    summary = {
        "status": "GEOMETRY_PREVIEW_GENERATED",
        "shape_xyz": list(shape),
        "qsgs_method": "random four-parameter face-growth",
        "parameters": {"P0": args.p0, "Px": args.px, "Py": args.py, "Pz": args.pz},
        "target_skeleton_fraction": args.skeleton_fraction,
        "actual_skeleton_fraction": float(skeleton.mean()),
        "actual_ag_fraction": float(ag.mean()),
        "random_seed": args.seed,
        "growth_generations": generations,
        "skeleton_connectivity_6face": sc_metrics,
        "ag_connectivity_6face": ag_metrics,
        "skeleton_local_diameter_voxels_p50_p90_p99": np.quantile(
            thickness_diameters, [0.50, 0.90, 0.99]
        ).round(3).tolist(),
        "both_phases_span_all_three_axes": all(sc_metrics["spans_xyz"]) and all(ag_metrics["spans_xyz"]),
        "thermal_or_flow_simulation_run": False,
    }
    summary_path = args.output_dir / "geometry_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    render_slices(skeleton, args.output_dir / "geometry_slices.png")
    print(json.dumps({"geometry": str(npz_path), "summary": str(summary_path), **summary}, indent=2))


if __name__ == "__main__":
    main()
