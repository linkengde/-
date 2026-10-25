"""Exact source-geometry Ag contact fractions on projected finite-volume faces."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent


def _breaks(lo: float, hi: float, pitch: float) -> np.ndarray:
    first = int(np.floor(lo / pitch)) + 1
    last = int(np.ceil(hi / pitch))
    internal = np.arange(first, last, dtype=float) * pitch
    internal = internal[(internal > lo) & (internal < hi)]
    return np.concatenate(([lo], internal, [hi]))


def _source_index(position: float, pitch: float, count: int) -> int:
    return min(max(int(np.floor(position / pitch)), 0), count - 1)


def exact_face_contact_fractions(
    source_ag: np.ndarray,
    coarse_ag_fraction: np.ndarray,
    size_A: np.ndarray,
    spacing_A: np.ndarray,
    positive_tolerance: float = 1.0e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return area fractions of true Ag contact on each positive-direction face.

    For each projected face, tangential subareas are split at source voxel
    boundaries. The returned face field is therefore exact for the saved
    piecewise-constant binary QSGS mask. x/y wrap periodically; z is open.
    Array entry ``axis[i,j,k]`` describes the face from that cell toward its
    positive neighbor (including the x/y periodic seam).
    """
    source_ag = np.asarray(source_ag, dtype=bool)
    coarse_ag = np.asarray(coarse_ag_fraction, dtype=float)
    size = np.asarray(size_A, dtype=float)
    spacing = np.asarray(spacing_A, dtype=float)
    if source_ag.ndim != 3 or coarse_ag.ndim != 3:
        raise ValueError("Source and coarse Ag fields must be three-dimensional")
    if size.shape != (3,) or spacing.shape != (3,) or np.any(size <= 0.0) or np.any(spacing <= 0.0):
        raise ValueError("size_A and spacing_A must contain three positive values")
    if np.any(~np.isfinite(coarse_ag)) or np.any(coarse_ag < 0.0) or np.any(coarse_ag > 1.0):
        raise ValueError("Coarse Ag fractions must be finite and lie in [0,1]")
    if not np.allclose(np.asarray(coarse_ag.shape) * spacing, size, rtol=0.0, atol=1.0e-8):
        raise ValueError("Coarse grid dimensions do not match size_A")

    source_pitch = size / np.asarray(source_ag.shape, dtype=float)
    faces = [np.zeros_like(coarse_ag, dtype=float) for _ in range(3)]
    for axis in range(3):
        tangential = [d for d in range(3) if d != axis]
        cell_count = coarse_ag.shape[axis]
        face_ids = list(range(1, cell_count))
        if axis < 2:
            face_ids.append(0)
        for face_id in face_ids:
            left_axis = (face_id - 1) % cell_count
            right_axis = face_id % cell_count
            if face_id == 0:
                left_position = size[axis] - source_pitch[axis] * 1.0e-9
                right_position = source_pitch[axis] * 1.0e-9
            else:
                face_position = face_id * spacing[axis]
                left_position = face_position - source_pitch[axis] * 1.0e-9
                right_position = face_position + source_pitch[axis] * 1.0e-9
            source_left_axis = _source_index(left_position, source_pitch[axis], source_ag.shape[axis])
            source_right_axis = _source_index(right_position, source_pitch[axis], source_ag.shape[axis])

            for first_tangent in range(coarse_ag.shape[tangential[0]]):
                lo0 = first_tangent * spacing[tangential[0]]
                hi0 = (first_tangent + 1) * spacing[tangential[0]]
                breaks0 = _breaks(lo0, hi0, source_pitch[tangential[0]])
                mid0 = 0.5 * (breaks0[:-1] + breaks0[1:])
                source_index0 = np.minimum(
                    (mid0 / source_pitch[tangential[0]]).astype(int),
                    source_ag.shape[tangential[0]] - 1,
                )
                for second_tangent in range(coarse_ag.shape[tangential[1]]):
                    lo1 = second_tangent * spacing[tangential[1]]
                    hi1 = (second_tangent + 1) * spacing[tangential[1]]
                    breaks1 = _breaks(lo1, hi1, source_pitch[tangential[1]])
                    mid1 = 0.5 * (breaks1[:-1] + breaks1[1:])
                    source_index1 = np.minimum(
                        (mid1 / source_pitch[tangential[1]]).astype(int),
                        source_ag.shape[tangential[1]] - 1,
                    )

                    coarse_left = [0, 0, 0]
                    coarse_right = [0, 0, 0]
                    coarse_left[axis] = left_axis
                    coarse_right[axis] = right_axis
                    coarse_left[tangential[0]] = coarse_right[tangential[0]] = first_tangent
                    coarse_left[tangential[1]] = coarse_right[tangential[1]] = second_tangent
                    if coarse_ag[tuple(coarse_left)] <= positive_tolerance or coarse_ag[tuple(coarse_right)] <= positive_tolerance:
                        continue

                    source_left = [0, 0, 0]
                    source_right = [0, 0, 0]
                    source_left[axis] = source_left_axis
                    source_right[axis] = source_right_axis
                    contact_area = 0.0
                    for i0, index0 in enumerate(source_index0):
                        source_left[tangential[0]] = source_right[tangential[0]] = int(index0)
                        width0 = breaks0[i0 + 1] - breaks0[i0]
                        for i1, index1 in enumerate(source_index1):
                            source_left[tangential[1]] = source_right[tangential[1]] = int(index1)
                            if source_ag[tuple(source_left)] and source_ag[tuple(source_right)]:
                                contact_area += width0 * (breaks1[i1 + 1] - breaks1[i1])
                    face_area = (hi0 - lo0) * (hi1 - lo1)
                    faces[axis][tuple(coarse_left)] = contact_area / face_area
    return tuple(faces)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    with np.load(MODEL_DIR / "qsgs_morphology.npz", allow_pickle=False) as data:
        source_ag = ~data["skeleton"].astype(bool)
        size_A = data["size_A"].astype(float)
    with np.load(HERE / f"qsgs_phase_fractions_{args.spacing}.npz", allow_pickle=False) as data:
        ag_fraction = data["ag_fraction"].astype(float)
        spacing_A = data["spacing_A"].astype(float)
    contact_x, contact_y, contact_z = exact_face_contact_fractions(
        source_ag, ag_fraction, size_A, spacing_A
    )
    output = args.output or HERE / f"qsgs_ag_face_contacts_{args.spacing}.npz"
    np.savez_compressed(
        output,
        contact_x=contact_x,
        contact_y=contact_y,
        contact_z=contact_z,
        ag_fraction=ag_fraction,
        spacing_A=spacing_A,
        size_A=size_A,
    )
    report = {
        "status": "EXACT_QSGS_AG_FACE_CONTACT_FRACTIONS",
        "spacing_label": args.spacing,
        "grid_shape_xyz": list(ag_fraction.shape),
        "source_ag_fraction": float(source_ag.mean()),
        "coarse_ag_fraction": float(ag_fraction.mean()),
        "contact_face_fraction_mean": {
            axis: float(np.mean(values[values > 0.0])) if np.any(values > 0.0) else 0.0
            for axis, values in zip("xyz", (contact_x, contact_y, contact_z))
        },
        "nonzero_face_count": {
            axis: int(np.count_nonzero(values))
            for axis, values in zip("xyz", (contact_x, contact_y, contact_z))
        },
        "output_npz": output.name,
        "method": "exact area integration of piecewise-constant source Ag voxels over each coarse face",
        "limitations": [
            "This corrects local electron conduction links only; it does not infer interface resistance.",
            "The correction is a finite-volume geometry approximation, not atomistic transport.",
            "No z-boundary contact factor is needed for the no-sink adiabatic case used in the pilot.",
        ],
    }
    report_path = output.with_suffix(".json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
