#!/usr/bin/env python3
"""Generate and audit exactly nested 10 Å and 5 Å QSGS thermal grids."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from qsgs_face_contacts import exact_face_contact_fractions
from qsgs_projection import project_phase_mask
from qsgs_source_integration import integrate_qsgs_source_power_W


HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent
SOURCE = MODEL_DIR / "qsgs_morphology.npz"
COARSE_SHAPE = (30, 31, 16)
FINE_SHAPE = (60, 62, 32)
PEAK_FLUX_W_M2 = 5.0 * 1.602176634e-19 / 1.0e-12 / 1.0e-20
GAUSSIAN_RADIUS_A = 30.0
ABSORPTION_DEPTH_A = 20.0


def block_average_2x(field: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=float)
    if any(length % 2 for length in field.shape):
        raise ValueError("Fine-grid dimensions must be even for 2× coarsening")
    nx, ny, nz = (length // 2 for length in field.shape)
    return field.reshape(nx, 2, ny, 2, nz, 2).mean(axis=(1, 3, 5))


def block_sum_2x(field: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=float)
    if any(length % 2 for length in field.shape):
        raise ValueError("Fine-grid dimensions must be even for 2× coarsening")
    nx, ny, nz = (length // 2 for length in field.shape)
    return field.reshape(nx, 2, ny, 2, nz, 2).sum(axis=(1, 3, 5))


def coarsen_face_contacts_2x(
    fine_contacts: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Area-average fine positive-face contacts onto aligned coarse faces."""
    cx_f, cy_f, cz_f = (np.asarray(value, dtype=float) for value in fine_contacts)
    nx, ny, nz = (length // 2 for length in cx_f.shape)
    cx = cx_f[1::2].reshape(nx, ny, 2, nz, 2).mean(axis=(2, 4))
    cy = cy_f[:, 1::2, :].reshape(nx, 2, ny, nz, 2).mean(axis=(1, 4))
    cz = cz_f[:, :, 1::2].reshape(nx, 2, ny, 2, nz).mean(axis=(1, 3))
    return cx, cy, cz


def relative_l1(actual: np.ndarray, expected: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    denominator = float(np.sum(np.abs(expected)))
    return float(np.sum(np.abs(actual - expected)) / denominator) if denominator else 0.0


def save_projection(path: Path, result, requested_spacing_A: float) -> None:
    np.savez_compressed(
        path,
        ag_fraction=result.ag_fraction,
        skeleton_fraction=result.skeleton_fraction,
        size_A=result.size_A,
        spacing_A=result.spacing_A,
        requested_spacing_A=requested_spacing_A,
    )


def save_contacts(path: Path, fractions: np.ndarray, size_A: np.ndarray, spacing_A: np.ndarray,
                  contacts: tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    np.savez_compressed(
        path,
        contact_x=contacts[0],
        contact_y=contacts[1],
        contact_z=contacts[2],
        ag_fraction=fractions,
        spacing_A=spacing_A,
        size_A=size_A,
    )


def main() -> None:
    with np.load(SOURCE, allow_pickle=False) as source:
        skeleton_mask = np.asarray(source["skeleton"], dtype=bool)
        size_A = np.asarray(source["size_A"], dtype=float)
    results = {
        "10A_nested": project_phase_mask(skeleton_mask, size_A, 10.0, COARSE_SHAPE),
        "5A_nested": project_phase_mask(skeleton_mask, size_A, 5.0, FINE_SHAPE),
    }
    contacts: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for label, result in results.items():
        phase_path = HERE / f"qsgs_phase_fractions_{label}.npz"
        face_path = HERE / f"qsgs_ag_face_contacts_{label}.npz"
        save_projection(phase_path, result, 10.0 if label.startswith("10A") else 5.0)
        contacts[label] = exact_face_contact_fractions(
            ~skeleton_mask,
            result.ag_fraction,
            size_A,
            result.spacing_A,
        )
        save_contacts(face_path, result.ag_fraction, size_A, result.spacing_A, contacts[label])

    coarse, fine = results["10A_nested"], results["5A_nested"]
    coarse_contacts_from_fine = coarsen_face_contacts_2x(contacts["5A_nested"])
    source_power_coarse = integrate_qsgs_source_power_W(
        ~skeleton_mask, COARSE_SHAPE, size_A, PEAK_FLUX_W_M2,
        GAUSSIAN_RADIUS_A, ABSORPTION_DEPTH_A,
    )
    source_power_fine = integrate_qsgs_source_power_W(
        ~skeleton_mask, FINE_SHAPE, size_A, PEAK_FLUX_W_M2,
        GAUSSIAN_RADIUS_A, ABSORPTION_DEPTH_A,
    )
    source_power_fine_coarsened = block_sum_2x(source_power_fine)
    ag_fine_coarsened = block_average_2x(fine.ag_fraction)
    face_metrics = {}
    for axis, coarse_face, fine_face in zip(
        "xyz", contacts["10A_nested"], coarse_contacts_from_fine
    ):
        face_metrics[axis] = {
            "max_abs_area_fraction_difference": float(np.max(np.abs(coarse_face - fine_face))),
            "relative_L1_area_weighted_difference": relative_l1(fine_face, coarse_face),
            "coarse_nonzero_face_count": int(np.count_nonzero(coarse_face > 1.0e-12)),
            "nested_fine_aggregated_nonzero_face_count": int(np.count_nonzero(fine_face > 1.0e-12)),
        }
    coarse_volume = float(np.prod(coarse.spacing_A))
    report = {
        "status": "STRICTLY_NESTED_QSGS_GRID_GEOMETRY_AUDIT",
        "source_qsgs_shape_xyz": list(skeleton_mask.shape),
        "physical_size_A": size_A.tolist(),
        "coarse_shape_xyz": list(COARSE_SHAPE),
        "fine_shape_xyz": list(FINE_SHAPE),
        "strict_2x_shape_relation": [FINE_SHAPE[i] == 2 * COARSE_SHAPE[i] for i in range(3)],
        "coarse_spacing_A_xyz": coarse.spacing_A.tolist(),
        "fine_spacing_A_xyz": fine.spacing_A.tolist(),
        "ag_volume_fraction_source": float(np.mean(~skeleton_mask)),
        "ag_volume_fraction_coarse": float(coarse.ag_fraction.mean()),
        "ag_volume_fraction_fine": float(fine.ag_fraction.mean()),
        "fine_ag_fraction_block_average_relative_L1_vs_coarse": relative_l1(
            ag_fine_coarsened, coarse.ag_fraction
        ),
        "fine_ag_fraction_block_average_max_abs_difference": float(
            np.max(np.abs(ag_fine_coarsened - coarse.ag_fraction))
        ),
        "coarse_source_power_W": float(source_power_coarse.sum()),
        "fine_source_power_W": float(source_power_fine.sum()),
        "fine_source_power_block_sum_relative_L1_vs_coarse": relative_l1(
            source_power_fine_coarsened, source_power_coarse
        ),
        "face_contact_area_audit": face_metrics,
        "phase_files": {
            "coarse": "qsgs_phase_fractions_10A_nested.npz",
            "fine_nested": "qsgs_phase_fractions_5A_nested.npz",
        },
        "face_contact_files": {
            "coarse": "qsgs_ag_face_contacts_10A_nested.npz",
            "fine_nested": "qsgs_ag_face_contacts_5A_nested.npz",
        },
        "limitations": [
            "This is a geometry and source-integral audit only; it does not establish transport convergence.",
            "The nested fine grid has slightly anisotropic cell widths because the physical box lengths are not integer multiples of 5 Å.",
        ],
    }
    path = HERE / "qsgs_nested_grid_geometry_audit.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
