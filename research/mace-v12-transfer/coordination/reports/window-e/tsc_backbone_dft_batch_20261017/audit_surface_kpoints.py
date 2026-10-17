#!/usr/bin/env python3
"""Offline audit of the proposed slab Monkhorst-Pack coordinates; no GPAW calculation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase import __version__ as ase_version
from ase.dft.kpoints import monkhorst_pack


def bytes_sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(a, dtype="<f8").tobytes(order="C")).hexdigest()


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    root = next((p for p in script_dir.parents if (p / ".git").exists()), None)
    if root is None:
        raise RuntimeError("Could not locate repository root")
    parent_log = root / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ti3SiC2_baseline_k6x6x2_v19/gpaw.log"
    parent_runner = root / "research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/run_single_point_msr1.py"
    from pathlib import Path as P
    def sha(path: P) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    bulk = monkhorst_pack((6, 6, 2))
    slab = monkhorst_pack((6, 6, 1))
    bulk_xy = np.unique(np.round(bulk[:, :2], 14), axis=0)
    slab_xy = np.unique(np.round(slab[:, :2], 14), axis=0)
    if bulk_xy.shape != slab_xy.shape or not np.array_equal(bulk_xy, slab_xy):
        raise RuntimeError("The surface in-plane MP coordinates do not equal the parent in-plane grid")
    if not np.allclose(slab[:, 2], 0.0, atol=1e-14):
        raise RuntimeError("A 6x6x1 standard MP mesh should have only z=0")
    log = parent_log.read_text(errors="replace")
    if "Monkhorst-Pack size: [6, 6, 2]" not in log or "Monkhorst-Pack shift: [0, 0, 0]" not in log:
        raise RuntimeError("Parent GPAW output does not confirm the expected MP mesh/shift")
    runner_text = parent_runner.read_text()
    if "kpts=tuple(method['kpoints']['size'])" not in runner_text:
        raise RuntimeError("Could not verify the production tuple-kpts convention in the archived runner")
    payload = {
        "status": "PASS_OFFLINE_COORDINATE_AUDIT; SURFACE_MESH_STILL_REQUIRES_C_APPROVAL",
        "ASE_version": ase_version,
        "parent_geometry": {
            "source_path": "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz",
            "mesh": [6, 6, 2],
            "mesh_type": "standard Monkhorst-Pack via GPAW kpts tuple",
            "GPAW_log_path": str(parent_log.relative_to(root)),
            "GPAW_log_sha256": sha(parent_log),
            "GPAW_runner_path": str(parent_runner.relative_to(root)),
            "GPAW_runner_sha256": sha(parent_runner),
            "coordinate_count": int(len(bulk)),
            "coordinate_sha256_float64_le": bytes_sha(bulk),
            "unique_in_plane_coordinate_count": int(len(bulk_xy)),
            "z_coordinates": sorted(set(float(x) for x in bulk[:, 2])),
        },
        "surface_candidate": {
            "mesh": [6, 6, 1],
            "mesh_type": "standard Monkhorst-Pack via GPAW kpts tuple; no explicit gamma flag or extra offset",
            "coordinate_count": int(len(slab)),
            "coordinate_sha256_float64_le": bytes_sha(slab),
            "unique_in_plane_coordinate_count": int(len(slab_xy)),
            "z_coordinates": sorted(set(float(x) for x in slab[:, 2])),
            "in_plane_grid_exactly_matches_parent": True,
            "parent_full_mesh_has_two_z_planes": True,
            "surface_mesh_has_gamma_z_only": True,
            "coordinate_list_fractional": slab.tolist(),
        },
        "interpretation": "Parent 6x6x2 numerical convergence is not a k-point convergence test for the new slab geometries. The proposed 6x6x1 preserves the parent's in-plane grid because the slab in-plane cell vectors equal the parent's vectors; z=1 removes vacuum-direction sampling.",
    }
    (script_dir / "KPOINT_REVIEW.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "surface_coordinate_sha256": payload["surface_candidate"]["coordinate_sha256_float64_le"], "points": len(slab)}, indent=2))


if __name__ == "__main__":
    main()
