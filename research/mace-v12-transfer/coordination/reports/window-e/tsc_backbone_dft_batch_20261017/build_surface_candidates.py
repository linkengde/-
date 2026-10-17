#!/usr/bin/env python3
"""Build four frozen Ti3SiC2 basal-surface DFT candidates from the archived COD CIF.

This script only creates geometries. It does not construct a calculator or run DFT.
Run from the repository root with the already-installed project ASE environment.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write


EXPECTED_ASE = "3.29.0"
EXPECTED_CIF_SHA256 = "6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3"
EXPECTED_PARENT_INPUT_SHA256 = "157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4"
VACUUM_GAP_A = 18.0
DISPLACEMENT_A = 0.02


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def coordinate_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    # ASE extxyz writes Cartesian positions to 8 decimal places by default.
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def minimum_distance(atoms) -> tuple[float, int, int]:
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return float(d[i, j]), int(i + 1), int(j + 1)


def symmetry_error(atoms, operation: str) -> float:
    symbols = atoms.get_chemical_symbols()
    frac = atoms.get_scaled_positions(wrap=True)
    zlo = float(atoms.positions[:, 2].min())
    zhi = float(atoms.positions[:, 2].max())
    errs = []
    for i, symbol in enumerate(symbols):
        if operation == "inversion_about_slab_center":
            xy_target = (-frac[i, :2]) % 1.0
        elif operation == "horizontal_mirror_about_slab_center":
            xy_target = frac[i, :2]
        else:
            raise ValueError(operation)
        z_target = zlo + zhi - float(atoms.positions[i, 2])
        candidates = []
        for j, other in enumerate(symbols):
            if other != symbol:
                continue
            dxy = (frac[j, :2] - xy_target + 0.5) % 1.0 - 0.5
            dz = float(atoms.positions[j, 2]) - z_target
            candidates.append(float(np.sqrt(np.dot(dxy @ atoms.cell[:2], dxy @ atoms.cell[:2]) + dz * dz)))
        errs.append(min(candidates))
    return max(errs)


def terminal_pair(atoms, element: str) -> tuple[int, int]:
    z = atoms.positions[:, 2]
    low = np.flatnonzero(np.abs(z - z.min()) < 1e-7)
    high = np.flatnonzero(np.abs(z - z.max()) < 1e-7)
    symbols = atoms.get_chemical_symbols()
    frac = atoms.get_scaled_positions(wrap=True)
    pairs = []
    for i in low:
        if symbols[i] != element:
            continue
        for j in high:
            if symbols[j] != element:
                continue
            delta = (frac[j, :2] - frac[i, :2] + 0.5) % 1.0 - 0.5
            if np.linalg.norm(delta) < 1e-8:
                pairs.append((int(i + 1), int(j + 1)))
    if not pairs:
        raise RuntimeError(f"No matching top/bottom {element} terminal site pair")
    return min(pairs)


def build_slab(bulk, lo: float, hi: float, termination: str, output_dir: Path):
    c0 = float(bulk.cell[2, 2])
    large = bulk.repeat((2, 2, 3))
    z_in_bulk_cells = large.positions[:, 2] / c0
    mask = (z_in_bulk_cells >= lo - 1e-8) & (z_in_bulk_cells <= hi + 1e-8)
    slab = large[mask]
    zlo = float(slab.positions[:, 2].min())
    zhi = float(slab.positions[:, 2].max())
    span = zhi - zlo
    if not np.isclose(span, 2 * c0, atol=1e-6):
        raise RuntimeError(f"Unexpected slab thickness {span}; expected 2c={2*c0}")

    cell = np.asarray(slab.cell).copy()
    cell[2] = [0.0, 0.0, span + VACUUM_GAP_A]
    slab.set_cell(cell, scale_atoms=False)
    slab.positions[:, 2] += VACUUM_GAP_A / 2 - zlo
    slab.set_pbc((True, True, True))
    slab.new_array("source_atom_id", np.arange(1, len(slab) + 1, dtype=np.int64))

    expected_counts = (
        {"Ti": 52, "Si": 16, "C": 32}
        if termination == "Ti"
        else {"Ti": 48, "Si": 20, "C": 32}
    )
    counts = {s: slab.get_chemical_symbols().count(s) for s in ("Ti", "Si", "C")}
    if counts != expected_counts or len(slab) != 100:
        raise RuntimeError(f"Unexpected {termination}-slab composition/count: {counts}, N={len(slab)}")
    ends = []
    for z in (float(slab.positions[:, 2].min()), float(slab.positions[:, 2].max())):
        ids = np.flatnonzero(np.abs(slab.positions[:, 2] - z) < 1e-7)
        ends.append(sorted(set(slab.get_chemical_symbols()[i] for i in ids)))
    if ends != [[termination], [termination]]:
        raise RuntimeError(f"Requested {termination} terminations, got {ends}")

    operation = "inversion_about_slab_center" if termination == "Ti" else "horizontal_mirror_about_slab_center"
    base_symmetry_error = symmetry_error(slab, operation)
    if base_symmetry_error > 1e-4:
        raise RuntimeError(f"Slab does not pass symmetry check: max error {base_symmetry_error} A")
    min_d, ai, aj = minimum_distance(slab)
    if min_d < 1.5:
        raise RuntimeError(f"Unphysical short contact {min_d} A between IDs {ai},{aj}")

    root_id = f"TSC100_{termination.upper()}0001_COD9009647_2x2x2"
    slab.info.update({
        "structure_id": f"{root_id}_REF",
        "source_family": f"surface_{termination.lower()}0001_cod9009647",
        "candidate_role": "TRAIN_CANDIDATE_PENDING_C_REVIEW",
        "source_cif": "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif",
        "source_cif_sha256": EXPECTED_CIF_SHA256,
        "geometry_status": "as-cleaved experimental-CIF geometry; fixed-cell; not relaxed",
        "surface_orientation": "(0001)",
        "surface_termination": termination,
        "vacuum_gap_angstrom": VACUUM_GAP_A,
        "symmetry_check": operation,
        "symmetry_max_error_angstrom": base_symmetry_error,
    })
    ref_path = output_dir / f"{root_id}_REF.extxyz"
    write(ref_path, slab, format="extxyz")

    bottom_id, top_id = terminal_pair(slab, termination)
    pert = slab.copy()
    pert.info = dict(slab.info)
    pert.info.update({
        "structure_id": f"{root_id}_RESP_Z020",
        "candidate_role": "TRAIN_CANDIDATE_PENDING_C_REVIEW",
        "displacement_mode": "one inversion/mirror-related terminal-site pair; bottom -z and top +z",
        "displacement_magnitude_angstrom": DISPLACEMENT_A,
        "displaced_source_atom_ids": [bottom_id, top_id],
    })
    pert.positions[bottom_id - 1, 2] -= DISPLACEMENT_A
    pert.positions[top_id - 1, 2] += DISPLACEMENT_A
    pert_symmetry_error = symmetry_error(pert, operation)
    if pert_symmetry_error > 1e-4:
        raise RuntimeError(f"Perturbed slab does not preserve symmetry: {pert_symmetry_error} A")
    pmin, pai, paj = minimum_distance(pert)
    if pmin < 1.5:
        raise RuntimeError(f"Perturbation creates short contact {pmin} A between IDs {pai},{paj}")
    pert_path = output_dir / f"{root_id}_RESP_Z020.extxyz"
    write(pert_path, pert, format="extxyz")

    rows = []
    for role, atoms, path, md, i, j, symerr in [
        ("reference", slab, ref_path, min_d, ai, aj, base_symmetry_error),
        ("response", pert, pert_path, pmin, pai, paj, pert_symmetry_error),
    ]:
        rows.append({
            "structure_id": atoms.info["structure_id"],
            "role_candidate": role,
            "source_family": atoms.info["source_family"],
            "path": path.relative_to(output_dir.parent).as_posix(),
            "file_sha256": sha256(path),
            "coordinate_sha256": coordinate_sha(atoms),
            "atom_count": len(atoms),
            "composition": counts,
            "cell_angstrom": np.asarray(atoms.cell).tolist(),
            "pbc": atoms.pbc.tolist(),
            "surface_span_angstrom": span,
            "vacuum_gap_angstrom": float(atoms.cell[2, 2] - span),
            "top_bottom_termination": ends,
            "symmetry_operation": operation,
            "symmetry_max_error_angstrom": symerr,
            "minimum_distance_angstrom": md,
            "minimum_distance_atom_ids": [i, j],
            "displaced_source_atom_ids": atoms.info.get("displaced_source_atom_ids", []),
            "dft_label_present": False,
        })
    return rows


def main() -> None:
    import ase

    if ase.__version__ != EXPECTED_ASE:
        raise RuntimeError(f"Expected ASE {EXPECTED_ASE}, got {ase.__version__}")
    script_dir = Path(__file__).resolve().parent
    root = next((p for p in script_dir.parents if (p / ".git").exists()), None)
    if root is None:
        raise RuntimeError("Could not locate repository root")
    cif_rel = Path("research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif")
    parent_rel = Path("research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz")
    cif_path, parent_path = root / cif_rel, root / parent_rel
    if sha256(cif_path) != EXPECTED_CIF_SHA256:
        raise RuntimeError("Source CIF SHA mismatch")
    if sha256(parent_path) != EXPECTED_PARENT_INPUT_SHA256:
        raise RuntimeError("Parent proposal input SHA mismatch")
    bulk = read(cif_path)
    if len(bulk) != 12 or bulk.get_chemical_formula() != "C4Si2Ti6":
        raise RuntimeError("Unexpected COD unit-cell chemistry")
    if not np.allclose(np.asarray(bulk.cell)[:2], np.asarray(read(parent_path).cell)[:2] / 2, atol=1e-8):
        raise RuntimeError("COD parent in-plane cell does not reproduce the audited 48-atom parent")

    out = script_dir / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    entries = []
    entries += build_slab(bulk, 0.0, 2.0, "Ti", out)
    entries += build_slab(bulk, 0.25, 2.25, "Si", out)
    manifest = {
        "status": "FOUR_GEOMETRY_CANDIDATES_GENERATED; NO_DFT_LABELS; C_REVIEW_REQUIRED",
        "generator": Path(__file__).name,
        "generator_sha256": sha256(Path(__file__).resolve()),
        "ASE_version": ase.__version__,
        "source_cif_path": cif_rel.as_posix(),
        "source_cif_sha256": EXPECTED_CIF_SHA256,
        "parent_proposal_path": parent_rel.as_posix(),
        "parent_proposal_sha256": EXPECTED_PARENT_INPUT_SHA256,
        "construction": "COD primitive cell -> 2x2 in-plane, two c repeats plus matching terminal plane; no relaxation or DFT",
        "entries": entries,
    }
    (script_dir / "SURFACE_CANDIDATE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"], "entries": len(entries), "manifest": "SURFACE_CANDIDATE_MANIFEST.json"}, indent=2))


if __name__ == "__main__":
    main()
