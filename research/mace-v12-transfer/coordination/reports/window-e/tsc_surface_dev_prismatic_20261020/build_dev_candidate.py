#!/workspace/.venvs/gpaw-mpi/bin/python
"""Construct a frozen prismatic Ti3SiC2 DEV-candidate parent and response pair."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import ase
from ase.build import surface
from ase.io import read, write


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
CIF_REL = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
)
EXPECTED_CIF_SHA = "6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3"
VACUUM_A = 18.0
LAYERS = 8
DISPLACEMENT_A = 0.02


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def inversion_map_error(atoms) -> tuple[float, list[tuple[int, int]]]:
    """Match each atom to its same-species inversion image about slab center."""
    symbols = atoms.get_chemical_symbols()
    positions = np.asarray(atoms.positions)
    frac = atoms.get_scaled_positions(wrap=True)
    zlo, zhi = float(positions[:, 2].min()), float(positions[:, 2].max())
    errors, mapping = [], []
    for i, sym in enumerate(symbols):
        target_xy = (-frac[i, :2]) % 1.0
        target_z = zlo + zhi - positions[i, 2]
        candidates = []
        for j, other in enumerate(symbols):
            if other != sym:
                continue
            dxy = (frac[j, :2] - target_xy + 0.5) % 1.0 - 0.5
            drxy = dxy @ atoms.cell[:2]
            dz = positions[j, 2] - target_z
            candidates.append((float(np.sqrt(np.dot(drxy, drxy) + dz * dz)), j))
        error, j = min(candidates)
        errors.append(error)
        if i < j:
            mapping.append((i + 1, j + 1))
    return max(errors), mapping


def min_distance(atoms) -> tuple[float, int, int]:
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return float(d[i, j]), int(i + 1), int(j + 1)


def terminal_plane(atoms, which: str) -> dict:
    z = atoms.positions[:, 2]
    edge = float(z.min() if which == "bottom" else z.max())
    ids = np.flatnonzero(np.abs(z - edge) < 1e-7)
    symbols = atoms.get_chemical_symbols()
    return {
        "side": which,
        "z_angstrom": edge,
        "atom_ids_1based": [int(i + 1) for i in ids],
        "composition": dict(sorted(Counter(symbols[i] for i in ids).items())),
    }


def main() -> None:
    if ase.__version__ != "3.29.0":
        raise RuntimeError(f"Expected ASE 3.29.0, got {ase.__version__}")
    cif = REPO / CIF_REL
    cif_sha = sha256(cif)
    if cif_sha != EXPECTED_CIF_SHA:
        raise RuntimeError(f"COD CIF SHA mismatch: {cif_sha}")
    bulk = read(cif, format="cif")
    if len(bulk) != 12 or bulk.get_chemical_formula() != "C4Si2Ti6":
        raise RuntimeError("Unexpected COD unit-cell composition")

    # Translating by half a bulk a-vector selects the symmetric mixed Ti/C
    # termination for the (100) cut; it does not strain or change bulk geometry.
    shifted = bulk.copy()
    shifted.positions += 0.5 * shifted.cell[0]
    slab = surface(shifted, (1, 0, 0), LAYERS, vacuum=None, periodic=False)
    slab.center(vacuum=VACUUM_A / 2, axis=2)
    slab.set_pbc((True, True, True))
    slab.new_array("source_atom_id", np.arange(1, len(slab) + 1, dtype=np.int64))

    counts = dict(sorted(Counter(slab.get_chemical_symbols()).items()))
    expected_counts = {"C": 32, "Si": 16, "Ti": 48}
    if len(slab) != 96 or counts != expected_counts:
        raise RuntimeError(f"Unexpected slab size/composition: N={len(slab)}, {counts}")
    base_symmetry_error, inversion_pairs = inversion_map_error(slab)
    if base_symmetry_error > 1e-5:
        raise RuntimeError(f"Base slab inversion symmetry failed: {base_symmetry_error} A")
    min_d, min_i, min_j = min_distance(slab)
    if min_d < 1.5:
        raise RuntimeError(f"Unphysical short contact {min_d} A at IDs {min_i},{min_j}")

    z = slab.positions[:, 2]
    lo = np.flatnonzero(np.abs(z - z.min()) < 1e-7)
    hi = np.flatnonzero(np.abs(z - z.max()) < 1e-7)
    symbols = slab.get_chemical_symbols()
    ti_pairs = []
    for i in lo:
        if symbols[i] != "Ti":
            continue
        target_xy = (-slab.get_scaled_positions(wrap=True)[i, :2]) % 1.0
        target_z = z.min() + z.max() - z[i]
        for j in hi:
            if symbols[j] != "Ti":
                continue
            frac = slab.get_scaled_positions(wrap=True)
            dxy = (frac[j, :2] - target_xy + 0.5) % 1.0 - 0.5
            drxy = dxy @ slab.cell[:2]
            dz = z[j] - target_z
            err = float(np.sqrt(np.dot(drxy, drxy) + dz * dz))
            if err < 1e-5:
                ti_pairs.append((int(i + 1), int(j + 1)))
    if not ti_pairs:
        raise RuntimeError("No inversion-related Ti terminal pair for response point")
    bottom_id, top_id = min(ti_pairs)

    family = "surface_prismatic_10m10_cod9009647_8layer"
    root_id = "TSC96_DEV_PRISM_10m10_COD9009647_8L"
    for atoms, structure_id, role in [
        (slab, f"{root_id}_REF", "DEV_GEOMETRY_CANDIDATE"),
    ]:
        atoms.info.update({
            "structure_id": structure_id,
            "source_family": family,
            "candidate_role": role,
            "material": "Ti3SiC2",
            "surface_orientation_ASE_Miller": [1, 0, 0],
            "surface_orientation_hexagonal_4index_equivalent": "(10-10)",
            "surface_termination": "symmetric mixed Ti-C / Ti-C prismatic faces",
            "source_cif": CIF_REL.as_posix(),
            "source_cif_sha256": cif_sha,
            "construction": "ASE surface; conventional bulk cell translated by 0.5*a before (100) cleavage; 8 layers; no relaxation",
            "shared_crystallographic_root_with_TRAIN": True,
            "independence_scope": "new surface-orientation geometry family; not an independent crystallographic-source blind test",
            "vacuum_gap_angstrom": VACUUM_A,
            "fixed_geometry": True,
        })
    ref = slab.copy()
    response = slab.copy()
    response.info = dict(slab.info)
    response.info.update({
        "structure_id": f"{root_id}_RESP_TIZ020",
        "candidate_role": "DEV_GEOMETRY_RESPONSE_CANDIDATE",
        "displacement_definition": "inversion-related terminal Ti pair; bottom -surface-normal and top +surface-normal",
        "displacement_magnitude_angstrom": DISPLACEMENT_A,
        "displaced_source_atom_ids": [bottom_id, top_id],
    })
    response.positions[bottom_id - 1, 2] -= DISPLACEMENT_A
    response.positions[top_id - 1, 2] += DISPLACEMENT_A
    response_symmetry_error, _ = inversion_map_error(response)
    if response_symmetry_error > 1e-5:
        raise RuntimeError(f"Response inversion symmetry failed: {response_symmetry_error} A")
    pert_min, pert_i, pert_j = min_distance(response)
    if pert_min < 1.5:
        raise RuntimeError(f"Response creates short contact {pert_min} A at IDs {pert_i},{pert_j}")

    out = ROOT / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    ref_path = out / f"{root_id}_REF.extxyz"
    response_path = out / f"{root_id}_RESP_TIZ020.extxyz"
    write(ref_path, ref, format="extxyz")
    write(response_path, response, format="extxyz")

    entries = []
    for atoms, path, role, md, mi, mj, symerr in [
        (ref, ref_path, "reference", min_d, min_i, min_j, base_symmetry_error),
        (response, response_path, "displacement_response", pert_min, pert_i, pert_j, response_symmetry_error),
    ]:
        entries.append({
            "structure_id": atoms.info["structure_id"],
            "role_candidate": role,
            "source_family": family,
            "path": path.relative_to(ROOT).as_posix(),
            "file_sha256": sha256(path),
            "geometry_sha256": geometry_sha(atoms),
            "atom_count": len(atoms),
            "composition": counts,
            "cell_angstrom": np.asarray(atoms.cell).tolist(),
            "pbc": atoms.pbc.tolist(),
            "surface_span_angstrom": float(z.max() - z.min()),
            "vacuum_gap_angstrom": float(atoms.cell[2, 2] - (z.max() - z.min())),
            "surface_termination_planes": [terminal_plane(atoms, "bottom"), terminal_plane(atoms, "top")],
            "minimum_distance_angstrom": md,
            "minimum_distance_atom_ids": [mi, mj],
            "inversion_symmetry_max_error_angstrom": symerr,
            "displaced_source_atom_ids": atoms.info.get("displaced_source_atom_ids", []),
            "DFT_label_present": False,
        })

    manifest = {
        "status": "DEV_CANDIDATE_GEOMETRY_PAIR_GENERATED_NO_ROLE_CHANGE",
        "role_assignment": "candidate only; TRAIN/DEV/TEST manifest unchanged",
        "candidate_family": family,
        "source_cif": CIF_REL.as_posix(),
        "source_cif_sha256": cif_sha,
        "ASE_version": ase.__version__,
        "builder_sha256": sha256(Path(__file__).resolve()),
        "construction": {
            "Miller_indices_conventional_hex_cell": [1, 0, 0],
            "four_index_equivalent": "(10-10)",
            "layers": LAYERS,
            "bulk_translation_before_cleavage": "0.5*a; selects equivalent Ti-C mixed prismatic terminations and inversion symmetry",
            "relaxation": False,
            "cell_strain": False,
            "vacuum_angstrom": VACUUM_A,
            "pbc": [True, True, True],
            "displacement_angstrom": DISPLACEMENT_A,
            "displaced_pair_source_atom_ids": [bottom_id, top_id],
            "displacement_direction": "bottom -z, top +z; z is surface normal",
        },
        "family_independence_note": "Distinct surface orientation and geometry-generation family from TRAIN basal (0001) slabs, but shares their COD 9009647 crystallographic root; not a source-independent blind test.",
        "entries": entries,
    }
    (ROOT / "STRUCTURE_AUDIT.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"], "entries": entries}, indent=2))


if __name__ == "__main__":
    main()
