#!/usr/bin/env python3
"""Generate and audit one symmetric Ti3SiC2 (10-11) geometry-only candidate.

No calculator is instantiated. No force/energy labels or data roles are created.
"""
from __future__ import annotations

import hashlib
import json
from importlib.metadata import version as package_version
from collections import Counter
from pathlib import Path

import ase
import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read, write
from scipy.optimize import linear_sum_assignment

EXPECTED_CIF_SHA256 = "6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3"
MILLER_INDEX = (1, 0, 1)  # hexagonal four-index equivalent (10-11)
MIN_SLAB_THICKNESS_A = 15.0
MIN_VACUUM_A = 18.0
DISPLACEMENT_A = 0.02
EXPECTED_COMPOSITION = {"Ti": 33, "Si": 12, "C": 24}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def coordinate_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def min_distance(atoms):
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return float(d[i, j]), [int(i + 1), int(j + 1)]


def slab_projection(atoms, normal):
    return np.asarray(atoms.positions) @ np.asarray(normal)


def inversion_map(atoms, normal):
    syms = atoms.get_chemical_symbols()
    frac = atoms.get_scaled_positions(wrap=True)
    proj = slab_projection(atoms, normal)
    lo, hi = float(proj.min()), float(proj.max())
    errors = []
    pairs = []
    for i, sym in enumerate(syms):
        target_xy = (-frac[i, :2]) % 1.0
        target_n = lo + hi - float(proj[i])
        best = (float("inf"), None)
        for j, other in enumerate(syms):
            if other != sym:
                continue
            dxy = (frac[j, :2] - target_xy + 0.5) % 1.0 - 0.5
            dcart = dxy[0] * atoms.cell[0] + dxy[1] * atoms.cell[1]
            dn = float(proj[j]) - target_n
            error = float(np.sqrt(np.dot(dcart, dcart) + dn * dn))
            if error < best[0]:
                best = (error, j)
        errors.append(best[0])
        if best[1] is not None and i < best[1]:
            pairs.append((i + 1, best[1] + 1, best[0]))
    return max(errors), pairs


def terminal_plane(atoms, normal, value, tol=1e-6):
    proj = slab_projection(atoms, normal)
    ids = np.flatnonzero(np.abs(proj - value) < tol)
    counts = dict(sorted(Counter(atoms[i].symbol for i in ids).items()))
    return ids, counts


def source_parser_crosscheck(cif: Path, parent):
    ase_parent = read(cif, format="cif")
    pmg_parent = AseAtomsAdaptor.get_atoms(parent)
    if ase_parent.get_chemical_symbols() != pmg_parent.get_chemical_symbols():
        raise RuntimeError("ASE/Pymatgen CIF parse atom order or species differ")
    if not np.allclose(ase_parent.cell.array, pmg_parent.cell.array, atol=1e-10, rtol=0):
        raise RuntimeError("ASE/Pymatgen CIF parse lattice differs")
    fa = ase_parent.get_scaled_positions(wrap=True)
    fb = pmg_parent.get_scaled_positions(wrap=True)
    all_delta = []
    for symbol in sorted(set(ase_parent.get_chemical_symbols())):
        ia = np.array([i for i, s in enumerate(ase_parent.get_chemical_symbols()) if s == symbol])
        ib = np.array([i for i, s in enumerate(pmg_parent.get_chemical_symbols()) if s == symbol])
        delta = fb[ib][None, :, :] - fa[ia][:, None, :]
        delta -= np.round(delta)
        cart = delta @ ase_parent.cell.array
        costs = np.linalg.norm(cart, axis=2)
        rows, cols = linear_sum_assignment(costs)
        all_delta.extend(costs[rows, cols].tolist())
    return {
        "ASE_version": ase.__version__,
        "Pymatgen_version": package_version("pymatgen"),
        "same_cell_and_species_order": True,
        "max_matched_parent_site_difference_angstrom": float(max(all_delta)),
        "rms_matched_parent_site_difference_angstrom": float(np.sqrt(np.mean(np.square(all_delta)))),
        "note": "Pymatgen emitted a CIF normalization warning for four fractional coordinates; matched against direct ASE CIF parsing.",
    }


def main():
    root = Path(__file__).resolve().parent
    repo = next(p for p in root.parents if (p / ".git").exists())
    cif = repo / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
    source_hash = sha256(cif)
    if source_hash != EXPECTED_CIF_SHA256:
        raise RuntimeError(f"Unexpected CIF SHA-256: {source_hash}")

    bulk = Structure.from_file(cif)
    if len(bulk) != 12 or bulk.composition.reduced_formula != "Ti3SiC2":
        raise RuntimeError(f"Unexpected parent structure: {len(bulk)} sites, {bulk.composition.reduced_formula}")
    parser_crosscheck = source_parser_crosscheck(cif, bulk)
    generator = SlabGenerator(
        bulk,
        MILLER_INDEX,
        min_slab_size=MIN_SLAB_THICKNESS_A,
        min_vacuum_size=MIN_VACUUM_A,
        center_slab=True,
        primitive=True,
        reorient_lattice=True,
    )
    slabs = generator.get_slabs(symmetrize=True, tol=1e-3)
    candidates = []
    for slab in slabs:
        composition = dict(sorted(Counter(site.specie.symbol for site in slab).items()))
        if composition != EXPECTED_COMPOSITION:
            continue
        atoms = AseAtomsAdaptor.get_atoms(slab)
        atoms.set_pbc((True, True, True))
        normal = np.asarray(slab.normal, dtype=float)
        normal /= np.linalg.norm(normal)
        proj = slab_projection(atoms, normal)
        lo, hi = float(proj.min()), float(proj.max())
        low_ids, low_comp = terminal_plane(atoms, normal, lo)
        high_ids, high_comp = terminal_plane(atoms, normal, hi)
        if low_comp != {"Ti": 1} or high_comp != {"Ti": 1}:
            continue
        sym_error, pair_map = inversion_map(atoms, normal)
        if sym_error > 1e-5:
            continue
        candidates.append((atoms, slab, normal, sym_error, pair_map, low_ids, high_ids, lo, hi))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one unique symmetric Ti/Ti slab with {EXPECTED_COMPOSITION}; found {len(candidates)}")

    atoms, slab, normal, sym_error, pair_map, low_ids, high_ids, lo, hi = candidates[0]
    atoms.new_array("source_atom_id", np.arange(1, len(atoms) + 1, dtype=np.int64))
    if len(atoms) != sum(EXPECTED_COMPOSITION.values()):
        raise RuntimeError("Unexpected atom count")
    id_pairs = [
        (a, b)
        for a, b, _error in pair_map
        if (a - 1) in low_ids and (b - 1) in high_ids
    ]
    if len(id_pairs) != 1:
        # Match the single inversion-related outermost Ti pair directly.
        id_pairs = []
        frac = atoms.get_scaled_positions(wrap=True)
        for i in low_ids:
            for j in high_ids:
                dxy = (frac[j, :2] + frac[i, :2] + 0.5) % 1.0 - 0.5
                dcart = dxy[0] * atoms.cell[0] + dxy[1] * atoms.cell[1]
                dn = float((atoms.positions[j] - atoms.positions[i]) @ normal)
                if np.linalg.norm(dcart) < 1e-5 and abs(dn - (hi - lo)) < 1e-5:
                    id_pairs.append((int(i + 1), int(j + 1)))
    if len(id_pairs) != 1:
        raise RuntimeError(f"Could not identify a unique inversion-paired terminal Ti pair: {id_pairs}")
    bottom_id, top_id = id_pairs[0]

    composition = dict(sorted(Counter(atoms.get_chemical_symbols()).items()))
    if composition != EXPECTED_COMPOSITION:
        raise RuntimeError(f"Composition mismatch after conversion: {composition}")
    cell_normal_gap = float(np.dot(atoms.cell[2], normal) - (hi - lo))
    if cell_normal_gap < MIN_VACUUM_A - 1e-6:
        raise RuntimeError(f"Insufficient vacuum along surface normal: {cell_normal_gap} A")
    dmin, dmin_ids = min_distance(atoms)
    if dmin < 1.5:
        raise RuntimeError(f"Unphysical short contact: {dmin} A at {dmin_ids}")

    family = "surface_pyramidal_10m11_cod9009647_symm_ti_69atom"
    base_info = {
        "source_family": family,
        "role_assignment": "NONE",
        "candidate_use": "unscored geometry-family DEV candidate only",
        "independence_scope": "new facet family, same COD 9009647 crystallographic root; not source-independent validation",
        "source_cif_path": cif.relative_to(repo).as_posix(),
        "source_cif_sha256": source_hash,
        "source_parser_crosscheck": parser_crosscheck,
        "surface_miller_ase_3index": list(MILLER_INDEX),
        "surface_miller_bravais_equivalent": "(10-11)",
        "surface_termination": "symmetric Ti/Ti",
        "min_slab_thickness_angstrom": MIN_SLAB_THICKNESS_A,
        "minimum_vacuum_along_normal_angstrom": MIN_VACUUM_A,
        "geometry_status": "as-cleaved; unrelaxed; no strain; no DFT labels",
    }
    ref = atoms.copy()
    ref.info.update(base_info)
    ref.info["structure_id"] = "TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF"
    response = ref.copy()
    response.info = dict(ref.info)
    response.info.update({
        "structure_id": "TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_RESP_TIZ020",
        "displacement_definition": "inversion-paired terminal Ti atoms displaced outward along the surface normal",
        "displacement_angstrom": DISPLACEMENT_A,
        "displaced_source_atom_ids": [bottom_id, top_id],
    })
    response.positions[bottom_id - 1] -= DISPLACEMENT_A * normal
    response.positions[top_id - 1] += DISPLACEMENT_A * normal
    resp_error, _ = inversion_map(response, normal)
    if resp_error > 1e-5:
        raise RuntimeError(f"Response does not preserve inversion symmetry: {resp_error}")
    resp_min, resp_min_ids = min_distance(response)
    if resp_min < 1.5:
        raise RuntimeError(f"Response creates unphysical contact: {resp_min} A at {resp_min_ids}")

    out = root / "inputs"
    out.mkdir(exist_ok=True)
    paths = [
        out / f"{ref.info['structure_id']}.extxyz",
        out / f"{response.info['structure_id']}.extxyz",
    ]
    write(paths[0], ref, format="extxyz")
    write(paths[1], response, format="extxyz")

    serialized = []
    for expected, path in zip((ref, response), paths):
        loaded = read(path, format="extxyz")
        if len(loaded) != len(expected) or loaded.get_chemical_symbols() != expected.get_chemical_symbols():
            raise RuntimeError(f"Serialized structure identity changed: {path.name}")
        if not np.array_equal(loaded.pbc, expected.pbc):
            raise RuntimeError(f"Serialized PBC changed: {path.name}")
        if not np.allclose(loaded.cell.array, expected.cell.array, atol=1e-12, rtol=0):
            raise RuntimeError(f"Serialized cell changed: {path.name}")
        if not np.array_equal(loaded.arrays["source_atom_id"], expected.arrays["source_atom_id"]):
            raise RuntimeError(f"Serialized source_atom_id changed: {path.name}")
        if coordinate_sha(loaded) != coordinate_sha(expected):
            raise RuntimeError(f"Serialized coordinate SHA changed: {path.name}")
        serialized.append(loaded)
    serialized_delta = np.linalg.norm(serialized[1].positions - serialized[0].positions, axis=1)
    expected_delta = np.zeros(len(ref))
    expected_delta[bottom_id - 1] = DISPLACEMENT_A
    expected_delta[top_id - 1] = DISPLACEMENT_A
    if not np.allclose(serialized_delta, expected_delta, atol=1e-7, rtol=0):
        raise RuntimeError("Serialized response differs from reference outside the selected terminal pair")

    entries = []
    for image, path, d, ids, inv_err in [
        (serialized[0], paths[0], dmin, dmin_ids, sym_error),
        (serialized[1], paths[1], resp_min, resp_min_ids, resp_error),
    ]:
        entries.append({
            "structure_id": image.info["structure_id"],
            "role_candidate": "UNASSIGNED_GEOMETRY_HOLDOUT_CANDIDATE",
            "source_family": family,
            "path": path.relative_to(root).as_posix(),
            "file_sha256": sha256(path),
            "coordinate_sha256": coordinate_sha(image),
            "atom_count": len(image),
            "composition": composition,
            "cell_angstrom": np.asarray(image.cell).tolist(),
            "pbc": image.pbc.tolist(),
            "surface_normal_cartesian": normal.tolist(),
            "surface_span_along_normal_angstrom": hi - lo,
            "vacuum_along_normal_angstrom": cell_normal_gap,
            "top_termination": {"Ti": 1},
            "bottom_termination": {"Ti": 1},
            "inversion_symmetry_max_error_angstrom": inv_err,
            "minimum_distance_angstrom": d,
            "minimum_distance_source_atom_ids": ids,
            "displaced_source_atom_ids": image.info.get("displaced_source_atom_ids", []),
            "dft_label_present": False,
            "role_assignment_changed": False,
            "serialized_readback_pass": True,
        })

    manifest = {
        "status": "GEOMETRY_AUDIT_PASS_NO_DFT_NO_ROLE_ASSIGNMENT_SAME_SOURCE",
        "generator": Path(__file__).name,
        "generator_sha256": sha256(Path(__file__).resolve()),
        "pymatgen_version": package_version("pymatgen"),
        "ase_version": ase.__version__,
        "source_cif_sha256": source_hash,
        "source_parser_crosscheck": parser_crosscheck,
        "generation": {
            "miller_ase_3index": list(MILLER_INDEX),
            "miller_bravais_equivalent": "(10-11)",
            "min_slab_thickness_angstrom": MIN_SLAB_THICKNESS_A,
            "min_vacuum_angstrom": MIN_VACUUM_A,
            "pymatgen_symmetrize": True,
            "selected_composition": EXPECTED_COMPOSITION,
            "selected_termination": "symmetric Ti/Ti",
            "relaxation": False,
            "strain": False,
            "displacement_angstrom": DISPLACEMENT_A,
            "displaced_source_atom_ids": [bottom_id, top_id],
        },
        "family_independence": "a new, previously unscored surface-orientation family, but derived from COD 9009647; cannot be called independent-source blind validation",
        "entries": entries,
    }
    (root / "GEOMETRY_CANDIDATE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=json_default) + "\n")
    cell_lengths = np.linalg.norm(atoms.cell.array, axis=1).tolist()
    report = (
        "# Ti3SiC2 pyramidal `(10-11)` geometry candidate\n\n"
        "Geometry checks passed. This package contains no DFT labels, no potential scoring, and no data-role assignment.\n\n"
        f"The source is COD 9009647 (SHA-256 `{source_hash}`). Pymatgen {package_version('pymatgen')} generated a symmetric Ti/Ti slab using ASE Miller indices `(1,0,1)`, equivalent to hexagonal `(10-11)`, with minimum slab thickness {MIN_SLAB_THICKNESS_A} A and requested minimum vacuum {MIN_VACUUM_A} A. ASE {ase.__version__} wrote the extxyz files.\n\n"
        f"Pymatgen reported normalizing four near-ideal fractional coordinates in the CIF. Direct ASE-versus-Pymatgen parent matching gives a maximum site displacement of {parser_crosscheck['max_matched_parent_site_difference_angstrom']:.3e} A and RMS {parser_crosscheck['rms_matched_parent_site_difference_angstrom']:.3e} A; this is recorded as a parser-level difference, not a new source.\n\n"
        f"The selected parent contains {len(ref)} atoms with composition Ti33Si12C24. Its cell lengths are `{cell_lengths}` A, PBC is `[true,true,true]`, actual surface-normal span is {hi-lo:.6f} A (requested minimum was {MIN_SLAB_THICKNESS_A:.1f} A), and periodic vacuum projected along the surface normal is {cell_normal_gap:.6f} A. This non-bulk-stoichiometric symmetric termination is an as-cleaved neutral-atom slab; no surface energy or thermodynamic stability is claimed.\n\n"
        f"The response geometry moves the inversion-paired terminal Ti source IDs {bottom_id} and {top_id} outward by {DISPLACEMENT_A} A along the surface normal. Inversion mapping errors are {sym_error:.3e} A (reference) and {resp_error:.3e} A (response). Minimum distances are {dmin:.6f} A (reference) and {resp_min:.6f} A (response).\n\n"
        "This is a new, not-yet-scored geometry family, but it shares the COD 9009647 crystallographic root with existing datasets. It is only a candidate for facet-held-out development diagnostics, not an independent-source blind test. Both candidates are unassigned; TEST remains outside E's package.\n"
    )
    (root / "GEOMETRY_AUDIT.md").write_text(report)
    print(json.dumps({"status": manifest["status"], "entries": entries}, indent=2, default=json_default))


if __name__ == "__main__":
    main()
