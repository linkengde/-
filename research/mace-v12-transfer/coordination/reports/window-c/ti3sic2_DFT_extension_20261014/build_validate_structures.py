#!/usr/bin/env python3
"""Build and validate the ten source-backed Stage C candidate geometries.

This script only edits extxyz input geometry files and a geometry audit JSON.
It never constructs a calculator, launches GPAW, evaluates a potential, or
changes any source/archive file. Run from the repository root with the GPAW
environment's Python (ASE 3.29.0).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write


ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
SOURCE = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "pbe_interface_v19_pure_phase_controls/inputs"
)
BASELINE = SOURCE / "Ti3SiC2_baseline_PROPOSAL.extxyz"
EXPECTED_BASELINE_SHA = "157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4"
EXPECTED_COUNTS = {"Ti": 24, "Si": 8, "C": 16}
DISPLACEMENT_A = 0.02
VACUUM_A = 20.0
SLAB_REPEAT_C = 2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def minimum_distance(atoms: Atoms) -> tuple[float, int, int]:
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    i, j = np.unravel_index(np.argmin(distances), distances.shape)
    return float(distances[i, j]), int(i), int(j)


def formula_counts(atoms: Atoms) -> dict[str, int]:
    return {s: int(sum(atoms.get_chemical_symbols().count(s) for _ in [0]))
            for s in ("Ag", "Ti", "Si", "C")
            if s in atoms.get_chemical_symbols()}


def validate(atoms: Atoms, label: str) -> dict:
    dmin, i, j = minimum_distance(atoms)
    counts = formula_counts(atoms)
    assert len(atoms) > 0, label
    assert counts.get("Ag", 0) == 0, (label, counts)
    assert dmin >= 1.5, (label, dmin, i, j)
    assert np.isfinite(atoms.positions).all(), label
    assert np.isfinite(atoms.cell.array).all(), label
    assert abs(float(atoms.get_volume())) > 1.0, label
    assert atoms.pbc.tolist() == [True, True, True], (label, atoms.pbc)
    return {
        "id": label,
        "natoms": len(atoms),
        "composition": counts,
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "minimum_distance_A": dmin,
        "minimum_distance_pair_1based": [i + 1, j + 1],
        "minimum_distance_species": [atoms[i].symbol, atoms[j].symbol],
        "geometry_pass": True,
    }


def slab_from_parent(base: Atoms, termination: str) -> tuple[Atoms, dict]:
    c = float(np.linalg.norm(base.cell[2]))
    base_pos = base.positions
    lo = 0.0 if termination == "Ti" else 4.405875
    hi = lo + SLAB_REPEAT_C * c
    tol = 1e-7
    local_ids = base.arrays["local_new_id"].astype(int)
    site_kinds = base.arrays["spacegroup_kinds"].astype(int)
    source_rows: list[tuple[float, int, int, int, str, np.ndarray, int]] = []
    for atom_idx, atom in enumerate(base):
        for n in range(-3, 5):
            z = float(base_pos[atom_idx, 2] + n * c)
            if lo - tol <= z <= hi + tol:
                source_rows.append((z, int(local_ids[atom_idx]), n, atom_idx,
                                    atom.symbol, base_pos[atom_idx].copy(),
                                    int(site_kinds[atom_idx])))
    # Stable atom order: height, parent-local ID, then periodic image index.
    source_rows.sort(key=lambda row: (round(row[0], 8), row[1], row[2]))
    expected_planes = 2 * SLAB_REPEAT_C * 6 + 1
    assert len(source_rows) == 100, (termination, len(source_rows))
    assert len({round(row[0], 6) for row in source_rows}) == expected_planes

    cell = base.cell.array.copy()
    cell[2] = np.array([0.0, 0.0, SLAB_REPEAT_C * c + VACUUM_A])
    positions = []
    symbols = []
    parent_ids = []
    repeat_indices = []
    kinds = []
    for z, parent_id, n, atom_idx, symbol, parent_pos, kind in source_rows:
        p = parent_pos.copy()
        p[2] = z - lo + VACUUM_A / 2.0
        positions.append(p)
        symbols.append(symbol)
        parent_ids.append(parent_id)
        repeat_indices.append(n)
        kinds.append(kind)

    atoms = Atoms(symbols=symbols, positions=np.asarray(positions), cell=cell,
                  pbc=[True, True, True])
    atoms.new_array("local_new_id", np.arange(1, len(atoms) + 1, dtype=int))
    atoms.new_array("parent_local_new_id", np.asarray(parent_ids, dtype=int))
    atoms.new_array("parent_repeat_index_c", np.asarray(repeat_indices, dtype=int))
    atoms.new_array("spacegroup_kinds", np.asarray(kinds, dtype=int))
    top_z = float(max(atoms.positions[:, 2]))
    top_candidates = [i for i, z in enumerate(atoms.positions[:, 2])
                      if abs(float(z) - top_z) < 1e-7]
    # Deterministic top-layer site chosen by parent-local ID then x/y.
    chosen = min(top_candidates,
                 key=lambda i: (int(atoms.arrays["parent_local_new_id"][i]),
                                float(atoms.positions[i, 0]),
                                float(atoms.positions[i, 1])))
    atoms.info.update({
        "phase": "Ti3SiC2_clean_symmetric_basal_slab",
        "source_baseline_sha256": EXPECTED_BASELINE_SHA,
        "source_parent_path": str(BASELINE),
        "termination": termination,
        "slab_repeat_c": SLAB_REPEAT_C,
        "vacuum_A": VACUUM_A,
        "geometry_role": "TRAIN_CANDIDATE_PENDING_ROLE_APPROVAL",
        "surface_response_atom_id_1based": chosen + 1,
        "surface_response_parent_local_id": int(atoms.arrays["parent_local_new_id"][chosen]),
        "surface_response_repeat_index_c": int(atoms.arrays["parent_repeat_index_c"][chosen]),
    })
    meta = {
        "termination": termination,
        "source_plane_z_A": lo,
        "material_thickness_A": SLAB_REPEAT_C * c,
        "vacuum_A": VACUUM_A,
        "cell_z_A": float(cell[2, 2]),
        "plane_count": expected_planes,
        "atom_count": len(atoms),
        "top_layer_response_atom_1based": chosen + 1,
        "top_layer_parent_local_id": int(atoms.arrays["parent_local_new_id"][chosen]),
    }
    return atoms, meta


def main() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    assert sha256(BASELINE) == EXPECTED_BASELINE_SHA
    baseline = read(BASELINE, format="extxyz")
    assert len(baseline) == 48
    assert formula_counts(baseline) == EXPECTED_COUNTS
    assert baseline.pbc.tolist() == [True, True, True]
    assert "local_new_id" in baseline.arrays
    assert sorted(baseline.arrays["local_new_id"].astype(int).tolist()) == list(range(1, 49))

    source_files = {
        "A_VOL_M02": SOURCE / "Ti3SiC2_volume_m02_PROPOSAL.extxyz",
        "A_VOL_P02": SOURCE / "Ti3SiC2_volume_p02_PROPOSAL.extxyz",
        "A_SHEAR_M005": SOURCE / "Ti3SiC2_shear_m005_PROPOSAL.extxyz",
        "A_SHEAR_P005": SOURCE / "Ti3SiC2_shear_p005_PROPOSAL.extxyz",
    }
    expected_source_sha = {
        "A_VOL_M02": "6be6f7f69262d9c71efa7c7cb3ab0f3b099e44cac9f5611b9659e8336394596e",
        "A_VOL_P02": "e24ef50128af4dfa20fca0964bc7d29d5d2f247d33ae98968ba92055e1f64d88",
        "A_SHEAR_M005": "9b715a257bf9dcde0f304d6c7f70e3a66326e0a3c57fc2a3f36706c73f2a60af",
        "A_SHEAR_P005": "889967da2059ec56b0f62ff5252214228be1bb7730c878a970685018c37be8de",
    }
    records: list[dict] = []
    for key, source in source_files.items():
        assert sha256(source) == expected_source_sha[key], key
        atoms = read(source, format="extxyz")
        assert len(atoms) == 48 and formula_counts(atoms) == EXPECTED_COUNTS
        out = INPUTS / f"{key}.extxyz"
        # Preserve source bytes exactly for the four existing proposal geometries.
        out.write_bytes(source.read_bytes())
        item = validate(read(out, format="extxyz"), key)
        item.update({
            "file": str(out.relative_to(ROOT)),
            "file_sha256": sha256(out),
            "source_file": str(source),
            "source_sha256": sha256(source),
        })
        records.append(item)

    for sign, suffix in [(-1, "M020"), (1, "P020")]:
        atoms = baseline.copy()
        idx = int(np.where(atoms.arrays["local_new_id"] == 7)[0][0])
        assert atoms[idx].symbol == "Si"
        atoms.positions[idx, 2] += sign * DISPLACEMENT_A
        atoms.info.update({
            "case": f"Si_local_id0007_z_{suffix}",
            "source_baseline_sha256": EXPECTED_BASELINE_SHA,
            "source_parent_path": str(BASELINE),
            "distortion": "single-Si fixed-cell z displacement",
            "displacement_A": sign * DISPLACEMENT_A,
            "displaced_local_new_id": 7,
            "geometry_role": "TRAIN_CANDIDATE_PENDING_ROLE_APPROVAL",
        })
        out = INPUTS / f"B_SI_ID0007_Z_{suffix}.extxyz"
        write(out, atoms, format="extxyz")
        item = validate(read(out, format="extxyz"), f"B_SI_ID0007_Z_{suffix}")
        item.update({
            "file": str(out.relative_to(ROOT)),
            "file_sha256": sha256(out),
            "source_file": str(BASELINE),
            "source_sha256": EXPECTED_BASELINE_SHA,
            "displaced_local_new_id": 7,
            "displacement_A": sign * DISPLACEMENT_A,
        })
        records.append(item)

    slab_meta = []
    for termination in ("Ti", "Si"):
        ref, meta = slab_from_parent(baseline, termination)
        ref.info["case"] = f"C_{termination.upper()}_TERM_REF"
        ref_path = INPUTS / f"C_{termination.upper()}_TERM_REF.extxyz"
        write(ref_path, ref, format="extxyz")
        item = validate(read(ref_path, format="extxyz"), f"C_{termination.upper()}_TERM_REF")
        item.update({
            "file": str(ref_path.relative_to(ROOT)),
            "file_sha256": sha256(ref_path),
            "source_file": str(BASELINE),
            "source_sha256": EXPECTED_BASELINE_SHA,
            "termination": termination,
            **meta,
        })
        records.append(item)

        resp = ref.copy()
        response_id = int(ref.info["surface_response_atom_id_1based"])
        ridx = response_id - 1
        resp.positions[ridx, 2] += DISPLACEMENT_A
        resp.info.update({
            "case": f"C_{termination.upper()}_TERM_TOP_ID{response_id:03d}_Z_P020",
            "distortion": "one top-layer atom fixed-cell outward z displacement",
            "displacement_A": DISPLACEMENT_A,
            "geometry_role": "TRAIN_CANDIDATE_PENDING_ROLE_APPROVAL",
        })
        resp_path = INPUTS / f"C_{termination.upper()}_TERM_TOP_Z_P020.extxyz"
        write(resp_path, resp, format="extxyz")
        item = validate(read(resp_path, format="extxyz"), f"C_{termination.upper()}_TERM_TOP_Z_P020")
        item.update({
            "file": str(resp_path.relative_to(ROOT)),
            "file_sha256": sha256(resp_path),
            "source_file": str(BASELINE),
            "source_sha256": EXPECTED_BASELINE_SHA,
            "termination": termination,
            "response_atom_1based": response_id,
            "displacement_A": DISPLACEMENT_A,
            **meta,
        })
        records.append(item)
        slab_meta.append(meta)

    assert len(records) == 10, len(records)
    report = {
        "scope": "geometry generation and static preflight only; no calculator/DFT/MD invoked",
        "ase_version": __import__("ase").__version__,
        "baseline_source": str(BASELINE),
        "baseline_source_sha256": EXPECTED_BASELINE_SHA,
        "bulk_composition": EXPECTED_COUNTS,
        "bulk_cell_A": baseline.cell.array.tolist(),
        "bulk_points": 6,
        "clean_slab_points": 4,
        "clean_slab_construction": {
            "cleavage": "symmetric basal slabs cut between periodic planes, retaining identical termination at both faces",
            "parent_c_repeat_count": SLAB_REPEAT_C,
            "vacuum_A": VACUUM_A,
            "pbc": [True, True, True],
            "slab_layers": slab_meta,
            "warning": "Ti- and Si-terminated slabs have different stoichiometries; their raw total energies are not directly comparable surface energies without chemical-potential/reference corrections.",
        },
        "structures": records,
        "all_source_backed_geometries_validated": all(x["geometry_pass"] for x in records),
        "independent_dev_test_structures": "not created: independent parent-family source evidence is not present in the reviewed E6 request/source bundle",
    }
    (ROOT / "geometry_validation.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
