#!/usr/bin/env python3
"""Create and geometry-screen label-free V17 validation/test proposals.

This reader extracts only symbols, positions, IDs, cell and PBC from extxyz.
It does not interpret energy/force fields and never traverses calculation folders.
"""
import argparse
import hashlib
import json
import math
import re
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from scipy.optimize import linear_sum_assignment

PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
WINDOW = Path("research/mace-v12-transfer/coordination/reports/window-b")
OUTPUT_REL = WINDOW / "v17_split_geometry_design"
NEAR_AG_RMS_A = 0.15
NEAR_FRAMEWORK_RMS_A = 0.15
EXACT_MAX_A = 1e-5
CELL_TOL_A = 1e-5
PAIR_FLOOR_ALLOWANCE_A = 0.05
TRANSVERSE_SHIFT_A = [0.25, 0.25, 0.0]
VALIDATION_SHIFTS_A = {
    "AgC": [-0.25, 0.35, 0.0],
    "AgSi": [-0.25, 0.35, 0.0],
    "AgTi": [0.35, -0.25, 0.0],
}
LOCAL_FRAMEWORK_SHIFT_A = 0.15
IFACES = ["AgC", "AgSi", "AgTi"]
PARENT_INPUTS = {
    "AgC": PROJECT / "pbe_interface_v14_parallel_acquisition/inputs/AgC_registry_strain_acq_v14_01.extxyz",
    "AgSi": PROJECT / "pbe_interface_v14_parallel_acquisition/inputs/AgSi_registry_strain_acq_v14_01.extxyz",
    "AgTi": PROJECT / "pbe_interface_v14_parallel_acquisition/inputs/AgTi_registry_probe_v13_01.extxyz",
}
WIDE_TESTS = {
    "AgC": "AgC_wide_distance_m0p20",
    "AgSi": "AgSi_wide_distance_p0p35",
    "AgTi": "AgTi_wide_distance_p0p20",
}
V16_CANDIDATE_DIR = WINDOW / "v16_unseen_validation_geometry_design/final_candidates"
V15_WIDE_DIR = WINDOW / "v15_wide_span_sampling_design"
V17_B_DIR = PROJECT / "pbe_interface_v17_parallel_targeted_acquisition"


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def geometry_sha(atoms):
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(atoms.positions, 10).tolist(),
        "cell_A": np.round(atoms.cell.array, 10).tolist(),
        "pbc": atoms.pbc.astype(bool).tolist(),
        "lammps_id": atoms.arrays.get("lammps_id", np.array([], dtype=int)).tolist(),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def extxyz_geometry_frames(path):
    """Read only identity and geometry fields from one or more extxyz frames."""
    result = []
    with Path(path).open("r", encoding="utf-8") as f:
        frame = 0
        while True:
            first = f.readline()
            if not first:
                break
            if not first.strip():
                continue
            n = int(first.strip())
            header = f.readline().rstrip("\n")
            match = re.search(r"Properties=([^ ]+)", header)
            need(match is not None, f"missing extxyz Properties in {path} frame {frame}")
            schema = match.group(1).split(":")
            fields = {}
            offset = 0
            for i in range(0, len(schema), 3):
                name, kind, width = schema[i], schema[i + 1], int(schema[i + 2])
                fields[name] = (offset, width, kind)
                offset += width
            need("species" in fields and "pos" in fields, f"missing geometry columns in {path} frame {frame}")
            lattice = re.search(r'Lattice="([^"]+)"', header)
            need(lattice is not None, f"missing Lattice in {path} frame {frame}")
            cell = np.asarray([float(x) for x in lattice.group(1).split()], dtype=float).reshape(3, 3)
            pbcm = re.search(r'pbc="([^"]+)"', header)
            pbc = [x.upper() in {"T", "TRUE", "1"} for x in pbcm.group(1).split()] if pbcm else [True, True, True]
            symbols, positions, ids, types, marked = [], [], [], [], []
            for _ in range(n):
                row = f.readline().split()
                s0, _, _ = fields["species"]
                p0, _, _ = fields["pos"]
                symbols.append(row[s0])
                positions.append([float(x) for x in row[p0:p0 + 3]])
                if "lammps_id" in fields:
                    q, _, _ = fields["lammps_id"]
                    ids.append(int(row[q]))
                if "lammps_type" in fields:
                    q, _, _ = fields["lammps_type"]
                    types.append(int(row[q]))
                if "central_pair" in fields:
                    q, _, _ = fields["central_pair"]
                    marked.append(int(row[q]))
            atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)
            if ids:
                atoms.new_array("lammps_id", np.asarray(ids, dtype=int))
            if types:
                atoms.new_array("lammps_type", np.asarray(types, dtype=int))
            if marked:
                atoms.new_array("central_pair", np.asarray(marked, dtype=int))
            result.append(atoms)
            frame += 1
    return result


def formula_key(atoms):
    return tuple(sorted(Counter(atoms.get_chemical_symbols()).items()))


def atom_id_map(atoms):
    ids = atoms.arrays.get("lammps_id")
    if ids is None or len(ids) != len(atoms) or len(set(map(int, ids))) != len(atoms):
        return None
    return {int(value): i for i, value in enumerate(ids)}


def compare_geometry(a, b):
    if len(a) != len(b) or formula_key(a) != formula_key(b) or not np.array_equal(a.pbc, b.pbc):
        return None
    cell_delta = float(np.max(np.abs(a.cell.array - b.cell.array)))
    if cell_delta > CELL_TOL_A:
        return None
    amap, bmap = atom_id_map(a), atom_id_map(b)
    if amap is not None and bmap is not None and set(amap) == set(bmap):
        ids = sorted(amap)
        ia = np.asarray([amap[x] for x in ids], dtype=int)
        ib = np.asarray([bmap[x] for x in ids], dtype=int)
        method = "persistent_lammps_id"
    elif a.get_chemical_symbols() == b.get_chemical_symbols():
        ia = np.arange(len(a), dtype=int)
        ib = np.arange(len(b), dtype=int)
        method = "same_element_order"
    else:
        pairs = []
        aa = a.positions - np.mean(a.positions, axis=0)
        bb = b.positions - np.mean(b.positions, axis=0)
        for element in sorted(set(a.get_chemical_symbols())):
            ai = np.asarray([i for i, x in enumerate(a.get_chemical_symbols()) if x == element], dtype=int)
            bi = np.asarray([i for i, x in enumerate(b.get_chemical_symbols()) if x == element], dtype=int)
            _, distances = find_mic((aa[ai, None, :] - bb[None, bi, :]).reshape(-1, 3), a.cell, pbc=a.pbc)
            rows, cols = linear_sum_assignment(np.asarray(distances).reshape(len(ai), len(bi)))
            pairs.extend((int(ai[r]), int(bi[c])) for r, c in zip(rows, cols))
        ia = np.asarray([x[0] for x in pairs], dtype=int)
        ib = np.asarray([x[1] for x in pairs], dtype=int)
        method = "species_assignment"
    delta, _ = find_mic(a.positions[ia] - b.positions[ib], a.cell, pbc=a.pbc)
    symbols = np.asarray(a.get_chemical_symbols())[ia]
    framework = symbols != "Ag"
    ag = symbols == "Ag"
    common = np.median(delta[framework], axis=0) if framework.any() else np.median(delta, axis=0)
    delta = delta - common
    magnitudes = np.linalg.norm(delta, axis=1)
    ag_rms = float(np.sqrt(np.mean(magnitudes[ag] ** 2))) if ag.any() else 0.0
    fw_rms = float(np.sqrt(np.mean(magnitudes[framework] ** 2))) if framework.any() else 0.0
    fw_max = float(np.max(magnitudes[framework])) if framework.any() else 0.0
    all_rms = float(np.sqrt(np.mean(magnitudes ** 2)))
    maximum = float(np.max(magnitudes))
    exact = bool(cell_delta <= CELL_TOL_A and maximum <= EXACT_MAX_A)
    near = bool(ag_rms <= NEAR_AG_RMS_A and fw_rms <= NEAR_FRAMEWORK_RMS_A)
    return {
        "matching_method": method,
        "cell_max_delta_A": cell_delta,
        "Ag_registry_RMS_A": ag_rms,
        "framework_RMS_A": fw_rms,
        "framework_max_displacement_A": fw_max,
        "all_atom_RMS_A": all_rms,
        "max_atom_displacement_A": maximum,
        "exact_duplicate": exact,
        "near_duplicate": near,
    }


def pair_minima(atoms):
    distances = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    ids = atoms.arrays.get("lammps_id", np.arange(len(atoms), dtype=int))
    minima = {}
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            pair = "-".join(sorted((symbols[i], symbols[j])))
            value = float(distances[i, j])
            if pair not in minima or value < minima[pair]["distance_A"]:
                minima[pair] = {
                    "distance_A": value,
                    "atom_ids": [int(ids[i]), int(ids[j])],
                    "elements_in_file_order": [symbols[i], symbols[j]],
                }
    return minima


def marked_pair(atoms):
    need("central_pair" in atoms.arrays and "lammps_id" in atoms.arrays, "candidate lacks central_pair or persistent IDs")
    marked = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=int) > 0)
    need(len(marked) == 2, "expected exactly two marked central-pair atoms")
    ag = [int(i) for i in marked if atoms[i].symbol == "Ag"]
    other = [int(i) for i in marked if atoms[i].symbol != "Ag"]
    need(len(ag) == len(other) == 1, "central pair must have one Ag and one framework atom")
    return ag[0], other[0]


def write_geometry(path, atoms, metadata):
    ids = atoms.arrays.get("lammps_id")
    types = atoms.arrays.get("lammps_type")
    marks = atoms.arrays.get("central_pair")
    need(ids is not None and types is not None and marks is not None, "candidate must preserve IDs, atom types and central-pair marks")
    cell = " ".join(f"{x:.12f}" for x in atoms.cell.array.reshape(-1))
    pbc = " ".join("T" if x else "F" for x in atoms.pbc)
    fields = [
        f'config_type="{metadata["candidate_id"]}"',
        f'candidate_role="{metadata["role"]}"',
        f'source_path="{metadata["direct_source_path"]}"',
        f'source_sha256={metadata["direct_source_sha256"]}',
        f'parent_input_path="{metadata["parent_input_path"]}"',
        f'parent_input_sha256={metadata["parent_input_sha256"]}',
        "proposal_only_no_DFT_label=T",
    ]
    header = (
        f'{len(atoms)}\nLattice="{cell}" Properties=species:S:1:pos:R:3:lammps_id:I:1:lammps_type:I:1:central_pair:I:1 '
        + " ".join(fields)
        + f' pbc="{pbc}"\n'
    )
    with Path(path).open("w", encoding="utf-8", newline="\n") as f:
        f.write(header)
        for i, atom in enumerate(atoms):
            x, y, z = atom.position
            f.write(
                f"{atom.symbol} {x:.12f} {y:.12f} {z:.12f} "
                f"{int(ids[i])} {int(types[i])} {int(marks[i])}\n"
            )


def nonzero_displacements(parent, candidate):
    pids, cids = atom_id_map(parent), atom_id_map(candidate)
    need(pids is not None and cids is not None and set(pids) == set(cids), "source/candidate persistent ID mismatch")
    out = []
    for atom_id in sorted(pids):
        i, j = pids[atom_id], cids[atom_id]
        delta, _ = find_mic(candidate.positions[j] - parent.positions[i], parent.cell, pbc=parent.pbc)
        if float(np.linalg.norm(delta)) > 1e-8:
            out.append({
                "id": int(atom_id),
                "element": parent[i].symbol,
                "delta_xyz_A": [float(x) for x in delta],
                "magnitude_A": float(np.linalg.norm(delta)),
            })
    return out


def inventory_geometry_files(project):
    split_files = sorted(
        p for p in project.glob("mace_periodic_v*/data/*.extxyz")
        if p.name in {"train.extxyz", "valid.extxyz", "test.extxyz"}
    )
    input_files = sorted(project.glob("pbe_*/inputs/**/*.extxyz"))
    return split_files, input_files


def create_candidates(repo, output):
    project = repo / PROJECT
    v16_manifest_path = repo / V16_CANDIDATE_DIR / "candidate_manifest.json"
    v15_manifest_path = repo / V15_WIDE_DIR / "candidate_manifest.json"
    v17_manifest_path = repo / V17_B_DIR / "input_manifest.json"
    v16_manifest = json.loads(v16_manifest_path.read_text())
    v15_manifest = json.loads(v15_manifest_path.read_text())
    v17_manifest = json.loads(v17_manifest_path.read_text())
    v16_by_iface = {row["interface"]: row for row in v16_manifest["records"]}
    v15_by_id = {row["candidate_id"]: row for row in v15_manifest["candidates"]}
    v17_by_label = {row["label"]: row for row in v17_manifest["records"]}

    split_files, input_files = inventory_geometry_files(project)
    need(len(split_files) > 0 and len(input_files) > 0, "missing historical geometry inventory")
    inventory = []
    source_file_records = []
    split_summary = []
    for path in split_files:
        frames = extxyz_geometry_frames(path)
        digest = sha256(path)
        source_file_records.append({
            "path": str(path.relative_to(repo)),
            "sha256": digest,
            "role": "historical_model_split_geometry_only",
            "frames": len(frames),
        })
        split_summary.append({
            "path": str(path.relative_to(repo)),
            "sha256": digest,
            "frames": len(frames),
        })
        inventory.extend({
            "path": str(path.relative_to(repo)),
            "frame_index_zero_based": i,
            "sha256": digest,
            "role": "historical_model_split_geometry_only",
            "atoms": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "geometry_sha256": geometry_sha(atoms),
            "atoms_object": atoms,
        } for i, atoms in enumerate(frames))
    for path in input_files:
        frames = extxyz_geometry_frames(path)
        digest = sha256(path)
        source_file_records.append({
            "path": str(path.relative_to(repo)),
            "sha256": digest,
            "role": "published_dft_input_geometry_only",
            "frames": len(frames),
        })
        inventory.extend({
            "path": str(path.relative_to(repo)),
            "frame_index_zero_based": i,
            "sha256": digest,
            "role": "published_dft_input_geometry_only",
            "atoms": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "geometry_sha256": geometry_sha(atoms),
            "atoms_object": atoms,
        } for i, atoms in enumerate(frames))

    v17_inputs = []
    for path in (project / "pbe_interface_v17_main_targeted_acquisition/inputs").glob("*.extxyz"):
        v17_inputs.append(path)
    for path in (project / "pbe_interface_v17_parallel_targeted_acquisition/inputs").glob("*.extxyz"):
        v17_inputs.append(path)
    need(len(set(v17_inputs)) == 8, f"expected eight planned V17 training inputs; found {len(set(v17_inputs))}")
    v17_input_checks = []
    for label, row in v17_by_label.items():
        path = repo / V17_B_DIR / row["input"]
        need(sha256(path) == row["input_sha256"], f"V17 B input hash mismatch: {label}")
        v17_input_checks.append({"label": label, "path": str(path.relative_to(repo)), "sha256": row["input_sha256"]})
    v17_main_manifest_path = repo / PROJECT / "pbe_interface_v17_main_targeted_acquisition/input_manifest.json"
    v17_main_manifest = json.loads(v17_main_manifest_path.read_text())
    for row in v17_main_manifest["records"]:
        path = repo / PROJECT / "pbe_interface_v17_main_targeted_acquisition" / row["input"]
        need(sha256(path) == row["input_sha256"], f"V17 A input hash mismatch: {row['label']}")
        v17_input_checks.append({"label": row["label"], "path": str(path.relative_to(repo)), "sha256": row["input_sha256"]})
    need(len(v17_input_checks) == 8, "V17 manifest set does not have eight distinct planned inputs")

    selected = []
    for interface in IFACES:
        parent_path = repo / PARENT_INPUTS[interface]
        parent_sha = sha256(parent_path)
        parent = extxyz_geometry_frames(parent_path)[0]
        source_label = {
            "AgC": "AgC_framework_neighbor_m_acq_v17_01",
            "AgSi": "AgSi_framework_neighbor_m_acq_v17_01",
            "AgTi": "AgTi_framework_neighbor_m_acq_v17_01",
        }[interface]
        parent_record = v17_by_label[source_label]
        ag_index, framework_index = marked_pair(parent)
        outward, _ = find_mic(parent.positions[framework_index] - parent.positions[ag_index], parent.cell, pbc=parent.pbc)
        outward = outward / np.linalg.norm(outward)

        # Role 1: fresh development-validation geometry. The previous V16 proposals
        # were rejected because they exactly matched the frozen V16 blind inputs.
        validation = parent.copy()
        validation_ag = np.flatnonzero(np.asarray(validation.get_chemical_symbols()) == "Ag")
        validation_shift = np.asarray(VALIDATION_SHIFTS_A[interface], dtype=float)
        validation.positions[validation_ag] += validation_shift
        validation.positions[framework_index] += outward * LOCAL_FRAMEWORK_SHIFT_A
        selected.append({
            "interface": interface,
            "candidate_id": f"{interface}_v17_dev_validation_01",
            "role": "development_validation",
            "atoms": validation,
            "parent_atoms": parent,
            "direct_source_path": str(parent_path.relative_to(repo)),
            "direct_source_sha256": parent_sha,
            "parent_input_path": str(parent_path.relative_to(repo)),
            "parent_input_sha256": parent_sha,
            "lineage": {
                "source_v17_diagnostic_parent_record": source_label,
                "source_dataset_path": parent_record["source_parent"]["dataset_path"],
                "source_dataset_sha256": parent_record["source_parent"]["dataset_sha256"],
                "source_frame_index_zero_based": parent_record["source_parent"]["frame_index_zero_based"],
                "source_config_type": parent_record["source_parent"]["config_type"],
                "source_frame_geometry_sha256": parent_record["source_parent"]["frame_geometry_sha256"],
                "Ag_sublattice_transverse_translation_A": VALIDATION_SHIFTS_A[interface],
                "local_framework_atom_id": int(validation.arrays["lammps_id"][framework_index]),
                "local_framework_element": validation[framework_index].symbol,
                "local_framework_displacement_from_marked_Ag_A": [float(x) for x in outward * LOCAL_FRAMEWORK_SHIFT_A],
                "local_framework_displacement_magnitude_A": LOCAL_FRAMEWORK_SHIFT_A,
                "fresh_geometry_generated_from_parent_input": True,
            },
        })

        source_label = {
            "AgC": "AgC_framework_neighbor_m_acq_v17_01",
            "AgSi": "AgSi_framework_neighbor_m_acq_v17_01",
            "AgTi": "AgTi_framework_neighbor_m_acq_v17_01",
        }[interface]
        v17row = v17_by_label[source_label]
        combo = parent.copy()
        ag_index, framework_index = marked_pair(combo)
        delta_vec = combo.positions[framework_index] - combo.positions[ag_index]
        outward, _ = find_mic(delta_vec, combo.cell, pbc=combo.pbc)
        outward = outward / np.linalg.norm(outward)
        ag_indices = np.flatnonzero(np.asarray(combo.get_chemical_symbols()) == "Ag")
        combo.positions[ag_indices] += np.asarray(TRANSVERSE_SHIFT_A, dtype=float)
        combo.positions[framework_index] += outward * LOCAL_FRAMEWORK_SHIFT_A
        selected.append({
            "interface": interface,
            "candidate_id": f"{interface}_v17_test_transverse_framework_01",
            "role": "withheld_test",
            "atoms": combo,
            "parent_atoms": parent,
            "direct_source_path": str(parent_path.relative_to(repo)),
            "direct_source_sha256": parent_sha,
            "parent_input_path": str(parent_path.relative_to(repo)),
            "parent_input_sha256": parent_sha,
            "lineage": {
                "source_v17_diagnostic_parent_record": source_label,
                "source_dataset_path": v17row["source_parent"]["dataset_path"],
                "source_dataset_sha256": v17row["source_parent"]["dataset_sha256"],
                "source_frame_index_zero_based": v17row["source_parent"]["frame_index_zero_based"],
                "source_config_type": v17row["source_parent"]["config_type"],
                "source_frame_geometry_sha256": v17row["source_parent"]["frame_geometry_sha256"],
                "Ag_sublattice_transverse_translation_A": TRANSVERSE_SHIFT_A,
                "local_framework_atom_id": int(combo.arrays["lammps_id"][framework_index]),
                "local_framework_element": combo[framework_index].symbol,
                "local_framework_displacement_from_marked_Ag_A": [float(x) for x in outward * LOCAL_FRAMEWORK_SHIFT_A],
                "local_framework_displacement_magnitude_A": LOCAL_FRAMEWORK_SHIFT_A,
                "candidate_uses_only_parent_geometry_and_declared_transforms": True,
            },
        })

        wide_id = WIDE_TESTS[interface]
        wide_row = v15_by_id[wide_id]
        wide_path = repo / wide_row["path"]
        need(sha256(wide_path) == wide_row["sha256"], f"wide-span candidate hash mismatch: {wide_id}")
        wide_atoms = extxyz_geometry_frames(wide_path)[0]
        wide_parent_input_path = repo / PARENT_INPUTS[interface]
        wide_parent = extxyz_geometry_frames(wide_parent_input_path)[0]
        selected.append({
            "interface": interface,
            "candidate_id": f"{interface}_v17_test_distance_01",
            "role": "withheld_test",
            "atoms": wide_atoms,
            "parent_atoms": wide_parent,
            "direct_source_path": str(wide_path.relative_to(repo)),
            "direct_source_sha256": sha256(wide_path),
            "parent_input_path": str(wide_parent_input_path.relative_to(repo)),
            "parent_input_sha256": sha256(wide_parent_input_path),
            "lineage": {
                "reuse_of_prior_unlabeled_wide_span_proposal": True,
                "prior_proposal_id": wide_id,
                "prior_proposal_path": str(wide_path.relative_to(repo)),
                "prior_proposal_sha256": sha256(wide_path),
                "proposal_geometry_sha256": wide_row["candidate_geometry_sha256"],
                "parent_dataset_path": wide_row["parent"]["path"],
                "parent_dataset_sha256": wide_row["parent"]["file_sha256"],
                "parent_frame_index_zero_based": wide_row["parent"]["frame"],
                "parent_config_type": wide_row["parent"]["config_type"],
                "parent_frame_geometry_sha256": wide_row["parent"]["geometry_sha256"],
                "axis": wide_row["axis"],
                "offset_A": wide_row["offset_A"],
                "candidate_method": wide_row["parent"]["source_method"],
            },
        })

    need(len(selected) == 9, f"expected nine role-assigned candidates; got {len(selected)}")
    output.mkdir(parents=True, exist_ok=False)
    candidates_dir = output / "geometries"
    candidates_dir.mkdir()
    for row in selected:
        atoms = row["atoms"]
        parent = row["parent_atoms"]
        need(len(atoms) == len(parent), f"source atom count mismatch: {row['candidate_id']}")
        need(atom_id_map(atoms) is not None and atom_id_map(parent) is not None, f"missing unique IDs: {row['candidate_id']}")
        need(set(atom_id_map(atoms)) == set(atom_id_map(parent)), f"candidate ID set changed: {row['candidate_id']}")
        need(np.array_equal(atoms.pbc, parent.pbc), f"PBC changed: {row['candidate_id']}")
        need(np.max(np.abs(atoms.cell.array - parent.cell.array)) <= 1e-10, f"cell changed: {row['candidate_id']}")
        ag, x = marked_pair(atoms)
        fname = row["candidate_id"] + ".extxyz"
        row["path"] = str((Path("geometries") / fname).as_posix())
        write_geometry(candidates_dir / fname, atoms, row)
        written = extxyz_geometry_frames(candidates_dir / fname)[0]
        need(len(written) == len(atoms), f"written atom count changed: {row['candidate_id']}")
        need(np.array_equal(written.arrays["lammps_id"], atoms.arrays["lammps_id"]), f"written IDs changed: {row['candidate_id']}")
        need(np.array_equal(written.arrays["lammps_type"], atoms.arrays["lammps_type"]), f"written atom types changed: {row['candidate_id']}")
        need(np.array_equal(written.arrays["central_pair"], atoms.arrays["central_pair"]), f"written marker changed: {row['candidate_id']}")
        need(np.max(np.abs(written.cell.array - atoms.cell.array)) <= 1e-10, f"written cell changed: {row['candidate_id']}")
        need(np.array_equal(written.pbc, atoms.pbc), f"written PBC changed: {row['candidate_id']}")
        need(np.max(np.abs(written.positions - atoms.positions)) <= 1e-10, f"written coordinates changed: {row['candidate_id']}")
        row["written_atoms"] = written
        row["file_sha256"] = sha256(candidates_dir / fname)
        row["geometry_sha256"] = geometry_sha(written)
        row["atom_displacements_from_parent"] = nonzero_displacements(parent, written)
        row["atom_count"] = len(written)
        row["formula"] = written.get_chemical_formula()
        row["cell_A"] = written.cell.array.tolist()
        row["pbc"] = written.pbc.astype(bool).tolist()
        row["marked_pair"] = {
            "Ag_id": int(written.arrays["lammps_id"][ag]),
            "framework_id": int(written.arrays["lammps_id"][x]),
            "framework_element": written[x].symbol,
            "distance_A": float(written.get_distance(ag, x, mic=True)),
        }
        row["species_pair_minima_A"] = pair_minima(written)

    floors_by_formula = {}
    for item in inventory:
        fkey = json.dumps(formula_key(item["atoms_object"]))
        floors = floors_by_formula.setdefault(fkey, {})
        for pair, val in pair_minima(item["atoms_object"]).items():
            floors[pair] = min(floors.get(pair, math.inf), val["distance_A"])

    screen_rows = []
    for row in selected:
        atoms = row["written_atoms"]
        hits = []
        comparable = []
        for item in inventory:
            cmp = compare_geometry(atoms, item["atoms_object"])
            if cmp is None:
                continue
            comparable.append((cmp["all_atom_RMS_A"], item, cmp))
            if cmp["exact_duplicate"] or cmp["near_duplicate"]:
                hits.append({"path": item["path"], "frame_index_zero_based": item["frame_index_zero_based"], **cmp})
        nearest = None
        if comparable:
            _, item, cmp = min(comparable, key=lambda x: (x[0], x[1]["path"], x[1]["frame_index_zero_based"]))
            nearest = {"path": item["path"], "frame_index_zero_based": item["frame_index_zero_based"], **cmp}
        floors = floors_by_formula.get(json.dumps(formula_key(atoms)), {})
        contact_screen = []
        for pair, observed in row["species_pair_minima_A"].items():
            floor = floors.get(pair)
            if floor is None:
                continue
            margin = float(observed["distance_A"] - floor)
            contact_screen.append({
                "species_pair": pair,
                "candidate_minimum_A": float(observed["distance_A"]),
                "candidate_ids": observed["atom_ids"],
                "observed_inventory_floor_A": float(floor),
                "margin_above_floor_A": margin,
                "minimum_allowed_A": float(floor - PAIR_FLOOR_ALLOWANCE_A),
                "pass": bool(observed["distance_A"] >= floor - PAIR_FLOOR_ALLOWANCE_A),
            })
        row["screening"] = {
            "comparable_inventory_frames": len(comparable),
            "exact_or_near_hits": hits,
            "nearest_inventory_geometry": nearest,
            "contact_floor_rule": "same-composition observed minimum per element pair; candidate may not be more than 0.05 A below it",
            "species_contact_checks": contact_screen,
            "exact_near_pass": len(hits) == 0,
            "species_contact_pass": all(x["pass"] for x in contact_screen),
        }
        row["decision"] = "PASS_GEOMETRY_TRIAGE" if row["screening"]["exact_near_pass"] and row["screening"]["species_contact_pass"] else "FAIL_GEOMETRY_TRIAGE"
        screen_rows.append(row)

    cross_role = []
    for i, left in enumerate(selected):
        for right in selected[i + 1:]:
            if formula_key(left["written_atoms"]) != formula_key(right["written_atoms"]):
                continue
            cmp = compare_geometry(left["written_atoms"], right["written_atoms"])
            if cmp is None:
                continue
            cross_role.append({"candidate_a": left["candidate_id"], "candidate_b": right["candidate_id"], **cmp})
    pair_hits = [x for x in cross_role if x["exact_duplicate"] or x["near_duplicate"]]
    need(not pair_hits, f"role candidate sets contain exact/near overlap: {pair_hits[:3]}")
    failures = [x["candidate_id"] for x in screen_rows if x["decision"] != "PASS_GEOMETRY_TRIAGE"]
    need(not failures, f"geometry/contact screen failed for: {failures}")

    manifest_paths = [
        repo / V16_CANDIDATE_DIR / "candidate_manifest.json",
        repo / V15_WIDE_DIR / "candidate_manifest.json",
        repo / V17_B_DIR / "input_manifest.json",
        repo / PROJECT / "pbe_interface_v17_main_targeted_acquisition/input_manifest.json",
        repo / PROJECT / "pbe_interface_v14_parallel_acquisition/input_manifest.json",
    ]
    source_manifest_records = [
        {"path": str(p.relative_to(repo)), "sha256": sha256(p), "role": "geometry/source manifest"}
        for p in manifest_paths
    ]
    summary = {
        "scope": "geometry-only V17 role proposal; no DFT output, unseen labels, model inference/training, MD or TTM read/run",
        "reader_policy": "extxyz parser converts only species, positions, lammps_id, lammps_type, central_pair, Lattice and PBC; energy/force fields are ignored",
        "comparison_inventory": {
            "historical_model_split_files": len(split_files),
            "historical_model_split_frames": sum(x["frames"] for x in split_summary),
            "published_pbe_input_files": len(input_files),
            "published_pbe_input_frames": sum(x["frames"] for x in source_file_records if x["role"] == "published_dft_input_geometry_only"),
            "total_geometry_frames": len(inventory),
            "V16_train_rows": len(extxyz_geometry_frames(repo / PROJECT / "mace_periodic_v16_interface_energy/data/train.extxyz")),
            "V17_planned_training_inputs": len(v17_input_checks),
            "V17_planned_training_input_records": v17_input_checks,
        },
        "duplicate_policy": {
            "exact": f"same composition/order or persistent IDs, same cell/PBC, maximum aligned minimum-image displacement <= {EXACT_MAX_A} A",
            "near": f"same persistent IDs where present, minimum-image align; Ag registry RMS <= {NEAR_AG_RMS_A} A AND framework RMS <= {NEAR_FRAMEWORK_RMS_A} A",
            "roles_pairwise_checked": True,
        },
        "species_contact_policy": f"for each element pair, candidate minimum distance must be >= the smallest observed same-composition split/input minimum minus {PAIR_FLOOR_ALLOWANCE_A} A; empirical geometry screen, not a universal bond cutoff",
        "candidate_count": len(selected),
        "development_validation_count": sum(x["role"] == "development_validation" for x in selected),
        "withheld_test_count": sum(x["role"] == "withheld_test" for x in selected),
        "pairwise_role_overlap_hits": pair_hits,
        "candidate_results": [],
        "limitations": [
            "All proposals inherit the existing Ti3SiC2-based small-cluster parent motifs; geometry separation does not establish independent morphology or thermal coverage.",
            "Three development-validation geometries are generated from existing parent inputs with interface-specific transverse Ag translations and a 0.15 A local framework displacement.",
            "Three withheld-test geometries reuse unlabelled wide-span proposals; their parent lineage remains correlated with existing train motifs.",
            "Three additional withheld-test geometries combine a lateral Ag translation with a 0.15 A local marked-framework displacement; they are controlled design probes, not relaxed structures.",
            "A must accept/freeze the role manifest before V17 training; no labels exist, so DFT labels need a separate owner assignment."
        ],
        "source_files": source_file_records,
        "source_manifests": source_manifest_records,
        "historical_split_files": split_summary,
    }
    role_rows = []
    for row in selected:
        s = row["screening"]
        role_rows.append({
            "candidate_id": row["candidate_id"],
            "interface": row["interface"],
            "proposed_role": row["role"],
            "role_status": "frozen_by_window_b_before_v17_training; A review required before use",
            "file": row["path"],
            "file_sha256": row["file_sha256"],
            "geometry_sha256": row["geometry_sha256"],
            "atom_count": row["atom_count"],
            "formula": row["formula"],
            "cell_A": row["cell_A"],
            "pbc": row["pbc"],
            "marked_pair": row["marked_pair"],
            "direct_source_path": row["direct_source_path"],
            "direct_source_sha256": row["direct_source_sha256"],
            "parent_input_path": row["parent_input_path"],
            "parent_input_sha256": row["parent_input_sha256"],
            "lineage": row["lineage"],
            "atom_displacements_from_parent": row["atom_displacements_from_parent"],
            "species_pair_minima_A": row["species_pair_minima_A"],
            "screening": s,
            "decision": row["decision"],
            "label_state": "geometry only; no energy or force fields",
            "allowed_use": "geometry proposal only; do not train or tune on withheld-test labels",
        })
        summary["candidate_results"].append({
            "candidate_id": row["candidate_id"],
            "interface": row["interface"],
            "role": row["role"],
            "file_sha256": row["file_sha256"],
            "geometry_sha256": row["geometry_sha256"],
            "decision": row["decision"],
            "nearest_inventory_geometry": s["nearest_inventory_geometry"],
            "exact_or_near_hit_count": len(s["exact_or_near_hits"]),
            "contact_check_count": len(s["species_contact_checks"]),
            "species_contact_pass": s["species_contact_pass"],
        })
    inventory_json = {
        "comparison_scope": summary["comparison_inventory"],
        "geometry_only_sources": source_file_records,
        "source_manifests": source_manifest_records,
        "note": "Input and model-split coordinates only; calculated DFT folders and output labels were not traversed.",
    }
    (output / "source_inventory.json").write_text(json.dumps(inventory_json, indent=2, allow_nan=False) + "\n")
    (output / "role_manifest.json").write_text(json.dumps({
        "task": "v17_split_geometry_design",
        "owner_instance": "c5035b48-f43b-4da0-b8e4-e2862817f86a",
        "roles_frozen_before_v17_training": True,
        "role_policy": {
            "development_validation": "three structures, one per interface; usable only for development validation after A review",
            "withheld_test": "six structures, two per interface; keep labels sealed from checkpoint/parameter selection",
            "geometry_sets_pairwise_separated": True,
            "parent_motif_independent": False,
        },
        "candidates": role_rows,
        "limitations": summary["limitations"],
    }, indent=2, allow_nan=False) + "\n")
    (output / "screening_manifest.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")

    candidate_lines = []
    for row in selected:
        s = row["screening"]
        d = row["marked_pair"]["distance_A"]
        candidate_lines.append(
            f'| {row["candidate_id"]} | {row["role"]} | {row["atom_count"]} | '
            f'{row["marked_pair"]["Ag_id"]}–{row["marked_pair"]["framework_id"]} ({d:.4f}) | '
            f'{len(s["exact_or_near_hits"])} | {len(row["species_pair_minima_A"])} | {row["file_sha256"][:12]}… |'
        )
    report = [
        "# V17 validation/test geometry design (label-free)",
        "",
        "## Current decision",
        "",
        "Window B prepared three development-validation proposals and six withheld-test geometry proposals, all before V17 training. Every structure is geometry-only and preserves atom count, atom order, persistent IDs, cell and PBC. The role manifest is proposed as frozen; A must review and accept it before training. Withheld-test labels must not be used for checkpoint or parameter selection.",
        "",
        "## Candidate screen",
        "",
        "| Candidate | Proposed role | Atoms | Marked pair (ID; distance Å) | Exact/near inventory hits | Species-pair checks | File SHA256 |",
        "|---|---|---:|---|---:|---:|---|",
        *candidate_lines,
        "",
        f'The comparison inventory contains {summary["comparison_inventory"]["historical_model_split_frames"]} frames from {summary["comparison_inventory"]["historical_model_split_files"]} prior model-split files and {summary["comparison_inventory"]["published_pbe_input_frames"]} geometry frames from {summary["comparison_inventory"]["published_pbe_input_files"]} published DFT input files. This includes all {summary["comparison_inventory"]["V16_train_rows"]} V16 training rows and the eight hash-verified planned V17 training inputs. Historical probe inputs are included via their input paths only.',
        "",
        "Exact duplicates use persistent IDs when available, equal PBC/cell and maximum minimum-image displacement ≤1e-5 Å. Near duplicates use Ag registry RMS ≤0.15 Å and non-Ag framework RMS ≤0.15 Å after a common framework translation is removed. Candidate roles were checked pairwise using the same rules. Contact checks compare each species-pair minimum with the smallest same-composition minimum in existing split/input geometries, allowing at most 0.05 Å compression; this is an empirical geometry guard, not a universal bond criterion.",
        "",
        "## Role and lineage limits",
        "",
        "The three development-validation candidates are newly generated from the original parent inputs with interface-specific lateral Ag shifts plus a +0.15 Å displacement of the marked framework atom away from its marked Ag partner. For each interface, one withheld-test geometry combines a +0.25/+0.25 Å lateral Ag translation with the same local framework mode. The second withheld-test geometry reuses one unlabelled distance-axis proposal from the earlier wide-span candidate set. Per-ID displacement vectors, parent input hashes, original parent frame hashes, geometry hashes, contact minima and nearest inventory hits are in role_manifest.json.",
        "",
        "These candidates are isolated from committed train/validation/test geometries and DFT input geometries under the stated thresholds, but they all inherit existing small-cluster Ti3SiC2 parent motifs. They do not establish morphology independence or thermal coverage. The six combined transverse/framework geometries are unrelaxed controlled designs. Geometry screening cannot establish DFT convergence or physical stability.",
        "",
        "No DFT output, unseen reference label, MACE inference/training, MD or TTM was read or run. DFT labels require a separate owner assignment. Keep withheld-test labels sealed until model/checkpoint selection is frozen.",
        "",
        "## Reproduction",
        "",
        "Run the geometry-only generator with ASE and SciPy available:",
        "",
        "    cd /workspace/-",
        "    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /workspace/.venvs/gpaw-mpi/bin/python -B research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py --repo-root /workspace/- --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/reproduction_01",
        "",
        "The generator refuses an existing output directory, verifies source/manifest hashes and rechecks all nine geometry files, IDs, cell/PBC, duplicates and contact minima.",
    ]
    (output / "audit_report.md").write_text("\n".join(report) + "\n")
    (output / "README.md").write_text(
        "# V17 split geometry design\n\n"
        "This bundle proposes one development-validation geometry and two withheld-test geometries per interface. "
        "All files contain geometry and identity only; no energy/force labels or DFT results are included. "
        "The V16 proposals were rejected after exact matches to V16 blind inputs were found; see the parent report. "
        "See audit_report.md for screening rules, limitations and role boundaries. "
        "A must review the role_manifest.json before training.\n\n"
        "Candidates live in geometries/. source_inventory.json lists hashed geometry-only inputs and historical split files. "
        "screening_manifest.json contains machine-readable duplicate/contact results. "
        "SHA256SUMS.txt covers all bundle files except itself.\n"
    )
    shutil.copy2(Path(__file__), output / "generate_split_candidates.py")
    hash_lines = []
    for path in sorted(p for p in output.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
        hash_lines.append(f"{sha256(path)}  {path.relative_to(output).as_posix()}")
    (output / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    output = args.output_dir.resolve()
    allowed = (repo / OUTPUT_REL).resolve()
    need(output.is_relative_to(allowed) and output != allowed.parent, "output must be under the v17_split_geometry_design report directory")
    need(not output.exists(), f"refusing existing output directory: {output}")
    create_candidates(repo, output)
    print(json.dumps({
        "output_dir": str(output),
        "candidate_count": 9,
        "development_validation": 3,
        "withheld_test": 6,
        "result": "PASS_GEOMETRY_TRIAGE",
    }, indent=2))


if __name__ == "__main__":
    main()
