#!/usr/bin/env python3
"""Read-only independent audit of the three published V19 dev labels.

Only the three designated V19 development archives are parsed for labels.
Training, residual, V17 development, and V18 development files are parsed
through selected geometry columns (symbols, positions, atom IDs, cell, PBC).
No test/sealed path is constructed or read.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[6]
OUT = Path(__file__).resolve().parent
BASE = ROOT / "research/mace-v12-transfer"
P = BASE / "periodic_interface_v4"
ENTRY = P / "pbe_interface_v19_independent_dev"
MANIFEST_PATH = ENTRY / "input_manifest.json"
GEOM_AUDIT_PATH = ENTRY / "geometry_audit.json"
TRAIN = P / "mace_periodic_v19_residual_response/training_candidate/train.extxyz"
RESIDUAL_DIR = P / "pbe_interface_v19_residual_acquisition/inputs"
V17_DIR = BASE / "coordination/reports/window-b/v17_split_geometry_design/final_set/geometries"
V18_DEV = P / "mace_periodic_v18_clean_core/data/valid.extxyz"
V18_DEV_SHA256 = "985045431b3b69242ebbb9966393794a18a39fad37e029915e9d75f3b1e06f84"
SOURCE_CIF = P / "v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
METHOD = {"xc": "PBE", "basis": "PW", "cutoff_eV": 500, "kpts": [6, 6, 1], "smearing_eV": 0.1}
PBC = [True, True, False]
CELL_TOL_A = 1e-5
EXACT_MAX_A = 1e-5
NEAR_AG_RMS_A = 0.15
NEAR_FRAMEWORK_RMS_A = 0.15
REGISTRY_SHIFT_FRAC = (-0.27, 0.19)
LABELS = [
    "AgC_v19_independent_dev_01",
    "AgSi_v19_independent_dev_01",
    "AgTi_v19_independent_dev_01",
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def header_value(header: str, key: str) -> str | None:
    match = re.search(r"(?:^|\s)" + re.escape(key) + r'=(?:"([^"]*)"|(\S+))', header)
    return (match.group(1) if match.group(1) is not None else match.group(2)) if match else None


def parse_extxyz(path: Path, with_labels: bool = False) -> list[dict]:
    """Read selected geometry columns; label columns are parsed only on dev outputs."""
    frames: list[dict] = []
    with path.open("r", encoding="utf-8") as stream:
        frame_index = 0
        while True:
            first = stream.readline()
            if not first:
                break
            if not first.strip():
                continue
            n = int(first.strip())
            header = stream.readline().rstrip("\n")
            prop_text = header_value(header, "Properties")
            require(prop_text is not None, f"{path}: frame {frame_index} missing Properties")
            schema = prop_text.split(":")
            require(len(schema) % 3 == 0, f"{path}: malformed Properties")
            fields: dict[str, tuple[int, int, str]] = {}
            offset = 0
            for i in range(0, len(schema), 3):
                name, kind, width_text = schema[i : i + 3]
                width = int(width_text)
                fields[name] = (offset, width, kind)
                offset += width
            require("species" in fields and "pos" in fields, f"{path}: missing geometry columns")
            lat = header_value(header, "Lattice")
            require(lat is not None, f"{path}: missing Lattice")
            cell = np.asarray([float(x) for x in lat.split()], dtype=float).reshape(3, 3)
            pbc_text = header_value(header, "pbc")
            pbc = [x.upper() in {"T", "TRUE", "1"} for x in pbc_text.split()] if pbc_text else [True] * 3
            require(len(pbc) == 3, f"{path}: malformed pbc")
            syms: list[str] = []
            positions: list[list[float]] = []
            ids: list[int] = []
            types: list[int] = []
            marked: list[int] = []
            forces: list[list[float]] = []
            s0 = fields["species"][0]
            p0 = fields["pos"][0]
            id0 = fields.get("lammps_id", (-1, 0, ""))[0]
            type0 = fields.get("lammps_type", (-1, 0, ""))[0]
            marker0 = fields.get("central_pair", (-1, 0, ""))[0]
            force_field = next((name for name in ("PW_PBE_forces", "forces") if name in fields), None)
            force0 = fields[force_field][0] if force_field else -1
            for _ in range(n):
                row = stream.readline().split()
                require(len(row) >= offset, f"{path}: truncated atomic row")
                syms.append(row[s0])
                positions.append([float(x) for x in row[p0 : p0 + 3]])
                if id0 >= 0:
                    ids.append(int(row[id0]))
                if type0 >= 0:
                    types.append(int(row[type0]))
                if marker0 >= 0:
                    marked.append(int(row[marker0]))
                if with_labels and force0 >= 0:
                    forces.append([float(x) for x in row[force0 : force0 + 3]])
            atoms = Atoms(symbols=syms, positions=positions, cell=cell, pbc=pbc)
            if ids:
                atoms.new_array("lammps_id", np.asarray(ids, dtype=int))
            if types:
                atoms.new_array("lammps_type", np.asarray(types, dtype=int))
            if marked:
                atoms.new_array("central_pair", np.asarray(marked, dtype=int))
            if with_labels and forces:
                atoms.new_array("PW_PBE_forces", np.asarray(forces, dtype=float))
            energy = None
            if with_labels:
                raw_energy = header_value(header, "PW_PBE_energy_eV")
                energy = float(raw_energy) if raw_energy is not None else None
            frames.append({"atoms": atoms, "header": header, "fields": fields, "energy": energy})
            frame_index += 1
    return frames


def geometry_sha(atoms: Atoms) -> str:
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(atoms.positions, 10).tolist(),
        "cell_A": np.round(atoms.cell.array, 10).tolist(),
        "pbc": atoms.pbc.astype(bool).tolist(),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def formula_key(atoms: Atoms) -> tuple:
    return tuple(sorted(Counter(atoms.get_chemical_symbols()).items()))


def id_map(atoms: Atoms) -> dict[int, int] | None:
    values = atoms.arrays.get("lammps_id")
    if values is None or len(values) != len(atoms) or len(set(map(int, values))) != len(atoms):
        return None
    return {int(value): i for i, value in enumerate(values)}


def geometry_compare(a: Atoms, b: Atoms) -> dict | None:
    """Periodic species-aware screen; no rotation or cell scaling is applied."""
    if len(a) != len(b) or formula_key(a) != formula_key(b) or not np.array_equal(a.pbc, b.pbc):
        return None
    cell_delta = float(np.max(np.abs(a.cell.array - b.cell.array)))
    if cell_delta > CELL_TOL_A:
        return None
    sa = np.asarray(a.get_chemical_symbols())
    sb = np.asarray(b.get_chemical_symbols())
    amap, bmap = id_map(a), id_map(b)
    if amap is not None and bmap is not None and set(amap) == set(bmap):
        ids = sorted(amap)
        ia = np.asarray([amap[x] for x in ids], dtype=int)
        ib = np.asarray([bmap[x] for x in ids], dtype=int)
        method = "persistent_lammps_id"
    elif list(sa) == list(sb):
        ia = np.arange(len(a), dtype=int)
        ib = np.arange(len(b), dtype=int)
        method = "same_element_order"
    else:
        # Species-preserving Hungarian assignment, iterated with a framework
        # median translation to handle different IDs/order and periodic images.
        shift = np.zeros(3)
        for _ in range(5):
            pairs_a: list[int] = []
            pairs_b: list[int] = []
            for element in sorted(set(sa)):
                xa, xb = np.flatnonzero(sa == element), np.flatnonzero(sb == element)
                delta = (a.positions[xa, None, :] - b.positions[None, xb, :] - shift).reshape(-1, 3)
                _, distances = find_mic(delta, a.cell, pbc=a.pbc)
                rows, cols = linear_sum_assignment(np.asarray(distances).reshape(len(xa), len(xb)) ** 2)
                pairs_a.extend(xa[rows].tolist())
                pairs_b.extend(xb[cols].tolist())
            ia, ib = np.asarray(pairs_a), np.asarray(pairs_b)
            delta, _ = find_mic(a.positions[ia] - b.positions[ib], a.cell, pbc=a.pbc)
            framework = sa[ia] != "Ag"
            new_shift = np.median(delta[framework] if framework.any() else delta, axis=0)
            if float(np.linalg.norm(new_shift - shift)) < 1e-10:
                shift = new_shift
                break
            shift = new_shift
        method = "species_hungarian_periodic_mic"
    delta, _ = find_mic(a.positions[ia] - b.positions[ib], a.cell, pbc=a.pbc)
    framework = sa[ia] != "Ag"
    delta = delta - np.median(delta[framework] if framework.any() else delta, axis=0)
    delta, _ = find_mic(delta, a.cell, pbc=a.pbc)
    magnitudes = np.linalg.norm(delta, axis=1)
    ag = sa[ia] == "Ag"
    framework = ~ag
    ag_rms = float(np.sqrt(np.mean(magnitudes[ag] ** 2))) if ag.any() else 0.0
    framework_rms = float(np.sqrt(np.mean(magnitudes[framework] ** 2))) if framework.any() else 0.0
    return {
        "matching_method": method,
        "cell_max_delta_A": cell_delta,
        "Ag_RMS_A": ag_rms,
        "framework_RMS_A": framework_rms,
        "all_atom_RMS_A": float(np.sqrt(np.mean(magnitudes**2))),
        "max_atom_displacement_A": float(np.max(magnitudes)),
        "exact_duplicate": bool(float(np.max(magnitudes)) <= EXACT_MAX_A),
        "near_duplicate": bool(ag_rms <= NEAR_AG_RMS_A and framework_rms <= NEAR_FRAMEWORK_RMS_A),
    }


def rdf_fingerprint(atoms: Atoms) -> np.ndarray:
    distances = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    bins = np.linspace(0.0, 8.0, 81)
    values: list[float] = []
    elements = ("Ag", "C", "Si", "Ti")
    for e1 in elements:
        for e2 in elements:
            pair = [
                float(distances[i, j])
                for i in range(len(atoms))
                for j in range(i + 1, len(atoms))
                if tuple(sorted((symbols[i], symbols[j]))) == tuple(sorted((e1, e2)))
            ]
            hist = np.histogram(pair, bins=bins)[0].astype(float)
            if pair:
                hist /= len(pair)
            values.extend(hist.tolist())
    return np.asarray(values, dtype=float)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denom) if denom else 0.0


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_record(path: Path, expected: str | None = None) -> dict:
    actual = sha(path)
    require(expected is None or actual == expected, f"SHA256 mismatch: {relative(path)} expected {expected}, got {actual}")
    return {"path": relative(path), "sha256": actual, "expected_sha256": expected}


def audit() -> dict:
    manifest = json.loads(MANIFEST_PATH.read_text())
    geometry_audit = json.loads(GEOM_AUDIT_PATH.read_text())
    require(manifest["status"] == "READY_LAUNCH_ENABLED", "Unexpected input manifest state")
    require(manifest["dataset_role"] == "independent_development" and manifest["roles_frozen_before_DFT"] is True,
            "Role/freeze metadata mismatch")
    require(manifest["method"] == METHOD and manifest["pbc"] == PBC and manifest["poissonsolver"] == {"dipolelayer": "xy"},
            "Common method/PBC metadata mismatch")
    require(geometry_audit["status"] == "PASS_GEOMETRY_AND_ROLE_SCREEN_PENDING_EXTERNAL_B_AUDIT",
            "Unexpected A geometry audit status")

    evidence: dict[str, dict] = {}
    evidence["input_manifest"] = file_record(MANIFEST_PATH)
    evidence["geometry_audit"] = file_record(GEOM_AUDIT_PATH, manifest["geometry_audit_sha256"])
    evidence["source_CIF"] = file_record(ROOT / manifest["source_CIF"], manifest["source_CIF_sha256"])
    runner_path = ENTRY / "run_pw_reference.py"
    evidence["reference_runner"] = file_record(runner_path)
    runner_text = runner_path.read_text()
    require("world.size != 4" in runner_text and "get_potential_energy(force_consistent=False)" in runner_text,
            "Registered runner does not document four ranks/native extrapolated energy")
    require('calc.results["free_energy"]' in runner_text, "Registered runner does not record free energy separately")
    require(evidence["source_CIF"]["sha256"] == "6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3",
            "Source CIF hash differs from the registered COD 9009647 source")

    candidate_records = {row["label"]: row for row in manifest["records"]}
    audit_records = {row["label"]: row for row in geometry_audit["candidate_records"]}
    require(set(candidate_records) == set(LABELS) and set(audit_records) == set(LABELS), "Unexpected label inventory")
    require(len({row["source_CIF_sha256"] for row in manifest["records"]}) == 1,
            "Entries do not share the declared source CIF")

    candidates: dict[str, dict] = {}
    archive_results: dict[str, dict] = {}
    source_lineage: dict[str, dict] = {}
    for label in LABELS:
        row = candidate_records[label]
        geom_row = audit_records[label]
        require(row["role"] == "independent_development" and row["A_lineage_review"] == "APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY",
                f"{label}: role metadata mismatch")
        require(row["source_CIF_sha256"] == manifest["source_CIF_sha256"], f"{label}: source hash metadata mismatch")
        input_path = ROOT / row["input"]
        require(row["input"] == geom_row["input"], f"{label}: input path differs in geometry audit")
        input_rec = file_record(input_path, row["input_sha256"])
        input_frames = parse_extxyz(input_path)
        require(len(input_frames) == 1, f"{label}: input is not single-frame")
        inp = input_frames[0]
        input_atoms: Atoms = inp["atoms"]
        require(len(input_atoms) == row["atoms"] == 52, f"{label}: atom count mismatch")
        require(formula_key(input_atoms) == tuple(sorted(Counter({"Ag": 4, "C": 16, "Si": 8, "Ti": 24}).items())),
                f"{label}: composition mismatch")
        require(input_atoms.pbc.tolist() == PBC, f"{label}: input PBC mismatch")
        require(header_value(inp["header"], "dataset_role") == "independent_development", f"{label}: missing input role")
        require(header_value(inp["header"], "candidate_role") == "independent_development", f"{label}: missing candidate role")
        require(header_value(inp["header"], "role_frozen_before_labeling") == "T", f"{label}: role freeze metadata mismatch")
        require(header_value(inp["header"], "proposal_only_no_DFT_label") == "T", f"{label}: input is not label-free")
        require(header_value(inp["header"], "config_type") == label, f"{label}: config_type metadata mismatch")
        require(header_value(inp["header"], "source_path") == manifest["source_CIF"], f"{label}: source path metadata mismatch")
        require(header_value(inp["header"], "source_sha256") == manifest["source_CIF_sha256"], f"{label}: source hash metadata mismatch")
        require(header_value(inp["header"], "parent_geometry_path") == geom_row["parent_geometry"], f"{label}: parent path metadata mismatch")
        require(header_value(inp["header"], "parent_geometry_sha256") == geom_row["parent_geometry_sha256"], f"{label}: parent hash metadata mismatch")
        require(not any(name.lower() in {"forces", "pw_pbe_forces", "ref_forces", "energy", "ref_energy"} for name in inp["fields"]),
                f"{label}: input unexpectedly contains labels")
        geom_hash = geometry_sha(input_atoms)
        require(geom_hash == row["geometry_sha256"] == geom_row["geometry_sha256"], f"{label}: geometry hash mismatch")
        require(input_rec["sha256"] == geom_row["input_sha256"], f"{label}: input hash differs in geometry audit")
        require(np.allclose(input_atoms.cell.array, np.asarray(geom_row["cell_A"]), atol=1e-10, rtol=0.0),
                f"{label}: cell differs from geometry audit")

        parent_path = ROOT / geom_row["parent_geometry"]
        parent_rec = file_record(parent_path, geom_row["parent_geometry_sha256"])
        parent_frames = parse_extxyz(parent_path)
        require(len(parent_frames) == 1, f"{label}: parent geometry is not single-frame")
        parent = parent_frames[0]["atoms"]
        imap, pmap = id_map(input_atoms), id_map(parent)
        require(imap is not None and pmap is not None and set(imap) == set(pmap), f"{label}: parent/input atom IDs differ")
        ids = sorted(imap)
        in_ix = np.asarray([imap[x] for x in ids])
        par_ix = np.asarray([pmap[x] for x in ids])
        require(np.array_equal(np.asarray(input_atoms.get_chemical_symbols())[in_ix], np.asarray(parent.get_chemical_symbols())[par_ix]),
                f"{label}: parent/input ID species mismatch")
        displacement, _ = find_mic(input_atoms.positions[in_ix] - parent.positions[par_ix], input_atoms.cell, pbc=input_atoms.pbc)
        sym = np.asarray(input_atoms.get_chemical_symbols())[in_ix]
        framework_disp = np.linalg.norm(displacement[sym != "Ag"], axis=1)
        shift = REGISTRY_SHIFT_FRAC[0] * input_atoms.cell[0] + REGISTRY_SHIFT_FRAC[1] * input_atoms.cell[1]
        ag_disp, _ = find_mic(displacement[sym == "Ag"] - shift, input_atoms.cell, pbc=input_atoms.pbc)
        require(float(np.max(framework_disp)) <= 1e-8, f"{label}: framework changed relative to registered parent")
        require(float(np.max(np.linalg.norm(ag_disp, axis=1))) <= 1e-8, f"{label}: Ag registry shift differs from recorded transformation")
        source_lineage[label] = {
            "source_CIF": evidence["source_CIF"],
            "parent_geometry": parent_rec,
            "parent_matches_registered_CIF_lineage": geom_row["source_CIF_sha256"] == manifest["source_CIF_sha256"],
            "framework_max_delta_from_parent_A": float(np.max(framework_disp)),
            "Ag_registry_shift_fractional": list(REGISTRY_SHIFT_FRAC),
            "Ag_max_delta_from_registered_shift_A": float(np.max(np.linalg.norm(ag_disp, axis=1))),
        }

        archive_dir = ENTRY / "calculations" / label
        checksum_path = archive_dir / "SHA256SUMS.txt"
        sums: dict[str, str] = {}
        for line in checksum_path.read_text().splitlines():
            if not line.strip():
                continue
            digest, filename = line.split(None, 1)
            filename = filename.strip().lstrip("*")
            require(Path(filename).name == filename, f"{label}: unsafe archive path in checksum list")
            require(filename not in sums, f"{label}: duplicate archive checksum entry")
            sums[filename] = digest
        output_name = f"{label}_PW_PBE.extxyz"
        required_archive = {output_name, "gpaw.log", "progress.json", "summary.json", "verification.json"}
        require(set(sums) == required_archive, f"{label}: archive checksum inventory mismatch")
        archive_hashes = {}
        for filename, expected in sums.items():
            archive_hashes[filename] = file_record(archive_dir / filename, expected)["sha256"]
        archive_results[label] = {"SHA256SUMS_path": relative(checksum_path), "checksums": archive_hashes}

        summary_path = archive_dir / "summary.json"
        verification_path = archive_dir / "verification.json"
        progress_path = archive_dir / "progress.json"
        log_path = archive_dir / "gpaw.log"
        summary = json.loads(summary_path.read_text())
        verification = json.loads(verification_path.read_text())
        progress = json.loads(progress_path.read_text())
        log_text = log_path.read_text(errors="replace")
        require(verification["status"] == "PASS", f"{label}: archive verifier did not pass")
        require(verification["withheld_test_paths_opened"] is False, f"{label}: archive metadata reports test access")
        require(verification["role"] == "independent_development", f"{label}: archive role mismatch")
        require(verification["input_sha256"] == input_rec["sha256"], f"{label}: verification input hash mismatch")
        require(verification["geometry_sha256"] == geom_hash, f"{label}: verification geometry hash mismatch")
        require(verification["fixed_geometry_exactly_preserved"] is True, f"{label}: verifier does not confirm fixed geometry")
        require(verification["archive_sha256"] == {k: v for k, v in archive_hashes.items() if k != "verification.json"},
                f"{label}: verification archive hash map mismatch")
        require(progress["status"] == "complete" and progress["scf_converged"] is True, f"{label}: progress not complete")
        require(progress["input_sha256"] == input_rec["sha256"] and progress["mpi_ranks"] == 4,
                f"{label}: progress hash/MPI metadata mismatch")
        require(summary["label"] == label and summary["dataset_role"] == "independent_development", f"{label}: summary role mismatch")
        require(summary["input_sha256"] == input_rec["sha256"] and summary["geometry_sha256"] == geom_hash,
                f"{label}: summary lineage hash mismatch")
        require(summary["method"] == METHOD and summary["pbc"] == PBC and summary["poissonsolver"] == manifest["poissonsolver"],
                f"{label}: summary method/PBC mismatch")
        require(summary["energy_convention"] == "GPAW native extrapolated energy, force_consistent=False; free energy separately recorded",
                f"{label}: energy convention mismatch")
        require(summary["frozen_geometry"] is True and summary["scf_converged"] is True and summary["mpi_ranks"] == 4,
                f"{label}: convergence/fixed geometry/MPI metadata mismatch")
        require(summary["scf_iterations"] == progress["iteration"], f"{label}: SCF iteration metadata mismatch")
        require("mode=PW(ecut=500.0)" in log_text and "size=array([6, 6, 1])" in log_text,
                f"{label}: GPAW log method evidence mismatch")
        require("width=0.1" in log_text and "dipolelayer='xy'" in log_text and "xc=XC(name='PBE')" in log_text,
                f"{label}: GPAW log smearing/dipole/XC evidence mismatch")
        converged_match = re.search(r"Converged in (\d+) steps", log_text)
        require(converged_match is not None and int(converged_match.group(1)) == summary["scf_iterations"],
                f"{label}: GPAW log convergence count mismatch")
        extrapolated = re.findall(r"^\s*Extrapolated:\s*([-+0-9.eE]+)", log_text, re.MULTILINE)
        free = re.findall(r"^\s*Free energy:\s*([-+0-9.eE]+)", log_text, re.MULTILINE)
        require(extrapolated and free, f"{label}: GPAW log lacks separated extrapolated/free energies")
        require(abs(float(extrapolated[-1]) - summary["energy_eV_cell"]) <= 5e-6,
                f"{label}: extrapolated energy differs from summary")
        require(abs(float(free[-1]) - summary["free_energy_eV_cell"]) <= 5e-6,
                f"{label}: free energy differs from summary")

        output_path = archive_dir / output_name
        output_frames = parse_extxyz(output_path, with_labels=True)
        require(len(output_frames) == 1, f"{label}: output is not single-frame")
        out = output_frames[0]
        out_atoms: Atoms = out["atoms"]
        require(out["energy"] is not None and math.isfinite(out["energy"]), f"{label}: missing/non-finite output energy")
        require(out_atoms.arrays.get("PW_PBE_forces") is not None, f"{label}: output force property missing")
        forces = out_atoms.arrays["PW_PBE_forces"]
        require(forces.shape == (len(input_atoms), 3) and np.isfinite(forces).all(), f"{label}: invalid/non-finite forces")
        require(np.array_equal(input_atoms.get_chemical_symbols(), out_atoms.get_chemical_symbols()), f"{label}: atom order/species changed")
        require(np.array_equal(input_atoms.arrays["lammps_id"], out_atoms.arrays["lammps_id"]), f"{label}: atom IDs/order changed")
        for array_name in ("lammps_type", "central_pair"):
            require(np.array_equal(input_atoms.arrays[array_name], out_atoms.arrays[array_name]), f"{label}: {array_name} changed")
        pos_delta = float(np.max(np.abs(input_atoms.positions - out_atoms.positions)))
        cell_delta = float(np.max(np.abs(input_atoms.cell.array - out_atoms.cell.array)))
        require(pos_delta <= 1e-10 and cell_delta <= 1e-10 and np.array_equal(input_atoms.pbc, out_atoms.pbc),
                f"{label}: fixed positions/cell/PBC differ")
        require(abs(out["energy"] - summary["energy_eV_cell"]) <= 1e-10, f"{label}: output energy differs from summary")
        fmax = float(np.max(np.linalg.norm(forces, axis=1)))
        rms_force = float(np.sqrt(np.mean(np.sum(forces**2, axis=1))))
        require(abs(fmax - summary["fmax_eV_A"]) <= 1e-8, f"{label}: force maximum differs from summary")
        require(abs(rms_force - summary["rms_force_eV_A"]) <= 1e-8, f"{label}: force RMS differs from summary")
        contact_indices = np.flatnonzero(input_atoms.arrays["central_pair"] == 1)
        require(len(contact_indices) == 2, f"{label}: marked contact pair count mismatch")
        target = {"AgC": "C", "AgSi": "Si", "AgTi": "Ti"}[row["interface"]]
        target_meta = geom_row["target_contact_pair"]
        marked_symbols = input_atoms.get_chemical_symbols()
        require(sorted(marked_symbols[i] for i in contact_indices) == sorted(["Ag", target]),
                f"{label}: marked contact species mismatch")
        marked_ids = {marked_symbols[i]: int(input_atoms.arrays["lammps_id"][i]) for i in contact_indices}
        require(marked_ids == {
            "Ag": int(target_meta["ag_lammps_id"]),
            target: int(target_meta["framework_lammps_id"]),
        }, f"{label}: marked contact IDs differ from geometry audit")
        contact_delta, contact_dist = find_mic(
            input_atoms.positions[contact_indices[0]] - input_atoms.positions[contact_indices[1]],
            input_atoms.cell,
            pbc=input_atoms.pbc,
        )
        del contact_delta
        contact_distance = float(contact_dist)
        require(abs(contact_distance - target_meta["distance_A"]) <= 1e-8, f"{label}: contact distance metadata mismatch")
        candidates[label] = {
            "interface": row["interface"],
            "role": row["role"],
            "input": input_rec,
            "output": file_record(output_path, archive_hashes[output_name]),
            "geometry_sha256": geom_hash,
            "atoms": len(input_atoms),
            "formula": "C16Ag4Si8Ti24",
            "pbc": input_atoms.pbc.tolist(),
            "input_label_fields_absent": True,
            "fixed_geometry": {
                "species_and_order_equal": True,
                "lammps_ids_and_order_equal": True,
                "positions_max_abs_delta_A": pos_delta,
                "cell_max_abs_delta_A": cell_delta,
                "pbc_equal": True,
            },
            "archive": {
                "verification_status": verification["status"],
                "scf_converged": summary["scf_converged"],
                "scf_iterations": summary["scf_iterations"],
                "method": summary["method"],
                "poissonsolver": summary["poissonsolver"],
                "energy_convention": summary["energy_convention"],
                "energy_eV_cell_extrapolated": summary["energy_eV_cell"],
                "free_energy_eV_cell": summary["free_energy_eV_cell"],
                "output_energy_matches_summary": True,
                "log_extrapolated_and_free_energy_match_summary": True,
                "forces_finite": True,
                "force_max_eV_A_recomputed": fmax,
                "force_rms_eV_A_recomputed": rms_force,
                "mpi_ranks": summary["mpi_ranks"],
                "archive_checksums": archive_results[label]["checksums"],
            },
            "marked_contact": {
                "species": f"Ag-{target}",
                "distance_A": contact_distance,
                "metadata_match": True,
            },
        }

    # Hash- and geometry-only screening inventory. Do not parse label columns.
    geom_reference = geometry_audit["reference_geometry_sets"]
    references: list[dict] = []
    train_rec = geom_reference["train36"]
    train_path = ROOT / train_rec["path"]
    train_file = file_record(train_path, train_rec["sha256"])
    train_frames = parse_extxyz(train_path)
    require(len(train_frames) == train_rec["frames"] == 36, "train36 frame count mismatch")
    evidence["train36"] = train_file
    for i, frame in enumerate(train_frames):
        references.append({"group": "train36", "id": f"train36_frame_{i:02d}", "atoms": frame["atoms"]})

    for record in geom_reference["six_v19_residual_inputs"]:
        path = ROOT / record["path"]
        rec = file_record(path, record["sha256"])
        evidence[f"residual::{record['label']}"] = rec
        frames = parse_extxyz(path)
        require(len(frames) == 1, f"{record['label']}: residual geometry not single-frame")
        references.append({"group": "six_v19_residual_inputs", "id": record["label"], "atoms": frames[0]["atoms"]})

    for record in geom_reference["three_public_v17_development_parents"]:
        path = ROOT / record["path"]
        rec = file_record(path, record["sha256"])
        evidence[f"v17_dev_parent::{record['label']}"] = rec
        frames = parse_extxyz(path)
        require(len(frames) == 1, f"{record['label']}: V17 dev geometry not single-frame")
        references.append({"group": "three_public_v17_development_parents", "id": record["label"], "atoms": frames[0]["atoms"]})

    v18_file = file_record(V18_DEV, V18_DEV_SHA256)
    evidence["v18_dev_valid_only"] = v18_file
    v18_frames = parse_extxyz(V18_DEV)
    require(len(v18_frames) == 3, "V18 development valid.extxyz is not the expected three-frame development set")
    for i, frame in enumerate(v18_frames):
        references.append({"group": "v18_development_valid_only", "id": f"v18_dev_frame_{i:02d}", "atoms": frame["atoms"]})

    comparisons: dict[str, dict] = {}
    rdf_vectors = {r["id"]: rdf_fingerprint(r["atoms"]) for r in references}
    for label, candidate in candidates.items():
        atoms = parse_extxyz(ROOT / candidate["input"]["path"])[0]["atoms"]
        groups: dict[str, dict] = {}
        for group in sorted({r["group"] for r in references}):
            selected = [r for r in references if r["group"] == group]
            comparable: list[dict] = []
            for ref in selected:
                result = geometry_compare(atoms, ref["atoms"])
                if result is not None:
                    comparable.append({"reference": ref["id"], **result})
            groups[group] = {
                "n_geometries_checked": len(selected),
                "n_comparable_same_composition_cell_pbc": len(comparable),
                "exact_duplicates": sum(x["exact_duplicate"] for x in comparable),
                "near_duplicates": sum(x["near_duplicate"] for x in comparable),
                "comparable_results": comparable,
                "max_species_pair_RDF_cosine_similarity_descriptive_only": max(
                    cosine(rdf_fingerprint(atoms), rdf_vectors[r["id"]]) for r in selected
                ),
            }
        comparisons[label] = groups

    candidate_pairwise = {}
    for left, right in itertools.combinations(LABELS, 2):
        a = parse_extxyz(ROOT / candidates[left]["input"]["path"])[0]["atoms"]
        b = parse_extxyz(ROOT / candidates[right]["input"]["path"])[0]["atoms"]
        candidate_pairwise[f"{left}__{right}"] = geometry_compare(a, b)

    all_geometry_hits = [
        {"candidate": label, "reference_group": group, "reference": item["reference"],
         "exact_duplicate": item["exact_duplicate"], "near_duplicate": item["near_duplicate"]}
        for label, group_map in comparisons.items()
        for group, result in group_map.items()
        for item in result["comparable_results"]
        if item["exact_duplicate"] or item["near_duplicate"]
    ]
    require(not all_geometry_hits, "Found an exact or near-duplicate cross-role geometry")
    pairwise_hits = [k for k, v in candidate_pairwise.items() if v and (v["exact_duplicate"] or v["near_duplicate"])]
    require(not pairwise_hits, "Found exact/near duplicates among the three V19 development entries")

    checks = [
        {"check": "manifest_and_registered_geometry_audit_hashes", "status": "PASS"},
        {"check": "three_input_hashes_roles_and_label_free_inputs", "status": "PASS"},
        {"check": "CIF_parent_lineage_and_registered_registry_transform", "status": "PASS"},
        {"check": "three_archive_SHA256SUMS_and_verification_records", "status": "PASS"},
        {"check": "converged_fixed_geometry_finite_label_method_and_MPI_evidence", "status": "PASS"},
        {"check": "train36_residual_V17_and_V18_dev_exact_near_geometry_screen", "status": "PASS"},
    ]
    audit_result = {
        "status": "PASS_WITH_SHARED_SOURCE_LIMITATION",
        "audit_kind": "read_only_independent_dataset_entry_audit",
        "audit_commit": "860bf6cb7aca2fa1fc57024a28c154d2f08c6228",
        "checks": checks,
        "decision_boundary": {
            "dataset_approved": False,
            "preflight_run": False,
            "training_run": False,
            "test_or_sealed_data_opened": False,
            "scope_note": "This report audits entries and split geometry; only A may record A_APPROVED_V19_DATASET after review.",
        },
        "geometry_thresholds": {
            "cell_tolerance_A": CELL_TOL_A,
            "exact_max_displacement_A": EXACT_MAX_A,
            "near_Ag_RMS_A": NEAR_AG_RMS_A,
            "near_framework_RMS_A": NEAR_FRAMEWORK_RMS_A,
            "alignment": "persistent IDs when the full ID sets match; otherwise same element order or species-preserving periodic MIC Hungarian assignment; framework-median translation; no rotation or cell scaling",
            "limitations": "Only equal-composition, equal-PBC, equal-cell (within tolerance) pairs are exact/near comparable. This local screen does not prove separation under all crystallographic symmetries or across different cells.",
        },
        "candidate_entries": candidates,
        "source_lineage": {
            "source_CIF": evidence["source_CIF"],
            "source_citation_claim": "COD 9009647 / AMCSD0013252; Kisi et al. (1998), DOI 10.1016/S0022-3697(98)00226-1",
            "remote_source_limit": manifest["source_evidence_limit"],
            "one_shared_source_group": True,
            "per_label": source_lineage,
            "interpretation": "The three entries have separately registered contact geometries and do not match the screened train/residual/V17/V18 development structures. They all inherit one uploaded Ti3SiC2 CIF parent, so the three rows are not three independent crystallographic source families. This supports split separation from the audited structures but does not establish source-diverse generalization or blind test performance.",
        },
        "geometry_screen": {
            "scope": ["train36 geometry columns only", "six V19 residual geometry inputs only", "three public V17 development parent geometries", "three V18 development frames from valid.extxyz only"],
            "training_or_residual_label_values_parsed": False,
            "test_or_sealed_paths_constructed_or_opened": False,
            "comparison_summary": comparisons,
            "candidate_pairwise": candidate_pairwise,
            "exact_or_near_hits": all_geometry_hits,
        },
        "evidence_files": evidence,
        "limitations": [
            "All three development entries share the same uploaded COD 9009647 Ti3SiC2 source CIF; treat them as one source family for independence claims.",
            "The COD source bytes/checksum were not externally verified; the repository records the earlier HTTP 403 and verifies only the uploaded CIF bytes and internal rebuild lineage.",
            "The fixed PBE/PW 500 eV, 6x6x1, sigma 0.10 eV, dipole-corrected slab labels are consistent; slab k-point, smearing, and vacuum convergence remain unresolved.",
            "A geometry screen with no exact/near hits is not a proof against all symmetry-equivalent structures and does not establish model accuracy or generalization.",
        ],
    }

    readme = [
        "# V19 independent dataset-entry audit",
        "",
        "## Result",
        "",
        "**PASS with a shared-source limitation.** All three fixed-geometry development entries pass the independent input, archive, method, label, and geometry checks listed below. This report does not approve a dataset; A retains that decision.",
        "",
        "- Input and archive SHA256 inventories match their pinned manifests and archive checksums.",
        "- GPAW logs show SCF convergence (AgC 59, AgSi 57, AgTi 76 iterations); output coordinates, cell, PBC, atom order, IDs, markers, energies, and finite forces match their registered inputs and summaries.",
        "- Native extrapolated energies match the output and log; free energies are separately recorded and cross-checked.",
        "- Exact/near screen found no hits against train36, six V19 residual input geometries, three public V17 development parents, or the three V18 development frames in `valid.extxyz`.",
        "",
        "## Independence boundary and caveats",
        "",
        "The three entries use distinct Ag contact configurations and are separated from the audited train/residual and V17/V18 development geometries under the stated periodic comparator. They all derive from the same uploaded Ti3SiC2 COD 9009647 CIF, so they represent one crystallographic source family, not three independent source families. This supports role separation from the checked geometries but does not establish source-diverse generalization or blind test performance. The repository cites COD 9009647 / AMCSD0013252 and Kisi et al. (1998), DOI 10.1016/S0022-3697(98)00226-1; the upstream bytes remain unverified after the recorded HTTP 403.",
        "",
        "The fixed PBE/PW 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, xy dipole-corrected labels are method-consistent. Slab k-point, smearing, and vacuum convergence remain unresolved; finite SCF labels do not establish production accuracy. The near-duplicate test only compares equal-composition structures with equal PBC and cells within 1e-5 A, uses 1e-5 A maximum displacement for exact geometry and Ag/framework RMS <=0.15 A for near geometry, and applies no rotation or cell scaling. It is a local screen, not a complete crystallographic symmetry proof.",
        "",
        "## Scope and decision state",
        "",
        "Only train36 geometry columns, the six residual input geometry files, public V17 development geometry files, and V18 `data/valid.extxyz` geometry columns were screened. The three designated V19 development archives were checked for labels. No test/sealed path was constructed or opened. No dataset was changed or approved, and no preflight, DFT/MPI, MACE inference, or training was run.",
        "",
        "`audit.json` contains per-file hashes, per-label evidence, geometry comparisons, and source lineage. `audit.py` reproduces the read-only checks; run it from the repository root with Python plus NumPy, ASE, and SciPy installed.",
        "",
    ]
    (OUT / "audit.json").write_text(json.dumps(audit_result, indent=2, sort_keys=True) + "\n")
    (OUT / "README.md").write_text("\n".join(readme))
    sums = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            sums.append(f"{sha(path)}  {path.name}\n")
    (OUT / "SHA256SUMS.txt").write_text("".join(sums))
    return audit_result


if __name__ == "__main__":
    result = audit()
    print(json.dumps({"status": result["status"], "checks": result["checks"]}, indent=2))
