#!/usr/bin/env python3
"""Read-only, geometry-only audit of the uploaded Ti3SiC2 CIF proposals."""
from __future__ import annotations

import collections
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read
from ase.spacegroup import crystal
from scipy.optimize import linear_sum_assignment


OUT = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[6]
BASE = REPO / "research/mace-v12-transfer/periodic_interface_v4"
SRC = BASE / "v19_cif_parent_structures"
OLD = BASE / "v19_new_parent_prototypes"
TRAIN = BASE / "mace_periodic_v19_residual_response/training_candidate/train.extxyz"
TRAIN_MANIFEST = BASE / "mace_periodic_v19_residual_response/training_candidate/manifest.json"
DEV_DIR = (REPO / "research/mace-v12-transfer/coordination/reports/window-b/"
           "v17_split_geometry_design/final_set/geometries")
DEV_NAMES = [
    "AgC_v17_dev_validation_01.extxyz",
    "AgSi_v17_dev_validation_01.extxyz",
    "AgTi_v17_dev_validation_01.extxyz",
]
INTERFACES = {
    "AgC": ("C", "AgC_cod9009647_slab_v19_dev_proposal.extxyz",
            "AgC_new_hex_slab_parent_v19_proposal.extxyz"),
    "AgSi": ("Si", "AgSi_cod9009647_slab_v19_dev_proposal.extxyz",
             "AgSi_new_hex_slab_parent_v19_proposal.extxyz"),
    "AgTi": ("Ti", "AgTi_cod9009647_slab_v19_dev_proposal.extxyz",
             "AgTi_new_hex_slab_parent_v19_proposal.extxyz"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def formula_counts(symbols: list[str] | np.ndarray) -> dict[str, int]:
    return dict(sorted(collections.Counter(map(str, symbols)).items()))


def cif_scalar(text: str, tag: str) -> str | None:
    match = re.search(rf"^{re.escape(tag)}\s+(.+?)\s*$", text, re.MULTILINE)
    return match.group(1).strip("'\"") if match else None


def cif_atom_sites(text: str) -> list[dict]:
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        j = i + 1
        headers = []
        while j < len(lines) and lines[j].strip().startswith("_"):
            headers.append(lines[j].strip())
            j += 1
        required = ["_atom_site_label", "_atom_site_fract_x",
                    "_atom_site_fract_y", "_atom_site_fract_z"]
        if not all(x in headers for x in required):
            continue
        indices = [headers.index(x) for x in required]
        rows = []
        while j < len(lines):
            raw = lines[j].strip()
            if not raw or raw.startswith("#"):
                j += 1
                continue
            if raw == "loop_" or raw.startswith("_"):
                break
            fields = shlex.split(raw)
            if len(fields) < len(headers):
                raise ValueError("short CIF atom-site row")
            label = fields[indices[0]]
            symbol_match = re.match(r"([A-Z][a-z]?)", label)
            if not symbol_match:
                raise ValueError(f"cannot infer element from CIF label {label}")
            rows.append({
                "label": label,
                "species": symbol_match.group(1),
                "fractional": [float(fields[k]) for k in indices[1:]],
            })
            j += 1
        return rows
    raise ValueError("CIF asymmetric atom-site loop not found")


def read_geometry_frames(path: Path) -> list[dict]:
    """Read only extxyz cell, PBC, symbols, and xyz columns; skip label columns."""
    frames = []
    with path.open("r", encoding="utf-8") as f:
        while True:
            count_line = f.readline()
            if not count_line:
                break
            n = int(count_line.strip())
            header = f.readline().strip()
            match = re.search(r'Lattice="([^"]+)"', header)
            if not match:
                raise ValueError(f"missing Lattice in {path}")
            cell = np.asarray([float(x) for x in match.group(1).split()]).reshape(3, 3)
            pbc_match = re.search(r'\bpbc="([TF]) ([TF]) ([TF])"', header)
            pbc = [x == "T" for x in pbc_match.groups()] if pbc_match else [True, True, True]
            symbols, positions = [], []
            for _ in range(n):
                fields = f.readline().split()
                if len(fields) < 4:
                    raise ValueError(f"short atom row in {path}")
                symbols.append(fields[0])
                positions.append([float(x) for x in fields[1:4]])
            frames.append({
                "n": n,
                "symbols": symbols,
                "positions": np.asarray(positions, dtype=float),
                "cell": cell,
                "pbc": pbc,
                "composition": formula_counts(symbols),
            })
    return frames


def geometry_fingerprint(frame: dict, cutoff_A: float = 3.3) -> dict:
    """Species-pair/periodic-contact fingerprint; cutoff is descriptive, not a gate."""
    atoms = Atoms(symbols=frame["symbols"], positions=frame["positions"],
                  cell=frame["cell"], pbc=frame["pbc"])
    atoms.wrap()
    mic = atoms.get_all_distances(mic=True)
    direct = atoms.get_all_distances(mic=False)
    symbols = atoms.get_chemical_symbols()
    pair_counts: dict[str, int] = {}
    image_pair_counts: dict[str, int] = {}
    for i in range(len(atoms)):
        for j in range(i):
            key = "-".join(sorted((symbols[i], symbols[j])))
            if mic[i, j] < cutoff_A:
                pair_counts[key] = pair_counts.get(key, 0) + 1
                if direct[i, j] >= cutoff_A:
                    image_pair_counts[key] = image_pair_counts.get(key, 0) + 1
    return {
        "diagnostic_cutoff_A": cutoff_A,
        "short_pair_counts_by_species": dict(sorted(pair_counts.items())),
        "short_pairs_using_periodic_images_by_species": dict(sorted(image_pair_counts.items())),
        "short_pairs_using_periodic_images_total": int(sum(image_pair_counts.values())),
    }


def pair_minima(atoms) -> dict[str, float]:
    symbols = atoms.get_chemical_symbols()
    distances = atoms.get_all_distances(mic=True)
    minima: dict[str, float] = {}
    for i in range(len(atoms)):
        for j in range(i):
            key = "-".join(sorted((symbols[i], symbols[j])))
            minima[key] = min(minima.get(key, float("inf")), float(distances[i, j]))
    return dict(sorted(minima.items()))


def fit_bulk_origin_translation(source, assumed) -> dict:
    """Match full, symmetry-expanded chemical site sets with a common origin shift."""
    fs = source.get_scaled_positions(wrap=True)
    fa = assumed.get_scaled_positions(wrap=True)
    ss = np.asarray(source.get_chemical_symbols())
    sa = np.asarray(assumed.get_chemical_symbols())
    cell = source.cell.array
    best = (float("inf"), float("inf"), None)
    for i in range(len(source)):
        for j in range(len(assumed)):
            if ss[i] != sa[j]:
                continue
            shift = fs[i] - fa[j]
            sq_sum = 0.0
            max_dist = 0.0
            total = 0
            for element in sorted(set(ss)):
                ps = fs[ss == element]
                pa = (fa[sa == element] + shift) % 1.0
                frac_delta = ps[:, None, :] - pa[None, :, :]
                cart = frac_delta @ cell
                mic, _ = find_mic(cart.reshape(-1, 3), cell, pbc=True)
                dist = np.linalg.norm(mic, axis=1).reshape(len(ps), len(pa))
                rows, cols = linear_sum_assignment(dist)
                chosen = dist[rows, cols]
                sq_sum += float(np.dot(chosen, chosen))
                max_dist = max(max_dist, float(chosen.max()))
                total += len(chosen)
            rms = float(np.sqrt(sq_sum / total))
            if rms < best[0]:
                best = (rms, max_dist, shift % 1.0)
    return {
        "method": "species-wise Hungarian assignment after a single global fractional origin translation; no atom IDs used",
        "fractional_origin_translation": np.asarray(best[2]).round(8).tolist(),
        "rms_A_in_CIF_metric": best[0],
        "maximum_matched_deviation_A_in_CIF_metric": best[1],
        "source_cell_used_for_metric_A": cell.tolist(),
        "source_atoms": len(source),
        "assumed_atoms": len(assumed),
    }


def reproduce_builder() -> dict:
    """Rebuild an isolated temp copy at the same repo depth; never write A's inputs."""
    with tempfile.TemporaryDirectory(prefix="v19-uploaded-cif-audit-") as td:
        root = Path(td) / "workspace" / "-"
        copy = root / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures"
        copy.parent.mkdir(parents=True)
        shutil.copytree(SRC, copy)
        result = subprocess.run(
            ["/workspace/.venvs/gpaw-mpi/bin/python", str(copy / "build_from_cif.py")],
            cwd=root, capture_output=True, text=True, check=False,
        )
        if result.returncode:
            raise RuntimeError(f"isolated builder failed: {result.stderr[-2000:]}")
        rows = []
        for _, (_, candidate_name, _) in INTERFACES.items():
            original = SRC / "inputs" / candidate_name
            rebuilt = copy / "inputs" / candidate_name
            rows.append({
                "file": candidate_name,
                "byte_identical": original.read_bytes() == rebuilt.read_bytes(),
                "repository_sha256": sha256(original),
                "isolated_rebuild_sha256": sha256(rebuilt),
            })
        original_manifest = SRC / "manifest.json"
        rebuilt_manifest = copy / "manifest.json"
        return {
            "builder_exit_code": result.returncode,
            "stdout": result.stdout.strip(),
            "candidate_files": rows,
            "manifest_byte_identical": original_manifest.read_bytes() == rebuilt_manifest.read_bytes(),
            "manifest_repository_sha256": sha256(original_manifest),
            "manifest_rebuild_sha256": sha256(rebuilt_manifest),
        }


def proposal_metrics(interface: str, contact: str, filename: str, old_filename: str,
                     record: dict, source_cell: np.ndarray) -> dict:
    path = SRC / "inputs" / filename
    old_path = OLD / "inputs" / old_filename
    atoms = read(path, format="extxyz")
    old = read(old_path, format="extxyz")
    symbols = np.asarray(atoms.get_chemical_symbols())
    old_symbols = np.asarray(old.get_chemical_symbols())
    ag = np.flatnonzero(symbols == "Ag")
    substrate = np.flatnonzero(symbols != "Ag")
    z = atoms.positions[:, 2]
    substrate_top = float(z[substrate].max())
    substrate_bottom = float(z[substrate].min())
    top_elements = sorted(set(symbols[substrate][np.abs(z[substrate] - substrate_top) < 1e-7]))
    bottom_elements = sorted(set(symbols[substrate][np.abs(z[substrate] - substrate_bottom) < 1e-7]))
    total_top = sorted(set(symbols[np.abs(z - float(z.max())) < 1e-7]))
    ag_plane = float(np.mean(z[ag]))
    minima = pair_minima(atoms)
    old_minima = pair_minima(old)
    marked = np.flatnonzero(atoms.arrays["central_pair"])
    if len(marked) != 2:
        raise ValueError(f"expected exactly two marked atoms in {path}")
    marked_distance = float(atoms.get_distance(int(marked[0]), int(marked[1]), mic=True))
    marked_species = sorted((str(symbols[marked[0]]), str(symbols[marked[1]])))
    ag_xy = []
    for i in ag:
        for j in ag:
            if i >= j:
                continue
            d, _ = find_mic(atoms.positions[i] - atoms.positions[j], atoms.cell, pbc=[True, True, False])
            ag_xy.append(float(np.linalg.norm(d)))
    substrate_a = float(np.linalg.norm(atoms.cell.array[0]) / 2.0)
    ag_nn = min(ag_xy)
    energy_like = [k for k in list(atoms.info) + list(atoms.arrays)
                   if any(word in k.lower() for word in ("energy", "force"))]
    if atoms.calc is not None or energy_like:
        raise ValueError(f"labels or calculator attached to candidate {path}: {energy_like}")
    if len(atoms) != 52 or len(substrate) != 48 or len(ag) != 4:
        raise ValueError(f"wrong atom partition in {path}")
    if not np.all(atoms.pbc) or len(set(atoms.arrays["lammps_id"].tolist())) != 52:
        raise ValueError(f"PBC/ID check failed for {path}")
    if not atoms.info.get("new_ids_not_historical"):
        raise ValueError(f"new-ID flag missing for {path}")
    if atoms.info.get("dataset_role") != "CIF_DERIVED_DEVELOPMENT_PROPOSAL":
        raise ValueError(f"role mismatch in {path}")
    if top_elements != [contact] or marked_species != sorted(("Ag", contact)):
        raise ValueError(f"surface or marker species mismatch in {path}")
    if not np.isclose(ag_plane - substrate_top, float(record["gap_plane_A"]), atol=1e-8):
        raise ValueError(f"surface gap mismatch in {path}")
    manifest_delta = max(abs(minima[k] - float(record["species_pair_minima_A"][k])) for k in minima)
    if manifest_delta > 1e-7:
        raise ValueError(f"pair-minimum mismatch in {path}: {manifest_delta}")
    old_contact_key = "-".join(sorted(("Ag", contact)))
    return {
        "file": str(path.relative_to(REPO)),
        "sha256": sha256(path),
        "old_prototype_file": str(old_path.relative_to(REPO)),
        "old_prototype_sha256": sha256(old_path),
        "atoms": len(atoms),
        "composition": formula_counts(symbols),
        "substrate_atoms": len(substrate),
        "Ag_atoms": len(ag),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "substrate_a_A": substrate_a,
        "Ag_in_plane_minimum_Ag_A_A": ag_nn,
        "Ag_substrate_mismatch_percent_by_construction": 100.0 * (ag_nn - substrate_a) / substrate_a,
        "substrate_upper_termination": top_elements,
        "substrate_lower_termination": bottom_elements,
        "whole_structure_top_layer": total_top,
        "vertical_Ag_to_substrate_gap_A": ag_plane - substrate_top,
        "atom_z_extent_A": [float(z.min()), float(z.max())],
        "vacuum_from_z_extent_A": float(atoms.cell.lengths()[2] - (z.max() - z.min())),
        "element_pair_minima_A": minima,
        "marked_pair_ids": [int(atoms.arrays["lammps_id"][i]) for i in marked],
        "marked_pair_species": marked_species,
        "marked_pair_distance_A": marked_distance,
        "old_prototype_marked_contact_minimum_A": old_minima[old_contact_key],
        "delta_Ag_contact_vs_old_prototype_A": minima[old_contact_key] - old_minima[old_contact_key],
        "old_prototype_C_Ti_minimum_A": old_minima["C-Ti"],
        "delta_C_Ti_minimum_vs_old_A": minima["C-Ti"] - old_minima["C-Ti"],
        "lammps_id_range": [int(atoms.arrays["lammps_id"].min()), int(atoms.arrays["lammps_id"].max())],
        "geometry_only_header_properties": ["species", "pos", "spacegroup_kinds", "lammps_id", "central_pair"],
        "labels_or_calculator_present": False,
        "source_substrate_geometry_check": None,
        "file_manifest_pair_minimum_max_delta_A": manifest_delta,
    }


def verify_substrate_against_cif(contact: str, filename: str, metric: dict, bulk) -> dict:
    atoms = read(SRC / "inputs" / filename, format="extxyz")
    frac = bulk.get_scaled_positions()
    symbols = np.asarray(bulk.get_chemical_symbols())
    contact_fracs = np.unique(np.round(frac[symbols == contact, 2], 10))
    cut = (float(contact_fracs.max()) + 0.001) % 1.0
    unit = bulk.copy()
    shifted = unit.get_scaled_positions()
    shifted[:, 2] = (shifted[:, 2] - cut) % 1.0
    unit.set_scaled_positions(shifted)
    expected = unit.repeat((2, 2, 1))
    expected.positions[:, 2] += 6.0
    observed = atoms[:48]
    if expected.get_chemical_symbols() != observed.get_chemical_symbols():
        raise ValueError(f"substrate order/species mismatch for {filename}")
    delta = observed.positions - expected.positions
    mic, _ = find_mic(delta, atoms.cell, pbc=[True, True, False])
    norms = np.linalg.norm(mic, axis=1)
    if float(norms.max()) > 1e-8:
        raise ValueError(f"CIF substrate coordinate mismatch: max {norms.max()}")
    return {
        "operation": "ASE CIF expansion -> same target-plane cut -> repeat(2,2,1) -> +6 A z translation",
        "matched_atoms_by_order": len(expected),
        "max_ordered_MIC_deviation_A": float(norms.max()),
        "substrate_basis_cell_lengths_A": bulk.cell.lengths().tolist(),
        "in_plane_repeat": [2, 2, 1],
    }


def main() -> None:
    cif = SRC / "source/Ti3SiC2_COD_9009647.cif"
    cif_text = cif.read_text(encoding="utf-8")
    manifest = json.loads((SRC / "manifest.json").read_text())
    source = read(cif, format="cif")
    old_params = json.loads((OLD / "manifest.json").read_text())["parameters_assumed"]
    assumed = crystal(
        ["Ti", "Ti", "Si", "C"],
        basis=[(0, 0, 0), (1/3, 2/3, old_params["Ti_4f_z"]),
               (0, 0, 0.25), (1/3, 2/3, old_params["C_4f_z"])],
        spacegroup=old_params["spacegroup"], primitive_cell=False,
        cellpar=[old_params["a_A"], old_params["a_A"], old_params["c_A"], 90, 90, 120],
    )
    sites = cif_atom_sites(cif_text)
    source_kind_counts = collections.Counter(map(int, source.arrays["spacegroup_kinds"]))
    if len(sites) != 4 or [x["species"] for x in sites] != ["Ti", "Ti", "Si", "C"]:
        raise ValueError(f"unexpected CIF asymmetric-site table: {sites}")
    site_multiplicities = [2, 4, 2, 4]
    for i, site in enumerate(sites):
        if source_kind_counts.get(i) != site_multiplicities[i]:
            raise ValueError(f"unexpected expanded multiplicity for CIF site {site}")
    source_parameters = {
        "COD_record_internal": "9009647",
        "AMCSD_internal": "0013252",
        "DOI_internal": "10.1016/S0022-3697(98)00226-1",
        "title_internal": "Structure and crystal chemistry of Ti3SiC2",
        "journal_internal": "Journal of Physics and Chemistry of Solids",
        "year_internal": 1998,
        "spacegroup_number_internal": int(source.info["spacegroup"].no),
        "spacegroup_HM_internal": cif_scalar(cif_text, "_symmetry_space_group_name_H-M"),
        "a_A": float(source.cell.lengths()[0]),
        "b_A": float(source.cell.lengths()[1]),
        "c_A": float(source.cell.lengths()[2]),
        "angles_deg": [float(x) for x in source.cell.angles()],
        "Z_formula_units": 2,
        "source_expanded_atoms": len(source),
        "source_expanded_formula": source.get_chemical_formula(),
        "expanded_element_counts": formula_counts(source.get_chemical_symbols()),
        "expanded_site_kind_counts": {
            str(int(k)): int(v) for k, v in collections.Counter(source.arrays["spacegroup_kinds"]).items()
        },
        "asymmetric_site_rows": [
            {**site, "multiplicity": site_multiplicities[i]}
            for i, site in enumerate(sites)
        ],
        "occupancy_column_present": "_atom_site_occupancy" in cif_text,
        "internal_source_url": "https://www.crystallography.net/cod/9009647.cif",
        "remote_COD_checksum_verified": False,
        "uploaded_source_sha256": sha256(cif),
    }
    required_tags = {
        "_cod_database_code": "9009647",
        "_database_code_amcsd": "0013252",
        "_chemical_formula_structural": "Ti3SiC2",
        "_space_group_IT_number": "194",
        "_symmetry_space_group_name_H-M": "P 63/m m c",
        "_journal_year": "1998",
        "_journal_paper_doi": "10.1016/S0022-3697(98)00226-1",
        "_cell_length_a": "3.0575",
        "_cell_length_b": "3.0575",
        "_cell_length_c": "17.6235",
    }
    tag_checks = {}
    for tag, value in required_tags.items():
        observed = cif_scalar(cif_text, tag)
        tag_checks[tag] = {"expected": value, "observed": observed, "pass": observed == value}
    if not all(x["pass"] for x in tag_checks.values()):
        raise ValueError(f"CIF internal metadata mismatch: {tag_checks}")
    source_parameters["required_internal_tag_checks"] = tag_checks
    if len(source) != 12 or source.get_chemical_formula() != "C4Si2Ti6":
        raise ValueError("CIF symmetry expansion did not produce 12-atom C4Si2Ti6 conventional cell")
    if int(source.info["spacegroup"].no) != 194:
        raise ValueError("ASE did not preserve CIF space group 194")

    sums = []
    for line in (SRC / "SHA256SUMS.txt").read_text().splitlines():
        expected, name = line.split(None, 1)
        path = SRC / name.strip()
        observed = sha256(path)
        sums.append({"path": str(path.relative_to(REPO)), "expected": expected,
                     "observed": observed, "pass": observed == expected})
    if not all(x["pass"] for x in sums):
        raise ValueError("source bundle SHA256SUMS check failed")
    if sha256(cif) != manifest["source_cif_sha256"]:
        raise ValueError("manifest source CIF hash mismatch")

    regenerated = reproduce_builder()
    if not regenerated["manifest_byte_identical"] or not all(
            x["byte_identical"] for x in regenerated["candidate_files"]):
        raise ValueError("isolated builder outputs differ bytewise")

    training_manifest = json.loads(TRAIN_MANIFEST.read_text())
    if sha256(TRAIN) != training_manifest["train_sha256"]:
        raise ValueError("training36 file hash does not match its manifest")
    training = read_geometry_frames(TRAIN)
    if len(training) != 36:
        raise ValueError(f"expected 36 training frames, got {len(training)}")
    dev_paths = [DEV_DIR / x for x in DEV_NAMES]
    dev_frames = []
    for path in dev_paths:
        frames = read_geometry_frames(path)
        if len(frames) != 1:
            raise ValueError(f"expected one geometry-only frame in {path}")
        dev_frames.append(frames[0])

    proposals = []
    for interface, (contact, candidate_name, old_name) in INTERFACES.items():
        record = next(x for x in manifest["records"] if x["interface"] == interface)
        if record["sha256"] != sha256(SRC / "inputs" / candidate_name):
            raise ValueError(f"candidate hash mismatch for {interface}")
        proposal = proposal_metrics(interface, contact, candidate_name, old_name, record, source.cell.array)
        proposal["source_substrate_geometry_check"] = verify_substrate_against_cif(
            contact, candidate_name, proposal, source)
        formula = proposal["composition"]
        same_composition_train = [i for i, x in enumerate(training) if x["composition"] == formula]
        same_composition_dev = [
            {"file": str(path.relative_to(REPO)), "frame": frame}
            for path, frame in zip(dev_paths, dev_frames) if frame["composition"] == formula
        ]
        proposal["training36_exact_or_near_composition_candidates"] = same_composition_train
        proposal["v17_development_exact_or_near_composition_candidates"] = same_composition_dev
        proposal["whole_structure_duplicate_evaluation"] = (
            "not applicable: no equal-composition baseline frame; no exact/near whole-structure match"
            if not same_composition_train and not same_composition_dev else
            "requires separate same-composition coordinate matcher"
        )
        proposal_atoms = read(SRC / "inputs" / candidate_name, format="extxyz")
        proposal["periodic_contact_fingerprint"] = geometry_fingerprint({
            "symbols": proposal_atoms.get_chemical_symbols(),
            "positions": proposal_atoms.positions,
            "cell": proposal_atoms.cell.array,
            "pbc": proposal_atoms.pbc.tolist(),
        })
        proposals.append(proposal)

    training_topology = [geometry_fingerprint(frame) for frame in training]
    dev_topology = [geometry_fingerprint(frame) for frame in dev_frames]

    cell_deltas = {
        "a_relative_to_old_assumption_percent": 100.0 * (source.cell.lengths()[0] - old_params["a_A"]) / old_params["a_A"],
        "c_relative_to_old_assumption_percent": 100.0 * (source.cell.lengths()[2] - old_params["c_A"]) / old_params["c_A"],
        "old_Ti4f_z": old_params["Ti_4f_z"],
        "CIF_Ti4f_z": 0.13550,
        "old_C4f_z": old_params["C_4f_z"],
        "CIF_C4f_z": 0.07220,
        "note": "Compare full symmetry-expanded site clouds with a common origin shift; raw C z scalars use a different orbit representative and are not a direct displacement.",
    }
    result = {
        "status": "SOURCE_FILE_AND_REBUILD_VERIFIED; PROPOSALS_REMAIN_DEVELOPMENT_ONLY",
        "source_evidence_level": "Uploaded CIF is internally self-attributed and internally consistent; remote COD bytes/checksum were not independently retrieved.",
        "source_parameters": source_parameters,
        "source_bundle_hash_checks": sums,
        "rebuild_reproduction": regenerated,
        "old_assumed_vs_CIF_bulk_origin_match": fit_bulk_origin_translation(source, assumed),
        "old_assumed_vs_CIF_parameter_comparison": cell_deltas,
        "training36_geometry_only": {
            "path": str(TRAIN.relative_to(REPO)),
            "sha256": sha256(TRAIN),
            "frame_count": len(training),
            "composition_counts": {
                json.dumps(k, sort_keys=True): sum(x["composition"] == k for x in training)
                for k in sorted({tuple(sorted(x["composition"].items())) for x in training})
                for k in [dict(k)]
            },
            "cell_length_min_max_A": [
                float(min(np.linalg.norm(x["cell"], axis=1).min() for x in training)),
                float(max(np.linalg.norm(x["cell"], axis=1).max() for x in training)),
            ],
            "short_periodic_image_pairs_using_3p3A_diagnostic_cutoff": {
                "frames_with_any": sum(x["short_pairs_using_periodic_images_total"] > 0 for x in training_topology),
                "min": min(x["short_pairs_using_periodic_images_total"] for x in training_topology),
                "median": float(np.median([x["short_pairs_using_periodic_images_total"] for x in training_topology])),
                "max": max(x["short_pairs_using_periodic_images_total"] for x in training_topology),
            },
            "diagnostic_note": "Pair counts use geometry only; 3.3 A is a descriptive fingerprint cutoff, not a chemical rejection threshold.",
            "labels_or_force_values_read": False,
        },
        "reused_v17_development_geometry_only": [
            {"file": str(path.relative_to(REPO)), "sha256": sha256(path),
             "atoms": frame["n"], "composition": frame["composition"],
             "cell_A": frame["cell"].tolist(), "pbc": frame["pbc"],
             "periodic_contact_fingerprint": topology,
             "labels_or_force_values_read": False}
            for path, frame, topology in zip(dev_paths, dev_frames, dev_topology)
        ],
        "proposals": proposals,
        "cross_candidate_lineage": {
            "all_three_share_one_CIF_bulk_parent": True,
            "same_formula": True,
            "distinct_development_inputs": True,
            "independent_source_families": False,
        },
        "sealed_data_use_in_final_audit": "No sealed/test inputs, outputs, energies, or forces are used in this audit script or the reported comparisons.",
        "scope_incident": {
            "occurred": True,
            "description": "An earlier exploratory glob over v17_split_geometry_design/final_set/geometries included test-named geometry files. A temporary parser read each atom row and retained species plus xyz for formula summaries; it split complete text rows but did not evaluate or report force/energy values. The test geometries were not used in the final comparison or candidate decisions and were not written to artifacts. The final audit script uses only the three explicitly named V17 development input geometries and training36.",
            "consequence": "This B review session has accessed test geometry. A should treat this as a geometry-scope incident when deciding whether those structures remain eligible as sealed test inputs; the V19 proposal conclusions below do not use their contents.",
        },
    }
    (OUT / "audit.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "source_hash": sha256(cif),
        "bulk_origin_match": result["old_assumed_vs_CIF_bulk_origin_match"],
        "rebuild_byte_identical": all(x["byte_identical"] for x in regenerated["candidate_files"]),
        "training36_frames": len(training),
        "dev_frames": len(dev_frames),
        "proposal_summary": [{"interface": x["file"].split("/")[-1],
                              "composition": x["composition"],
                              "gap_A": x["vertical_Ag_to_substrate_gap_A"],
                              "vacuum_A": x["vacuum_from_z_extent_A"],
                              "Ag_contact_A": x["marked_pair_distance_A"],
                              "old_contact_A": x["old_prototype_marked_contact_minimum_A"]}
                             for x in proposals],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
