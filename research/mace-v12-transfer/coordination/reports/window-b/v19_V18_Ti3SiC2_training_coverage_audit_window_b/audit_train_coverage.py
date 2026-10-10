#!/usr/bin/env python3
"""Read-only V18 training-set coverage audit; no model inference or label acquisition."""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list


def find_repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("Cannot locate repository root")


REPO = find_repo_root(Path(__file__).resolve())
MODEL_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core")
TRAIN_REL = MODEL_ROOT / "data/train.extxyz"
SELECTION_REL = MODEL_ROOT / "training_clean_core_01/selection_record.json"
MANIFEST_REL = MODEL_ROOT / "data/dataset_manifest.json"
CONTROL_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs")
SCREEN_METRICS_REL = Path(
    "research/mace-v12-transfer/coordination/reports/window-a/"
    "v19_V18_pure_phase_force_screen/metrics.json"
)
CONTROLS = {
    "Ag_baseline_k8x8x8_v19": CONTROL_ROOT / "Ag_baseline_PROPOSAL.extxyz",
    "Ti3SiC2_baseline_k10x10x2_v19": CONTROL_ROOT / "Ti3SiC2_baseline_PROPOSAL.extxyz",
}
TIC_FIRST_SHELL_CUTOFF_A = 2.50
PAIR_DISTANCE_CUTOFF_A = 4.00
COORDINATION_CUTOFFS_A = (2.50, 3.50)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_file(relative: Path) -> Path:
    path = (REPO / relative).resolve()
    if REPO.resolve() not in path.parents or not path.is_file():
        raise RuntimeError(f"Missing or out-of-repository input: {relative}")
    return path


def geometry_sha256(atoms) -> str:
    payload = b"".join(
        (
            np.asarray(atoms.numbers, dtype="<i4").tobytes(),
            np.asarray(atoms.positions, dtype="<f8").tobytes(),
            np.asarray(atoms.cell.array, dtype="<f8").tobytes(),
            np.asarray(atoms.pbc, dtype="u1").tobytes(),
        )
    )
    return hashlib.sha256(payload).hexdigest()


def json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def family_from_config(config_type: str) -> str:
    if config_type.startswith("AgC"):
        return "Ag-C"
    if config_type.startswith("AgSi"):
        return "Ag-Si"
    if config_type.startswith("AgTi"):
        return "Ag-Ti"
    if config_type.startswith("Ti_gap"):
        return "Ti-gap"
    return "other"


def series_from_config(config_type: str) -> str:
    value = config_type.lower()
    if "ti_gap" in value:
        return "short-range Ti gap"
    if "framework_neighbor" in value:
        return "framework-neighbor displacement pair"
    if "transverse_ag" in value:
        return "transverse Ag displacement pair"
    if "widespan_distance" in value:
        return "wide-span distance"
    if "widespan_registry" in value:
        return "wide-span registry"
    if "residual_shell" in value:
        return "targeted residual shell"
    if "registry_strain" in value:
        return "registry strain"
    if "registry_probe" in value:
        return "registry probe"
    if "lateral_registry" in value:
        return "lateral registry"
    if "validation" in value:
        return "validation-derived frame in train.extxyz"
    if "_d2p" in value:
        return "contact-distance scan"
    return "other"


def quantiles(values: list[float]) -> dict:
    if not values:
        return {"count": 0, "min": None, "p10": None, "median": None, "p90": None, "max": None}
    arr = np.asarray(values, dtype=float)
    return {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "p10": float(np.quantile(arr, 0.10)),
        "median": float(np.median(arr)),
        "p90": float(np.quantile(arr, 0.90)),
        "max": float(np.max(arr)),
    }


def pair_distance_map(atoms, max_distance: float = PAIR_DISTANCE_CUTOFF_A) -> dict:
    symbols = np.asarray(atoms.get_chemical_symbols())
    distances = atoms.get_all_distances(mic=True)
    pair_types = (("Ti", "C"), ("Ti", "Si"), ("Ti", "Ti"), ("C", "Ag"), ("Ti", "Ag"), ("Ag", "Si"), ("Ag", "C"))
    output = {}
    for first, second in pair_types:
        first_indices = np.where(symbols == first)[0]
        second_indices = np.where(symbols == second)[0]
        if first == second:
            block = distances[np.ix_(first_indices, second_indices)]
            vals = block[np.triu_indices_from(block, k=1)]
        elif first_indices.size and second_indices.size:
            vals = distances[np.ix_(first_indices, second_indices)].ravel()
        else:
            vals = np.asarray([], dtype=float)
        vals = vals[(vals > 1.0e-8) & (vals <= max_distance)]
        output[f"{first}-{second}"] = quantiles(vals.tolist())
    return output


def aggregate_pair_distances(atoms_list: list) -> dict:
    pair_types = (("Ti", "C"), ("Ti", "Si"), ("Ti", "Ti"), ("C", "Ag"), ("Ti", "Ag"), ("Ag", "Si"), ("Ag", "C"))
    collected = {f"{first}-{second}": [] for first, second in pair_types}
    for atoms in atoms_list:
        symbols = np.asarray(atoms.get_chemical_symbols())
        distances = atoms.get_all_distances(mic=True)
        for first, second in pair_types:
            first_indices = np.where(symbols == first)[0]
            second_indices = np.where(symbols == second)[0]
            if first == second:
                block = distances[np.ix_(first_indices, second_indices)]
                vals = block[np.triu_indices_from(block, k=1)]
            elif first_indices.size and second_indices.size:
                vals = distances[np.ix_(first_indices, second_indices)].ravel()
            else:
                vals = np.asarray([], dtype=float)
            vals = vals[(vals > 1.0e-8) & (vals <= PAIR_DISTANCE_CUTOFF_A)]
            collected[f"{first}-{second}"].extend(vals.tolist())
    return {key: quantiles(values) for key, values in collected.items()}


def ti_c_local(atoms) -> dict:
    symbols = np.asarray(atoms.get_chemical_symbols())
    distances = atoms.get_all_distances(mic=True)
    ti_indices = np.where(symbols == "Ti")[0]
    c_indices = np.where(symbols == "C")[0]
    if not ti_indices.size or not c_indices.size:
        return {
            "ti_count": int(ti_indices.size),
            "c_count": int(c_indices.size),
            "ti_nearest_C_A": quantiles([]),
            "c_nearest_Ti_A": quantiles([]),
            "min_image_pair_distances_Ti_C_le4A": quantiles([]),
            "first_shell_coordination": {},
        }

    block = distances[np.ix_(ti_indices, c_indices)]
    ti_nearest = np.min(block, axis=1).tolist()
    c_nearest = np.min(block, axis=0).tolist()
    pair_values = block[(block > 1.0e-8) & (block <= PAIR_DISTANCE_CUTOFF_A)].tolist()

    # ASE's neighbor list includes distinct periodic images, which matters for
    # the one short-cell Ti-gap training frame.
    i_idx, j_idx, dists = neighbor_list(
        "ijd", atoms, cutoff=max(COORDINATION_CUTOFFS_A), self_interaction=False
    )
    ti_coord = Counter()
    c_coord = Counter()
    ti_seen = {int(idx): {cut: 0 for cut in COORDINATION_CUTOFFS_A} for idx in ti_indices}
    c_seen = {int(idx): {cut: 0 for cut in COORDINATION_CUTOFFS_A} for idx in c_indices}
    tic_pairs = {cut: 0 for cut in COORDINATION_CUTOFFS_A}
    for i, j, distance in zip(i_idx, j_idx, dists):
        if symbols[i] == "Ti" and symbols[j] == "C":
            for cut in COORDINATION_CUTOFFS_A:
                if distance <= cut:
                    ti_seen[int(i)][cut] += 1
                    c_seen[int(j)][cut] += 1
                    tic_pairs[cut] += 1
    for cut in COORDINATION_CUTOFFS_A:
        ti_coord[cut] = Counter(row[cut] for row in ti_seen.values())
        c_coord[cut] = Counter(row[cut] for row in c_seen.values())
    return {
        "ti_count": int(ti_indices.size),
        "c_count": int(c_indices.size),
        "ti_nearest_C_A": quantiles(ti_nearest),
        "c_nearest_Ti_A": quantiles(c_nearest),
        "min_image_pair_distances_Ti_C_le4A": quantiles(pair_values),
        "first_shell_coordination": {
            str(cut): {
                "cutoff_A": cut,
                "Ti_neighbor_count_histogram": {str(k): int(v) for k, v in sorted(ti_coord[cut].items())},
                "C_neighbor_count_histogram": {str(k): int(v) for k, v in sorted(c_coord[cut].items())},
                "periodic_Ti_C_pair_count": int(tic_pairs[cut]),
            }
            for cut in COORDINATION_CUTOFFS_A
        },
    }


def force_rows(groups: dict[str, list], frames: list) -> list[dict]:
    rows = []
    for group, indices in groups.items():
        grouped_frames = [frames[index] for index in indices]
        for element in ("Ag", "Ti", "Si", "C"):
            samples = []
            for atoms in grouped_frames:
                symbols = atoms.get_chemical_symbols()
                element_indices = [i for i, symbol in enumerate(symbols) if symbol == element]
                if element_indices:
                    samples.append(np.asarray(atoms.arrays["REF_forces"], dtype=float)[element_indices])
            if not samples:
                continue
            values = np.concatenate(samples, axis=0)
            vector_norms = np.linalg.norm(values, axis=1)
            for component_index, component in enumerate("xyz"):
                comp = values[:, component_index]
                rows.append(
                    {
                        "group": group,
                        "species": element,
                        "component": component,
                        "atom_samples": int(len(values)),
                        "frame_count": len(grouped_frames),
                        "mean_eV_A": float(np.mean(comp)),
                        "RMS_eV_A": float(np.sqrt(np.mean(comp**2))),
                        "median_abs_eV_A": float(np.median(np.abs(comp))),
                        "p95_abs_eV_A": float(np.quantile(np.abs(comp), 0.95)),
                        "max_abs_eV_A": float(np.max(np.abs(comp))),
                        "positive_count": int(np.sum(comp > 0.0)),
                        "negative_count": int(np.sum(comp < 0.0)),
                        "abs_Fz_gt_0p1_count": int(np.sum(np.abs(values[:, 2]) > 0.1)) if component == "z" else "",
                        "abs_Fz_gt_0p5_count": int(np.sum(np.abs(values[:, 2]) > 0.5)) if component == "z" else "",
                        "vector_RMS_eV_A": float(np.sqrt(np.mean(vector_norms**2))),
                    }
                )
    return rows


def build_run() -> dict:
    selection_path = repo_file(SELECTION_REL)
    manifest_path = repo_file(MANIFEST_REL)
    train_path = repo_file(TRAIN_REL)
    selection = json.loads(selection_path.read_text())
    manifest = json.loads(manifest_path.read_text())

    # Validate lineage before ASE reads/parses the training frames.
    expected_sha = selection.get("training_input_sha256")
    manifest_sha = manifest.get("train_input_sha256")
    actual_sha = file_sha256(train_path)
    hash_checks = {
        "selected_model_record_matches_train": expected_sha == actual_sha,
        "dataset_manifest_matches_train": manifest_sha == actual_sha,
        "selected_epoch_is_79": selection.get("selected_epoch") == 79,
        "selection_blind_labels_unused": selection.get("blind_labels_used_for_selection") is False,
        "selection_test_file_not_passed_to_training": selection.get("test_file_passed_to_training") is False,
        "dataset_manifest_withheld_labels_not_opened": manifest.get("withheld_labels_opened") is False,
    }
    if not all(hash_checks.values()):
        raise RuntimeError(f"Training input lineage check failed before parsing: {hash_checks}")

    # The only training structures loaded by this audit.
    frames = read(train_path, index=":")
    if not frames:
        raise RuntimeError("V18 training file contains no frames")
    for index, atoms in enumerate(frames):
        if "REF_forces" not in atoms.arrays:
            raise RuntimeError(f"Missing REF_forces in V18 train frame {index}")
        force = np.asarray(atoms.arrays["REF_forces"], dtype=float)
        if force.shape != (len(atoms), 3) or not np.isfinite(force).all():
            raise RuntimeError(f"Invalid REF_forces in V18 train frame {index}")

    screen_metrics_path = repo_file(SCREEN_METRICS_REL)
    screen_metrics = json.loads(screen_metrics_path.read_text())
    screen_by_label = {row["label"]: row for row in screen_metrics["references"]}
    control_rows = {}
    control_atoms = {}
    for label, relative in CONTROLS.items():
        source = repo_file(relative)
        if label not in screen_by_label or file_sha256(source) != screen_by_label[label].get("DFT_source_sha256"):
            raise RuntimeError(f"Public control geometry hash mismatch for {label}")
        if str(relative) != screen_by_label[label].get("DFT_source_path"):
            raise RuntimeError(f"Unexpected public control geometry path for {label}")
        atoms = read(source, index=0)
        if geometry_sha256(atoms) != screen_by_label[label].get("result_geometry_sha256"):
            raise RuntimeError(f"Public control geometry differs from A-screened geometry for {label}")
        control_atoms[label] = atoms
        control_rows[label] = {
            "label": label,
            "source_path": str(relative),
            "source_sha256": file_sha256(source),
            "geometry_sha256": geometry_sha256(atoms),
            "formula": atoms.get_chemical_formula(),
            "atom_count": len(atoms),
            "species_counts": {key: int(value) for key, value in Counter(atoms.get_chemical_symbols()).items()},
            "pbc": atoms.pbc.tolist(),
            "cell_A": atoms.cell.array.tolist(),
            "cell_lengths_A": atoms.cell.lengths().tolist(),
            "cell_angles_deg": atoms.cell.angles().tolist(),
            "screen_force_RMSE_eV_A": screen_by_label[label]["all_atom_force_vector_RMSE_eV_A"],
            "screen_force_component_RMSE_eV_A": screen_by_label[label]["force_component_RMSE_eV_A"],
            "screen_per_species_force_RMSE_eV_A": {
                species: values["force_vector_RMSE_eV_A"]
                for species, values in screen_by_label[label]["per_species"].items()
            },
            "geometry_only_read": True,
        }

    formula_counts = Counter(atoms.get_chemical_formula() for atoms in frames)
    total_species = Counter(symbol for atoms in frames for symbol in atoms.get_chemical_symbols())
    family_indices = defaultdict(list)
    series_counts = Counter()
    frame_rows = []
    provenance_rows = []
    unique_configs = set()
    source_dataset_paths = Counter()
    source_refs = Counter()
    parent_geometry_hashes = Counter()
    source_geometry_hashes = Counter()
    role_counts = Counter()
    pbc_counts = Counter()
    cell_class_counts = Counter()

    for index, atoms in enumerate(frames):
        info = atoms.info
        config_type = str(info.get("config_type", "unknown"))
        unique_configs.add(config_type)
        family = family_from_config(config_type)
        series = series_from_config(config_type)
        family_indices[family].append(index)
        series_counts[series] += 1
        role_counts[str(info.get("dataset_role", "not_recorded"))] += 1
        pbc_counts[str(tuple(bool(x) for x in atoms.pbc))] += 1
        lengths = np.asarray(atoms.cell.lengths(), dtype=float)
        angles = np.asarray(atoms.cell.angles(), dtype=float)
        if np.all(np.abs(lengths - 28.0) <= 0.15) and np.all(np.abs(angles - 90.0) <= 0.2):
            cell_class = "large approximately cubic 28A box"
        else:
            cell_class = "compact/noncubic cell"
        cell_class_counts[cell_class] += 1
        counts = Counter(atoms.get_chemical_symbols())
        non_ag = {element: int(counts.get(element, 0)) for element in ("Ti", "Si", "C")}
        if non_ag["Si"]:
            ti_ratio_error = non_ag["Ti"] / (3.0 * non_ag["Si"]) - 1.0
            c_ratio_error = non_ag["C"] / (2.0 * non_ag["Si"]) - 1.0
        else:
            ti_ratio_error = None
            c_ratio_error = None
        local = ti_c_local(atoms)
        first_shell = local["first_shell_coordination"][str(TIC_FIRST_SHELL_CUTOFF_A)]
        row = {
            "frame_index": index,
            "config_type": config_type,
            "interface_family": family,
            "geometry_series": series,
            "formula": atoms.get_chemical_formula(),
            "atom_count": len(atoms),
            "species_counts": {key: int(value) for key, value in sorted(counts.items())},
            "contains_Ag_Ti_Si_C": all(counts.get(element, 0) > 0 for element in ("Ag", "Ti", "Si", "C")),
            "exact_pure_Ti3SiC2_formula": counts.get("Ag", 0) == 0 and non_ag["Ti"] == 3 * non_ag["Si"] and non_ag["C"] == 2 * non_ag["Si"],
            "Ti_per_3Si_relative_error": ti_ratio_error,
            "C_per_2Si_relative_error": c_ratio_error,
            "pbc": atoms.pbc.tolist(),
            "cell_lengths_A": lengths.tolist(),
            "cell_angles_deg": angles.tolist(),
            "cell_class": cell_class,
            "dataset_role": info.get("dataset_role"),
            "proposal_role": info.get("proposal_role"),
            "proposal_provenance": info.get("proposal_provenance"),
            "geometry_family_metadata": info.get("geometry_family"),
            "REF_forces_present": "REF_forces" in atoms.arrays,
            "PW_PBE_forces_present": "PW_PBE_forces" in atoms.arrays,
            "REF_forces_vector_RMS_eV_A": float(np.sqrt(np.mean(np.sum(np.asarray(atoms.arrays["REF_forces"]) ** 2, axis=1)))),
            "Ti_C_nearest_per_Ti_A": local["ti_nearest_C_A"],
            "C_Ti_nearest_per_C_A": local["c_nearest_Ti_A"],
            "Ti_C_first_shell_Ti_coordination_hist": first_shell["Ti_neighbor_count_histogram"],
            "Ti_C_first_shell_C_coordination_hist": first_shell["C_neighbor_count_histogram"],
            "Ti_C_periodic_pair_count_le2p5A": first_shell["periodic_Ti_C_pair_count"],
        }
        frame_rows.append(row)

        source_type = info.get("source_config_type")
        if isinstance(source_type, str) and source_type.startswith("source_structure_sha256="):
            source_type = None
        source_key = info.get("source_label") or source_type or info.get("parent_config_type") or config_type
        provenance = {
            "frame_index": index,
            "config_type": config_type,
            "interface_family": family,
            "source_record_key_from_train_metadata": source_key,
            "source_label": info.get("source_label"),
            "source_config_type": info.get("source_config_type"),
            "source_dataset_path": info.get("source_dataset_path"),
            "source_dataset_sha256": info.get("source_dataset_sha256"),
            "source_structure_sha256": info.get("source_structure_sha256"),
            "source_geometry_sha256": info.get("source_geometry_sha256"),
            "source_frame_geometry_sha256": info.get("source_frame_geometry_sha256"),
            "parent_config_type": info.get("parent_config_type"),
            "parent_file_sha256": info.get("parent_file_sha256"),
            "parent_geometry_sha256": info.get("parent_geometry_sha256"),
            "previous_role_metadata": info.get("previous_role"),
            "validation_role_metadata": info.get("validation_role"),
            "proposal_provenance": info.get("proposal_provenance"),
            "proposal_role": info.get("proposal_role"),
        }
        provenance_rows.append(provenance)
        source_dataset_paths[str(info.get("source_dataset_path", "not_recorded"))] += 1
        source_refs[str(source_key)] += 1
        if info.get("parent_geometry_sha256"):
            parent_geometry_hashes[str(info["parent_geometry_sha256"])] += 1
        for key in ("source_geometry_sha256", "source_structure_sha256", "source_frame_geometry_sha256"):
            if info.get(key):
                source_geometry_hashes[str(info[key])] += 1

    groups = {"all_train": list(range(len(frames)))}
    groups.update(family_indices)
    forces = force_rows(groups, frames)

    coordination_rows = []
    coordination_data = {}
    training_pair_distances = {}
    for group, indices in groups.items():
        result = aggregate_ti_c([frames[index] for index in indices])
        coordination_data[group] = result
        training_pair_distances[group] = aggregate_pair_distances([frames[index] for index in indices])
        for center in ("Ti", "C"):
            coord = result["first_shell_coordination"][str(TIC_FIRST_SHELL_CUTOFF_A)]
            histogram = coord["Ti_neighbor_count_histogram"] if center == "Ti" else coord["C_neighbor_count_histogram"]
            coordination_rows.append(
                {
                    "group": group,
                    "center_species": center,
                    "frame_count": result["frame_count"],
                    "center_atom_count": result["ti_count"] if center == "Ti" else result["c_count"],
                    "neighbor_species": "C" if center == "Ti" else "Ti",
                    "first_shell_cutoff_A": TIC_FIRST_SHELL_CUTOFF_A,
                    "coordination_histogram": histogram,
                    "Ti_C_periodic_pair_count_le2p5A": coord["periodic_Ti_C_pair_count"],
                    "nearest_opposite_species_A": result["ti_nearest_C_A"] if center == "Ti" else result["c_nearest_Ti_A"],
                    "Ti_C_pair_distance_min_image_A_le4": result["min_image_pair_distances_Ti_C_le4A"],
                }
            )

    control_coordination = {}
    control_distance_pairs = {}
    for label, atoms in control_atoms.items():
        if "Ti" in atoms.get_chemical_symbols() and "C" in atoms.get_chemical_symbols():
            control_coordination[label] = aggregate_ti_c([atoms])
            control_distance_pairs[label] = pair_distance_map(atoms)
        else:
            control_coordination[label] = {"frame_count": 1, "ti_count": 0, "c_count": 0}
            control_distance_pairs[label] = pair_distance_map(atoms)

    formula_counts = Counter(atoms.get_chemical_formula() for atoms in frames)
    all_four = sum(row["contains_Ag_Ti_Si_C"] for row in frame_rows)
    pbc_all = sum(bool(np.all(atoms.pbc)) for atoms in frames)
    force_array_counts = Counter(
        "REF_and_PW_PBE" if "PW_PBE_forces" in atoms.arrays else "REF_only"
        for atoms in frames
    )
    source_dataset_paths.pop("not_recorded", None)
    stoich_candidates = [
        {
            "frame_index": row["frame_index"],
            "config_type": row["config_type"],
            "formula": row["formula"],
            "Ti_per_3Si_relative_error": row["Ti_per_3Si_relative_error"],
            "C_per_2Si_relative_error": row["C_per_2Si_relative_error"],
        }
        for row in frame_rows
        if row["Ti_per_3Si_relative_error"] is not None
    ]
    stoich_candidates.sort(
        key=lambda row: max(abs(row["Ti_per_3Si_relative_error"]), abs(row["C_per_2Si_relative_error"]))
    )

    return {
        "task": "v19_V18_Ti3SiC2_training_coverage_audit_window_b",
        "status": "READ_ONLY_COVERAGE_AUDIT_COMPLETE",
        "scope": "Read-only structural/label coverage inventory of the selected V18 train.extxyz, its selection/dataset hash records, its own embedded per-frame provenance metadata, and the two A-screened public control geometries. No MACE inference, training, DFT/MPI, dataset edits, development/test/holdout/sealed structures or labels were accessed.",
        "lineage": {
            "train_path": str(TRAIN_REL),
            "train_sha256": actual_sha,
            "selection_record_path": str(SELECTION_REL),
            "selection_record_sha256": file_sha256(selection_path),
            "dataset_manifest_path": str(MANIFEST_REL),
            "dataset_manifest_sha256": file_sha256(manifest_path),
            "selected_epoch": selection.get("selected_epoch"),
            "selected_model_sha256": selection.get("model_sha256"),
            "expected_train_sha256_in_selection": expected_sha,
            "expected_train_sha256_in_dataset_manifest": manifest_sha,
            "checks": hash_checks,
            "train_file_parsed_after_hash_pass": True,
        },
        "dataset_summary": {
            "frame_count": len(frames),
            "total_atom_samples": int(sum(total_species.values())),
            "formula_counts": dict(sorted(formula_counts.items())),
            "species_atom_counts": {key: int(value) for key, value in sorted(total_species.items())},
            "frame_counts_by_interface_family": {key: len(value) for key, value in sorted(family_indices.items())},
            "frame_counts_by_geometry_series": dict(sorted(series_counts.items())),
            "unique_config_type_count": len(unique_configs),
            "all_frames_have_Ag_Ti_Si_C": all_four == len(frames),
            "exact_pure_Ti3SiC2_formula_frame_count": sum(row["exact_pure_Ti3SiC2_formula"] for row in frame_rows),
            "frames_with_Ti3SiC2_nonAg_ratios_within_10pct": sum(
                max(abs(row["Ti_per_3Si_relative_error"]), abs(row["C_per_2Si_relative_error"])) <= 0.10
                for row in stoich_candidates
            ),
            "closest_Ti3SiC2_nonAg_composition_ratios": stoich_candidates[:5],
            "pbc_counts": dict(sorted(pbc_counts.items())),
            "fully_periodic_frame_count": int(pbc_all),
            "cell_class_counts": dict(sorted(cell_class_counts.items())),
            "cell_vector_length_minmax_A": {
                f"v{axis}": [float(min(np.asarray(a.cell.lengths())[axis] for a in frames)), float(max(np.asarray(a.cell.lengths())[axis] for a in frames))]
                for axis in range(3)
            },
            "REF_forces_frame_count": len(frames),
            "force_array_presence_counts": dict(force_array_counts),
            "frames_with_frame_metadata_dataset_role": int(sum("dataset_role" in a.info for a in frames)),
            "frames_with_source_dataset_path": int(sum("source_dataset_path" in a.info for a in frames)),
            "frames_with_source_label_or_type": int(sum(any(k in a.info for k in ("source_label", "source_config_type")) for a in frames)),
            "frames_with_parent_geometry_hash": int(sum("parent_geometry_sha256" in a.info for a in frames)),
            "frames_marked_direct_archive_verified": int(sum(a.info.get("proposal_provenance") == "directly_archive_verified" for a in frames)),
            "frames_marked_proposal_role_train": int(sum(a.info.get("proposal_role") == "train" for a in frames)),
            "source_record_key_counts": dict(sorted(source_refs.items())),
            "source_dataset_path_counts": dict(sorted(source_dataset_paths.items())),
            "repeated_parent_geometry_hash_counts": {key: int(value) for key, value in parent_geometry_hashes.items() if value > 1},
            "source_geometry_hash_occurrence_counts": {key: int(value) for key, value in source_geometry_hashes.items()},
            "previous_role_metadata_counts_within_train_only": dict(Counter(str(a.info.get("previous_role", "not_recorded")) for a in frames)),
        },
        "reference_geometry_comparison": {
            "control_geometries": control_rows,
            "training_all_ti_c_coordination": coordination_data["all_train"],
            "training_interface_group_ti_c_coordination": {key: coordination_data[key] for key in sorted(family_indices)},
            "training_pair_distances_within4A_by_group": training_pair_distances,
            "control_ti_c_coordination": control_coordination,
            "control_pair_distances_within4A": control_distance_pairs,
            "first_shell_definition": {
                "Ti_C_cutoff_A": TIC_FIRST_SHELL_CUTOFF_A,
                "basis": "The A-screened Ti3SiC2 control has its first Ti-C shell at 2.088-2.176 A; the next Ti-C shell starts above 3.6 A. 2.50 A is inside the observed gap.",
            },
        },
        "training_REF_forces_distribution": {
            "label_key": "REF_forces",
            "force_distribution_rows": forces,
            "interpretation_limit": "Large z force components exist in the interface training labels; their presence does not establish matched pristine Ti3SiC2 backbone-response coverage or explain the V18 residual.",
        },
        "public_control_screen_metrics_from_A": {
            "metrics_path": str(SCREEN_METRICS_REL),
            "metrics_sha256": file_sha256(screen_metrics_path),
            "controls": [
                {
                    "label": row["label"],
                    "all_atom_force_vector_RMSE_eV_A": row["screen_force_RMSE_eV_A"],
                    "force_component_RMSE_eV_A": row["screen_force_component_RMSE_eV_A"],
                    "per_species_force_vector_RMSE_eV_A": row["screen_per_species_force_RMSE_eV_A"],
                }
                for row in control_rows.values()
            ],
        },
        "coverage_diagnosis": "The selected training set contains many Ti-C contacts and substantial REF_forces z components, but all 30 frames are Ag-containing interface/contact structures, no frame has the pure Ti3SiC2 formula, and first-shell Ti/C coordination is less complete and more dispersed than in the screened bulk control. This supports a coverage gap hypothesis; counts alone do not prove the cause of the model residual.",
        "frame_inventory": frame_rows,
        "training_provenance": provenance_rows,
        "coordination_csv_rows": coordination_rows,
        "force_csv_rows": forces,
        "prohibited_data_access": {
            "development_structures_or_labels_read": False,
            "test_structures_or_labels_read": False,
            "holdout_or_sealed_structures_or_labels_read": False,
            "other_training_dataset_frames_opened": False,
            "DFT_or_MPI_run": False,
            "MACE_inference_run": False,
            "training_run": False,
        },
    }


def aggregate_ti_c(atoms_list: list) -> dict:
    ti_nearest = []
    c_nearest = []
    pair_values = []
    pair_counts = Counter({cut: 0 for cut in COORDINATION_CUTOFFS_A})
    ti_count = 0
    c_count = 0
    ti_hist = {cut: Counter() for cut in COORDINATION_CUTOFFS_A}
    c_hist = {cut: Counter() for cut in COORDINATION_CUTOFFS_A}
    for atoms in atoms_list:
        symbols = np.asarray(atoms.get_chemical_symbols())
        dmat = atoms.get_all_distances(mic=True)
        tis = np.where(symbols == "Ti")[0]
        carbons = np.where(symbols == "C")[0]
        ti_count += len(tis)
        c_count += len(carbons)
        if not len(tis) or not len(carbons):
            continue
        block = dmat[np.ix_(tis, carbons)]
        ti_nearest.extend(np.min(block, axis=1).tolist())
        c_nearest.extend(np.min(block, axis=0).tolist())
        pair_values.extend(block[(block > 1.0e-8) & (block <= PAIR_DISTANCE_CUTOFF_A)].tolist())
        ii, jj, dd = neighbor_list("ijd", atoms, cutoff=max(COORDINATION_CUTOFFS_A), self_interaction=False)
        ti_per_atom = {int(idx): {cut: 0 for cut in COORDINATION_CUTOFFS_A} for idx in tis}
        c_per_atom = {int(idx): {cut: 0 for cut in COORDINATION_CUTOFFS_A} for idx in carbons}
        for i, j, distance in zip(ii, jj, dd):
            if symbols[i] != "Ti" or symbols[j] != "C":
                continue
            for cut in COORDINATION_CUTOFFS_A:
                if distance <= cut:
                    ti_per_atom[int(i)][cut] += 1
                    c_per_atom[int(j)][cut] += 1
                    pair_counts[cut] += 1
        for cut in COORDINATION_CUTOFFS_A:
            ti_hist[cut].update(row[cut] for row in ti_per_atom.values())
            c_hist[cut].update(row[cut] for row in c_per_atom.values())
    return {
        "frame_count": len(atoms_list),
        "ti_count": ti_count,
        "c_count": c_count,
        "ti_nearest_C_A": quantiles(ti_nearest),
        "c_nearest_Ti_A": quantiles(c_nearest),
        "min_image_pair_distances_Ti_C_le4A": quantiles(pair_values),
        "first_shell_coordination": {
            str(cut): {
                "cutoff_A": cut,
                "Ti_neighbor_count_histogram": {str(k): int(v) for k, v in sorted(ti_hist[cut].items())},
                "C_neighbor_count_histogram": {str(k): int(v) for k, v in sorted(c_hist[cut].items())},
                "periodic_Ti_C_pair_count": int(pair_counts[cut]),
            }
            for cut in COORDINATION_CUTOFFS_A
        },
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n")
        return
    columns = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else value for key, value in row.items()})


def main() -> None:
    inventory = build_run()
    output_dir = Path(__file__).resolve().parent
    (output_dir / "inventory.json").write_text(json.dumps(inventory, indent=2, ensure_ascii=False) + "\n")
    write_csv(output_dir / "frame_inventory.csv", inventory["frame_inventory"])
    write_csv(output_dir / "training_provenance.csv", inventory["training_provenance"])
    write_csv(output_dir / "coordination_summary.csv", inventory["coordination_csv_rows"])
    write_csv(output_dir / "force_distribution.csv", inventory["force_csv_rows"])
    summary = inventory["dataset_summary"]
    print(
        json.dumps(
            {
                "status": inventory["status"],
                "train_frames": summary["frame_count"],
                "formulas": summary["formula_counts"],
                "pure_Ti3SiC2_frames": summary["exact_pure_Ti3SiC2_formula_frame_count"],
                "interface_families": summary["frame_counts_by_interface_family"],
                "coordination_2p5A": inventory["reference_geometry_comparison"]["training_all_ti_c_coordination"]["first_shell_coordination"][str(TIC_FIRST_SHELL_CUTOFF_A)],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
