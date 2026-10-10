#!/usr/bin/env python3
"""Reproduce the allowed frozen-model diagnostics; never trains or runs DFT/MD."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from ase.geometry import find_mic
from ase.io import read
from mace.calculators import MACECalculator


ROOT = Path(os.environ.get("CLOUD_C_REPO_ROOT", Path(__file__).resolve().parents[6])).resolve()
BASE = ROOT / "research/mace-v12-transfer"
OUT = Path(__file__).resolve().parent
AUDIT = BASE / "coordination/reports/window-c/v18_unified_potential_audit_20261010"
PURE = BASE / "periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations"
RESIDUAL = BASE / "periodic_interface_v4/pbe_interface_v19_residual_acquisition/calculations"
GAP = BASE / "periodic_interface_v4/pbe_interface_v19_gap_scan/calculations"

FOUNDATION_PATH = BASE / "periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
V18_PATH = BASE / "periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
FOUNDATION_SHA = "2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012"
V18_SHA = "757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62"

REUSED_LABELS = {
    "Ag_baseline_k8x8x8_v19": "solid_Ag_baseline",
    "Ti3SiC2_baseline_k10x10x2_v19": "solid_Ti3SiC2_baseline_kz2",
    "Ti3SiC2_baseline_k10x10x4_v19": "solid_Ti3SiC2_baseline_kz4",
    "Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19": "Ti3SiC2_Ti_response_minus",
    "Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19": "Ti3SiC2_Ti_response_plus",
    "Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19": "Ti3SiC2_C_response_minus_only",
    "AgSi_COD9009647_pilot_k6x6": "AgSi_COD_slab_original_vacuum",
    "AgSi_COD9009647_vacuum26_k6x6_sigma0p10": "AgSi_COD_slab_vacuum26",
}

AG_GROUPS = {
    "displacement_k6": [
        "Ag_baseline_k6x6x6_v19",
        "Ag_displacement_m01_k6x6x6_v19",
        "Ag_displacement_p01_k6x6x6_v19",
    ],
    "shear_k8": [
        "Ag_baseline_k8x8x8_v19",
        "Ag_shear_m005_k8x8x8_v19",
        "Ag_shear_p005_k8x8x8_v19",
    ],
    "volume_k8": [
        "Ag_baseline_k8x8x8_v19",
        "Ag_volume_m02_k8x8x8_v19",
        "Ag_volume_p02_k8x8x8_v19",
    ],
}
AG_LABELS = sorted({label for labels in AG_GROUPS.values() for label in labels})
INTERFACE_LABELS = [
    "AgC_residual_mode_m_06_v19_acq",
    "AgC_residual_mode_p_06_v19_acq",
    "AgSi_residual_mode_m_03_v19_acq",
    "AgSi_residual_mode_p_03_v19_acq",
    "AgTi_residual_mode_m_03_v19_acq",
    "AgTi_residual_mode_p_03_v19_acq",
]
NEW_SLAB_LABEL = "AgSi_rigid_Ag_z_m020_k6x6x1_sigma0p10_dipoleXY_v19"

torch.set_num_threads(4)
torch.set_num_interop_threads(1)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def verify_sums(directory: Path) -> list[dict[str, str]]:
    manifest = directory / "SHA256SUMS.txt"
    if not manifest.is_file():
        raise RuntimeError(f"No SHA256SUMS.txt in {directory.name}")
    checked = []
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        item = (directory / name).resolve()
        if directory.resolve() not in item.parents or not item.is_file():
            raise RuntimeError(f"Invalid or missing checksum entry: {directory.name}/{name}")
        actual = sha(item)
        if actual != expected:
            raise RuntimeError(f"Checksum mismatch: {directory.name}/{name}")
        checked.append({"file": name, "sha256": actual})
    if not checked:
        raise RuntimeError(f"Empty checksum manifest: {directory.name}")
    return checked


def safe_label(label: str) -> None:
    lowered = label.casefold()
    forbidden = ("sealed", "blind", "holdout", "withheld", "test")
    if any(word in lowered for word in forbidden):
        raise RuntimeError(f"Refusing non-public or sealed label: {label}")


def load_verified_archive(label: str, folder: Path, expected_source: str | None = None,
                          expected_result: str | None = None) -> tuple[object, np.ndarray, dict]:
    safe_label(label)
    checksums = verify_sums(folder)
    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"SCF/label check failed: {label}")
    if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
        raise RuntimeError(f"Archive verifier failed: {label}")
    result = folder / f"{label}_PW_PBE.extxyz"
    result_hash = sha(result)
    if verification.get("sha256", {}).get(result.name) != result_hash:
        raise RuntimeError(f"Result hash disagrees with verification: {label}")
    if expected_source and summary.get("source_sha256") != expected_source:
        raise RuntimeError(f"Source geometry hash disagrees with public ledger: {label}")
    if expected_result and result_hash != expected_result:
        raise RuntimeError(f"Result hash disagrees with published audit: {label}")
    atoms = read(result, index=0)
    if "PW_PBE_forces" not in atoms.arrays:
        raise RuntimeError(f"No published DFT forces in archive: {label}")
    forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    if forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
        raise RuntimeError(f"Invalid DFT force array: {label}")
    if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all():
        raise RuntimeError(f"Invalid structure geometry: {label}")
    record = {
        "label": label,
        "summary": summary,
        "verification": verification,
        "source_sha256": summary.get("source_sha256", ""),
        "result_sha256": result_hash,
        "summary_sha256": sha(folder / "summary.json"),
        "dataset_role": summary.get("dataset_role", ""),
        "scf_converged": True,
        "verification_status": "PASS",
        "checksums": checksums,
        "result_path": str(result.relative_to(ROOT)),
    }
    return atoms, forces, record


def summarize(label: str, category: str, family: str, model_name: str,
              symbols: np.ndarray, reference: np.ndarray, predicted: np.ndarray,
              origin: str, role: str) -> tuple[dict, list[dict], list[dict]]:
    error = predicted - reference
    norms = np.linalg.norm(error, axis=1)
    metrics = {
        "label": label,
        "category": category,
        "parent_family_id": family,
        "model": model_name,
        "prediction_origin": origin,
        "dataset_role": role,
        "n_atoms": int(len(symbols)),
        "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(norms**2))),
        "x_RMSE_eV_A": float(np.sqrt(np.mean(error[:, 0]**2))),
        "y_RMSE_eV_A": float(np.sqrt(np.mean(error[:, 1]**2))),
        "z_RMSE_eV_A": float(np.sqrt(np.mean(error[:, 2]**2))),
        "max_atom_force_vector_error_eV_A": float(norms.max()),
        "z_squared_error_fraction": float(np.sum(error[:, 2]**2) / max(float(np.sum(error**2)), 1e-30)),
    }
    element_rows = []
    atom_rows = []
    for element in sorted(set(symbols.tolist())):
        mask = symbols == element
        el_norms = norms[mask]
        element_rows.append({
            "label": label,
            "category": category,
            "parent_family_id": family,
            "model": model_name,
            "element": element,
            "n_atoms": int(mask.sum()),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(el_norms**2))),
            "x_RMSE_eV_A": float(np.sqrt(np.mean(error[mask, 0]**2))),
            "y_RMSE_eV_A": float(np.sqrt(np.mean(error[mask, 1]**2))),
            "z_RMSE_eV_A": float(np.sqrt(np.mean(error[mask, 2]**2))),
            "max_atom_force_vector_error_eV_A": float(el_norms.max()),
        })
    return metrics, element_rows, atom_rows


def add_atom_rows(label: str, category: str, family: str, model_name: str,
                  symbols: np.ndarray, atom_ids: np.ndarray, positions: np.ndarray,
                  reference: np.ndarray, predicted: np.ndarray, origin: str) -> list[dict]:
    error = predicted - reference
    rows = []
    for i in range(len(symbols)):
        row = {
            "label": label,
            "category": category,
            "parent_family_id": family,
            "model": model_name,
            "prediction_origin": origin,
            "index_zero_based": i,
            "atom_id": int(atom_ids[i]),
            "element": str(symbols[i]),
            "x_A": float(positions[i, 0]) if np.isfinite(positions[i, 0]) else "",
            "y_A": float(positions[i, 1]) if np.isfinite(positions[i, 1]) else "",
            "z_A": float(positions[i, 2]) if np.isfinite(positions[i, 2]) else "",
            "DFT_Fx_eV_A": float(reference[i, 0]),
            "DFT_Fy_eV_A": float(reference[i, 1]),
            "DFT_Fz_eV_A": float(reference[i, 2]),
            "model_Fx_eV_A": float(predicted[i, 0]),
            "model_Fy_eV_A": float(predicted[i, 1]),
            "model_Fz_eV_A": float(predicted[i, 2]),
            "error_Fx_eV_A": float(error[i, 0]),
            "error_Fy_eV_A": float(error[i, 1]),
            "error_Fz_eV_A": float(error[i, 2]),
            "error_norm_eV_A": float(np.linalg.norm(error[i])),
        }
        rows.append(row)
    return rows


def csv_write(name: str, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {name}")
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def public_repo_path(value: str) -> str:
    """Reduce archived path metadata to a repository-relative path, never publish host paths."""
    parts = Path(value).parts
    try:
        index = parts.index("research")
    except ValueError:
        return "see published DFT coverage ledger"
    return Path(*parts[index:]).as_posix()


def model_pair_force(atoms, forces: np.ndarray) -> tuple[float, str]:
    markers = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    indices = np.flatnonzero(markers)
    if len(indices) != 2:
        raise RuntimeError(f"Expected exactly two central_pair atoms; saw {len(indices)}")
    ag = [int(i) for i in indices if atoms[int(i)].symbol == "Ag"]
    x = [int(i) for i in indices if atoms[int(i)].symbol != "Ag"]
    if len(ag) != 1 or len(x) != 1:
        raise RuntimeError("Marked central_pair must contain exactly one Ag and one contact atom")
    i, j = ag[0], x[0]
    vector, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    unit = vector / np.linalg.norm(vector)
    return float(np.dot(forces[j] - forces[i], unit)), f"Ag_id={atoms.arrays['lammps_id'][i]};X_id={atoms.arrays['lammps_id'][j]};X={atoms[j].symbol}"


def write_outputs() -> None:
    audit_checks = verify_sums(AUDIT)
    audit_by_file = {x["file"]: x["sha256"] for x in audit_checks}
    for fname in ("per_atom.csv", "per_element.csv", "ag_response.csv", "dft_coverage.csv", "same_geometry_old_potential_comparison.csv"):
        if fname not in audit_by_file:
            raise RuntimeError(f"Required published audit file missing from hash inventory: {fname}")

    pure_eligibility = read_csv(BASE / "coordination/reports/window-c/v19_controlled_repair_design_20261010/pure_phase_source_eligibility.csv")
    pure_records = {row["label"]: row for row in pure_eligibility}
    if not pure_records or any(row["train_eligible_now"] != "False" for row in pure_records.values()):
        raise RuntimeError("Pure-phase records no longer match the published control-only role review")

    dft_coverage_rows = read_csv(AUDIT / "dft_coverage.csv")
    dft_coverage = {row["label"]: row for row in dft_coverage_rows}
    interface_audit = read_csv(AUDIT / "same_geometry_old_potential_comparison.csv")
    interface_by_label = defaultdict(list)
    for row in interface_audit:
        interface_by_label[row["label"]].append(row)
    ag_prior = {row["label"]: row for row in read_csv(AUDIT / "ag_response.csv")}

    records: dict[str, dict] = {}
    structures: dict[str, object] = {}
    references: dict[str, np.ndarray] = {}
    summaries: dict[str, dict] = {}
    archive_rows = []

    def load_input(label: str, category: str, family: str, folder: Path,
                   expected_source: str | None = None, expected_result: str | None = None,
                   role: str | None = None, mesh: str = "") -> None:
        atoms, forces, rec = load_verified_archive(label, folder, expected_source, expected_result)
        if role is not None and role not in rec["dataset_role"]:
            raise RuntimeError(f"Data role mismatch for {label}: {rec['dataset_role']}")
        rec.update({"category": category, "parent_family_id": family, "mesh": mesh,
                    "atom_count": len(atoms), "formula": atoms.get_chemical_formula(),
                    "pbc": atoms.pbc.tolist(), "cell_A": atoms.cell.array.tolist()})
        records[label] = rec
        structures[label] = atoms
        references[label] = forces
        summaries[label] = rec["summary"]

    # All seven unique Ag reference geometries used by the published same-mesh response panel.
    ag_meshes = {
        "Ag_baseline_k6x6x6_v19": "6x6x6",
        "Ag_displacement_m01_k6x6x6_v19": "6x6x6",
        "Ag_displacement_p01_k6x6x6_v19": "6x6x6",
        "Ag_baseline_k8x8x8_v19": "8x8x8",
        "Ag_shear_m005_k8x8x8_v19": "8x8x8",
        "Ag_shear_p005_k8x8x8_v19": "8x8x8",
        "Ag_volume_m02_k8x8x8_v19": "8x8x8",
        "Ag_volume_p02_k8x8x8_v19": "8x8x8",
    }
    for label, mesh in ag_meshes.items():
        item = pure_records.get(label)
        if not item or item["phase"] != "solid_Ag" or item["verification_status"] != "PASS" or item["scf_converged"] != "True":
            raise RuntimeError(f"Ag control not fully verified: {label}")
        if item["dataset_role"] != "numerical_pure_phase_control_not_training_or_independent_validation":
            raise RuntimeError(f"Unexpected Ag data role: {label}")
        load_input(label, "solid_Ag", "Ag_fcc_bulk_parent_v19", PURE / label,
                   item["source_sha256"], item["result_sha256"], "numerical_pure_phase_control", mesh)

    # Read only the already completed, verified Ti3SiC2 bulk controls needed for geometry-identity checks.
    tsc_labels = [
        "Ti3SiC2_baseline_k10x10x2_v19",
        "Ti3SiC2_baseline_k10x10x4_v19",
        "Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19",
        "Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19",
        "Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19",
    ]
    for label in tsc_labels:
        item = pure_records.get(label)
        if not item or item["phase"] != "solid_Ti3SiC2" or item["verification_status"] != "PASS" or item["scf_converged"] != "True":
            raise RuntimeError(f"Ti3SiC2 control is not fully verified: {label}")
        load_input(label, "solid_Ti3SiC2", "Ti3SiC2_bulk_parent_v19", PURE / label,
                   item["source_sha256"], item["result_sha256"], "numerical_pure_phase_control", item["kpts"])

    for label in INTERFACE_LABELS:
        rows = interface_by_label.get(label, [])
        if not rows or {row["model"] for row in rows} != {"foundation", "V18"}:
            raise RuntimeError(f"Interface reference is absent from verified audit: {label}")
        if any(row["independent_validation"] != "False" for row in rows):
            raise RuntimeError(f"Interface parent-family must remain development diagnostic: {label}")
        item = dft_coverage.get(label)
        if not item or item["scf_converged"] != "True" or item["verification_PASS"] != "True" or item["sha256_all_pass"] != "True":
            raise RuntimeError(f"Interface DFT archive does not pass the published verification ledger: {label}")
        atom = "AgC" if "AgC_" in label else "AgSi" if "AgSi_" in label else "AgTi"
        source_row = rows[0]
        folder = RESIDUAL / label
        load_input(label, f"interface_{atom}", f"{atom}_residual_acquisition_parent_v19", folder,
                   source_row["DFT_source_sha256"], source_row["DFT_result_sha256"], "training_acquisition_reused_development_family")

    load_input(NEW_SLAB_LABEL, "AgSi_slab_rigid_Ag_minus020", "AgSi_rigid_gap_parent_A_v19",
               GAP / NEW_SLAB_LABEL)
    slab_verification = records[NEW_SLAB_LABEL]["verification"]
    if slab_verification.get("owner_task") != "window-a":
        raise RuntimeError("The new AgSi endpoint is not owned/published by A")
    if "diagnostic only" not in records[NEW_SLAB_LABEL]["dataset_role"]:
        raise RuntimeError("The new AgSi endpoint role is not a public diagnostic")

    # Reuse all previously published two-model per-atom predictions for the exact audited structures.
    old_atom_rows = read_csv(AUDIT / "per_atom.csv")
    for row in old_atom_rows:
        if row["label"] not in REUSED_LABELS:
            continue
        safe_label(row["label"])
    old_groups = defaultdict(list)
    for row in old_atom_rows:
        if row["label"] in REUSED_LABELS:
            old_groups[(row["label"], row["model"])].append(row)
    expected_model_keys = {(label, model) for label in REUSED_LABELS for model in ("foundation", "V18")}
    if set(old_groups) != expected_model_keys:
        raise RuntimeError("Published two-model atom predictions are incomplete for the selected audit structures")

    models = {
        "MACE-MP-0b3-medium": {"path": FOUNDATION_PATH, "sha256": FOUNDATION_SHA, "source": "ACEsuit/mace-foundations release mace_mp_0b3; repository NOTICE.md", "license": "MIT (official MACE 0.3.16 catalog)", "version": "0b3 medium", "prediction_origin": "new frozen CPU float64 inference where listed"},
        "V18": {"path": V18_PATH, "sha256": V18_SHA, "source": "project-frozen V18 epoch79 derivative; selected model record", "license": "No standalone model-weight license was found in the archived project record", "version": "frozen V18 / MACE-Torch 0.3.16 runtime", "prediction_origin": "reused published inference or new frozen CPU float64 inference"},
    }
    calcs = {}
    model_info = []
    for name, meta in models.items():
        model_path = meta["path"]
        if sha(model_path) != meta["sha256"]:
            raise RuntimeError(f"Model weight hash mismatch: {name}")
        calc = MACECalculator(model_paths=str(model_path), device="cpu", default_dtype="float64")
        model = calc.models[0]
        z_table = [int(z) for z in model.atomic_numbers.detach().cpu().tolist()]
        required = {6, 14, 22, 47}
        if not required.issubset(z_table):
            raise RuntimeError(f"Model element mapping lacks Ag/Ti/Si/C: {name}")
        calcs[name] = calc
        model_info.append({
            "model": name, "version": meta["version"], "source": meta["source"], "license": meta["license"],
            "repository_relative_weight_path": str(model_path.relative_to(ROOT)), "sha256": meta["sha256"],
            "atomic_numbers_order": json.dumps(z_table), "Ag_Ti_Si_C_supported": True,
            "r_max_A": float(model.r_max), "num_interactions": int(model.num_interactions),
            "inference_runtime": "Python 3.12.14; PyTorch 2.5.1+cpu; MACE-Torch 0.3.16; ASE 3.29.0; float64; CPU",
            "inference_compatible_here": True,
        })

    predictions: dict[tuple[str, str], np.ndarray] = {}
    energies: dict[tuple[str, str], float] = {}
    metric_rows: list[dict] = []
    element_rows: list[dict] = []
    atom_rows: list[dict] = []

    def append_metrics(label: str, model_name: str, symbols: np.ndarray, positions: np.ndarray,
                       atom_ids: np.ndarray, reference: np.ndarray, predicted: np.ndarray,
                       category: str, family: str, origin: str, role: str) -> None:
        m, els, _ = summarize(label, category, family, model_name, symbols, reference, predicted, origin, role)
        metric_rows.append(m)
        element_rows.extend(els)
        atom_rows.extend(add_atom_rows(label, category, family, model_name, symbols, atom_ids, positions,
                                       reference, predicted, origin))

    # Values already in the audit are reused without recalculating the models.
    for (label, stored_model), rows in old_groups.items():
        rows.sort(key=lambda row: int(row["index_zero_based"]))
        if stored_model not in ("foundation", "V18"):
            continue
        model_name = "MACE-MP-0b3-medium" if stored_model == "foundation" else "V18"
        symbols = np.asarray([row["element"] for row in rows])
        ref = np.asarray([[float(row[f"DFT_F{axis}_eV_A"]) for axis in "xyz"] for row in rows])
        pred = np.asarray([[float(row[f"model_F{axis}_eV_A"]) for axis in "xyz"] for row in rows])
        atom_ids = np.asarray([int(row["atom_id"]) for row in rows])
        if label in references:
            if not np.array_equal(atom_ids, np.asarray(structures[label].arrays.get("lammps_id", atom_ids), dtype=int)):
                raise RuntimeError(f"Audit/per-archive atom ID ordering mismatch: {label}")
            if not np.allclose(ref, references[label], rtol=0.0, atol=1e-12):
                raise RuntimeError(f"Audit/per-archive reference force mismatch: {label}")
            if not np.array_equal(symbols, np.asarray(structures[label].get_chemical_symbols())):
                raise RuntimeError(f"Audit/per-archive element order mismatch: {label}")
        positions = np.asarray([[np.nan, np.nan, float(row["z_A"])] for row in rows])
        category = REUSED_LABELS[label]
        family = ("Ti3SiC2_bulk_parent_v19" if label.startswith("Ti3SiC2_") else
                  "AgSi_COD9009647_slab_family" if label.startswith("AgSi_COD") else
                  "Ag_fcc_bulk_parent_v19")
        role = dft_coverage.get(label, {}).get("role", "published public diagnostic")
        append_metrics(label, model_name, symbols, positions, atom_ids, ref, pred, category, family,
                       "reused_v18_audit_per_atom.csv", role)
        predictions[(label, model_name)] = pred
        predictions[(label, "MACE-MP-0b3-medium" if stored_model == "foundation" else "V18")] = pred

    # New inference is limited to fully verified solid-Ag controls, six public local interfaces, and A's completed minus endpoint.
    new_labels = AG_LABELS + INTERFACE_LABELS + [NEW_SLAB_LABEL]
    for label in new_labels:
        atoms = structures[label]
        symbols = np.asarray(atoms.get_chemical_symbols())
        atom_ids = np.asarray(atoms.arrays.get("lammps_id", np.arange(1, len(atoms) + 1)), dtype=int)
        category = records[label]["category"]
        family = records[label]["parent_family_id"]
        role = records[label]["dataset_role"]
        for model_name, calc in calcs.items():
            work = atoms.copy()
            work.calc = calc
            predicted = np.asarray(work.get_forces(), dtype=float)
            energy = float(work.get_potential_energy())
            if predicted.shape != references[label].shape or not np.isfinite(predicted).all() or not math.isfinite(energy):
                raise RuntimeError(f"Nonfinite or malformed frozen inference: {label}/{model_name}")
            if label in REUSED_LABELS:
                prior = predictions[(label, model_name)]
                difference = float(np.max(np.abs(predicted - prior)))
                if difference > 1e-8:
                    raise RuntimeError(f"Reused published inference disagrees with this deterministic check: {label}/{model_name}")
            else:
                predictions[(label, model_name)] = predicted.copy()
            energies[(label, model_name)] = energy
            if label not in REUSED_LABELS:
                append_metrics(label, model_name, symbols, atoms.positions.copy(), atom_ids,
                               references[label], predicted, category, family,
                               "new_frozen_CPU_float64_inference", role)

    # Reproduction checks against the previously published same-input aggregate scores.
    metric_lookup = {(row["label"], row["model"]): row for row in metric_rows}
    repro = []
    for label in AG_LABELS:
        old = ag_prior[label]
        for model_name, old_prefix in (("MACE-MP-0b3-medium", "foundation"), ("V18", "V18")):
            actual = metric_lookup[(label, model_name)]["force_vector_RMSE_eV_A"]
            expected = float(old[f"{old_prefix}_force_RMSE_eV_A"])
            if abs(actual - expected) > 2e-7:
                raise RuntimeError(f"Ag response RMSE does not reproduce the published audit: {label}/{model_name}")
            repro.append({"label": label, "model": model_name, "metric": "force_vector_RMSE_eV_A", "prior": expected, "recomputed": actual, "abs_difference": abs(actual - expected), "status": "PASS"})
    for label in INTERFACE_LABELS:
        for prior_row in interface_by_label[label]:
            model_name = "MACE-MP-0b3-medium" if prior_row["model"] == "foundation" else "V18"
            actual = metric_lookup[(label, model_name)]["force_vector_RMSE_eV_A"]
            expected = float(prior_row["model_force_vector_RMSE_eV_A"])
            if abs(actual - expected) > 2e-8:
                raise RuntimeError(f"Interface RMSE does not reproduce the published audit: {label}/{model_name}")
            repro.append({"label": label, "model": model_name, "metric": "force_vector_RMSE_eV_A", "prior": expected, "recomputed": actual, "abs_difference": abs(actual - expected), "status": "PASS"})

    # Relative energy/force response on matched-mesh Ag parent families; never compare absolute model zeros.
    ag_response_rows = []
    for group, labels in AG_GROUPS.items():
        baseline = labels[0]
        ref_base = summaries[baseline]
        for label in labels:
            summary = summaries[label]
            atoms = structures[label]
            record = {
                "response_family": group,
                "parent_family_id": "Ag_fcc_bulk_parent_v19",
                "label": label,
                "kpts": json.dumps(summary["method"]["kpts"]),
                "n_atoms": len(atoms),
                "DFT_relative_native_meV_atom": (float(summary["energy_eV_cell"]) - float(ref_base["energy_eV_cell"])) * 1000.0 / len(atoms),
                "DFT_relative_force_consistent_meV_atom": (float(summary["free_energy_eV_cell"]) - float(ref_base["free_energy_eV_cell"])) * 1000.0 / len(atoms),
                "dataset_role": records[label]["dataset_role"],
            }
            for model_name in calcs:
                row = metric_lookup[(label, model_name)]
                base_e = energies[(baseline, model_name)]
                ag_response_rows.append({**record,
                    "model": model_name,
                    "model_relative_energy_meV_atom": (energies[(label, model_name)] - base_e) * 1000.0 / len(atoms),
                    "force_vector_RMSE_eV_A": row["force_vector_RMSE_eV_A"],
                    "max_atom_force_vector_error_eV_A": row["max_atom_force_vector_error_eV_A"],
                })
    # Validate each selected Ag relative energy and force against the previous published table.
    for row in ag_response_rows:
        old = ag_prior[row["label"]]
        key = "foundation" if row["model"] == "MACE-MP-0b3-medium" else "V18"
        expected_e = float(old[f"{key}_relative_energy_meV_atom"])
        if abs(row["model_relative_energy_meV_atom"] - expected_e) > 2e-3:
            raise RuntimeError(f"Ag relative energy does not reproduce published audit: {row['label']}/{row['model']}")

    response_rows = []
    # The symmetric one-atom Ag displacement response.
    bm, bp = AG_GROUPS["displacement_k6"][1:]
    am, ap = structures[bm], structures[bp]
    changed = np.flatnonzero(np.linalg.norm(ap.positions - am.positions, axis=1) > 1e-9)
    if len(changed) != 1:
        raise RuntimeError("Ag displacement pair must change exactly one atom")
    i = int(changed[0])
    displacement = ap.positions[i] - am.positions[i]
    span = float(np.linalg.norm(displacement))
    unit = displacement / span
    response = {"response": "Ag_single_atom_displacement", "parent_family_id": "Ag_fcc_bulk_parent_v19",
                "status": "complete_symmetric_pair", "atom_index_zero_based": i,
                "atom_id": int(ap.arrays.get("lammps_id", np.arange(1, len(ap)+1))[i]),
                "element": "Ag", "direction_xyz": json.dumps(unit.tolist()), "span_A": span}
    response["DFT_restoring_slope_eV_A2"] = float(-np.dot(references[bp][i] - references[bm][i], unit) / span)
    for model_name in calcs:
        response[f"{model_name}_restoring_slope_eV_A2"] = float(-np.dot(predictions[(bp, model_name)][i] - predictions[(bm, model_name)][i], unit) / span)
    response_rows.append(response)

    # Matched-k Ag response derivatives use only relative energies within each parent family.
    def deformation(parent, target):
        return target.cell.array.T @ np.linalg.inv(parent.cell.array.T)

    for group, measure, param_name in (("shear_k8", "shear", "engineering_shear"), ("volume_k8", "volume", "log_volume")):
        labels = AG_GROUPS[group]
        parent = structures[labels[0]]
        minus, plus = structures[labels[1]], structures[labels[2]]
        if measure == "shear":
            fm, fp = deformation(parent, minus), deformation(parent, plus)
            offdiag = [(a, b) for a in range(3) for b in range(3) if a != b]
            pair = max(offdiag, key=lambda ij: max(abs(fm[ij]), abs(fp[ij])))
            pm, pp = float(fm[pair]), float(fp[pair])
            if abs(pm + pp) > 2e-7 or abs(pp - pm) < 1e-8:
                raise RuntimeError("Ag shear probes are not a symmetric matched strain pair")
            delta = pp - pm
        else:
            vm = float(minus.get_volume() / parent.get_volume())
            vp = float(plus.get_volume() / parent.get_volume())
            pm, pp = math.log(vm), math.log(vp)
            if not np.isclose(vm, 0.98, atol=2e-6) or not np.isclose(vp, 1.02, atol=2e-6):
                raise RuntimeError("Ag volume probes are not the published -2%/+2% volume bracket")
            delta = pp - pm
        row = {"response": f"Ag_{measure}", "parent_family_id": "Ag_fcc_bulk_parent_v19",
               "status": "complete_symmetric_pair" if measure == "shear" else "complete_paired_bracket_secant",
               "parameter": param_name,
               "parameter_minus": pm, "parameter_plus": pp, "parameter_span": delta,
               "DFT_native_dE_dparameter_eV_cell": (float(summaries[labels[2]]["energy_eV_cell"]) - float(summaries[labels[1]]["energy_eV_cell"])) / delta,
               "DFT_force_consistent_dE_dparameter_eV_cell": (float(summaries[labels[2]]["free_energy_eV_cell"]) - float(summaries[labels[1]]["free_energy_eV_cell"])) / delta}
        if measure == "shear":
            row["deformation_gradient_offdiagonal_index"] = f"{pair[0]},{pair[1]}"
        for model_name in calcs:
            row[f"{model_name}_dE_dparameter_eV_cell"] = (energies[(labels[2], model_name)] - energies[(labels[1], model_name)]) / delta
        response_rows.append(row)

    # Ti restoring response is complete; C+ is still running and is never opened or synthesized here.
    for element, minus_label, plus_label, atom_index, axis in (
        ("Ti", "Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19", "Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19", 2, 2),
        ("C", "Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19", "Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19", 8, 2),
    ):
        if (plus_label not in structures or minus_label not in structures or
                any((label, model_name) not in predictions
                    for label in (minus_label, plus_label) for model_name in calcs)):
            response_rows.append({"response": f"Ti3SiC2_{element}_z_restoring_slope", "parent_family_id": "Ti3SiC2_bulk_parent_v19",
                                  "status": "PAIR_INCOMPLETE_PLUS_NOT_ARCHIVED", "minus_label": minus_label,
                                  "plus_label": plus_label, "scope": "minus-only public pointwise force errors remain in per-atom/element CSV"})
            continue
        minus_atoms = structures[minus_label]
        plus_atoms = structures[plus_label]
        changed_indices = np.flatnonzero(np.linalg.norm(plus_atoms.positions - minus_atoms.positions, axis=1) > 1e-9)
        if len(changed_indices) != 1 or int(changed_indices[0]) != atom_index:
            raise RuntimeError(f"Expected only the target {element} atom to move in the symmetric pair")
        if not np.array_equal(plus_atoms.arrays.get("lammps_id"), minus_atoms.arrays.get("lammps_id")):
            raise RuntimeError(f"Atom identity/order changed in the {element} response pair")
        actual_span = float(plus_atoms.positions[atom_index, axis] - minus_atoms.positions[atom_index, axis])
        if not np.isclose(abs(actual_span), 0.04, atol=1e-8):
            raise RuntimeError(f"Unexpected Ti/C displacement span for {element}")
        row = {"response": f"Ti3SiC2_{element}_z_restoring_slope", "parent_family_id": "Ti3SiC2_bulk_parent_v19",
               "status": "complete_symmetric_pair", "minus_label": minus_label, "plus_label": plus_label,
               "atom_index_zero_based": atom_index, "atom_id": int(old_groups[(minus_label, "foundation")][0]["atom_id"]),
               "span_A": actual_span, "DFT_restoring_slope_eV_A2": float(-(references[plus_label][atom_index, axis] - references[minus_label][atom_index, axis]) / actual_span)}
        row["atom_id"] = int(next(r["atom_id"] for r in old_groups[(minus_label, "foundation")] if int(r["index_zero_based"]) == atom_index))
        for model_name in calcs:
            row[f"{model_name}_restoring_slope_eV_A2"] = float(-(predictions[(plus_label, model_name)][atom_index, axis] - predictions[(minus_label, model_name)][atom_index, axis]) / actual_span)
        response_rows.append(row)

    # Contact separation force and Ag-layer generalized force for local references and the completed A endpoint.
    interface_rows = []
    for label in INTERFACE_LABELS + [NEW_SLAB_LABEL]:
        atoms = structures[label]
        ag_mask = np.asarray(atoms.get_chemical_symbols()) == "Ag"
        ag_z = atoms.positions[ag_mask, 2]
        planes = sorted(set(np.round(ag_z, 3).tolist()))
        rec = records[label]
        contact = "Ag-C" if "AgC" in label else "Ag-Si" if "AgSi" in label else "Ag-Ti"
        for model_name in calcs:
            model_force, pair_description = model_pair_force(atoms, predictions[(label, model_name)])
            dft_force, _ = model_pair_force(atoms, references[label])
            model_gz = float(predictions[(label, model_name)][ag_mask, 2].sum())
            dft_gz = float(references[label][ag_mask, 2].sum())
            interface_rows.append({
                "label": label, "category": rec["category"], "contact": contact,
                "parent_family_id": rec["parent_family_id"], "model": model_name,
                "central_pair_ids": pair_description,
                "DFT_marked_pair_separating_force_eV_A": dft_force,
                "model_marked_pair_separating_force_eV_A": model_force,
                "marked_pair_signed_error_eV_A": model_force - dft_force,
                "marked_pair_abs_error_eV_A": abs(model_force - dft_force),
                "DFT_all_Ag_generalized_Fz_eV_A": dft_gz,
                "model_all_Ag_generalized_Fz_eV_A": model_gz,
                "all_Ag_generalized_Fz_error_eV_A": model_gz - dft_gz,
                "Ag_z_plane_count_rounded_0p001A": len(planes),
                "Ag_single_layer_mask": len(planes) == 1,
                "paired_Ag_z_endpoint_complete": label != NEW_SLAB_LABEL,
                "dataset_role": rec["dataset_role"],
            })

    # Three-point curvature for Ag shear/volume response. Energies are recentered
    # within each model/reference family so cross-model absolute zeros are never used.
    elastic_rows = []
    for group, measure, param_name in (("shear_k8", "shear", "engineering_shear"),
                                       ("volume_k8", "volume", "log_volume")):
        labels = AG_GROUPS[group]
        parent, minus, plus = (structures[x] for x in labels)
        if measure == "shear":
            pm = float(deformation(parent, minus)[0, 1])
            pp = float(deformation(parent, plus)[0, 1])
        else:
            pm = math.log(minus.get_volume() / parent.get_volume())
            pp = math.log(plus.get_volume() / parent.get_volume())
        parameters = np.asarray([pm, 0.0, pp], dtype=float)
        if abs(parameters[0]) < 1e-10 or abs(parameters[2]) < 1e-10 or not parameters[0] < 0 < parameters[2]:
            raise RuntimeError(f"Invalid Ag {measure} response bracket")
        row = {"response": f"Ag_{measure}_three_point_energy_curvature",
               "parent_family_id": "Ag_fcc_bulk_parent_v19", "parameter": param_name,
               "parameter_minus": pm, "parameter_parent": 0.0, "parameter_plus": pp,
               "minus_label": labels[1], "parent_label": labels[0], "plus_label": labels[2],
               "energy_zero_policy": "subtract the same-model/reference parent energy before fitting"}
        def add_fit(prefix: str, values: list[float]) -> None:
            relative = np.asarray(values, dtype=float) - float(values[1])
            quadratic, linear, _ = np.polyfit(parameters, relative, 2)
            row[f"{prefix}_dE_dparameter_at_parent_eV_cell"] = float(linear)
            row[f"{prefix}_d2E_dparameter2_at_parent_eV_cell"] = float(2.0 * quadratic)
        fit_labels = [labels[1], labels[0], labels[2]]  # minus, parent, plus
        add_fit("DFT_native", [float(summaries[x]["energy_eV_cell"]) for x in fit_labels])
        add_fit("DFT_force_consistent", [float(summaries[x]["free_energy_eV_cell"]) for x in fit_labels])
        for model_name in calcs:
            add_fit(model_name, [energies[(x, model_name)] for x in fit_labels])
        elastic_rows.append(row)

    # Copy old hybrid-potential results into their own source-labeled table; never mix with MACE rows.
    hybrid_rows = []
    for label in INTERFACE_LABELS:
        source = interface_by_label[label][0]
        hybrid_rows.append({
            "label": label, "contact": source["contact"], "role": source["role"],
            "model_source": "published_candidate02_EAM_TSCBOP_Morse_hybrid_LAMMPS_run0",
            "old_hybrid_force_vector_RMSE_eV_A": source["old_candidate02_hybrid_force_vector_RMSE_eV_A"],
            "old_hybrid_marked_pair_signed_error_eV_A": source["old_candidate02_hybrid_marked_pair_signed_error_eV_A"],
            "old_hybrid_marked_pair_abs_error_eV_A": source["old_candidate02_hybrid_abs_pair_error_eV_A"],
            "DFT_source_sha256": source["DFT_source_sha256"],
            "DFT_result_sha256": source["DFT_result_sha256"],
            "independent_validation": "False",
            "source_report": "coordination/reports/window-c/v18_unified_potential_audit_20261010/same_geometry_old_potential_comparison.csv",
            "note": "Historical LAMMPS output; not run by this task and not a MACE prediction.",
        })

    for label, rec in records.items():
        archive_rows.append({
            "label": label, "category": rec["category"], "parent_family_id": rec["parent_family_id"],
            "formula": rec["formula"], "n_atoms": rec["atom_count"], "pbc": json.dumps(rec["pbc"]),
            "dataset_role": rec["dataset_role"], "scf_converged": rec["scf_converged"],
            "archive_verification": rec["verification_status"], "source_sha256": rec["source_sha256"],
            "result_sha256": rec["result_sha256"], "summary_sha256": rec["summary_sha256"],
            "result_repo_path": rec["result_path"], "mesh": rec["mesh"],
            "independent_validation": False, "training_use": False,
        })
    # Add re-used audit references with their archived hash/SCF ledger metadata.
    for label, category in REUSED_LABELS.items():
        if label in records:
            continue
        source = dft_coverage[label]
        family = ("Ti3SiC2_bulk_parent_v19" if label.startswith("Ti3SiC2_") else
                  "AgSi_COD9009647_slab_family" if label.startswith("AgSi_COD") else
                  "Ag_fcc_bulk_parent_v19")
        archive_rows.append({
            "label": label, "category": category, "parent_family_id": family,
            "formula": source["formula"], "n_atoms": source["atoms"], "pbc": "published audit input",
            "dataset_role": source["role"], "scf_converged": source["scf_converged"],
            "archive_verification": "PASS" if source["verification_PASS"] == "True" else "FAIL",
            "source_sha256": source["source_sha256"], "result_sha256": source["result_sha256"],
            "summary_sha256": "covered by published audit manifest", "result_repo_path": public_repo_path(source["path"]),
            "mesh": "from published DFT coverage ledger", "independent_validation": False, "training_use": False,
        })

    # Model records, including MATPES's license/provenance blocker; no MATPES artifact was fetched or loaded.
    model_info.append({
        "model": "MACE-MATPES-PBE-0", "version": "official catalog identifier; artifact version/hash unavailable",
        "source": "ACEsuit/mace v0.3.16 model catalog, upstream commit 4d2da09413ac1407f37cdbb6b81fa28e4c15655e; catalog asset MACE-matpes-pbe-omat-ft.model",
        "license": "ASL in official catalog; project permission/scope not confirmed",
        "repository_relative_weight_path": "not present", "sha256": "not available; model was not downloaded",
        "atomic_numbers_order": "catalog declares 89 elements; exact artifact z_table not inspected",
        "Ag_Ti_Si_C_supported": "catalog element count only; exact map unverified",
        "r_max_A": "unknown", "num_interactions": "unknown",
        "inference_runtime": "catalog requires MACE >=0.3.10; available runtime 0.3.16 meets minimum",
        "inference_compatible_here": "API minimum passes; weight-level compatibility unverified",
    })

    # Emit machine-readable outputs.
    csv_write("structure_metrics.csv", metric_rows)
    csv_write("per_element_errors.csv", element_rows)
    csv_write("per_atom_errors.csv", atom_rows)
    csv_write("ag_response.csv", ag_response_rows)
    csv_write("backbone_interface_responses.csv", response_rows)
    csv_write("ag_elastic_response.csv", elastic_rows)
    csv_write("interface_separation_and_layer_forces.csv", interface_rows)
    csv_write("old_hybrid_vs_dft.csv", hybrid_rows)
    csv_write("reference_inventory.csv", archive_rows)
    csv_write("model_inventory.csv", model_info)
    csv_write("reproduction_checks.csv", repro)

    proof = {
        "status": "completed_with_MATPES_blocker",
        "baseline_audit_commit": "e455b71f6ef91a0955495001765e659aaee12174",
        "comparison_start_commit": "5857001",
        "runtime": {"python": "3.12.14", "pytorch": "2.5.1+cpu", "mace_torch": "0.3.16", "ase": "3.29.0", "numpy": "1.26.4", "scipy": "1.17.1", "device": "CPU", "dtype": "float64", "torch_threads": 4},
        "models": {k: {"sha256": v["sha256"], "repo_relative_path": str(v["path"].relative_to(ROOT))} for k, v in models.items()},
        "audit_source_sha256": audit_by_file,
        "reused_reference_labels": list(REUSED_LABELS),
        "scored_reference_labels": new_labels,
        "new_metric_labels": [label for label in new_labels if label not in REUSED_LABELS],
        "deterministic_recheck_labels": [label for label in new_labels if label in REUSED_LABELS],
        "MATPES_status": {"license": "ASL; no project scope confirmation", "metadata_endpoint_http": 403,
                          "model_download": "not attempted after license blocker", "inference": "not run", "bypass": False},
        "prohibited_work_started": {"DFT": False, "MACE_training": False, "MD": False, "LAMMPS": False, "TTM": False},
        "unique_new_metric_label_count": len(set(new_labels) - set(REUSED_LABELS)),
        "deterministic_recheck_label_count": len(set(new_labels) & set(REUSED_LABELS)),
        "reused_prediction_count": len(REUSED_LABELS) * 2,
        "predicted_model_names": list(models),
        "note": "All scored geometries are public development diagnostics or numerical controls, not an independent test set. No absolute cross-model energy comparison is made.",
    }
    (OUT / "provenance.json").write_text(json.dumps(proof, indent=2) + "\n")
    print(json.dumps({"structures_new_metric_labels": len(set(new_labels) - set(REUSED_LABELS)),
                      "structures_deterministically_rechecked": len(set(new_labels) & set(REUSED_LABELS)),
                      "structures_reused_from_published_audit": len(REUSED_LABELS),
                      "metrics": len(metric_rows), "elements": len(element_rows), "atoms": len(atom_rows),
                      "model_info": model_info}, indent=2))


if __name__ == "__main__":
    start = time.perf_counter()
    write_outputs()
    print(f"elapsed_s={time.perf_counter() - start:.2f}")
