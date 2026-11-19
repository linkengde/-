#!/usr/bin/env python3
"""Build an exploratory ACE fitting set from verified v24 and v11 labels.

The v24 periodic-interface labels use their native extrapolated PBE energies
and finite-smearing forces. The v11 hotspot crops contribute only their
registered center-Si force; crop energies and all edge forces are excluded by
the fitting weights written into the accompanying manifest.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "v24_c_response_run"
V11 = ROOT / "v11_hotspot_reference/archive/round7-v11-si-hotspot-labels-20261118"
BASE = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/candidate_data")
OUT = ROOT / "fit_candidate_acquisition"
SELECTED_V24 = [
    "AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24",
    "AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24",
    "AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24",
]
V11_TRAIN = [0, 2, 3]
V11_DEV = 1


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    base_manifest = json.loads((BASE / "candidate_manifest.json").read_text())
    base_train_path = BASE / "train_kconv.extxyz"
    base_dev_path = BASE / "dev.extxyz"
    if sha(base_train_path) != base_manifest["candidate_train_sha256"]:
        raise SystemExit("base training dataset hash mismatch")
    if sha(base_dev_path) != base_manifest["candidate_dev_sha256"]:
        raise SystemExit("base development dataset hash mismatch")

    train = read(base_train_path, index=":")
    dev = read(base_dev_path, index=":")
    train_roles = [{"role": "existing_training_candidate", "config_type": a.info.get("config_type")}
                   for a in train]

    v24_manifest = json.loads((RUN / "input_manifest.json").read_text())
    v24_by_label = {job["label"]: job for job in v24_manifest["preflight_jobs"]}
    for label in SELECTED_V24:
        folder = RUN / "calculations" / label
        summary_path = folder / "summary.json"
        verification_path = folder / "verification.json"
        if not summary_path.is_file() or not verification_path.is_file():
            raise SystemExit(f"verified v24 result is missing: {label}")
        summary = json.loads(summary_path.read_text())
        verification = json.loads(verification_path.read_text())
        job = v24_by_label[label]
        if verification.get("status") != "PASS" or summary.get("scf_converged") is not True:
            raise SystemExit(f"v24 result did not pass the registered SCF/archive audit: {label}")
        if summary.get("source_sha256") != job["input_sha256"]:
            raise SystemExit(f"v24 source hash mismatch: {label}")
        atoms = read(folder / f"{label}_PW_PBE.extxyz")
        if atoms.info.get("PW_PBE_energy_eV") != summary.get("native_extrapolated_energy_eV_cell"):
            raise SystemExit(f"v24 native energy mismatch: {label}")
        forces = np.asarray(atoms.arrays.get("PW_PBE_forces"), dtype=float)
        if forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise SystemExit(f"v24 force array invalid: {label}")
        atoms.info["REF_energy"] = float(atoms.info["PW_PBE_energy_eV"])
        atoms.arrays["REF_forces"] = forces.copy()
        atoms.info["dataset_role"] = "training_acquisition_candidate"
        atoms.info["training_label_policy"] = (
            "native_extrapolated_PBE_energy; all-atom forces at Fermi-Dirac 0.10 eV; "
            "periodic interface acquisition; not validation"
        )
        train.append(atoms)
        train_roles.append({"role": "v24_periodic_interface_train_acquisition",
                            "config_type": label, "input_sha256": job["input_sha256"],
                            "result_sha256": sha(folder / f"{label}_PW_PBE.extxyz"),
                            "energy_convention": summary["energy_convention"]})

    split_path = V11 / "v12_label_split.json"
    split = json.loads(split_path.read_text())
    if split.get("audit_status") != "PASS" or split.get("training_indices") != V11_TRAIN:
        raise SystemExit("v11 hotspot split audit or locked training indices mismatch")
    if split.get("development_validation_indices") != [V11_DEV]:
        raise SystemExit("v11 hotspot development holdout mismatch")
    hotspot_root = V11 / "round7_pod_snap/models/v11_candidate/active_si_hotspots/labels"
    hotspot_records = {int(row["index"]): row for row in split["labels"]}
    for i in V11_TRAIN + [V11_DEV]:
        record = hotspot_records[i]
        path = hotspot_root / f"si_hotspot_{i:02d}_pw350_gamma.extxyz"
        if sha(path) != record["files"]["extxyz"]:
            raise SystemExit(f"v11 hotspot label hash mismatch: {path.name}")
        atoms = read(path)
        center_id = int(record["center_atom_id"])
        ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
        center = np.flatnonzero(ids == center_id)
        if len(center) != 1 or atoms.get_chemical_symbols()[int(center[0])] != "Si":
            raise SystemExit(f"v11 center atom is ambiguous or not Si: {path.name}")
        idx = int(center[0])
        expected_force = np.asarray(record["center_force_DFT_eV_A"], dtype=float)
        if not np.allclose(atoms.arrays["REF_forces"][idx], expected_force, rtol=0, atol=1e-7):
            raise SystemExit(f"v11 center force differs from the locked audit: {path.name}")
        pos = np.asarray(atoms.positions[idx], dtype=float)
        lengths = np.asarray(atoms.cell.lengths(), dtype=float)
        if np.any(pos < 6.0) or np.any(lengths - pos < 6.0):
            raise SystemExit(f"v11 center lacks the registered 6 Å local shell: {path.name}")
        atoms.info["force_target_center_atom_id"] = center_id
        atoms.info["training_label_policy"] = (
            "center_Si_force_only; crop energy and all noncenter forces have zero fit weight"
        )
        atoms.info["dataset_role"] = (
            "training_acquisition_candidate" if i in V11_TRAIN else "development_hotspot_holdout"
        )
        if i in V11_TRAIN:
            train.append(atoms)
            train_roles.append({"role": "v11_center_Si_force_only_training_candidate",
                                "config_type": atoms.info.get("config_type"),
                                "center_atom_id": center_id,
                                "label_sha256": sha(path),
                                "crop_energy_weight": 0.0,
                                "noncenter_force_weight": 0.0})
        else:
            dev_hotspot = atoms

    OUT.mkdir(parents=True, exist_ok=True)
    train_out = OUT / "train_acquisition.extxyz"
    dev_out = OUT / "dev_periodic.extxyz"
    hotspot_dev_out = OUT / "dev_hotspot01_center_force.extxyz"
    write(train_out, train, format="extxyz")
    write(dev_out, dev, format="extxyz")
    write(hotspot_dev_out, dev_hotspot, format="extxyz")
    manifest = {
        "status": "exploratory_fit_candidate_pending_model_fit",
        "base_train_sha256": sha(base_train_path),
        "base_dev_sha256": sha(base_dev_path),
        "base_training_frames": len(read(base_train_path, index=":")),
        "v24_acquisition_labels": SELECTED_V24,
        "v11_training_hotspot_indices": V11_TRAIN,
        "v11_development_hotspot_index": V11_DEV,
        "v11_split_is_blind_test": False,
        "training_frames": len(train),
        "periodic_development_frames": len(dev),
        "separate_hotspot_development_frames": 1,
        "crop_energy_and_noncenter_force_targets": "zero fit weight",
        "v24_energy_convention": "native extrapolated energy stored in REF_energy; finite-smearing free energy retained separately",
        "training_roles": train_roles,
        "files": {p.name: sha(p) for p in [train_out, dev_out, hotspot_dev_out]},
        "source_commit": "2c68ef6a3ce37a339770c74c3a2cccb15ba9b5ba",
        "v11_source_commit": "beb4c8794c3edbb754d0efc8c00090184ca852b4",
        "blind_test_opened": False,
    }
    (OUT / "candidate_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
