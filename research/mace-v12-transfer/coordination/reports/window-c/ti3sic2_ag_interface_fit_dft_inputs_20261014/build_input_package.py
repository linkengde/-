#!/usr/bin/env python3
"""Create geometry-only DFT input candidates; never launches GPAW."""
from __future__ import annotations

import hashlib
import csv
import json
import shutil
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
STAGE_C = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_DFT_extension_20261014"
SOURCE_217 = ROOT / "source/Ag_Ti3SiC2_12p5198deg_217atom_dft_lite.extxyz"
INTERFACE_SOURCE_SHA = "945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2"
INTERFACE_SOURCE_COMMIT = "f1a080e14cd5ffd9739dcd2a3df2d188c853a8a8"
INTERFACE_SOURCE_PATH = "research/mace-v12-transfer/coordination/reports/window-c/low_mismatch_217_dft_lite_20261013/Ag_Ti3SiC2_12p5198deg_217atom_dft_lite.extxyz"
TSC_FAMILY = "TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY"
IFACE_FAMILY = "AG_TSC_SI_TERM_TWIST12P51983_217_PARENT_FAMILY"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def min_distance(atoms, idx_a, idx_b):
    pos = atoms.positions
    pairs = []
    for i in idx_a:
        for j in idx_b:
            dvec, _ = find_mic(pos[j] - pos[i], atoms.cell, atoms.pbc)
            pairs.append((float(np.linalg.norm(dvec)), int(i), int(j)))
    return sorted(pairs)


def method_for(kind: str) -> dict:
    if kind == "bulk_tsc":
        return {
            "code": "GPAW 26.7.0 / ASE 3.29.0 / gpaw-data 1.2.1",
            "resources": {"mpi_ranks": 4, "omp_threads": 1},
            "paw_setup_sha256": {
                "Ti.PBE.gz": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
                "Si.PBE.gz": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
                "C.PBE.gz": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
            },
            "xc": "PBE",
            "mode": "PW",
            "cutoff_eV": 500,
            "kpoints": [10, 10, 4],
            "kpoint_convention": "GPAW default Monkhorst-Pack integer mesh; no explicit gamma override or extra offset",
            "occupations": {"type": "FermiDirac", "width_eV": 0.10},
            "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
            "convergence": {"energy_eV_per_electron": 1e-5, "density": 1e-5, "eigenstates_eV": 1e-5},
            "maxiter": 160,
            "geometry": "fixed ions and fixed cell; single point",
            "outputs": ["native_energy_eV", "force_consistent_free_energy_eV", "forces_eV_per_A", "stress_eV_per_A3_full_tensor"],
        }
    return {
        "code": "GPAW 26.7.0 / ASE 3.29.0 / gpaw-data 1.2.1",
        "resources": {"mpi_ranks": 4, "omp_threads": 1},
        "paw_setup_sha256": {
            "Ag.PBE.gz": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
            "Ti.PBE.gz": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
            "Si.PBE.gz": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
            "C.PBE.gz": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
        },
        "xc": "PBE",
        "mode": "PW",
        "cutoff_eV": 500,
        "kpoints": [2, 2, 1],
        "kpoint_convention": "Gamma-centered mesh as proposed for the 217-atom DFT-lite cell; convergence NOT demonstrated",
        "occupations": {"type": "FermiDirac", "width_eV": 0.10},
        "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
        "convergence": {"energy_eV_per_electron": 1e-5, "density": 1e-5, "eigenstates_eV": 1e-5},
        "maxiter": 160,
        "geometry": "fixed ions and fixed cell; rigid Ag z translation only between input candidates",
        "outputs": ["native_energy_eV", "force_consistent_free_energy_eV", "forces_eV_per_A", "stress_eV_per_A3_full_tensor_diagnostic_only"],
        "stress_note": "3D stress of a slab with vacuum is not a bulk stress target; retain only as diagnostic unless a slab convention is approved.",
    }


def main() -> None:
    if sha(SOURCE_217) != INTERFACE_SOURCE_SHA:
        raise SystemExit("217-atom source SHA mismatch")
    bulk_specs = [
        ("TSC_VOL_M02", "A_VOL_M02.extxyz", "C-A-VOL-M02", "isotropic volume -2%"),
        ("TSC_VOL_P02", "A_VOL_P02.extxyz", "C-A-VOL-P02", "isotropic volume +2%"),
        ("TSC_SHEAR_M005", "A_SHEAR_M005.extxyz", "C-A-SHEAR-M005", "shear -0.5%"),
        ("TSC_SHEAR_P005", "A_SHEAR_P005.extxyz", "C-A-SHEAR-P005", "shear +0.5%"),
        ("TSC_SI_ID0007_Z_M020", "B_SI_ID0007_Z_M020.extxyz", "C-B-SI-ID0007-Z-M020", "Si local_new_id=7 z -0.02 A at fixed cell"),
        ("TSC_SI_ID0007_Z_P020", "B_SI_ID0007_Z_P020.extxyz", "C-B-SI-ID0007-Z-P020", "Si local_new_id=7 z +0.02 A at fixed cell"),
    ]
    records = []
    for structure_id, filename, task_id, perturbation in bulk_specs:
        src = STAGE_C / "inputs" / filename
        dst = ROOT / "inputs/backbone" / filename
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        atoms = read(dst)
        if len(atoms) != 48 or atoms.get_chemical_formula() not in ("C16Si8Ti24", "Ti24Si8C16"):
            raise SystemExit(f"Unexpected TSC geometry identity: {filename}")
        config = {
            "task_id": task_id,
            "structure_id": structure_id,
            "structure_path": f"inputs/backbone/{filename}",
            "structure_sha256": sha(dst),
            "natoms": len(atoms),
            "composition": {s: int((atoms.symbols == s).sum()) for s in ("Ti", "Si", "C")},
            "parent_family": TSC_FAMILY,
            "source_geometry_sha256": sha(dst),
            "source_e7_commit": "5e3b3f7c31fd9bf6156cf339428095dea54c1de0",
            "perturbation": perturbation,
            "proposed_role": "TRAIN_CANDIDATE_PENDING_ROLE_FREEZE",
            "method": method_for("bulk_tsc"),
            "execution_authorized": False,
            "dft_label_status": "NOT_GENERATED",
        }
        config_path = ROOT / "configs" / f"{structure_id}.json"
        json_write(config_path, config)
        records.append({**config, "config_path": str(config_path.relative_to(ROOT)), "config_sha256": sha(config_path)})

    source_atoms = read(SOURCE_217)
    if len(source_atoms) != 217 or {s: int((source_atoms.symbols == s).sum()) for s in ("Ag", "Ti", "Si", "C")} != {"Ag": 84, "Ti": 57, "Si": 38, "C": 38}:
        raise SystemExit("217-atom interface source identity mismatch")
    source_info = dict(source_atoms.info)
    source_symbols = np.asarray(source_atoms.get_chemical_symbols())
    ag_idx = np.where(source_symbols == "Ag")[0]
    tsc_idx = np.where(source_symbols != "Ag")[0]
    top_tsc_z = float(source_atoms.positions[tsc_idx, 2].max())
    ag_bottom_z = float(source_atoms.positions[ag_idx, 2].min())
    observed_gap = ag_bottom_z - top_tsc_z
    if abs(observed_gap - 2.0) > 1e-6:
        raise SystemExit(f"Unexpected source plane gap: {observed_gap}")

    interface_records = []
    # Coarse rigid-layer scan. It is a geometry proposal, not an optimized d0.
    for gap in np.round(np.arange(1.8, 3.401, 0.2), 1):
        atoms = source_atoms.copy()
        delta = float(gap - observed_gap)
        atoms.positions[ag_idx, 2] += delta
        structure_id = f"AGTSC_SI_GAP_{gap:.1f}A".replace(".", "p")
        task_id = f"C-IFACE-SI-GAP-{gap:.1f}".replace(".", "p")
        atoms.info.update(source_info)
        atoms.info.update({
            "structure_id": structure_id,
            "source_family": IFACE_FAMILY,
            "role_proposal": "TRAIN_CANDIDATE_FOR_FUTURE_UNIFIED_FOUR_ELEMENT_MODEL",
            "target_interface_plane_gap_A": float(gap),
            "rigid_Ag_translation_z_A": delta,
            "source_structure_sha256": INTERFACE_SOURCE_SHA,
            "geometry_operation": "translate all Ag atoms rigidly along +z; no relaxation, no lateral shift, no cell change",
        })
        out = ROOT / "inputs/interface_scan" / f"{structure_id}.extxyz"
        out.parent.mkdir(parents=True, exist_ok=True)
        write(out, atoms, format="extxyz")
        # Re-read and ensure the serialization did not change atom order or identities.
        roundtrip = read(out)
        if len(roundtrip) != 217 or not np.array_equal(roundtrip.symbols.numbers, atoms.symbols.numbers):
            raise SystemExit(f"Round-trip structure identity failure: {structure_id}")
        roundtrip_symbols = np.asarray(roundtrip.get_chemical_symbols())
        cross = min_distance(roundtrip, np.where(roundtrip_symbols == "Ag")[0], np.where(roundtrip_symbols != "Ag")[0])
        minima_by_el = {}
        for element in ("Si", "Ti", "C"):
            idx = np.where(roundtrip_symbols == element)[0]
            pairs = min_distance(roundtrip, np.where(roundtrip_symbols == "Ag")[0], idx)
            minima_by_el[f"Ag-{element}"] = {"distance_A": pairs[0][0], "atom_i_0based": pairs[0][1], "atom_j_0based": pairs[0][2]}
        config = {
            "task_id": task_id,
            "structure_id": structure_id,
            "structure_path": str(out.relative_to(ROOT)),
            "structure_sha256": sha(out),
            "natoms": len(roundtrip),
            "composition": {s: int((roundtrip_symbols == s).sum()) for s in ("Ag", "Ti", "Si", "C")},
            "parent_family": IFACE_FAMILY,
            "source_geometry": {"commit": INTERFACE_SOURCE_COMMIT, "path": INTERFACE_SOURCE_PATH, "sha256": INTERFACE_SOURCE_SHA},
            "parent_394_geometry_sha256": source_info.get("parent_structure_sha256", "NOT_PRESENT"),
            "geometry_operation": {"target_gap_A": float(gap), "rigid_Ag_z_translation_A": delta, "cell_and_pbc_unchanged": True, "ions_unrelaxed": True},
            "minimum_cross_interface_distances": minima_by_el,
            "minimum_any_Ag_to_TSC_pair_A": cross[0][0],
            "proposed_role": "TRAIN_CANDIDATE_PENDING_ROLE_FREEZE",
            "fit_scope_note": "Interface labels cannot constrain the three-element Ti-Si-C TSCBOP; intended only for a future unified Ag-Ti-Si-C model.",
            "method": method_for("interface"),
            "execution_authorized": False,
            "dft_label_status": "NOT_GENERATED",
        }
        config_path = ROOT / "configs" / f"{structure_id}.json"
        json_write(config_path, config)
        interface_records.append({**config, "config_path": str(config_path.relative_to(ROOT)), "config_sha256": sha(config_path)})

    json_write(ROOT / "INPUT_MANIFEST.json", {
        "package_status": "GEOMETRY_AND_DFT_INPUTS_PREPARED_ONLY",
        "execution_authorized": False,
        "dft_started": False,
        "fitting_started": False,
        "tscbop_modified": False,
        "backbone_tasks": records,
        "interface_scan_tasks": interface_records,
    })
    with (ROOT / "STRUCTURE_FAMILY_LABELS.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["task_id", "structure_id", "parent_family", "source", "source_sha256", "geometry_operation", "natoms", "composition", "proposed_split", "current_role_or_status", "required_outputs", "execution_status", "notes"])
        legacy = [
            ("E5_TSC_BASELINE", "Ti3SiC2_baseline", "48-atom Kisi/COD baseline", "157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4", "baseline", "NUMERICAL_CONTROL_DIAGNOSTIC_NOT_TRAIN", "energy/free_energy/forces; no stress"),
            ("E5_TSC_TI4F_M020", "Ti4f_z_m020", "B archived fixed-cell response", "5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf", "Ti4f z -0.02 A", "NUMERICAL_CONTROL_DIAGNOSTIC_NOT_TRAIN", "energy/free_energy/forces; no stress"),
            ("E5_TSC_TI4F_P020", "Ti4f_z_p020", "B archived fixed-cell response", "5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b", "Ti4f z +0.02 A", "NUMERICAL_CONTROL_DIAGNOSTIC_NOT_TRAIN", "energy/free_energy/forces; no stress"),
            ("E5_TSC_C4F_M020", "C4f_z_m020", "B archived fixed-cell response", "19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f", "C4f z -0.02 A", "NUMERICAL_CONTROL_DIAGNOSTIC_NOT_TRAIN", "energy/free_energy/forces; no stress"),
            ("E5_TSC_C4F_P020", "C4f_z_p020", "B archived fixed-cell response", "6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72", "C4f z +0.02 A", "NUMERICAL_CONTROL_DIAGNOSTIC_NOT_TRAIN", "energy/free_energy/forces; no stress"),
            ("E5_TSC_INTERFACE_DERIVED_SLAB", "TSC_only_114", "E5 Ag-interface-derived TSC slab", "0cd573c024b7f792d68fb6bf9e46ae4b01a56ae65f29d86aec1726444ffa4682", "Ag removed from interface-derived slab; ancestry not independent", "FORCE_DIAGNOSTIC_KEEP_CURRENT_ROLE", "force label only; energy/free_energy/stress unverified"),
        ]
        for task, sid, source, geom_sha, operation, role_status, outputs in legacy:
            family = "TSC_INTERFACE_DERIVED_SLAB_FAMILY" if task.endswith("SLAB") else TSC_FAMILY
            natoms = "114" if task.endswith("SLAB") else "48"
            composition = "Ti57Si19C38" if task.endswith("SLAB") else "Ti24Si8C16"
            writer.writerow([task, sid, family, f"Historical E5/B archive; {source}", geom_sha, operation, natoms, composition, "EXCLUDE_FROM_DEV_TEST; possible future TRAIN only by separate review", role_status, outputs, "EXISTING_LABELS_NOT_COPIED", "Retain historical numerical-control/diagnostic role; previously used in E5 analysis; do not use as independent validation."])
        for row in records + interface_records:
            src = row["source_geometry_sha256"] if "source_geometry_sha256" in row else row["source_geometry"]["sha256"]
            if row["parent_family"] == TSC_FAMILY:
                source = "Stage C geometry-only candidate; E7 commit 5e3b3f7c31fd9bf6156cf339428095dea54c1de0"
                operation = row["perturbation"]
                notes = "Baseline-derived family; geometry only; keep entire family together. Existing numerical-control labels remain unchanged."
                outputs = "native energy; force-consistent free energy; forces; full stress"
            else:
                source = f"217-atom low-mismatch parent, commit {INTERFACE_SOURCE_COMMIT}, {INTERFACE_SOURCE_PATH}"
                operation = f"Rigid Ag z translation to plane gap {row['geometry_operation']['target_gap_A']:.1f} A; no relaxation"
                notes = "One interface family; TRAIN candidate only for future unified four-element potential, not TSCBOP; k-point convergence pending."
                outputs = "native energy; force-consistent free energy; forces; full stress diagnostic only"
            writer.writerow([
                row["task_id"], row["structure_id"], row["parent_family"], source, src,
                operation, row["natoms"], ";".join(f"{k}{v}" for k, v in row["composition"].items()),
                "TRAIN_CANDIDATE", "GEOMETRY_ONLY; pending role freeze", outputs, "NOT_RUN", notes,
            ])
    print(f"Prepared {len(records)} backbone and {len(interface_records)} interface geometry/config pairs.")


if __name__ == "__main__":
    main()
