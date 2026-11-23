#!/usr/bin/env python3
"""从仓库已归档的 Ag 纯相数值控制中抽取训练候选，不改动原记录或既有划分。"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from ase.io import read, write

REPO = Path(__file__).resolve().parents[3]
CONTROL = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls"
OUTDIR = Path(__file__).resolve().parent
SELECTED = [
    "Ag_baseline_k8x8x8_v19",
    "Ag_displacement_m01_k6x6x6_v19",
    "Ag_displacement_p01_k6x6x6_v19",
    "Ag_volume_m02_k8x8x8_v19",
    "Ag_volume_p02_k8x8x8_v19",
    "Ag_shear_m005_k8x8x8_v19",
    "Ag_shear_p005_k8x8x8_v19",
]

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def verify_manifest(path: Path) -> int:
    """逐项核验结果目录自身的完整性清单。"""
    checked = 0
    for line in (path / "SHA256SUMS.txt").read_text().splitlines():
        if not line.strip():
            continue
        expected, rel = line.split(maxsplit=1)
        rel = rel.strip().lstrip("*")
        target = path / rel
        if not target.is_file() or sha(target) != expected:
            raise RuntimeError(f"SHA-256 mismatch: {target}")
        checked += 1
    if checked == 0:
        raise RuntimeError(f"empty checksum manifest: {path}")
    return checked

manifest = json.loads((CONTROL / "input_manifest.json").read_text())
record_by_label = {row["label"]: row for row in manifest["records"]}
candidate_frames = []
audit_rows = []
source_families = set()

for label in SELECTED:
    result_dir = CONTROL / "calculations" / label
    output = result_dir / f"{label}_PW_PBE.extxyz"
    summary_path = result_dir / "summary.json"
    verify_path = result_dir / "verification.json"
    checked_files = verify_manifest(result_dir)
    summary = json.loads(summary_path.read_text())
    verification = json.loads(verify_path.read_text())
    record = record_by_label[label]
    source = CONTROL / record["input"]
    output_hash = sha(output)
    if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
        raise RuntimeError(f"verification did not pass for {label}")
    if not summary.get("scf_converged"):
        raise RuntimeError(f"SCF is not converged for {label}")
    if sha(source) != summary["source_sha256"] or summary["source_sha256"] != record["input_sha256"]:
        raise RuntimeError(f"input hash mismatch for {label}")
    if output_hash != verification["sha256"][output.name]:
        raise RuntimeError(f"output hash mismatch for {label}")
    if not re.search(rf"Converged in\s+{summary['scf_iterations']}\s+steps", (result_dir / "gpaw.log").read_text()):
        raise RuntimeError(f"GPAW convergence record missing for {label}")

    atoms = read(output)
    forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    energy = float(atoms.info["PW_PBE_energy_eV"])
    if (len(atoms) != 32 or atoms.get_chemical_symbols() != ["Ag"] * 32
            or atoms.pbc.tolist() != [True, True, True]
            or forces.shape != (32, 3) or not np.isfinite(forces).all()
            or not np.isfinite(energy)):
        raise RuntimeError(f"structure or label values invalid for {label}")
    atoms.info.update({
        "REF_energy": energy,
        "source_DFT_label": label,
        "training_label_policy": "native extrapolated PBE/PW500 total energy; finite-smearing (0.10 eV) GPAW forces",
        "model_data_role": "exploratory_train_candidate_unassigned",
        "data_family": "Ag32_COD9009647_pure_solid_controls",
        "source_result_sha256": output_hash,
        "source_result_role": summary.get("dataset_role", "numerical control; split review pending"),
    })
    atoms.arrays["REF_forces"] = forces.copy()
    # 保留原有 PW_PBE_* 目标字段，便于与原始计算结果逐项对照。
    candidate_frames.append(atoms)
    source_families.add("Ag32_COD9009647_pure_solid_controls")
    audit_rows.append({
        "label": label, "source_file": str(output.relative_to(REPO)),
        "source_sha256": output_hash, "input_sha256": summary["source_sha256"],
        "verified_result_files": checked_files, "scf_iterations": summary["scf_iterations"],
        "kpts": summary["method"]["kpts"], "energy_eV_cell": energy,
        "Fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "role": "train candidate only; not assigned or merged",
    })

OUTDIR.mkdir(exist_ok=True)
out = OUTDIR / "Ag_solid_controls_train_candidate_unassigned.extxyz"
write(out, candidate_frames, format="extxyz")
report = {
    "status": "PASS_SOURCE_AND_DFT_AUDIT_SPLIT_PENDING",
    "frames": len(candidate_frames), "atoms_per_frame": 32,
    "source_family_count": len(source_families),
    "source_families": sorted(source_families),
    "selected_controls": audit_rows,
    "output": out.name, "output_sha256": sha(out),
    "limitations": [
        "all seven structures belong to one Ag32 crystalline source family; never split across train and development/test",
        "static pristine/displacement/volume/shear controls do not provide liquid or solid-liquid coexistence coverage",
        "these controls were originally registered as numerical controls; this curated copy is only a training candidate pending project-level split review",
    ],
}
(OUTDIR / "curation_report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
