#!/workspace/.venvs/gpaw-diagnosis/bin/python
"""将已通过标签审查的单帧标签映射为训练数据字段，但不并入任何数据划分。"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "dft_run"
audit = json.loads((RUN / "label_audit.json").read_text())
if audit["status"] != "PASS_LABEL_INTEGRITY_TRAIN_SPLIT_PENDING":
    raise SystemExit("标签完整性审查未通过，拒绝生成训练候选")
a = read(RUN / "audited_candidate.extxyz")
forces = np.asarray(a.arrays["PW_PBE_forces"], dtype=float)
energy = float(a.info["PW_PBE_energy_eV"])
if forces.shape != (len(a), 3) or not np.isfinite(forces).all() or not np.isfinite(energy):
    raise SystemExit("能量或力不完整，拒绝生成训练候选")
a.info["REF_energy"] = energy
a.arrays["REF_forces"] = forces.copy()
a.info["training_label_policy"] = "native extrapolated PBE/PW500 total energy; finite-smearing (0.10 eV) GPAW forces"
a.info["proposal_role"] = "train_candidate_unassigned"
a.info["model_data_role"] = "exploratory_training_candidate_not_yet_merged"
output = RUN / "training_candidate_unassigned.extxyz"
write(output, a, format="extxyz")
print(json.dumps({
    "output": output.name,
    "frames": 1,
    "natoms": len(a),
    "REF_energy_eV_cell": energy,
    "force_key": "REF_forces",
    "status": "candidate only; no train/dev assignment or merge performed",
}, indent=2))
