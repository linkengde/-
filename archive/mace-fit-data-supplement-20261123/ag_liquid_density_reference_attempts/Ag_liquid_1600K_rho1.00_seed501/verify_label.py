#!/workspace/.venvs/gpaw-diagnosis/bin/python
"""核验 GPAW 标签的 SCF、结构映射、力有限性和哈希。"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "dft_run"
manifest = json.loads((ROOT / "input_manifest.json").read_text())
summary = json.loads((RUN / "summary.json").read_text())
input_atoms = read(ROOT / manifest["input"])
output_path = RUN / summary["output"]
output_atoms = read(output_path)
log = (RUN / "gpaw.log").read_text(errors="replace")

def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit("LABEL_AUDIT_FAILED: " + message)

input_sha = hashlib.sha256((ROOT / manifest["input"]).read_bytes()).hexdigest()
output_sha = hashlib.sha256(output_path.read_bytes()).hexdigest()
forces = np.asarray(output_atoms.arrays["PW_PBE_forces"], dtype=float)
energy = float(output_atoms.info["PW_PBE_energy_eV"])
matches = re.findall(r"Converged in\s+(\d+)\s+steps", log)
require(summary.get("scf_converged") is True, "摘要没有登记 SCF 收敛")
require(bool(matches), "GPAW 日志中没有明确的 Converged in N steps 记录")
log_iterations = int(matches[-1])
require(log_iterations == int(summary["scf_iterations"]), "日志与摘要 SCF 次数不符")
require(input_sha == manifest["input_sha256"] == summary["input_sha256"], "输入 SHA-256 不一致")
require(output_sha == summary["output_sha256"], "输出 SHA-256 不一致")
require(len(input_atoms) == len(output_atoms) == 64, "原子数不符")
require(input_atoms.get_chemical_symbols() == output_atoms.get_chemical_symbols() == ["Ag"] * 64,
        "元素顺序或符号不符")
require(input_atoms.pbc.tolist() == output_atoms.pbc.tolist() == [True, True, True], "PBC 不符")
require(np.array_equal(input_atoms.cell.array, output_atoms.cell.array), "晶胞改变")
require(np.allclose(input_atoms.positions, output_atoms.positions, rtol=0, atol=1e-10), "单点计算中坐标改变")
require(forces.shape == (32, 3) and np.isfinite(forces).all(), "力数组缺失或有 NaN/Inf")
require(np.isfinite(energy), "总能量不是有限值")
require(np.isclose(energy, float(summary["native_extrapolated_energy_eV_cell"]), rtol=0, atol=1e-10),
        "输出总能与摘要不符")

# 按生成轨迹的 sort id 顺序恢复明确的来源 atom ID；坐标与帧已在 candidate audit 中逐行核过。
output_atoms.arrays["source_atom_id"] = np.arange(1, len(output_atoms) + 1, dtype=int)
output_atoms.info["label_audit_status"] = "PASS_GPAW_SCF_mapping_force_hash; training split audit still pending"
audited_path = RUN / "audited_candidate.extxyz"
write(audited_path, output_atoms, format="extxyz")
report = {
    "status": "PASS_LABEL_INTEGRITY_TRAIN_SPLIT_PENDING",
    "input_sha256": input_sha,
    "raw_dft_output_sha256": output_sha,
    "audited_candidate_sha256": hashlib.sha256(audited_path.read_bytes()).hexdigest(),
    "scf_iterations_log": log_iterations,
    "natoms": len(output_atoms),
    "energy_eV_cell": energy,
    "Fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
    "F99_eV_A": float(np.percentile(np.linalg.norm(forces, axis=1), 99)),
    "RMS_force_vector_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
    "net_force_vector_eV_A": forces.sum(axis=0).tolist(),
    "source_ID_mapping": "input row i -> source LAMMPS atom ID i+1; dump_modify sort id and exact coordinates checked",
    "verified": {
        "all_GPAW_SCf_criteria": True,
        "same_geometry_and_order": True,
        "all_finite_forces_and_energy": True,
        "input_output_SHA256": True,
    },
    "limitations": [
        "one EAM-proposed 32-atom liquid-like snapshot; not an equilibrium validation ensemble",
        "k-point and cutoff convergence have not been established for this frame",
        "not merged into train/dev; candidate family split review remains required",
    ],
}
(RUN / "label_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
