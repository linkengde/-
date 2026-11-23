#!/workspace/.venvs/gpaw-diagnosis/bin/python
"""审查纯 Ag 液态候选的来源、邻居距离、位移和径向分布。"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from ase.io import read
from ase.geometry import find_mic

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
source = ROOT / manifest["input"]
assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest["input_sha256"]
input_atoms = read(source)
frames = read(ROOT / "ag_liquid_sample.lammpstrj", index=":")
assert len(input_atoms) == 32 and len(frames) >= 3
for atoms in frames:
    atoms.set_chemical_symbols(["Ag"] * len(atoms))
assert input_atoms.get_chemical_symbols() == ["Ag"] * 32
assert np.allclose(input_atoms.positions, frames[-1].positions, rtol=0, atol=1e-8)

# LAMMPS dump 在生成脚本中按 atom ID 排序，故第 i 行对应 ID=i+1。
# 输入文件没有逐原子 ID 数组；此映射由 dump_modify sort id 和逐帧坐标核对确定。
source_ids = list(range(1, len(input_atoms) + 1))
D = input_atoms.get_all_distances(mic=True)
np.fill_diagonal(D, np.inf)
nearest = D.min(axis=1)
pairs = D[np.triu_indices(len(input_atoms), 1)]

# dump 每 1000 步保存一次，时间步长 0.5 fs；末尾相隔 2000 步即 1 ps。
prior = frames[-3]
last = frames[-1]
delta = last.positions - prior.positions
mic_delta, _ = find_mic(delta, last.cell, pbc=last.pbc)
displacements = np.linalg.norm(mic_delta, axis=1)

# 从 LAMMPS CNA 输出提取类别；5 表示 unknown，不能单独视为液态证明。
dump_lines = (ROOT / "candidate_structure.dump").read_text().splitlines()
atom_header = next(i for i, line in enumerate(dump_lines) if line.startswith("ITEM: ATOMS"))
columns = dump_lines[atom_header].split()[2:]
cna_i = columns.index("c_cna")
rows = [line.split() for line in dump_lines[atom_header + 1:] if line.strip()]
cna = [int(float(row[cna_i])) for row in rows]
unique, counts = np.unique(cna, return_counts=True)

# LAMMPS ASE reader不保留这个自定义dump中的力列；直接从原始行读取EAM力范数。
trajectory_lines = (ROOT / "ag_liquid_sample.lammpstrj").read_text().splitlines()
trajectory_force_max = []
i = 0
while i < len(trajectory_lines):
    if trajectory_lines[i].startswith("ITEM: TIMESTEP"):
        step = int(trajectory_lines[i + 1])
        natoms = int(trajectory_lines[i + 3])
        header = trajectory_lines[i + 8].split()[2:]
        index = {name: j for j, name in enumerate(header)}
        data = np.asarray([[float(x) for x in line.split()]
                           for line in trajectory_lines[i + 9:i + 9 + natoms]])
        force = data[:, [index["fx"], index["fy"], index["fz"]]]
        trajectory_force_max.append({
            "step": step,
            "EAM_Fmax_eV_A": float(np.linalg.norm(force, axis=1).max()),
        })
        i += 9 + natoms
    else:
        i += 1

volume = float(input_atoms.get_volume())
density = len(input_atoms) * 107.8682 / 0.602214076 / volume
# 立方体周期小胞的 RDF 仅统计最小镜像距离不超过半盒长的区间。
width = 0.05
rmax = min(input_atoms.cell.lengths()) / 2
edges = np.arange(0.0, rmax + width, width)
hist, edges = np.histogram(pairs, bins=edges)
centers = 0.5 * (edges[:-1] + edges[1:])
shell = 4 * np.pi / 3 * (edges[1:] ** 3 - edges[:-1] ** 3)
ideal = len(input_atoms) * density / 2 * shell
rdf = np.divide(hist, ideal, out=np.zeros_like(hist, dtype=float), where=ideal > 0)
np.savetxt(ROOT / "rdf.csv", np.c_[centers, rdf], delimiter=",", header="r_A,g_r", comments="")

summary = {
    "status": "liquid_like_geometry_candidate_only",
    "DFT_label_status": "pending",
    "atom_count": len(input_atoms),
    "source_snapshot_step": int(manifest["source_snapshot"]["timestep"]),
    "source_atom_id_mapping": "input row i maps to LAMMPS atom ID i+1; source dump was sorted by ID and coordinates were checked against selected frame",
    "source_atom_ids": source_ids,
    "density_g_cm3": density,
    "minimum_pair_distance_A": float(pairs.min()),
    "nearest_neighbor_percentiles_A": {str(p): float(x) for p, x in zip([0, 10, 50, 90, 100], np.percentile(nearest, [0, 10, 50, 90, 100]))},
    "displacement_1ps_percentiles_A": {str(p): float(x) for p, x in zip([50, 90, 100], np.percentile(displacements, [50, 90, 100]))},
    "cna_counts": {str(int(k)): int(v) for k, v in zip(unique, counts)},
    "sampled_trajectory_EAM_force_max": trajectory_force_max,
    "thermo_note": "generate_liquid.stdout records a one-frame instantaneous T=8002.44 K at step 40000; the selected step 44000 is 2 ps later at T=1384.28 K. This is a small-cell exploratory snapshot, not a claim of a fully equilibrated liquid ensemble.",
    "phase_screen_note": "CNA unknown class and 1 ps mobility support a liquid-like EAM-sampled geometry; EAM is only a structure proposer, not a DFT reference. Use several independent snapshots and DFT trajectories for final phase coverage.",
    "input_sha256": manifest["input_sha256"],
    "rdf_file": "rdf.csv",
}
(ROOT / "candidate_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
