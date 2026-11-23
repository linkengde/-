#!/workspace/.venvs/gpaw-diagnosis/bin/python
"""为单个纯 Ag 液态候选构型生成可审查的 GPAW DFT 标签。"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "dft_run"
MANIFEST = json.loads((ROOT / "input_manifest.json").read_text())
SOURCE = ROOT / MANIFEST["input"]
LABEL = MANIFEST["label"]
assert world.size == int(MANIFEST["method"]["mpi_ranks"]), (
    f"需要 {MANIFEST['method']['mpi_ranks']} 个 MPI rank，当前为 {world.size}"
)
source_sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
assert source_sha == MANIFEST["input_sha256"], "输入结构哈希与登记清单不符"

# 所有 rank 共享同一锁；防止误启动第二个相同 DFT 任务。
lock_file = open(ROOT / "dft_label.lock", "a+")
if world.rank == 0:
    try:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("已有 Ag 液态 DFT 标签任务持有锁") from exc
    RUN.mkdir(exist_ok=True)
    if (RUN / "summary.json").exists():
        raise SystemExit("检测到已有正式标签摘要；拒绝重复覆盖")
world.barrier()

# 设置来自当前 GPAW 安装的标准 Ag PAW 数据路径。
import gpaw_data
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())

from ase.io import read, write
from gpaw import FermiDirac, GPAW, PW
from gpaw.mixer import Mixer

atoms = read(SOURCE)
assert len(atoms) == 64
assert atoms.get_chemical_symbols() == ["Ag"] * 64
assert atoms.pbc.tolist() == [True, True, True]
assert not any(k in atoms.info for k in ("energy", "REF_energy", "PW_PBE_energy_eV"))
assert not any(k in atoms.arrays for k in ("forces", "REF_forces", "PW_PBE_forces"))

started = time.time()
progress_path = RUN / "progress.json"
converged = {"value": False}

calc = GPAW(
    mode=PW(500), xc="PBE", kpts=(2, 2, 2), symmetry="off",
    occupations=FermiDirac(0.10),
    mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    maxiter=160,
    parallel={"kpt": 4, "domain": 1, "band": 1, "sl_auto": False},
    txt=str(RUN / "gpaw.log"),
)
atoms.calc = calc

if world.rank == 0:
    progress_path.write_text(json.dumps({
        "label": LABEL, "status": "running", "iteration": 0,
        "input_sha256": source_sha, "mpi_ranks": world.size,
        "started_epoch_s": started,
    }, indent=2) + "\n")

def scf_step(context):
    """每轮更新进度，每 10 轮通过临时文件安全替换检查点。"""
    iteration = int(context.niter)
    if world.rank == 0:
        progress_path.write_text(json.dumps({
            "label": LABEL, "status": "running", "iteration": iteration,
            "elapsed_s": time.time() - started, "input_sha256": source_sha,
            "mpi_ranks": world.size,
        }, indent=2) + "\n")
    if iteration and iteration % 10 == 0:
        temporary = RUN / "state.next.gpw"
        stable = RUN / "state.gpw"
        if world.rank == 0:
            temporary.unlink(missing_ok=True)
        world.barrier()
        calc.write(str(temporary), mode="all")
        world.barrier()
        if world.rank == 0:
            os.replace(temporary, stable)
            checkpoint_hash = hashlib.sha256(stable.read_bytes()).hexdigest()
            (RUN / "checkpoint.json").write_text(json.dumps({
                "label": LABEL, "input_sha256": source_sha,
                "checkpoint_sha256": checkpoint_hash,
                "checkpoint_bytes": stable.stat().st_size,
                "last_scf_iteration": iteration,
            }, indent=2) + "\n")
        world.barrier()

calc.hooks["scf_step"] = scf_step
calc.hooks["converged"] = lambda: converged.__setitem__("value", True)

# 只有 GPAW 明确报告达到全部 SCF 收敛准则后，才读取并保存能量、力。
native_energy = float(atoms.get_potential_energy(force_consistent=False))
free_energy = float(calc.results["free_energy"])
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(calc.get_number_of_iterations())
if not converged["value"]:
    if world.rank == 0:
        progress_path.write_text(json.dumps({
            "label": LABEL, "status": "not_converged", "iteration": iterations,
            "elapsed_s": time.time() - started, "input_sha256": source_sha,
        }, indent=2) + "\n")
    raise SystemExit(f"GPAW 未在 {iterations} 次迭代内收敛；没有写出 DFT 标签")
if not np.isfinite(native_energy) or not np.isfinite(free_energy) or not np.isfinite(forces).all():
    raise SystemExit("DFT 返回非有限能量或力；没有写出正式标签")

if world.rank == 0:
    output_atoms = atoms.copy()
    output_atoms.calc = None
    output_atoms.info.update({
        "config_type": LABEL,
        "PW_PBE_energy_eV": native_energy,
        "PW_PBE_free_energy_eV": free_energy,
        "label_status": "GPAW PBE/PW500 k=3x3x3 converged; independent audit pending",
        "model_data_role": "exploratory_training_candidate_not_yet_merged",
        "source_input_sha256": source_sha,
        "DFT_method": "GPAW PBE/PW500, Gamma-centered 3x3x3, Fermi-Dirac 0.10 eV",
    })
    output_atoms.arrays["PW_PBE_forces"] = forces.copy()
    output_path = RUN / f"{LABEL}_PW_PBE.extxyz"
    write(output_path, output_atoms, format="extxyz")
    summary = {
        "label": LABEL, "status": "computed_independent_audit_pending",
        "input": str(SOURCE.relative_to(ROOT)), "input_sha256": source_sha,
        "output": output_path.name,
        "output_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "atom_count": len(atoms), "formula": atoms.get_chemical_formula(),
        "cell_angstrom": atoms.cell.array.tolist(), "pbc": atoms.pbc.tolist(),
        "native_extrapolated_energy_eV_cell": native_energy,
        "free_energy_eV_cell": free_energy,
        "scf_converged": True, "scf_iterations": iterations,
        "max_force_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_vector_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "mean_force_vector_eV_A": forces.mean(axis=0).tolist(),
        "method": MANIFEST["method"],
        "elapsed_s": time.time() - started,
        "target_policy": "periodic full-cell energy and all-atom PBE forces; candidate only until label and split audit",
        "training_integration": "not performed",
    }
    (RUN / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    progress_path.write_text(json.dumps({
        "label": LABEL, "status": "complete", "iteration": iterations,
        "elapsed_s": summary["elapsed_s"], "scf_converged": True,
        "input_sha256": source_sha, "mpi_ranks": world.size,
    }, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
world.barrier()
