#!/usr/bin/env python3
"""从 8 ps 完整原子+TTM 检查点继续进行同一热源协议的加热。"""

import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
# 首段驱动器把输入文件哈希记录为相对于势函数模型目录的路径，续算时用同一根目录核验。
MODEL_DIR = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
B_DIR = MODEL_DIR / "potential_validation/geometry_constraint_B"
CHECKPOINT_DIR = Path(os.environ.get("AG_TTM_CHECKPOINT_DIR", str(HERE / "checkpoint_latest")))
OUTPUT_DIR = Path(os.environ.get("AG_TTM_CONTINUE_OUTPUT_DIR", str(HERE / "continued_from_8ps")))
PARAMETERS_CSV = MODEL_DIR / "ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv"

sys.path.insert(0, str(MODEL_DIR))
from mpi4py import MPI  # noqa: E402
from lammps import lammps  # noqa: E402
from ag_ttm_reference import (  # noqa: E402
    AgParameters,
    AgTTMGrid,
    AtomGridMap,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
    KB_EV_K,
    M_VV2E_METAL,
)


def sha256(path):
    """计算文件 SHA256，用于拒绝损坏或不匹配的检查点。"""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    meta_path = CHECKPOINT_DIR / "CHECKPOINT.json"
    state_path = CHECKPOINT_DIR / "ttm_state.npz"
    restart_path = CHECKPOINT_DIR / "lammps.restart"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if sha256(state_path) != meta["ttm_state_sha256"]:
        raise RuntimeError("TTM 检查点 SHA256 不匹配")
    if sha256(restart_path) != meta["lammps_restart_sha256"]:
        raise RuntimeError("LAMMPS restart SHA256 不匹配")
    for relative_path, expected_hash in meta.get("input_sha256", {}).items():
        input_path = MODEL_DIR / relative_path
        if not input_path.exists() or sha256(input_path) != expected_hash:
            raise RuntimeError(f"续算输入文件缺失或哈希不匹配：{relative_path}")

    nsteps = int(os.environ.get("AG_TTM_CONTINUE_STEPS", "2000"))
    if nsteps <= 0:
        raise ValueError("AG_TTM_CONTINUE_STEPS 必须为正整数")
    with np.load(state_path, allow_pickle=False) as state:
        old_step = int(state["completed_step"])
        dt_ps = float(state["dt_ps"])
        source_scale = float(state["source_scale"])
        source_total = float(state["source_total_eV"])
        transfer_total = float(state["electron_to_lattice_total_eV"])
        e_md0 = float(state["initial_md_energy_eV"])
        e_electron0 = float(state["initial_electron_energy_eV"])
        ttm_time_ps = float(state["ttm_time_ps"])
        spacing = tuple(float(q) for q in state["grid_spacing_A"])
        fraction_ag = np.array(state["fraction_ag"], copy=True)
        source_power = np.array(state["source_power_eV_ps"], copy=True)
        temperature = np.array(state["grid_temperature_K"], copy=True)
        atom_ids_ttm = np.array(state["atom_ids_ttm"], dtype=np.int64, copy=True)
        atom_types_ttm = np.array(state["atom_types_ttm"], dtype=np.int32, copy=True)
        bottom_ag_ids = np.array(state["bottom_ag_ids"], dtype=np.int64, copy=True)
        map_cell_ids = np.array(state["map_cell_ids"], dtype=np.int64, copy=True)
        map_weights = np.array(state["map_weights"], dtype=float, copy=True)
        force_max_seen = float(state["force_max_seen_eV_A"])

    if not np.isclose(dt_ps, 0.0005, rtol=0.0, atol=1.0e-15):
        raise RuntimeError(f"该续算脚本只接受原有 0.5 fs 步长，检查点为 {dt_ps} ps")

    # 读入二进制 restart 后重新施加势、分组和积分 fix；LAMMPS restart 不保存 fix 定义。
    # 势文件在模型目录下以相对路径引用，续算工作目录需与原始驱动一致。
    os.chdir(MODEL_DIR)
    lmp = lammps(cmdargs=["-screen", "none", "-log", "none", "-sf", "omp"], comm=MPI.COMM_WORLD)
    for command in [
        f"read_restart {restart_path}",
        f"include {B_DIR}/potential.geometry_wca_fixed_skeleton.inc",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "group Ag type 4",
        "group skeleton type 1 2 3",
    ]:
        lmp.command(command)
    # 恢复原始被冻结的底部 Ag 原子 ID，不按已经移动后的坐标重新分组。
    lmp.command("group bottom_Ag id " + " ".join(str(int(q)) for q in bottom_ag_ids))
    for command in [
        "group mobile_Ag subtract Ag bottom_Ag",
        "fix freeze_skeleton skeleton setforce 0.0 0.0 0.0",
        "fix freeze_bottom_Ag bottom_Ag setforce 0.0 0.0 0.0",
        "mass 4 107.8682",
        f"timestep {dt_ps}",
        "fix integrate mobile_Ag nve",
        "thermo 1",
        "thermo_style custom step pe ke etotal temp",
        "thermo_modify lost error flush yes",
        "run 0",
    ]:
        lmp.command(command)

    nlocal = int(lmp.extract_setting("nlocal"))
    ids = np.asarray(lmp.numpy.extract_atom("id")[:nlocal], dtype=np.int64)
    types = np.asarray(lmp.numpy.extract_atom("type")[:nlocal], dtype=np.int32)
    x = lmp.numpy.extract_atom("x")
    v = lmp.numpy.extract_atom("v")
    f = lmp.numpy.extract_atom("f")
    id_to_index = {int(atom_id): i for i, atom_id in enumerate(ids)}
    missing = [int(atom_id) for atom_id in atom_ids_ttm if int(atom_id) not in id_to_index]
    if missing:
        raise RuntimeError(f"restart 缺少 {len(missing)} 个 TTM Ag 原子 ID")
    newton_idx = np.asarray([id_to_index[int(atom_id)] for atom_id in atom_ids_ttm], dtype=int)
    if not np.array_equal(types[newton_idx], atom_types_ttm) or np.any(types[newton_idx] != 4):
        raise RuntimeError("restart 中 TTM 原子元素映射不匹配")
    if len(atom_ids_ttm) != len(map_cell_ids):
        raise RuntimeError("TTM 原子 ID 与 atom-to-grid 映射长度不一致")
    bottom_id_set = set(int(q) for q in bottom_ag_ids)
    mobile_idx = np.asarray(
        [i for i, (atom_id, atom_type) in enumerate(zip(ids, types))
         if atom_type == 4 and int(atom_id) not in bottom_id_set], dtype=int
    )

    parameters = AgParameters.from_csv(PARAMETERS_CSV, allow_low_temperature_extension=True)
    grid = AgTTMGrid(parameters, fraction_ag, spacing, source_power, initial_temperature_K=300.0)
    if grid.temperature_K.shape != temperature.shape:
        raise RuntimeError("TTM 温度网格形状不匹配")
    grid.temperature_K = temperature
    grid.time_ps = ttm_time_ps
    atom_map = AtomGridMap(map_cell_ids, map_weights, tuple(int(q) for q in fraction_ag.shape))
    masses = np.full(len(newton_idx), 107.8682)
    active_cells = fraction_ag.ravel() > 0.0
    hot_cells = source_power.ravel() >= 0.01 * source_power.max()
    source_weights = source_power.ravel()[hot_cells]
    start_step = int(lmp.get_thermo("step"))
    if start_step != old_step:
        raise RuntimeError(f"LAMMPS step={start_step} 与 TTM checkpoint step={old_step} 不一致")
    if not np.isclose(float(lmp.get_thermo("etotal")),
                      float(meta.get("md_energy_at_checkpoint_eV", lmp.get_thermo("etotal"))),
                      rtol=1.0e-10, atol=1.0e-7):
        raise RuntimeError("LAMMPS restart 能量与检查点报告不一致")

    output_metrics = OUTPUT_DIR / "continuation_metrics.csv"
    dt_ps = 0.0005
    started = time.monotonic()
    with output_metrics.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step", "time_ps", "T_newton_Ag_K", "Tl_source_weighted_K",
                         "Te_max_K", "Fmax_mobile_Ag_eV_A", "source_total_eV",
                         "electron_to_lattice_eV", "energy_ledger_residual_eV"])
        for offset in range(1, nsteps + 1):
            lmp.command("run 1 pre no post no")
            step = start_step + offset
            e_md = float(lmp.get_thermo("etotal"))
            velocity = np.array(v[:nlocal][newton_idx], copy=True)
            fmax = float(np.linalg.norm(np.asarray(f[:nlocal][mobile_idx]), axis=1).max())
            force_max_seen = max(force_max_seen, fmax)
            if not np.isfinite(fmax) or fmax > 50.0:
                raise RuntimeError(f"续算安全停止：step={step}, Fmax={fmax:.6g} eV/A")
            tl = lattice_temperature_by_cell(velocity, masses, atom_map, smoothing_radius_cells=1)
            increment = grid.advance(dt_ps, tl, heat_on=True)
            velocity_new, applied = rescale_lattice_kinetic(
                velocity, masses, atom_map, increment.lattice_energy_per_cell_eV
            )
            v[:nlocal][newton_idx] = velocity_new
            source_total += increment.source_to_electron_eV
            transfer_total += float(applied.sum())
            te_max = float(grid.temperature_K.max())
            if te_max > 20000.0:
                raise RuntimeError(f"TTM 温度超过表格上限：{te_max:.6g} K")
            tl_flat = lattice_temperature_by_cell(
                velocity_new, masses, atom_map, smoothing_radius_cells=1
            ).ravel()
            tl_source = float(np.average(tl_flat[hot_cells], weights=source_weights))
            kinetic = 0.5 * M_VV2E_METAL * np.sum(masses[:, None] * velocity_new**2)
            temp_ag = float(2.0 * kinetic / (3.0 * len(newton_idx) * KB_EV_K))
            combined = e_md + float(applied.sum()) + float(grid.electron_energy_eV.sum())
            residual = combined - (e_md0 + e_electron0) - source_total
            if offset % 20 == 0 or offset == nsteps:
                writer.writerow([step, f"{step * dt_ps:.6f}", f"{temp_ag:.8f}",
                                 f"{tl_source:.8f}", f"{te_max:.8f}", f"{fmax:.8f}",
                                 f"{source_total:.10f}", f"{transfer_total:.10f}", f"{residual:.12g}"])
                stream.flush()
                print(f"step={step} t={step*dt_ps:.5f}ps T_Ag={temp_ag:.2f}K "
                      f"Tl_source={tl_source:.2f}K Te_max={te_max:.2f}K "
                      f"Fmax={fmax:.4f}eV/A source={source_total:.3f}eV "
                      f"ledger={residual:.6g}eV", flush=True)

    final_step = start_step + nsteps
    lmp.command("run 0 pre no post no")
    frame_name = f"frame_{final_step * dt_ps:06.3f}ps_step{final_step:05d}.dump.gz"
    lmp.command(f"write_dump all custom/gz {OUTPUT_DIR / frame_name} "
                "id type x y z vx vy vz fx fy fz modify sort id")
    next_restart = OUTPUT_DIR / "lammps.restart"
    lmp.command(f"write_restart {next_restart}")
    new_state_path = OUTPUT_DIR / "ttm_state.npz"
    with new_state_path.open("wb") as stream:
        np.savez_compressed(
            stream, completed_step=np.asarray(final_step), dt_ps=np.asarray(dt_ps),
            source_scale=np.asarray(source_scale), source_total_eV=np.asarray(source_total),
            electron_to_lattice_total_eV=np.asarray(transfer_total),
            initial_md_energy_eV=np.asarray(e_md0), initial_electron_energy_eV=np.asarray(e_electron0),
            ttm_time_ps=np.asarray(grid.time_ps), grid_temperature_K=grid.temperature_K,
            fraction_ag=grid.f_ag, source_power_eV_ps=grid.source_power,
            grid_spacing_A=grid.spacing_A, atom_ids_ttm=atom_ids_ttm,
            atom_types_ttm=atom_types_ttm, map_cell_ids=atom_map.cell_ids,
            map_weights=atom_map.weights, bottom_ag_ids=bottom_ag_ids,
            force_max_seen_eV_A=np.asarray(force_max_seen),
        )
    next_meta = {
        "checkpoint_time_ps": final_step * dt_ps,
        "completed_step": final_step,
        "dt_ps": dt_ps,
        "atom_count": int(lmp.get_natoms()),
        "ttm_atom_count": int(len(atom_ids_ttm)),
        "grid_shape": list(grid.temperature_K.shape),
        "grid_spacing_A": [float(q) for q in grid.spacing_A],
        "ttm_time_ps": float(grid.time_ps),
        "source_total_eV": float(source_total),
        "electron_to_lattice_total_eV": float(transfer_total),
        "initial_md_energy_eV": float(e_md0),
        "initial_electron_energy_eV": float(e_electron0),
        "md_energy_at_checkpoint_eV": float(lmp.get_thermo("etotal")),
        "mobile_ag_force_max_seen_eV_A": float(force_max_seen),
        "input_sha256": meta.get("input_sha256", {}),
        "lammps_restart_sha256": sha256(next_restart),
        "ttm_state_sha256": sha256(new_state_path),
        "atomic_frame": frame_name,
        "restart_note": "LAMMPS 原子状态与 TTM 电子网格状态成对保存，可再次用续算脚本读取。",
    }
    (OUTPUT_DIR / "CHECKPOINT.json").write_text(
        json.dumps(next_meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    result = {
        "status": "COMPLETED",
        "start_step": start_step,
        "final_step": final_step,
        "final_time_ps": final_step * dt_ps,
        "source_total_eV": source_total,
        "electron_to_lattice_total_eV": transfer_total,
        "energy_ledger_residual_eV": residual,
        "mobile_ag_force_max_seen_eV_A": force_max_seen,
        "runtime_s": time.monotonic() - started,
        "restart_file": str(next_restart),
        "ttm_state_file": str(new_state_path),
    }
    (OUTPUT_DIR / "continuation_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("FINAL " + json.dumps(result, ensure_ascii=False), flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
