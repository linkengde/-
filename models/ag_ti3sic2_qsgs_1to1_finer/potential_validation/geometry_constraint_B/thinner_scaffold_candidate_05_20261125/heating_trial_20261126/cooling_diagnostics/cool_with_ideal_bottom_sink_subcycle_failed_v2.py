#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""停热后用 Ag TTM 网格底部 Dirichlet 热沉冷却；不改变刚性骨架理想化。"""

import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI
from lammps import lammps

HERE = Path(__file__).resolve().parent
ROOT = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
B_DIR = ROOT / "potential_validation/geometry_constraint_B"
CHECKPOINT_DIR = Path(os.environ["AG_TTM_CHECKPOINT_DIR"])
OUTPUT_DIR = Path(os.environ.get("AG_TTM_COOL_OUTPUT_DIR", str(HERE / "cooling_sink")))
BOTTOM_SINK_Z_MAX_A = float(os.environ.get("AG_TTM_BOTTOM_SINK_Z_MAX_A", "50.0"))
BOTTOM_SINK_TEMPERATURE_K = float(os.environ.get("AG_TTM_BOTTOM_SINK_TEMPERATURE_K", "300.0"))
PARAMETERS_CSV = ROOT / "ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv"

sys.path.insert(0, str(ROOT))
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
    """计算文件哈希，校验 LAMMPS 原子态和 TTM 网格检查点。"""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    meta_path = CHECKPOINT_DIR / "CHECKPOINT.json"
    state_path = CHECKPOINT_DIR / "ttm_state.npz"
    restart_path = CHECKPOINT_DIR / "lammps.restart"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if sha256(state_path) != meta["ttm_state_sha256"]:
        raise RuntimeError("输入 TTM 检查点 SHA256 不匹配")
    if sha256(restart_path) != meta["lammps_restart_sha256"]:
        raise RuntimeError("输入 LAMMPS 检查点 SHA256 不匹配")
    for relative_path, expected_hash in meta.get("input_sha256", {}).items():
        input_path = ROOT / relative_path
        if not input_path.exists() or sha256(input_path) != expected_hash:
            raise RuntimeError(f"输入依赖缺失或哈希不匹配：{relative_path}")

    with np.load(state_path, allow_pickle=False) as state:
        old_step = int(state["completed_step"])
        dt_ps = float(state["dt_ps"])
        source_scale = float(state["source_scale"])
        source_total = float(state["source_total_eV"])
        transfer_total = float(state["electron_to_lattice_total_eV"])
        sink_total = float(state["bottom_sink_removed_eV"]) if "bottom_sink_removed_eV" in state else 0.0
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
        raise RuntimeError(f"冷却要求原有 0.5 fs 步长，检查点为 {dt_ps} ps")
    if BOTTOM_SINK_Z_MAX_A <= 0.0:
        raise ValueError("底部热沉厚度必须为正")

    # LAMMPS 续算：固定骨架和最低 Ag 层；其他 Ag 自由运动，不再加激光热源。
    os.chdir(ROOT)
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

    nsteps = int(os.environ.get("AG_TTM_COOL_STEPS", "1000"))
    if nsteps <= 0:
        raise ValueError("AG_TTM_COOL_STEPS 必须为正整数")
    coupling_substeps = int(os.environ.get("AG_TTM_COOL_COUPLING_SUBSTEPS", "20"))
    if coupling_substeps <= 0:
        raise ValueError("AG_TTM_COOL_COUPLING_SUBSTEPS 必须为正整数")
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
        raise RuntimeError("restart 中 TTM 原子映射不匹配")
    bottom_id_set = set(int(q) for q in bottom_ag_ids)
    mobile_idx = np.asarray([
        i for i, (atom_id, atom_type) in enumerate(zip(ids, types))
        if atom_type == 4 and int(atom_id) not in bottom_id_set
    ], dtype=int)

    parameters = AgParameters.from_csv(PARAMETERS_CSV, allow_low_temperature_extension=True)
    grid = AgTTMGrid(parameters, fraction_ag, spacing, source_power, initial_temperature_K=300.0)
    if grid.temperature_K.shape != temperature.shape:
        raise RuntimeError("TTM 温度网格形状与检查点不匹配")
    grid.temperature_K = temperature
    grid.time_ps = ttm_time_ps
    atom_map = AtomGridMap(map_cell_ids, map_weights, tuple(int(q) for q in fraction_ag.shape))
    masses = np.full(len(newton_idx), 107.8682)

    # z=35 和 45 Å 是当前 Ag TTM 网格最靠底的两个有效单元，构成 300 K 理想 Dirichlet 热沉。
    nz = fraction_ag.shape[2]
    z_centers_A = (np.arange(nz, dtype=float) + 0.5) * spacing[2]
    sink_mask = np.broadcast_to(
        z_centers_A[None, None, :] <= BOTTOM_SINK_Z_MAX_A, fraction_ag.shape
    ).copy()
    sink_mask &= fraction_ag > 0.0
    if not np.any(sink_mask):
        raise RuntimeError("底部理想热沉区域没有 Ag 电子网格单元")
    start_step = int(lmp.get_thermo("step"))
    if start_step != old_step:
        raise RuntimeError(f"LAMMPS step={start_step} 与检查点 step={old_step} 不一致")
    if not np.isclose(float(lmp.get_thermo("etotal")),
                      float(meta.get("md_energy_at_checkpoint_eV", lmp.get_thermo("etotal"))),
                      rtol=1.0e-10, atol=1.0e-7):
        raise RuntimeError("LAMMPS restart 能量与检查点报告不一致")

    output_metrics = OUTPUT_DIR / "cooling_metrics.csv"
    started = time.monotonic()
    residual = float(meta.get("energy_ledger_residual_eV", 0.0))
    with output_metrics.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step", "time_ps", "T_ttm_Ag_K", "Te_grid_max_K",
                         "Te_bottom_mean_K", "Fmax_mobile_Ag_eV_A", "source_total_eV",
                         "electron_to_lattice_total_eV", "bottom_sink_removed_total_eV",
                         "energy_ledger_residual_eV"])
        for offset in range(1, nsteps + 1):
            lmp.command("run 1 pre no post no")
            step = start_step + offset
            e_md = float(lmp.get_thermo("etotal"))
            velocity = np.array(v[:nlocal][newton_idx], copy=True)
            fmax = float(np.linalg.norm(np.asarray(f[:nlocal][mobile_idx]), axis=1).max())
            force_max_seen = max(force_max_seen, fmax)
            if not np.isfinite(fmax) or fmax > 50.0:
                raise RuntimeError(f"冷却安全停止：step={step}, Fmax={fmax:.6g} eV/A")
            # 对 0.5 fs 的 MD 步把 TTM 耦合分成 20 次小交换，
            # 每次都更新原子速度和局部晶格温度，避免一次性抽走某个原子的全部动能。
            # 每个子步都保留 300 K 理想底部边界，热源保持关闭。
            coupling_dt_ps = dt_ps / coupling_substeps
            for _substep in range(coupling_substeps):
                tl = lattice_temperature_by_cell(
                    velocity, masses, atom_map, smoothing_radius_cells=1
                )
                increment = grid.advance(
                    coupling_dt_ps, tl, heat_on=False, coupling_on=True,
                    bottom_sink_mask=sink_mask,
                    bottom_sink_temperature_K=BOTTOM_SINK_TEMPERATURE_K,
                )
                velocity, applied = rescale_lattice_kinetic(
                    velocity, masses, atom_map,
                    increment.lattice_energy_per_cell_eV,
                )
                transfer_total += float(applied.sum())
                sink_total += float(increment.bottom_sink_removed_eV)
            velocity_new = velocity
            v[:nlocal][newton_idx] = velocity_new
            te_max = float(grid.temperature_K.max())
            if te_max > 20000.0:
                raise RuntimeError(f"冷却 TTM 温度超过表格上限：{te_max:.6g} K")
            kinetic = 0.5 * M_VV2E_METAL * np.sum(masses[:, None] * velocity_new**2)
            temp_ag = float(2.0 * kinetic / (3.0 * len(newton_idx) * KB_EV_K))
            te_bottom = float(np.mean(grid.temperature_K[sink_mask]))
            combined = e_md + float(applied.sum()) + float(grid.electron_energy_eV.sum())
            residual = combined - (e_md0 + e_electron0) - source_total + sink_total
            if offset % 20 == 0 or offset == nsteps:
                writer.writerow([
                    step, f"{step * dt_ps:.6f}", f"{temp_ag:.8f}", f"{te_max:.8f}",
                    f"{te_bottom:.8f}", f"{fmax:.8f}", f"{source_total:.10f}",
                    f"{transfer_total:.10f}", f"{sink_total:.10f}", f"{residual:.12g}",
                ])
                stream.flush()
                print(f"step={step} t={step*dt_ps:.5f}ps T_Ag={temp_ag:.2f}K "
                      f"Te_max={te_max:.2f}K Te_bottom={te_bottom:.2f}K "
                      f"Fmax={fmax:.4f}eV/A sink_removed={sink_total:.3f}eV "
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
            bottom_sink_removed_eV=np.asarray(sink_total),
            initial_md_energy_eV=np.asarray(e_md0), initial_electron_energy_eV=np.asarray(e_electron0),
            ttm_time_ps=np.asarray(grid.time_ps), grid_temperature_K=grid.temperature_K,
            fraction_ag=grid.f_ag, source_power_eV_ps=grid.source_power,
            grid_spacing_A=grid.spacing_A, atom_ids_ttm=atom_ids_ttm,
            atom_types_ttm=atom_types_ttm, map_cell_ids=atom_map.cell_ids,
            map_weights=atom_map.weights, bottom_ag_ids=bottom_ag_ids,
            force_max_seen_eV_A=np.asarray(force_max_seen),
        )
    next_meta = {
        "status": "COMPLETED",
        "cooling_mode": "source off; 300 K ideal bottom Dirichlet sink on Ag TTM electron cells",
        "bottom_sink_z_max_A": BOTTOM_SINK_Z_MAX_A,
        "bottom_sink_temperature_K": BOTTOM_SINK_TEMPERATURE_K,
        "coupling_substeps_per_md_step": coupling_substeps,
        "bottom_sink_removed_total_eV": sink_total,
        "checkpoint_time_ps": final_step * dt_ps,
        "completed_step": final_step,
        "dt_ps": dt_ps,
        "atom_count": int(lmp.get_natoms()),
        "ttm_atom_count": int(len(atom_ids_ttm)),
        "grid_shape": list(grid.temperature_K.shape),
        "grid_spacing_A": [float(q) for q in grid.spacing_A],
        "ttm_time_ps": float(grid.time_ps),
        "source_total_eV": source_total,
        "electron_to_lattice_total_eV": float(transfer_total),
        "initial_md_energy_eV": float(e_md0),
        "initial_electron_energy_eV": float(e_electron0),
        "md_energy_at_checkpoint_eV": float(lmp.get_thermo("etotal")),
        "mobile_ag_force_max_seen_eV_A": float(force_max_seen),
        "input_sha256": meta.get("input_sha256", {}),
        "lammps_restart_sha256": sha256(next_restart),
        "ttm_state_sha256": sha256(new_state_path),
        "atomic_frame": frame_name,
        "restart_note": "冷却检查点保留完整 LAMMPS/TTM 原子和电子态，可继续关闭热源冷却。",
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
        "bottom_sink_removed_total_eV": sink_total,
        "energy_ledger_residual_eV": residual,
        "mobile_ag_force_max_seen_eV_A": force_max_seen,
        "bottom_sink_z_max_A": BOTTOM_SINK_Z_MAX_A,
        "bottom_sink_temperature_K": BOTTOM_SINK_TEMPERATURE_K,
        "runtime_s": time.monotonic() - started,
        "restart_file": str(next_restart),
        "ttm_state_file": str(new_state_path),
    }
    (OUTPUT_DIR / "cooling_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("FINAL " + json.dumps(result, ensure_ascii=False), flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
