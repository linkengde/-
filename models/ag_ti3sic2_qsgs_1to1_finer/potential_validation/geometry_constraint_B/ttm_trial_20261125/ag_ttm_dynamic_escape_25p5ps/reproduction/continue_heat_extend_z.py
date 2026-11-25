#!/usr/bin/env python3
"""从检查点继续加热，并允许在独立副本中扩大纯真空 z 上边界。"""

import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent.parent
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
        # 只扩展空真空，不移动原子、不改势能或 TTM 材料网格；默认不改变边界。
        *([f"change_box all z final 0.0 {float(os.environ['AG_TTM_ZHI_A']):.6f} units box"]
          if os.environ.get("AG_TTM_ZHI_A") else []),
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

    # 显式解除已经逸出凝聚相的 Ag 原子与其初始电子网格单元的耦合。
    # 仅接受高于材料表面、具有向外速度且远离其他 Ag 的原子，避免误解除表面热振动原子。
    detach_text = os.environ.get("AG_TTM_DETACHED_IDS", "").strip()
    requested_detach_ids = sorted({int(item.strip()) for item in detach_text.split(",") if item.strip()})
    inherited_detached_ids = [int(q) for q in meta.get("detached_atom_ids_from_ttm", [])]
    detached_ids = sorted(set(inherited_detached_ids) | set(requested_detach_ids))
    if requested_detach_ids:
        missing_detach = sorted(set(requested_detach_ids) - set(int(q) for q in atom_ids_ttm))
        if missing_detach:
            raise RuntimeError(f"请求解除的 ID 不在 TTM 映射中：{missing_detach}")
        min_z = float(os.environ.get("AG_TTM_DETACH_MIN_Z_A", "160.0"))
        min_nearest = float(os.environ.get("AG_TTM_DETACH_MIN_NEAREST_AG_A", "8.0"))
        ag_idx = np.flatnonzero(types[:nlocal] == 4)
        lx, ly = 200.0, 200.0
        for atom_id in requested_detach_ids:
            idx = id_to_index[atom_id]
            if types[idx] != 4:
                raise RuntimeError(f"解除映射的 ID {atom_id} 不是 Ag")
            if float(x[idx, 2]) <= min_z or float(v[idx, 2]) <= 0.0:
                raise RuntimeError(f"ID {atom_id} 未满足高于表面且向外运动的条件")
            others = ag_idx[ag_idx != idx]
            dx = np.abs(np.mod(x[others, 0], lx) - np.mod(x[idx, 0], lx))
            dy = np.abs(np.mod(x[others, 1], ly) - np.mod(x[idx, 1], ly))
            dx = np.minimum(dx, lx - dx)
            dy = np.minimum(dy, ly - dy)
            dz = x[others, 2] - x[idx, 2]
            nearest = float(np.sqrt(dx * dx + dy * dy + dz * dz).min())
            if nearest <= min_nearest:
                raise RuntimeError(
                    f"ID {atom_id} 最近 Ag 距离 {nearest:.3f} A，未达到逸出解除阈值"
                )
            print(f"TTM_DETACH id={atom_id} z={float(x[idx,2]):.3f}A "
                  f"vz={float(v[idx,2]):.4f}A/ps nearest_Ag={nearest:.3f}A", flush=True)
        keep = ~np.isin(atom_ids_ttm, np.asarray(requested_detach_ids, dtype=np.int64))
        atom_ids_ttm = atom_ids_ttm[keep]
        atom_types_ttm = atom_types_ttm[keep]
        map_cell_ids = map_cell_ids[keep]
        map_weights = map_weights[keep]
        if np.any(np.isin(atom_ids_ttm, np.asarray(detached_ids, dtype=np.int64))):
            raise RuntimeError("逸出原子仍残留在 TTM 映射中")

    # `newton_idx` 仅对应仍与电子网格耦合的 Ag；脱离原映射的原子继续按 NVE 自由运动。
    newton_idx = np.asarray([id_to_index[int(atom_id)] for atom_id in atom_ids_ttm], dtype=int)
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
    # 自动脱耦判据：高于初始材料表面 10 Å、vz>0，且最近 Ag 距离大于 8 Å。
    geometry = json.loads((HERE.parent / "geometry_report.json").read_text(encoding="utf-8"))
    surface_z_A = float(geometry["material_dimensions_A"][2])
    escape_z_A = float(os.environ.get("AG_TTM_ESCAPE_Z_A", str(surface_z_A + 10.0)))
    escape_nearest_A = float(os.environ.get("AG_TTM_ESCAPE_NEAREST_AG_A", "8.0"))
    detached_events = []
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
            newly_detached = []
            # 每 0.01 ps 搜索一次真正离开 Ag 凝聚相且向外运动的 TTM 原子。
            if offset % 20 == 0:
                candidates = newton_idx[(x[newton_idx, 2] > escape_z_A) & (v[newton_idx, 2] > 0.0)]
                if len(candidates):
                    ag_idx = np.flatnonzero(types[:nlocal] == 4)
                    lx, ly = 200.0, 200.0
                    ag_pos = np.column_stack((np.mod(x[ag_idx, 0], lx),
                                              np.mod(x[ag_idx, 1], ly),
                                              x[ag_idx, 2] + 100.0))
                    spatial = cKDTree(ag_pos, boxsize=[lx, ly, 600.0])
                    query_pos = np.column_stack((np.mod(x[candidates, 0], lx),
                                                 np.mod(x[candidates, 1], ly),
                                                 x[candidates, 2] + 100.0))
                    nearest = spatial.query(query_pos, k=2)[0][:, 1]
                    newly_detached = [int(ids[idx]) for idx, dist in zip(candidates, nearest)
                                      if float(dist) > escape_nearest_A]
                    if newly_detached:
                        keep = ~np.isin(atom_ids_ttm, np.asarray(newly_detached, dtype=np.int64))
                        if int(keep.sum()) == len(atom_ids_ttm):
                            raise RuntimeError("自动逸出筛查未能从 TTM 映射移除候选原子")
                        for idx, dist in zip(candidates, nearest):
                            if int(ids[idx]) in newly_detached:
                                event = {"step": step, "time_ps": step * dt_ps,
                                         "atom_id": int(ids[idx]), "z_A": float(x[idx, 2]),
                                         "vz_A_ps": float(v[idx, 2]),
                                         "nearest_Ag_A": float(dist)}
                                detached_events.append(event)
                                print("TTM_AUTO_DETACH " + json.dumps(event, ensure_ascii=False), flush=True)
                        atom_ids_ttm = atom_ids_ttm[keep]
                        atom_types_ttm = atom_types_ttm[keep]
                        map_cell_ids = map_cell_ids[keep]
                        map_weights = map_weights[keep]
                        newton_idx = np.asarray([id_to_index[int(atom_id)] for atom_id in atom_ids_ttm], dtype=int)
                        atom_map = AtomGridMap(map_cell_ids, map_weights, tuple(int(q) for q in fraction_ag.shape))
                        masses = np.full(len(newton_idx), 107.8682)
                        detached_ids = sorted(set(detached_ids) | set(newly_detached))
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
        "restart_note": "LAMMPS 原子状态与 TTM 电子网格状态成对保存；已解除映射的逸出原子继续按 NVE 运动。",
        "detached_atom_ids_from_ttm": detached_ids,
        "detachment_events_this_segment": detached_events,
        "escape_detection": {"surface_z_A": surface_z_A, "z_threshold_A": escape_z_A,
                             "minimum_nearest_ag_A": escape_nearest_A,
                             "checked_every_steps": 20},
        "continuation_script_sha256": sha256(Path(__file__).resolve()),
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
        "detached_atom_ids_from_ttm": detached_ids,
        "detachment_events_this_segment": detached_events,
        "escape_detection": {"surface_z_A": surface_z_A, "z_threshold_A": escape_z_A,
                             "minimum_nearest_ag_A": escape_nearest_A,
                             "checked_every_steps": 20},
        "ttm_atom_count": int(len(atom_ids_ttm)),
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
