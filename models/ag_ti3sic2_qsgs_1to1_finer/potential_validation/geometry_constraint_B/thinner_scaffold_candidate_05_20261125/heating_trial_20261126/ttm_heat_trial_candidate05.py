#!/usr/bin/env python3
"""QSGS 更细骨架候选：刚性骨架 + Ag EAM + TTM 参考耦合探索试跑。

此处的 Ti/Si/C 骨架只作刚性几何支撑；Ag–骨架为排斥 WCA 理想壁。
Python TTM 是仓库数值参考耦合器，不是生产级原生 LAMMPS TTM fix。
"""

import csv
import hashlib
import json
import math
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
B_DIR = ROOT / "potential_validation/geometry_constraint_B"
TRIAL_DIR = Path(os.environ.get(
    "AG_TTM_TRIAL_OUTPUT_DIR",
    str(B_DIR / "thinner_scaffold_candidate_05_20261125/heating_trial_20261126/ttm_initial_0p5ps"),
))
DATA = Path(os.environ.get(
    "AG_TTM_MODEL_DATA",
    str(B_DIR / "thinner_scaffold_candidate_05_20261125/ag_ti3sic2_qsgs.data"),
))
sys.path.insert(0, str(ROOT))

from mpi4py import MPI  # noqa: E402
from lammps import lammps  # noqa: E402
from ag_ttm_reference import (  # noqa: E402
    AgParameters,
    AgTTMGrid,
    build_ag_grid_fields,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
    KB_EV_K,
    M_VV2E_METAL,
)


def save_json(path, payload):
    """原子替换写报告，防止异常退出留下截断 JSON。"""
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    temp.replace(path)


def write_frame(lmp, name):
    """保存 OVITO 可读的完整压缩原子帧，含位置、速度和受力。"""
    path = TRIAL_DIR / f"{name}.dump.gz"
    lmp.command(
        f"write_dump all custom/gz {path} id type x y z vx vy vz fx fy fz modify sort id"
    )
    return path


def save_restart_checkpoint(lmp, grid, fields, ids, types, newton_idx, step,
                     source_total, transfer_total, e_md0, e_electron0,
                     dt_ps, source_scale, force_max_seen, bottom_ag_ids, checkpoint_dir):
    """保存可续算的 LAMMPS 原子态与完整 TTM 电子网格状态。"""
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    restart_path = checkpoint_dir / "lammps.restart"
    lmp.command(f"write_restart {restart_path}")

    # 原子映射按 atom ID 保存，续算时不依赖 LAMMPS 内部数组顺序。
    atom_ids_ttm = np.asarray(ids[:int(lmp.extract_setting("nlocal"))][newton_idx], dtype=np.int64)
    atom_types_ttm = np.asarray(types[:int(lmp.extract_setting("nlocal"))][newton_idx], dtype=np.int32)
    state_path = checkpoint_dir / "ttm_state.npz"
    tmp_path = checkpoint_dir / "ttm_state.npz.tmp"
    with tmp_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            completed_step=np.asarray(step, dtype=np.int64),
            dt_ps=np.asarray(dt_ps, dtype=np.float64),
            source_scale=np.asarray(source_scale, dtype=np.float64),
            source_total_eV=np.asarray(source_total, dtype=np.float64),
            electron_to_lattice_total_eV=np.asarray(transfer_total, dtype=np.float64),
            initial_md_energy_eV=np.asarray(e_md0, dtype=np.float64),
            initial_electron_energy_eV=np.asarray(e_electron0, dtype=np.float64),
            ttm_time_ps=np.asarray(grid.time_ps, dtype=np.float64),
            grid_temperature_K=grid.temperature_K,
            fraction_ag=grid.f_ag,
            source_power_eV_ps=grid.source_power,
            grid_spacing_A=grid.spacing_A,
            atom_ids_ttm=atom_ids_ttm,
            atom_types_ttm=atom_types_ttm,
            bottom_ag_ids=np.asarray(bottom_ag_ids, dtype=np.int64),
            map_cell_ids=fields.atom_map.cell_ids,
            map_weights=fields.atom_map.weights,
            force_max_seen_eV_A=np.asarray(force_max_seen, dtype=np.float64),
        )
    tmp_path.replace(state_path)

    def sha256(path):
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    provenance_inputs = [
        DATA,
        B_DIR / "potential.geometry_wca_fixed_skeleton.inc",
        ROOT / "potentials/Ag_u3.eam",
        ROOT / "ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv",
    ]

    metadata = {
        "checkpoint_time_ps": float(step * dt_ps),
        "completed_step": int(step),
        "atom_count": int(lmp.get_natoms()),
        "md_energy_at_checkpoint_eV": float(lmp.get_thermo("etotal")),
        "ttm_atom_count": int(len(atom_ids_ttm)),
        "grid_shape": list(grid.temperature_K.shape),
        "grid_spacing_A": [float(v) for v in grid.spacing_A],
        "ttm_time_ps": float(grid.time_ps),
        "source_total_eV": float(source_total),
        "electron_to_lattice_total_eV": float(transfer_total),
        "mobile_ag_force_max_seen_eV_A": float(force_max_seen),
        "initial_md_energy_eV": float(e_md0),
        "initial_electron_energy_eV": float(e_electron0),
        "lammps_restart_sha256": sha256(restart_path),
        "ttm_state_sha256": sha256(state_path),
        "atomic_frame": f"frame_{step * dt_ps:06.3f}ps_step{step:05d}.dump.gz",
        "input_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in provenance_inputs
        },
        "restart_note": "The binary LAMMPS restart contains atom coordinates, velocities, IDs, types and timestep; ttm_state.npz contains the electron grid and cumulative TTM accounting.",
    }
    save_json(checkpoint_dir / "CHECKPOINT.json", metadata)


def main():
    os.chdir(ROOT)
    nsteps = int(os.environ.get("AG_TTM_TRIAL_STEPS", "2000"))
    source_scale = float(os.environ.get("AG_TTM_SOURCE_SCALE", "1.0"))
    dt_ps = 0.0005  # 0.5 fs
    force_abort_eVA = 50.0  # 数值安全停止线，不是物理验收值
    if nsteps <= 0 or not math.isfinite(source_scale) or source_scale < 0.0:
        raise ValueError("模拟步数或热源倍率不合法")

    geometry = json.loads((TRIAL_DIR / "geometry_report.json").read_text())
    lx, ly, lz_box = map(float, geometry["box_A"])
    z_material_lo, z_material_hi = map(float, geometry["material_z_A"])

    # 1. 装入归档原始几何；所有 Ti/Si/C 固定，底部 Ag 薄层固定。
    #    Ag-Ag 使用 Foiles Ag_u3 EAM；Ag–骨架使用 WCA 纯排斥几何边界。
    lmp = lammps(cmdargs=["-screen", "none", "-log", "none", "-sf", "omp"], comm=MPI.COMM_WORLD)
    commands = [
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        f"read_data {DATA}",
        f"include {B_DIR}/potential.geometry_wca_fixed_skeleton.inc",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "region r_bottom block INF INF INF INF 0.0 20.0 units box",
        "group bottom region r_bottom",
        "group Ag type 4",
        "group skeleton type 1 2 3",
        "group bottom_Ag intersect Ag bottom",
        "group mobile_Ag subtract Ag bottom_Ag",
        "fix freeze_skeleton skeleton setforce 0.0 0.0 0.0",
        "fix freeze_bottom_Ag bottom_Ag setforce 0.0 0.0 0.0",
        "mass 4 107.8682",
        "velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian",
        "timestep 0.0005",
        "fix integrate mobile_Ag nve",
        "thermo 1",
        "thermo_style custom step pe ke etotal temp",
        "thermo_modify lost error flush yes",
        "run 0",
    ]
    for command in commands:
        lmp.command(command)

    nlocal = int(lmp.extract_setting("nlocal"))
    x = lmp.numpy.extract_atom("x")
    v = lmp.numpy.extract_atom("v")
    f = lmp.numpy.extract_atom("f")
    types = lmp.numpy.extract_atom("type")
    ids = lmp.numpy.extract_atom("id")
    mobile_idx = np.flatnonzero((types[:nlocal] == 4) & (x[:nlocal, 2] >= 20.0))
    bottom_ag_ids = np.asarray(
        ids[:nlocal][(types[:nlocal] == 4) & (x[:nlocal, 2] < 20.0)], dtype=np.int64
    )
    newton_idx = np.flatnonzero(
        (types[:nlocal] == 4)
        & (x[:nlocal, 2] >= 40.0)
        & (x[:nlocal, 2] <= z_material_hi)
    )
    if len(mobile_idx) < 1000 or len(newton_idx) < 1000:
        raise RuntimeError(
            f"移动 Ag 或 TTM 层数量异常：mobile={len(mobile_idx)}, TTM={len(newton_idx)}"
        )

    mobile_force0 = float(np.linalg.norm(np.asarray(f[:nlocal][mobile_idx]), axis=1).max())
    if mobile_force0 > force_abort_eVA:
        raise RuntimeError(f"初始移动 Ag Fmax={mobile_force0:.6g} eV/A，超过安全停止线")

    # 2. 只在实际材料厚度 0–150 Å 建电子网格；真空区不进入 TTM。
    xyz_newton = np.array(x[:nlocal][newton_idx], copy=True)
    masses = np.full(len(newton_idx), 107.8682)
    shape = np.ceil(np.array([lx, ly, z_material_hi - z_material_lo]) / 10.0).astype(int)
    spacing = np.array([lx, ly, z_material_hi - z_material_lo]) / shape
    r2 = (xyz_newton[:, 0] - lx / 2.0) ** 2 + (xyz_newton[:, 1] - ly / 2.0) ** 2
    source_per_atom = (5.0 * source_scale / 20.0) * np.exp(
        -r2 / (30.0**2) - (z_material_hi - xyz_newton[:, 2]) / 20.0
    )
    fields = build_ag_grid_fields(
        xyz_newton,
        source_per_atom,
        tuple(shape),
        tuple(spacing),
        4.09**3 / 4.0,
    )
    params = AgParameters.from_csv(
        ROOT / "ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv",
        allow_low_temperature_extension=True,
    )
    grid = AgTTMGrid(
        params,
        fields.fraction_ag,
        tuple(spacing),
        fields.source_power_eV_ps,
        initial_temperature_K=300.0,
    )

    e_md0 = float(lmp.get_thermo("etotal"))
    e_electron0 = float(grid.electron_energy_eV.sum())
    source_total = 0.0
    transfer_total = 0.0
    force_max_seen = mobile_force0
    electron_temperature_max = float(grid.temperature_K.max())
    active_cells = fields.fraction_ag.ravel() > 0.0
    hot_cells = fields.source_power_eV_ps.ravel() >= 0.01 * fields.source_power_eV_ps.max()
    source_weights = fields.source_power_eV_ps.ravel()[hot_cells]
    started = time.monotonic()
    status = "COMPLETED"
    completed_step = 0

    metrics_file = TRIAL_DIR / "step_metrics.csv"
    with metrics_file.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "step", "time_ps", "T_newton_Ag_K", "Tl_source_weighted_K", "Tl_active_max_K",
            "Te_active_min_K", "Te_max_K", "Fmax_mobile_Ag_eV_A", "source_total_eV",
            "electron_to_lattice_eV", "energy_ledger_residual_eV", "ttm_substeps",
        ])
        print(
            f"INIT atoms={lmp.get_natoms()} Ag_mobile={len(mobile_idx)} Ag_TTM={len(newton_idx)} "
            f"grid={tuple(int(q) for q in shape)} source_power={fields.source_power_eV_ps.sum():.9f}eV/ps "
            f"Fmax={mobile_force0:.6f}eV/A", flush=True
        )
        write_frame(lmp, "frame_00000_step0")

        for step in range(1, nsteps + 1):
            lmp.command("run 1 pre no post no")
            completed_step = step
            e_md = float(lmp.get_thermo("etotal"))
            velocity = np.array(v[:nlocal][newton_idx], copy=True)
            fmax = float(np.linalg.norm(np.asarray(f[:nlocal][mobile_idx]), axis=1).max())
            force_max_seen = max(force_max_seen, fmax)
            if not math.isfinite(fmax) or fmax > force_abort_eVA:
                status = f"FORCE_SAFETY_STOP_FMAX_{fmax:.6g}_EV_A"
                write_frame(lmp, f"safety_stop_step{step:05d}")
                print(status, flush=True)
                break

            tl = lattice_temperature_by_cell(
                velocity, masses, fields.atom_map, smoothing_radius_cells=1
            )
            increment = grid.advance(dt_ps, tl, heat_on=True)
            velocity_new, applied = rescale_lattice_kinetic(
                velocity, masses, fields.atom_map, increment.lattice_energy_per_cell_eV
            )
            v[:nlocal][newton_idx] = velocity_new
            source_total += increment.source_to_electron_eV
            transfer_total += float(applied.sum())

            tl_flat = lattice_temperature_by_cell(
                velocity_new, masses, fields.atom_map, smoothing_radius_cells=1
            ).ravel()
            te_active = grid.temperature_K.ravel()[active_cells]
            te_min = float(te_active.min())
            te_max = float(grid.temperature_K.max())
            electron_temperature_max = max(electron_temperature_max, te_max)
            if te_max > 20000.0:
                status = f"TTM_TABLE_RANGE_STOP_TE_{te_max:.6g}_K"
                write_frame(lmp, f"safety_stop_step{step:05d}")
                print(status, flush=True)
                break

            tl_source = float(np.average(tl_flat[hot_cells], weights=source_weights))
            tl_max = float(np.nanmax(tl_flat[active_cells]))
            kinetic = 0.5 * M_VV2E_METAL * np.sum(masses[:, None] * velocity_new**2)
            temp_ag = float(2.0 * kinetic / (3.0 * len(newton_idx) * KB_EV_K))
            combined = e_md + float(applied.sum()) + float(grid.electron_energy_eV.sum())
            residual = combined - (e_md0 + e_electron0) - source_total

            if step % 20 == 0 or step == 1 or step == nsteps:
                writer.writerow([
                    step, f"{step * dt_ps:.6f}", f"{temp_ag:.8f}", f"{tl_source:.8f}",
                    f"{tl_max:.8f}", f"{te_min:.8f}", f"{te_max:.8f}", f"{fmax:.8f}",
                    f"{source_total:.10f}", f"{transfer_total:.10f}", f"{residual:.12g}",
                    increment.substeps,
                ])
                stream.flush()
                elapsed = time.monotonic() - started
                peak_rss_gib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2)
                print(
                    f"step={step} t={step*dt_ps:.5f}ps T_Ag={temp_ag:.2f}K "
                    f"Tl_source={tl_source:.2f}K Tl_max={tl_max:.2f}K "
                    f"Te_minmax={te_min:.2f},{te_max:.2f}K Fmax={fmax:.4f}eV/A "
                    f"source={source_total:.3f}eV e2l={transfer_total:.3f}eV "
                    f"ledger={residual:.6g}eV substeps={increment.substeps} "
                    f"elapsed={elapsed:.1f}s peak_RSS={peak_rss_gib:.3f}GiB",
                    flush=True,
                )

            if step % 500 == 0 or step == nsteps:
                lmp.command("run 0 pre no post no")
                write_frame(lmp, f"frame_{step*dt_ps:06.3f}ps_step{step:05d}")

    lmp.command("run 0 pre no post no")
    e_mdf = float(lmp.get_thermo("etotal"))
    e_electronf = float(grid.electron_energy_eV.sum())
    residual_final = (e_mdf + e_electronf) - (e_md0 + e_electron0) - source_total
    v_final = np.array(v[:nlocal][newton_idx], copy=True)
    ke_final = 0.5 * M_VV2E_METAL * np.sum(masses[:, None] * v_final**2)
    elapsed = time.monotonic() - started
    summary = {
        "status": status,
        "model_atoms": int(lmp.get_natoms()),
        "Ag_atoms": int(np.count_nonzero(types[:nlocal] == 4)),
        "Ti3SiC2_atoms": int(np.count_nonzero(types[:nlocal] != 4)),
        "geometry_data": str(DATA),
        "Ag_mobile_atoms": int(len(mobile_idx)),
        "Ag_TTM_atoms": int(len(newton_idx)),
        "box_A": [lx, ly, lz_box],
        "material_z_A": [z_material_lo, z_material_hi],
        "skeleton_handling": "all Ti/Si/C atoms fixed as an ideal rigid support",
        "bottom_boundary": "Ag with initial z<20 A fixed; no separate thermostat during this source-on trial",
        "Ag_Ag_potential": "Foiles Ag_u3.eam",
        "Ag_skeleton_boundary": "repulsive WCA geometric support; not a fitted physical interface",
        "TTM_implementation": "repository Python reference operator coupled through LAMMPS Python API; not native LAMMPS fix ttm",
        "duration_requested_ps": nsteps * dt_ps,
        "duration_completed_ps": completed_step * dt_ps,
        "timestep_fs": dt_ps * 1000.0,
        "source_peak_eV_ps_A2": 5.0 * source_scale,
        "source_scale": source_scale,
        "source_profile": "XY circular Gaussian (R=30 A) times z surface exponential (skin=20 A)",
        "electron_grid_shape": [int(q) for q in shape],
        "electron_grid_spacing_A": [float(q) for q in spacing],
        "source_total_eV": source_total,
        "electron_to_lattice_energy_eV": transfer_total,
        "combined_energy_residual_eV": residual_final,
        "newton_Ag_temperature_final_K": float(2.0 * ke_final / (3.0 * len(newton_idx) * KB_EV_K)),
        "mobile_Ag_Fmax_max_eV_A": force_max_seen,
        "electron_temperature_max_K": electron_temperature_max,
        "force_abort_threshold_eV_A": force_abort_eVA,
        "atom_to_grid_map": "fixed at starting coordinates for this short exploratory run",
        "runtime_s": elapsed,
        "process_peak_RSS_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2),
        "trajectory_frames": sorted(p.name for p in TRIAL_DIR.glob("frame_*.dump.gz")),
        "limitations": [
            "The source model was built as geometry only and has not passed force minimization or thermal equilibration.",
            "All Ti3SiC2 atoms are immobile, so the scaffold does not vibrationally conduct heat.",
            "The WCA wall provides only geometric exclusion and has no calibrated Ag/MAX-phase binding or adhesion.",
            "The TTM operator is a numerical reference prototype with deterministic mean-field energy handoff.",
            "The atom-to-grid map remains fixed as Ag moves; large structural changes invalidate a production interpretation.",
            "This run contains only the heat-on stage; it does not establish a stable melt pool, heat-off response, or resolidification.",
        ],
    }
    save_json(TRIAL_DIR / "trial_summary.json", summary)
    if status == "COMPLETED" and completed_step == nsteps:
        save_restart_checkpoint(
            lmp, grid, fields, ids, types, newton_idx, completed_step,
            source_total, transfer_total, e_md0, e_electron0, dt_ps,
            source_scale, force_max_seen, bottom_ag_ids,
            TRIAL_DIR / "checkpoint_latest",
        )
    print("FINAL " + json.dumps(summary, ensure_ascii=False), flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
