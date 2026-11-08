#!/usr/bin/env python3
"""Resume the Round-5 MTP + Ag-only TTM reference run from a saved state."""
from __future__ import annotations

import csv
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
REFERENCE_MODEL = MODEL.parent / "ag_ti3sic2_qsgs_1to1_finer"
BASE_RUN = MODEL / "runs" / "ttm_mtp_round5_dt0p1fs_0p05ps_20261107"
SOURCE_RUN = Path(os.environ.get("TTM_RESTART_FROM", str(BASE_RUN))).resolve()
PREVIOUS_SUMMARY = Path(os.environ.get("TTM_PREVIOUS_SUMMARY", str(SOURCE_RUN / "summary.json"))).resolve()
sys.path.insert(0, str(REFERENCE_MODEL))

from lammps import lammps
from ag_ttm_reference import (
    AgGridFields,
    AgParameters,
    AgTTMGrid,
    AtomGridMap,
    KB_EV_K,
    M_VV2E_METAL,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)

DT_PS = 0.0001
TARGET_PS = float(os.environ.get("TTM_TARGET_PS", "0.1"))
CHECKPOINT_INTERVAL = int(os.environ.get("TTM_CHECKPOINT_INTERVAL", "500"))
STEP_LIMIT = int(os.environ.get("TTM_STEP_LIMIT", "0"))
RESTART_TIME_PS = float(os.environ.get("TTM_RESTART_TIME_PS", "0.05"))
ATOMIC_MASSES = {1: 47.867, 2: 28.085, 3: 12.011, 4: 107.8682}
SEED = 190427


def read_initial_model(path: Path):
    with gzip.open(path, "rt") as stream:
        lines = stream.readlines()
    natoms = int(next(line.split()[0] for line in lines if "atoms" in line))
    start = next(i for i, line in enumerate(lines) if line.startswith("Atoms #")) + 2
    atoms = np.loadtxt(lines[start:start + natoms])
    atoms = atoms[np.argsort(atoms[:, 0])]
    ids = atoms[:, 0].astype(np.int64)
    types = atoms[:, 1].astype(np.int32)
    xyz = atoms[:, 2:5]
    if len(ids) != natoms or not np.array_equal(ids, np.arange(1, natoms + 1)):
        raise RuntimeError("Initial model atom IDs are not a contiguous 1..N sequence")
    return ids, types, xyz


def contiguous_ranges(values: np.ndarray) -> list[str]:
    values = np.asarray(values, dtype=np.int64)
    if values.size == 0:
        return []
    values = np.unique(values)
    starts = np.r_[0, np.flatnonzero(np.diff(values) > 1) + 1]
    ends = np.r_[starts[1:] - 1, len(values) - 1]
    result = []
    for start, end in zip(starts, ends):
        a, b = int(values[start]), int(values[end])
        result.append(str(a) if a == b else f"{a}:{b}")
    return result


def add_id_group(lmp, name: str, ids: np.ndarray) -> None:
    lmp.command(f"group {name} empty")
    tokens = contiguous_ranges(ids)
    for start in range(0, len(tokens), 400):
        lmp.command(f"group {name} id " + " ".join(tokens[start:start + 400]))


def kinetic_eV(velocity: np.ndarray, mass: np.ndarray) -> float:
    return float(0.5 * M_VV2E_METAL * np.sum(mass[:, None] * velocity**2))


def temperature_K(velocity: np.ndarray, mass: np.ndarray) -> float:
    if len(mass) == 0:
        return float("nan")
    return 2.0 * kinetic_eV(velocity, mass) / (3.0 * len(mass) * KB_EV_K)


def main() -> None:
    os.chdir(MODEL)
    report = json.loads((MODEL / "geometry_report.json").read_text())
    previous = json.loads(PREVIOUS_SUMMARY.read_text())
    restart_tag = f"{RESTART_TIME_PS:.4f}".rstrip("0").rstrip(".").replace(".", "p")
    restart_path = SOURCE_RUN / "checkpoints" / f"ttm_{restart_tag}ps.restart"
    state_path = SOURCE_RUN / "checkpoints" / f"electrons_{restart_tag}ps.npz"
    if not restart_path.is_file() or not state_path.is_file():
        raise FileNotFoundError("The 0.05 ps LAMMPS/TTM restart pair is incomplete")

    model_ids, model_type, model_xyz = read_initial_model(MODEL / "ag_ti3sic2_qsgs.data.gz")
    initial_z = model_xyz[:, 2]
    initial_fixed = initial_z < 20.0
    initial_bath = (initial_z >= 20.0) & (initial_z < 40.0)
    initial_newton = (initial_z >= 40.0) & (initial_z < 150.0)

    state = np.load(state_path)
    te_saved = np.array(state["electron_temperature_K"], copy=True)
    fraction_ag = np.array(state["fraction_ag"], copy=True)
    source_power = np.array(state["source_power_eV_ps"], copy=True)
    spacing = tuple(float(v) for v in state["spacing_A"])
    grid_shape = tuple(int(v) for v in te_saved.shape)
    ag_ids = np.array(state["atom_ids_newton_ag"], dtype=np.int64, copy=True)
    atom_map = AtomGridMap(
        np.array(state["atom_cell_ids"], dtype=np.int64, copy=True),
        np.array(state["atom_cell_weights"], dtype=float, copy=True),
        grid_shape,
    )
    fields = AgGridFields(
        fraction_ag,
        source_power,
        atom_map,
        int(previous["ag_grid_fill_diagnostics"]["raw_overfilled_cells"]),
        float(previous["ag_grid_fill_diagnostics"]["raw_max_fill_fraction"]),
        float(previous["ag_grid_fill_diagnostics"]["redistributed_ag_volume_A3"]),
    )
    parameters = AgParameters.from_csv(
        REFERENCE_MODEL / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv",
        allow_low_temperature_extension=True,
    )
    grid = AgTTMGrid(
        parameters,
        fraction_ag,
        spacing,
        source_power,
        initial_temperature_K=300.0,
    )
    grid.temperature_K[:] = te_saved
    grid.time_ps = float(state["time_ps"])
    previous_source = float(state["cumulative_source_eV"])
    previous_transfer = float(state["cumulative_electron_to_lattice_eV"])
    previous_grid_residual = float(state["cumulative_ttm_residual_eV"])

    lx, ly, lz = (float(v) for v in report["material_dimensions_A"])
    MTP_DIR = (
        REFERENCE_MODEL
        / "potential_validation"
        / "mtp_thermal_acquisition_pilot_20261106"
        / "mtp4_ag108_ti3sic2_tscbop_thermal_interface_weighted_round5"
    )
    MTP_PATH = MTP_DIR / "MTP1m_AgTiSiC_thermal_interface_weighted_round5_DEVELOPMENT_ONLY.mtp"
    MTP_INI = HERE / "mlip.ini"
    MTP_INI.write_text(f"mtp-filename {MTP_PATH}\nselect FALSE\n")
    OUT_LOG = HERE / "lammps.log"

    lmp = lammps(name="mtp", cmdargs=["-screen", "none", "-log", str(OUT_LOG)])
    for command in (
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        f"read_restart {restart_path}",
        f"pair_style mlip {MTP_INI}",
        "pair_coeff * *",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "timestep 0.0001",
        "region r_fixed block INF INF INF INF 0.0 20.0 units box",
        "region r_bath block INF INF INF INF 20.0 40.0 units box",
        "region r_newton block INF INF INF INF 40.0 150.0 units box",
        "region r_sputter block INF INF INF INF 150.0 200.0 units box",
        "region r_evaporation block INF INF INF INF 200.0 250.0 units box",
        "group fixed_region region r_fixed",
        "group bath_region region r_bath",
        "group newton_region region r_newton",
        "run 0",
    ):
        lmp.command(command)

    nlocal = int(lmp.extract_setting("nlocal"))
    if nlocal != int(report["atoms"]):
        raise RuntimeError(f"Expected all atoms on one rank, got {nlocal} local of {lmp.get_natoms()}")
    x = lmp.numpy.extract_atom("x")
    v = lmp.numpy.extract_atom("v")
    ids_ptr = lmp.numpy.extract_atom("id")
    type_ptr = lmp.numpy.extract_atom("type")
    ids = np.asarray(ids_ptr[:nlocal], dtype=np.int64).copy()
    atom_type = np.asarray(type_ptr[:nlocal], dtype=np.int32).copy()
    xyz = np.array(x[:nlocal], copy=True)
    velocity = np.array(v[:nlocal], copy=True)
    masses_all = np.array([ATOMIC_MASSES[int(t)] for t in atom_type])
    if len(np.unique(ids)) != nlocal:
        raise RuntimeError("Restart contains duplicate atom IDs")
    initial_index = ids - 1
    if np.any(initial_index < 0) or np.any(initial_index >= len(model_ids)):
        raise RuntimeError("Restart atom IDs do not map to the original model")
    if not np.array_equal(atom_type, model_type[initial_index]):
        raise RuntimeError("Restart atom types do not match the original model")
    current_z_region = (
        (xyz[:, 2] < 20.0),
        (xyz[:, 2] >= 20.0) & (xyz[:, 2] < 40.0),
        (xyz[:, 2] >= 40.0) & (xyz[:, 2] < 150.0),
    )
    original_region = (initial_fixed[initial_index], initial_bath[initial_index], initial_newton[initial_index])
    missing_ids = [ids[was & ~now] for was, now in zip(original_region, current_z_region)]
    wrong_ids = [ids[now & ~was] for was, now in zip(original_region, current_z_region)]

    for group_name, values in zip(("missing_fixed", "wrong_fixed", "missing_bath", "wrong_bath", "missing_newton", "wrong_newton"),
                                  (missing_ids[0], wrong_ids[0], missing_ids[1], wrong_ids[1], missing_ids[2], wrong_ids[2])):
        if len(values):
            add_id_group(lmp, group_name, values)
        else:
            lmp.command(f"group {group_name} empty")
    lmp.command("group fixed_corrected union fixed_region missing_fixed")
    lmp.command("group fixed subtract fixed_corrected wrong_fixed")
    lmp.command("group bath_cross union missing_bath wrong_bath")
    lmp.command("group bath_plus union bath_region bath_cross")
    lmp.command("group bath subtract bath_plus wrong_bath")
    lmp.command("group newton_cross union missing_newton wrong_newton")
    lmp.command("group newton_plus union newton_region newton_cross")
    lmp.command("group newton subtract newton_plus wrong_newton")
    lmp.command("group Ag type 4")
    lmp.command("group newton_Ag intersect newton Ag")
    lmp.command("group newton_skeleton type 1 2 3")
    lmp.command("group mobile union bath newton")
    lmp.command("group sputter dynamic all region r_sputter every 10")
    lmp.command("group evaporation dynamic all region r_evaporation every 10")
    lmp.command("fix freeze fixed setforce 0.0 0.0 0.0")
    lmp.command("fix bath_control bath nvt temp 300.0 300.0 0.1")
    lmp.command("fix integrate_newton newton nve")
    lmp.command("compute T_bath bath temp")
    lmp.command("compute T_newton newton temp")
    lmp.command("thermo 20")
    lmp.command("thermo_style custom step time atoms temp c_T_bath c_T_newton pe ke etotal ecouple econserve fmax")
    lmp.command("thermo_modify flush yes")
    lmp.command("run 0")

    id_to_index = {int(atom_id): index for index, atom_id in enumerate(ids)}
    try:
        ag_indices = np.array([id_to_index[int(atom_id)] for atom_id in ag_ids], dtype=np.int64)
    except KeyError as exc:
        raise RuntimeError(f"Saved Ag-to-TTM mapping references missing atom {exc}") from exc
    if not np.all(atom_type[ag_indices] == 4):
        raise RuntimeError("Saved TTM map contains a non-Ag atom")
    group_masks = {
        "fixed": initial_fixed[initial_index],
        "bath": initial_bath[initial_index],
        "newton": initial_newton[initial_index],
        "newton_ag": (atom_type == 4) & initial_newton[initial_index],
        "newton_skeleton": np.isin(atom_type, (1, 2, 3)) & initial_newton[initial_index],
    }
    if int(group_masks["newton_ag"].sum()) != len(ag_ids):
        raise RuntimeError("Static Newton-layer Ag set differs from saved TTM map")

    start_step = int(round(float(state["time_ps"]) / DT_PS))
    actual_start_step = int(lmp.get_thermo("step"))
    if actual_start_step != start_step:
        raise RuntimeError(f"Restart step {actual_start_step} does not match TTM state step {start_step}")
    restart_atom_energy = float(lmp.get_thermo("pe")) + kinetic_eV(np.array(v[:nlocal], copy=True), masses_all)
    restart_electron_energy = float(grid.electron_energy_eV.sum())
    if not np.isclose(restart_atom_energy, float(previous["final_atom_energy_eV"]), rtol=0.0, atol=2.0e-4):
        raise RuntimeError(
            f"LAMMPS restart atom energy {restart_atom_energy:.9f} eV does not match saved "
            f"checkpoint {previous['final_atom_energy_eV']:.9f} eV"
        )
    if not np.isclose(restart_electron_energy, float(previous["final_electron_energy_eV"]), rtol=0.0, atol=2.0e-7):
        raise RuntimeError(
            f"Restored TTM energy {restart_electron_energy:.9f} eV does not match saved "
            f"checkpoint {previous['final_electron_energy_eV']:.9f} eV"
        )
    target_step = int(round(TARGET_PS / DT_PS))
    requested_steps = max(0, target_step - start_step)
    steps_to_run = min(requested_steps, STEP_LIMIT) if STEP_LIMIT > 0 else requested_steps
    if steps_to_run <= 0:
        raise ValueError("Target duration must be later than the saved restart")

    ag_masses = np.full(len(ag_ids), ATOMIC_MASSES[4])
    segment_atom_start = float(lmp.get_thermo("pe")) + kinetic_eV(velocity, masses_all)
    segment_electron_start = float(grid.electron_energy_eV.sum())
    segment_ecouple_start = float(lmp.get_thermo("ecouple"))
    segment_source = 0.0
    segment_transfer = 0.0
    segment_grid_residual = 0.0
    segment_substeps = 0
    max_substeps = 0

    old_total_ecouple = float(previous["thermostat_ecouple_change_eV"])
    old_total_residual = float(previous["source_corrected_combined_energy_residual_eV"])
    old_total_source = previous_source
    old_total_transfer = previous_transfer
    old_total_grid_residual = previous_grid_residual
    old_total_substeps = previous.get("electron_substeps_total")
    old_max_substeps = int(previous["electron_substeps_max_per_md_step"])

    print(
        f"RESUME step={start_step} time_ps={start_step*DT_PS:.4f} target_step={start_step+steps_to_run} "
        f"atoms={nlocal} fixed/bath/newton={int(group_masks['fixed'].sum())}/"
        f"{int(group_masks['bath'].sum())}/{int(group_masks['newton'].sum())} "
        f"Ag_Newton={len(ag_ids)} source_power={source_power.sum():.9f} eV/ps "
        f"Te_mean={np.average(grid.temperature_K[fraction_ag>0],weights=fraction_ag[fraction_ag>0]):.3f}K",
        flush=True,
    )

    history = list(previous.get("checkpoints", []))
    final_row = None
    for offset in range(1, steps_to_run + 1):
        lmp.command("run 1 pre no post no")
        velocity_ag = np.array(v[ag_indices], copy=True)
        lattice_temp = lattice_temperature_by_cell(
            velocity_ag,
            ag_masses,
            fields.atom_map,
            smoothing_radius_cells=1,
        )
        accounting = grid.advance(DT_PS, lattice_temp, heat_on=True, coupling_on=True)
        velocity_new, applied = rescale_lattice_kinetic(
            velocity_ag,
            ag_masses,
            fields.atom_map,
            accounting.lattice_energy_per_cell_eV,
        )
        v[ag_indices] = velocity_new
        segment_source += accounting.source_to_electron_eV
        segment_transfer += float(applied.sum())
        segment_grid_residual += accounting.conservation_residual_eV
        segment_substeps += accounting.substeps
        max_substeps = max(max_substeps, accounting.substeps)
        abs_step = start_step + offset

        if offset % 20 == 0 or offset == 1 or offset == steps_to_run:
            atom_now = float(lmp.get_thermo("pe")) + kinetic_eV(np.array(v[:nlocal], copy=True), masses_all)
            electron_now = float(grid.electron_energy_eV.sum())
            ecouple_now = float(lmp.get_thermo("ecouple"))
            segment_residual = (
                atom_now + electron_now - segment_atom_start - segment_electron_start
                - segment_source + (ecouple_now - segment_ecouple_start)
            )
            active = fraction_ag > 0.0
            print(
                f"step={abs_step} time_ps={abs_step*DT_PS:.5f} "
                f"Tbath={temperature_K(np.array(v[:nlocal][group_masks['bath']],copy=True),masses_all[group_masks['bath']]):.2f}K "
                f"Tnewton={temperature_K(np.array(v[:nlocal][group_masks['newton']],copy=True),masses_all[group_masks['newton']]):.2f}K "
                f"Tnewton_Ag={temperature_K(np.array(v[ag_indices],copy=True),ag_masses):.2f}K "
                f"Te_mean={np.average(grid.temperature_K[active],weights=fraction_ag[active]):.2f}K "
                f"source_total={old_total_source+segment_source:.6f}eV "
                f"e2l_total={old_total_transfer+segment_transfer:.6f}eV "
                f"segment_residual={segment_residual:.6e}eV "
                f"ttm_residual={old_total_grid_residual+segment_grid_residual:.3e}eV "
                f"substeps={accounting.substeps}",
                flush=True,
            )

        if abs_step % CHECKPOINT_INTERVAL == 0 or offset == steps_to_run:
            lmp.command("run 0 pre no post no")
            time_ps = abs_step * DT_PS
            tag = f"{time_ps:.4f}".rstrip("0").rstrip(".").replace(".", "p")
            out_ckpt = HERE / "checkpoints"
            out_ckpt.mkdir(parents=True, exist_ok=True)
            restart_out = out_ckpt / f"ttm_{tag}ps.restart"
            dump_out = out_ckpt / f"atoms_{tag}ps.dump"
            state_out = out_ckpt / f"electrons_{tag}ps.npz"
            lmp.command(f"write_restart {restart_out}")
            lmp.command(f"write_dump all custom {dump_out} id type x y z vx vy vz modify sort id")

            atom_now = float(lmp.get_thermo("pe")) + kinetic_eV(np.array(v[:nlocal], copy=True), masses_all)
            electron_now = float(grid.electron_energy_eV.sum())
            ecouple_now = float(lmp.get_thermo("ecouple"))
            segment_residual = (
                atom_now + electron_now - segment_atom_start - segment_electron_start
                - segment_source + (ecouple_now - segment_ecouple_start)
            )
            total_source = old_total_source + segment_source
            total_transfer = old_total_transfer + segment_transfer
            total_grid_residual = old_total_grid_residual + segment_grid_residual
            total_ecouple = old_total_ecouple + (ecouple_now - segment_ecouple_start)
            total_residual = old_total_residual + segment_residual
            active = fraction_ag > 0.0
            hot = source_power >= 0.01 * source_power.max()
            xyz_now = np.array(x[:nlocal], copy=True)
            type_now = atom_type
            tl = lattice_temperature_by_cell(
                np.array(v[ag_indices], copy=True),
                ag_masses,
                fields.atom_map,
                smoothing_radius_cells=1,
            )
            hot_tl = float(np.average(tl[hot], weights=source_power[hot]))
            for name in ("bath", "newton", "newton_ag", "newton_skeleton"):
                group_masks[name] = group_masks[name]
            mobile_momentum = np.sum(masses_all[:, None] * np.array(v[:nlocal], copy=True), axis=0)
            row = {
                "step": abs_step,
                "time_ps": time_ps,
                "electron_temperature_mean_K": float(np.average(grid.temperature_K[active], weights=fraction_ag[active])),
                "electron_temperature_max_K": float(np.max(grid.temperature_K[active])),
                "source_weighted_Ag_lattice_temperature_K": hot_tl,
                "source_energy_eV": total_source,
                "electron_to_lattice_eV": total_transfer,
                "electron_grid_energy_eV": electron_now,
                "atom_energy_eV": atom_now,
                "thermostat_ecouple_change_eV": total_ecouple,
                "combined_source_corrected_residual_eV": total_residual,
                "cumulative_ttm_grid_residual_eV": total_grid_residual,
                "ttm_checkpoint": restart_out.name,
                "dump": dump_out.name,
                "electron_state": state_out.name,
                "sputter_region_atoms_z150_200": int(np.count_nonzero((xyz_now[:, 2] >= 150.0) & (xyz_now[:, 2] < 200.0))),
                "evaporation_region_atoms_z200_250": int(np.count_nonzero((xyz_now[:, 2] >= 200.0) & (xyz_now[:, 2] < 250.0))),
                "max_z_A": float(np.max(xyz_now[:, 2])),
                "bath_temperature_K": temperature_K(np.array(v[:nlocal][group_masks["bath"]], copy=True), masses_all[group_masks["bath"]]),
                "newton_temperature_K": temperature_K(np.array(v[:nlocal][group_masks["newton"]], copy=True), masses_all[group_masks["newton"]]),
                "newton_ag_temperature_K": temperature_K(np.array(v[ag_indices], copy=True), ag_masses),
                "newton_skeleton_temperature_K": temperature_K(np.array(v[:nlocal][group_masks["newton_skeleton"]], copy=True), masses_all[group_masks["newton_skeleton"]]),
                "mobile_momentum_delta_g_mol_A_ps": mobile_momentum.tolist(),
                "atom_count": int(lmp.get_natoms()),
                "max_force_eV_A": float(lmp.get_thermo("fmax")),
                "continuation_segment_residual_eV": segment_residual,
            }
            history.append(row)
            np.savez_compressed(
                state_out,
                electron_temperature_K=grid.temperature_K,
                fraction_ag=fields.fraction_ag,
                source_power_eV_ps=fields.source_power_eV_ps,
                spacing_A=grid.spacing_A,
                atom_cell_ids=fields.atom_map.cell_ids,
                atom_cell_weights=fields.atom_map.weights,
                atom_ids_newton_ag=ag_ids,
                time_ps=np.array(time_ps),
                cumulative_source_eV=np.array(total_source),
                cumulative_electron_to_lattice_eV=np.array(total_transfer),
                cumulative_ttm_residual_eV=np.array(total_grid_residual),
            )
            final_row = row
            print(json.dumps(row), flush=True)

    if final_row is None:
        raise RuntimeError("No checkpoint was written")
    total_step = int(final_row["step"])
    result = dict(previous)
    result.update({
        "status": "EXPLORATORY_ROUND5_MTP_AG_ONLY_TTM_REFERENCE_COMPOSITE_CONTINUED_RUN",
        "duration_ps": total_step * DT_PS,
        "steps": total_step,
        "continuation_from": str(restart_path),
        "continuation_start_ps": start_step * DT_PS,
        "continuation_steps": total_step - start_step,
        "source": dict(previous["source"], total_source_energy_eV=float(final_row["source_energy_eV"])),
        "electron_to_lattice_energy_eV": float(final_row["electron_to_lattice_eV"]),
        "electron_temperature_final_mean_K": float(final_row["electron_temperature_mean_K"]),
        "electron_temperature_final_max_K": float(final_row["electron_temperature_max_K"]),
        "final_atom_energy_eV": float(final_row["atom_energy_eV"]),
        "final_electron_energy_eV": float(final_row["electron_grid_energy_eV"]),
        "thermostat_ecouple_change_eV": float(final_row["thermostat_ecouple_change_eV"]),
        "source_corrected_combined_energy_residual_eV": float(final_row["combined_source_corrected_residual_eV"]),
        "cumulative_ttm_grid_residual_eV": float(final_row["cumulative_ttm_grid_residual_eV"]),
        "electron_substeps_total": None if old_total_substeps is None else int(old_total_substeps) + segment_substeps,
        "electron_substeps_total_complete": old_total_substeps is not None and previous.get("electron_substeps_total_complete", True),
        "electron_substeps_since_resume_start": segment_substeps,
        "electron_substeps_max_per_md_step": max(old_max_substeps, max_substeps),
        "checkpoints": history,
    })
    result["limitations"] = list(previous["limitations"]) + [
        "This continuation preserves the static layer memberships from the original configuration, including atoms that crossed z-layer bounds before the restart.",
        "Continuation energy closure is accumulated from the previous checkpoint; the MTP energy remains uncalibrated.",
    ]
    (HERE / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    with (HERE / "checkpoints.csv").open("w", newline="") as stream:
        fieldnames = list(dict.fromkeys(key for item in history for key in item))
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(history)
    print(json.dumps(result, indent=2), flush=True)
    lmp.close()


if __name__ == "__main__":
    main()
