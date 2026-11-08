#!/usr/bin/env python3
"""Exploratory Round-5 MTP + Ag-only TTM reference coupling on finer QSGS."""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
REFERENCE_MODEL = MODEL.parent / 'ag_ti3sic2_qsgs_1to1_finer'
sys.path.insert(0, str(REFERENCE_MODEL))

from lammps import lammps
from ag_ttm_reference import (
    AgParameters,
    AgTTMGrid,
    KB_EV_K,
    M_VV2E_METAL,
    build_ag_grid_fields,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)


STEPS = int(os.environ.get('TTM_HEATING_STEPS', '1000'))
CHECKPOINT_STEPS = int(os.environ.get('TTM_CHECKPOINT_STEPS', '500'))
SOURCE_SCALE = float(os.environ.get('TTM_SOURCE_SCALE', '1.0'))
DT_PS = 0.0001
GRID_SPACING_A = 10.0
R_A = 30.0
SKIN_A = 20.0
I0_EVP_S_A2 = 5.0 * SOURCE_SCALE
SEED = 190427
ATOMIC_MASSES = {1: 47.867, 2: 28.085, 3: 12.011, 4: 107.8682}
AG_ATOMIC_VOLUME_A3 = 4.09**3 / 4.0

if STEPS <= 0 or CHECKPOINT_STEPS <= 0 or not np.isfinite(SOURCE_SCALE) or SOURCE_SCALE < 0.0:
    raise ValueError('Invalid step count, checkpoint interval, or source scale')


def kinetic_eV(velocity: np.ndarray, mass: np.ndarray) -> float:
    return float(0.5 * M_VV2E_METAL * np.sum(mass[:, None] * velocity**2))


def temperature_K(velocity: np.ndarray, mass: np.ndarray) -> float:
    if len(mass) == 0:
        return float('nan')
    return 2.0 * kinetic_eV(velocity, mass) / (3.0 * len(mass) * KB_EV_K)


def checkpoint(
    lmp,
    step: int,
    grid: AgTTMGrid,
    fields,
    ag_ids: np.ndarray,
    cumulative_source: float,
    cumulative_transfer: float,
    cumulative_grid_residual: float,
    baseline_total_eV: float,
    baseline_electron_eV: float,
    baseline_ecouple_eV: float,
    atom_ids: np.ndarray,
    atom_type: np.ndarray,
    x: np.ndarray,
    v: np.ndarray,
    masses_all: np.ndarray,
    group_masks: dict[str, np.ndarray],
    initial_momentum: np.ndarray,
    history: list[dict],
) -> dict:
    lmp.command('run 0 pre no post no')
    time_ps = step * DT_PS
    tag = f'{time_ps:.3f}'.rstrip('0').rstrip('.').replace('.', 'p')
    restart_path = HERE / 'checkpoints' / f'ttm_{tag}ps.restart'
    dump_path = HERE / 'checkpoints' / f'atoms_{tag}ps.dump'
    state_path = HERE / 'checkpoints' / f'electrons_{tag}ps.npz'
    lmp.command(f'write_restart {restart_path}')
    lmp.command(f'write_dump all custom {dump_path} id type x y z vx vy vz modify sort id')

    v_all = np.array(v[:len(atom_ids)], copy=True)
    pe = float(lmp.get_thermo('pe'))
    atom_total = pe + kinetic_eV(v_all, masses_all)
    electron_total = float(grid.electron_energy_eV.sum())
    ecouple = float(lmp.get_thermo('ecouple'))
    residual = (
        atom_total + electron_total - baseline_total_eV - baseline_electron_eV
        - cumulative_source + (ecouple - baseline_ecouple_eV)
    )
    active = fields.fraction_ag > 0.0
    hot = fields.source_power_eV_ps >= 0.01 * fields.source_power_eV_ps.max()
    tl = lattice_temperature_by_cell(
        np.array(v[:len(atom_ids)][group_masks['newton_ag']], copy=True),
        np.full(int(group_masks['newton_ag'].sum()), ATOMIC_MASSES[4]),
        fields.atom_map,
        smoothing_radius_cells=1,
    )
    hot_tl = float(np.average(
        tl[hot], weights=fields.source_power_eV_ps[hot]
    ))
    active_te = grid.temperature_K[active]
    type_now = atom_type[:len(atom_ids)]
    xyz_now = x[:len(atom_ids)]
    counts = {
        'sputter_region_atoms_z150_200': int(np.count_nonzero((xyz_now[:, 2] >= 150.0) & (xyz_now[:, 2] < 200.0))),
        'evaporation_region_atoms_z200_250': int(np.count_nonzero((xyz_now[:, 2] >= 200.0) & (xyz_now[:, 2] < 250.0))),
        'max_z_A': float(np.max(xyz_now[:, 2])),
    }
    for name in ('bath', 'newton', 'newton_ag', 'newton_skeleton'):
        mask = group_masks[name]
        counts[f'{name}_temperature_K'] = temperature_K(v_all[mask], masses_all[mask])
    momentum = np.sum(masses_all[:, None] * v_all, axis=0)
    counts['mobile_momentum_delta_g_mol_A_ps'] = (momentum - initial_momentum).tolist()
    counts['atom_count'] = int(lmp.get_natoms())
    counts['max_force_eV_A'] = float(lmp.get_thermo('fmax'))

    np.savez_compressed(
        state_path,
        electron_temperature_K=grid.temperature_K,
        fraction_ag=fields.fraction_ag,
        source_power_eV_ps=fields.source_power_eV_ps,
        spacing_A=grid.spacing_A,
        atom_cell_ids=fields.atom_map.cell_ids,
        atom_cell_weights=fields.atom_map.weights,
        atom_ids_newton_ag=ag_ids,
        time_ps=np.array(time_ps),
        cumulative_source_eV=np.array(cumulative_source),
        cumulative_electron_to_lattice_eV=np.array(cumulative_transfer),
        cumulative_ttm_residual_eV=np.array(cumulative_grid_residual),
    )
    row = {
        'step': step,
        'time_ps': time_ps,
        'electron_temperature_mean_K': float(np.average(active_te, weights=fields.fraction_ag[active])),
        'electron_temperature_max_K': float(np.max(active_te)),
        'source_weighted_Ag_lattice_temperature_K': hot_tl,
        'source_energy_eV': cumulative_source,
        'electron_to_lattice_eV': cumulative_transfer,
        'electron_grid_energy_eV': electron_total,
        'atom_energy_eV': atom_total,
        'thermostat_ecouple_change_eV': ecouple - baseline_ecouple_eV,
        'combined_source_corrected_residual_eV': residual,
        'cumulative_ttm_grid_residual_eV': cumulative_grid_residual,
        'ttm_checkpoint': restart_path.name,
        'dump': dump_path.name,
        'electron_state': state_path.name,
        **counts,
    }
    history.append(row)
    return row


def main() -> None:
    os.chdir(MODEL)
    report = json.loads((MODEL / 'geometry_report.json').read_text())
    lx, ly, lz = report['material_dimensions_A']
    lengths = np.array([lx, ly, lz], dtype=float)
    grid_shape = np.ceil(lengths / GRID_SPACING_A).astype(int)
    spacing = lengths / grid_shape
    OUT_LOG = HERE / 'lammps.log'
    MODEL_DATA = MODEL / 'ag_ti3sic2_qsgs.data.gz'
    OUT_CKPT = HERE / 'checkpoints'
    MTP_DIR = (REFERENCE_MODEL / 'potential_validation' / 'mtp_thermal_acquisition_pilot_20261106'
               / 'mtp4_ag108_ti3sic2_tscbop_thermal_interface_weighted_round5')
    MTP_PATH = MTP_DIR / 'MTP1m_AgTiSiC_thermal_interface_weighted_round5_DEVELOPMENT_ONLY.mtp'
    MTP_INI = HERE / 'mlip.ini'
    MTP_INI.write_text(f'mtp-filename {MTP_PATH}\nselect FALSE\n')
    OUT_CKPT.mkdir(parents=True, exist_ok=True)

    lmp = lammps(name='mtp', cmdargs=['-screen', 'none', '-log', str(OUT_LOG)])
    commands = [
        'units metal',
        'atom_style atomic',
        'boundary p p f',
        'atom_modify sort 0 0',
        f'read_data {MODEL_DATA}',
        f'pair_style mlip {MTP_INI}',
        'pair_coeff * *',
        'neighbor 2.0 bin',
        'neigh_modify delay 0 every 1 check yes',
        'region r_fixed block INF INF INF INF 0.0 20.0 units box',
        'region r_bath block INF INF INF INF 20.0 40.0 units box',
        'region r_newton block INF INF INF INF 40.0 150.0 units box',
        'region r_sputter block INF INF INF INF 150.0 200.0 units box',
        'region r_evaporation block INF INF INF INF 200.0 250.0 units box',
        'group fixed region r_fixed',
        'group bath region r_bath',
        'group newton region r_newton',
        'group Ag type 4',
        'group newton_Ag intersect newton Ag',
        'group newton_skeleton type 1 2 3',
        'group mobile union bath newton',
        'group sputter dynamic all region r_sputter every 10',
        'group evaporation dynamic all region r_evaporation every 10',
        f'velocity mobile create 300.0 {SEED} mom yes rot no dist gaussian',
        'velocity fixed set 0.0 0.0 0.0',
        f'timestep {DT_PS}',
        'fix freeze fixed setforce 0.0 0.0 0.0',
        'fix bath_control bath nvt temp 300.0 300.0 0.1',
        'fix integrate_newton newton nve',
        'compute T_bath bath temp',
        'compute T_newton newton temp',
        'thermo 20',
        'thermo_style custom step time atoms temp c_T_bath c_T_newton pe ke etotal ecouple econserve fmax',
        'thermo_modify flush yes',
        'run 0',
    ]
    for command in commands:
        lmp.command(command)

    nlocal = int(lmp.extract_setting('nlocal'))
    if nlocal != int(report['atoms']):
        raise RuntimeError(f'This reference pilot expects one rank with all atoms local; got {nlocal}')
    x = lmp.numpy.extract_atom('x')
    v = lmp.numpy.extract_atom('v')
    ids_ptr = lmp.numpy.extract_atom('id')
    type_ptr = lmp.numpy.extract_atom('type')
    ids = np.asarray(ids_ptr[:nlocal], dtype=np.int64).copy()
    atom_type = np.asarray(type_ptr[:nlocal], dtype=np.int32).copy()
    xyz = np.array(x[:nlocal], copy=True)
    velocity = np.array(v[:nlocal], copy=True)
    masses_all = np.array([ATOMIC_MASSES[int(t)] for t in atom_type])
    z0 = xyz[:, 2]
    masks = {
        'fixed': z0 < 20.0,
        'bath': (z0 >= 20.0) & (z0 < 40.0),
        'newton': (z0 >= 40.0) & (z0 < 150.0),
        'newton_ag': (atom_type == 4) & (z0 >= 40.0) & (z0 < 150.0),
        'newton_skeleton': np.isin(atom_type, (1, 2, 3)) & (z0 >= 40.0) & (z0 < 150.0),
    }
    if int(masks['newton_ag'].sum()) == 0:
        raise RuntimeError('No Newton-layer Ag atoms mapped to TTM')
    ag_ids = ids[masks['newton_ag']]
    ag_xyz = xyz[masks['newton_ag']]

    dx = ag_xyz[:, 0] - lx / 2.0
    dy = ag_xyz[:, 1] - ly / 2.0
    dx = (dx + lx / 2.0) % lx - lx / 2.0
    dy = (dy + ly / 2.0) % ly - ly / 2.0
    depth = lz - ag_xyz[:, 2]
    q_atom = (I0_EVP_S_A2 / SKIN_A) * np.exp(
        -(dx * dx + dy * dy) / (R_A * R_A) - depth / SKIN_A
    )
    fields = build_ag_grid_fields(
        ag_xyz,
        q_atom,
        tuple(int(v) for v in grid_shape),
        tuple(float(v) for v in spacing),
        AG_ATOMIC_VOLUME_A3,
    )
    parameters = AgParameters.from_csv(
        REFERENCE_MODEL / 'ttm_data' / 'Ag_TTM_hybrid_figure_fit_300_20000K.csv',
        allow_low_temperature_extension=True,
    )
    grid = AgTTMGrid(
        parameters,
        fields.fraction_ag,
        tuple(float(v) for v in spacing),
        fields.source_power_eV_ps,
        initial_temperature_K=300.0,
    )
    e0 = float(grid.electron_energy_eV.sum())
    atom_e0 = float(lmp.get_thermo('etotal'))
    ecouple0 = float(lmp.get_thermo('ecouple'))
    initial_momentum = np.sum(masses_all[:, None] * velocity, axis=0)
    cumulative_source = 0.0
    cumulative_transfer = 0.0
    cumulative_grid_residual = 0.0
    total_substeps = 0
    max_substeps = 0
    history: list[dict] = []
    print(
        f'INITIAL atoms={nlocal} Ag_Newton={int(masks["newton_ag"].sum())} '
        f'grid={tuple(int(v) for v in grid_shape)} spacing_A={spacing.tolist()} '
        f'I0={I0_EVP_S_A2:g} eV/(ps A^2) R={R_A:g} A skin={SKIN_A:g} A '
        f'source_power={fields.source_power_eV_ps.sum():.9f} eV/ps '
        f'initial_atom_etotal={atom_e0:.9f} eV electron_energy={e0:.9f} eV',
        flush=True,
    )

    for step in range(1, STEPS + 1):
        lmp.command('run 1 pre no post no')
        velocity_md = np.array(v[:nlocal][masks['newton_ag']], copy=True)
        lattice_temp = lattice_temperature_by_cell(
            velocity_md,
            np.full(len(ag_ids), ATOMIC_MASSES[4]),
            fields.atom_map,
            smoothing_radius_cells=1,
        )
        accounting = grid.advance(DT_PS, lattice_temp, heat_on=True, coupling_on=True)
        velocity_new, applied = rescale_lattice_kinetic(
            velocity_md,
            np.full(len(ag_ids), ATOMIC_MASSES[4]),
            fields.atom_map,
            accounting.lattice_energy_per_cell_eV,
        )
        v[:nlocal][masks['newton_ag']] = velocity_new
        cumulative_source += accounting.source_to_electron_eV
        cumulative_transfer += float(applied.sum())
        cumulative_grid_residual += accounting.conservation_residual_eV
        total_substeps += accounting.substeps
        max_substeps = max(max_substeps, accounting.substeps)

        if step % 20 == 0 or step == 1 or step == STEPS:
            pe = float(lmp.get_thermo('pe'))
            atom_e = pe + kinetic_eV(np.array(v[:nlocal], copy=True), masses_all)
            electron_e = float(grid.electron_energy_eV.sum())
            ecouple = float(lmp.get_thermo('ecouple'))
            residual = atom_e + electron_e - atom_e0 - e0 - cumulative_source + (ecouple - ecouple0)
            active = fields.fraction_ag > 0.0
            print(
                f'step={step} time_ps={step*DT_PS:.5f} '
                f'Tbath={temperature_K(np.array(v[:nlocal][masks["bath"]],copy=True),masses_all[masks["bath"]]):.2f}K '
                f'Tnewton={temperature_K(np.array(v[:nlocal][masks["newton"]],copy=True),masses_all[masks["newton"]]):.2f}K '
                f'Tnewton_Ag={temperature_K(np.array(v[:nlocal][masks["newton_ag"]],copy=True),np.full(len(ag_ids),ATOMIC_MASSES[4])):.2f}K '
                f'Te_mean={np.average(grid.temperature_K[active],weights=fields.fraction_ag[active]):.2f}K '
                f'Te_max={grid.temperature_K[active].max():.2f}K '
                f'source={cumulative_source:.6f}eV e2l={cumulative_transfer:.6f}eV '
                f'ecouple_delta={ecouple-ecouple0:.6f}eV residual={residual:.6e}eV '
                f'ttm_residual={cumulative_grid_residual:.3e}eV substeps={accounting.substeps}',
                flush=True,
            )

        if step % CHECKPOINT_STEPS == 0 or step == STEPS:
            row = checkpoint(
                lmp, step, grid, fields, ag_ids, cumulative_source,
                cumulative_transfer, cumulative_grid_residual, atom_e0, e0,
                ecouple0, ids, atom_type, x, v, masses_all, masks,
                initial_momentum, history,
            )
            print(json.dumps(row), flush=True)

    lmp.command('run 0 pre no post no')
    final_atom_e = float(lmp.get_thermo('pe')) + kinetic_eV(np.array(v[:nlocal], copy=True), masses_all)
    final_electron_e = float(grid.electron_energy_eV.sum())
    final_ecouple = float(lmp.get_thermo('ecouple'))
    final_residual = final_atom_e + final_electron_e - atom_e0 - e0 - cumulative_source + (final_ecouple - ecouple0)
    result = {
        'status': 'EXPLORATORY_ROUND5_MTP_AG_ONLY_TTM_REFERENCE_COMPOSITE_RUN',
        'model': str(MODEL),
        'reference_operator': str(REFERENCE_MODEL / 'ag_ttm_reference.py'),
        'potential': str(MTP_PATH),
        'potential_classification': 'force-only development candidate; energies and virials are uncalibrated',
        'duration_ps': STEPS * DT_PS,
        'steps': STEPS,
        'dt_fs': DT_PS * 1000.0,
        'atoms': int(lmp.get_natoms()),
        'newton_layer_ag_atoms_initial': int(len(ag_ids)),
        'electron_grid_shape': [int(v) for v in grid_shape],
        'electron_grid_spacing_A': spacing.tolist(),
        'ag_grid_fill_diagnostics': {
            'raw_overfilled_cells': int(fields.raw_overfilled_cells),
            'raw_max_fill_fraction': float(fields.raw_max_fill_fraction),
            'redistributed_ag_volume_A3': float(fields.redistributed_ag_volume_A3),
        },
        'source': {
            'I0_eV_ps_A2': I0_EVP_S_A2,
            'R_A': R_A,
            'skin_A': SKIN_A,
            'surface_z_A': lz,
            'center_xy_A': [lx / 2.0, ly / 2.0],
            'integrated_source_power_eV_ps': float(fields.source_power_eV_ps.sum()),
            'total_source_energy_eV': cumulative_source,
            'source_scale': SOURCE_SCALE,
        },
        'electron_to_lattice_energy_eV': cumulative_transfer,
        'electron_temperature_final_mean_K': float(np.average(
            grid.temperature_K[fields.fraction_ag > 0.0],
            weights=fields.fraction_ag[fields.fraction_ag > 0.0],
        )),
        'electron_temperature_final_max_K': float(grid.temperature_K[fields.fraction_ag > 0.0].max()),
        'final_atom_energy_eV': final_atom_e,
        'final_electron_energy_eV': final_electron_e,
        'thermostat_ecouple_change_eV': final_ecouple - ecouple0,
        'source_corrected_combined_energy_residual_eV': final_residual,
        'cumulative_ttm_grid_residual_eV': cumulative_grid_residual,
        'electron_substeps_total': total_substeps,
        'electron_substeps_max_per_md_step': max_substeps,
        'checkpoints': history,
        'limitations': [
            'TTM electronic subsystem, source, and electron-lattice exchange are Ag-only; Ti3SiC2 receives heat through atomistic interactions.',
            'Reference coupling is deterministic mean-field velocity scaling, not a compiled LAMMPS fix and not a production-validated TTM-MD method.',
            'The initial geometry is unrelaxed and the TTM atom-to-cell map is held fixed during this short run.',
            'The candidate Ag Ce/G/ke table and its sub-300 K continuation are exploratory fits.',
            'The round-5 MTP is force-only; its energy and virial are not calibrated, so energy ledgers are numerical diagnostics only.',
            'This run does not by itself establish melting, a stable melt pool, or physical interfacial heat transfer.',
        ],
    }
    (HERE / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    with (HERE / 'checkpoints.csv').open('w', newline='') as stream:
        if history:
            writer = csv.DictWriter(stream, fieldnames=list(history[0]))
            writer.writeheader()
            writer.writerows(history)
    print(json.dumps(result, indent=2), flush=True)
    lmp.close()


if __name__ == '__main__':
    main()
