#!/usr/bin/env bash
set -euo pipefail

LMP_BIN="${LMP_BIN:-/workspace/venvs/podfit/bin/lmp}"
OVITO_PYTHON="${OVITO_PYTHON:-/workspace/software/ovito-lavapipe/ovito-python.sh}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

python prepare_gaussian_pulse.py
"${LMP_BIN}" -in in.local_gaussian_pulse_1ps_nve.lmp -log log.local_gaussian_pulse_1ps_nve.lammps -screen screen.local_gaussian_pulse_1ps_nve.out
python prepare_gaussian_pulse_9000K.py
"${LMP_BIN}" -in in.local_gaussian_pulse_9000K_nve.lmp -log log.local_gaussian_pulse_9000K_nve.lammps -screen screen.local_gaussian_pulse_9000K_nve.out
python prepare_gaussian_pulse_22500K.py
"${LMP_BIN}" -in in.local_gaussian_pulse_22500K_nve.lmp -log log.local_gaussian_pulse_22500K_nve.lammps -screen screen.local_gaussian_pulse_22500K_nve.out
python prepare_gaussian_pulse_surface_local.py
"${LMP_BIN}" -in in.local_gaussian_pulse_surface_local_2ps.lmp -log log.local_gaussian_pulse_surface_local_2ps.lammps -screen screen.local_gaussian_pulse_surface_local_2ps.out
python prepare_gaussian_pulse_surface_local_33750.py
"${LMP_BIN}" -in in.local_gaussian_pulse_surface_local_33750_2ps.lmp -log log.local_gaussian_pulse_surface_local_33750_2ps.lammps -screen screen.local_gaussian_pulse_surface_local_33750_2ps.out
python analyze_pulse.py
python analyze_surface_local_33750.py
python analyze_ejection_timeseries.py
python analyze_surface_local.py
python analyze_escape_clusters.py
"${OVITO_PYTHON}" render_center_slice_ovito.py
"${OVITO_PYTHON}" render_center_slice_22500K_ovito.py
"${OVITO_PYTHON}" render_center_slice_surface_local_2ps_ovito.py
"${OVITO_PYTHON}" render_center_slice_surface_local_33750_2ps_ovito.py
