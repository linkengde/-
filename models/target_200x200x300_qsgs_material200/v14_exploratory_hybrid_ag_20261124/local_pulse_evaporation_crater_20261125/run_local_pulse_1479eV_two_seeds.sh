#!/usr/bin/env bash
set -euo pipefail
LMP_BIN="${LMP_BIN:-/workspace/venvs/podfit/bin/lmp}"
OVITO_PYTHON="${OVITO_PYTHON:-/workspace/software/ovito-lavapipe/ovito-python.sh}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
python prepare_local_pulse_1479eV_two_seeds.py
mkdir -p trajectory_surface_local_1479eV_seed20261133_2ps trajectory_surface_local_1479eV_seed20261134_2ps
"${LMP_BIN}" -in in.local_gaussian_pulse_1479eV_seed20261133_2ps.lmp -log log.local_gaussian_pulse_1479eV_seed20261133_2ps.lammps
"${LMP_BIN}" -in in.local_gaussian_pulse_1479eV_seed20261134_2ps.lmp -log log.local_gaussian_pulse_1479eV_seed20261134_2ps.lammps
python analyze_1479eV_two_seeds.py
"${OVITO_PYTHON}" render_1479eV_two_seeds_ovito.py
