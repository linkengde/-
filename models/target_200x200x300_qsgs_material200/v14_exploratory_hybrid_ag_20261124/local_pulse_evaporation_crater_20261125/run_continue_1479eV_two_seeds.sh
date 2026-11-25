#!/usr/bin/env bash
set -euo pipefail
LMP_BIN="${LMP_BIN:-/workspace/venvs/podfit/bin/lmp}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
"${LMP_BIN}" -in in.continue_1479eV_seed20261133_2ps.lmp -log log.continue_1479eV_seed20261133_2ps.lammps -screen screen.continue_1479eV_seed20261133_2ps.out
"${LMP_BIN}" -in in.continue_1479eV_seed20261134_2ps.lmp -log log.continue_1479eV_seed20261134_2ps.lammps -screen screen.continue_1479eV_seed20261134_2ps.out
