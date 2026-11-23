#!/workspace/.venvs/gpaw-diagnosis/bin/python
"""GPAW 500 eV/3x3x3 单点计算的初始化预检，不进入 SCF。"""
from pathlib import Path
import json
import gpaw
from ase.io import read
from gpaw import GPAW, PW, FermiDirac
from gpaw.mixer import Mixer

ROOT = Path(__file__).resolve().parent
# 该脚本只允许通过 GPAW dry-run 初始化，不得意外进入昂贵的 SCF。
gpaw.dry_run = 1
atoms = read(ROOT / 'input.extxyz')
assert len(atoms) == 64 and atoms.pbc.tolist() == [True, True, True]
assert set(atoms.get_chemical_symbols()) == {'Ag'}
calc = GPAW(mode=PW(500), xc='PBE', kpts=(2, 2, 2), symmetry='off',
            occupations=FermiDirac(0.10),
            mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
            convergence={'energy': 1e-5, 'density': 1e-5, 'eigenstates': 1e-5},
            maxiter=160,
            parallel={'kpt': 4, 'domain': 1, 'band': 1, 'sl_auto': False},
            txt=str(ROOT / 'preflight_gpaw.log'))
atoms.calc = calc
calc.initialize(atoms)
raise SystemExit('PRECHECK_ERROR: dry-run initialization returned unexpectedly')
