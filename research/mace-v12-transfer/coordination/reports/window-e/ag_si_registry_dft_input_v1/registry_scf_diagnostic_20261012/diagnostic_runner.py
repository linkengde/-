#!/usr/bin/env python3
"""Isolated R1 SCF diagnostic runner; never modifies the frozen screening run."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read
from gpaw import FermiDirac, GPAW, PW
from gpaw.mixer import Mixer
from gpaw.mpi import world

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "R1_diagnostic_case"
STRUCTURE = CASE / "structure.xyz"
EXPECTED_STRUCTURE_SHA256 = "b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="read-only environment/geometry preflight")
    parser.add_argument("--execute", action="store_true", help="run this diagnostic copy only")
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--mixer-beta", type=float, default=0.10)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.execute == args.check:
        parser.error("choose exactly one of --check or --execute")
    if args.maxiter <= 160 or args.maxiter > 240:
        parser.error("diagnostic maxiter must be 161..240")
    if not (0.05 <= args.mixer_beta <= 0.12):
        parser.error("diagnostic mixer beta must be 0.05..0.12")

    import ase
    import gpaw
    from importlib.metadata import version
    from gpaw.setup_data import search_for_file

    if gpaw.__version__ != "26.7.0" or ase.__version__ != "3.29.0" or version("gpaw-data") != "1.2.1":
        raise RuntimeError("diagnostic runtime version mismatch")
    actual_structure_sha = sha256(STRUCTURE)
    if actual_structure_sha != EXPECTED_STRUCTURE_SHA256:
        raise RuntimeError("copied R1 structure SHA-256 mismatch")
    paw_hashes = {}
    expected_paw = {
        "Ag": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
        "C": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
        "Si": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
        "Ti": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
    }
    for symbol, expected in expected_paw.items():
        filename, _ = search_for_file(f"{symbol}.PBE")
        actual = sha256(Path(filename))
        if actual != expected:
            raise RuntimeError(f"PAW hash mismatch for {symbol}")
        paw_hashes[symbol] = actual

    if args.check:
        atoms = read(STRUCTURE, format="extxyz")
        print(json.dumps({
            "status": "DIAGNOSTIC_PREFLIGHT_PASS",
            "dft_started": False,
            "registry_id": "R1_DIAGNOSTIC_ONLY",
            "structure_sha256": actual_structure_sha,
            "n_atoms": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "pbc": atoms.pbc.tolist(),
            "kpoint_mesh": [6, 6, 1],
            "gpaw": gpaw.__version__, "ase": ase.__version__, "gpaw_data": version("gpaw-data"),
            "mpi_ranks": world.size,
            "maxiter": args.maxiter, "mixer_beta": args.mixer_beta,
            "paw_sha256": paw_hashes,
        }, indent=2))
        return 0

    if world.size != 4:
        raise RuntimeError("diagnostic run requires the validated four-rank allocation")
    if args.output_dir is None:
        parser.error("--output-dir is required with --execute")
    out = args.output_dir.resolve()
    if world.rank == 0:
        out.mkdir(parents=True, exist_ok=False)
    world.barrier()

    atoms = read(STRUCTURE, format="extxyz")
    before_positions = atoms.positions.copy()
    before_cell = atoms.cell.array.copy()
    before_pbc = atoms.pbc.copy()
    log_path = out / "gpaw.log"
    result_path = out / "result.json"
    calc = GPAW(
        mode=PW(500), xc="PBE", kpts=(6, 6, 1),
        occupations=FermiDirac(0.1),
        convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
        maxiter=args.maxiter, mixer=Mixer(args.mixer_beta, 8, 100),
        symmetry={"point_group": True, "time_reversal": True},
        poissonsolver={}, parallel={"sl_auto": True}, txt=str(log_path),
    )
    atoms.calc = calc
    try:
        energy = float(atoms.get_potential_energy(force_consistent=False))
        free_energy = float(atoms.get_potential_energy(force_consistent=True))
        forces = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
        stress = np.asarray(atoms.get_stress(voigt=False, apply_constraint=False), dtype=float)
        if not np.array_equal(atoms.positions, before_positions):
            raise RuntimeError("diagnostic mutated atomic coordinates")
        if not np.array_equal(atoms.cell.array, before_cell) or not np.array_equal(atoms.pbc, before_pbc):
            raise RuntimeError("diagnostic mutated cell/PBC")
        if world.rank == 0:
            log_text = log_path.read_text(errors="replace") if log_path.exists() else ""
            match = re.search(r"Converged in (\d+) steps", log_text)
            result = {
                "status": "DIAGNOSTIC_COMPLETE",
                "registry_id": "R1_DIAGNOSTIC_ONLY",
                "source_structure_sha256": actual_structure_sha,
                "kpoint_mesh": [6, 6, 1], "cutoff_eV": 500, "xc": "PBE",
                "smearing_eV": 0.1, "maxiter": args.maxiter,
                "mixer": {"class": "Mixer", "beta": args.mixer_beta, "nmaxold": 8, "weight": 100},
                "scf_converged": bool(match), "scf_steps": int(match.group(1)) if match else None,
                "total_energy_eV": energy, "force_consistent_free_energy_eV": free_energy,
                "forces_eV_A": forces.tolist(), "stress_eV_A3": stress.tolist(),
                "geometry_modified": False, "paw_sha256": paw_hashes,
            }
            result_path.write_text(json.dumps(result, indent=2) + "\n")
            print(json.dumps({"status": result["status"], "scf_converged": result["scf_converged"], "scf_steps": result["scf_steps"]}, indent=2))
    except BaseException as exc:
        if world.rank == 0:
            (out / "failure.json").write_text(json.dumps({
                "status": "DIAGNOSTIC_FAILED", "exception_type": type(exc).__name__, "message": str(exc),
                "maxiter": args.maxiter, "mixer_beta": args.mixer_beta,
            }, indent=2) + "\n")
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
