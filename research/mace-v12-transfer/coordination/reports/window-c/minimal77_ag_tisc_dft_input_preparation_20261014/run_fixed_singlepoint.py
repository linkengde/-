#!/usr/bin/env python3
"""Fixed-geometry GPAW input runner. Cannot run without a separate bound release JSON."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_fixed_singlepoint.py CONFIG.json EXECUTION_RELEASE.json")
    config_path = Path(sys.argv[1]).resolve()
    release_path = Path(sys.argv[2]).resolve()
    config = json.loads(config_path.read_text())
    structure_path = (ROOT / config["structure_path"]).resolve()
    config_sha = sha(config_path)
    structure_sha = sha(structure_path)
    runner_sha = sha(Path(__file__).resolve())
    release = json.loads(release_path.read_text()) if release_path.is_file() else {}
    required = {
        "execution_authorized": True,
        "task_id": config["task_id"],
        "config_sha256": config_sha,
        "structure_sha256": structure_sha,
        "runner_sha256": runner_sha,
    }
    mismatches = [key for key, expected in required.items() if release.get(key) != expected]
    if config.get("execution_authorized") is not False or config.get("launch_enabled") is not False:
        raise SystemExit("Input configs are immutable and must remain execution_authorized=false / launch_enabled=false")
    if mismatches:
        raise SystemExit("DFT blocked; missing or mismatched release fields: " + ", ".join(mismatches))
    if not structure_path.is_file() or structure_sha != config["structure_sha256"]:
        raise SystemExit("Structure missing or SHA mismatch")
    if os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("Set OMP_NUM_THREADS=1 before execution")

    atoms = read(structure_path)
    method = config["method"]
    symbols = atoms.get_chemical_symbols()
    symbols_sha = hashlib.sha256(" ".join(symbols).encode()).hexdigest()
    composition = {symbol: symbols.count(symbol) for symbol in sorted(set(symbols))}
    atom_ids_ok = ("atom_id" in atoms.arrays and
                   np.array_equal(atoms.arrays["atom_id"], np.arange(1, int(config["natoms"]) + 1)))
    if (len(atoms) != int(config["natoms"])
            or composition != config["composition"]
            or symbols_sha != config["ordered_symbols_sha256"]
            or not atom_ids_ok
            or atoms.pbc.tolist() != config["pbc"]
            or not np.allclose(atoms.cell.array, config["cell_A"], atol=1e-10)
            or not np.isfinite(atoms.positions).all()):
        raise SystemExit("Structure identity, atom order, cell, PBC, or coordinate check failed")

    import ase
    import gpaw
    import gpaw_data
    if ase.__version__ != "3.29.0" or gpaw.__version__ != "26.7.0":
        raise SystemExit(f"Runtime mismatch ASE={ase.__version__}, GPAW={gpaw.__version__}")
    if importlib.metadata.version("gpaw-data") != "1.2.1":
        raise SystemExit("gpaw-data version mismatch")
    setup_root = Path(os.environ.get("GPAW_SETUP_PATH", str(gpaw_data.datapath()))).resolve()
    os.environ["GPAW_SETUP_PATH"] = str(setup_root)
    for name, expected in method["paw_setup_sha256"].items():
        p = setup_root / name
        if not p.is_file() and (setup_root / "setups" / name).is_file():
            p = setup_root / "setups" / name
        if not p.is_file() or sha(p) != expected:
            raise SystemExit(f"PAW setup missing or SHA mismatch: {name}")

    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer
    from gpaw.mpi import world
    from gpaw.occupations import FermiDirac

    if world.size != 4:
        raise SystemExit(f"Expected four MPI ranks, got {world.size}")
    out = Path(release["output_dir"]).resolve()
    if out.exists():
        raise SystemExit(f"Refusing to overwrite existing result directory: {out}")
    if world.rank == 0:
        out.mkdir(parents=True, exist_ok=False)
    world.barrier()
    calc = GPAW(
        mode=PW(float(method["cutoff_eV"])), xc=method["xc"],
        kpts=tuple(method["kpoints"]),
        occupations=FermiDirac(float(method["occupations"]["width_eV"])),
        mixer=Mixer(float(method["mixer"]["beta"]), int(method["mixer"]["nmaxold"]), int(method["mixer"]["weight"])),
        convergence={
            "energy": float(method["convergence"]["energy_eV_per_electron"]),
            "density": float(method["convergence"]["density"]),
            "eigenstates": float(method["convergence"]["eigenstates_eV"]),
        },
        maxiter=int(method["maxiter"]), txt=str(out / "gpaw.log"),
    )
    atoms.calc = calc
    converged = {"value": False}
    calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
    start = time.time()
    energy = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if not (np.isfinite(energy) and np.isfinite(free_energy) and np.isfinite(forces).all() and np.isfinite(stress).all()):
        raise SystemExit("Nonfinite energy/force/stress result")
    atoms.info["native_energy_eV"] = energy
    atoms.info["force_consistent_free_energy_eV"] = free_energy
    atoms.info["scf_converged"] = bool(converged["value"])
    atoms.arrays["dft_forces_eV_per_A"] = forces
    write(out / "result.extxyz", atoms, format="extxyz")
    summary = {
        "task_id": config["task_id"], "config_sha256": config_sha,
        "structure_sha256": structure_sha, "runner_sha256": runner_sha,
        "native_energy_eV": energy, "force_consistent_free_energy_eV": free_energy,
        "forces_eV_per_A": forces.tolist(), "stress_eV_per_A3_diagnostic_only": stress.tolist(),
        "scf_converged": bool(converged["value"]), "scf_iterations": int(calc.get_number_of_iterations()),
        "elapsed_s": time.time() - start,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not converged["value"]:
        raise SystemExit("SCF did not converge; retain diagnostics and mark FAILED_SCF")


if __name__ == "__main__":
    main()
