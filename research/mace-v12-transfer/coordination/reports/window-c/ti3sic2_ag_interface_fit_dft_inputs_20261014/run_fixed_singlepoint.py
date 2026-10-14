#!/usr/bin/env python3
"""Authorization-gated GPAW single-point runner; this package contains no authorization."""
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


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_fixed_singlepoint.py CONFIG.json AUTHORIZATION.json")
    config_path = Path(sys.argv[1]).resolve()
    auth_path = Path(sys.argv[2]).resolve()
    config = json.loads(config_path.read_text())
    structure_path = (ROOT / config["structure_path"]).resolve()
    config_sha = digest(config_path)
    structure_sha = digest(structure_path)
    runner_sha = digest(Path(__file__).resolve())
    auth = json.loads(auth_path.read_text()) if auth_path.is_file() else {}
    required = {
        "execution_authorized": True,
        "task_id": config["task_id"],
        "config_sha256": config_sha,
        "structure_sha256": structure_sha,
        "runner_sha256": runner_sha,
    }
    mismatches = [key for key, expected in required.items() if auth.get(key) != expected]
    if config.get("execution_authorized") is not False:
        raise SystemExit("Package configs are input-only and must remain execution_authorized=false")
    if mismatches:
        raise SystemExit("Execution blocked; missing/mismatched authorization fields: " + ", ".join(mismatches))
    method = config["method"]
    if os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("Set OMP_NUM_THREADS=1 before launch")

    atoms = read(structure_path)
    if len(atoms) != int(config["natoms"]):
        raise SystemExit("Atom count mismatch")
    if digest(structure_path) != config["structure_sha256"]:
        raise SystemExit("Structure SHA mismatch")
    if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all():
        raise SystemExit("Nonfinite geometry")
    if any(key in atoms.info for key in ("energy", "REF_energy", "PW_PBE_energy_eV")):
        raise SystemExit("Unexpected energy label in input geometry")
    if any(key in atoms.arrays for key in ("forces", "REF_forces", "PW_PBE_forces")):
        raise SystemExit("Unexpected force label in input geometry")

    import ase
    import gpaw
    import gpaw_data
    if ase.__version__ != "3.29.0" or gpaw.__version__ != "26.7.0":
        raise SystemExit(f"Runtime version mismatch: ASE={ase.__version__}, GPAW={gpaw.__version__}; expected ASE 3.29.0/GPAW 26.7.0")
    try:
        paw_pkg_version = importlib.metadata.version("gpaw-data")
    except importlib.metadata.PackageNotFoundError:
        paw_pkg_version = "unknown"
    if paw_pkg_version != "1.2.1":
        raise SystemExit(f"gpaw-data version mismatch: {paw_pkg_version}; expected 1.2.1")
    setup_root = Path(os.environ.get("GPAW_SETUP_PATH", str(gpaw_data.datapath()))).resolve()
    os.environ["GPAW_SETUP_PATH"] = str(setup_root)
    for setup_name, expected_sha in method["paw_setup_sha256"].items():
        setup_path = setup_root / setup_name
        if not setup_path.is_file() and (setup_root / "setups" / setup_name).is_file():
            setup_path = setup_root / "setups" / setup_name
        if not setup_path.is_file() or digest(setup_path) != expected_sha:
            raise SystemExit(f"PAW setup missing or SHA mismatch: {setup_name}")

    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer
    from gpaw.mpi import world
    from gpaw.occupations import FermiDirac

    if world.size != 4:
        raise SystemExit(f"Expected exactly 4 MPI ranks, got {world.size}")
    out = Path(auth["output_dir"]).resolve()
    if out.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {out}")
    if world.rank == 0:
        out.mkdir(parents=True, exist_ok=False)
    world.barrier()
    if "Gamma-centered" in method["kpoint_convention"]:
        kpts = {"size": tuple(method["kpoints"]), "gamma": True}
    else:
        kpts = tuple(method["kpoints"])
    calc = GPAW(
        mode=PW(float(method["cutoff_eV"])),
        xc=method["xc"],
        kpts=kpts,
        occupations=FermiDirac(float(method["occupations"]["width_eV"])),
        mixer=Mixer(float(method["mixer"]["beta"]), int(method["mixer"]["nmaxold"]), int(method["mixer"]["weight"])),
        parallel={"sl_auto": True},
        convergence={
            "energy": float(method["convergence"]["energy_eV_per_electron"]),
            "density": float(method["convergence"]["density"]),
            "eigenstates": float(method["convergence"]["eigenstates_eV"]),
        },
        maxiter=int(method["maxiter"]),
        txt=str(out / "gpaw.log"),
    )
    atoms.calc = calc
    converged = {"value": False}
    calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
    start = time.time()
    native_energy = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(calc.results["free_energy"]) if "free_energy" in calc.results else None
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if not (np.isfinite(native_energy) and np.isfinite(forces).all() and np.isfinite(stress).all() and (free_energy is None or np.isfinite(free_energy))):
        raise SystemExit("Nonfinite DFT result; retain log and mark failed")
    atoms.info["native_energy_eV"] = native_energy
    atoms.info["force_consistent_free_energy_eV"] = free_energy
    atoms.info["scf_converged"] = bool(converged["value"])
    atoms.info["task_id"] = config["task_id"]
    atoms.arrays["dft_forces_eV_per_A"] = forces
    write(out / "result.extxyz", atoms, format="extxyz")
    summary = {
        "task_id": config["task_id"],
        "structure_sha256": structure_sha,
        "config_sha256": config_sha,
        "runner_sha256": runner_sha,
        "native_energy_eV": native_energy,
        "force_consistent_free_energy_eV": free_energy,
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3": stress.tolist(),
        "stress_usage": "bulk target" if config["parent_family"] == "TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY" else "slab diagnostic only; vacuum-diluted tensor not a bulk fitting target",
        "scf_converged": bool(converged["value"]),
        "scf_iterations": int(calc.get_number_of_iterations()),
        "elapsed_s": time.time() - start,
        "dataset_role": config["proposed_role"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not converged["value"]:
        raise SystemExit("SCF did not converge; result must not enter the labeled dataset")


if __name__ == "__main__":
    main()
