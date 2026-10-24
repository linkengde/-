#!/usr/bin/env python3
"""Run one fixed-geometry GPAW diagnostic; never writes to a training dataset."""

from __future__ import annotations

import hashlib
import gzip
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read, write

from gpaw.mpi import world


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "Ti3SiC2_Bai2010_DFT_candidate.cif"
OUT = ROOT / "exploratory_dft" / "run_01"
METHOD = {
    "code": "GPAW 26.7.0",
    "xc": "PBE",
    "basis": "plane-wave",
    "cutoff_eV": 500,
    "kpts": [12, 12, 2],
    "smearing_eV": 0.1,
    "mixer": "Pulay beta=0.05, nmaxold=8, weight=100",
    "convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    "maxiter": 160,
    "mpi_ranks": 4,
    "geometry": "fixed; single point only",
    "kpoint_note": "12x12x2 approximately matches the reciprocal spacing of the existing 2x2x1 Ti3SiC2 control at 6x6x2; it is not a convergence claim.",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if world.size != METHOD["mpi_ranks"]:
        raise SystemExit(f"Expected {METHOD['mpi_ranks']} MPI ranks, received {world.size}")
    launch_status = np.zeros(1, dtype=np.int32)
    if world.rank == 0:
        if OUT.exists():
            launch_status[0] = 1
        else:
            OUT.mkdir(parents=True, exist_ok=False)
    world.broadcast(launch_status, 0)
    if launch_status[0]:
        raise SystemExit(f"Refusing to overwrite existing exploratory run directory: {OUT}")
    world.barrier()

    import _gpaw
    import gpaw_data
    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer
    from gpaw.occupations import FermiDirac

    if not hasattr(_gpaw, "scalapack_diagonalize_dc"):
        raise SystemExit("Required ScaLAPACK diagonalizer is unavailable")
    os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())

    input_hash = sha256(SOURCE)
    atoms = read(SOURCE, format="cif")
    if len(atoms) != 12 or atoms.get_chemical_formula() != "C4Si2Ti6":
        raise SystemExit("Unexpected source-cell contents")
    if not atoms.pbc.all():
        raise SystemExit("Expected a fully periodic bulk cell")
    if any(key in atoms.info for key in ("REF_energy", "PW_PBE_energy_eV", "energy")):
        raise SystemExit("Source geometry unexpectedly contains an energy label")

    started = time.time()
    progress_path = OUT / "progress.json"
    if world.rank == 0:
        progress_path.write_text(
            json.dumps(
                {
                    "status": "running",
                    "input_sha256": input_hash,
                    "method": METHOD,
                    "dataset_role": "exploratory_only_not_training_or_validation",
                    "started_utc_epoch_s": started,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    calc = GPAW(
        mode=PW(METHOD["cutoff_eV"]),
        xc=METHOD["xc"],
        kpts=tuple(METHOD["kpts"]),
        occupations=FermiDirac(METHOD["smearing_eV"]),
        mixer=Mixer(0.05, 8, 100),
        parallel={"sl_auto": True},
        convergence=METHOD["convergence"],
        maxiter=METHOD["maxiter"],
        txt=str(OUT / "gpaw.log"),
    )
    atoms.calc = calc
    scf_converged = {"value": False}

    def mark_converged() -> None:
        scf_converged["value"] = True

    calc.hooks["converged"] = mark_converged

    def update_progress(ctx) -> None:
        niter = int(ctx.niter)
        if world.rank == 0 and niter % 10 == 0:
            progress_path.write_text(
                json.dumps(
                    {
                        "status": "running",
                        "input_sha256": input_hash,
                        "iteration": niter,
                        "elapsed_s": time.time() - started,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )

    calc.hooks["scf_step"] = update_progress
    energy_native = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(calc.results["free_energy"])
    forces = np.asarray(atoms.get_forces(), dtype=float)
    iterations = int(calc.get_number_of_iterations())
    elapsed_s = time.time() - started
    if not np.isfinite(energy_native) or not np.isfinite(free_energy) or not np.isfinite(forces).all():
        raise SystemExit("Nonfinite GPAW energy or force")

    if world.rank == 0:
        output_atoms = atoms.copy()
        output_atoms.calc = None
        output_atoms.info.update(
            {
                "config_type": "Bai2010_Ti3SiC2_12atom_exploratory_single_point",
                "source_doi": "10.1016/j.solidstatesciences.2010.03.007",
                "source_cif_sha256": input_hash,
                "dataset_role": "exploratory_only_not_training_or_validation",
                "training_label_eligible": False,
                "DFT_label_present": True,
                "energy_native_eV": energy_native,
                "free_energy_eV": free_energy,
                "energy_convention": "GPAW native extrapolated energy and free_energy both retained separately",
                "scf_converged": bool(scf_converged["value"]),
            }
        )
        output_atoms.arrays["DFT_forces_eV_A"] = forces.copy()
        write(OUT / "result.extxyz", output_atoms, format="extxyz")
        summary = {
            "status": "SCF_CONVERGED_EXPLORATORY_ONLY" if scf_converged["value"] else "SCF_NOT_CONVERGED",
            "input_file": SOURCE.name,
            "input_sha256": input_hash,
            "atoms": len(atoms),
            "formula": atoms.get_chemical_formula(),
            "cell_lengths_A": [float(x) for x in atoms.cell.lengths()],
            "cell_angles_deg": [float(x) for x in atoms.cell.angles()],
            "volume_A3": float(atoms.get_volume()),
            "energy_native_eV_cell": energy_native,
            "free_energy_eV_cell": free_energy,
            "energy_per_atom_native_eV": energy_native / len(atoms),
            "free_energy_per_atom_eV": free_energy / len(atoms),
            "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
            "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
            "scf_iterations": iterations,
            "scf_converged": bool(scf_converged["value"]),
            "elapsed_s": elapsed_s,
            "gpaw_log_raw_sha256": sha256(OUT / "gpaw.log"),
            "method": METHOD,
            "role": "exploratory_only_not_training_or_validation",
            "limitations": [
                "The Bai2010 numerical parameters are not yet independently checked against the paper full text.",
                "One fixed-geometry point does not establish k-point or cutoff convergence.",
                "This result is not included in any train/dev/test dataset or B1 role ledger.",
            ],
        }
        (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        raw_log_path = OUT / "gpaw.log"
        raw_log_sha = sha256(raw_log_path)
        with raw_log_path.open("rb") as source_log, (OUT / "gpaw.log.gz").open("wb") as compressed_log:
            with gzip.GzipFile(fileobj=compressed_log, mode="wb", compresslevel=9, mtime=0, filename="") as gzip_log:
                shutil.copyfileobj(source_log, gzip_log)
        raw_log_path.unlink()
        (OUT / "gpaw.log.raw.sha256").write_text(f"{raw_log_sha}  gpaw.log\n", encoding="utf-8")
        progress_path.write_text(
            json.dumps(
                {
                    "status": summary["status"],
                    "input_sha256": input_hash,
                    "iterations": iterations,
                    "elapsed_s": elapsed_s,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        files = sorted(path for path in OUT.iterdir() if path.is_file() and path.name != "SHA256SUMS.txt")
        (OUT / "SHA256SUMS.txt").write_text(
            "".join(f"{sha256(path)}  {path.name}\n" for path in files), encoding="utf-8"
        )
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    world.barrier()
    if not scf_converged["value"]:
        raise SystemExit("SCF did not converge; see gpaw.log and progress.json")


if __name__ == "__main__":
    main()
