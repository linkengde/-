#!/usr/bin/env python3
"""Guarded GPAW/PBE fixed-geometry registry single-point input.

Default behavior only validates the frozen geometry and prints the locked plan.
A DFT calculation requires --run and a matching C-approved release JSON.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter
from importlib.metadata import version
from pathlib import Path
import sys
import numpy as np
from ase.io import read

REGISTRY_ID = "R1"
SOURCE_COMMIT = "9a5482d93a8756850ab457a605cf75dfa1f8a3d4"
SOURCE_PATH = "research/mace-v12-transfer/coordination/reports/window-c/ag_si_dev_v2_structure_generation/Ag_Ti3SiC2_Si_R1.xyz"
STRUCTURE_SHA256 = "b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890"
EXPECTED_CELL_A = np.array([[3.071, 0.0, 0.0], [-1.5355, 2.6595640150220112, 0.0], [0.0, 0.0, 71.26826600017428]])
EXPECTED_PBC = [True, True, False]
EXPECTED_COUNTS = {"Ag": 7, "C": 4, "Si": 3, "Ti": 6}
PAW_SHA256 = {
    "Ag": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
    "C": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
    "Si": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
    "Ti": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
}
MESHES = {"coarse": (6, 6, 1), "medium": (8, 8, 1)}
KPOINT_CONVENTION = "GPAW default Monkhorst-Pack integer mesh; no explicit gamma override or extra offset"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit("HOLD: " + message)


def load_frozen_structure():
    path = Path(__file__).resolve().parent / "structure.xyz"
    raw = path.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    require(sha == STRUCTURE_SHA256, "frozen structure SHA-256 mismatch")
    atoms = read(path, format="extxyz")
    require(len(atoms) == 20, "expected exactly 20 atoms")
    require(dict(sorted(Counter(atoms.get_chemical_symbols()).items())) == EXPECTED_COUNTS,
            "element composition mismatch")
    require(np.array_equal(atoms.pbc, EXPECTED_PBC), "PBC mismatch")
    require(np.allclose(atoms.cell.array, EXPECTED_CELL_A, rtol=0.0, atol=5e-12),
            "cell mismatch")
    require(np.isfinite(atoms.positions).all(), "non-finite coordinates")
    zmin = float(atoms.positions[:, 2].min())
    zmax = float(atoms.positions[:, 2].max())
    vacuum_bottom = zmin
    vacuum_top = float(atoms.cell.array[2, 2] - zmax)
    return atoms, raw, {
        "registry_id": REGISTRY_ID,
        "source_commit": SOURCE_COMMIT,
        "source_path": SOURCE_PATH,
        "structure_path": str(path),
        "structure_sha256": sha,
        "n_atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "element_counts": dict(sorted(Counter(atoms.get_chemical_symbols()).items())),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "z_vacuum_bottom_A": vacuum_bottom,
        "z_vacuum_top_A": vacuum_top,
        "vacuum_warning": min(vacuum_bottom, vacuum_top) < 18.0,
    }


def check_release(release_path: Path, expected_sha: str, mesh_name: str):
    raw = release_path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected_sha,
            "release config SHA-256 mismatch")
    cfg = json.loads(raw)
    mesh = list(MESHES[mesh_name])
    require(cfg.get("execution_authorized") is True, "C release does not authorize DFT execution")
    approval = cfg.get("C_approval", {})
    require(approval.get("signed_by") == "C" and bool(approval.get("approval_id")),
            "C approval identity is missing")
    expected = {
        "registry_id": REGISTRY_ID,
        "source_commit": SOURCE_COMMIT,
        "structure_sha256": STRUCTURE_SHA256,
        "xc": "PBE",
        "cutoff_eV": 500,
        "smearing_eV": 0.1,
        "kpoint_mesh": mesh,
        "kpoint_convention": KPOINT_CONVENTION,
        "fixed_geometry": True,
        "relax_atoms": False,
        "relax_cell": False,
        "save_energy": True,
        "save_forces": True,
        "save_stress": True,
    }
    for key, value in expected.items():
        require(cfg.get(key) == value, f"C release field mismatch: {key}")
    require(cfg.get("mpi_ranks", 0) > 0, "C release must lock a positive MPI rank count")
    return cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", choices=tuple(MESHES),
                        help="coarse=6x6x1, medium=8x8x1; required only with --run")
    parser.add_argument("--run", action="store_true",
                        help="execute GPAW only after a matching C release is supplied")
    parser.add_argument("--release-config", type=Path)
    parser.add_argument("--release-sha256")
    args = parser.parse_args()
    atoms, raw, preflight = load_frozen_structure()
    preflight["prepared_kpoint_meshes"] = {name: list(mesh) for name, mesh in MESHES.items()}
    preflight["kpoint_convention"] = KPOINT_CONVENTION
    preflight["status"] = "INPUT_VALIDATED_DFT_NOT_EXECUTED"
    preflight["method"] = {
        "xc": "PBE", "basis": "GPAW PAW plane waves", "cutoff_eV": 500,
        "occupations": "FermiDirac", "smearing_eV": 0.1,
        "scf_convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
        "max_scf_iterations": 160, "relax_atoms": False, "relax_cell": False,
    }
    if not args.run:
        print(json.dumps(preflight, indent=2))
        print("DFT_NOT_RUN")
        return
    require(args.mesh in MESHES, "select one mesh to execute")
    require(args.release_config is not None and args.release_sha256,
            "--run requires --release-config and --release-sha256")
    cfg = check_release(args.release_config, args.release_sha256, args.mesh)

    import ase
    import gpaw
    from gpaw import GPAW, PW, FermiDirac
    from gpaw.mixer import Mixer
    from gpaw.mpi import world
    from gpaw.setup_data import search_for_file
    require(gpaw.__version__ == cfg.get("gpaw_version") == "26.7.0", "GPAW version mismatch")
    require(ase.__version__ == cfg.get("ase_version") == "3.29.0", "ASE version mismatch")
    require(version("gpaw-data") == cfg.get("gpaw_data_version") == "1.2.1", "PAW data version mismatch")
    require(world.size == cfg["mpi_ranks"], "MPI rank count mismatch")
    paw_files = {}
    for symbol in sorted(EXPECTED_COUNTS):
        filename, _ = search_for_file(f"{symbol}.PBE")
        actual = hashlib.sha256(Path(filename).read_bytes()).hexdigest()
        require(actual == PAW_SHA256[symbol], f"PAW SHA-256 mismatch for {symbol}")
        paw_files[symbol] = {"path": filename, "sha256": actual}

    before_positions = atoms.positions.copy()
    before_cell = atoms.cell.array.copy()
    before_pbc = atoms.pbc.copy()
    mesh = MESHES[args.mesh]
    result_dir = Path(__file__).resolve().parent / "results" / args.mesh
    result_dir.mkdir(parents=True, exist_ok=False)
    log_path = result_dir / "gpaw.log"
    result_path = result_dir / "result.json"
    calc = GPAW(
        mode=PW(500), xc="PBE", kpts=mesh,
        occupations=FermiDirac(0.1),
        convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
        maxiter=160, mixer=Mixer(0.05, 8, 100),
        symmetry={"point_group": True, "time_reversal": True},
        poissonsolver={}, parallel={"sl_auto": True}, txt=str(log_path),
    )
    atoms.calc = calc
    energy = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False, apply_constraint=False), dtype=float)
    require(np.array_equal(atoms.positions, before_positions), "atomic coordinates changed")
    require(np.array_equal(atoms.cell.array, before_cell), "cell changed")
    require(np.array_equal(atoms.pbc, before_pbc), "PBC changed")
    require(np.isfinite(energy) and np.isfinite(free_energy) and np.isfinite(forces).all()
            and np.isfinite(stress).all(), "non-finite DFT labels")
    require(np.allclose(stress, stress.T, rtol=0.0, atol=1e-10), "stress tensor is not symmetric")
    require(bool(calc.scf.converged), "SCF did not converge")
    output = {
        "registry_id": REGISTRY_ID,
        "source_commit": SOURCE_COMMIT,
        "source_path": SOURCE_PATH,
        "structure_sha256": STRUCTURE_SHA256,
        "cell_A": atoms.cell.array.tolist(), "pbc": atoms.pbc.tolist(),
        "total_energy_eV": energy, "free_energy_eV": free_energy,
        "forces_eV_A": forces.tolist(), "stress_eV_A3": stress.tolist(),
        "stress_convention": "ASE tensile-positive Cartesian full 3x3",
        "stress_supervision_allowed": False,
        "kpoint_mesh": list(mesh), "scf_converged": bool(calc.scf.converged),
        "paw_files": paw_files, "C_release_sha256": args.release_sha256,
        "geometry_mutated": False,
        "stress_note": "Slab stress is normalized by the full vacuum-containing cell; do not use as bulk stress without C approval.",
    }
    with result_path.open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"status": "DFT_SINGLE_POINT_COMPLETE", "result": str(result_path)}, indent=2))


if __name__ == "__main__":
    main()
