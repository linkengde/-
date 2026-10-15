#!/usr/bin/env python3
"""Compute isolated elemental solid reference frames for v10 using GPAW PW-PBE."""
from pathlib import Path
import json, os, time, re
import numpy as np
import gpaw_data
from ase.build import bulk
from ase.io import write

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "frames"
OUT.mkdir(exist_ok=True)
os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
from gpaw import GPAW, PW
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

specs = [
    ("Ag_fcc", "Ag", "fcc", {"a": 4.152, "cubic": True}, (6, 6, 6)),
    ("Ti_hcp", "Ti", "hcp", {"a": 2.94, "c": 4.67}, (12, 12, 8)),
    ("Si_diamond", "Si", "diamond", {"a": 5.47, "cubic": True}, (8, 8, 8)),
    ("C_diamond", "C", "diamond", {"a": 3.57, "cubic": True}, (10, 10, 10)),
]
records = []
for name, symbol, crystal, kwargs, kpts in specs:
    atoms = bulk(symbol, crystal, **kwargs)
    atoms.pbc = True
    start = time.time()
    calc = GPAW(
        mode=PW(500), xc="PBE", kpts=kpts,
        occupations=FermiDirac(0.1), mixer=Mixer(0.05, 8, 100),
        convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
        maxiter=160, txt=str(ROOT / f"{name}.gpaw.log"),
    )
    atoms.calc = calc
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    atoms.info["REF_energy"] = energy
    atoms.arrays["REF_forces"] = forces
    atoms.info["config_type"] = f"v10_element_reference_{name}"
    atoms.info["dft_method"] = "GPAW PW-PBE 500 eV; Fermi smearing 0.1 eV"
    atoms.info["energy_reference_note"] = "pure elemental solid training anchor; not an isolated-atom chemical potential"
    atoms.calc = None
    write(OUT / f"{name}.extxyz", atoms, format="extxyz")
    log_text = (ROOT / f"{name}.gpaw.log").read_text()
    scf_match = re.findall(r"Converged in (\d+) steps", log_text)
    rec = {
        "name": name, "element": symbol, "crystal": crystal,
        "atoms": len(atoms), "formula": atoms.get_chemical_formula(),
        "cell_A": atoms.cell.array.tolist(), "volume_A3": float(atoms.get_volume()),
        "kpts": list(kpts), "cutoff_eV": 500, "xc": "PBE", "smearing_eV": 0.1,
        "energy_eV_cell": energy, "energy_eV_atom": energy / len(atoms),
        "scf_converged": bool(scf_match), "scf_iterations": int(scf_match[-1]) if scf_match else None,
        "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
        "rms_force_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
        "stress_eV_A3": stress.tolist(), "wall_seconds": time.time() - start,
        "reference_role": "training-only independent composition anchor for four-species energy gauge",
    }
    records.append(rec)
    (ROOT / "progress.json").write_text(json.dumps({"completed": records}, indent=2) + "\n")
    print(json.dumps(rec), flush=True)
    if not rec["scf_converged"]:
        raise RuntimeError(f"SCF did not converge for {name}")

summary = {
    "status": "COMPLETE",
    "method": "GPAW 26.7.0, PAW-PBE, plane waves 500 eV, Fermi smearing 0.1 eV, convergence energy/density/eigenstates 1e-5",
    "energy_convention": "ASE atoms.get_potential_energy(), matching the existing v9 interface labels",
    "note": "Crystalline elemental anchors span the four species composition directions. They make the energy-composition design matrix identifiable; they do not by themselves validate interface binding or claim isolated-atom formation energies.",
    "records": records,
}
(ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"status": summary["status"], "anchors": len(records)}, indent=2), flush=True)
