#!/workspace/venvs/dft-labels/bin/python
"""Evaluate POD/SNAP v12l against the held-out Si01 center DFT force."""
from __future__ import annotations

import json
import subprocess
import numpy as np
from pathlib import Path
from ase.io import read

ROOT = Path(__file__).resolve().parent
MODEL = ROOT
HOLDOUT = MODEL / "holdout_si_hotspot_104359.extxyz"
DATA = MODEL / "holdout_label_evaluation/POD_v12f_si01/local.data"
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
ELEMENTS = "Ag C Si Ti"
MASS_LINES = "mass 1 107.8682\nmass 2 12.011\nmass 3 28.085\nmass 4 47.867"

atoms = read(HOLDOUT)
center_id = int(atoms.info["center_atom_id"])
source_ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
idx = np.flatnonzero(source_ids == center_id)
if len(idx) != 1:
    raise RuntimeError("Could not map unique center atom in held-out crop")
lmp_center_id = int(idx[0]) + 1
dft_force = np.asarray(atoms.info["DFT_center_force_eV_A"], dtype=float)
if dft_force.shape != (3,) or not np.isfinite(dft_force).all():
    raise RuntimeError("Invalid held-out DFT center force")

models = {
    "POD_v12l": (
        f"pair_style pod\npair_coeff * * {ROOT/'POD/AgTiSiC_POD_v12l.pod'} "
        f"{ROOT/'POD/pod_v12l_coefficients.pod'} {ELEMENTS}"
    ),
    "SNAP_v12l": (
        f"pair_style snap\npair_coeff * * {ROOT/'SNAP/AgTiSiC_SNAP_v12l.snapcoeff'} "
        f"{ROOT/'SNAP/AgTiSiC_SNAP_v12l.snapparam'} {ELEMENTS}"
    ),
}
results = []
for name, pair in models.items():
    out = ROOT / "holdout_si01" / name
    out.mkdir(parents=True, exist_ok=True)
    local_data = out / "local.data"
    if not local_data.exists():
        local_data.write_bytes(DATA.read_bytes())
    dump = out / "forces.dump"
    inp = out / "run.in"
    inp.write_text(f"""clear
units metal
atom_style atomic
boundary p p p
read_data {local_data}
{MASS_LINES}
{pair}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 1
thermo_style custom step atoms pe
run 0
write_dump all custom {dump} id type fx fy fz modify sort id
""")
    subprocess.run([str(LMP), "-in", str(inp), "-log", str(out/"log.lammps"),
                    "-screen", "none"], cwd=out, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    lines = dump.read_text().splitlines()
    head = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS "))
    cols = lines[head].split()[2:]
    table = {}
    for row in lines[head+1:]:
        if row.startswith("ITEM:"):
            break
        vals = row.split()
        table[int(vals[cols.index("id")])] = np.array([float(vals[cols.index(k)]) for k in ("fx", "fy", "fz")])
    if lmp_center_id not in table:
        raise RuntimeError(f"Center LAMMPS ID {lmp_center_id} missing from {name} dump")
    force = table[lmp_center_id]
    err = float(np.linalg.norm(force - dft_force))
    cosine = float(np.dot(force, dft_force) / (np.linalg.norm(force)*np.linalg.norm(dft_force)))
    results.append({"model": name, "center_atom_id": center_id,
                    "lammps_atom_id": lmp_center_id,
                    "DFT_force_eV_A": dft_force.tolist(),
                    "model_force_eV_A": force.tolist(),
                    "vector_error_eV_A": err,
                    "direction_cosine": cosine,
                    "role": "development holdout; excluded from training"})

summary = {"status": "EVALUATED", "holdout_center_atom_id": center_id,
           "holdout_crop_atoms": len(atoms), "models": results}
(ROOT / "holdout_si01/summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
