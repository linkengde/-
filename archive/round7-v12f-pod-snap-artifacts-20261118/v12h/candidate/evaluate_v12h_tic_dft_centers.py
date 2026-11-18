#!/workspace/venvs/dft-labels/bin/python
"""Compare v12h POD/SNAP with the three converged Ti-C DFT center labels."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
ACTIVE = Path(
    "/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/active_tic_hotspot_expansion"
)
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
ELEMENTS = ["Ag", "C", "Si", "Ti"]
MASS_LINES = "mass 1 107.8682\nmass 2 12.011\nmass 3 28.085\nmass 4 47.867"

CASES = [
    ("initial_pair1", "tic_pair_initial.extxyz", "tic_pair_initial_pw350_gamma.extxyz"),
    ("cold_10fs_pair1", "tic_pair_cold_10fs.extxyz", "tic_pair_cold_10fs_pw350_gamma.extxyz"),
    ("cold_10fs_pair2", "tic_pair2_cold_10fs.extxyz", "tic_pair2_cold_10fs_pw350_gamma.extxyz"),
]
PAIR_FILES = {
    "POD_v12h": (
        "pair_style pod\n"
        f"pair_coeff * * {ROOT/'POD/AgTiSiC_POD_v12h.pod'} "
        f"{ROOT/'POD/pod_v12h_coefficients.pod'} Ag C Si Ti"
    ),
    "SNAP_v12h": (
        "pair_style snap\n"
        f"pair_coeff * * {ROOT/'SNAP/AgTiSiC_SNAP_v12h.snapcoeff'} "
        f"{ROOT/'SNAP/AgTiSiC_SNAP_v12h.snapparam'} Ag C Si Ti"
    ),
}

results = []
for case, candidate_name, label_name in CASES:
    candidate = read(ACTIVE / "candidates" / candidate_name)
    dft = read(ACTIVE / "labels" / label_name)
    source_ids = np.asarray(candidate.arrays["source_atom_id"], dtype=int)
    label_ids = np.asarray(dft.arrays["source_atom_id"], dtype=int)
    if not np.array_equal(source_ids, label_ids):
        raise RuntimeError(f"atom order/source IDs differ for {case}")
    if not np.allclose(candidate.positions, dft.positions, atol=1.0e-8, rtol=0.0):
        raise RuntimeError(f"coordinates differ between candidate and DFT label for {case}")
    center_ids = [int(candidate.info["center_atom_id"]), int(candidate.info["partner_atom_id"])]
    idxs = [int(np.flatnonzero(source_ids == atom_id)[0]) for atom_id in center_ids]
    dft_forces = np.asarray(dft.arrays["REF_forces"], dtype=float)
    unit = candidate.copy()
    unit.set_pbc(True)
    out = ROOT / "tic_center_validation" / case
    out.mkdir(parents=True, exist_ok=True)
    data = out / "local.data"
    write(data, unit, format="lammps-data", atom_style="atomic", specorder=ELEMENTS)
    for model, pair in PAIR_FILES.items():
        run_dir = out / model
        run_dir.mkdir(exist_ok=True)
        dump = run_dir / "forces.dump"
        inp = run_dir / "run.in"
        inp.write_text(
            f"""clear
units metal
atom_style atomic
boundary p p p
read_data {data}
{MASS_LINES}
{pair}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 1
thermo_style custom step atoms pe
run 0
write_dump all custom {dump} id type fx fy fz modify sort id
"""
        )
        subprocess.run(
            [str(LMP), "-in", str(inp), "-log", str(run_dir / "log.lammps"), "-screen", "none"],
            cwd=run_dir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        rows = dump.read_text().splitlines()
        h = next(i for i, row in enumerate(rows) if row.startswith("ITEM: ATOMS "))
        cols = rows[h].split()[2:]
        force_cols = [cols.index(k) for k in ("fx", "fy", "fz")]
        model_forces = np.asarray(
            [[float(x) for x in row.split()] for row in rows[h + 1 : h + 1 + len(unit)]]
        )[:, force_cols]
        center_results = []
        for atom_id, idx in zip(center_ids, idxs):
            ref = dft_forces[idx]
            pred = model_forces[idx]
            center_results.append(
                {
                    "source_atom_id": atom_id,
                    "symbol": candidate[idx].symbol,
                    "DFT_force_eV_A": ref.tolist(),
                    "model_force_eV_A": pred.tolist(),
                    "vector_error_eV_A": float(np.linalg.norm(pred - ref)),
                }
            )
        # Relative mass-weighted net-force projection along Ti -> C.
        ti_idx = next(i for i in idxs if candidate[i].symbol == "Ti")
        c_idx = next(i for i in idxs if candidate[i].symbol == "C")
        direction = candidate.positions[c_idx] - candidate.positions[ti_idx]
        direction /= np.linalg.norm(direction)
        masses = candidate.get_masses()

        def relative_projection(force_array: np.ndarray) -> float:
            return float(
                np.dot(force_array[c_idx] / masses[c_idx] - force_array[ti_idx] / masses[ti_idx], direction)
            )

        results.append(
            {
                "case": case,
                "model": model,
                "pair_distance_A": float(np.linalg.norm(candidate.positions[c_idx] - candidate.positions[ti_idx])),
                "center_forces": center_results,
                "DFT_relative_projection_eV_A_amu": relative_projection(dft_forces),
                "model_relative_projection_eV_A_amu": relative_projection(model_forces),
            }
        )

summary = {
    "status": "EVALUATED",
    "role": "targeted local validation; three DFT-labeled Ti-C crops used in v12h training; not independent blind test",
    "models": results,
}
out_path = ROOT / "tic_center_validation" / "summary.json"
out_path.parent.mkdir(exist_ok=True)
out_path.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
