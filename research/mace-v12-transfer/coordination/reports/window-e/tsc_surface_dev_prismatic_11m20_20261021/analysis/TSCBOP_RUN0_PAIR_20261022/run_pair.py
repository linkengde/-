#!/usr/bin/env python3
"""Frozen-coordinate TSCBOP/LAMMPS run-0 comparison for a DFT pair."""
from pathlib import Path
import csv
import hashlib
import json
import os
import subprocess
import time

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
POT = Path("/workspace/-/research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e5_force_jacobian_20261014/potentials/TSCBOP.source.tersoff")
POT_SHA = "729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7"
LMP = Path("/workspace/.venvs/lammps-run0/bin/lmp")
TASKS = [
    ("TSC120_DEV_PRISM_11m20_COD9009647_10L_REF", "REF", "attempt-0003"),
    ("TSC120_DEV_PRISM_11m20_COD9009647_10L_RESP_TIZ020", "RESP_TIZ020", "attempt-0001"),
]
TYPE = {"Ti": 1, "Si": 2, "C": 3}
MASS = {"Ti": 47.867, "Si": 28.085, "C": 12.011}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def geometry_sha(atoms):
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes())
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes())
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def make_lammps_input(atoms, run_dir):
    run_dir.mkdir(parents=True, exist_ok=True)
    cell = np.asarray(atoms.cell, float)
    if atoms.pbc.tolist() != [True, True, True]:
        raise ValueError("expected p p p slab box")
    if max(abs(cell[0, 1]), abs(cell[0, 2]), abs(cell[1, 2])) > 1e-12:
        raise ValueError("cell is not in the restricted lower-triangular convention")
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    frac = np.asarray(atoms.positions) @ np.linalg.inv(cell)
    images = np.floor(frac + 1e-12).astype(int)
    wrapped = (frac - images) @ cell
    ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
    symbols = atoms.get_chemical_symbols()
    if len(set(ids.tolist())) != len(atoms) or ids.min() < 1:
        raise ValueError("source atom IDs must be unique positive integers")
    if not set(symbols).issubset(TYPE):
        raise ValueError(f"unsupported element(s): {sorted(set(symbols)-set(TYPE))}")

    lines = [
        "Frozen-coordinate TSCBOP diagnostic; source IDs retained", "",
        f"{len(atoms)} atoms", "3 atom types", "",
        f"0 {lx:.16g} xlo xhi", f"0 {ly:.16g} ylo yhi", f"0 {lz:.16g} zlo zhi",
        f"{xy:.16g} {xz:.16g} {yz:.16g} xy xz yz", "", "Masses", "",
    ]
    for element, atom_type in TYPE.items():
        lines.append(f"{atom_type} {MASS[element]}")
    lines += ["", "Atoms # atomic", ""]
    with open(run_dir / "atoms.data", "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
        for atom_id, element, xyz in zip(ids, symbols, wrapped):
            f.write(f"{atom_id} {TYPE[element]} {xyz[0]:.16g} {xyz[1]:.16g} {xyz[2]:.16g}\n")
    with open(run_dir / "coordinate_map.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["source_atom_id", "element", "x_source_A", "y_source_A", "z_source_A", "x_wrapped_A", "y_wrapped_A", "z_wrapped_A"])
        for atom_id, element, pos, xyz in zip(ids, symbols, atoms.positions, wrapped):
            w.writerow([int(atom_id), element, *map(float, pos), *map(float, xyz)])
    (run_dir / "TSCBOP.tersoff").write_bytes(POT.read_bytes())
    (run_dir / "run.in").write_text(
        "clear\nunits metal\natom_style atomic\nboundary p p p\n"
        f"read_data {(run_dir / 'atoms.data').resolve()}\n"
        "pair_style tersoff\n"
        f"pair_coeff * * {(run_dir / 'TSCBOP.tersoff').resolve()} Ti Si C\n"
        "neighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\n"
        "thermo 1\nthermo_style custom step atoms pe\n"
        "dump forces all custom 1 forces.dump id type x y z fx fy fz\n"
        'dump_modify forces sort id format line "%d %d %.16g %.16g %.16g %.16g %.16g %.16g"\n'
        "run 0 post no\nundump forces\n"
    )
    return ids, symbols, wrapped


def parse_forces_dump(path):
    lines = path.read_text().splitlines()
    section = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    cols = lines[section].split()[2:]
    col = {name: i for i, name in enumerate(cols)}
    inv_type = {v: k for k, v in TYPE.items()}
    records = {}
    for line in lines[section + 1:]:
        if not line.strip():
            continue
        fields = line.split()
        atom_id = int(fields[col["id"]])
        coord_force = np.array([float(fields[col[key]]) for key in ("x", "y", "z", "fx", "fy", "fz")])
        if not np.isfinite(coord_force).all():
            raise ValueError(f"non-finite LAMMPS record at source ID {atom_id}")
        records[atom_id] = (inv_type[int(fields[col["type"]])], coord_force[:3], coord_force[3:])
    return records


def metrics(error, ids, mask):
    selected = error[mask]
    norms = np.linalg.norm(selected, axis=1)
    selected_ids = ids[mask]
    return {
        "atom_count": int(mask.sum()),
        "component_rmse_eV_A": float(np.sqrt(np.mean(selected ** 2))),
        "vector_rmse_eV_A": float(np.sqrt(np.mean(np.sum(selected ** 2, axis=1)))),
        "p95_vector_error_eV_A": float(np.percentile(norms, 95)),
        "max_vector_error_eV_A": float(norms.max()),
        "max_error_source_atom_id": int(selected_ids[np.argmax(norms)]),
    }


def main():
    if sha(POT) != POT_SHA:
        raise ValueError("frozen TSCBOP source SHA mismatch")
    if not LMP.is_file():
        raise FileNotFoundError(LMP)
    version_lines = subprocess.run([str(LMP), "-h"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout.splitlines()
    version = next(line.strip() for line in version_lines if line.strip())
    results = {}
    for structure_id, suffix, attempt in TASKS:
        inp = ROOT / "inputs" / f"{structure_id}.extxyz"
        result_dir = ROOT / "results" / structure_id / attempt
        result = json.loads((result_dir / "RESULT.json").read_text())
        status = json.loads((result_dir / "RUN_STATUS.json").read_text())
        if status.get("status") != "SUCCESS" or status.get("scf_converged") is not True:
            raise ValueError(f"{structure_id}: DFT is not converged")
        atoms = read(inp)
        labelled = read(result_dir / "labelled.extxyz")
        if not all(np.array_equal(getattr(atoms, key), getattr(labelled, key)) for key in ("numbers", "positions", "cell", "pbc")):
            raise ValueError(f"{structure_id}: labelled frame geometry differs from input")
        ids = np.asarray(atoms.arrays.get("source_atom_id"))
        labelled_ids = np.asarray(labelled.arrays.get("source_atom_id"))
        if ids.shape != (len(atoms),) or not np.array_equal(ids, labelled_ids) or not np.array_equal(atoms.numbers, labelled.numbers):
            raise ValueError(f"{structure_id}: source ID mapping missing or changed")
        if sha(inp) != result["input_sha256"] or geometry_sha(atoms) != result["geometry_sha256"]:
            raise ValueError(f"{structure_id}: input or geometry SHA mismatch")
        dft_forces = np.asarray(result["forces_eV_per_angstrom"], dtype=float)
        if dft_forces.shape != (len(atoms), 3) or not np.isfinite(dft_forces).all():
            raise ValueError(f"{structure_id}: DFT force labels invalid or array order mismatch")

        run_dir = OUT / "runs" / suffix
        run_ids, symbols, wrapped = make_lammps_input(atoms, run_dir)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OMP_DYNAMIC"] = "FALSE"
        started = time.monotonic()
        proc = subprocess.run([str(LMP), "-in", "run.in", "-log", "lammps.log", "-screen", "screen.log", "-nocite"], cwd=run_dir, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        elapsed = time.monotonic() - started
        (run_dir / "launcher.log").write_text(proc.stdout)
        if proc.returncode != 0 or not (run_dir / "forces.dump").is_file():
            raise RuntimeError(f"{structure_id}: LAMMPS run 0 failed, exit {proc.returncode}")
        lmp = parse_forces_dump(run_dir / "forces.dump")
        if set(lmp) != set(run_ids.tolist()):
            raise ValueError(f"{structure_id}: LAMMPS atom IDs differ")
        by_id = {int(atom_id): idx for idx, atom_id in enumerate(run_ids)}
        lmp_forces = np.zeros_like(dft_forces)
        max_coordinate_error = 0.0
        for atom_id in run_ids:
            idx = by_id[int(atom_id)]
            element, xyz, force = lmp[int(atom_id)]
            if element != symbols[idx]:
                raise ValueError(f"{structure_id}: element mismatch at source ID {atom_id}")
            err = float(np.linalg.norm(xyz - wrapped[idx]))
            max_coordinate_error = max(max_coordinate_error, err)
            if err > 1e-8:
                raise ValueError(f"{structure_id}: coordinate mismatch at source ID {atom_id}: {err}")
            lmp_forces[idx] = force
        error = lmp_forces - dft_forces
        with open(OUT / f"{suffix}_FORCE_COMPARISON.csv", "w", newline="") as f:
            cols = ["source_atom_id", "element", "x_A", "y_A", "z_A", "dft_fx_eV_A", "dft_fy_eV_A", "dft_fz_eV_A", "tscbop_fx_eV_A", "tscbop_fy_eV_A", "tscbop_fz_eV_A", "error_fx_eV_A", "error_fy_eV_A", "error_fz_eV_A", "force_error_norm_eV_A"]
            w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
            w.writeheader()
            for i, atom_id in enumerate(run_ids):
                row = {"source_atom_id": int(atom_id), "element": symbols[i], "x_A": atoms.positions[i, 0], "y_A": atoms.positions[i, 1], "z_A": atoms.positions[i, 2], "force_error_norm_eV_A": float(np.linalg.norm(error[i]))}
                row.update({f"dft_f{axis}_eV_A": dft_forces[i, k] for k, axis in enumerate("xyz")})
                row.update({f"tscbop_f{axis}_eV_A": lmp_forces[i, k] for k, axis in enumerate("xyz")})
                row.update({f"error_f{axis}_eV_A": error[i, k] for k, axis in enumerate("xyz")})
                w.writerow(row)
        results[suffix] = {
            "structure_id": structure_id,
            "atoms": atoms,
            "ids": run_ids,
            "symbols": np.asarray(symbols),
            "dft_forces": dft_forces,
            "lammps_forces": lmp_forces,
            "error": error,
            "input_sha256": sha(inp),
            "geometry_sha256": result["geometry_sha256"],
            "config_sha256": result["config_sha256"],
            "scf_steps": result["scf_steps"],
            "run_input_sha256": sha(run_dir / "run.in"),
            "data_sha256": sha(run_dir / "atoms.data"),
            "coordinate_map_sha256": sha(run_dir / "coordinate_map.csv"),
            "max_coordinate_error_A": max_coordinate_error,
            "lammps_runtime_s": elapsed,
        }

    ref, resp = results["REF"], results["RESP_TIZ020"]
    if not np.array_equal(ref["ids"], resp["ids"]) or not np.array_equal(ref["symbols"], resp["symbols"]):
        raise ValueError("reference and response ID/element sequences differ")
    ref_atoms, resp_atoms = ref["atoms"], resp["atoms"]
    if not np.array_equal(ref_atoms.cell, resp_atoms.cell) or not np.array_equal(ref_atoms.pbc, resp_atoms.pbc):
        raise ValueError("reference and response cell/PBC differ")
    moved = np.linalg.norm(resp_atoms.positions - ref_atoms.positions, axis=1) > 1e-7
    moved_ids = ref["ids"][moved]
    if moved_ids.tolist() != [4, 113]:
        raise ValueError(f"expected only source IDs [4, 113] to move; found {moved_ids.tolist()}")
    displacements = resp_atoms.positions[moved] - ref_atoms.positions[moved]
    if not np.allclose(np.sort(displacements[:, 2]), [-0.02, 0.02], atol=1e-8) or not np.allclose(displacements[:, :2], 0, atol=1e-10):
        raise ValueError("displacement is not the prescribed +/-0.02 A z-pair")

    response_error = (resp["lammps_forces"] - ref["lammps_forces"]) - (resp["dft_forces"] - ref["dft_forces"])
    dft_delta_forces = resp["dft_forces"] - ref["dft_forces"]
    lammps_delta_forces = resp["lammps_forces"] - ref["lammps_forces"]
    # A local response subset uses the 3.3 A first-neighbour neighborhood around either moved terminal Ti.
    moved_indices = np.flatnonzero(moved)
    dists = ref_atoms.get_distances(moved_indices[0], np.arange(len(ref_atoms)), mic=True)
    dists2 = ref_atoms.get_distances(moved_indices[1], np.arange(len(ref_atoms)), mic=True)
    near = np.minimum(dists, dists2) < 3.3
    with open(OUT / "RESPONSE_DIFFERENCE.csv", "w", newline="") as f:
        cols = ["source_atom_id", "element", "dft_delta_fx_eV_A", "dft_delta_fy_eV_A", "dft_delta_fz_eV_A", "tscbop_delta_fx_eV_A", "tscbop_delta_fy_eV_A", "tscbop_delta_fz_eV_A", "response_error_norm_eV_A", "within_3p3A_of_displaced_Ti"]
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for i, atom_id in enumerate(ref["ids"]):
            dft_delta = resp["dft_forces"][i] - ref["dft_forces"][i]
            lmp_delta = resp["lammps_forces"][i] - ref["lammps_forces"][i]
            row = {"source_atom_id": int(atom_id), "element": ref["symbols"][i], "response_error_norm_eV_A": float(np.linalg.norm(response_error[i])), "within_3p3A_of_displaced_Ti": bool(near[i])}
            row.update({f"dft_delta_f{axis}_eV_A": dft_delta[k] for k, axis in enumerate("xyz")})
            row.update({f"tscbop_delta_f{axis}_eV_A": lmp_delta[k] for k, axis in enumerate("xyz")})
            w.writerow(row)

    summary = {
        "scope": "Two frozen-coordinate LAMMPS run 0 calculations only; no minimization, MD, DFT, or parameter changes.",
        "role": "same COD 9009647 parent family; geometry-family diagnostic only; no TRAIN/DEV/TEST role assignment",
        "potential": {"source_path": str(POT), "sha256": POT_SHA, "pair_style": "tersoff", "pair_coeff_element_order": ["Ti", "Si", "C"], "units": "metal; energies eV and forces eV/Angstrom"},
        "lammps": {"executable": str(LMP), "version": version, "OMP_NUM_THREADS": 1, "steps": 0},
        "structures": {},
        "response_error": {},
        "response_definition": "(F_TSCBOP_RESP-F_TSCBOP_REF) - (F_DFT_RESP-F_DFT_REF)",
        "local_neighborhood_definition": "reference-geometry minimum-image distance <3.3 A from either displaced terminal Ti source ID",
    }
    for suffix, row in results.items():
        all_atoms = np.ones(len(row["ids"]), dtype=bool)
        summary["structures"][suffix] = {
            "structure_id": row["structure_id"], "input_sha256": row["input_sha256"], "geometry_sha256": row["geometry_sha256"], "config_sha256": row["config_sha256"],
            "scf_steps": row["scf_steps"], "atom_count": len(row["ids"]), "max_lammps_vs_wrapped_coordinate_error_A": row["max_coordinate_error_A"],
            "lammps_run0_runtime_s": row["lammps_runtime_s"], "run_input_sha256": row["run_input_sha256"], "lammps_data_sha256": row["data_sha256"], "coordinate_map_sha256": row["coordinate_map_sha256"],
            "all": metrics(row["error"], row["ids"], all_atoms),
            "by_element": {element: metrics(row["error"], row["ids"], row["symbols"] == element) for element in ("Ti", "Si", "C")},
        }
    summary["response_error"] = {
        "displaced_source_atom_ids": moved_ids.tolist(),
        "displacements_A": {str(int(i)): (resp_atoms.positions[j] - ref_atoms.positions[j]).tolist() for j, i in enumerate(ref["ids"]) if moved[j]},
        "all": metrics(response_error, ref["ids"], np.ones(len(ref["ids"]), dtype=bool)),
        "within_3p3A_of_displaced_Ti": metrics(response_error, ref["ids"], near),
        "near_atom_count": int(near.sum()),
        "dft_response_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(dft_delta_forces ** 2, axis=1)))),
        "tscbop_response_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(lammps_delta_forces ** 2, axis=1)))),
        "dft_response_vector_rms_within_3p3A_eV_A": float(np.sqrt(np.mean(np.sum(dft_delta_forces[near] ** 2, axis=1)))),
        "tscbop_response_vector_rms_within_3p3A_eV_A": float(np.sqrt(np.mean(np.sum(lammps_delta_forces[near] ** 2, axis=1)))),
        "by_element": {element: metrics(response_error, ref["ids"], ref["symbols"] == element) for element in ("Ti", "Si", "C")},
        "dft_response_vector_rms_by_element_eV_A": {element: float(np.sqrt(np.mean(np.sum(dft_delta_forces[ref["symbols"] == element] ** 2, axis=1)))) for element in ("Ti", "Si", "C")},
        "tscbop_response_vector_rms_by_element_eV_A": {element: float(np.sqrt(np.mean(np.sum(lammps_delta_forces[ref["symbols"] == element] ** 2, axis=1)))) for element in ("Ti", "Si", "C")},
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
