#!/usr/bin/env python3
"""Create temporary same-phase FCC/BCC springs for the Ag-W v5 interface band.

These bonds are geometric preconditioning only. The final structure must be
minimized again with Dai EFS and no springs, and the final atomic-only data
file is exported separately.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
PACK = HERE.parents[2]
MODEL = PACK / "models/ag_w_qsgs_200A"


def read_atomic(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        lines = stream.readlines()
    section = next(i for i, line in enumerate(lines) if line.strip().startswith("Atoms")) + 1
    rows = []
    for line in lines[section:]:
        fields = line.split()
        if not fields:
            if rows:
                break
            continue
        if len(fields) < 5:
            if rows:
                break
            continue
        rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    array = np.asarray(rows, dtype=float)
    return array[:, 0].astype(np.int64), array[:, 1].astype(np.int8), array[:, 2:5].copy()


def write_atomic_data(path: Path, ids: np.ndarray, types: np.ndarray, pos: np.ndarray) -> None:
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as stream:
        stream.write("Ag-W QSGS v5 temporary same-phase spring network\n\n")
        stream.write(f"{len(pos)} atoms\n\n2 atom types\n\n")
        stream.write("0 200 xlo xhi\n0 200 ylo yhi\n0 300 zlo zhi\n\n")
        stream.write("Masses\n\n1 183.84 # W\n2 107.8682 # Ag\n\n")
        stream.write("Atoms # atomic\n\n")
        for atom_id, atom_type, xyz in zip(ids, types, pos):
            stream.write(f"{int(atom_id)} {int(atom_type)} {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interface-band-A", type=float, default=12.0)
    parser.add_argument("--spring-K-eV-A2", type=float, default=8.0)
    parser.add_argument("--cross-soft-eV", type=float, default=50.0)
    args = parser.parse_args()
    if min(args.interface_band_A, args.spring_K_eV_A2, args.cross_soft_eV) <= 0.0:
        parser.error("interface band, spring stiffness, and soft repulsion must be positive")
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory {out}")
    out.mkdir(parents=True, exist_ok=True)
    candidate, report_path = args.candidate.resolve(), args.report.resolve()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    candidate_hash = hashlib.sha256(candidate.read_bytes()).hexdigest()
    if candidate_hash != report["candidate"]["sha256"]:
        raise SystemExit("Candidate SHA-256 does not match its geometry report")
    ids, types, positions = read_atomic(candidate)

    sys.path.insert(0, str(PACK / "scripts"))
    sys.path.insert(0, str(HERE))
    import build_ag_w_qsgs as builder
    from reconstruct_interface_v5 import interpolate_grid, signed_distance_and_normal

    source_mask = np.load(MODEL / "qsgs_morphology.npz")["W_backbone"].astype(bool)
    sigma = float(report["parameters"]["mask_smoothing_vox"])
    if sigma:
        from scipy.ndimage import gaussian_filter
        occupancy = gaussian_filter(source_mask.astype(np.float32), sigma=sigma,
                                    mode=("wrap", "wrap", "nearest"))
        nvox = int(source_mask.sum())
        take = np.argsort(occupancy.ravel(), kind="stable")[-nvox:]
        mask = np.zeros(occupancy.size, dtype=bool)
        mask[take] = True
        mask = mask.reshape(source_mask.shape)
    else:
        mask = source_mask.copy()
    phi, _ = signed_distance_and_normal(mask, float(report["parameters"]["normal_smoothing_vox"]))
    w_base, w_strain = builder.make_lattice_sites("W")
    ag_base, ag_strain = builder.make_lattice_sites("Ag")
    w_base = w_base[builder.in_mask(mask, w_base)]
    ag_base = ag_base[~builder.in_mask(mask, ag_base)]
    if len(w_base) != int(np.count_nonzero(types == 1)) or len(ag_base) != int(np.count_nonzero(types == 2)):
        raise RuntimeError("Candidate lattice population does not match its report")
    base_positions = np.vstack((w_base, ag_base))
    phi_atom = interpolate_grid(phi, base_positions)
    in_band = np.abs(phi_atom) <= args.interface_band_A

    # Temporary same-phase Morse wells supply a short-range elastic network.
    # Pair interactions honor the periodic boundaries without explicit bonds.
    a_w = float(builder.LATTICE["W"][1])
    a_ag = float(builder.LATTICE["Ag"][1])
    sx_w, sy_w = (1.0 + np.asarray(w_strain, dtype=float)).tolist()
    sx_ag, sy_ag = (1.0 + np.asarray(ag_strain, dtype=float)).tolist()
    r0_w = float(a_w / 2.0 * np.sqrt(sx_w**2 + sy_w**2 + 1.0))
    r_xy = float(a_ag / 2.0 * np.sqrt(sx_ag**2 + sy_ag**2))
    r_xz = float(a_ag / 2.0 * np.sqrt(sx_ag**2 + 1.0))
    r_yz = float(a_ag / 2.0 * np.sqrt(sy_ag**2 + 1.0))
    r0_ag = float((4.0 * r_xy + 4.0 * r_xz + 4.0 * r_yz) / 12.0)
    morse_alpha = float(np.sqrt(args.spring_K_eV_A2 / 2.0))

    # Sort the temporary atom IDs into interface and far-field ranges.
    order = np.concatenate((np.flatnonzero(in_band), np.flatnonzero(~in_band)))
    new_ids = np.arange(1, len(order) + 1, dtype=np.int64)
    new_types, new_positions = types[order], positions[order]
    # Keep the minimizer seed infinitesimally inside the temporary z caps;
    # LAMMPS wall/harmonic rejects particles exactly on a wall plane.
    new_positions[:, 2] = np.clip(new_positions[:, 2], 1.0e-3, 200.0 - 1.0e-3)
    n_interface = int(in_band.sum())

    data_path = out / "interface_precondition.data.gz"
    write_atomic_data(data_path, new_ids, new_types, new_positions)
    input_path = out / "in.spring_precondition"
    input_text = f"""# Temporary collision-aware Ag-W interface precondition
units metal
atom_style atomic
boundary p p f
atom_modify sort 0 0
read_data {data_path}
mass 1 183.84
mass 2 107.8682
pair_style hybrid/overlay soft 3.0 morse 3.0
pair_coeff 1 1 soft {args.cross_soft_eV:.12g} 2.70
pair_coeff 1 1 morse 1.0 {morse_alpha:.12g} {r0_w:.12g} 3.0
pair_coeff 1 2 soft {args.cross_soft_eV:.12g} 2.75
pair_coeff 1 2 morse 0.0 {morse_alpha:.12g} 2.8 3.0
pair_coeff 2 2 soft {args.cross_soft_eV:.12g} 2.48
pair_coeff 2 2 morse 1.0 {morse_alpha:.12g} {r0_ag:.12g} 3.0
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
group interface id 1:{n_interface}
group far id {n_interface + 1}:{len(order)}
group mobile union all
fix far_hold far setforce 0.0 0.0 0.0
fix bottomcap all wall/harmonic zlo 0.0 100.0 2.0 2.0
fix topcap all wall/harmonic zhi 200.0 100.0 2.0 2.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
variable f2 atom fx*fx+fy*fy+fz*fz
compute fmax mobile reduce max v_fmag
compute f2mean mobile reduce ave v_f2
variable frms equal sqrt(c_f2mean)
thermo 100
thermo_style custom step pe c_fmax v_frms
thermo_modify flush yes
min_style fire
min_modify dmax 0.05
minimize 0.0 0.1 10000 50000
unfix far_hold
minimize 0.0 1.0 10000 50000
write_data geometry_precondition.data
write_dump all custom geometry_precondition.dump id type x y z fx fy fz modify sort id
"""
    input_path.write_text(input_text, encoding="utf-8")
    manifest = {
        "status": "READY_FOR_TEMPORARY_SAME_PHASE_MORSE_SOFTCORE_PRECONDITION",
        "candidate_sha256": candidate_hash,
        "input_data": data_path.name,
        "input_data_sha256": hashlib.sha256(data_path.read_bytes()).hexdigest(),
        "lammps_input": input_path.name,
        "total_atoms": int(len(order)),
        "interface_band_A": args.interface_band_A,
        "interface_atom_count": n_interface,
        "far_atom_count": int(len(order) - n_interface),
        "same_phase_elasticity": "Temporary same-phase Morse pair wells with cutoff below the next-nearest shell.",
        "same_phase_morse_r0_A": {"W": r0_w, "Ag": r0_ag},
        "same_phase_morse_alpha_A-1": morse_alpha,
        "spring_K_eV_A2": args.spring_K_eV_A2,
        "cross_soft_repulsion_eV": args.cross_soft_eV,
        "temporary_soft_repulsion_cutoffs_A": {"W_W": 2.70, "Ag_W": 2.75, "Ag_Ag": 2.48},
        "temporary_constraints": "Temporary same-phase Morse elasticity and Ag-W soft repulsion; far-field atoms are held in stage one and released in stage two. Temporary z caps keep atoms inside the material slab. All temporary terms are removed before spring-free EFS convergence.",
    }
    (out / "spring_precondition_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
