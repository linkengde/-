#!/usr/bin/env python3
"""Prepare frozen, labeled ACE/MTP candidate datasets after k-point review."""
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write


HERE = Path(__file__).resolve().parent
FIT_ROOT = HERE.parent
LABEL = HERE / "k1x6x3" / "ti_gap_1x6x3.extxyz"
FRAME = "Ti_gap_2p3_rank_recovery_v17_01"
TYPE = {"Ag": 0, "C": 1, "Si": 2, "Ti": 3}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def as_cfg(atoms_list, path):
    with Path(path).open("w") as out:
        for source_atoms in atoms_list:
            atoms = source_atoms.copy()
            if atoms.pbc.any() and not atoms.pbc.all():
                old_distances = atoms.get_all_distances(mic=True)
                cell = np.asarray(atoms.cell).copy()
                for axis in range(3):
                    if not atoms.pbc[axis]:
                        height = float(np.linalg.norm(cell[axis]))
                        span = float(np.ptp(atoms.get_scaled_positions(wrap=False)[:, axis]) * height)
                        if height - span < 8.0:
                            cell[axis] *= (span + 8.0) / height
                atoms.set_cell(cell, scale_atoms=False)
                atoms.pbc = True
                new_distances = atoms.get_all_distances(mic=True)
                mask = (old_distances < 6.0) | (new_distances < 6.0)
                if not np.allclose(old_distances[mask], new_distances[mask], atol=1e-8):
                    raise ValueError(f"Vacuum padding changed local neighbors: {atoms.info.get('config_type')}")
            elif not atoms.pbc.all():
                raise ValueError(f"Unperiodic frame cannot be serialized to MLIP CFG: {atoms.info.get('config_type')}")
            energy = atoms.info.get("REF_energy")
            forces = atoms.arrays.get("REF_forces")
            if energy is None or forces is None or np.asarray(forces).shape != (len(atoms), 3):
                raise ValueError(f"Missing energy/forces: {atoms.info.get('config_type')}")
            out.write("BEGIN_CFG\n Size\n    %d\n Supercell\n" % len(atoms))
            for row in np.asarray(atoms.cell, dtype=float):
                out.write("    %.16f %.16f %.16f\n" % tuple(row))
            out.write(" AtomData:  id type cartes_x cartes_y cartes_z fx fy fz\n")
            for i, (symbol, position, force) in enumerate(zip(atoms.get_chemical_symbols(), atoms.positions, forces), 1):
                if symbol not in TYPE:
                    raise ValueError(f"Unexpected species {symbol}")
                values = [*position, *force]
                out.write("    %d %d %s\n" % (i, TYPE[symbol], " ".join("%.16f" % float(x) for x in values)))
            out.write(" Energy\n    %.16f\n" % float(energy))
            out.write(" Feature   config_type %s\nEND_CFG\n\n" % atoms.info.get("config_type", "unknown"))


def main():
    if not LABEL.is_file():
        raise FileNotFoundError(LABEL)
    dft = read(LABEL)
    if dft.info.get("config_type") != FRAME:
        raise ValueError(f"Wrong relabel frame: {dft.info.get('config_type')}")
    train = read(FIT_ROOT / "train.extxyz", ":")
    dev = read(FIT_ROOT / "dev.extxyz", ":")
    targets = [a for a in train if a.info.get("config_type") == FRAME]
    if len(targets) != 1:
        raise ValueError(f"Expected one training frame {FRAME}, found {len(targets)}")
    old = targets[0]
    if (old.get_chemical_symbols() != dft.get_chemical_symbols() or
            not np.allclose(old.positions, dft.positions, atol=1e-10, rtol=0) or
            not np.allclose(old.cell, dft.cell, atol=1e-10, rtol=0)):
        raise ValueError("DFT label geometry does not match the frozen training frame")
    old.info["REF_energy"] = float(dft.info["REF_energy"])
    old.arrays["REF_forces"] = np.asarray(dft.arrays["REF_forces"], dtype=float).copy()
    old.info["source_method"] = dft.info["source_method"]
    old.info["kpoint_mesh_reviewed"] = "1x6x3; see result_resume.json and comparison.json"
    outdir = HERE / "fit_candidate"
    outdir.mkdir(exist_ok=True)
    train_xyz, dev_xyz = outdir / "train_kconv.extxyz", outdir / "dev.extxyz"
    write(train_xyz, train, format="extxyz")
    write(dev_xyz, dev, format="extxyz")
    as_cfg(train, outdir / "train.cfg")
    as_cfg(dev, outdir / "dev.cfg")
    summary = {
        "status": "candidate datasets prepared; exploratory only",
        "replaced_frame": FRAME,
        "mlip_type_map": TYPE,
        "training_frames": len(train),
        "development_frames": len(dev),
        "source_train_sha256": sha256(FIT_ROOT / "train.extxyz"),
        "kpoint_label_sha256": sha256(LABEL),
        "candidate_train_sha256": sha256(train_xyz),
        "candidate_dev_sha256": sha256(dev_xyz),
        "candidate_train_cfg_sha256": sha256(outdir / "train.cfg"),
        "candidate_dev_cfg_sha256": sha256(outdir / "dev.cfg"),
    }
    (outdir / "candidate_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
