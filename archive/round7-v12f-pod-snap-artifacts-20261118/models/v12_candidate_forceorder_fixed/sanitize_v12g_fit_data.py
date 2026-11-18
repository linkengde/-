#!/workspace/venvs/dft-labels/bin/python
"""Create verified minimal local extxyz files for FitPOD and FitSNAP."""
from __future__ import annotations

import json
import numpy as np
from ase import Atoms
from ase.io import read, write
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE if (HERE / "dataset/local_train/local_dft_targets.xyz").is_file() else HERE / "v12g_si04_candidate"
SOURCE = ROOT / "dataset/local_train/local_dft_targets.xyz"
FIT_FILES = [ROOT / "dataset/train/local_dft_targets.xyz", ROOT / "SNAP/data/Local.xyz"]

source_frames = read(SOURCE, index=":")
if len(source_frames) != 14:
    raise RuntimeError(f"Expected 14 source local frames; found {len(source_frames)}")

fit_frames = []
for src in source_frames:
    forces = np.asarray(src.get_forces(apply_constraint=False), dtype=float)
    if forces.shape != (len(src), 3) or not np.isfinite(forces).all():
        raise RuntimeError("Local target frame has invalid forces")
    required = ("source_atom_id", "DFT_forces", "POD_v10_forces", "is_DFT_force_target")
    if any(key not in src.arrays for key in required):
        raise RuntimeError("Rich source frame lacks force provenance arrays")
    target = np.asarray(src.arrays["is_DFT_force_target"], dtype=bool)
    ids = np.asarray(src.arrays["source_atom_id"], dtype=int)
    dft_forces = np.asarray(src.arrays["DFT_forces"], dtype=float)
    baseline_forces = np.asarray(src.arrays["POD_v10_forces"], dtype=float)
    if target.shape != (len(src),) or int(target.sum()) < 1:
        raise RuntimeError("Expected at least one DFT force target in each local frame")
    if not np.array_equal(forces[target], dft_forces[target]):
        raise RuntimeError("Selected DFT force target was changed")
    if not np.array_equal(forces[~target], baseline_forces[~target]):
        raise RuntimeError("Local baseline force anchors were changed")
    if src.info.get("config_type") == "POD_v12g_Si_hotspot_04" and int(ids[target][0]) != 244669:
        raise RuntimeError("Si04 center ID does not match the declared training candidate")
    clean = Atoms(symbols=src.get_chemical_symbols(), positions=src.positions,
                  cell=src.cell, pbc=src.pbc)
    clean.new_array("forces", forces.copy())
    clean.info["energy"] = float(src.info.get("energy", src.get_potential_energy()))
    clean.info["config_type"] = str(src.info.get("config_type", "local_force_target"))
    fit_frames.append(clean)

for path in FIT_FILES:
    write(path, fit_frames, format="extxyz")
    check = read(path, index=":")
    if len(check) != len(source_frames):
        raise RuntimeError(f"Frame count mismatch in {path}")
    for i, (expected, actual) in enumerate(zip(fit_frames, check)):
        if not np.allclose(expected.positions, actual.positions, atol=0, rtol=0):
            raise RuntimeError(f"Position changed in {path}, frame {i}")
        if not np.allclose(expected.arrays["forces"], actual.get_forces(), atol=1e-12, rtol=0):
            raise RuntimeError(f"Force changed in {path}, frame {i}")
        if float(expected.info["energy"]) != float(actual.get_potential_energy()):
            raise RuntimeError(f"Energy changed in {path}, frame {i}")
    header = path.read_text().splitlines()[1]
    if "Properties=species:S:1:pos:R:3:forces:R:3" not in header:
        raise RuntimeError(f"Unexpected property order in {path}: {header}")

report = {
    "status": "PASS",
    "source": str(SOURCE),
    "fit_files": [str(p) for p in FIT_FILES],
    "frames": len(fit_frames),
    "property_order": "species, position[3], forces[3]",
    "note": "Rich source IDs and audit arrays are retained in dataset/local_train. Marked rows use DFT force targets; other rows use POD v10 force anchors and receive the configured local force weight.",
}
(ROOT / "local_fit_input_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
