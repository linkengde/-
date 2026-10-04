#!/usr/bin/env python3
"""Create reference-free pair-distance variants for independent PW-PBE checks."""
from pathlib import Path
import hashlib
import json

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
DEST = ROOT / "inputs/validation_distance_scans"
if DEST.exists() and any(DEST.iterdir()):
    raise SystemExit(f"Refusing to overwrite non-empty input directory: {DEST}")
DEST.mkdir(parents=True, exist_ok=True)

SPECS = [
    {
        "label": "AgC_validation_d2p30",
        "source": PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz",
        "contact": "C",
        "target_distance_A": 2.30,
    },
    {
        "label": "AgSi_validation_d2p60",
        "source": PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz",
        "contact": "Si",
        "target_distance_A": 2.60,
    },
    {
        "label": "AgTi_validation_d2p60",
        "source": ROOT.parent / "mace_periodic_v12_interface_energy/data/test.extxyz",
        "source_index": 2,
        "contact": "Ti",
        "target_distance_A": 2.60,
    },
]

records = []
for spec in SPECS:
    source = Path(spec["source"])
    source_atoms = read(source, index=spec.get("source_index", 0))
    atoms = source_atoms.copy()
    atoms.calc = None
    central = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    indices = np.flatnonzero(central)
    if len(indices) != 2:
        raise SystemExit(f"Expected two central-pair atoms in {source}")
    symbols = atoms.get_chemical_symbols()
    ag = [int(i) for i in indices if symbols[i] == "Ag"]
    contact = [int(i) for i in indices if symbols[i] == spec["contact"]]
    if len(ag) != 1 or len(contact) != 1:
        raise SystemExit(f"Unexpected central-pair elements in {source}: {[symbols[i] for i in indices]}")
    ia, ib = ag[0], contact[0]
    vector = atoms.positions[ib] - atoms.positions[ia]
    distance = float(np.linalg.norm(vector))
    unit = vector / distance
    displacement = float(spec["target_distance_A"]) - distance
    atoms.positions[ia] -= 0.5 * displacement * unit
    atoms.positions[ib] += 0.5 * displacement * unit
    actual_distance = float(np.linalg.norm(atoms.positions[ib] - atoms.positions[ia]))
    if not np.isclose(actual_distance, spec["target_distance_A"], rtol=0, atol=1e-10):
        raise SystemExit(f"Failed to create requested central-pair distance for {spec['label']}")

    # Keep only geometry identity arrays. Reference labels from the source must
    # never leak into these validation inputs.
    for key in ("REF_forces", "PW_PBE_forces", "forces", "stress", "virials"):
        atoms.arrays.pop(key, None)
    atoms.info = {
        "config_type": spec["label"],
        "source_config_type": str(source_atoms.info.get("config_type", "")),
        "source_structure_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "validation_target_distance_A": float(spec["target_distance_A"]),
        "validation_source_distance_A": distance,
        "central_pair_elements": f"Ag-{spec['contact']}",
    }
    atoms.pbc = True
    input_path = DEST / f"{spec['label']}.extxyz"
    write(input_path, atoms, format="extxyz")
    reread = read(input_path)
    saved_distance = float(np.linalg.norm(reread.positions[ib] - reread.positions[ia]))
    if not np.isclose(saved_distance, spec["target_distance_A"], rtol=0, atol=1e-8):
        raise SystemExit(f"Serialized input missed requested pair distance for {spec['label']}")
    if "energy" in reread.info or "REF_energy" in reread.info or "PW_PBE_energy_eV" in reread.info:
        raise SystemExit(f"Reference energy leaked into validation input: {input_path}")
    if any(key in reread.arrays for key in ("forces", "REF_forces", "PW_PBE_forces")):
        raise SystemExit(f"Reference forces leaked into validation input: {input_path}")
    records.append({
        "label": spec["label"],
        "input": str(input_path.relative_to(ROOT)),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source": str(source.relative_to(PROJECT)),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "central_pair": {"elements": f"Ag-{spec['contact']}", "Ag_id": int(atoms.arrays["lammps_id"][ia]), "contact_id": int(atoms.arrays["lammps_id"][ib]), "source_distance_A": distance, "target_distance_A": float(spec["target_distance_A"]), "input_distance_A": saved_distance},
        "labels_removed_from_source": True,
    })

manifest = {"method": "Fixed-geometry central-pair displacement; all other atoms and cell preserved.", "records": records}
(ROOT / "validation_distance_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
