#!/usr/bin/env python3
"""Make new interface-registry holdouts by lateral Ag translation and small local displacements."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
INPUTS = ROOT / "inputs"
if INPUTS.exists() and any(INPUTS.iterdir()):
    raise SystemExit(f"Refusing to overwrite existing holdout inputs in {INPUTS}")
INPUTS.mkdir(parents=True, exist_ok=True)

SOURCES = [
    {
        "label": "AgC_lateral_registry_holdout_v13",
        "source": PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz",
        "contact": "C",
        "seed": 13001,
    },
    {
        "label": "AgSi_lateral_registry_holdout_v13",
        "source": PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz",
        "contact": "Si",
        "seed": 13002,
    },
]

records = []
for item in SOURCES:
    source_path = item["source"]
    source_hash = sha256(source_path.read_bytes()).hexdigest()
    atoms = read(source_path)
    if not np.all(atoms.pbc):
        raise SystemExit(f"Source must be fully periodic: {source_path}")
    if "central_pair" not in atoms.arrays or int(np.sum(atoms.arrays["central_pair"])) != 2:
        raise SystemExit(f"Expected exactly one tagged central pair in {source_path}")
    symbols = np.asarray(atoms.get_chemical_symbols())
    central = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
    ag_index = int(central[symbols[central] == "Ag"][0])
    contact_index = int(central[symbols[central] == item["contact"]][0])
    old_vector = atoms.positions[contact_index] - atoms.positions[ag_index]
    old_distance = float(np.linalg.norm(old_vector))

    # Translate the full Ag sublattice laterally, changing registry while preserving Ag geometry.
    lateral = np.cross(old_vector, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(lateral) < 1e-8:
        lateral = np.cross(old_vector, np.array([0.0, 1.0, 0.0]))
    lateral /= np.linalg.norm(lateral)
    lateral *= 0.18
    ag_indices = np.flatnonzero(symbols == "Ag")
    atoms.positions[ag_indices] += lateral

    # Add a small, deterministic thermal-like displacement to the nearby carbide/silicide environment.
    distances = atoms.get_all_distances(mic=True)
    interface_atoms = np.flatnonzero(np.min(distances[:, central], axis=1) < 4.2)
    movable = interface_atoms[~np.isin(interface_atoms, ag_indices)]
    rng = np.random.default_rng(item["seed"])
    atoms.positions[movable] += rng.normal(0.0, 0.015, size=(len(movable), 3))

    # Strip every reference label so these remain blind DFT checks.
    for key in ("REF_energy", "PW_PBE_energy_eV", "energy"):
        atoms.info.pop(key, None)
    for key in ("REF_forces", "PW_PBE_forces", "forces"):
        atoms.arrays.pop(key, None)
    atoms.calc = None
    atoms.info["config_type"] = item["label"]
    atoms.info["validation_role"] = "independent_v13_geometry_holdout"
    atoms.info["source_geometry_sha256"] = source_hash
    atoms.info["registry_translation_A"] = " ".join(f"{x:.10f}" for x in lateral)
    atoms.info["local_noise_sigma_A"] = 0.015
    output = INPUTS / f"{item['label']}.extxyz"
    write(output, atoms, format="extxyz")

    new_pair_vector = atoms.positions[contact_index] - atoms.positions[ag_index]
    new_distance = float(np.linalg.norm(new_pair_vector))
    records.append({
        "label": item["label"],
        "input": f"inputs/{output.name}",
        "input_sha256": sha256(output.read_bytes()).hexdigest(),
        "source": str(source_path.relative_to(PROJECT)),
        "source_sha256": source_hash,
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "central_pair": {
            "elements": f"Ag-{item['contact']}",
            "Ag_id": int(atoms.arrays["lammps_id"][ag_index]),
            "contact_id": int(atoms.arrays["lammps_id"][contact_index]),
            "source_distance_A": old_distance,
            "input_distance_A": new_distance,
            "registry_translation_A": lateral.tolist(),
        },
        "interface_atoms_randomly_displaced": int(len(movable)),
        "local_noise_sigma_A": 0.015,
        "random_seed": item["seed"],
        "labels_removed_from_input": True,
    })

manifest = {
    "purpose": "Independent v13 checks using a lateral Ag-registry shift and small interface-local displacements, rather than central-pair distance-only scans.",
    "method": "Translate all Ag atoms by 0.18 A perpendicular to the tagged contact vector; apply deterministic 0.015 A Gaussian displacements to non-Ag atoms within 4.2 A of the tagged pair; fixed-cell single-point PBE labels remain outside v13 training. GPAW uses compiled ScaLAPACK through parallel.sl_auto=true.",
    "records": records,
}
(ROOT / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
