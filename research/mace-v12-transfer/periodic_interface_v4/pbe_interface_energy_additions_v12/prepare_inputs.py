#!/usr/bin/env python3
"""Extract frozen AgSi/AgC interface scan frames for same-method PW-PBE labels."""
from pathlib import Path
import hashlib
import json
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "mace_periodic_v11_energyfocus/data/train.extxyz"
OUT = ROOT / "inputs"
OUT.mkdir(parents=True, exist_ok=True)
if any(OUT.iterdir()):
    raise SystemExit(f"Refusing to overwrite existing input files in {OUT}")

wanted = {
    "AgSi_d2p2": "AgSi_2p2A_LCAO",
    "AgSi_d2p8": "AgSi_2p8A_LCAO",
    "AgC_d2p0": "AgC_secondary_2.0A_LCAO",
    "AgC_d2p4": "AgC_secondary_2.4A_LCAO",
}
frames = read(SOURCE, index=":")
records = []
for out_name, config_type in wanted.items():
    matches = [a.copy() for a in frames if a.info.get("config_type") == config_type]
    if len(matches) != 1:
        raise SystemExit(f"Expected one {config_type}, found {len(matches)}")
    atoms = matches[0]
    atoms.pbc = True
    path = OUT / f"{out_name}.extxyz"
    write(path, atoms, format="extxyz")
    central = atoms.arrays["central_pair"].astype(bool)
    ids = atoms.arrays["lammps_id"].astype(int)
    symbols = atoms.get_chemical_symbols()
    pair = {symbols[i]: int(i) for i in central.nonzero()[0]}
    distance = float(atoms.get_distance(pair["Ag"], pair["Si"] if "Si" in pair else pair["C"], mic=True))
    records.append({
        "label": out_name,
        "source_config_type": config_type,
        "source_file": str(SOURCE),
        "input_file": str(path),
        "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "contact_distance_A": distance,
        "contact_ids": {e: int(ids[i]) for e, i in pair.items()},
        "cell_A": atoms.cell.lengths().tolist(),
        "pbc": atoms.pbc.tolist(),
        "prior_energy_label": "LCAO/force-scan labels excluded; new PW-PBE single-point energy and forces will replace these frames in v12 training",
    })

manifest = {
    "purpose": "Add same-method GPAW PW-PBE energy and force labels for AgSi and AgC contacts.",
    "method": {"xc": "PBE", "basis": "plane wave", "cutoff_eV": 500, "kpts": [1, 1, 1], "smearing_eV": 0.1, "convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5}, "mixer": {"beta": 0.05, "nmaxold": 8, "weight": 100}},
    "holdouts_preserved": ["AgSi external d=2.4541 A (51 atoms)", "AgC external d=2.2026 A (48 atoms)", "AgTi d=2.70 A frozen test"],
    "records": records,
}
(ROOT / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
