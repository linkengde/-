#!/usr/bin/env python3
"""Build new, residual-targeted Ag-C/Ag-Si/Ag-Ti PW-PBE single-point inputs."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
SOURCE_ROOT = PROJECT / "pbe_interface_v14_parallel_acquisition"
REPORT = PROJECT.parent / "coordination/reports/window-b/v13_force_localization.json"
V14 = PROJECT / "mace_periodic_v14_interface_energy/data"
MAIN_HOLDOUTS = PROJECT / "pbe_interface_v14_main_holdouts/calculations"
B_HOLDOUTS = SOURCE_ROOT / "calculations"
INPUTS = ROOT / "inputs"

if any(INPUTS.iterdir()):
    raise SystemExit("v15 inputs already exist; inspect their hashes instead of regenerating them.")

report = json.loads(REPORT.read_text())
sources = {row["label"]: row for row in report["sources"]}
specs = [
    ("AgC_residual_shell_v15_01", "AgC_registry_strain_acq_v14_01", "AgC", 0.020, 15114),
    ("AgSi_residual_shell_v15_01", "AgSi_registry_strain_acq_v14_01", "AgSi", 0.035, 15102),
    ("AgTi_residual_shell_v15_01", "AgTi_registry_probe_v13_01", "AgTi", 0.030, 15103),
]

protected = [read(p, index=":") for p in (V14 / "train.extxyz", V14 / "valid.extxyz", V14 / "test.extxyz")]
for folder in (MAIN_HOLDOUTS, B_HOLDOUTS):
    for output in sorted(folder.glob("*/*_PW_PBE.extxyz")):
        protected.extend(read(output, index=":"))

records = []
for label, source_label, interface, sigma, seed in specs:
    source = SOURCE_ROOT / "inputs" / f"{source_label}.extxyz"
    atoms = read(source)
    row = sources[source_label]
    by_id = {int(atom_id): i for i, atom_id in enumerate(atoms.arrays["lammps_id"])}
    pair_ids = [int(value) for value in row["marked_pair"]["ids"]]
    top = row["top_10_atom_errors"][:3]
    focus_ids = list(dict.fromkeys(pair_ids + [int(atom["id"]) for atom in top]))
    missing = sorted(set(focus_ids) - set(by_id))
    if missing:
        raise SystemExit(f"Residual target IDs absent from {source_label}: {missing}")
    focus = np.asarray([by_id[value] for value in focus_ids], dtype=int)
    distances = atoms.get_all_distances(mic=True)
    local = np.flatnonzero(np.min(distances[:, focus], axis=1) <= 4.5)
    rng = np.random.default_rng(seed)
    atoms.positions[local] += rng.normal(0.0, sigma, (len(local), 3))

    # Move the largest-residual atom a small, documented distance along its v13 force residual.
    top_atom = top[0]
    idx = by_id[int(top_atom["id"])]
    residual = np.asarray(top_atom["force_error_vector_eV_A"], dtype=float)
    atoms.positions[idx] += 0.05 * residual / np.linalg.norm(residual)

    for key in ("REF_energy", "PW_PBE_energy_eV", "energy", "free_energy", "source_method"):
        atoms.info.pop(key, None)
    for key in ("REF_forces", "PW_PBE_forces", "forces", "magmoms"):
        atoms.arrays.pop(key, None)
    atoms.calc = None
    atoms.info.update(
        config_type=label,
        dataset_role="v15_acquisition_candidate",
        geometry_family=f"{interface}_residual_targeted_local_perturbation",
        source_label=source_label,
        targeted_lammps_ids=",".join(map(str, focus_ids)),
        correlation_note="Targeted local perturbation of an existing small-cell motif; not an independent morphology or thermal snapshot.",
    )

    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    minimum = float(d.min())
    symbol = np.asarray(atoms.get_chemical_symbols())
    ag = np.flatnonzero(symbol == "Ag")
    other = np.flatnonzero(symbol != "Ag")
    min_ag_other = float(d[np.ix_(ag, other)].min())
    if minimum < 1.75 or min_ag_other < 1.75:
        raise SystemExit(f"Unexpected close contact for {label}: min={minimum:.6f} Å, Ag-other={min_ag_other:.6f} Å")
    for prior in protected:
        if (len(prior) == len(atoms) and prior.get_chemical_symbols() == atoms.get_chemical_symbols()
                and np.allclose(prior.cell.array, atoms.cell.array, atol=1e-12, rtol=0)
                and np.allclose(prior.positions, atoms.positions, atol=1e-10, rtol=0)):
            raise SystemExit("Generated geometry duplicates a protected v14/train/test geometry: " + label)
    out = INPUTS / f"{label}.extxyz"
    write(out, atoms, format="extxyz")
    pair = [by_id[value] for value in pair_ids]
    records.append({
        "label": label,
        "role": "v15_acquisition_candidate",
        "input": f"inputs/{out.name}",
        "input_sha256": sha256(out.read_bytes()).hexdigest(),
        "source": f"pbe_interface_v14_parallel_acquisition/inputs/{source.name}",
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "source_label": source_label,
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "targeted_lammps_ids": focus_ids,
        "perturbation": {"local_gaussian_sigma_A": sigma, "seed": seed,
                         "top_residual_atom_id": int(top_atom["id"]),
                         "top_residual_displacement_A": 0.05},
        "central_pair": {
            "elements": f"Ag-{interface[2:]}",
            "Ag_id": pair_ids[0] if symbol[pair[0]] == "Ag" else pair_ids[1],
            "contact_id": pair_ids[1] if symbol[pair[0]] == "Ag" else pair_ids[0],
            "input_distance_A": float(atoms.get_distance(pair[0], pair[1], mic=True)),
        },
        "minimum_pair_distance_A": minimum,
        "minimum_Ag_nonAg_distance_A": min_ag_other,
        "purpose": "Add local-force coverage for v13 residual hotspots if v14 screening identifies a need; excluded from v14 training and all frozen/reserved tests.",
    })

manifest = {
    "owner_task": "window-a",
    "purpose": "Residual-targeted DFT acquisition candidates prepared in parallel with the B v14 model run.",
    "method": "GPAW 26.7.0 PW-PBE, 500 eV, Gamma, four MPI ranks, fixed-geometry single points.",
    "v14_isolation": "These labels are not part of the frozen v14 train/valid/test or its three reserved holdouts. Do not use them to change the v14 model before scoring it.",
    "source_report_sha256": sha256(REPORT.read_bytes()).hexdigest(),
    "records": records,
}
(ROOT / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"prepared": len(records), "records": [
    {k: r[k] for k in ("label", "role", "formula", "central_pair", "minimum_pair_distance_A", "minimum_Ag_nonAg_distance_A")}
    for r in records
]}, indent=2))
