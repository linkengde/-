#!/usr/bin/env python3
"""Package DFT inputs for the three frozen M77 geometries; never launches GPAW."""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
SOURCE_DIR = REPO / "research/mace-v12-transfer/coordination/reports/window-c/minimal77_d4j_model_readiness_20261014"
SOURCE_MANIFEST = SOURCE_DIR / "MODEL_MANIFEST.json"
INPUT_DIR = ROOT / "DFT_INPUTS"
SOURCE_217 = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_fit_dft_inputs_20261014/inputs/interface_scan/AGTSC_SI_GAP_2p4A.extxyz"

PAW_SHA = {
    "Ag.PBE.gz": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
    "C.PBE.gz": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
    "Si.PBE.gz": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
    "Ti.PBE.gz": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
}
METHOD = {
    "code": "GPAW 26.7.0 / ASE 3.29.0 / gpaw-data 1.2.1",
    "xc": "PBE",
    "paw": True,
    "mode": "PW",
    "cutoff_eV": 500,
    "occupations": {"type": "FermiDirac", "width_eV": 0.10},
    "convergence": {"energy_eV_per_electron": 1e-5, "density": 1e-5, "eigenstates_eV": 1e-5},
    "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
    "maxiter": 160,
    "kpoints": [4, 4, 1],
    "kpoint_convention": "GPAW default Monkhorst-Pack integer mesh; no explicit gamma override or extra offset",
    "kpoint_status": "PROPOSED_FIRST_MESH_NOT_CONVERGENCE_VALIDATED; prior M77 design recommends a 5x5x1 denser central-gap comparison",
    "geometry": "fixed ions and fixed cell; no relaxation",
    "resources": {"mpi_ranks": 4, "omp_threads": 1},
    "outputs": ["native_energy_eV", "force_consistent_free_energy_eV", "forces_eV_per_A", "stress_eV_per_A3_diagnostic_only"],
    "paw_setup_sha256": PAW_SHA,
    "stress_note": "Slab 3D stress includes vacuum dilution; diagnostic only, not a fitting target unless a slab stress convention is approved.",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main() -> None:
    if INPUT_DIR.exists():
        raise SystemExit(f"Refusing to overwrite existing input directory: {INPUT_DIR}")
    source_manifest = json.loads(SOURCE_MANIFEST.read_text())
    if source_manifest["natoms"] != 77 or len(source_manifest["structures"]) != 3:
        raise SystemExit("Unexpected D5 model manifest identity")
    if sha(SOURCE_217) != "f776237b5b86ce9292d2c667adc6055f74389743f967a64481d5f5e1a981ced6":
        raise SystemExit("Archived 217-atom comparison input SHA mismatch")

    struct_dir = INPUT_DIR / "structures"
    config_dir = INPUT_DIR / "configs"
    struct_dir.mkdir(parents=True)
    config_dir.mkdir(parents=True)
    records = []
    for s in source_manifest["structures"]:
        src = SOURCE_DIR / s["path"]
        dst = struct_dir / Path(s["path"]).name
        if sha(src) != s["sha256"]:
            raise SystemExit(f"D5 source structure SHA mismatch: {src}")
        shutil.copyfile(src, dst)
        atoms = read(dst)
        symbols = atoms.get_chemical_symbols()
        composition = dict(sorted(Counter(symbols).items()))
        if len(atoms) != 77 or composition != {"Ag": 28, "C": 14, "Si": 14, "Ti": 21}:
            raise SystemExit(f"Unexpected structure identity: {dst.name}")
        if not np.array_equal(atoms.arrays["atom_id"], np.arange(1, 78)):
            raise SystemExit(f"Atom ID/order array mismatch: {dst.name}")
        if atoms.pbc.tolist() != [True, True, True]:
            raise SystemExit(f"PBC mismatch: {dst.name}")
        task_id = "C-D6-" + s["structure_id"].replace("M77_SHARED_STRAIN_", "").replace("_", "-")
        config = {
            "schema": "ag_ti3sic2_fixed_geometry_dft_input_v1",
            "task_id": task_id,
            "structure_id": s["structure_id"],
            "structure_path": f"DFT_INPUTS/structures/{dst.name}",
            "structure_sha256": sha(dst),
            "natoms": 77,
            "composition": composition,
            "ordered_symbols_sha256": hashlib.sha256(" ".join(symbols).encode()).hexdigest(),
            "atom_id_policy": "extxyz atom_id=1..77; substrate 1..49 in bottom-to-top layer/site order, Ag 50..77 in bottom-to-top ABC-A layer/site order",
            "cell_A": atoms.cell.array.tolist(),
            "pbc": [True, True, True],
            "interface_plane_gap_A": float(s["gap_A"]),
            "biaxial_strain_percent": {"Ti3SiC2": -3.0044589063391136, "Ag": 3.0975227030672903},
            "twist_deg": 0.0,
            "parent_family": "M77_AREA7_SHARED_STRAIN_DERIVED_FROM_COD9009647_AND_ARCHIVED_217_PARENT",
            "data_role": "SCREENING_CANDIDATE_ROLE_UNASSIGNED",
            "dft_label_status": "NOT_GENERATED",
            "execution_authorized": False,
            "launch_enabled": False,
            "geometry_policy": "fixed positions and fixed cell; no optimization, no atom motion",
            "method": METHOD,
        }
        config_path = config_dir / f"{s['structure_id']}.json"
        json_write(config_path, config)
        records.append({
            "task_id": task_id,
            "structure_id": s["structure_id"],
            "structure_path": f"DFT_INPUTS/structures/{dst.name}",
            "structure_sha256": sha(dst),
            "input_config_path": f"DFT_INPUTS/configs/{config_path.name}",
            "input_config_sha256": sha(config_path),
            "natoms": len(atoms),
            "composition": composition,
            "cell_a1_A": atoms.cell.array[0].tolist(),
            "cell_a2_A": atoms.cell.array[1].tolist(),
            "cell_a3_A": atoms.cell.array[2].tolist(),
            "cell_volume_A3": float(atoms.cell.volume),
            "interface_gap_A": float(s["gap_A"]),
            "TSC_strain_percent": -3.0044589063391136,
            "Ag_strain_percent": 3.0975227030672903,
            "twist_deg": 0.0,
            "PBC": "T T T",
            "execution_authorized": False,
        })

    with (ROOT / "STRUCTURE_LIST.csv").open("w", newline="") as f:
        fields = ["task_id", "structure_id", "structure_path", "structure_sha256", "input_config_path",
                  "input_config_sha256", "natoms", "composition", "cell_a1_A", "cell_a2_A", "cell_a3_A",
                  "cell_volume_A3", "interface_gap_A", "TSC_strain_percent", "Ag_strain_percent",
                  "twist_deg", "PBC", "execution_authorized"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in records:
            row["composition"] = json.dumps(row["composition"], sort_keys=True)
            for key in ("cell_a1_A", "cell_a2_A", "cell_a3_A"):
                row[key] = json.dumps(row[key])
            w.writerow(row)

    manifest = {
        "package_status": "INPUT_FILES_PREPARED_NOT_RELEASED",
        "DFT_STARTED": False,
        "execution_authorized": False,
        "D5_source_manifest_sha256": sha(SOURCE_MANIFEST),
        "builder_sha256": sha(Path(__file__)),
        "runner_path": "run_fixed_singlepoint.py",
        "runner_sha256": sha(ROOT / "run_fixed_singlepoint.py"),
        "structures": records,
        "method": METHOD,
        "notes": [
            "4x4x1 is the first proposed comparable-density mesh from the archived M77 resource plan, not a demonstrated converged mesh.",
            "A 5x5x1 central-gap comparison is proposed before labels are released; no k-point test was run.",
            "Peak memory and wall time are unmeasured for these exact geometries.",
        ],
    }
    json_write(ROOT / "INPUT_MANIFEST.json", manifest)


if __name__ == "__main__":
    main()
