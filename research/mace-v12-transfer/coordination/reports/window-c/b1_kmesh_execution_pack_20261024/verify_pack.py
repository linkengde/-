#!/usr/bin/env python3
"""Read-only validation of the B1 mesh pack; deliberately contains no DFT runner."""
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

from ase.io import read

ROOT = Path(__file__).resolve().parents[6]
PACK = Path(__file__).resolve().parent
BASE = ROOT / "research/mace-v12-transfer"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest_path = PACK / "execution_manifest.json"
manifest = json.loads(manifest_path.read_text())
assert manifest["schema"] == "B1_PHASE_KMESH_EXECUTION_PACK_V1"
assert manifest["global_launch_enabled"] is False
assert manifest["global_dft_authorized"] is False
assert manifest["training_enabled"] is False
assert manifest["dataset_integration_enabled"] is False
assert manifest["record_count"] == 12 == len(manifest["records"])

plan_path = BASE / manifest["source_mesh_plan"]
ledger_path = BASE / manifest["source_role_ledger"]
assert sha(plan_path) == manifest["source_mesh_plan_sha256"]
assert sha(ledger_path) == manifest["source_role_ledger_sha256"]
roles = {r["config_type"]: r for r in csv.DictReader(ledger_path.open())}

pairs = {}
for record in manifest["records"]:
    assert record["launch_enabled"] is False
    assert record["dft_authorized"] is False
    assert record["dataset_integration_enabled"] is False
    assert record["energy_target"] == "free_energy_eV"
    assert record["native_energy_retained_separately"] is True
    assert record["candidate_role"] == "TRAINING_CANDIDATE_ONLY"
    assert record["pair_member"] in ("coarse", "fine")
    assert len(record["kpts"]) == 3 and all(type(k) is int and k >= 1 for k in record["kpts"])

    source = (ROOT / record["input"]).resolve()
    assert source.is_relative_to(ROOT)
    assert sha(source) == record["input_sha256"]
    role = roles[record["config_type"]]
    assert role["file_sha256"] == record["input_sha256"]
    assert role["candidate_role"] == "TRAINING_CANDIDATE_ONLY"
    assert role["dft_label_present"] == "False"
    assert role["dft_launch_enabled"] == "False"
    assert role["dataset_integrated"] == "False"

    atoms = read(source)
    assert len(atoms) == record["natoms"]
    assert atoms.pbc.tolist() == [True, True, True] == record["pbc"]
    composition = Counter(atoms.get_chemical_symbols())
    if record["phase"] in ("solid_Ag", "liquid_Ag"):
        assert set(composition) == {"Ag"}
    else:
        assert set(composition) == {"Ti", "Si", "C"}
        assert composition["Ti"] == 3 * composition["Si"]
        assert composition["C"] == 2 * composition["Si"]
    assert not any(k in atoms.info for k in ("REF_energy", "PW_PBE_energy_eV", "free_energy_eV", "energy"))
    assert not any(k in atoms.arrays for k in ("REF_forces", "PW_PBE_forces", "forces"))

    method = record["method"]
    assert method["kpts"] == record["kpts"] and method["pbc"] == record["pbc"]
    assert method["poissonsolver"] == record["poissonsolver"] == {}
    assert method["GPAW_version"] == "26.7.0"
    assert method["ASE_version"] == "3.29.0"
    assert method["PAW_data_version"] == "1.2.1"
    assert method["xc"] == "PBE" and method["cutoff_eV"] == 500
    assert method["smearing_eV"] == 0.1
    pairs.setdefault(record["pair_id"], []).append(record)

assert len(pairs) == 6
for pair_id, members in pairs.items():
    assert len(members) == 2, pair_id
    assert {m["pair_member"] for m in members} == {"coarse", "fine"}
    assert len({m["input_sha256"] for m in members}) == 1, pair_id
    for key in ("input", "input_sha256", "config_type", "method"):
        if key == "method":
            normalized = []
            for member in members:
                value = dict(member[key])
                value.pop("kpts")
                normalized.append(value)
            assert normalized[0] == normalized[1], pair_id
        else:
            assert members[0][key] == members[1][key], pair_id

# Verify the already-installed scientific stack and exact PAW payloads. No install or calculation.
import ase
import gpaw
import gpaw_data
assert ase.__version__ == "3.29.0"
assert gpaw.__version__ == "26.7.0"
setup_dir = Path(gpaw_data.datapath())
for name, expected in manifest["provisional_method"]["PAW_setup_sha256"].items():
    assert sha(setup_dir / name) == expected, name
gpaw_info = subprocess.check_output(
    [str(Path(sys.executable).with_name("gpaw")), "info"], text=True
)
assert re.search(r"MPI enabled\s+yes", gpaw_info)
assert re.search(r"scalapack\s+yes", gpaw_info)

print(
    json.dumps(
        {
            "status": "PASS_PACK_VALIDATED_LAUNCH_DISABLED",
            "records": len(manifest["records"]),
            "same_geometry_mesh_pairs": len(pairs),
            "candidate_inputs_label_free": True,
            "software": {"GPAW": gpaw.__version__, "ASE": ase.__version__, "MPI": True, "ScaLAPACK": True},
            "paw_setup_hashes": "PASS",
            "dft_launched": False,
            "training_or_integration": False,
        },
        indent=2,
    )
)
