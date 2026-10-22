#!/usr/bin/env python3
"""Verify one compact, converged V19 independent-development archive."""
import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

from entry import HERE, REPO, load, sha

if len(sys.argv) != 3:
    raise SystemExit("usage: verify_result.py LABEL ARCHIVE_DIR")
label, archive_arg = sys.argv[1], Path(sys.argv[2]).resolve()
manifest, record = load(label, launch=False)
archive = (REPO / record["result_directory"] / "calculations" / label).resolve()
assert archive_arg == archive and archive.is_dir()
input_path = (REPO / record["input"]).resolve()
output_path = archive / f"{label}_PW_PBE.extxyz"
summary_path = archive / "summary.json"
assert output_path.is_file() and summary_path.is_file() and (archive / "gpaw.log").is_file()
source = read(input_path)
result = read(output_path)
summary = json.loads(summary_path.read_text())
assert len(source) == len(result) == record["atoms"] == 52
assert source.get_chemical_symbols() == result.get_chemical_symbols()
assert np.max(np.abs(source.cell.array - result.cell.array)) < 1e-10
assert np.array_equal(source.pbc, result.pbc)
assert np.max(np.abs(source.positions - result.positions)) < 1e-10
assert np.array_equal(source.arrays["lammps_id"], result.arrays["lammps_id"])
assert np.array_equal(source.arrays["lammps_type"], result.arrays["lammps_type"])
assert np.array_equal(source.arrays["central_pair"], result.arrays["central_pair"])
assert "PW_PBE_forces" in result.arrays and result.arrays["PW_PBE_forces"].shape == (len(result), 3)
assert np.isfinite(result.arrays["PW_PBE_forces"]).all()
assert np.isfinite(float(result.info["PW_PBE_energy_eV"]))
assert summary["label"] == label and summary["input_sha256"] == record["input_sha256"]
assert summary["scf_converged"] is True and summary["scf_iterations"] > 0
assert summary["dataset_role"] == "independent_development"
assert summary["method"] == record["method"] and summary["pbc"] == record["pbc"]
assert summary["poissonsolver"] == record["poissonsolver"]
assert summary["atoms"] == 52 and summary["formula"] == record["formula"]
assert all(np.isfinite(float(summary[key])) for key in ("energy_eV_cell", "free_energy_eV_cell", "elapsed_s", "fmax_eV_A", "rms_force_eV_A"))
result_record = {
    "status": "PASS",
    "label": label,
    "archive_sha256": {
        p.name: sha(p) for p in sorted(archive.iterdir()) if p.is_file() and p.name not in {"verification.json", "SHA256SUMS.txt"}
    },
    "input_sha256": record["input_sha256"],
    "geometry_sha256": record["geometry_sha256"],
    "scf_converged": True,
    "fixed_geometry_exactly_preserved": True,
    "role": "independent_development",
    "withheld_test_paths_opened": False,
}
(archive / "verification.json").write_text(json.dumps(result_record, indent=2) + "\n")
print(json.dumps(result_record, indent=2))
