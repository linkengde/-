#!/usr/bin/env python3
"""Evaluate the registered baseline ACE on v11 Si hotspot development center 01."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from pyace import ACEBBasisSet, PyACECalculator

ROOT = Path(__file__).resolve().parent
MODEL = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/ACE/AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight1.ace.yaml")
BASE_METRICS = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/ACE/metrics_weight1.json")
HOTSPOT = ROOT / "v11_hotspot_reference/archive/round7-v11-si-hotspot-labels-20261118/round7_pod_snap/models/v11_candidate/active_si_hotspots/labels/si_hotspot_01_pw350_gamma.extxyz"
OUT = ROOT / "baseline_ace_hotspot01_metrics.json"
GATE = 0.5


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    base = json.loads(BASE_METRICS.read_text())
    atoms = read(HOTSPOT)
    center_id = int(atoms.info["center_atom_id"])
    ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
    center = np.flatnonzero(ids == center_id)
    if len(center) != 1:
        raise SystemExit("registered development center is not unique")
    center_idx = int(center[0])
    reference = np.asarray(atoms.arrays["REF_forces"][center_idx], dtype=float)
    atoms.calc = PyACECalculator(ACEBBasisSet(str(MODEL)))
    predicted = np.asarray(atoms.get_forces()[center_idx], dtype=float)
    error = float(np.linalg.norm(predicted - reference))
    result = {
        "status": "baseline direct-center evaluation",
        "baseline_model": str(MODEL),
        "baseline_model_sha256": sha(MODEL),
        "baseline_fit_metrics_sha256": sha(BASE_METRICS),
        "baseline_train": {key: base["train"][key] for key in ["energy_rmse_meV_atom", "force_component_rmse_eV_A"]},
        "baseline_periodic_development": {key: base["development"][key] for key in ["energy_rmse_meV_atom", "force_component_rmse_eV_A"]},
        "hotspot_label": str(HOTSPOT),
        "hotspot_label_sha256": sha(HOTSPOT),
        "hotspot_center_atom_id": center_id,
        "hotspot_center_reference_force_eV_A": reference.tolist(),
        "hotspot_center_ACE_force_eV_A": predicted.tolist(),
        "hotspot_center_vector_error_eV_A": error,
        "hotspot_center_gate_eV_A": GATE,
        "train_dev_split": "v11 indices 00/02/03 training, 01 development; not blind test",
        "energy_and_force_gate_pass": bool(base["development"]["energy_rmse_meV_atom"] <= 10.0 and base["development"]["force_component_rmse_eV_A"] <= 0.1),
        "hotspot_center_gate_pass": error <= GATE,
        "all_gates_pass": bool(base["development"]["energy_rmse_meV_atom"] <= 10.0 and base["development"]["force_component_rmse_eV_A"] <= 0.1 and error <= GATE),
        "interpretation": "The existing aggregate periodic development gates pass; the held-out Si hotspot center gate fails and needs center-force training data.",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
