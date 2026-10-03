#!/usr/bin/env python3
"""Verify the per-run persistent archives and write a compact manifest."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
DEST = ROOT / "calculations"
labels = ["AgSi_d2p2", "AgSi_d2p8", "AgC_d2p0", "AgC_d2p4"]
records = []
for label in labels:
    path = DEST / label
    summary_path = path / "summary.json"
    summary = json.loads(summary_path.read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise SystemExit(f"Invalid or non-converged persistent result: {label}")
    names = [f"{label}_PW_PBE.extxyz", "summary.json", "gpaw.log", "progress.json"]
    if any(not (path / name).is_file() or (path / name).stat().st_size == 0 for name in names):
        raise SystemExit(f"Incomplete persistent archive: {label}")
    records.append({
        "label": label,
        "scf_iterations": summary["scf_iterations"],
        "energy_eV_cell": summary["energy_eV_cell"],
        "source_sha256": summary["source_sha256"],
        "extxyz_sha256": hashlib.sha256((path / f"{label}_PW_PBE.extxyz").read_bytes()).hexdigest(),
        "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        "archived_files": names,
        "state_gpw_archived": False,
    })
(ROOT / "archive_manifest.json").write_text(json.dumps({"records": records}, indent=2) + "\n")
print(json.dumps(records, indent=2))
