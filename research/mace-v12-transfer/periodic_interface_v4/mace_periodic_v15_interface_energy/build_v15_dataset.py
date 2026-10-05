#!/usr/bin/env python3
"""Build A-approved 29-frame V15 numerical-screening data without blind labels."""
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
DRAFT = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py"
EXPECTED_DRAFT_SHA256 = "47515d2dad01341faaa217fa59d264a88ca82411510360e1b99897a55984237c"
OUT = ROOT / "data"
if OUT.exists() and any(OUT.iterdir()):
    raise SystemExit("V15 data exists; preserve it and review its manifest")
if hashlib.sha256(DRAFT.read_bytes()).hexdigest() != EXPECTED_DRAFT_SHA256:
    raise SystemExit("Reviewed B builder hash changed; A review required")
spec = importlib.util.spec_from_file_location("reviewed_v15_builder", DRAFT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
OUT.mkdir(parents=True, exist_ok=True)
try:
    result = module.build(REPO, OUT, {}, exclude_unresolved=True)
    if result["split_sizes"] != {"train": 29, "valid": 2, "test": 3} or result["composition_rank"] != 4:
        raise RuntimeError("Unexpected approved V15 split or rank")
    manifest_path = OUT / "dataset_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["version"] = "V15 A-integrated 29-frame numerical screening dataset"
    manifest.pop("not_production_or_training_approval", None)
    manifest["screening_authorization"] = "A-approved limited numerical training/evaluation; physical-method consistency remains UNKNOWN; no production MD/TTM authorization"
    manifest["reviewed_builder_sha256"] = EXPECTED_DRAFT_SHA256
    manifest["integration_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    report_path = OUT / "verification_report.json"
    report = json.loads(report_path.read_text())
    report["production_dataset_modified"] = True
    report["integration_scope"] = "New V15 screening data only; historical datasets unchanged"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    (OUT / "SHA256SUMS.txt").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in sorted(OUT.iterdir()) if p.is_file() and p.name != "SHA256SUMS.txt"))
    print(json.dumps(result, indent=2))
except Exception as exc:
    (OUT / "diagnostics_failure.json").write_text(json.dumps({"error": str(exc)}, indent=2) + "\n")
    raise
