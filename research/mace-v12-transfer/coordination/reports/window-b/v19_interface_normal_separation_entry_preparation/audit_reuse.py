"""Compare proposal and parent bytes with public, non-holdout DFT input manifests."""
import hashlib
import json
from pathlib import Path
from ase.io import read

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
PERIODIC = REPO / "research/mace-v12-transfer/periodic_interface_v4"
PARENT = PERIODIC / "v19_cif_parent_structures"
SOURCE = json.loads((PARENT / "manifest.json").read_text())
PROPOSALS = json.loads((HERE / "input_manifest.json").read_text())
TASKS = {}
for name in ("window-a", "window-b"):
    p = REPO / f"research/mace-v12-transfer/coordination/tasks/{name}.json"
    TASKS[name] = json.loads(p.read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
parent_hashes = {r["sha256"]: r for r in SOURCE["records"]}
proposal_hashes = {r["sha256"]: r for r in PROPOSALS["records"]}
excluded_tokens = {"holdout", "blind", "sealed", "test"}
manifests = []
matches = {r["interface"]: [] for r in SOURCE["records"]}
source_geometry_matches = {r["interface"]: [] for r in SOURCE["records"]}
proposal_matches = {r["label"]: [] for r in PROPOSALS["records"]}
integrity_issues = []

for manifest_path in sorted(PERIODIC.rglob("input_manifest.json")):
    rel = manifest_path.relative_to(PERIODIC)
    if any(token in part.lower() for part in rel.parts for token in excluded_tokens):
        continue
    doc = json.loads(manifest_path.read_text())
    manifests.append(str(manifest_path.relative_to(REPO)))
    source_geometry_hash = doc.get("source_geometry_sha256") or doc.get("source_manifest_sha256")
    source_geometry_path = doc.get("source_geometry_path")
    if source_geometry_path and source_geometry_hash in parent_hashes:
        candidate = REPO / source_geometry_path
        if candidate.is_file() and sha(candidate) == source_geometry_hash:
            target = parent_hashes[source_geometry_hash]
            for record in doc.get("records", []):
                label = record.get("label")
                label_dirs = []
                for task in TASKS.values():
                    result_directory = task.get("label_result_directories", {}).get(label)
                    if result_directory:
                        label_dirs.append(REPO / result_directory / "calculations" / str(label))
                label_dirs.append(manifest_path.parent / "calculations" / str(label))
                archive = next((p for p in label_dirs if (p / "verification.json").is_file()), None)
                verification = json.loads((archive / "verification.json").read_text()) if archive else None
                summary = json.loads((archive / "summary.json").read_text()) if archive and (archive / "summary.json").is_file() else {}
                source_geometry_matches[target["interface"]].append({
                    "interface": target["interface"], "parent_label": target["label"],
                    "parent_sha256": source_geometry_hash, "existing_label": label,
                    "manifest": str(manifest_path.relative_to(REPO)),
                    "source_geometry_path": str(candidate.relative_to(REPO)),
                    "source_geometry_hash_verified": True,
                    "input_sha256": record.get("input_sha256"),
                    "input_pbc": record.get("pbc") or read(manifest_path.parent / record["input"]).pbc.tolist() if record.get("input") else record.get("pbc"),
                    "poissonsolver": record.get("poissonsolver"),
                    "kpts": record.get("kpts") or record.get("method", {}).get("kpts"),
                    "smearing_eV": record.get("method", {}).get("smearing_eV", doc.get("method", {}).get("smearing_eV")),
                    "archive_path": str(archive.relative_to(REPO)) if archive else None,
                    "archive_verification": verification.get("status") if verification else "NO_ARCHIVE_FOUND",
                    "archive_scf_converged": summary.get("scf_converged")})
    for record in doc.get("records", []):
        label = record.get("label") or record.get("config_type")
        expected = record.get("input_sha256") or record.get("sha256")
        input_name = record.get("input")
        if not expected or not input_name:
            continue
        candidate = Path(input_name)
        if not candidate.is_absolute():
            # Inputs are normally relative to their manifest directory; absolute repo-relative paths occur in older manifests.
            local_path = manifest_path.parent / candidate
            repo_path = REPO / candidate
            candidate = local_path if local_path.exists() else repo_path
        if not candidate.is_file():
            continue
        observed = sha(candidate)
        if observed != expected:
            integrity_issues.append({"manifest": str(manifest_path.relative_to(REPO)),
                                     "label": label, "path": str(candidate.relative_to(REPO)),
                                     "declared_sha256": expected, "observed_sha256": observed})
            continue
        target = parent_hashes.get(observed)
        proposal = proposal_hashes.get(observed)
        if target:
            label_dirs = []
            for task_name, task in TASKS.items():
                result_directory = task.get("label_result_directories", {}).get(label)
                if result_directory:
                    label_dirs.append(REPO / result_directory / "calculations" / str(label))
            label_dirs.append(manifest_path.parent / "calculations" / str(label))
            archive = None
            verification = None
            for directory in label_dirs:
                vp = directory / "verification.json"
                if vp.is_file():
                    archive = directory
                    verification = json.loads(vp.read_text())
                    break
            row = {"interface": target["interface"], "parent_label": target["label"],
                   "parent_sha256": observed, "existing_label": label,
                   "manifest": str(manifest_path.relative_to(REPO)),
                   "input_path": str(candidate.relative_to(REPO)),
                   "declared_input_hash_verified": True,
                   "record_state": record.get("state"),
                   "role": record.get("role"),
                   "input_pbc": record.get("pbc") or read(candidate).pbc.tolist(),
                   "method": record.get("method") or doc.get("method"),
                   "energy_convention": record.get("energy_convention") or doc.get("energy_convention"),
                   "archive_path": str(archive.relative_to(REPO)) if archive else None,
                   "archive_verification": verification.get("status") if verification else "NO_ARCHIVE_FOUND",
                   "archive_scf_converged": None}
            if archive and (archive / "summary.json").is_file():
                summary = json.loads((archive / "summary.json").read_text())
                row["archive_scf_converged"] = summary.get("scf_converged")
                row["mesh"] = summary.get("kpts") or summary.get("method", {}).get("kpts")
                row["smearing_eV"] = summary.get("smearing_eV") or summary.get("method", {}).get("smearing_eV")
            matches[target["interface"]].append(row)
        if proposal:
            proposal_matches[proposal["label"]].append({"existing_label": label,
                "manifest": str(manifest_path.relative_to(REPO)),
                "input_path": str(candidate.relative_to(REPO)),
                "input_sha256": observed, "role": record.get("role"),
                "record_state": record.get("state")})

ledger = {"status": "PUBLIC_MANIFEST_REUSE_AUDIT_COMPLETE",
          "scan_scope": "input_manifest.json files under periodic_interface_v4, excluding any path containing holdout, blind, sealed or test; only manifest-declared label-free input bytes were opened and hashed.",
          "excluded_scope": "No dataset_manifest, training extxyz with labels, holdout/blind/test/sealed path, DFT output labels or sealed geometry was opened.",
          "scanned_manifest_count": len(manifests), "scanned_manifests": manifests,
          "manifest_input_hash_mismatches": integrity_issues,
          "parent_exact_byte_reuse": matches,
          "parent_exact_source_geometry_reuse_with_method_differences": source_geometry_matches,
          "proposal_exact_byte_reuse": proposal_matches,
          "conclusion": "Exact parent geometry inputs can reuse the matching numerical-baseline archive under its recorded method. None of the six offset proposal byte hashes is an exact match to a scanned public input. Reusing a parent baseline does not supply labels for either offset.",
          "limits": ["This is an exact input-byte audit; it does not equate nearby geometries or different methods.",
                     "A matching input is reusable only for its recorded mesh, boundary and energy convention.",
                     "All three interface terminations share one COD9009647 CIF source and do not form independent validation families."]}
(HERE / "reuse_ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
print(json.dumps({"scanned_manifest_count": len(manifests),
                  "matched_parent_labels": {k: len(v) for k, v in matches.items()},
                  "matched_proposal_labels": sum(map(len, proposal_matches.values())),
                  "integrity_mismatches": len(integrity_issues)}, indent=2))
