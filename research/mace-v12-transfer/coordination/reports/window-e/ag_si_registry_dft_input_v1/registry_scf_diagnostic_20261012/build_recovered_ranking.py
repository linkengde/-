#!/usr/bin/env python3
"""Assemble the four recovered 6x6 registry single-point records."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "RANKING_6x6_RECOVERED.json"

source_hashes = json.loads((ROOT / "SOURCE_SHA256.json").read_text())
diag_results = {
    "R1": ROOT / "R1_diagnostic_case/run_maxiter200_mixer_beta010/result.json",
    "R3": ROOT / "R3_diagnostic_case/run_maxiter200_mixer_beta010/result.json",
}
recovered_logs = json.loads((ROOT / "R2_R4_LOG_RECOVERY.json").read_text())
records = []

for registry in ("R1", "R2", "R3", "R4"):
    input_sha = source_hashes[f"{registry}/input.py"]
    structure_sha = source_hashes[f"{registry}/structure.xyz"]
    if registry in diag_results:
        result = json.loads(diag_results[registry].read_text())
        if not result.get("scf_converged") or result.get("geometry_modified"):
            raise SystemExit(f"{registry} diagnostic result is not converged with fixed geometry")
        record = {
            "registry_id": registry,
            "total_energy_eV": result["total_energy_eV"],
            "force_consistent_free_energy_eV": result["force_consistent_free_energy_eV"],
            "forces_eV_A": result["forces_eV_A"],
            "force_count": len(result["forces_eV_A"]),
            "forces_atom_order": "zero-based atom order in the original structure.xyz",
            "scf_steps": result["scf_steps"],
            "scf_converged": True,
            "diagnostic_recovered": True,
            "recovery_method": "INDEPENDENT_SCF_DIAGNOSTIC_COPY",
            "diagnostic_configuration": {
                "maxiter": result["maxiter"],
                "mixer": result["mixer"],
            },
            "original_input_sha256": input_sha,
            "original_structure_sha256": structure_sha,
            "original_input_modified": False,
            "original_structure_modified": False,
            "stress_eV_A3": result["stress_eV_A3"],
            "stress_status": "AVAILABLE_FROM_DIAGNOSTIC_COPY",
        }
    else:
        result = next(x for x in recovered_logs if x["registry_id"] == registry)
        if not result.get("scf_converged"):
            raise SystemExit(f"{registry} log does not establish SCF convergence")
        forces = result["forces_eV_A"]
        record = {
            "registry_id": registry,
            "total_energy_eV": result["total_energy_eV"],
            "force_consistent_free_energy_eV": result["force_consistent_free_energy_eV"],
            "forces_eV_A": [[f["fx_eV_A"], f["fy_eV_A"], f["fz_eV_A"]] for f in forces],
            "force_count": result["force_count"],
            "forces_atom_order": "zero-based atom order in the original structure.xyz",
            "scf_steps": result["scf_steps"],
            "scf_converged": True,
            "diagnostic_recovered": True,
            "recovery_method": "EXISTING_GPAW_LOG_RECOVERY_NO_RECALCULATION",
            "original_input_sha256": input_sha,
            "original_structure_sha256": structure_sha,
            "original_input_modified": False,
            "original_structure_modified": False,
            "stress_eV_A3": None,
            "stress_status": "STRESS_UNAVAILABLE",
        }
    if record["force_count"] != 20 or len(record["forces_eV_A"]) != 20:
        raise SystemExit(f"{registry} recovered force count is not 20")
    records.append(record)

ranking = sorted(records, key=lambda x: x["total_energy_eV"])
document = {
    "document_id": "RANKING_6x6_RECOVERED",
    "source_batch_commit": "ab364698dca879284be61c66f0340b02a6aaccd7",
    "scientific_settings": {
        "xc": "PBE", "basis": "PAW plane waves", "cutoff_eV": 500,
        "kpoint_mesh": [6, 6, 1], "smearing_eV": 0.1,
        "geometry": "fixed", "relax_atoms": False, "relax_cell": False,
    },
    "diagnostic_recovered": True,
    "diagnostic_note": "R1 and R3 use independent SCF diagnostic copies with maxiter=200 and Mixer beta=0.10; this does not modify or replace their original screening inputs or registry definitions. R2 and R4 are recovered from existing GPAW logs without recalculation.",
    "ranking_method": "ascending total_energy_eV; records retain registry ID and recovery method",
    "ranked_registry_ids_by_total_energy_eV": [x["registry_id"] for x in ranking],
    "results_in_registry_order": records,
    "eight_by_eight_run": False,
    "d0_scan_run": False,
    "geometry_or_registry_modified": False,
    "data_role_assigned": False,
}
RESULTS.write_text(json.dumps(document, indent=2) + "\n")
print(json.dumps({"output": str(RESULTS), "ranked_registry_ids": document["ranked_registry_ids_by_total_energy_eV"], "records": len(records)}, indent=2))
