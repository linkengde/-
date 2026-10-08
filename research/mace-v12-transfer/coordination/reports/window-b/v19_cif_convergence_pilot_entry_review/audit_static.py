#!/usr/bin/env python3
"""Read-only static audit of the disabled V19 matched-kmesh pilot entry."""
from __future__ import annotations

import ast
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = next(
    parent
    for parent in HERE.parents
    if (parent / "research/mace-v12-transfer/coordination").is_dir()
)
COORD = REPO / "research/mace-v12-transfer/coordination"
PILOT = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot"
SOURCE = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures"
sha256 = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

checks: dict[str, object] = {}
inventory = {}
for line in (PILOT / "SHA256SUMS.txt").read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    inventory[name.strip()] = digest
inventory_results = {}
for name, expected in inventory.items():
    path = PILOT / name
    observed = sha256(path) if path.is_file() else None
    inventory_results[name] = {"expected": expected, "observed": observed, "pass": expected == observed}
checks["pilot_inventory"] = inventory_results

source_manifest = SOURCE / "manifest.json"
source_manifest_sha = sha256(source_manifest)
manifest = json.loads((PILOT / "input_manifest.json").read_text())
task = json.loads((COORD / "tasks/window-b.json").read_text())
task_cards = {
    owner: json.loads((COORD / "tasks" / f"{owner}.json").read_text())
    for owner in ("window-a", "window-b")
}
identity = subprocess.check_output(
    [sys.executable, str(COORD / "sync_tasks.py"), "identity"], text=True
).strip()
records = {record["owner_task"]: record for record in manifest["records"]}
checks["source_manifest_hash"] = {
    "declared": manifest["source_manifest_sha256"],
    "observed": source_manifest_sha,
    "pass": manifest["source_manifest_sha256"] == source_manifest_sha,
}
checks["owner_and_input_records"] = {
    owner: {
        "label": record["label"],
        "owner_instance_matches_task_card": record["owner_instance"] == task_cards[owner]["owner_instance"],
        "owner_instance_matches_this_runtime": owner == "window-b" and record["owner_instance"] == identity,
        "label_registered_with_owner_task": record["label"] in task_cards[owner].get("labels", {}),
        "input_hash_matches_record": sha256(PILOT / record["input"]) == record["input_sha256"],
        "input_hash_matches_CIF_AgSi_proposal": record["input_sha256"]
        == sha256(SOURCE / "inputs/AgSi_cod9009647_slab_v19_dev_proposal.extxyz"),
        "kpts": record["kpts"],
        "blocked_state": record["state"],
    }
    for owner, record in records.items()
}
checks["launch_disabled"] = manifest["launch_enabled"] is False
checks["pilot_method_and_geometry_control"] = {
    "cutoff_eV": manifest["cutoff_eV"],
    "fermi_smearing_eV": manifest["fermi_smearing_eV"],
    "MPI_ranks": manifest["MPI_ranks"],
    "matched_geometry": manifest["matched_geometry"],
    "pilot_only": manifest["pilot_only"],
    "same_input_hash_for_both_meshes": len({r["input_sha256"] for r in manifest["records"]}) == 1,
}

queue_text = (PILOT / "run_queue.py").read_text()
calc_text = (PILOT / "run_pw_reference.py").read_text()
verify_text = (PILOT / "verify_result.py").read_text()
compare_text = (PILOT / "compare_mesh_results.py").read_text()
python_files = [PILOT / name for name in (
    "run_queue.py", "run_pw_reference.py", "verify_result.py", "compare_mesh_results.py"
)]
syntax = {}
for path in python_files:
    try:
        ast.parse(path.read_text(), filename=str(path))
        syntax[path.name] = True
    except SyntaxError as error:
        syntax[path.name] = {"line": error.lineno, "message": error.msg}
checks["python_syntax"] = syntax

disabled_test = "if not manifest[\"launch_enabled\"]:" in calc_text
checks["calculator_disabled_guard_precedes_GPAW_import"] = {
    "pass": disabled_test and calc_text.index("if not manifest[\"launch_enabled\"]:") < calc_text.index("import gpaw_data"),
    "evidence": "run_pw_reference.py checks launch_enabled before importing gpaw_data/GPAW",
}
claim_pos = queue_text.index("sync.claim(TASK)")
manifest_pos = queue_text.index('manifest=json.loads((ROOT/\'input_manifest.json\')')
enabled_pos = queue_text.index("manifest['launch_enabled'] is True")
checks["queue_checks_disabled_before_claim"] = {
    "pass": enabled_pos < claim_pos,
    "observed_order": ["sync.claim(TASK)", "load manifest", "assert launch_enabled"],
    "impact": "A disabled queue invocation mutates/claims the task card before refusing launch.",
}
checks["queue_checks_current_owner"] = "sync.ensure_owner(TASK,task)" in queue_text
checks["direct_calculator_checks_current_owner"] = bool(
    re.search(r"owner_task|owner_instance|sync\.identity", calc_text)
)
checks["record_bound_kpoints_in_calculator"] = "kpts=tuple(record[\"kpts\"])" in calc_text
checks["verifier_checks_recorded_kpoints"] = "s['method']['kpts']==r['kpts']" in verify_text
checks["archive_inventory_validated_before_rewrite"] = bool(
    re.search(r"SHA256SUMS\.txt.*(check|verify)|sha256sum\s+-c", queue_text, re.I)
)
checks["stale_output_refused_when_archive_exists"] = "if out.exists()" in queue_text[:queue_text.index("if archive.exists()")]
state_path_comment = re.search(r"state_path\s*=\s*OUT\s*/\s*[\"']state\.gpw[\"']\s*#([^\n]+)", calc_text)
checks["checkpoint_comment_matches_path"] = {
    "pass": not state_path_comment or "/tmp" not in state_path_comment.group(1),
    "observed_path": "OUT/state.gpw" if state_path_comment else None,
    "observed_comment": state_path_comment.group(1).strip() if state_path_comment else None,
}
checks["runtime_disk_monitor_present"] = bool(
    re.search(r"while child\.poll\(\) is None:[\s\S]{0,500}disk_usage", queue_text)
)
checks["comparison_reports_required_metrics"] = {
    key: key in compare_text
    for key in (
        "native_energy_delta_meV_atom", "free_energy_delta_meV_atom",
        "force_vector_difference_RMSE_eV_A", "maximum_atom_force_difference_eV_A",
        "pair_signed_projection_delta_eV_A",
    )
}

workspace = Path("/workspace")
workspace_free = shutil.disk_usage(workspace).free
out = workspace / "mace_v19_cif_pilot_window-b"
archive = PILOT / "calculations" / records["window-b"]["label"]
lock_path = workspace / ".mace-v19-dft.lock"
process_lines = subprocess.check_output(["ps", "-eo", "args="], text=True).splitlines()
active = [
    line.strip() for line in process_lines
    if any(token in line for token in ("run_queue.py window-b", "run_pw_reference.py", "mpirun", "prterun", "gpaw"))
]
lock_entries = []
if lock_path.exists():
    stat = lock_path.stat()
    device = f"{os.major(stat.st_dev):02x}:{os.minor(stat.st_dev):02x}"
    lock_entries = [
        line.strip() for line in Path("/proc/locks").read_text().splitlines()
        if device in line and str(stat.st_ino) in line.split()
    ]
input_header = (PILOT / records["window-b"]["input"]).read_text().splitlines()[1]
cell_values = [float(value) for value in re.search(r'Lattice="([^"]+)"', input_header).group(1).split()]
cell = [cell_values[index:index + 3] for index in (0, 3, 6)]
volume = abs(
    cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
    - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
    + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0])
)
bohr_A = 0.529177210903
cutoff_Ha = manifest["cutoff_eV"] / 27.211386245988
kmax_bohr = math.sqrt(2 * cutoff_Ha)
plane_waves = (volume / bohr_A**3) * kmax_bohr**3 / (6 * math.pi**2)
raw_coeff_MiB = {
    str(nbands): plane_waves * nbands * 4 * 16 / (1024**2)
    for nbands in (118, 260)
}
checks["workspace_and_old_run_snapshot"] = {
    "free_bytes": workspace_free,
    "B_run_output_exists": out.exists(),
    "B_archive_exists": archive.exists(),
    "persistent_lock_file_exists": lock_path.exists(),
    "persistent_lock_file_bytes": lock_path.stat().st_size if lock_path.exists() else None,
    "persistent_lock_file_held_entries": lock_entries,
    "matching_DFT_processes": active,
}
checks["disk_estimate_only"] = {
    "cell_volume_A3": volume,
    "estimated_plane_waves_per_kpoint_at_declared_cutoff": plane_waves,
    "raw_complex128_wavefunction_coefficients_MiB_at_4_kpoints_for_band_count": raw_coeff_MiB,
    "caveat": "Coefficient payload estimate only; excludes GPAW metadata, densities, restart/write overhead and filesystem headroom. Band count is a bracket, not measured from a GPAW run.",
}

result = {
    "audit": "v19_cif_convergence_pilot_entry_review",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "runtime_identity": identity,
    "registered_window_b_owner": task.get("owner_instance"),
    "pilot_launch_enabled": manifest["launch_enabled"],
    "checks": checks,
    "decision": "BLOCKED_PENDING_A_REVIEW_AND_JOB_REGISTRATION",
    "calculations_launched": False,
}
(HERE / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
