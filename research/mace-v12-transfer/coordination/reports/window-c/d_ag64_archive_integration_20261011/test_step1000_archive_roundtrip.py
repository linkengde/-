from __future__ import annotations

import argparse
import sys
import copy
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "core"))
import github_archive_backend as backend_module
import segment_core as core

TASK_ID = "C-D-AGL-1250K-rho097-seed01"
OWNER = "D"


def h_array(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Non-scientific step1000 copy/archive/restore check")
    parser.add_argument("--remote-url", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--expected-branch-head", required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    source_archive_sha_before = backend_module.sha256_file(args.source_archive)
    checkpoint_copy_tmp = tempfile.TemporaryDirectory(prefix="stage-d-step1000-checkpoint-copy-")
    checkpoint_copy = Path(checkpoint_copy_tmp.name) / "step1000-source-copy.tar.gz"
    shutil.copyfile(args.source_archive, checkpoint_copy)
    if backend_module.sha256_file(checkpoint_copy) != source_archive_sha_before:
        raise RuntimeError("STEP1000_CHECKPOINT_COPY_SHA_MISMATCH")
    source = core.verify_source_archive(checkpoint_copy, args.source_manifest)
    if source["archive_sha256"] != core.EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("SOURCE_STEP1000_ARCHIVE_BINDING_FAIL")
    restored = core.restore_source_from_archive(checkpoint_copy)
    if restored.global_step != 1000 or restored.atoms.calc is not None:
        raise RuntimeError("SOURCE_CHECKPOINT_NOT_CLEAN_STEP1000_COPY")

    source_positions = restored.atoms.get_positions().copy()
    source_velocities = restored.atoms.get_velocities().copy()
    source_rng = copy.deepcopy(restored.rng_state_all)

    # Synthetic state-only evolution. This is deliberately not ASE integration:
    # there is no calculator, no EAM load, no force/energy call and no MD step.
    fixture = np.random.Generator(np.random.PCG64(20261011))
    for _global_step in range(1001, 2001):
        draw = restored.thermostat_rng.standard_normal((64, 3))
        restored.atoms.set_positions(
            restored.atoms.get_positions() + 1.0e-8 * draw,
            apply_constraint=False)
        restored.atoms.set_momenta(
            restored.atoms.get_momenta() + 1.0e-8 * restored.atoms.get_masses()[:, None] * draw,
            apply_constraint=False)
        fixture.random(1)
    final_step = 2000
    final_positions = restored.atoms.get_positions().copy()
    final_velocities = restored.atoms.get_velocities().copy()
    final_rng_all = copy.deepcopy(restored.rng_state_all)
    final_rng_all["thermostat"] = copy.deepcopy(restored.thermostat_rng.bit_generator.state)
    if np.array_equal(final_positions, source_positions) or np.array_equal(final_velocities, source_velocities):
        raise RuntimeError("SYNTHETIC_OUTPUT_DID_NOT_CHANGE_FROM_SOURCE_COPY")

    provenance = {
        "test_type": "NON_SCIENTIFIC_SYNTHETIC_STATE_ONLY",
        "source_task_id": TASK_ID,
        "source_owner": OWNER,
        "source_step": 1000,
        "source_archive_sha256": source["archive_sha256"],
        "source_state_sha256": source["state_sha256"],
        "source_arrays_sha256": source["arrays_sha256"],
        "no_calculator_loaded": True,
        "no_energy_force_or_integrator_calls": True,
        "synthetic_state_steps": 1000,
    }
    snapshot = core.snapshot_state(restored.atoms, final_step, final_rng_all,
                                   restored.thermostat_rng, provenance)

    with tempfile.TemporaryDirectory(prefix="stage-d-archive-roundtrip-") as tmp:
        root = Path(tmp)
        local_output = root / "simulated-output"
        bundle = core.write_checkpoint_bundle(local_output, "synthetic-step1000-copy",
                                              snapshot,
                                              thermal_rows=[],
                                              trajectory_frames=[
                                                  {"step": 1000, "positions": source_positions},
                                                  {"step": 2000, "positions": final_positions},
                                              ])
        archive_path = root / "simulated-output-step2000.tar.gz"
        archive_sha = core.create_archive(bundle, archive_path)
        git_backend = backend_module.GitHubArchiveBackend(
            args.remote_url, args.branch, args.expected_branch_head,
            "research/mace-v12-transfer/coordination/reports/window-c/d_ag64_archive_integration_20261011/backend_roundtrip_tests",
            author_name="Codex Cloud C", author_email="codex-cloud-c@users.noreply.github.com")
        receipt = git_backend.upload(archive_path, archive_sha)
        verified = git_backend.verify(receipt)
        restore_dir = root / "fresh-github-restore"
        restore_result = git_backend.restore(receipt, restore_dir)
        restored_bundle = core.restore_segment_bundle(restore_dir)
        positions_ok = np.array_equal(restored_bundle.atoms.get_positions(), final_positions)
        velocities_ok = np.array_equal(restored_bundle.atoms.get_velocities(), final_velocities)
        rng_ok = (restored_bundle.rng_state_all == final_rng_all
                  and restored_bundle.thermostat_rng.bit_generator.state == final_rng_all["thermostat"])
        if not (positions_ok and velocities_ok and rng_ok):
            raise RuntimeError("FRESH_RESTORE_STATE_MISMATCH")
        if backend_module.sha256_file(args.source_archive) != source_archive_sha_before:
            raise RuntimeError("ORIGINAL_STEP1000_ARCHIVE_BYTES_CHANGED")
        report = {
            "schema_version": 1,
            "test_status": "ARCHIVE_BACKEND_PRODUCTION_PASS",
            "scientific_calculation_performed": False,
            "md_integrator_called": False,
            "calculator_loaded": False,
            "source_checkpoint_copy": {
                "task_id": TASK_ID, "owner": OWNER, "step": 1000,
                "copied_before_test": True,
                "source_archive_commit": "145ee8050cf95a3a1988e868c6b935e10538f723",
                "archive_backend_evidence_commit": "c1c9f866426e572e82ab88157654a015802707f8",
                "copy_sha256": source_archive_sha_before,
                "archive_sha256": source["archive_sha256"],
                "state_sha256": source["state_sha256"],
                "arrays_sha256": source["arrays_sha256"],
                "original_bytes_unchanged": True,
                "positions_sha256_before": h_array(source_positions),
                "velocities_sha256_before": h_array(source_velocities),
                "rng_state_before": source_rng,
            },
            "simulated_output": {
                "step": final_step,
                "state_label": "SYNTHETIC_NOT_SCIENTIFIC_NOT_MD",
                "positions_sha256_expected": h_array(final_positions),
                "velocities_sha256_expected": h_array(final_velocities),
                "rng_state_expected": final_rng_all,
            },
            "archive": {
                "sha256": archive_sha,
                "remote_branch": args.branch,
                "remote_commit": receipt["commit_sha"],
                "remote_path": receipt["archive_path"],
                "upload_pass": True,
                "fresh_fetch_verify": verified,
                "restore": restore_result,
                "positions_exact": positions_ok,
                "velocities_exact": velocities_ok,
                "rng_exact": rng_ok,
            },
            "code_sha256": {
                "production_adapter": backend_module.sha256_file(Path(__file__).with_name("production_segment_adapter.py")),
                "archive_backend": backend_module.sha256_file(Path(backend_module.__file__)),
                "segment_core": backend_module.sha256_file(Path(core.__file__)),
                "test_script": backend_module.sha256_file(Path(__file__)),
            },
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": report["test_status"],
                          "archive_sha256": archive_sha,
                          "remote_commit": receipt["commit_sha"],
                          "positions_exact": positions_ok,
                          "velocities_exact": velocities_ok,
                          "rng_exact": rng_ok,
                          "report": str(args.report)}, sort_keys=True))
    checkpoint_copy_tmp.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
