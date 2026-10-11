"""Hash-bound continuation entry point for D task step 1300 -> 2000.

The default path is read-only preflight. Scientific integration requires an
explicit --execute flag and a C-signed config whose SHA and all source pins
match. This module never initializes atoms, velocities, or RNG state.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath

import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "core"))

import github_archive_backend as archive_backend
import segment_core as core

TASK_ID = "C-D-AGL-1250K-rho097-seed01"
OWNER = "D"
START_STEP = 1300
STOP_STEP = 2000
MAX_ADDED_STEPS = STOP_STEP - START_STEP
REPOSITORY = "https://github.com/linkengde/-.git"
D_BRANCH = "cloud-d-b1-dft"
D_BRANCH_HEAD = "9aea23179dd73da749d67f7be0b3a3b86e5df085"
ARCHIVE_RECEIPT_BASE_HEAD = "8269f3f70b92a758507d71644c0c60f53db3cb59"
ARCHIVE_PATH = (
    "codex_segment_archives/D/C-D-AGL-1250K-rho097-seed01/"
    "step-000001300-e75322b9bdd4b6bc2d485225a76c889536219ac22b6fe299858f96c1411839ee.tar.gz"
)
ARCHIVE_SHA256 = "e75322b9bdd4b6bc2d485225a76c889536219ac22b6fe299858f96c1411839ee"
ARCHIVE_RECEIPT_SHA256 = "ede5dfee81846911e289540391432c5dfe324aac159cd89aa6dc08cafb93335f"
STATE_SHA256 = "3b3950729735765dc5c4d205d16a3654b4e9f8ba193c58bf8421685f36d66b9e"
ARRAYS_SHA256 = "6eb03c97a33f6c2f3b5d9521c98435a77e5ee9228286090602733c3d22a65a2a"
OLD_CORE_SHA256 = "e0bba3702c54bbdcb39f5f7e5641e1e1969a2ea5af8faceee8a0d50316f440dd"
OLD_RELEASE_SHA256 = "026be260ca65ab55b60c96a407c6181649e5f245bcf2d138993e246d2aa56d9b"
POTENTIAL_SHA256 = core.EXPECTED_POTENTIAL_SHA256
ARCHIVE_BACKEND_SHA256 = "8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346"
ARCHIVE_TEST_SHA256 = "1f9f5dc0536b2736cf8b05583631bf7c3dcd3f595f73db38a3534d97ff33c924"
SOURCE_EXECUTION_STATUS = "SYNTHETIC_OR_UNSIGNED_CANDIDATE_ONLY"


class ContinuationRefused(RuntimeError):
    pass


def _require(condition: bool, reason: str) -> None:
    if not condition:
        raise ContinuationRefused(reason)


def sha(path: str | Path) -> str:
    return archive_backend.sha256_file(path)


def validate_release_config(path: Path, expected_sha256: str) -> dict:
    _require(path.is_file(), "SIGNED_CONTINUATION_CONFIG_MISSING")
    _require(sha(path) == expected_sha256, "SIGNED_CONTINUATION_CONFIG_SHA_MISMATCH")
    adapter_sha = sha(Path(__file__))
    core_sha = sha(Path(core.__file__))
    backend_sha = sha(Path(archive_backend.__file__))
    cfg = json.loads(path.read_text(encoding="utf-8"))
    prior_release = HERE / "source_step1000_release_reference.json"
    _require(prior_release.is_file() and sha(prior_release) == OLD_RELEASE_SHA256,
             "PRIOR_STEP1000_SIGNED_RELEASE_REFERENCE_SHA_MISMATCH")
    prior = json.loads(prior_release.read_text(encoding="utf-8"))
    _require(prior.get("execution_authorized") is True
             and prior.get("task_id") == TASK_ID
             and prior.get("start_step") == 1000
             and prior.get("stop_at_step") == 2000,
             "PRIOR_STEP1000_RELEASE_REFERENCE_NOT_AUTHORIZED")
    pins = {
        "config_type": "SIGNED_CONTINUATION_CONFIG",
        "release_status": "C_SIGNED_FOR_CONTINUATION_1300_2000",
        "execution_authorized": True,
        "task_id": TASK_ID,
        "owner": OWNER,
        "start_step": START_STEP,
        "stop_at_step": STOP_STEP,
        "max_added_steps": MAX_ADDED_STEPS,
        "continuation_adapter_sha256": adapter_sha,
        "continuation_core_sha256": core_sha,
        "archive_backend_sha256": backend_sha,
        "source_archive_sha256": ARCHIVE_SHA256,
        "source_archive_receipt_sha256": ARCHIVE_RECEIPT_SHA256,
        "source_checkpoint_state_sha256": STATE_SHA256,
        "source_checkpoint_arrays_sha256": ARRAYS_SHA256,
        "source_executor_code_sha256": OLD_CORE_SHA256,
        "prior_step1000_release_config_sha256": OLD_RELEASE_SHA256,
        "potential_sha256": POTENTIAL_SHA256,
        "archive_backend_test_sha256": ARCHIVE_TEST_SHA256,
        "d_branch": D_BRANCH,
        "d_branch_head_at_release": D_BRANCH_HEAD,
        "source_archive_receipt_base_head": ARCHIVE_RECEIPT_BASE_HEAD,
        "remote_archive_verified": True,
        "remote_archive_method": "GITHUB_PUSH_FRESH_FETCH_SHA256",
        "prohibited_step_gt_2000": True,
        "automatic_continue_to_120000": False,
        "velocity_or_rng_reinitialization_allowed": False,
    }
    for key, expected in pins.items():
        _require(cfg.get(key) == expected, "RELEASE_PIN_MISMATCH:" + key)
    _require(cfg.get("sampling") == core.SIGNED_SAMPLING_PARAMETERS,
             "SCIENTIFIC_OR_RESOURCE_PARAMETERS_MISMATCH")
    _require(cfg.get("source_checkpoint_status") == SOURCE_EXECUTION_STATUS,
             "SOURCE_CHECKPOINT_STATUS_NOT_ACKNOWLEDGED")
    _require(cfg.get("source_checkpoint_status_disposition") ==
             "RETAIN_ORIGINAL_METADATA; TRUST_CHAIN_BINDS_PRIOR_C_SIGNED_RELEASE_AND_D_AUDIT",
             "SOURCE_CHECKPOINT_STATUS_DISPOSITION_MISMATCH")
    _require(cfg.get("max_wall_s") == 600,
             "WALLTIME_BUDGET_CHANGED_WITHOUT_RELEASE")
    _require(cfg.get("new_output_root") ==
             "/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step1300-2000",
             "OUTPUT_ROOT_PIN_MISMATCH")
    _require(cfg.get("release_signature", {}).get("signed_by_c") == "Cloud C (Codex)"
             and cfg.get("release_signature", {}).get("signature_semantics") ==
             "Procedural C release; immutable Git commit and config SHA provide audit integrity, not a cryptographic identity signature.",
             "C_RELEASE_SIGNATURE_METADATA_INVALID")
    return cfg


def _clone_pinned_d_branch(tmp: Path, cfg: dict) -> Path:
    remote_head = subprocess.run(
        ["git", "ls-remote", REPOSITORY, f"refs/heads/{D_BRANCH}"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env={**os.environ, "GIT_TERMINAL_PROMPT": "0", "GIT_OPTIONAL_LOCKS": "0"},
        timeout=120, check=False,
    )
    _require(remote_head.returncode == 0, "D_REMOTE_HEAD_READ_FAILED")
    rows = [line.split() for line in remote_head.stdout.splitlines() if line.strip()]
    _require(len(rows) == 1 and len(rows[0]) == 2 and rows[0][0] == cfg["d_branch_head_at_release"],
             "D_BRANCH_HEAD_MOVED_SINCE_RELEASE")
    clone = tmp / "d-source-fresh-clone"
    result = subprocess.run(
        ["git", "clone", "--quiet", "--depth", "1", "--single-branch", "--no-tags",
         "--branch", D_BRANCH, REPOSITORY, str(clone)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env={**os.environ, "GIT_TERMINAL_PROMPT": "0", "GIT_OPTIONAL_LOCKS": "0"},
        timeout=180, check=False,
    )
    _require(result.returncode == 0, "D_BRANCH_FRESH_CLONE_FAILED")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=clone, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    _require(head.returncode == 0 and head.stdout.strip() == cfg["d_branch_head_at_release"],
             "D_FRESH_CLONE_HEAD_MISMATCH")
    return clone


def _extract_bundle_verified(archive_path: Path, bundle_parent: Path) -> tuple[Path, dict]:
    manifest, files = archive_backend._read_archive_members(archive_path, ARCHIVE_SHA256)
    _require(manifest.get("task_id") == TASK_ID and manifest.get("owner") == OWNER
             and manifest.get("global_step") == START_STEP,
             "SOURCE_ARCHIVE_CONTENT_IDENTITY_MISMATCH")
    roots = {PurePosixPath(name).parts[0] for name in files}
    _require(len(roots) == 1, "SOURCE_ARCHIVE_ROOT_DIRECTORY_AMBIGUOUS")
    root_name = next(iter(roots))
    bundle = bundle_parent / root_name
    _require(not bundle.exists(), "SOURCE_BUNDLE_TEMP_PATH_EXISTS")
    bundle.mkdir(parents=True, exist_ok=False)
    for archive_name, content in files.items():
        rel = PurePosixPath(archive_name)
        _require(not rel.is_absolute() and ".." not in rel.parts and rel.parts[0] == root_name,
                 "SOURCE_ARCHIVE_PATH_UNSAFE")
        if len(rel.parts) == 1:
            continue
        target = bundle.joinpath(*rel.parts[1:])
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    state_path = bundle / "state.json"
    arrays_path = bundle / "state_arrays.npz"
    _require(sha(state_path) == STATE_SHA256, "SOURCE_STATE_JSON_SHA_MISMATCH")
    _require(sha(arrays_path) == ARRAYS_SHA256, "SOURCE_ARRAYS_SHA_MISMATCH")
    restored = core.restore_segment_bundle_with_executor(bundle, OLD_CORE_SHA256)
    _require(restored.global_step == START_STEP and len(restored.atoms) == 64,
             "SOURCE_CHECKPOINT_STEP_OR_ATOM_COUNT_MISMATCH")
    _require(set(restored.rng_state_all) == {"packing", "velocity", "thermostat"},
             "SOURCE_RNG_STREAM_SET_MISMATCH")
    state = restored.sampler_state
    _require(state.get("execution_status") == SOURCE_EXECUTION_STATUS,
             "SOURCE_EXECUTION_STATUS_CHANGED")
    provenance = state.get("source_provenance", {})
    _require(provenance.get("execution_release_config_sha256") == OLD_RELEASE_SHA256
             and provenance.get("synthetic") is False
             and provenance.get("potential_sha256") == POTENTIAL_SHA256
             and provenance.get("source_step") == 1000,
             "SOURCE_SIGNED_ANCESTRY_OR_POTENTIAL_PROVENANCE_MISMATCH")
    return bundle, {"restored": restored, "manifest": manifest}


def verify_source(path: Path, cfg: dict, tmp: Path) -> tuple[Path, dict, Path]:
    _require(sha(path) == cfg["source_archive_sha256"], "SOURCE_ARCHIVE_SHA_MISMATCH")
    _require(path.stat().st_size <= archive_backend.MAX_ARCHIVE_BYTES,
             "SOURCE_ARCHIVE_EXCEEDS_50_MIB")
    sidecar = path.with_suffix(path.suffix + ".sha256")
    receipt_path = path.with_suffix(path.suffix + ".receipt.json")
    _require(sidecar.is_file() and sidecar.read_text(encoding="utf-8").split()[0] == ARCHIVE_SHA256,
             "SOURCE_ARCHIVE_SIDECAR_MISMATCH")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    _require(sha(receipt_path) == ARCHIVE_RECEIPT_SHA256,
             "SOURCE_ARCHIVE_RECEIPT_SHA_MISMATCH")
    _require(receipt.get("archive_sha256") == ARCHIVE_SHA256
             and receipt.get("global_step") == START_STEP
             and receipt.get("task_id") == TASK_ID and receipt.get("owner") == OWNER
             and receipt.get("archive_path") == ARCHIVE_PATH
             and receipt.get("base_head") == ARCHIVE_RECEIPT_BASE_HEAD,
             "SOURCE_ARCHIVE_RECEIPT_IDENTITY_MISMATCH")
    bundle, restored_evidence = _extract_bundle_verified(path, tmp)
    return bundle, restored_evidence, receipt_path


def _ensure_new_output_root(output_root: Path) -> None:
    _require(not output_root.exists(), "OUTPUT_ROOT_ALREADY_EXISTS_NO_OVERWRITE")
    probe = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=HERE,
                           text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           check=False)
    if probe.returncode == 0:
        checkout = Path(probe.stdout.strip()).resolve()
        target = output_root.resolve()
        _require(target != checkout and checkout not in target.parents,
                 "OUTPUT_ROOT_MUST_BE_OUTSIDE_CHECKOUT")


def preflight(config_path: Path, config_sha256: str, potential_path: Path,
              output_root: Path) -> dict:
    cfg = validate_release_config(config_path, config_sha256)
    _require(output_root.resolve() == Path(cfg["new_output_root"]).resolve(),
             "OUTPUT_ROOT_DIFFERS_FROM_SIGNED_CONFIG")
    _require(potential_path.is_file() and sha(potential_path) == POTENTIAL_SHA256,
             "POTENTIAL_FILE_MISSING_OR_SHA_MISMATCH")
    _ensure_new_output_root(output_root)
    with tempfile.TemporaryDirectory(prefix="codex-c-step1300-preflight-") as temp:
        tmp = Path(temp)
        clone = _clone_pinned_d_branch(tmp, cfg)
        source = clone / ARCHIVE_PATH
        bundle, evidence, _ = verify_source(source, cfg, tmp)
        restored = evidence["restored"]
        return {
            "status": "PREFLIGHT_PASS_NO_SCIENCE",
            "task_id": TASK_ID,
            "start_step": START_STEP,
            "stop_step": STOP_STEP,
            "source_archive_sha256": sha(source),
            "source_state_sha256": restored.state_sha256,
            "source_arrays_sha256": restored.arrays_sha256,
            "restored_global_step": restored.global_step,
            "atom_count": len(restored.atoms),
            "rng_streams": sorted(restored.rng_state_all),
            "thermal_rows": len(restored.thermal_rows or []),
            "trajectory_frames": len(restored.trajectory_frames or []),
            "output_root_absent": True,
            "potential_sha256": sha(potential_path),
            "release_config_sha256": config_sha256,
            "verified_source_bundle": str(bundle),
        }


def run_continuation(config_path: Path, config_sha256: str, potential_path: Path,
                     output_root: Path) -> dict:
    cfg = validate_release_config(config_path, config_sha256)
    _require(output_root.resolve() == Path(cfg["new_output_root"]).resolve(),
             "OUTPUT_ROOT_DIFFERS_FROM_SIGNED_CONFIG")
    _require(potential_path.is_file() and sha(potential_path) == POTENTIAL_SHA256,
             "POTENTIAL_FILE_MISSING_OR_SHA_MISMATCH")
    _ensure_new_output_root(output_root)
    with tempfile.TemporaryDirectory(prefix="codex-c-step1300-execution-") as temp:
        tmp = Path(temp)
        clone = _clone_pinned_d_branch(tmp, cfg)
        source_archive = clone / ARCHIVE_PATH
        bundle_dir, evidence, _ = verify_source(source_archive, cfg, tmp)
        restored = evidence["restored"]
        _require(restored.global_step == START_STEP, "CONTINUATION_MUST_START_AT_1300")
        from ase.calculators.eam import EAM
        restored.atoms.calc = EAM(potential=str(potential_path))
        dyn = core.BoundedLangevin(
            restored.atoms, restored.thermostat_rng,
            start_step=1000, stop_at_step=STOP_STEP, initial_step=START_STEP,
            continuation_start_step=START_STEP,
            release_lock_path=config_path, release_lock_sha256=config_sha256,
            execution_code_sha256=sha(Path(core.__file__)),
            continuation_adapter_path=Path(__file__),
        )
        backend = archive_backend.GitHubArchiveBackend(
            remote_url=REPOSITORY, branch=D_BRANCH,
            expected_head=cfg["d_branch_head_at_release"],
            archive_prefix="codex_segment_archives",
            author_name="Codex Cloud D",
            author_email="codex-cloud-d@users.noreply.github.com",
        )
        _require(backend._remote_head() == cfg["d_branch_head_at_release"],
                 "D_BRANCH_HEAD_MOVED_BEFORE_EXECUTION")
        receipts: list[dict] = []

        def remote_roundtrip(archive_path: Path, archive_sha: str) -> None:
            receipt = backend.upload(archive_path, archive_sha)
            verified = backend.verify(receipt)
            restore_dir = output_root / "remote_restore" / f"step-{receipt['global_step']:09d}-{archive_sha[:16]}"
            restored_receipt = backend.restore(receipt, restore_dir)
            remote_bundle = restore_dir
            recovered = core.restore_segment_bundle(remote_bundle)
            local_bundle = archive_path.parent.parent / f"segment-step-{receipt['global_step']:09d}"
            local = core.restore_segment_bundle(local_bundle)
            local_manifest = json.loads((local_bundle / "MANIFEST.json").read_text(encoding="utf-8"))
            _require(recovered.global_step == receipt["global_step"], "REMOTE_RESTORE_STEP_MISMATCH")
            _require(recovered.arrays_sha256 == local_manifest["members_sha256"]["state_arrays.npz"],
                     "REMOTE_RESTORE_STATE_ARRAYS_MISMATCH")
            _require(recovered.state_sha256 == local_manifest["members_sha256"]["state.json"],
                     "REMOTE_RESTORE_STATE_JSON_MISMATCH")
            _require(np.array_equal(recovered.atoms.get_positions(), local.atoms.get_positions())
                     and np.array_equal(recovered.atoms.get_momenta(), local.atoms.get_momenta())
                     and np.array_equal(recovered.atoms.get_velocities(), local.atoms.get_velocities()),
                     "REMOTE_RESTORE_ATOM_STATE_MISMATCH")
            _require(recovered.rng_state_all == json.loads((remote_bundle / "rng_state.json").read_text(encoding="utf-8")),
                     "REMOTE_RESTORE_RNG_STATE_MISMATCH")
            _require(recovered.rng_state_all == local.rng_state_all,
                     "REMOTE_RESTORE_RNG_DIFFERS_FROM_LOCAL_CHECKPOINT")
            receipts.append({**receipt, "verify": verified, "restore": restored_receipt})

        output_root.mkdir(parents=True, exist_ok=False)
        provenance = {
            "source_archive_sha256": ARCHIVE_SHA256,
            "source_checkpoint_state_sha256": STATE_SHA256,
            "source_checkpoint_arrays_sha256": ARRAYS_SHA256,
            "source_step": START_STEP,
            "source_d_branch_head": D_BRANCH_HEAD,
            "source_executor_code_sha256": OLD_CORE_SHA256,
            "source_execution_release_config_sha256": OLD_RELEASE_SHA256,
            "production_adapter_sha256": sha(Path(__file__)),
            "segment_core_sha256": sha(Path(core.__file__)),
            "archive_backend_sha256": sha(Path(archive_backend.__file__)),
            "archive_backend_test_sha256": ARCHIVE_TEST_SHA256,
            "execution_release_config_sha256": config_sha256,
            "potential_sha256": POTENTIAL_SHA256,
            "execution_authorized": True,
            "synthetic": False,
        }
        outputs = core.advance_segment_with_outputs(
            dyn, rng_state_all=restored.rng_state_all,
            source_provenance=provenance, output_root=output_root,
            remote_roundtrip=remote_roundtrip,
            existing_thermal_rows=restored.thermal_rows,
            existing_trajectory_frames=restored.trajectory_frames,
            checkpoint_every=100, max_wall_s=cfg["max_wall_s"],
        )
        result = {
            "status": "CONTINUATION_COMPLETE_REMOTE_ARCHIVE_VERIFIED",
            "task_id": TASK_ID, "start_step": START_STEP,
            "stop_step": dyn.nsteps, "added_steps": dyn.nsteps - START_STEP,
            "config_sha256": config_sha256,
            "continuation_adapter_sha256": sha(Path(__file__)),
            "continuation_core_sha256": sha(Path(core.__file__)),
            "archive_backend_sha256": sha(Path(archive_backend.__file__)),
            "outputs": outputs, "remote_receipts": receipts,
        }
        result_path = output_root / "REMOTE_ARCHIVE_RECEIPTS.json"
        temp_path = result_path.with_suffix(".json.tmp")
        temp_path.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        os.replace(temp_path, result_path)
        return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-config", type=Path, required=True)
    parser.add_argument("--release-config-sha256", required=True)
    parser.add_argument("--potential", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path(
        "/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step1300-2000"))
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    try:
        _require(args.preflight != args.execute, "EXACTLY_ONE_OF_PREFLIGHT_OR_EXECUTE_REQUIRED")
        if args.preflight:
            print(json.dumps(preflight(args.release_config, args.release_config_sha256,
                                       args.potential, args.output_root), sort_keys=True))
            return 0
        result = run_continuation(args.release_config, args.release_config_sha256,
                                  args.potential, args.output_root)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ContinuationRefused, core.SegmentRefused,
            archive_backend.ArchiveBackendError, OSError, ValueError,
            KeyError, subprocess.SubprocessError) as exc:
        print(json.dumps({"status": "HOLD", "reason": str(exc)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
