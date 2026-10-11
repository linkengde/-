"""Hash-bound continuation entry point from the archived D step-2000 checkpoint.

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
REPOSITORY = "https://github.com/linkengde/-.git"
D_BRANCH = "cloud-d-b1-dft"
SOURCE_ARCHIVE_PATH = (
    "codex_segment_archives/D/C-D-AGL-1250K-rho097-seed01/"
    "step-000002000-c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965.tar.gz"
)
SOURCE_ARCHIVE_SHA256 = "c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965"
SOURCE_ARCHIVE_RECEIPT_SHA256 = "c0835b9210cab6752c2b2d0f1eddf6c2fec4075d5fb5957c3a45240ba1dc3550"
SOURCE_ARCHIVE_RECEIPT_BASE_HEAD = "3ae9aa1038235d38d19b02218081e8a3fdaaaf28"
SOURCE_ARCHIVE_COMMIT = "bc3ff0631db3b7d78dcc17556befa98c2f454480"
SOURCE_CHECKPOINT_STATE_SHA256 = "5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a"
SOURCE_CHECKPOINT_ARRAYS_SHA256 = "0f7830ebeba5f9a83da77b15e153546b24fb3efaf0299de767d9b4fbe4639f78"
SOURCE_EXECUTOR_CODE_SHA256 = "3033c3234c4074358f8ef54a6783c3c4c9495549c12c93641f4f777f3451ed76"
SOURCE_RELEASE_CONFIG_SHA256 = "70ca4190ecda29bb8212ba5ab144c36c6c99285b1ef9307393ef19dd409ed52c"
POTENTIAL_SHA256 = core.EXPECTED_POTENTIAL_SHA256
SOURCE_EXECUTION_STATUS = "EXECUTED_UNDER_SIGNED_CONTINUATION_RELEASE"


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
    pins = {
        "config_type": "SIGNED_CONTINUATION_CONFIG",
        "release_status": "C_SIGNED_FOR_STEP2000_CONTINUATION",
        "execution_authorized": True,
        "task_id": TASK_ID,
        "owner": OWNER,
        "continuation_adapter_sha256": adapter_sha,
        "continuation_core_sha256": core_sha,
        "archive_backend_sha256": backend_sha,
        "source_archive_path": SOURCE_ARCHIVE_PATH,
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "source_archive_receipt_sha256": SOURCE_ARCHIVE_RECEIPT_SHA256,
        "source_checkpoint_state_sha256": SOURCE_CHECKPOINT_STATE_SHA256,
        "source_checkpoint_arrays_sha256": SOURCE_CHECKPOINT_ARRAYS_SHA256,
        "source_executor_code_sha256": SOURCE_EXECUTOR_CODE_SHA256,
        "source_archive_commit": SOURCE_ARCHIVE_COMMIT,
        "d_branch": D_BRANCH,
        "d_branch_head_at_release": SOURCE_ARCHIVE_COMMIT,
        "source_archive_receipt_base_head": SOURCE_ARCHIVE_RECEIPT_BASE_HEAD,
        "potential_sha256": POTENTIAL_SHA256,
        "segment_checkpoint_every": core.CONTINUATION_CHECKPOINT_INTERVAL,
        "prohibited_step_gt_end": True,
        "automatic_continue_to_120000": False,
        "velocity_or_rng_reinitialization_allowed": False,
    }
    for key, expected in pins.items():
        _require(cfg.get(key) == expected, "RELEASE_PIN_MISMATCH:" + key)
    core.validate_segment_bounds(cfg.get("start_step"), cfg.get("end_step"),
                                 cfg.get("max_added_steps"))
    _require(cfg.get("sampling") == core.SIGNED_SAMPLING_PARAMETERS,
             "SCIENTIFIC_OR_RESOURCE_PARAMETERS_MISMATCH")
    _require(cfg.get("source_checkpoint_status") == SOURCE_EXECUTION_STATUS,
             "SOURCE_CHECKPOINT_STATUS_NOT_ACKNOWLEDGED")
    _require(cfg.get("max_wall_s") == 3600,
             "WALLTIME_BUDGET_CHANGED_WITHOUT_RELEASE")
    _require(isinstance(cfg.get("new_output_root"), str)
             and Path(cfg["new_output_root"]).is_absolute(),
             "OUTPUT_ROOT_MUST_BE_SIGNED_ABSOLUTE_PATH")
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


def _extract_bundle_verified(archive_path: Path, bundle_parent: Path,
                             cfg: dict) -> tuple[Path, dict]:
    manifest, files = archive_backend._read_archive_members(
        archive_path, cfg["source_archive_sha256"])
    _require(manifest.get("task_id") == TASK_ID and manifest.get("owner") == OWNER
             and manifest.get("global_step") == cfg["start_step"],
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
    _require(sha(state_path) == cfg["source_checkpoint_state_sha256"],
             "SOURCE_STATE_JSON_SHA_MISMATCH")
    _require(sha(arrays_path) == cfg["source_checkpoint_arrays_sha256"],
             "SOURCE_ARRAYS_SHA_MISMATCH")
    restored = core.restore_segment_bundle_with_executor(
        bundle, cfg["source_executor_code_sha256"],
        expected_global_step=cfg["start_step"])
    _require(restored.global_step == cfg["start_step"] and len(restored.atoms) == 64,
             "SOURCE_CHECKPOINT_STEP_OR_ATOM_COUNT_MISMATCH")
    _require(set(restored.rng_state_all) == {"packing", "velocity", "thermostat"},
             "SOURCE_RNG_STREAM_SET_MISMATCH")
    state = restored.sampler_state
    _require(state.get("execution_status") == SOURCE_EXECUTION_STATUS,
             "SOURCE_EXECUTION_STATUS_CHANGED")
    provenance = state.get("source_provenance", {})
    _require(provenance.get("execution_release_config_sha256") == SOURCE_RELEASE_CONFIG_SHA256
             and provenance.get("synthetic") is False
             and provenance.get("potential_sha256") == POTENTIAL_SHA256
             and provenance.get("source_step") == 1900,
             "SOURCE_SIGNED_ANCESTRY_OR_POTENTIAL_PROVENANCE_MISMATCH")
    return bundle, {"restored": restored, "manifest": manifest}


def verify_source(path: Path, cfg: dict, tmp: Path) -> tuple[Path, dict, Path]:
    _require(sha(path) == cfg["source_archive_sha256"], "SOURCE_ARCHIVE_SHA_MISMATCH")
    _require(path.stat().st_size <= archive_backend.MAX_ARCHIVE_BYTES,
             "SOURCE_ARCHIVE_EXCEEDS_50_MIB")
    sidecar = path.with_suffix(path.suffix + ".sha256")
    receipt_path = path.with_suffix(path.suffix + ".receipt.json")
    _require(sidecar.is_file()
             and sidecar.read_text(encoding="utf-8").split()[0] == cfg["source_archive_sha256"],
             "SOURCE_ARCHIVE_SIDECAR_MISMATCH")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    _require(sha(receipt_path) == cfg["source_archive_receipt_sha256"],
             "SOURCE_ARCHIVE_RECEIPT_SHA_MISMATCH")
    _require(receipt.get("archive_sha256") == cfg["source_archive_sha256"]
             and receipt.get("global_step") == cfg["start_step"]
             and receipt.get("task_id") == TASK_ID and receipt.get("owner") == OWNER
             and receipt.get("archive_path") == cfg["source_archive_path"]
             and receipt.get("base_head") == cfg["source_archive_receipt_base_head"],
             "SOURCE_ARCHIVE_RECEIPT_IDENTITY_MISMATCH")
    bundle, restored_evidence = _extract_bundle_verified(path, tmp, cfg)
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
    with tempfile.TemporaryDirectory(prefix="codex-c-step2000-preflight-") as temp:
        tmp = Path(temp)
        clone = _clone_pinned_d_branch(tmp, cfg)
        source = clone / cfg["source_archive_path"]
        bundle, evidence, _ = verify_source(source, cfg, tmp)
        restored = evidence["restored"]
        return {
            "status": "PREFLIGHT_PASS_NO_SCIENCE",
            "task_id": TASK_ID,
            "start_step": cfg["start_step"],
            "end_step": cfg["end_step"],
            "max_added_steps": cfg["max_added_steps"],
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
    with tempfile.TemporaryDirectory(prefix="codex-c-step2000-execution-") as temp:
        tmp = Path(temp)
        clone = _clone_pinned_d_branch(tmp, cfg)
        source_archive = clone / cfg["source_archive_path"]
        bundle_dir, evidence, _ = verify_source(source_archive, cfg, tmp)
        restored = evidence["restored"]
        _require(restored.global_step == cfg["start_step"],
                 "CONTINUATION_RESTORED_STEP_DIFFERS_FROM_SIGNED_START_STEP")
        from ase.calculators.eam import EAM
        restored.atoms.calc = EAM(potential=str(potential_path))
        dyn = core.BoundedLangevin(
            restored.atoms, restored.thermostat_rng,
            start_step=cfg["start_step"], end_step=cfg["end_step"],
            max_added_steps=cfg["max_added_steps"], initial_step=cfg["start_step"],
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
            receipt = backend.upload(
                archive_path, archive_sha, start_step=cfg["start_step"],
                end_step=cfg["end_step"], max_added_steps=cfg["max_added_steps"],
                source_checkpoint_sha256=cfg["source_checkpoint_state_sha256"])
            verified = backend.verify(receipt)
            restore_dir = output_root / "remote_restore" / f"step-{receipt['global_step']:09d}-{archive_sha[:16]}"
            restored_receipt = backend.restore(receipt, restore_dir)
            remote_bundle = restore_dir
            recovered = core.restore_segment_bundle(
                remote_bundle, expected_global_step=receipt["global_step"])
            local_bundle = archive_path.parent.parent / f"segment-step-{receipt['global_step']:09d}"
            local = core.restore_segment_bundle(
                local_bundle, expected_global_step=receipt["global_step"])
            local_manifest = json.loads((local_bundle / "MANIFEST.json").read_text(encoding="utf-8"))
            _require(recovered.global_step == receipt["global_step"], "REMOTE_RESTORE_STEP_MISMATCH")
            _require(receipt.get("start_step") == cfg["start_step"]
                     and receipt.get("end_step") == cfg["end_step"]
                     and receipt.get("source_checkpoint_sha256") == cfg["source_checkpoint_state_sha256"]
                     and receipt.get("checkpoint_sha256") == local_manifest["members_sha256"]["state.json"],
                     "REMOTE_RECEIPT_SEGMENT_OR_CHECKPOINT_BINDING_MISMATCH")
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
            "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
            "source_checkpoint_state_sha256": SOURCE_CHECKPOINT_STATE_SHA256,
            "source_checkpoint_arrays_sha256": SOURCE_CHECKPOINT_ARRAYS_SHA256,
            "source_step": cfg["start_step"],
            "source_d_branch_head": SOURCE_ARCHIVE_COMMIT,
            "source_executor_code_sha256": SOURCE_EXECUTOR_CODE_SHA256,
            "source_execution_release_config_sha256": SOURCE_RELEASE_CONFIG_SHA256,
            "production_adapter_sha256": sha(Path(__file__)),
            "segment_core_sha256": sha(Path(core.__file__)),
            "archive_backend_sha256": sha(Path(archive_backend.__file__)),
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
            checkpoint_every=cfg["segment_checkpoint_every"], max_wall_s=cfg["max_wall_s"],
        )
        _require(dyn.nsteps == cfg["end_step"],
                 "SEGMENT_STOPPED_BEFORE_OR_AFTER_SIGNED_END_STEP")
        result = {
            "status": "CONTINUATION_COMPLETE_REMOTE_ARCHIVE_VERIFIED",
            "task_id": TASK_ID, "start_step": cfg["start_step"],
            "end_step": dyn.nsteps, "added_steps": dyn.nsteps - cfg["start_step"],
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
        "/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step2000-extension"))
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
