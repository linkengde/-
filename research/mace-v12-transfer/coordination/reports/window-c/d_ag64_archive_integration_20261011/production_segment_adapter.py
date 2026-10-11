"""Production entry point for the one signed D Ag64 1000->2000 segment.

Importing this module performs no calculation. The preflight mode is read-only.
The --execute path is unavailable unless a hash-bound C execution release is
provided and the pinned D branch head, source checkpoint, potential and code
all verify before ASE's integrator is constructed.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "core"))

import github_archive_backend as archive_backend
import segment_core as core

REPOSITORY = "https://github.com/linkengde/-.git"
EXPECTED_TASK_ID = "C-D-AGL-1250K-rho097-seed01"
EXPECTED_OWNER = "D"
EXPECTED_START = 1000
EXPECTED_STOP = 2000
EXPECTED_POTENTIAL_SHA256 = core.EXPECTED_POTENTIAL_SHA256
EXPECTED_D_BRANCH = "cloud-d-b1-dft"
EXPECTED_D_BRANCH_HEAD_AT_RELEASE = "c1c9f866426e572e82ab88157654a015802707f8"
EXPECTED_D_ARCHIVE_TEST_SHA256 = "1f9f5dc0536b2736cf8b05583631bf7c3dcd3f595f73db38a3534d97ff33c924"


def _sha(path: Path) -> str:
    return archive_backend.sha256_file(path)


def validate_execution_release(path: Path, expected_sha256: str,
                                expected_adapter_sha256: str,
                                expected_backend_sha256: str) -> dict:
    if not path.is_file() or _sha(path) != expected_sha256:
        raise core.SegmentRefused("EXECUTION_RELEASE_CONFIG_SHA_MISMATCH")
    actual_adapter = _sha(Path(__file__))
    actual_backend = _sha(Path(archive_backend.__file__))
    actual_core = _sha(Path(core.__file__))
    if expected_adapter_sha256 != actual_adapter:
        raise core.SegmentRefused("PRODUCTION_ADAPTER_SHA_MISMATCH")
    if expected_backend_sha256 != actual_backend:
        raise core.SegmentRefused("ARCHIVE_BACKEND_SHA_MISMATCH")
    cfg = json.loads(path.read_text(encoding="utf-8"))
    if cfg.get("release_status") != "C_SIGNED_FOR_SEGMENT_1000_2000" or cfg.get("execution_authorized") is not True:
        raise core.SegmentRefused("EXECUTION_RELEASE_NOT_C_SIGNED_OR_AUTHORIZED")
    pins = {
        "owner": EXPECTED_OWNER,
        "task_id": EXPECTED_TASK_ID,
        "start_step": EXPECTED_START,
        "stop_at_step": EXPECTED_STOP,
        "max_added_steps": EXPECTED_STOP - EXPECTED_START,
        "production_adapter_sha256": actual_adapter,
        "archive_backend_sha256": actual_backend,
        "segment_core_sha256": actual_core,
        "execution_code_sha256": actual_core,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "d_branch": EXPECTED_D_BRANCH,
        "d_branch_head_at_release": EXPECTED_D_BRANCH_HEAD_AT_RELEASE,
        "archive_backend_test_sha256": EXPECTED_D_ARCHIVE_TEST_SHA256,
        "remote_archive_verified": True,
        "remote_archive_method": "GITHUB_PUSH_FRESH_FETCH_SHA256",
    }
    for key, expected in pins.items():
        if cfg.get(key) != expected:
            raise core.SegmentRefused("EXECUTION_RELEASE_PIN_MISMATCH:" + key)
    if cfg.get("source_archive_sha256") != core.EXPECTED_ARCHIVE_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_SOURCE_ARCHIVE_MISMATCH")
    if cfg.get("source_archive_manifest_sha256") != core.EXPECTED_ARCHIVE_MANIFEST_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_SOURCE_MANIFEST_MISMATCH")
    if cfg.get("source_checkpoint_state_sha256") != core.EXPECTED_STATE_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_STATE_SHA_MISMATCH")
    if cfg.get("source_state_arrays_sha256") != core.EXPECTED_ARRAYS_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_ARRAYS_SHA_MISMATCH")
    if cfg.get("source_executor_code_sha256") != core.EXPECTED_OLD_RUNNER_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_ORIGINAL_RUNNER_SHA_MISMATCH")
    if cfg.get("source_scientific_parameters_sha256") != core.EXPECTED_SCIENTIFIC_PARAMETERS_SHA256:
        raise core.SegmentRefused("EXECUTION_RELEASE_SOURCE_IDENTITY_MISMATCH")
    if cfg.get("sampling") != core.SIGNED_SAMPLING_PARAMETERS:
        raise core.SegmentRefused("EXECUTION_RELEASE_SAMPLING_PARAMETERS_MISMATCH")
    if cfg.get("output_archive_policy", {}).get("max_member_or_archive_bytes") != archive_backend.MAX_ARCHIVE_BYTES:
        raise core.SegmentRefused("EXECUTION_RELEASE_ARCHIVE_LIMIT_MISMATCH")
    return cfg


def _ensure_outside_checkout(output_root: Path) -> None:
    result = archive_backend.subprocess.run(
        ["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True,
        stdout=archive_backend.subprocess.PIPE, stderr=archive_backend.subprocess.PIPE,
        check=False)
    if result.returncode == 0:
        checkout = Path(result.stdout.strip()).resolve()
        target = output_root.resolve()
        if target == checkout or checkout in target.parents:
            raise core.SegmentRefused("OUTPUT_ROOT_MUST_BE_OUTSIDE_GIT_CHECKOUT")


def run_segment(*, release_path: Path, release_sha256: str,
                adapter_sha256: str, backend_sha256: str,
                source_archive: Path, source_manifest: Path,
                potential_path: Path, output_root: Path,
                remote_url: str = REPOSITORY) -> dict:
    cfg = validate_execution_release(release_path, release_sha256,
                                     adapter_sha256, backend_sha256)
    if output_root.exists():
        raise core.SegmentRefused("OUTPUT_ROOT_ALREADY_EXISTS_NO_OVERWRITE")
    _ensure_outside_checkout(output_root)
    if archive_backend.sha256_file(potential_path) != EXPECTED_POTENTIAL_SHA256:
        raise core.SegmentRefused("POTENTIAL_SHA_MISMATCH")
    evidence = core.verify_source_archive(source_archive, source_manifest)
    if evidence["archive_sha256"] != cfg["source_archive_sha256"]:
        raise core.SegmentRefused("SOURCE_ARCHIVE_CONFIG_MISMATCH")
    restored = core.restore_source_from_archive(source_archive)
    if restored.global_step != EXPECTED_START or restored.atoms.calc is not None:
        raise core.SegmentRefused("SOURCE_CHECKPOINT_NOT_STEP1000_OR_CALCULATOR_ATTACHED")
    backend = archive_backend.GitHubArchiveBackend(
        remote_url=remote_url, branch=cfg["d_branch"],
        expected_head=cfg["d_branch_head_at_release"],
        archive_prefix=cfg["output_archive_policy"]["remote_archive_prefix"],
        author_name="Codex Cloud D",
        author_email="codex-cloud-d@users.noreply.github.com")
    if backend._remote_head() != cfg["d_branch_head_at_release"]:
        raise core.SegmentRefused("D_BRANCH_HEAD_MOVED_SINCE_C_SIGNED_RELEASE")
    # The calculator is constructed only after every read-only identity/storage
    # gate passes. The preflight path never instantiates it.
    from ase.calculators.eam import EAM
    restored.atoms.calc = EAM(potential=str(potential_path))
    dyn = core.BoundedLangevin(
        restored.atoms, restored.thermostat_rng,
        start_step=EXPECTED_START, stop_at_step=EXPECTED_STOP,
        initial_step=EXPECTED_START,
        release_lock_path=release_path,
        release_lock_sha256=release_sha256,
        execution_code_sha256=_sha(Path(core.__file__)))
    source_provenance = {
        "source_archive_sha256": evidence["archive_sha256"],
        "source_manifest_sha256": evidence["source_manifest_sha256"],
        "source_step": EXPECTED_START,
        "source_d_branch_head": cfg["d_branch_head_at_release"],
        "production_adapter_sha256": _sha(Path(__file__)),
        "archive_backend_sha256": _sha(Path(archive_backend.__file__)),
        "segment_core_sha256": _sha(Path(core.__file__)),
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "archive_backend_test_sha256": EXPECTED_D_ARCHIVE_TEST_SHA256,
        "execution_release_config_sha256": release_sha256,
        "synthetic": False,
    }
    receipts: list[dict] = []

    def remote_roundtrip(archive_path: Path, archive_sha: str) -> None:
        receipt = backend.upload(archive_path, archive_sha)
        checked = backend.verify(receipt)
        restore_path = output_root / "remote_restore" / f"step-{receipt['global_step']:09d}-{archive_sha[:16]}"
        restored_receipt = backend.restore(receipt, restore_path)
        bundle = core.restore_segment_bundle(restore_path)
        if bundle.global_step != receipt["global_step"]:
            raise archive_backend.ArchiveBackendError("RESTORED_STEP_MISMATCH")
        receipts.append({**receipt, "verify": checked, "restore": restored_receipt})

    output_root.mkdir(parents=True, exist_ok=False)
    outputs = core.advance_segment_with_outputs(
        dyn, rng_state_all=restored.rng_state_all,
        source_provenance=source_provenance, output_root=output_root,
        remote_roundtrip=remote_roundtrip, checkpoint_every=100,
        max_wall_s=600.0)
    result = {"status": "SEGMENT_COMPLETE_REMOTE_ARCHIVE_VERIFIED",
              "task_id": EXPECTED_TASK_ID, "start_step": EXPECTED_START,
              "stop_step": dyn.nsteps, "added_steps": dyn.nsteps - EXPECTED_START,
              "production_adapter_sha256": _sha(Path(__file__)),
              "archive_backend_sha256": _sha(Path(archive_backend.__file__)),
              "segment_core_sha256": _sha(Path(core.__file__)),
              "outputs": outputs, "remote_receipts": receipts,
              "final_d_branch_head": backend.expected_head}
    result_path = output_root / "REMOTE_ARCHIVE_RECEIPTS.json"
    temp_path = result_path.with_suffix(".json.tmp")
    temp_path.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    os.replace(temp_path, result_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-config", type=Path, required=True)
    parser.add_argument("--release-config-sha256", required=True)
    parser.add_argument("--adapter-sha256", required=True)
    parser.add_argument("--backend-sha256", required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--stop-step", type=int, required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--potential", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    try:
        if (args.owner, args.task_id, args.start_step, args.stop_step) != (
                EXPECTED_OWNER, EXPECTED_TASK_ID, EXPECTED_START, EXPECTED_STOP):
            raise core.SegmentRefused("CLI_TASK_OR_STEP_BOUNDARY_MISMATCH")
        cfg = validate_execution_release(args.release_config, args.release_config_sha256,
                                         args.adapter_sha256, args.backend_sha256)
        source = core.verify_source_archive(args.source_archive, args.source_manifest)
        if source["archive_sha256"] != cfg["source_archive_sha256"]:
            raise core.SegmentRefused("CLI_SOURCE_ARCHIVE_BINDING_MISMATCH")
        if archive_backend.sha256_file(args.potential) != EXPECTED_POTENTIAL_SHA256:
            raise core.SegmentRefused("CLI_POTENTIAL_SHA_MISMATCH")
        if args.preflight:
            print(json.dumps({"status": "PREFLIGHT_PASS_NO_SCIENCE",
                              "task_id": EXPECTED_TASK_ID,
                              "start_step": EXPECTED_START,
                              "stop_step": EXPECTED_STOP,
                              "release_config_sha256": args.release_config_sha256}))
            return 0
        if not args.execute:
            raise core.SegmentRefused("EXPLICIT_PREFLIGHT_OR_EXECUTE_FLAG_REQUIRED")
        result = run_segment(release_path=args.release_config,
                             release_sha256=args.release_config_sha256,
                             adapter_sha256=args.adapter_sha256,
                             backend_sha256=args.backend_sha256,
                             source_archive=args.source_archive,
                             source_manifest=args.source_manifest,
                             potential_path=args.potential,
                             output_root=args.output_root)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (core.SegmentRefused, archive_backend.ArchiveBackendError,
            OSError, ValueError, KeyError) as exc:
        print(json.dumps({"status": "HOLD", "reason": str(exc)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
