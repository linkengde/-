# Ag64 step2000 extension review fix report

## Result

```text
ENVIRONMENT_EVIDENCE_PASS=YES
BACKEND_V2_REMOTE_RESTORE_PASS=YES
TIMEOUT_RECOVERY_PASS=YES
ASE_3_29_0_RUNTIME_AVAILABLE=NO
```

`ENVIRONMENT_EVIDENCE_PASS=YES` means the evidence records the measured runtime and SHA-256 values of the source files at the tested commit, and the evidence-file digest verifies. It does not mean the runtime is execution-ready: ASE 3.29.0 is absent, so adapter/core imports in the real Python environment are blocked. No ASE installation was attempted, consistent with the earlier instruction not to install it.

## Commits and evidence binding

- Branch: `ag64-step2000-extension-dev`
- Tested source commit: `148c9b1f29f43eac6991363cb4895797b597bebc`
- Runtime evidence: `ENVIRONMENT_EVIDENCE.json`
- SHA-256 of the exact evidence JSON bytes: `3131a14b92ca9e8ec314b9d1e8be92af0fee9875dd7674e7845d6ac2049c82de`
- Evidence hash file: `ENVIRONMENT_EVIDENCE_SHA256.txt`

The evidence binds these exact committed file contents:

| File | SHA-256 |
|---|---|
| `continuation_adapter.py` | `4dfc0485f90d67cb56ab9ee8c8721ee2bd749aa13027d15d6c0f99c7d2a92fbb` |
| `core/segment_core.py` | `16610a77a5ff486668f3be81a1747d5180a09f2e488b4ec8be43a580ed7f6972` |
| `github_archive_backend.py` | `0a8d8b71ee938e768f47bc9e819c5953f7f7bb0aa6415f248916b4a300a3928f` |

The adapter and archive-backend source files were not changed. The production potential was not opened or loaded. Its pinned SHA is recorded from the core source constant only.

## Changes

- `core/segment_core.py`: added scoped SIGINT/SIGTERM capture, graceful-stop polling, and checkpoint-on-stop handling for walltime expiry, `KeyboardInterrupt`, `TimeoutError`, and other step-loop exceptions. A stop saves the latest valid state in the existing checkpoint bundle format and raises `SegmentCheckpointSaved`; it never auto-resumes or crosses the configured end step. Emergency step-2500 checkpoints are local-only, so the established archive cadence remains every 1000 steps and at the configured end.
- `tests/test_step2000_preproduction.py`: added virtual timeout, signal, interrupt, graceful-stop, restore, and resume checks. The virtual runner checks state and RNG against an uninterrupted virtual control.
- `tests/test_backend_v2_local_remote.py`: exercises the existing production backend methods against a temporary local bare Git remote, including the historical unbound-step guard.
- Evidence/report files were regenerated under this extension directory. No sampling settings, adapter, backend, potential, or step1900→2000 package was changed.

## Environment evidence

- Python: `3.12.14`
- NumPy: `2.3.5`
- ASE required: `3.29.0`; found: **missing**
- Searched runtime and project environment locations, including Python package paths, venv/conda markers, and pip cache; no ASE installation or cached ASE wheel was found.
- Actual import checks:

```text
CONTINUATION_ADAPTER_IMPORT=BLOCKED ModuleNotFoundError: No module named 'ase'
SEGMENT_CORE_IMPORT=BLOCKED ModuleNotFoundError: No module named 'ase'
GITHUB_ARCHIVE_BACKEND_IMPORT=PASS
```

Offline tests import adapter/core using the test-only ASE stubs. They do not provide evidence that the installed ASE runtime is compatible. That remains a runtime readiness blocker.

## Backend v2 remote restore

The actual `GitHubArchiveBackend.upload()`, `verify()`, and `restore()` methods passed against a temporary local bare Git remote. This is a synthetic transport fixture, not a scientific step-3000 result; no mock backend was used. A fresh clone was made after upload, then the backend's own independent fresh-clone verification and restoration ran.

```text
UPLOAD_STATUS=PASS
FRESH_CLONE_STATUS=PASS
VERIFY_STATUS=PASS
RESTORE_STATUS=PASS
BACKEND_V2_REMOTE_RESTORE_PASS=YES
historical unbound step-3000 upload=BLOCKED (PASS)
```

Verified digests and bindings:

| Item | SHA-256 / value |
|---|---|
| Transport fixture archive | `3a2724344302d2201231946739b42db6602a616fb93e057e07f89deb3e704e8d` |
| Manifest | `1c7d933a3c48a7912927f6cc91784a591bbdda54bec0a0654cf876e2bb27e9e2` |
| Receipt | `b9dd63096b4a6c87748c555993036166eb6dcaea71888793e49a8998109b63a4` |
| Restored checkpoint `state.json` | `3efa8a6728b48196498e8fafe3bef1d1a5c4dffe45a446a18b0c313a5ed7ba07` |
| Restored arrays | `5d625bc34b42ac8a339efc302a31caafc78ff20111ca4cc1e366974c1116f2e4` |
| Restored RNG state | `9f31f80534f18bb25f441c3ca1474ba1aafaae5130fd584311e82aa78831cb91` |
| Receipt segment | `2000 → 3000`, `max_added_steps=1000` |
| Source checkpoint SHA bound in receipt | `5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a` |

Arrays and RNG state matched exactly after restore. The legacy upload path still rejected an unbound step-3000 archive, confirming its historical step restriction remains active. The detailed run output is in `BACKEND_V2_REMOTE_RESTORE_TEST.json`.

Full backend test output:

```json
{
  "RESTORE_STATUS": "PASS",
  "UPLOAD_STATUS": "PASS",
  "VERIFY_STATUS": "PASS",
  "archive_sha256": "3a2724344302d2201231946739b42db6602a616fb93e057e07f89deb3e704e8d",
  "arrays_equal_after_restore": true,
  "arrays_sha256": "5d625bc34b42ac8a339efc302a31caafc78ff20111ca4cc1e366974c1116f2e4",
  "backend_v2_remote_restore_pass": true,
  "base_head": "1578305a67c7f9468644cb9ef505f07424b6140b",
  "checkpoint_state_sha256": "3efa8a6728b48196498e8fafe3bef1d1a5c4dffe45a446a18b0c313a5ed7ba07",
  "fresh_clone_head": "b67fc1f041bb071f011696f36df9b6d04d21a6f9",
  "fresh_clone_status": "PASS",
  "historical_step_guard": "PASS_BLOCKED_UNBOUND_STEP_3000",
  "manifest_sha256": "1c7d933a3c48a7912927f6cc91784a591bbdda54bec0a0654cf876e2bb27e9e2",
  "manifest_step": 3000,
  "md_run": false,
  "potential_loaded": false,
  "production_backend_modified": false,
  "real_backend_methods": [
    "upload",
    "verify",
    "restore"
  ],
  "receipt_segment": {
    "checkpoint_sha256": "3efa8a6728b48196498e8fafe3bef1d1a5c4dffe45a446a18b0c313a5ed7ba07",
    "end_step": 3000,
    "max_added_steps": 1000,
    "source_checkpoint_sha256": "5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a",
    "start_step": 2000
  },
  "receipt_sha256": "b9dd63096b4a6c87748c555993036166eb6dcaea71888793e49a8998109b63a4",
  "restore_result": {
    "archive_sha256": "3a2724344302d2201231946739b42db6602a616fb93e057e07f89deb3e704e8d",
    "global_step": 3000,
    "member_count": 6,
    "owner": "D",
    "restore_path": "<temporary-local-test-directory>",
    "status": "REMOTE_RESTORE_VERIFIED",
    "task_id": "C-D-AGL-1250K-rho097-seed01"
  },
  "rng_equal_after_restore": true,
  "rng_state_sha256": "9f31f80534f18bb25f441c3ca1474ba1aafaae5130fd584311e82aa78831cb91",
  "source_checkpoint_sha256": "5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a",
  "synthetic_transport_fixture_only": true,
  "test_remote_type": "local_bare_git_file_remote",
  "upload_commit_sha": "b67fc1f041bb071f011696f36df9b6d04d21a6f9",
  "verify_result": {
    "archive_sha256": "3a2724344302d2201231946739b42db6602a616fb93e057e07f89deb3e704e8d",
    "checkpoint_sha256": "3efa8a6728b48196498e8fafe3bef1d1a5c4dffe45a446a18b0c313a5ed7ba07",
    "end_step": 3000,
    "fresh_head": "b67fc1f041bb071f011696f36df9b6d04d21a6f9",
    "member_count": 6,
    "start_step": 2000,
    "status": "REMOTE_FRESH_FETCH_VERIFIED",
    "upload_commit": "b67fc1f041bb071f011696f36df9b6d04d21a6f9"
  }
}
```

## Timeout recovery test

The runner was simulated with deterministic array-only virtual dynamics, not ASE. It stopped at step 2500 when the virtual walltime expired, created a local `checkpoint_step2500` bundle without creating an out-of-cadence archive, restored it, and continued from 2501 through 3000. Restored positions, velocities, and RNG matched the saved state; the resumed final positions, velocities, and RNG matched an uninterrupted virtual control. The visited sequence was exactly 2001–3000, with no duplicate and no step 3001. The step-3000 archive callback passed.

SIGTERM, `KeyboardInterrupt`, and the graceful-stop callback also produced recovery checkpoints in the virtual tests.

Full unit-test log:

```text
test_candidate_config_parses_and_2000_to_3000_boundary_passes ... ok
test_restored_global_step_is_2000 ... ok
test_rng_streams_restore_without_state_change ... ok
test_signal_interrupt_and_graceful_stop_each_save_current_checkpoint ... ok
test_step_3001_is_refused_without_integrator_call ... ok
test_walltime_checkpoint_restore_resume_matches_continuous_virtual_run ... ok
test_wrong_checkpoint_sha_is_rejected ... ok
test_wrong_potential_sha_is_rejected_before_source_access ... ok

Ran 8 tests in 0.123s
OK
```

Commands:

```text
python -m py_compile ag64_step2000_extension_dev/core/segment_core.py ag64_step2000_extension_dev/tests/test_step2000_preproduction.py ag64_step2000_extension_dev/tests/test_backend_v2_local_remote.py
AG64_STEP2000_ARCHIVE=/tmp/ag64-step2000.tar.gz python -m unittest discover -s ag64_step2000_extension_dev/tests -v
python ag64_step2000_extension_dev/tests/test_backend_v2_local_remote.py
```

The backend command output, including upload/verify/restore commit, fresh-clone head, receipt and member hashes, is saved in `BACKEND_V2_REMOTE_RESTORE_TEST.json`. Timeout details are saved in `TIMEOUT_RECOVERY_TEST.json`.

## Remaining blocker and execution state

- The evidence and offline transport/recovery checks pass.
- ASE 3.29.0 remains unavailable in this runtime. Real adapter/core import and any real ASE dry-run remain blocked until an approved existing ASE 3.29.0 environment is provided.
- No real MD was started; no production potential or ASE integrator was used.
- No `execution_authorized=true` configuration was created.
- The existing unsigned candidate configuration was not edited; because the core SHA changed, it must be reviewed and rebound before any future execution request.
