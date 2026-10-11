"""Offline release/checkpoint checks. No calculator or integration step runs."""
from __future__ import annotations

import hashlib
import json
import tempfile
import copy
from pathlib import Path

import continuation_adapter as adapter
import segment_core as core


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def must_refuse(call, expected: str) -> None:
    try:
        call()
    except Exception as exc:
        assert expected in str(exc), (expected, repr(exc))
    else:
        raise AssertionError(f"expected refusal: {expected}")


def main() -> None:
    here = Path(__file__).resolve().parent
    config = here / "SIGNED_CONTINUATION_CONFIG.json"
    config_sha = digest(config)
    sidecar_sha = (here / "SIGNED_CONTINUATION_CONFIG.sha256").read_text().split()[0]
    assert config_sha == sidecar_sha
    cfg = adapter.validate_release_config(config, config_sha)
    assert cfg["start_step"] == 1300 and cfg["stop_at_step"] == 2000
    assert adapter.sha(here / "continuation_adapter.py") == cfg["continuation_adapter_sha256"]
    assert adapter.sha(here / "core/segment_core.py") == cfg["continuation_core_sha256"]
    assert adapter.sha(here / "github_archive_backend.py") == cfg["archive_backend_sha256"]

    with tempfile.TemporaryDirectory(prefix="c-step1300-offline-test-") as tmp_name:
        tmp = Path(tmp_name)
        clone = adapter._clone_pinned_d_branch(tmp, cfg)
        archive = clone / adapter.ARCHIVE_PATH
        bundle_path, evidence, _ = adapter.verify_source(archive, cfg, tmp)
        restored = evidence["restored"]
        assert restored.global_step == 1300 and len(restored.atoms) == 64
        assert restored.state_sha256 == "3b3950729735765dc5c4d205d16a3654b4e9f8ba193c58bf8421685f36d66b9e"
        assert restored.arrays_sha256 == "6eb03c97a33f6c2f3b5d9521c98435a77e5ee9228286090602733c3d22a65a2a"
        assert set(restored.rng_state_all) == {"packing", "velocity", "thermostat"}
        assert restored.thermal_rows and restored.thermal_rows[-1]["step"] == 1300
        assert restored.trajectory_frames and restored.trajectory_frames[-1]["step"] == 1300

        # Constructor validates the signed lock. At the signed stop, the
        # integrator refuses before any calculator or force call can occur.
        positions_before = restored.atoms.get_positions().copy()
        momenta_before = restored.atoms.get_momenta().copy()
        rng_before = copy.deepcopy(restored.thermostat_rng.bit_generator.state)
        dyn = core.BoundedLangevin(
            restored.atoms, restored.thermostat_rng,
            start_step=1000, stop_at_step=2000, initial_step=2000,
            continuation_start_step=1300,
            release_lock_path=config, release_lock_sha256=config_sha,
            execution_code_sha256=adapter.sha(here / "core/segment_core.py"),
            continuation_adapter_path=here / "continuation_adapter.py",
        )
        assert restored.thermostat_rng.bit_generator.state == rng_before
        assert (restored.atoms.get_positions() == positions_before).all()
        assert (restored.atoms.get_momenta() == momenta_before).all()
        must_refuse(dyn.step, "SIGNED_STOP_AT_STEP_REACHED_OR_STEP_OUT_OF_RANGE")
        assert dyn.nsteps == 2000
        assert restored.thermostat_rng.bit_generator.state == rng_before
        assert (restored.atoms.get_positions() == positions_before).all()
        assert (restored.atoms.get_momenta() == momenta_before).all()

        # A changed task identity or signed boundary is rejected even when the
        # mutated file is re-hashed by the caller.
        mutated = dict(cfg)
        mutated["task_id"] = "wrong-task"
        bad = tmp / "wrong-task.json"
        bad.write_text(json.dumps(mutated, sort_keys=True, indent=2) + "\n")
        must_refuse(lambda: adapter.validate_release_config(bad, digest(bad)),
                    "RELEASE_PIN_MISMATCH:task_id")
        mutated = dict(cfg)
        mutated["start_step"] = 1000
        bad = tmp / "wrong-start.json"
        bad.write_text(json.dumps(mutated, sort_keys=True, indent=2) + "\n")
        must_refuse(lambda: adapter.validate_release_config(bad, digest(bad)),
                    "RELEASE_PIN_MISMATCH:start_step")
        must_refuse(lambda: adapter.validate_release_config(config, "0" * 64),
                    "SIGNED_CONTINUATION_CONFIG_SHA_MISMATCH")
        assert not (tmp / "science-output").exists()
        result = {
            "status": "PASS_NO_SCIENTIFIC_INTEGRATION",
            "source_step": restored.global_step,
            "atom_count": len(restored.atoms),
            "rng_streams": sorted(restored.rng_state_all),
            "thermal_rows": len(restored.thermal_rows or []),
            "trajectory_frames": len(restored.trajectory_frames or []),
            "checkpoint_state_sha256": restored.state_sha256,
            "checkpoint_arrays_sha256": restored.arrays_sha256,
            "boundary_2000_refusal": "PASS_BEFORE_CALCULATOR_OR_FORCE_ACCESS",
            "wrong_task_refusal": "PASS",
            "wrong_start_refusal": "PASS",
            "wrong_config_sha_refusal": "PASS",
            "fresh_d_branch_clone_and_receipt_check": "PASS",
            "md_steps_executed": 0,
            "dft_or_training_executed": False,
            "ase_langevin_constructor_rng_draws": 0,
            "checkpoint_state_unchanged_before_any_step": True,
        }
        (here / "OFFLINE_RELEASE_TEST_RESULTS.json").write_text(
            json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
