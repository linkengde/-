"""Read-only checkpoint restore regression for the step-48000 continuation entry."""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE))
import followup_runner_step48000_extension as runner


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    config_path = PACKAGE / "FOLLOWUP_EXECUTION_CONFIG_STEP48000_60000.json"
    config = json.loads(config_path.read_text())
    repo = runner._git_root(PACKAGE)
    restored = runner._load_resume_checkpoint(repo, config)
    checkpoint = restored["checkpoint"]
    atoms = restored["atoms"]
    if restored["checkpoint_sha256"] != "f21165280a3524612ab51702f7b28ee56643aace598d484f2966d0ce25291b57":
        raise AssertionError("wrong checkpoint digest")
    if checkpoint["global_step"] != 48000 or restored["state"]["global_step"] != 48000:
        raise AssertionError("wrong restored step")
    if restored["state"]["stage"] != "A_HOT_1600K":
        raise AssertionError("wrong restored stage")
    if restored["qa_monitor"].get("gate_status") != "PENDING":
        raise AssertionError("QA is not PENDING at the archived boundary")
    if list(atoms.get_chemical_symbols()) != ["Ag"] * 64:
        raise AssertionError("wrong composition")
    if atoms.get_positions().tolist() != checkpoint["positions_A"]:
        raise AssertionError("positions did not restore exactly")
    if (atoms.get_momenta() * runner.fs).tolist() != checkpoint["momenta_amu_A_per_fs"]:
        raise AssertionError("momenta did not restore exactly")
    for name in runner.RNG_NAMES:
        if restored["rngs"][name].bit_generator.state != checkpoint["rng_states"][name]:
            raise AssertionError(f"RNG stream {name} did not restore exactly")
    if config["step_policy"]["maximum_global_step"] != 60000:
        raise AssertionError("global step safety cap changed")
    if config["step_policy"]["maximum_added_steps"] != 12000:
        raise AssertionError("continuation step budget is wrong")
    output = Path(config["output_directory"])
    if output.exists():
        raise AssertionError("output path must be new and isolated before preflight")
    # The immutable archived copy and separately fresh-fetched copy must agree.
    fresh = Path("/tmp/d-step48000-fresh-fetch") / config["resume_checkpoint"]["archive_path"]
    if sha(fresh / "checkpoint.json") != sha(restored["archive_dir"] / "checkpoint.json"):
        raise AssertionError("fresh-fetch checkpoint differs")
    print("CHECKPOINT_RESTORE_PASS=YES")
    print("POSITIONS_MOMENTA_RNG_PASS=YES")
    print("STAGE_STEP_QA_BOUNDARY_PASS=YES")
    print("FRESH_FETCH_CHECKPOINT_MATCH_PASS=YES")
    print("OUTPUT_ISOLATION_PASS=YES")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
