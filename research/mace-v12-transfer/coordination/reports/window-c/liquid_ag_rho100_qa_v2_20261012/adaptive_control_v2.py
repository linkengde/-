#!/usr/bin/env python3
"""Pure control gates for the bounded Ag64 V2 continuation.

No calculator or integrator is imported here.  The production runner must call
these functions before opening the next chunk; failed archive verification can
never authorize a phase transition or another chunk.
"""
from __future__ import annotations


class ControlError(RuntimeError):
    pass


GLOBAL_STEP_CAP = 60000
CHUNK_STEP_CAP = 1000
RAMP_STEPS = 10000
CANDIDATE_STEPS = 10000
CHUNK_WALL_SECONDS = 3600


def next_chunk_end(current_step: int, phase_end: int, global_cap: int = GLOBAL_STEP_CAP) -> int:
    if current_step < 1000 or phase_end < current_step or global_cap < current_step:
        raise ControlError("invalid continuation boundary")
    end = min(current_step + CHUNK_STEP_CAP, phase_end, global_cap)
    if end <= current_step:
        raise ControlError("no authorized step remains")
    return end


def can_start_next_step(current_step: int, segment_end: int, elapsed_seconds: float,
                        wall_cap_seconds: int = CHUNK_WALL_SECONDS) -> bool:
    """Check limits before beginning a step; a started step may finish."""
    if elapsed_seconds < 0 or wall_cap_seconds <= 0:
        raise ControlError("invalid wall-clock value")
    return current_step < segment_end and elapsed_seconds < wall_cap_seconds


def archive_gate(result: dict | None) -> bool:
    required = ("archive_upload_pass", "archive_fresh_fetch_pass", "archive_restore_pass",
                "checkpoint_restore_pass", "sha256_match_pass")
    return isinstance(result, dict) and all(result.get(name) is True for name in required)


def advance_after_archive(state: dict, current_step: int, gate_status: str,
                          archive_result: dict | None, candidate_spacing: dict | None = None) -> str:
    """Apply a state change only after a verified remote checkpoint restore."""
    if not archive_gate(archive_result):
        state["stage"] = "HOLD_ARCHIVE_OR_RESTORE_FAILURE"
        return state["stage"]
    if current_step > GLOBAL_STEP_CAP:
        state["stage"] = "HOLD_GLOBAL_STEP_LIMIT"
        return state["stage"]
    stage = state.get("stage")
    if stage == "B_COOL_1600_TO_1250K" and current_step >= int(state["stage_b_start_step"]) + RAMP_STEPS:
        state["stage"] = "C_HOLD_1250K"
        state["stage_c_start_step"] = current_step
        return "CHUNK_COMPLETE"
    if stage == "D_CANDIDATE_1250K" and current_step >= int(state["stage_d_start_step"]) + CANDIDATE_STEPS:
        state["stage"] = "CANDIDATE_WINDOW_COMPLETE_FINAL_QA_REQUIRED"
        return state["stage"]
    if stage == "A_HOT_1600K" and current_step >= 20000 and gate_status != "PASS":
        state["stage"] = "HOLD_STAGE_A_CAP_WITHOUT_LIQUID_GATE"
        return state["stage"]
    if stage == "C_HOLD_1250K" and current_step >= min(60000, int(state["stage_c_start_step"]) + 20000) and gate_status != "PASS":
        state["stage"] = "HOLD_STAGE_C_CAP_WITHOUT_LIQUID_GATE"
        return state["stage"]
    if current_step >= GLOBAL_STEP_CAP:
        state["stage"] = "HOLD_GLOBAL_STEP_LIMIT"
        return state["stage"]
    if gate_status == "FAIL_WINDOW":
        state["stage"] = "HOLD_QA_FAIL_WINDOW"
        return state["stage"]
    if gate_status == "QA_HOLD":
        state["stage"] = "HOLD_QA"
        return state["stage"]
    if stage == "A_HOT_1600K" and gate_status == "PASS":
        state["liquid_onset_step"] = current_step
        state["stage_b_start_step"] = current_step
        state["stage"] = "B_COOL_1600_TO_1250K"
    elif stage == "C_HOLD_1250K" and gate_status == "PASS":
        if not isinstance(candidate_spacing, dict) or int(candidate_spacing.get("snapshot_spacing_steps", 0)) < 1000:
            state["stage"] = "HOLD_CANDIDATE_SPACING_UNVERIFIED"
            return state["stage"]
        state["candidate_snapshot_spacing_steps"] = int(candidate_spacing["snapshot_spacing_steps"])
        state["stage_d_start_step"] = current_step
        state["stage"] = "D_CANDIDATE_1250K"
    elif gate_status not in {"PASS", "PENDING", "FAIL_WINDOW"}:
        state["stage"] = "HOLD_UNKNOWN_QA_STATUS"
    if state.get("stage") == "CHUNK_COMPLETE":
        raise ControlError("CHUNK_COMPLETE is not a valid phase name")
    return "CHUNK_COMPLETE" if not str(state.get("stage", "")).startswith("HOLD_") else state["stage"]
