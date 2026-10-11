"""Non-production validator prototype for a signed 1000->2000 segment proposal.

This module does no MD. It is not wired into run_sampling_codex.py and must not
be treated as an execution release.
"""


class SegmentHold(ValueError):
    pass


def next_step_guard(completed_step, signed_end_step):
    """Pure loop guard prototype: refuses the first attempted step past end."""
    if signed_end_step != 2000:
        raise SegmentHold("ONLY_STEP_2000_ENDPOINT_SUPPORTED")
    if completed_step >= signed_end_step:
        raise SegmentHold("SIGNED_SEGMENT_END_REACHED")
    if completed_step < 1000:
        raise SegmentHold("CHECKPOINT_START_STEP_NOT_REACHED")
    return completed_step + 1


def validate_segment(config, state, *, owner, task_id, requested_stop_step,
                     checkpoint_state_sha256, checkpoint_archive_sha256,
                     actual_config_sha256, externally_pinned_config_sha256):
    if config.get("release_status") != "SIGNED_BY_C":
        raise SegmentHold("SEGMENT_CONFIG_NOT_SIGNED_BY_C")
    if not externally_pinned_config_sha256 or actual_config_sha256 != externally_pinned_config_sha256:
        raise SegmentHold("SIGNED_CONFIG_SHA_MISMATCH")
    if owner != "D" or owner != config.get("owner") or task_id != config.get("task_id"):
        raise SegmentHold("OWNER_OR_TASK_MISMATCH")
    start = config.get("checkpoint_step_required")
    end = config.get("authorized_end_step_if_later_signed")
    if state.get("task_id") != task_id or state.get("step") != start:
        raise SegmentHold("CHECKPOINT_TASK_OR_START_STEP_MISMATCH")
    if checkpoint_state_sha256 != config.get("checkpoint_state_sha256"):
        raise SegmentHold("CHECKPOINT_STATE_SHA_MISMATCH")
    if checkpoint_archive_sha256 != config.get("checkpoint_archive_sha256"):
        raise SegmentHold("CHECKPOINT_ARCHIVE_SHA_MISMATCH")
    if requested_stop_step is None or requested_stop_step != end:
        raise SegmentHold("EXPLICIT_STOP_STEP_MISSING_OR_MISMATCH")
    if not isinstance(end, int) or not isinstance(start, int) or (start, end) != (1000, 2000):
        raise SegmentHold("SEGMENT_LENGTH_NOT_EXACTLY_1000_STEPS")
    if config.get("maximum_added_steps_if_later_signed") != end - start:
        raise SegmentHold("SEGMENT_EXCEEDS_SIGNED_STEP_CAP")
    if state.get("scientific_parameters_sha256") != config.get("scientific_parameters_sha256_must_remain"):
        raise SegmentHold("SCIENTIFIC_PARAMETER_IDENTITY_CHANGED")
    if state.get("source_manifest_sha256") != config.get("source_manifest_sha256_must_remain"):
        raise SegmentHold("SOURCE_MANIFEST_IDENTITY_CHANGED")
    rng = state.get("rng_state")
    if not isinstance(rng, dict) or set(rng) != {"packing", "velocity", "thermostat"}:
        raise SegmentHold("RNG_STATE_MISSING_OR_INCOMPLETE")
    if any(rng[name].get("bit_generator") != "PCG64" for name in rng):
        raise SegmentHold("RNG_TYPE_CHANGED")
    if config.get("stage_continue_default_endpoint_warning") != 120000:
        raise SegmentHold("KNOWN_DEFAULT_CONTINUE_ENDPOINT_NOT_BOUND")
    return {"status": "VALIDATOR_PASS_ONLY_NOT_COMPUTE_AUTHORIZATION",
            "start_step": start, "end_step": end, "added_steps": end - start}
