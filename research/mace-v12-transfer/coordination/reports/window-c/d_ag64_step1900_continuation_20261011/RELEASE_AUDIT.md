# Release audit

{
  "task_id": "C-D-AGL-1250K-rho097-seed01",
  "start_step": 1900,
  "stop_at_step": 2000,
  "max_added_steps": 100,
  "d_branch_head_at_release": "3ae9aa1038235d38d19b02218081e8a3fdaaaf28",
  "source_checkpoint_state_sha256": "3b9918e46775a8b3e48d4ee4135118651240d9779d86813d861d66ea4abcc3f9",
  "source_checkpoint_arrays_sha256": "90d322c51442f0c94f8ba8cf513afbe1e6d7b4cb8505b501030cdb2d645de91d",
  "source_archive_sha256": "9081c7e5e55f0c05a8372b3d4795b0e0b8007346fd401c5de45baec5502510ec",
  "potential_sha256": "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447",
  "continuation_adapter_sha256": "0ddfc7db5d85b863ff7e16b25d0e4369350f44faf0c021206cd09a26eab9f71c",
  "continuation_core_sha256": "3033c3234c4074358f8ef54a6783c3c4c9495549c12c93641f4f777f3451ed76",
  "archive_backend_sha256": "8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346"
}

Source state identifies the prior step1600 signed config and executor; its original state is preserved. Restore validation compares positions, momenta, velocities and PCG64 state exactly without constructing a calculator. Source archive and all inner members are checked from a fresh remote clone. The inherited scientific_parameters_sha256 remains a source provenance identity; new adapter/core hashes are recorded separately.

Signing is procedural C authorization bound by Git commit and configuration SHA, not a cryptographic identity signature. Only this task and this 100-step interval are authorized. Other five tasks, DFT and training remain HOLD. No execution was dispatched.
