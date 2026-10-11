# F Ag64 warmup archive and recovery plan

## Scope and prerequisites

This plan targets only F's two assigned warmup tasks (`C-F-AGL-1250K-rho100-seed02` and `C-F-AGL-1250K-rho103-seed01`). The current Codex Cloud workspace is task-local temporary storage. It is suitable for live execution scratch, but it is not the cross-task recovery layer. The current release config still has `storage_release=false` and `warmup_release=false`; no sampler root, archive, structure, or trajectory is created by this plan.

C must publish and sign a Codex-specific release config and SHA that binds the F task roots, permits a temporary execution root with GitHub branch archival, and defines the archive/restore procedure. Until then, both tasks stay HOLD. A/B paths and checkpoints remain outside this workspace and untouched.

## Proposed staging and archive

After a signed release only:

1. Use a task-scoped sampling root under `/workspace`, outside `/workspace/-`, separate from the DFT root, and precreate it only after the signed config names its exact path. Run each assigned task serially under its own task ID.
2. At the completed 1000-step warmup boundary, package only that task's restart bundle and the small identity files needed to verify it: `sampler_state.json`, the selected `restart/step-*` bundle (`state.json`, `state_arrays.npz`, `state_atoms.traj`, `trajectory.traj`, `thermal.csv`), `source_manifest.json`, `stage_wall_budget.json`, and integrity log. Include the C config SHA, task ID, owner, sampler/source SHA identities, and a SHA-256 manifest in the archive envelope.
3. Do not include the EAM parameter file, software environments, A/B data, other owners' files, or unrelated checkout files. Do not upload an archive unless the repository's visibility and C's data-handling rules permit the trajectory/state bytes.
4. Enforce the sampler's 512 MiB per-task archive cap and 32 MiB checkpoint reserve. For normal Git blobs, stop and ask C for a different archive transport if any individual file reaches 50 MiB or the archive exceeds 50 MiB; this leaves margin below GitHub's 100 MB per-file limit. No LFS or other large-file service is configured or tested.
5. Commit the archive and manifest only to `cloud-f-b1-dft`, using a normal fast-forward push. Record the commit SHA and archive SHA in C's report. Never force-push or edit `main`.

## Restore procedure

In a later isolated task, fetch the F branch into a new temporary staging directory, verify the archive SHA-256 against the C-signed manifest, verify every inner-file SHA and owner/task/config identity, then extract to a new, empty, C-approved task root outside the checkout. Check path traversal and symlinks before extraction, confirm the expected checkpoint step/status, and retain the fetched commit SHA with the restored state. Do not merge restored data into an existing root or overwrite any old checkpoint. The sampler's own restart checks and C's task-specific resume release remain required before any continuation.

The recovery test in this report uses only the deterministic synthetic `RECOVERY_TEST.bin`; it proves Git transport, hashing, and byte restoration, not scientific archive behavior or future task persistence.
