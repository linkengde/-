# V19 CIF matched-kmesh pilot entry review

## Decision

**BLOCKED for execution pending A's review and job registration.** The package is suitable in scope for a fixed-ion numerical mesh-sensitivity probe, once the entry and archive safeguards below are corrected. It cannot support an equilibrium-surface, MACE-validation, or dataset-integration claim.

No DFT, model inference/training, MD, or TTM was started. `launch_enabled` remains `false`; both owner jobs remain absent from their task-card label maps. A must review the changes, register the two labels with the respective owner cards, and explicitly enable/authorize the jobs before any launch.

## Source, input and owner checks

All seven files listed in the pilot package's `SHA256SUMS.txt` match. The source manifest hash is `f543091f67ed1a85bf005869d95d5a9931ebf524defd0c3e3e81d93dde7f2342`. The pilot geometry is byte-identical to the CIF-derived AgSi proposal (`8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`). Both mesh records use that same input hash; Gamma is assigned to A and 2×2×1 to B. The runtime identity is B's registered instance, `c5035b48-f43b-4da0-b8e4-e2862817f86a`, and both manifest owner IDs match their corresponding task-card owners.

The separate CIF review found the uploaded source internally consistent with its COD/AMCSD metadata and reproduced the slab files. It did not verify the remote COD checksum. The pilot geometry is an asymmetric, one-sided Ag-covered, one-conventional-cell-thick slab with about 16 Å vacuum. It is intentionally a frozen stress geometry; these properties make a matched numerical comparison useful while ruling out an equilibrium or general validation interpretation.

The setup is controlled as intended: identical fixed geometry, PBE plane waves at 500 eV, Fermi smearing 0.1 eV, and four MPI ranks. The calculator takes k-points from each manifest record, and the verifier checks the recorded k-points. Native extrapolated energy and free energy are stored separately. Settings are spread across README, manifest and code, however: the manifest does not record XC, energy convention, GPAW version, mixer and SCF thresholds. Put the full method record in the manifest and have the verifier compare it with the archived summary before registration.

## Entry and runner findings

| Area | Result | Evidence and required correction |
|---|---|---|
| Launch disabled | PASS | The calculator checks `launch_enabled` before importing `gpaw_data` or GPAW. A prior negative guard probe for the B label exited with “Pilot not approved” and created no output directory. |
| Queue disabled guard | FAIL | `run_queue.py` calls `sync.claim(TASK)` before reading/asserting `launch_enabled`. Calling the disabled queue can mutate task state before refusing. Move all manifest/owner/launch checks ahead of the claim, or make a disabled claim a read-only refusal. |
| Owner registration | BLOCKED | Neither pilot label appears in either owner's current `labels` map. Register `AgSi_COD9009647_pilot_gamma` only to A and `AgSi_COD9009647_pilot_k2x2` only to B, with exact instance IDs, before enabling. |
| Direct calculator ownership | FAIL | `run_pw_reference.py` validates the input hash and global disabled flag but does not enforce the record's `owner_task`/`owner_instance`. Add an owner check against the current sync identity and registered task before GPAW import; retain the queue check as a second layer. |
| Existing-run handling | PARTIAL | Queue refuses an unarchived output directory only when no archive exists. If an archive and stale output coexist, it verifies the archive and leaves the stale run unchecked. Stop for same-job review whenever both paths exist. |
| Archive hashes | FAIL | `verify_result.py` reports hashes but does not validate an existing `SHA256SUMS.txt`. The queue rewrites `verification.json` and the inventory for an existing archive, which can mask altered files and makes the verification record refer to the previous inventory. For an existing archive, first validate its complete inventory and refuse any mismatch; never regenerate hashes over an existing archive without explicit review. |
| Result verification | PARTIAL | It checks input hash, symbols, exact positions/cell/PBC/IDs/pair markers, convergence text, four MPI ranks, finite labels, energy and key method fields. Add checks for method/version/energy-convention fields, summary formula/count, progress iteration, and force-derived `fmax`/RMS consistency. |
| Checkpoint path | FAIL | `state_path` resolves to `OUT/state.gpw`, and `OUT` is under `/workspace/mace_v19_cif_pilot_<owner>/...`; the comment says it is kept in `/tmp`. The checkpoint is not archived, but stays in that workspace run directory. Correct the comment and document retention, or use an explicit approved local checkpoint path. Do not delete checkpoints automatically. |
| Resource guard | FAIL | Queue checks only that free space exceeds 500 MiB before a run; it does not monitor disk during SCF/checkpoint writes. For the 1107.7 Å³ cell and 500 eV cutoff, the plane-wave sphere estimate is about 28,122 coefficients per k-point. Raw complex128 coefficients for four k-points are about 203 MiB at 118 bands or 446 MiB at 260 bands, before densities, GPAW metadata, file-write overhead, and safety margin. This is an estimate, not a measured `.gpw` size. The task's A-space snapshot was 994 MiB; I cannot refresh A's separate environment from B. B's current workspace has 6.24 GB free. Use a measured comparable checkpoint or conservative size estimate, require adequate headroom, and check free space during the run before each checkpoint. |

The persistent `/workspace/.mace-v19-dft.lock` file is zero bytes, has no matching active lock entry, and no matching GPAW/MPI process or B run/archive directory was present at audit time. The file was left untouched. This is a status snapshot, not permission to launch.

## Comparison logic and interpretation

`compare_mesh_results.py` first invokes result verification for both archives, requires identical positions/cell/species/IDs, and reports native and free-energy change per atom, atom-vector force-difference RMSE, maximum atom force-difference norm, and the signed marked-pair projection. The pair axis is the minimum-image Ag-to-Si vector; the projected quantity is the Si-minus-Ag force difference between meshes. These definitions are suitable for a numerical sensitivity diagnostic. The proposed 2 meV/atom, 0.01 eV/Å force-vector RMSE, and 0.02 eV/Å pair-projection budgets are provisional triggers, not proof of convergence.

The Gamma/2×2×1 pair is a sensible minimal first comparison on one frozen geometry. It does not establish k-point convergence by itself; if a budget is exceeded, extend the mesh and settle one common target convention before using labels. Vacuum and any supported dipole treatment are separate checks. None of these results should be used as equilibrium labels or MACE validation data.

## Required changes before A registers jobs

1. Move launch and owner preflight ahead of `sync.claim`; enforce owner identity inside the calculator entry too.
2. Register each exact label on its owner card while keeping launch disabled until A explicitly approves.
3. Make existing archive handling immutable: validate the existing complete SHA inventory, reject mismatches and archive/output coexistence, and never overwrite evidence silently.
4. Add the missing method and force-summary checks to verification; compare only after both complete archives pass them.
5. Correct the checkpoint-path documentation and add a conservative preflight plus runtime disk reserve check.

The report does not edit the production-intended pilot package or grant DFT authorization. The static evidence and current status snapshot are in `audit.json`; `audit_static.py` reproduces file hashes, task/owner checks, static entry findings and the local resource snapshot.
