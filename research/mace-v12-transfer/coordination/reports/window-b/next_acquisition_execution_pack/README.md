# Future acquisition execution pack — REVIEW ONLY, launch disabled

The six proposals are unchanged copies of the previous batch recommendation, not new candidates or blind tests. Every planned calculation is `blocked_pending_A_review_and_V15_screen`. There is no authorization flag, automatic six-job launcher or owner claim. A must publish explicit review, production integration and owner assignments before any calculation. B did not run/import GPAW, perform training/inference, DFT, MD or TTM, or inspect V15 blind calculations/summaries/labels.

## Proposed balanced ownership

| Proposed owner | Planned label | Candidate |
|---|---|---|
| B | AgC_widespan_distance_acq_next_01 | AgC_wide_distance_p0p35 |
| B | AgSi_widespan_distance_acq_next_01 | AgSi_wide_distance_p0p50 |
| B | AgTi_widespan_distance_acq_next_01 | AgTi_wide_distance_p0p35 |
| A | AgC_widespan_registry_acq_next_01 | AgC_wide_registry_p0p45 |
| A | AgSi_widespan_registry_acq_next_01 | AgSi_wide_registry_m0p30 |
| A | AgTi_widespan_registry_acq_next_01 | AgTi_wide_registry_p0p45 |

Names do not appear in any previously tracked project path. `owner_assigned` is false in every manifest record. These are proposed training acquisitions; do not change the existing blind split or treat these as blind validation.

## Rescreen and exact identity

Six input files preserve source bytes/SHA256, IDs/order, cell/PBC and energy/force label absence. Screening used 36 committed train/valid/test and acquisition/frozen-blind INPUT files containing 178 frames, including the integrated 29-frame V15 dataset. Same-composition ID mappings are resolved; no exact/near/other-PBC coordinate overlap is found. No exact/near pair among the six is found. Source inventory and per-candidate results are in screening_report.json. Original species minima and parent hashes are in input_manifest.json. Geometry PASS is triage, not chemical stability or transferable-model validation.

Policy retained from wide-span design: align framework translation using matching persistent IDs/elements; same cell within1e-8 Å. Exact max displacement1e-6 Å; near Ag and framework RMS≤0.15 Å, framework max≤0.25 Å. PBC identity is recorded separately. No rotated/permuted equivalence is claimed. Shared training-parent motif remains a correlation. A must rerun the screen if committed data or planned acquisitions change after this snapshot, without inspecting blind labels for selection.

## Draft behavior and future A integration

`run_pw_reference_draft.py` has an unconditional stop **before imports of GPAW or any output**. No argument or environment variable enables it. It is source for A review, not an executable production entry. Static compilation was performed; no calculator imports or runner execution occurred.

The disabled source describes explicit `--repo-root`, `--label`, `--output-dir` arguments, four MPI ranks, PW-PBE500/Gamma/FermiDirac0.1, original mixer/convergence settings, and a new `/tmp` work directory. Atomic rank0 output-directory creation refuses any existing run/checkpoint/active marker. It preserves `.gpw` checkpoints locally and failure diagnostics, requires finite native energy/forces and records free_energy separately if supported. Selected target is explicitly `get_potential_energy(force_consistent=False)`; do not substitute free energy. Successful calculation removes its active marker; failed/unconverged runs retain evidence. No automatic result publication or queue launch is included.

`archive_draft.py` is a verification/copy draft with explicit `--repo-root`, `--label`, `--run-dir`, `--archive-dir` arguments. It requires completed four-rank convergence, source linkage/hash, geometry/IDs/marked pair/PBC identity, finite labels and chosen/free-energy agreement, then stages compact files with verified copy hashes. It refuses existing archives/staging and active run markers. Checkpoints are not copied. A must integrate the archive location/role and authorization/ownership rules with the actual task cards; archive verification does not grant permission to start jobs.

Intended argument shapes for **A-integrated future versions only**, not commands to launch this pack:

```text
four-rank MPI runner: --repo-root REPO --label ASSIGNED_LABEL --output-dir NEW_TMP_RUN
archive verifier: --repo-root REPO --label ASSIGNED_LABEL --run-dir COMPLETED_RUN --archive-dir NEW_COMPACT_ARCHIVE
```

## Readiness checklist

- [x] Six copied inputs exactly match source hashes and contain no labels.
- [x] Current committed split/input inventory screened; zero hits/unresolved mappings.
- [x] Existing priority/near-pair review retained; balanced ownership proposed only.
- [x] Input manifest marks all work blocked; runner launch cannot be enabled by a flag.
- [x] Python syntax compiled without execution; no GPAW import/calculation.
- [x] Four ranks, convention, nonempty/active guards, hash/finite/convergence checks drafted.
- [ ] A reviews V15 numerical screen and explicitly assigns owners/production labels.
- [ ] A reviews/integrates runner/archive drafts and validates their runtime behavior.
- [ ] A confirms available independent environments/resources and absence of conflicting active work.

Runtime readiness is UNDETERMINED. Preparing source is not proof a calculation can run or converge. Hash checks:

```bash
cd /workspace/-
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/next_acquisition_execution_pack/SHA256SUMS.txt
```
