# Two fresh cloud windows: resume after completed four-mesh pilots

User requested restart only after both owners finished and uploaded. That condition is satisfied for Gamma/2x2/3x3/4x4; all compact archives and B ladder assessment were independently verified by A. No additional DFT is to start in the old environments. Live processes and virtual environments do not migrate through Git.

Primary read-only results entry: [V19 CIF convergence bundle](reports/window-b/v19_cif_convergence_bundle/README.md), published e4a65a2. A independently verified all44 root checksum entries. The bundle duplicates source archives without removing originals; no launch scripts are included.

## Completed; do not repeat

Six signed residual labels, training36 candidate/rank4 audit, uploaded COD9009647 CIF and source-derived52atom slabs, Gamma51/2x2SCF41/3x3SCF51/4x4SCF41 pilots and their verified compact archives are on main. V18 fixed model/checkpoint/logs/selection/evaluation are on main; V19 model is NOT trained. Do not read sealed reference labels. V17 test geometry exposure is documented in coordination/V17_TEST_GEOMETRY_EXPOSURE.json; these inputs are not strictly untouched geometry evidence.

Gamma/2x2 force difference1.332792eV/A;3x3/4x4 energy7.043195meV/atom and force0.074783eV/A still exceed numerical budgets, although pair difference0.010197 passes. Do not call any of these four meshes converged. Numerical pilot labels are excluded from model training/development integration. All slabs share one crystallographic parent; they are related surface probes, not three independently sourced bulk phases.

## New A

Select repository linkengde/-, branch main, NEW independent cloud environment. Read this document, README, WORK_LOG_2026-10.md, RESTART_V19_HANDOVER.json and the B mesh ladder report. Use existing checkout; no worktree needed. Then:

```bash
cd /workspace/-
python3 research/mace-v12-transfer/coordination/claim_restarted_window.py window-a
```

The helper requires a new local instance identity, verifies handover hashes, preserves old owner/progress history and pushes ordinary owner transfer; no computation starts. Never manually erase identity to bypass an old-instance refusal. Verify/install only missing dependencies with existing environment docs/script. GPAW26.7/ASE3.29/gpaw-data1.2.1/MPI/ScaLAPACK/OpenMPI readiness and four actual cores required. Check disk/memory first. Training requires torch2.5.1CPU/mace0.3.16 only when needed. If Git native auth fails but GH_TOKEN binding is already injected, use git config --local credential.helper '!gh auth git-credential'; do not print/copy tokens.

After BOTH new owners have published IDs, A prepares and pins same-geometry5x5x1(A) and6x6x1(B) records with those new IDs; old manifests retain historical IDs. Explicitly register new labels before launch and use exact label selection. Fix owner/source/method/inventory and disk guards as in the approved pilot pack. Authorize only those two numeric jobs, one per machine, one BLAS thread/rank; do not repeat old meshes. Compare4/5,4/6 and5/6 after both verify. If still over budgets, analyze smearing separately, then vacuum/dipole sensitivity; do not silently change targets or lower gates. Provide independent common-target development labels and entry audit before V19 fitting. Keep fixed V18 recipe for the controlled training comparison unless evidence supports a separately recorded experiment.

## New B

Use another NEW independent cloud environment, same repo/main. Read this document and B_CONTINUOUS_WORK.md, then:

```bash
cd /workspace/-
python3 research/mace-v12-transfer/coordination/claim_restarted_window.py window-b
```

Verify missing-only environment setup. Wait for A's owner-bound6x6 job registration, then execute exact B label. Verify/publish compact result and run mesh/source/target-consistency analysis; proceed through ready queue automatically. Do not retrain or inspect sealed labels. Pending official V19 entry audit is dependency-blocked, not permission to reuse old acquisition parents as fresh validation.

## Retained files versus local-only state

All key inputs/labels/logs/models/selected checkpoints/reports are on main. Large old state.gpw and unselected checkpoints are unnecessary for these completed jobs and not required to redo transfer. Do not upload multiGB checkpoints. No old live job needs restart. User authorized cleaning useless local files after verified compact backup; retain active/unfinished states and frozen artifacts.
