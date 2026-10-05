# Window B current assignment

Window B's read-only V15 dataset-integrity audit is complete and published. Window A is running the AgSi V15 registry blind DFT holdout. B now has a separate geometry-only task to prepare wider-span candidate structures for possible future **training acquisition**. The proposals are not DFT labels and are not blind tests.

## Assignment: wider Ag-X distance and registry sampling design

Use the V14 force errors, force-localization report, V15 force-separation plan, B's dataset audit, and committed PW-PBE training geometries to design a broader next acquisition. The user requested a wider exploratory span, followed by narrowing if the outer points prove excessive.

First synchronize and verify this is the registered B machine:

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py identity
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

The B identity must be `c5035b48-f43b-4da0-b8e4-e2862817f86a`. Confirm `v15_wide_span_sampling_design` is assigned in `window-b.json`. If identity differs, fast-forward fails, or the task is absent, stop and report it.

Publish the start event:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_wide_span_sampling_design --state running --iteration 0 --note "Designing broader label-free Ag-X distance/registry candidates from committed training references; no DFT or model run."
```

Design separate axes for each interface:

1. **Ag-X normal contact distance.** Starting from a committed, labeled PW-PBE training parent, propose a deliberately wider coarse bracket. Use offsets on the order of ±0.20 and ±0.35 Å from a documented anchor where the species-specific geometry supports them; choose interface-specific values if one endpoint is unsafe. Do not concentrate only on near-neighbor distances already represented in training.
2. **Lateral Ag registry.** Propose clearly distinct rigid translations, evaluated separately from the normal-distance scan so A can tell which change helps.
3. **Local surroundings.** Identify whether a small number of additional neighbors around the Ag-X contact need controlled perturbations to cover the Ag-C/Ag-Ti contact shells or the broader Ag-Si error. Keep these separate from the rigid registry proposals.

Generate label-free candidate extxyz frames only if the parent/IDs and cell allow an unambiguous edit. Preserve all atom IDs, cell and PBC; remove all energy/force labels. For every frame report its parent/hash, target and actual marked-pair distance, translation vector, species-resolved shortest contacts, affected framework neighbors, displacement/RMS, and exact/near-duplicate checks against existing train/valid/test, prior acquisitions and frozen holdout **inputs**. Never read or use any V15 holdout energy/force label, summary or DFT output while designing candidates. Do not inspect A's `/tmp`, active job files or checkpoints.

Use species- and framework-role-aware contact screening; do not apply a universal all-pair cutoff that rejects bonded C-Ti framework neighbors. Flag rather than force any endpoint with unsupported Ag-X contacts, severe framework damage, or an ambiguous marked pair. Recommend where to narrow the span if an outer proposal looks excessive. State clearly that the candidate structures are not validated training data until A selects and labels them.

Do not start DFT, MACE inference/training, MD or TTM. Do not change any dataset, builder, archive, or the active/frozen holdout inputs. Publish only under:

```text
research/mace-v12-transfer/coordination/reports/window-b/v15_wide_span_sampling_design/
```

Include a design report, candidate extxyz (if generated), JSON manifest with source/input hashes and screening details, and SHA256 list. Then finish:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_wide_span_sampling_design --state completed --iteration 1 --note "Published broader-span label-free Ag-X distance/registry proposals and geometry audit; no calculations."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```

Window A decides whether any proposal receives a DFT label after V15 blind evaluation. Do not call V15 passed or training-ready based on this geometry-only design.
