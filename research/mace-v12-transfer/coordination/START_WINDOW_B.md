# Window B current assignment

V14 is SCREEN_FAIL. The v14 per-atom force localization and source-geometry audit are complete. Window A's AgC v15 targeted DFT label passed archive verification; A is running the assigned AgSi label, then will run AgTi.

## Run now: propose new v15 registry holdout geometries

This is geometry generation and screening only. Do not start DFT, MACE inference/training, alter existing data/holdouts, or run MD/TTM.

First sync and verify B ownership:

    cd /workspace/-
    git fetch origin main
    git merge --ff-only origin/main
    python3 research/mace-v12-transfer/coordination/sync_tasks.py identity
    python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote

The B identity must be c5035b48-f43b-4da0-b8e4-e2862817f86a. If it differs or fast-forward fails, stop and report it.

Publish progress before work:

    python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_registry_candidate_design --state running --iteration 0 --note "Generating and screening label-free new Ag registry candidates for v15 holdouts."

Use these parent motifs from the V14 training file, not any currently running A output:
- Ag-C: AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train
- Ag-Si: AgSi_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train
- Ag-Ti: AgTi_registry_probe_v13_01_periodic_PW_PBE_energy_force_v14_train

Generate label-free alternatives by translating the full Ag sublattice rigidly relative to the carbide/silicide framework. Use reproducible lateral shifts (screen several magnitudes around 0.4-0.7 Angstrom and choose a safe, clearly distinct registry); do not add random noise in this first geometry-only comparison. Preserve atom IDs, cell, PBC and central-pair tags. Strip any energy/force labels from proposal inputs.

Write candidate extxyz files, a source/input SHA256 manifest and a report under:
coordination/reports/window-b/v15_registry_candidates/

For each interface, report the Ag registry vector, marked-pair distance, shortest total and Ag-nonAg distances, nearest neighbors, displacement/RMS from source, and exact/near-duplicate checks against v14 train/valid/test, all v13/v14 holdouts/acquisitions and the existing v15 targeted candidates. Reject candidates with pair distances below 1.75 Angstrom or exact duplicates. Select one candidate per interface if it passes. Explicitly label these as new local registries derived from existing small-cluster motifs; they do not establish independent morphology or extended-interface coverage.

After publishing the candidate files and report:

    python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_registry_candidate_design --state completed --iteration 1 --note "Published screened label-free registry candidates; no DFT or model run."
    python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b

A will review candidates and decide which receive DFT labels. If the geometry axes or distinctness are ambiguous, document the issue and stop; do not start DFT.
