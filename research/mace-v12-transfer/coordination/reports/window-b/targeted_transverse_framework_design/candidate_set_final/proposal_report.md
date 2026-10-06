# Paired transverse and framework V16 acquisition proposals

## Purpose and limits

This bundle separates two training-data factors suggested by the matched V14/V15 diagnosis: transverse Ag registry and local framework-contact response. It contains geometry only. It does not assign DFT ownership or a final dataset role, read unseen labels, or run DFT, model inference/training, MD or TTM. A decides whether/when to label these points.

For each interface, a fixed parent geometry is the control. The transverse pair translates every Ag atom by ±0.15 Å in the xy direction orthogonal to the marked Ag–X vector. The framework pair leaves Ag fixed and moves only the marked non-Ag contact by ±0.04 Å along that atom's historical V15 force-residual direction. Thus Ag registry and framework displacement are separate design arms.

## Priority subset

| Interface | Ag transverse pair | Framework contact pair | Pair projection behavior |
|---|---|---|---|
| AgC | Ag IDs shifted ±0.15 Å; `AgC_transverse_Ag_minus/plus.extxyz` | C#604804 shifted ±0.04 Å; `AgC_framework_neighbor_minus/plus.extxyz` | transverse +: projection Δ +0.000000 Å, direct distance Δ +0.004972 Å; framework +: projection Δ -0.039715 Å, direct distance Δ -0.039710 Å |
| AgSi | Ag IDs shifted ±0.15 Å; `AgSi_transverse_Ag_minus/plus.extxyz` | Si#366934 shifted ±0.04 Å; `AgSi_framework_neighbor_minus/plus.extxyz` | transverse +: projection Δ -0.000000 Å, direct distance Δ +0.004606 Å; framework +: projection Δ +0.008002 Å, direct distance Δ +0.008316 Å |
| AgTi | Ag IDs shifted ±0.15 Å; `AgTi_transverse_Ag_minus/plus.extxyz` | Ti#479027 shifted ±0.04 Å; `AgTi_framework_neighbor_minus/plus.extxyz` | transverse +: projection Δ +0.000000 Å, direct distance Δ +0.004343 Å; framework +: projection Δ -0.034298 Å, direct distance Δ -0.034215 Å |

The lateral Ag arm preserves the marked-pair longitudinal projection to numerical precision; its direct Euclidean distance changes slightly at second order. Framework displacement can change the longitudinal projection and direct distance, and those changes are recorded separately. Each atom ID and displacement vector is listed in `proposal_manifest.json`.

## Contact and split checks

Each proposal keeps atom order/IDs, cell and PBC. Species-pair minima are checked against the empirical lower envelope from current train/valid/test geometries and official input geometries for the same formula, with a 0.05 Å allowance. This is a comparative screen, not a chemical bond definition. The frozen withheld candidates, six planned V16 acquisition inputs, all current split geometries and prior B proposals are also included in exact/near comparison.

These are deliberately small perturbations around their fixed control, so the parent is an intentional near-geometry match. The manifest separates near hits to that source/control lineage from any other-lineage near hits; it does not hide or call the pair independent. The source motifs are existing small clusters, not new material morphologies. Thermal disorder and extended-interface coverage remain open.

## Exact and near-geometry screening results

The periodic exact-geometry screen found zero exact input-geometry hits for all 12 proposals among the screened train/valid/test frames, official acquisition/holdout inputs, frozen V16 validation proposals and prior B geometries. It matches persistent IDs where available and uses same-species assignment otherwise. A separate conservative near flag is raised when both Ag-registry RMS and non-Ag framework RMS are at most 0.15 Å after MIC matching. These are review flags, not duplicate declarations; repeated frames copied between model versions and input archives are counted as separate file/frame hits. See each candidate's complete hit list and cell deltas in `proposal_manifest.json`.

| Proposal | Exact hits | Other-lineage near hits | Nearest near hit |
|---|---:|---:|---|
| AgC_transverse_Ag_m | 0 | 2 | `train.extxyz` frame 26 (all-atom RMS 0.0410 Å) |
| AgC_transverse_Ag_p | 0 | 0 | none |
| AgC_framework_neighbor_m | 0 | 5 | `AgC_registry_strain_acq_v14_01.extxyz` frame 0 (all-atom RMS 0.0058 Å) |
| AgC_framework_neighbor_p | 0 | 5 | `AgC_registry_strain_acq_v14_01.extxyz` frame 0 (all-atom RMS 0.0058 Å) |
| AgSi_transverse_Ag_m | 0 | 1 | `AgSi_wide_registry_p0p20.extxyz` frame 0 (all-atom RMS 0.0289 Å) |
| AgSi_transverse_Ag_p | 0 | 2 | `AgSi_widespan_registry_acq_v16_01.extxyz` frame 0 (all-atom RMS 0.0866 Å) |
| AgSi_framework_neighbor_m | 0 | 5 | `AgSi_registry_strain_acq_v14_01.extxyz` frame 0 (all-atom RMS 0.0056 Å) |
| AgSi_framework_neighbor_p | 0 | 5 | `AgSi_registry_strain_acq_v14_01.extxyz` frame 0 (all-atom RMS 0.0056 Å) |
| AgTi_transverse_Ag_m | 0 | 17 | `train.extxyz` frame 25 (all-atom RMS 0.0565 Å) |
| AgTi_transverse_Ag_p | 0 | 5 | `train.extxyz` frame 28 (all-atom RMS 0.1128 Å) |
| AgTi_framework_neighbor_m | 0 | 5 | `AgTi_registry_probe_v13_01.extxyz` frame 0 (all-atom RMS 0.0057 Å) |
| AgTi_framework_neighbor_p | 0 | 5 | `AgTi_registry_probe_v13_01.extxyz` frame 0 (all-atom RMS 0.0057 Å) |

AgTi_transverse_Ag_m has the largest overlap flag count (17 file/frame hits, many mirrored historical split copies); its closest listed frame is 0.0565 Å all-atom RMS away. It is not an exact duplicate, but A should review this branch against the historical split before assigning it for labeling. No candidate is near a frozen V16 validation proposal under this screen. Near-hit counts and species-specific contact checks do not establish independent statistical coverage.

## A review

Priority is one paired example of each arm per interface (12 geometries total). First review the AgSi marked-Si transverse pair, then the AgC marked-C and AgTi marked-Ti framework pairs and their opposite Ag shifts. Treat AgTi_transverse_Ag_m as review-only until A checks the close historical split hits. The largest noncontact shell residuals (AgC Ti#602709; AgTi C#480975) can be a second-stage framework pair if A needs separate shell-atom coverage; they were not added to this compact first set.

The fixed source parent is the control; its existing label is not duplicated here. Keep these candidates outside any frozen V16 validation role. Do not label before A approves the acquisition role and timing.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /workspace/.venvs/mace-v14-assist/bin/python -B \
  research/mace-v12-transfer/coordination/reports/window-b/targeted_transverse_framework_design/generate_training_proposals.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/targeted_transverse_framework_design/reproduction_01
```

The script refuses nonempty output. Hashes cover the generator, all 12 geometry files, manifest, screening summary and this report.
