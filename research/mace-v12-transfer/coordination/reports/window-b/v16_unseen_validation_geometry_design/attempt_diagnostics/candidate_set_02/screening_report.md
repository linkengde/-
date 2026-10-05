# V16 withheld-validation geometry design (label-free)

## Scope

Three geometries were generated from the existing V14 training parent motifs, one each for AgC, AgSi and AgTi. Only coordinates/IDs/cell/PBC, input manifests and the historical V15 force-localization report were used. No DFT calculation output, unseen DFT label, MACE model/inference/training, MD or TTM was accessed. The six planned V16 acquisition input geometries were included in deduplication; their calculation folders were not opened.

The proposals move the Ag sublattice laterally while leaving the carbide/silicide framework fixed. This changes the local Ag-neighbor environment sampled by framework atoms identified in historical localization, while avoiding claims of independent morphology. Each file contains geometry/identity metadata only and remains unassigned until A freezes and assigns the split.

## Selected candidates

| Interface | Lateral shift (Å) | Marked pair distance (Å) | Closest known geometry RMS (Å) | Main framework-neighbor response target |
|---|---:|---:|---:|---|
| AgC | (+0.522, -0.671, +0.000) | 1.989 | 0.087 | C#604804 Δnearest-Ag=-0.271 Å, Ti#602709 Δnearest-Ag=-0.808 Å |
| AgSi | (+0.717, +0.456, +0.000) | 3.086 | 0.358 | Si#366934 Δnearest-Ag=+0.646 Å, C#368838 Δnearest-Ag=-0.041 Å |
| AgTi | (+0.499, +0.688, +0.000) | 3.101 | 0.344 | Ti#479027 Δnearest-Ag=-0.478 Å, C#480975 Δnearest-Ag=+0.690 Å |

The nearest-geometry, no-exact/no-near check covered every current train/validation/test frame found under all model-version data directories, prior V13–V15 frozen input geometries, the six V16 acquisition inputs, and the earlier B label-free registry/wide-span proposals. Matching uses persistent IDs where possible and same-species assignment otherwise. Full paths, hashes and nearest hits are in `candidate_manifest.json` and `screening_details.json`.

## Geometry and contact audit

Candidates preserve the source atom count, atomic ordering, persistent IDs, cell and PBC. The only coordinate changes are a rigid xy translation of the Ag sublattice. All pairwise element minima use minimum-image distances and are reported by species pair. The short-contact screen is comparative: candidate minima may not be more than 0.05 Å below the smallest observed distance in the screened same-composition set. It is not a universal bond-length judgment; the framework's existing C–Ti minima are preserved.

No exact or near duplicate passed into the selected set. The nearest known geometry is reported for each candidate. The candidates are intentionally motif-correlated proposals, not independent-material validation. They do not test framework relaxation or thermal disorder; a separate paired framework-perturbation set would be needed to test that factor directly.

## For A

Treat all three as geometry proposals only. A should review the source lineage, hash/contact audit and intended validation role, then freeze accepted inputs before V16 training. Do not use them for training if their purpose remains an unseen validation set. Keep the planned six V16 acquisitions and all scored V15 probes in a distinct acquisition/history role.

## Reproduction

Run from the repository root in the ASE environment:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /workspace/.venvs/mace-v14-assist/bin/python -B \
  research/mace-v12-transfer/coordination/reports/window-b/v16_unseen_validation_geometry_design/generate_candidates.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v16_unseen_validation_geometry_design/reproduction_01
```

The script refuses nonempty output directories. The SHA inventory covers this generator, the three label-free extxyz files, manifest, screening details and this report.
