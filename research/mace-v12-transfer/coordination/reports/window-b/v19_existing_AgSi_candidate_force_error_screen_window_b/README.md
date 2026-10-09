# Existing MACE candidate eligibility for AgSi force screening

## Decision

The V18 clean-core model and selected epoch-79 checkpoint have matching, published SHA256 records. A train-only comparison found no exact full-geometry match and no training frame with the same atom count, composition, and PBC as either CIF reference. There is local contact overlap: the target Ag–Si marked-pair distance is 2.779511 Å, and the nearest clean30 training pair is 2.800000 Å. Any future result is therefore a one-family development diagnostic, not independent validation.

**Inference was not run.** Neither PyTorch nor `mace-torch` is installed in the only available Python environments. The assigned task explicitly forbids package installation. `eligibility.json` records force, species, marked-pair, and Ag-layer normal-force metrics as `null` with status `NOT_RUN_RUNTIME_MISSING`; these values do not mean zero error.

## Model and training-data lineage

- Candidate: `MACE_periodic_v18_clean_core.model`, selected epoch 79 from the fixed 80-epoch clean30 run.
- Model SHA256: `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`.
- Published epoch-79 checkpoint SHA256: `d9faab45f32df7366b32f44a9f91c1aab9c72f2df422f24c9996fa3243b2c137`.
- Training input: 30 frames; SHA256 `20279d8aae53baf186546eb4a4dc396bb71e92de9d1c874d1c6fa0022313c045`.
- Dataset manifest SHA256: `36fb89632293cd5fcd494a4ad7665aeabb575ee1f6cc7633fa913ff7a1f1bdec`; it records 30 training frames, 3 development frames, no test split, and no opened withheld labels.
- Selection record says epoch 79 was recovered exactly from the checkpoint and reviewed by A. The original launcher exit receipt is unavailable (`training_exit_code: null`), which remains a lineage limitation.

The only newer V19 manifest found is `mace_periodic_v19_residual_response/training_candidate/manifest.json`. It is marked `TRAINING_CANDIDATE_ONLY_BLOCKED_PENDING_INDEPENDENT_DEVELOPMENT`, has `training_authorized=false`, and has no executable `.model`, `.pt`, or `.safetensors` file. It is not a model candidate for inference.

The audit read only the non-sealed `data/train.extxyz` and its manifests. It did not read development or test structure files or any sealed labels. Hash checks and structure metadata are reproducible from `audit_candidate_eligibility.py`.

## Reference archives

Both requested references have 52 atoms with formula `C16Ag4Si8Ti24`; their archive checksums, convergence flags, source hashes, and verifiers pass.

| Reference | Mesh | Smearing | PBC | Source SHA256 | Role |
|---|---|---|---|---|---|
| `AgSi_COD9009647_pilot_k6x6` | 6×6×1 | 0.10 eV | `[true,true,true]` | `8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c` | CIF numerical-convergence pilot |
| `AgSi_COD9009647_vacuum26_k6x6_sigma0p10` | 6×6×1 | 0.10 eV | `[true,true,true]` | `81e73a08a241634fc751ed51323177576f502391049ff9476c4c7acf5665f46d` | Same CIF family, expanded z cell |

These two structures are related controls from one uploaded-CIF family, not independent holdouts. Their DFT native and free energies remain separate in the archive and must not be mixed.

## Overlap assessment

The V18 training pool contains 30 frames. None has the target's `C16Ag4Si8Ti24` composition and 52-atom count, so there is no same-composition global geometry candidate for near-RMSD matching. Neither target's full geometry hash matches a training frame, and neither reference source hash appears in the training manifest.

Both targets have a marked Ag–Si distance of 2.779511 Å. Clean30 frame 1 (`AgSi_d2p8_periodic_PW_PBE_energy_force`) has a marked Ag–Si distance of 2.800000 Å, only 0.020489 Å away. This is local motif overlap even though the full formula, cell, atom count, and geometry differ. A future score must be described as a diagnostic on a related AgSi contact family, not a clean blind test.

## Metrics and limitations

If an existing MACE runtime becomes available under a separately authorized task, keep the references and energy conventions separate and use the project's provisional gates:

- All-atom force-vector RMSE: `sqrt(mean_i(sum_xyz((F_model[i] - F_DFT[i])**2)))` ≤0.05 eV/Å.
- Marked-pair separating force: `dot(F_X - F_Ag, unit_vector_Ag_to_X)`; report signed model and reference values and the absolute error, gated at ≤0.10 eV/Å.
- Also report errors by species, maximum atom-vector error, marked-pair identity/distance, and Ag-layer force projected along the slab normal (z here).

No force, species, pair-separation, or Ag-layer-normal model errors are available from this audit because predictions were not generated. A one-family pass would not establish independent validation, transfer across registries, or production readiness.

## Reproduce

From the repository root, run the read-only eligibility and training-overlap audit with:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_window_b/audit_candidate_eligibility.py
```

This script checks hashes and archive records and reads only the non-sealed V18 training frames plus the two public compact DFT references. It does not import MACE or PyTorch, run inference, train, or launch DFT/MPI. The output is `eligibility.json`; `SHA256SUMS.txt` covers the README, script, and output.
