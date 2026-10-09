# V18 AgSi public-reference diagnostic inference

## Scope and validity

This report screens the already selected V18 clean-core model against only the two task-authorized, public, archive-verified AgSi CIF references. It does not use sealed holdouts or test labels and does not train, fine-tune, edit model weights, launch DFT/MPI, or modify calculation inputs.

The two geometries come from one uploaded-CIF family. Their marked Ag–Si contact distances also overlap local V18 training coverage (reference about 2.7795 Å; nearest clean30 training motif about 2.8000 Å). This is a one-family diagnostic, not independent validation, a production pass, or evidence of transfer to other registries or compositions.

The selected epoch-79 `.model` is the inference artifact. The published epoch-79 checkpoint, training input, dataset manifest, selection record, and artifact manifest are hash-checked before inference. The script reads only `data/train.extxyz` for overlap screening; development and test structure files remain unopened.

## Runtime

Inference uses the isolated `/workspace/.venvs/mace-cpu` environment with Python 3.12, PyTorch 2.5.1 CPU, mace-torch 0.3.16, and ASE 3.29.0. No packages were installed in or changed within `/workspace/.venvs/gpaw-mpi`.

Post-install inventory confirmed `/workspace/.venvs/gpaw-mpi` still has Python 3.12.14, GPAW 26.7.0, and ASE 3.29.0, with no `torch` or `mace-torch`; the MACE environment has Python 3.12.14, PyTorch 2.5.1+cpu, mace-torch 0.3.16, and ASE 3.29.0. The installed MACE runtime also records NumPy 2.5.3 and SciPy 1.18.1 in `environment.json`.

## Reproduce

From the repository root:

```bash
/workspace/.venvs/mace-cpu/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/infer_and_evaluate.py
```

The script rechecks candidate lineage hashes, selection safeguards, the two archive checksum lists and stored verifier results, source geometry hashes, atom order/IDs, and V18 train-only geometry overlap before loading the model. It writes `metrics.json` and `environment.json` alongside itself.

## Metrics

- All-atom force-vector RMSE is `sqrt(mean_i(sum_xyz((F_model[i] - F_DFT[i])**2)))` in eV/Å; the maximum atom-vector error and per-species values are also recorded.
- The marked-pair signed separation force is `dot(F_X - F_Ag, minimum_image_unit_vector_Ag_to_X)`. Model, DFT, signed error, and absolute error are reported separately.
- The net and mean normal force of all Ag atoms are reported along +z, with signed model-minus-DFT errors.
- Model energy is compared separately with the DFT native extrapolated energy and the DFT free energy. No energy gate is applied, and the two DFT conventions are not mixed.
- Provisional diagnostic gates are force-vector RMSE ≤0.05 eV/Å and absolute marked-pair separation-force error ≤0.10 eV/Å for each reference. These thresholds do not make this related-CIF screen independent validation.

`metrics.json` contains per-reference values and model/input/reference hashes. `environment.json` records package versions and CPU runtime details. `SHA256SUMS.txt` covers this README, the inference script, metrics, and environment record.

## Result

The CPU smoke run passed on both 52-atom references, with finite energies and forces. Both references fail the provisional force and separation-force gates by a wide margin:

| Public reference | Force-vector RMSE (eV/Å) | Max atom-vector error (eV/Å) | DFT → model separation force (eV/Å) | Signed separation error (eV/Å) |
|---|---:|---:|---:|---:|
| `AgSi_COD9009647_pilot_k6x6` | 2.63375 | 4.55774 | −1.53713 → −3.25427 | −1.71714 |
| `AgSi_COD9009647_vacuum26_k6x6_sigma0p10` | 2.63351 | 4.55706 | −1.53874 → −3.25427 | −1.71553 |

The all-Ag net +z force is −4.91663 eV/Å from MACE versus −2.63378 and −2.63638 eV/Å from DFT, respectively. Native-extrapolated energy errors are about −652.25 meV/atom; free-energy errors are about −644.74 meV/atom. These energy comparisons remain separate and ungated. The force and separation failures make the V18 model unsuitable even for this provisional AgSi force screen; the result does not support independent validation or production use.
