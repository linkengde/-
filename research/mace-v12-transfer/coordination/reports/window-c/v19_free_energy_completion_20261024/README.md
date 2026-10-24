# V19 train36 free-energy candidate recovery — 2026-10-24

The V19 36-frame training candidate can carry a complete alternative force-consistent `free_energy_eV` key without running new DFT or changing the registered V19 dataset. The first 30 frames match the geometry and native-energy labels in the already audited V18 free-energy sidecar; force differences after extxyz serialization are below `5.0e-9 eV/Å`. The six V19 acquisitions already carry free energy and are reconciled here against their owner-verified output, summary, verification record, and GPAW log hashes. The native-minus-free energy differences range from 0.580422 to 1.166765 eV/cell, confirming these are distinct labels.

`v19_train36_free_energy_candidate.extxyz` preserves the existing `REF_energy` and `REF_forces` fields and ensures all 36 frames also have `free_energy_eV`. It remains a separate candidate sidecar. The V19 task manifest still says training is unauthorized; this recovery does not select the energy convention, approve training, replace `train.extxyz`, or provide free energies for V19 development frames.

The recovery script reads only the V19 `training_candidate/train.extxyz`, its manifest and six registered training archives, plus the previously audited V18 training sidecar and audit. It does not open validation, test, holdout, or sealed data. No DFT or training was run.

Reproduce from the repository root:

```bash
/workspace/.venvs/mace/bin/python research/mace-v12-transfer/coordination/reports/window-c/v19_free_energy_completion_20261024/recover_v19_free_energy_candidate.py
```
