# MACE v12 Ag/Ti3SiC2 handoff package

This lightweight package is intended to continue the Ag-Si/Ag-C DFT-label and MACE v12 workflow in another cloud window. It contains the four fixed interface input structures, GPAW scripts and settings, the exact v11 train/validation/test split, independent AgSi/AgC holdouts, and v9-v11 model files used by the comparison script.

## DFT scope and current handoff state

The four required PW-PBE single points are AgSi 2.2 A, AgSi 2.8 A, AgC 2.0 A, and AgC 2.4 A (GPAW 26.7, PBE, PW 500 eV, Gamma, 0.1 eV smearing). Their input structures and hashes are in `periodic_interface_v4/pbe_interface_energy_additions_v12/`.

The original cloud still owns any live GPAW processes. This repository does **not** include their rolling `state.gpw` checkpoints (about 3.9 GB total at the last check), because these are large and still being updated. The live jobs were launched from the original cloud workspace, so their automatic archives will land there first. After they finish, copy only the four compact result directories (each label `.extxyz`, `summary.json`, `gpaw.log`, and `progress.json`) into `periodic_interface_v4/pbe_interface_energy_additions_v12/calculations/`, then commit them here before starting dataset construction. If any calculation is missing after the original window ends, it can be rerun from its checked-in input structure.

## Continue after the labels are present

1. Use GPAW 26.7.0, ASE, and the matching GPAW PBE setups. Run `python archive_completed_pw.py` in the additions directory; it must verify all four converged results.
2. From `mace_periodic_v12_interface_energy/`, run `python build_v12_dataset.py`. It refuses to proceed if any reference is unconverged or the four-element energy-composition rank is not 4/4.
3. Install the MACE training environment and provide a compatible local foundation model with `FOUNDATION_MODEL=/path/to/mace-mp-0b3-medium.model bash run_v12_training.sh`. The training script is portable and no longer depends on the old cloud's absolute paths.
4. Run `python evaluate_v12_candidates.py` after training. It compares v9-v12 on the frozen test and independent AgSi/AgC holdouts.

The package does not include the external MACE foundation-model download or environment caches; those can be installed again in the new window. It also does not establish that v12 is validated for long MD or TTM.
