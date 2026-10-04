# Environment validation status

Validated in the original cloud on 2026-10-04 (Asia/Shanghai):

- GPAW import/version: PASS, GPAW 26.7.0; ASE 3.29.0; `gpaw_data` setup path resolved. Actual PW-PBE references PASS for AgC 2.0 A (61 SCF iterations), AgC 2.4 A (61), and AgSi 2.8 A (49): all converged, finite energy/forces, exact input hashes verified, and compact outputs are in the repository. AgSi 2.2 A is being calculated in the new cloud window.
- Saved MACE model inference: PASS for v11 checkpoint `MACE_periodic_v11_energyfocus_run-44.model` on the first frozen test frame (20 atoms). Energy `-121.7960428793 eV`, forces shape `(20, 3)`, all energy/force values finite, Fmax `0.9223657 eV/A`.
- MACE training command availability: PASS, `mace_run_train --help` exit code 0 with writable XDG/Matplotlib/Fontconfig caches. This verifies the command interface only; no v12 fit was started.
- MACE environment: MACE 0.3.16, PyTorch 2.5.1+cpu, ASE 3.29.0.
- Foundation model: bundled `mace-mp-0b3-medium.model`; SHA256 is in `file_manifest.json`.

## Still required

The new window must finish AgSi 2.2 A and add its converged archive. The other three labels are already present in `calculations/`. Only after all four `summary.json` files report `scf_converged: true` should it run `archive_completed_pw.py`, `build_v12_dataset.py`, v12 training, and frozen/external evaluation. No long MD or TTM is authorized by this handoff.

## New-window handoff action

1. Open this repository on branch `main`, follow `environment_setup.md`, and run `./run_AgSi_d2p2.sh` from `periodic_interface_v4/pbe_interface_energy_additions_v12/`.
2. Require `summary.json` to report `scf_converged: true`; verify its source hash matches `input_manifest.json`.
3. Pull the latest `main` (it now contains the other three converged archives).
4. Commit and push only the compact AgSi_d2p2 archive files (`AgSi_d2p2_PW_PBE.extxyz`, `summary.json`, `gpaw.log`, `progress.json`); do not add `state.gpw`. Then run `archive_completed_pw.py`, build the dataset, train v12, and run frozen/external evaluation only after all four summaries pass.
