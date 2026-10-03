# Environment validation status

Validated in the original cloud on 2026-10-04 (Asia/Shanghai):

- GPAW import/version: PASS, GPAW 26.7.0; ASE 3.29.0; `gpaw_data` setup path resolved. The three AgSi/AgC fixed-geometry PW-PBE jobs are still running; their convergence is the pending full plane-wave validation. AgSi 2.2 A is queued for the next cloud window.
- Saved MACE model inference: PASS for v11 checkpoint `MACE_periodic_v11_energyfocus_run-44.model` on the first frozen test frame (20 atoms). Energy `-121.7960428793 eV`, forces shape `(20, 3)`, all energy/force values finite, Fmax `0.9223657 eV/A`.
- MACE training command availability: PASS, `mace_run_train --help` exit code 0 with writable XDG/Matplotlib/Fontconfig caches. This verifies the command interface only; no v12 fit was started.
- MACE environment: MACE 0.3.16, PyTorch 2.5.1+cpu, ASE 3.29.0.
- Foundation model: bundled `mace-mp-0b3-medium.model`; SHA256 is in `file_manifest.json`.

## Still required

The three existing PW-PBE calculations must converge and their compact archives must be copied here. The new window must run the missing AgSi 2.2 A point and add its converged archive. Only after all four `summary.json` files report `scf_converged: true` should it run `archive_completed_pw.py`, `build_v12_dataset.py`, v12 training, and frozen/external evaluation. No long MD or TTM is authorized by this handoff.

## New-window handoff action

1. Open this repository on branch `main`, follow `environment_setup.md`, and run `./run_AgSi_d2p2.sh` from `periodic_interface_v4/pbe_interface_energy_additions_v12/`.
2. Require `summary.json` to report `scf_converged: true`; verify its source hash matches `input_manifest.json`.
3. Commit and push only the compact AgSi_d2p2 archive files (`AgSi_d2p2_PW_PBE.extxyz`, `summary.json`, `gpaw.log`, `progress.json`). Do not add `state.gpw`.
4. Wait for the original window to upload the other three converged archives. Pull `main`, then run `archive_completed_pw.py`, build the dataset, train v12, and run the frozen/external evaluation.
