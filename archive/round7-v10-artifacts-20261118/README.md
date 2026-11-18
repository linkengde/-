# Round7 POD/SNAP v10 artifact archive

Archive branch: `codex/round7-v10-artifacts-20261118`  
Base: repository `main` at `69ce2715d85d8ce8c0046437f35edc731b708ea6`  
Scope: existing exploratory v10 model files, fit-ready data, the exact 347,464-atom structure used by the v10 static check, and completed test/static-force evidence. No model was refit and no simulation was started for this archive.

## Restoring files

From a fresh checkout of this branch:

```bash
git clone --branch codex/round7-v10-artifacts-20261118 https://github.com/linkengde/-.git
cd -/archive/round7-v10-artifacts-20261118
sha256sum -c SHA256SUMS.txt
cp -a round7_pod_snap /workspace/potential-fit/
```

The final `cp` restores the archived paths under `/workspace/potential-fit/round7_pod_snap`, matching the absolute paths in the original fit and LAMMPS input files. `MANIFEST.json` lists every archived file, its original source path, byte size, and SHA-256.

## Contents and labels

- `round7_pod_snap/dataset_v10_v25_supplement/`: 128 training, 12 development-validation, and 2 historical challenge configurations. The challenge split was excluded from fitting. There is no fresh untouched test set.
- `round7_pod_snap/models/final_v10/POD/`: POD descriptor, fitted coefficients, exact fit inputs, logs, error tables, and evaluation outputs.
- `round7_pod_snap/models/final_v10/SNAP/`: SNAP parameters and coefficients, FitSNAP configuration/log/metrics, the exact aggregated `Train.xyz` and `Validation.xyz`, and evaluation outputs.
- `round7_pod_snap/models/final_v10/finite_difference_check/`: saved LAMMPS energy-force consistency inputs, outputs, script, and summary.
- `round7_pod_snap/models/final_v10/static_347k_preflight/`: the complete v10 `run 0` force dumps, inputs, logs, and summaries for both potentials.
- `round7_pod_snap/models/target_200x200x300_qsgs/`: canonical 200×200×300 Å, 347,464-atom model, QSGS build inputs/mask, geometry report, and the exact timestep-55 `start_wrapped.data` read by the v10 static preflight.
- `round7_pod_snap/models/final_v10/pod_minimize_347k*`: selected FIRE inputs, logs, audits, and stop/failure reports. Checkpoint and initial force dumps are excluded.

The v10 dataset manifest records native PBE total-energy and force conventions, source files/hashes, splits, and atom counts. Fit units are eV, eV/Å, and LAMMPS `metal`; model element order is Ag, C, Si, Ti. The full-model LAMMPS static inputs explicitly map source data types Ti/Si/C/Ag to model types Ag/C/Si/Ti.

## Results and limits

The saved v10 `run 0` used 347,464 atoms and did not integrate MD steps. It reports POD Fmax 51.070 eV/Å (F99 13.009 eV/Å) and SNAP Fmax 138.637 eV/Å (F99 31.919 eV/Å). The test summaries and complete force arrays are included.

The conservative POD FIRE attempt stopped at step 52 after force rebound and worsening local contacts. It did not produce a converged final structure; the archive does not mark the structure as equilibrated or heating-ready. SNAP v10 has no matching full-model FIRE pass. The historical `attempt_008` TTM trajectory used POD v9; it is not included or presented as a v10 trajectory. Only its timestep-55 structure file is archived because the v10 static input used that exact file.

## Deliberately excluded

- Four current v11 Si hotspot DFT jobs, their live logs, and all v11 candidate outputs.
- FIRE checkpoint dumps and redundant initial dumps.
- The `attempt_008` TTM trajectory/final trajectory and unrelated Round7/MTP files.
- Temporary build trees, Python caches, and unrelated model versions.

All individual files are below GitHub's 100 MB file limit; the coordinate files are text, so this archive does not require Git LFS.
