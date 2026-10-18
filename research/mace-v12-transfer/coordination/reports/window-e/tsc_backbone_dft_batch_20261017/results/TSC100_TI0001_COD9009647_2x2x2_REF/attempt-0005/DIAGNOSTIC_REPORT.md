# Independent SCF continuation: Ti-terminated reference

Status: **SCF_CONVERGED; result extraction and integrity checks passed.**

- Structure: `TSC100_TI0001_COD9009647_2x2x2_REF` (100 atoms, fixed geometry).
- Restart: attempt-0003 step-20 checkpoint; SHA-256 `38d99727d5b5b3f5152e918685e600790b4a36d819400b552b02519ab92eb2cb`. The supervisor verified that its hash was unchanged after the run.
- Input SHA-256: `b1948a44b89efe17ee7c2891b984167982067feeffa774b7f6fa77d9e9158d86`.
- GPAW 26.7.0 / ASE 3.29.0; PBE/PAW, plane-wave cutoff 500 eV, standard Monkhorst-Pack 6×6×1, Fermi–Dirac 0.10 eV, non-spin-polarized. Coordinates, cell, atom order, and PBC were verified unchanged.
- SCF converged at logical step 40 of 160 after 20 continuation steps; the flushed GPAW log has the matching convergence marker. Energy, all 100×3 force components, and the 3×3 stress tensor were saved after convergence and verified finite.
- Total energy: -798.275237007910 eV; force-consistent free energy: -799.090087579150 eV; maximum atomic-force norm: 1.419841 eV/Å. The execution lock defines no force-magnitude acceptance cutoff, so the force is reported as a measurement rather than pass/fail.
- Mixer: MSR1 beta=0.02. This is the sole change from the original locked beta=0.05, used to damp oscillation in the earlier failed run. All other mixer settings and physical calculation settings were preserved. The data role remains `TRAIN_FORCE_LABEL_CANDIDATE`.
- Runtime: 8988.9 s. cgroup peak memory: 17.23 GiB; OOM/kill events: zero. Intermediate GPW writes were disabled; the source checkpoint remains intact.

`FINAL_RESULT.json`, `FINAL_STRUCTURE.extxyz`, `CONVERGED_ENERGY_AND_EVIDENCE.json`, `FORCES_EVA.json`, and `STRESS_EVA3.json` contain the result. `RESULT_VALIDATION.json` records the integrity checks.
