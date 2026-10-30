# V19 train49 Ag/Ti₃SiC₂ interface response experiment

This versioned experiment appends the six verified v22 interface force labels to the byte-identical train43 candidate. The three-frame frozen dev3 is copied byte-for-byte. It is an experiment-only fit and does not modify official train36/dev3 or use test/sealed data.

## Added calculations

Six fixed-geometry GPAW labels cover AgTi, AgSi and AgC registry-B interfaces with rigid Ag-layer displacements of −0.20 Å and +0.20 Å along the surface normal. Each contains 52 atoms and uses PBE, a 500 eV plane-wave cutoff, 6×6×1 k points, 0.10 eV Fermi smearing, xy dipole correction and four MPI ranks. Each archived source result passed its archive verifier. Result, input and archive hashes are recorded in the generated dataset manifest.

The candidate set has 49 training frames and the original three frozen development frames. A byte-prefix check preserves all train43 data, and exact geometry duplicates against train43/dev3 are rejected during dataset construction. The six new configurations remain acquisition candidates from one source family; they are not official data assignments.

## Fixed training and screen

The fit repeats the train43 V19 recipe: seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, and REF_energy/REF_forces targets. The selected epoch-79 model is evaluated only against the frozen dev3. Existing gates remain per-interface force-vector RMSE ≤0.05 eV/Å and separation-force error ≤0.10 eV/Å.

Build and preflight the experiment data before fitting:

```bash
/workspace/.venvs/mace/bin/python prepare_train49.py
/workspace/.venvs/mace/bin/python preflight_train49.py
bash run_train49.sh
```

Dataset construction and preflight passed on 2026-10-30. The 80-epoch fit completed at 09:55 UTC, and the independent verifier confirmed full epoch-log/checkpoint coverage and tensor identity between the epoch-79 checkpoint and exported model. After this verification, intermediate epoch checkpoints were pruned; epoch 79 is retained.

## Frozen-dev3 result

The candidate fails the per-interface force-vector RMSE gate (≤0.05 eV/Å) for all three interfaces. All separation-force errors pass the ≤0.10 eV/Å gate. Aggregate force-vector RMSE is 0.09039 eV/Å, higher than train43's 0.06003 eV/Å. The fixed final epoch was used; no alternate epoch was selected using dev3.

| Interface | Force-vector RMSE (eV/Å) | Gate | Separation-force error (eV/Å) | Gate |
| --- | ---: | --- | ---: | --- |
| AgC | 0.08906 | Fail | 0.06556 | Pass |
| AgSi | 0.09129 | Fail | 0.07930 | Pass |
| AgTi | 0.09081 | Fail | 0.07230 | Pass |

Compared with train43, the force-vector RMSE rose by 0.03501 eV/Å for AgC, 0.02017 for AgSi and 0.03761 for AgTi. The AgTi separation-force error improved from 0.15869 to 0.07230 eV/Å; AgC and AgSi separation errors increased but remain below threshold. Frozen-dev3 residual analysis shows Ti force-vector RMSE increased for each interface: AgC 0.0315→0.1068, AgSi 0.0786→0.1118 and AgTi 0.0628→0.1176 eV/Å. The Ag component improved for AgC and AgTi. These are correlations on the development set, not proof that the added Ag labels caused the Ti errors. The complete atom- and distance-resolved comparison is in run/results/train49_residual_diagnostics.json.

The next acquisition will target the topmost Ti orbit in each registered training interface geometry with paired ±0.02 Å normal displacements, while retaining train43 as the base for a separate experiment. Those inputs will be derived from v22 training geometries; no dev3 coordinates or labels will be copied into training.

The result status is TRAIN49_EXPERIMENT_DEV3_SCREEN_FAIL. The model remains an experiment artifact and is not suitable for promotion to the official split, production use or molecular dynamics. Native extrapolated energy and finite-smearing forces follow the established candidate convention; exact energy-force gradient consistency is not claimed. These three frozen structures do not establish blind-test performance, and slab k-point, smearing and vacuum convergence remain incomplete. See run/results/train49_verification.json and run/SHA256SUMS.txt for the complete evaluation and artifact inventory.
