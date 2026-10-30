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

Dataset construction and read-only preflight passed on 2026-10-30. The train49 fit is now running with the fixed recipe; fit outputs and frozen-dev3 evaluation will be archived in a follow-up commit. Native extrapolated energy and finite-smearing forces follow the established candidate convention; exact energy-force gradient consistency is not claimed. Passing frozen dev3 would not establish blind-test performance or production/MD readiness.
