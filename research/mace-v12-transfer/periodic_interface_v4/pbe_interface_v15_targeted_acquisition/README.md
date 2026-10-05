# Residual-targeted v15 DFT acquisition candidates

This A-owned package prepares three fixed-geometry PW-PBE single points from the v13 force-localization report: one each for Ag-C, Ag-Si and Ag-Ti. Each input adds a small deterministic displacement to the local environment around the largest residual atoms and the marked contact pair. The inputs remain outside the v14 train/validation/test data and all three reserved v14 holdouts.

These labels are candidate v15 acquisition data. They are not a new validation set and do not establish transfer to extended interfaces, finite-temperature structures or production candidates. Use them for v15 only after reviewing the completed v14 evaluation; if needed for v15 training, reserve fresh blind geometries.

The generator checks protected v14 geometries and close contacts. The runner uses the already installed four-rank GPAW environment, publishes SCF progress, archives only compact results, and keeps large `state.gpw` files in `/tmp`.

```bash
cd /workspace/-
GPAW_PYTHON=/workspace/.venvs/gpaw-mpi/bin/python \
DFT_RUN_ROOT=/tmp/mace_v15_targeted_window_a \
bash research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v15_targeted_acquisition/run_parallel_labels.sh
```

The task remains owned by window A. Do not run the package from B or alter the frozen v14 dataset and holdouts.
