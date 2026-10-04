# MACE v14 readiness

The v14 construction, training and evaluation scripts are prepared. They have not run because the independent DFT labels have not all arrived.

Window A's two v13 checks and window B's three acquisition geometries are scored first with `evaluate_v13_acquisition.py`. Then `build_v14_dataset.py` imports exactly those five verified labels and keeps the three reserved v14 holdouts outside training. It preserves the v13 validation/test files byte-for-byte and requires rank 4/4.

After dataset construction, run `run_v14_training.sh`; it mirrors the v13 training settings with seed 45 and 80 epochs. The trainer refuses to overwrite an existing v14 model. Freeze that checkpoint, then run:

```bash
python evaluate_v14_cycle.py v14
```

The v14 evaluation includes the three reserved interface checks and frozen test set, records hashes, and reports energy, force-vector and separating-force errors per interface. Passing these small cluster-motif checks alone does not establish transfer to a periodic extended interface, high temperature, liquid Ag or pressure conditions. See the staged gates in `../../coordination/ITERATION_PROTOCOL.md`.
