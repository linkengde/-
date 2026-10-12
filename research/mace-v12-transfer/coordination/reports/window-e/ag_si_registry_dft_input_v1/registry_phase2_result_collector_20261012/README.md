# Ag–Si phase-2 result collector

`collect_registry_phase2_results.py` reads only already-written result files for the C-reviewed registry IDs `R1` and `R2`. It does not import GPAW/ASE, execute calculation inputs, change structures or parameters, or modify the input result tree. It writes `REGISTRY_PHASE2_RESULTS.json` beside the script unless `--output` is supplied.

The phase-2 mesh is locked to `8×8×1`. A result with another mesh is marked `KPOINT_MISMATCH`; its numerical labels are withheld. The collector accepts a `result.json` sidecar and/or a GPAW `gpaw.log`, plus optional `C_EXECUTION_RELEASE_RECORD.json`, `SCF_STATUS.json`, and `RUN_STATUS.json` metadata in each registry directory. It records energy, force-consistent free energy, force vectors and status, stress tensor and status, SCF convergence and step count. Multiple result files/logs are treated as ambiguous instead of selecting one silently.

Example, with an existing read-only phase-2 result tree:

```bash
python collect_registry_phase2_results.py \
  --phase2-root /path/to/phase2_8x8x1 \
  --output /path/to/REGISTRY_PHASE2_RESULTS.json
```

`phase2_runner.py` is a separate, fixed-scope runner for R1/R2 at `8×8×1`, fixed geometry, PBE/PAW, four MPI ranks, one OpenMP thread, and `REGISTRY_SCREENING_ONLY`. Its default command is a read-only hash preflight:

```bash
python phase2_runner.py --preflight
```

Execution requires a future C-issued `REGISTRY_PHASE2_EXEC_AUTH.json`, its matching SHA-256, an already-existing empty work root, and the exact SHA bindings in `PHASE2_EXEC_AUTH_REQUEST.json`. While the request remains `PREPARED_NOT_ISSUED`, `--execute` fails closed. The runner determines SCF convergence and iteration count from the GPAW log; it never reads `ASECalculator.scf`.

At preparation time, this checkout contains no phase-2 result directory or issued C execution authorization. The generated JSON therefore lists both reviewed IDs, preserves the known frozen input/structure hashes, leaves `config_sha` and result labels null, and reports `RESULTS_NOT_FOUND`. The phase-1 6×6×1 configuration hash is intentionally not substituted for the phase-2 parameter-lock SHA.
