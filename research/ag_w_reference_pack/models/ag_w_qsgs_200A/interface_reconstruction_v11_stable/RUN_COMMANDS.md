# V11 run commands

Run from the repository root in the existing Ag–W environment. The commands below are prepared but have **not** been executed in this chat session.

## 0. Quick code checks

```bash
/workspace/.venvs/ag-qsgs/bin/python -m py_compile \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_reference.py \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/scan_v11_padded_z_grid.py \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/test_v11_padded_z_smoke.py

/workspace/.venvs/ag-qsgs/bin/python \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/test_v11_padded_z_smoke.py
```

Both commands must pass before running the full structure scan.

## 1. Padded-z mapping phase scan

```bash
/workspace/.venvs/ag-qsgs/bin/python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/scan_v11_padded_z_grid.py
```

Expected output:
`research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/v11_padded_z_mapping_scan.json`

The scan uses the frozen V9 restart and checks all four 5 Å phases on a 40×40×34 mesh spanning z=35–205 Å. It performs no MD integration, no TTM step, and no laser source.

**Decision:** do not continue unless (a) all static Newton atoms are covered by the padded z domain and (b) at least one phase passes the same 20 Å direct production max-flow solver.

## 2. Actual LAMMPS mapping-only preflight

Use the best phase from the JSON. The following shows the (0,0) Å phase; replace the offsets if another phase is best.

```bash
LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/python \
  research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py \
  --data research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz \
  --audit research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/v9_exploratory_state_manifest.json \
  --exploratory \
  --restart-dump research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/v9_nve2_end.dump.gz \
  --restart-step 2200 \
  --output research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v11_stable/v11_mapping_only.json \
  --pulse-fwhm-ps 0.01 --pulse-center-ps 0.02 --cooling-ps 0.02 \
  --mapping-only --grid-A 5 \
  --map-z-min-A 35 --map-z-max-A 205 --allow-above-surface \
  --grid-offset-x-A 0 --grid-offset-y-A 0 \
  --box-xy-A 200 --max-map-transfer-A 20
```

This must return `AGW_MAPPING_PREFLIGHT_PASS`, report zero excluded static-Newton atoms, and report no cell above total fill 1 within the solver tolerance. It still executes zero TTM steps and applies no source.

## 3. Rebuild the 300 K cold baseline

Only after the padded-grid mapping scan is successful, run:

```bash
cd research/ag_w_reference_pack/models/ag_w_qsgs_200A
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/lmp -pk omp 4 \
  -in interface_reconstruction_v11_stable/in.v11_cold_reinitialize \
  -log interface_reconstruction_v11_stable/v11_cold_reinitialize.log
cd ../../../../
```

The input minimizes the V9 geometry with the bottom fixed, creates a fresh common 300 K Maxwell distribution for all mobile atoms, uses one shared NVT thermostat for 2 ps, and follows with 0.1 ps NVE. It does not force Ag and W to equal temperatures using separate thermostats.

If the final 0.5 ps NVT tail still has a strong potential-energy slope or a large Ag/W temperature split, extend the same NVT state in 1 ps increments up to 5 ps before making a cold-baseline decision.

## 4. TTM source/no-source pair — only after cold PASS

Create a new exploratory manifest for the exact V11 NVE dump. It must contain the new dump SHA-256 and actual timestep. Do not reuse the V9 restart hash.

Run the no-source control and the source case from the **same** V11 NVE restart, with identical mesh, seed, pulse envelope, and interface-conductance assumption. Only `--source-scale` differs:

- Control: `--source-scale 0`
- Gaussian source: `--source-scale 1` for the initial code-path smoke test.

Keep:
- grid 5 Å, z=35–205 Å;
- the best accepted x/y phase;
- max direct transfer 20 Å;
- nominal source plane z=200 Å with explicit above-plane clamp;
- I0=5 eV/(ps Å²), R=30 Å, skin depth=20 Å;
- FWHM=0.01 ps, pulse centre=0.02 ps, cooling=0.02 ps;
- provisional Ag/W electronic interface conductance 1 GW m⁻² K⁻¹, explicitly marked exploratory.

Do not describe this short smoke test as a laser-melting experiment or as quantitative heat-transfer validation.
