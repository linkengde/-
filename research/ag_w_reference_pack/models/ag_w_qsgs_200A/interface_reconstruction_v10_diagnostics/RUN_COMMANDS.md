# V10 reproduction commands

Run Python commands from the repository root. Run LAMMPS commands from `research/ag_w_reference_pack/models/ag_w_qsgs_200A/` so the relative paths in the inputs resolve.

```bash
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/scan_v10_5A_phases.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/scan_v10_transfer_distance.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/analyze_v10_local_structure.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/summarize_v10_resolution_comparison.py
```

The actual exploratory environment used `/workspace/.venvs/ag-qsgs/bin/python` and `/workspace/.venvs/ag-qsgs/bin/lmp` with four OpenMP threads. The equivalent LAMMPS calls are:

```bash
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/lmp \
  -var stateDump interface_reconstruction_v9_hybrid/v9_nve2_end.dump.gz \
  -var stateStep 2200 \
  -var output interface_reconstruction_v10_diagnostics/v10_voronoi_force_final.dump \
  -in interface_reconstruction_v10_diagnostics/in.v10_spatial_geometry \
  -log interface_reconstruction_v10_diagnostics/v10_voronoi_force_final.log

OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/lmp \
  -in interface_reconstruction_v10_diagnostics/in.v10_temperature_crosscheck \
  -log interface_reconstruction_v10_diagnostics/v10_temperature_crosscheck.log

OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/python ttm/agw_ttm_md_pilot.py \
  --data interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz \
  --audit interface_reconstruction_v10_diagnostics/v10_exploratory_state_manifest.json \
  --exploratory \
  --restart-dump interface_reconstruction_v9_hybrid/v9_nve2_end.dump.gz \
  --restart-step 2200 --output interface_reconstruction_v10_diagnostics/mapping_5A_best_final.json \
  --pulse-fwhm-ps 1 --pulse-center-ps 1 --cooling-ps 1 \
  --mapping-only --grid-A 5 --grid-offset-x-A 0 --grid-offset-y-A 0 \
  --max-map-transfer-A 20
```

Each LAMMPS spatial and temperature input uses `run 0`. The mapping-only pilot exits at the mapping gate; its JSON records zero TTM steps and no laser source. The two `.dump.gz` files are LAMMPS `run 0` snapshots containing forces and Voronoi cell values for independent inspection.
