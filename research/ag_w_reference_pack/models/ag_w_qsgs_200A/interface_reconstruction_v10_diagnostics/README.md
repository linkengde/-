# Ag–W V10 diagnostics

V10 freezes the V9 final NVE restart and tests the remaining 5 Å mapping phases, local density/neighbor/Voronoi statistics, transfer-distance sensitivity, and the Ag/W cold-state temperature split.

## Gate results

- **5 Å capacity mapping: FAIL.** All four phases fail the shared direct-flow gate. The best residual is 1,181.24 Å³; a real LAMMPS mapping-only run confirms failure.
- **Transfer distance: diagnostic only.** At 40 Å the residual falls to 0.000311 Å³, but the maximum-flow gate still fails and the required transport is strongly nonlocal (34.87 Å weighted mean distance).
- **Cold baseline: FAIL.** LAMMPS and Python temperatures agree within 0.001 K. The reported Ag/W split is partly affected by common z drift, but remains after removing each group's COM drift. Potential energy continues to fall during NVT.
- **TTM/source: not run.** No TTM time step or laser source was executed.

See [V10_DIAGNOSTIC_REPORT.md](V10_DIAGNOSTIC_REPORT.md) for the full findings and [SHA256SUMS](SHA256SUMS) for file checksums.

## Reproduction

Run from the repository root with the Ag QSGS virtual environment used for V9:

```bash
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/scan_v10_5A_phases.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/scan_v10_transfer_distance.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/analyze_v10_local_structure.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v10_diagnostics/summarize_v10_resolution_comparison.py
```

LAMMPS diagnostics use `in.v10_spatial_geometry` and `in.v10_temperature_crosscheck`; both execute `run 0`. `mapping_5A_best_final` is the mapping-only preflight and exits before any TTM integration or source application.

No V9 input, original QSGS file, atom count, atom ID, box dimension, or atomic-volume constant was changed.
