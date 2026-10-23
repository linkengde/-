# V19 remaining data-source and method review

## Current result

This read-only follow-up identifies the exact remaining blocker for adding independent data. The 12 earlier inherited rows are now source-reconciled in the [archive reconciliation report](../v19_inherited_train36_archive_reconciliation_20261023/README.md). A qualified second Ti₃SiC₂ parent has still not been found. External structure endpoints were blocked by the environment's HTTP 403 tunnel, and public structure fixtures searched locally contained no Ti₃SiC₂ reduced-formula match. No substitute structure was fabricated and no DFT, MPI, or training job was run.

## Existing slab controls

The registered AgSi parent controls are already archived and verified. In the fixed geometry comparisons, both pass the current numerical budgets:

| Pair | Native energy difference | Free energy difference | Force-vector RMS | Maximum atom-vector difference |
|---|---:|---:|---:|---:|
| 26 Å, dipole OFF vs ON | 0.000737 meV/atom | 0.000941 meV/atom | 0.000845 eV/Å | 0.002078 eV/Å |
| 26 Å vs 36 Å, dipole ON | 0.006762 meV/atom | 0.006551 meV/atom | 0.001056 eV/Å | 0.002063 eV/Å |

These results support the provisional 26 Å, 6×6×1, σ=0.10 eV, z-nonperiodic, xy-dipole reference **for that frozen AgSi parent**. They do not test slab k-mesh or smearing convergence. The 4×4×1/6×6×1/8×8×1 and σ=0.05/0.10/0.20 eV controls previously run with fully periodic z do not answer the nonperiodic-z slab question. The source metrics and archive hashes are in [the boundary/dipole comparison](../v19_boundary_dipole_comparison/report.json).

## V17 rank-recovery row clarification

V19 train36 uses the archived Gamma-point V17 relabel, and the V19 row matches that output exactly in ordered atoms, positions, cell, energy, and forces. The older V16 training frame has the same geometry but an unverified legacy label:

- Legacy V16 stored energy: **−207.4639683217 eV/cell**.
- Verified Gamma relabel used by V19: **−212.0015102184 eV/cell**.
- Difference: **−4.5375418967 eV/cell** (−141.798 meV/atom).
- Force-component RMSE: **1.590919 eV/Å**; per-atom force-vector RMSE: **2.755552 eV/Å**; maximum atom-vector difference: **5.596606 eV/Å**.

The V16 metadata declares 1×4×2, but the original DFT output identity, SCF log, exact getter, and full settings were not recovered. Therefore the discrepancy cannot be attributed to the k-point mesh alone. V19 uses the directly verified Gamma output; do not replace it with the legacy V16 values.

## Second-parent source search

The local project search still finds only copies of COD 9009647, which are the same crystallographic source family. This review also inspected public source fixtures at fixed Git revisions:

| Repository | Revision | CIF/MCIF files scanned | Reduced Ti₃SiC₂ matches |
|---|---|---:|---:|
| `materialsproject/pymatgen` | `c6507866179a6c37ee7d0df6cbb79896a2ebba52` | 88 | 0 |
| `materialsproject/pymatgen-test-files` | `85005108e99e40cbd59e78f067459b0ca664efc6` | 110 | 0 |
| `cod-developers/cod-tools` | `0060d4652a5a262fbda2276801f1244060d67128` | 1,516 | 0 across 2,065 formula fields |

Direct Materials Project OPTIMADE, OQMD OPTIMADE, NOMAD OPTIMADE, COD, and the UGR COD mirror requests each failed with `Tunnel connection failed: 403 Forbidden`. This records an environment access block, not evidence that those databases contain no structure.

## Next action

Keep the V19 model out of dynamics and do not train against dev3. Find a source-verified Ti₃SiC₂ parent with a distinct crystallographic source family before creating another independent development label. Once a source is available, screen exact and near-duplicate geometry against train36 and dev3, select fixed slab settings, register the labels, and have B audit the entry before A approves a new split or fit.

The current V19 candidate still fails the interface force/separation gates and Ti₃SiC₂ backbone force screen. Recovered train sources and numerical controls do not change those scores.

## Reproduction

The local-data comparisons are reproduced from the repository root with:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_remaining_data_source_and_method_review_20261023/review_remaining_v19_gaps.py
```

The script reads the approved V19 training file, the exact V16 source training frame, the verified V17 Gamma output, and the existing fixed-parent control report. It does not open dev/test/sealed data or start DFT/MPI.
