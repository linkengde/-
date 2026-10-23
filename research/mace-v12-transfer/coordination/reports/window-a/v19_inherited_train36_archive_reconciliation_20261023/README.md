# V19 train36 inherited archive reconciliation

> **Further follow-up:** The separate V17 rank-recovery row is clarified in [the remaining-gap review](../v19_remaining_data_source_and_method_review_20261023/README.md). V19 uses the archived Gamma relabel exactly; its historical V16 value differs substantially but lacks original SCF/output provenance, so the difference cannot be assigned to k-point sampling alone.

## Result

The 12 V19 train36 rows previously marked as lacking independently verified original settings were reconciled to their existing source inputs, GPAW outputs, summaries, logs, and archive manifests. **All 12 passed** input/output/summary SHA-256 checks, exact geometry and atom-ID checks, exact training-force comparison, energy comparison, and SCF convergence checks. No DFT or MPI job was run.

This updates the earlier lineage audit's unresolved-archive finding. The earlier report remains a record of the evidence available during that audit; this follow-up provides the recovered evidence.

| Train index | Label | Historical archive group | SCF iterations | GPAW backend | Energy (eV/cell) |
|---:|---|---|---:|---|---:|
| 0 | `AgSi_d2p2` | V12 contact | 47 | cgpaw, domain 4 | -206.9663089771 |
| 1 | `AgSi_d2p8` | V12 contact | 49 | serial, domain 1 | -207.4892339445 |
| 2 | `AgC_d2p0` | V12 contact | 61 | serial, domain 1 | -259.1583543152 |
| 3 | `AgC_d2p4` | V12 contact | 61 | serial, domain 1 | -258.9332034547 |
| 4 | `AgC_validation_d2p30` | V12 distance scan | 60 | cgpaw, domain 4 | -259.0564822555 |
| 5 | `AgSi_validation_d2p60` | V12 distance scan | 49 | cgpaw, domain 4 | -207.5280818326 |
| 6 | `AgTi_validation_d2p60` | V12 distance scan | 37 | cgpaw, domain 4 | -177.0164599820 |
| 7 | `AgC_lateral_registry_holdout_v13` | V13 lateral holdout, later included in train | 62 | cgpaw, domain 4 | -258.9555196097 |
| 8 | `AgSi_lateral_registry_holdout_v13` | V13 lateral holdout, later included in train | 46 | cgpaw, domain 4 | -207.5444919019 |
| 9 | `AgTi_registry_probe_v13_01` | V14 acquisition | 38 | cgpaw, domain 4 | -176.5798462749 |
| 10 | `AgC_registry_strain_acq_v14_01` | V14 acquisition | 58 | cgpaw, domain 4 | -257.5833982482 |
| 11 | `AgSi_registry_strain_acq_v14_01` | V14 acquisition | 49 | cgpaw, domain 4 | -208.6934581847 |

All twelve summaries record GPAW 26.7.0, PBE plane waves, 500 eV cutoff, Gamma `1×1×1`, and Fermi smearing 0.10 eV. The four historical runner scripts call ASE `atoms.get_potential_energy()` and store that value as `PW_PBE_energy_eV` / `energy_eV_cell`; they do not call `get_free_energy()`. The archived extxyz files do not contain a separate free-energy value.

Nine logs report the GPAW `cgpaw` backend with `domain=4`; three V12 contact logs report the serial backend with `domain=1`. The serial runs remain valid fixed-geometry PBE labels; this is a runtime execution difference, not a difference in the recorded DFT method. Rank count is not separately present in those three serial summaries.

For every row, the V19 `REF_energy` equals the archived summary and output energy to within `1e-10` eV, and V19 `REF_forces` are bitwise identical to archived `PW_PBE_forces`. Atomic numbers, IDs, positions, cell, and PBC also match exactly. All source inputs and archived outputs are included in the repository paths listed in [reconciliation.json](reconciliation.json).

## What this resolves and what remains

- The 12 rows' source archives, numerical settings, energy getter, forces, and energies are now traceable. No label needs to be re-run for source recovery.
- Train36 is still fully periodic and Gamma-sampled, while dev3 uses TTF slabs with `6×6×1`, 0.10 eV smearing, and xy dipole correction. That training/development target mismatch remains a material confounder.
- The V12 distance-scan and V13 holdout rows are members of the current train36 set; they are not independent V19 development examples.
- The separate V17 rank-recovery row still has an unresolved source-versus-replacement k-point convergence question.
- The V19 candidate still fails the held-out force and separation-force gates and the Ti₃SiC₂ Ti/C pure-phase force screen. Reconciliation of existing labels does not make that model suitable for dynamics.
- The three dev labels still share COD 9009647 as a source family. This audit does not establish source-independent development coverage.

## Reproduction

From the repository root, run in the pinned MACE environment:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_inherited_train36_archive_reconciliation_20261023/reconcile_inherited_train36.py
```

The script reads only the approved V19 `train.extxyz` and the twelve exact input/output archive folders. It does not read `valid.extxyz`, test files, or sealed labels. SHA-256 values and per-row checks are in [reconciliation.json](reconciliation.json).
