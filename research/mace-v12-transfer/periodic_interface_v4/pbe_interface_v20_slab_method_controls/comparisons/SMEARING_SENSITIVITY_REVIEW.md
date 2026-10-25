# V20 slab smearing sensitivity review

Review date: 2026-10-25

## Scope and integrity

This review compares the restarted, fixed-geometry `6x6x1` controls at Fermi smearing `0.05 eV` and `0.20 eV` against the retained V19 `0.10 eV` baseline. The three records use the same 52-atom AgSi slab geometry, PBE, 500 eV plane-wave cutoff, PAW setups, SCF thresholds and mixer, and xy dipole correction. Both V20 archives passed `verify_result.py`; the comparison script also verified the archives and exact geometry identity before computing the differences.

These controls use the V20 parent with `c = 44.2058765 Å` and about 26 Å of vacuum. The approved V19 AgSi development geometry uses `c = 34.2058765 Å` and about 16 Å of vacuum. The results below establish sensitivity for the V20 control parent; exact-development-geometry checks are registered separately.

The reported deltas are absolute differences from the `0.10 eV` baseline. The provisional limits are 2 meV/atom for each energy convention, 0.01 eV/Å for force-vector RMS, and 0.02 eV/Å for the marked Ag–Si pair-force projection.

| Smearing | Native energy delta (meV/atom) | Free-energy delta (meV/atom) | Force-vector RMS (eV/Å) | Maximum atom force-vector difference (eV/Å) | Absolute marked-pair projection delta (eV/Å) | Budget result |
|---:|---:|---:|---:|---:|---:|---|
| 0.05 eV | 0.0364 | 5.6536 | 0.01516 | 0.02757 | 0.02804 | Native energy passes; free energy, force RMS, and pair projection exceed limits |
| 0.20 eV | 0.5890 | 22.0606 | 0.03039 | 0.07325 | 0.04895 | Native energy passes; free energy, force RMS, and pair projection exceed limits |

Both SCF calculations converged: 59 iterations at `0.05 eV` and 54 at `0.20 eV`. These are completed calculations; `SENSITIVITY_EXCEEDS_BUDGET` describes the comparison outcome, not SCF failure.

## Interpretation and current method choice

The native extrapolated energy alone would make both smearing changes look small: its difference remains below 0.6 meV/atom. The force-consistent free-energy differences and force changes are larger and exceed their provisional limits at both tested smearings. The `0.20 eV` setting has the larger deviations. Therefore, the native-energy comparison by itself does not establish smearing insensitivity for this slab.

For the current workflow, retain `0.10 eV` as the provisional reference setting because it is the existing baseline and neither restarted alternative meets all comparison limits. This is a workflow choice supported by the current matched comparison, not a claim that `0.10 eV` is universally optimal or fully converged with respect to every numerical choice. Keep all three controls out of training and validation data. No further DFT runs or dataset edits are included in this assignment.

## Machine-readable evidence

- `comparison_AgSi_slab_k6x6x1_sigma0p05_v20__AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19.json`
- `comparison_AgSi_slab_k6x6x1_sigma0p20_v20__AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19.json`

`SHA256SUMS.txt` records checksums for this review and both comparison records.
