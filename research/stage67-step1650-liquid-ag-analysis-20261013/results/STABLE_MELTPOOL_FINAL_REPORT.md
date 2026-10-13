# Stage67 step1650 stable meltpool: final liquid Ag extraction report

## Result

- `FULL_LIQUID_LIKE_AG_EXTRACTION = PASS`
- Reconstructed Liquid-like Ag candidates: **8,728**
- Largest connected candidate cluster: **1,032 atoms**
- Connected components at 3.74625 Å: **204**
- MD/LAMMPS execution: **none**
- Original dump modified: **no**

**Liquid-like Ag classification was reconstructed from available CellTl/PTM-related criteria in the archived input package.**

The 8,728 count is a reconstructed candidate population, not an independent per-atom label file supplied by the archive. The archive lacks original per-atom PTM/Q6 labels and a complete LiquidCandidate ID table. The classification here was regenerated from the original step1650 dump with the archived reproducible mask; it was not inferred from the 1,032-atom cluster.

## Source and method

- Input package commit: `9fa29e3e9dc8c120303f12bcaf8978b701b03604`; prior maximum-cluster report commit: `5fe8b480aaf144ba39e6a9388c8f247e6fadec5a`.
- Frame: step 1650, t = 0.165 ps; source dump SHA-256 `91d04a6fa3cb1b2fd9eb668b392cad04ecb593620723beac3e12ea4e7b8e85fc`.
- Model: Stage67 Ag–Ti₃SiC₂ transverse model; 303,372 atoms; cell 150.1 × 200 × 200 Å; `pp pp pp` periodicity.
- Ag type: 4. The dump includes `CellTl` and `CellTl_valid` per atom.
- Reconstructed mask: Ag type 4, OVITO PTM `Structure Type == 0` (`Other`), `CellTl_valid == 1`, and `CellTl >= 975 K`.
- PTM implementation: OVITO 3.16.1, RMSD cutoff 0.1, selected Ag centers only. The archived Q6 summary records Q6 threshold 0.446766; its accompanying extraction report states this value was summary metadata and the per-atom candidate mask used PTM Other. No threshold was changed or added.

The CellTl-valid/hot subset contains 8,732 Ag atoms: 4 are PTM FCC and excluded, and 8,728 are PTM Other candidates. This exactly reproduces the archived summary count of 8,728.

## Cluster and dimensions

PBC-aware connectivity at the archived cutoff 3.74625 Å produces 204 components. The maximum component has 1,032 Ag atoms, and its complete ID set exactly matches the archived 1,032-ID list. Its center of mass is (68.920100, 88.326150, 89.191421) Å; its X/Y/Z dimensions are (49.297781, 59.589429, 79.495963) Å. Distances to the nearest computational-box faces are (50.324403, 60.367842, 60.215459) Å; these are not physical free surfaces because all three directions are periodic. The audited maximum-cluster flags are `PBC_contact=NO` and `PBC_spanning=NO`.

The `meltpool_size.json` also records a minimum periodic coordinate envelope for all 8,728 candidates. This aggregate contains disconnected components and must not be read as the size of one connected melt pool.

## Integrity and scope

All four uploaded ZIP parts, the reconstructed ZIP, and the package SHA256 list passed. The original dump SHA, 303,372 atom count, step, and CellTl fields passed validation. The output dump preserves selected source atom records byte-for-byte and in source order, changing only the atom-count header to 8,728. IDs are listed in ascending order. The data represent the archived step1650 frame and supplied 0.155–0.165 ps stability window; the package does not include a full time-resolved trajectory or original per-atom classification table. No longer-time or thermodynamic stability claim is made.

See `MANIFEST.json` and `SHA256SUMS.txt` for artifact hashes and analysis environment.
