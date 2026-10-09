# CIF-derived Ag/interface normal-separation preparation

This is a geometry-only preparation package for A to review. It contains six fixed-ion proposals: move the complete four-atom Ag layer by −0.20 Å toward the framework or +0.20 Å outward for each AgC, AgSi and AgTi termination. The script derives the outward normal from the surface lattice vectors, preserves the cell, PBC, framework coordinates, atom order, IDs and pair marker, and rejects Ag/framework contacts below 0.80 times the sum of the ASE covalent radii. That cutoff screens severe overlaps; it is not a bonding or stability criterion.

The source inputs are the three verified files in `periodic_interface_v4/v19_cif_parent_structures/inputs/`. The generator rechecks the full source package `SHA256SUMS.txt`. All structures have 52 atoms, formula C16Ag4Si8Ti24, four Ag atoms and CIF hash `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`. The AgC/AgSi/AgTi source input hashes and every proposal hash, ordered-identity digest, translation vector, pair distance, gap and vacuum are recorded in `input_manifest.json`. The six extxyz files contain no energies or forces; all are `launch_enabled=false`, `owner=null`, with role `train_or_numerical_diagnostic_only`. No common DFT method is assigned pending A's mesh, boundary and energy-target review.

The geometry audit produced these contact-plane gaps and nearest Ag/framework distances:

| Interface | Parent gap (Å) | Inward proposal gap / nearest contact (Å) | Outward proposal gap / nearest contact (Å) | Minimum screened distance (Å) |
|---|---:|---:|---:|---:|
| AgC | 2.40 | 2.20 / 2.337 | 2.60 / 2.717 | 2.337 |
| AgSi | 2.60 | 2.40 / 2.593 | 2.80 / 2.967 | 2.593 |
| AgTi | 2.60 | 2.40 / 2.807 | 2.80 / 3.156 | 2.807 |

The geometry-only input audit scanned 17 public `input_manifest.json` files under `periodic_interface_v4`, excluding any manifest path containing `holdout`, `blind`, `sealed` or `test`. It opened and hashed manifest-declared input files and read summary/verification metadata only to confirm existing public archives. It found 11 verified AgSi numerical-pilot archives using exactly the AgSi parent input bytes: Gamma and 2×2×1, 3×3×1, 4×4×1, 5×5×1, 6×6×1 and 8×8×1 at 0.10 eV smearing, plus 5×5×1 and 6×6×1 at 0.05 eV and 0.20 eV. The source-geometry ledger also finds two converged 6×6×1 AgSi dipole controls, but those use PBC `[true,true,false]` whereas the parent proposal uses full PBC; treat them as boundary controls, not interchangeable labels for a full-PBC target. No exact parent input/archive was found for AgC or AgTi in this scan, and none of the six offset proposal hashes matches an existing input. The ledger records per-label archive checks, methods, input hashes and roles. Parent results do not provide labels for the proposed offsets. No dataset manifests, training-label files or sealed paths were read.

## Joint energy and force analysis for a future registered batch

For each termination, keep the parent (`s=0`) and the two rigid translations (`s=±0.20 Å`) at identical cell, PBC, k mesh, cutoff, smearing, PAW setups, convergence thresholds and boundary treatment. Save GPAW's native extrapolated energy and force-consistent free energy separately. Define the parent-referenced changes as

`ΔE_native(s) = E_native(s) − E_native(0)` and `ΔF(s) = F_free(s) − F_free(0)`.

Report both per cell and divided once by the projected interface area `|a × b|` in eV/Å². The cell has one Ag/framework contact plane. Do not compare raw total-energy zeros between terminations. An absolute interface/adhesion energy needs a matched separated-framework and Ag reference with the same atom counts, in-plane cell, PBC/vacuum and numerical method.

For the outward unit normal `n`, the separation-conjugate force on the rigid Ag subsystem is `G_Ag(s) = Σ(i in Ag) F_i · n`. Use GPAW forces from the free-energy calculator convention for force consistency. Also report the marked-pair diagnostic `G_pair = 0.5 (F_Ag − F_X) · n`, where the persistent pair marker selects Ag and the termination atom, and retain the full force vector for every atom in ID order. The pair projection is a local diagnostic; `G_Ag` is conjugate to translating all four Ag atoms.

At the finite step `h=0.20 Å`, compare the free-energy central secant `−[F_free(+h) − F_free(−h)]/(2h)` with `[G_Ag(+h)+G_Ag(−h)]/2`; report their difference and treat it as a finite-step consistency check with an O(h²) truncation term. Report the ± energy asymmetry `ΔE(+h)−ΔE(−h)` for native and free energies. Do not compare native-energy derivatives directly with forces at finite smearing. Include marked-pair distance change, net Ag force, pair force and all-atom force-vector RMS/max changes for each gap. These three related terminations share one CIF parent and cannot serve as three independent validation families.

## Reproduction and checks

From this directory, run `python prepare.py` to regenerate the six label-free geometries and geometry manifest, then `python audit_reuse.py` to rebuild the public-input reuse ledger. These scripts create no calculator and launch no DFT, inference, training, MD or heating. `SHA256SUMS.txt` inventories the final package.

The uploaded CIF was internally consistent and the builder reproduced its three parent structures byte-for-byte in the prior review. The live COD bytes were not independently fetched. A should review/freeze the surface registry, target settings and any required additional controls, then explicitly register any DFT jobs before execution.
