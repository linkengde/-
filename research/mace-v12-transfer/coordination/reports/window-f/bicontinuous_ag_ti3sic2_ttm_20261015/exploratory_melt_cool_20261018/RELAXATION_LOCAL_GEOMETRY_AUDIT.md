# Follow-up local geometry audit of bounded FIRE relaxation

This audit reads the original source and its archived LAMMPS force dumps only. It does not modify the source structure or run another simulation. Distances use minimum-image Y/Z and nonperiodic X. See `MINIMIZATION_SHORT_CONTACTS.csv`, `FORCE_SPIKE_LOCAL_NEIGHBORS.csv`, and `CONNECTIVITY_CUTOFF_SENSITIVITY.csv` for atom-level records.

## Species-resolved displacement from the unrelaxed model

| phase | atoms | median Å | P95 Å | P99 Å | max Å | >1 Å | >3 Å |
|---|---:|---:|---:|---:|---:|---:|---:|
| Ti | 13576 | 0.360 | 0.866 | 1.260 | 2.328 | 366 | 0 |
| Si | 4768 | 0.369 | 0.959 | 1.574 | 3.265 | 212 | 4 |
| C | 8800 | 0.351 | 0.881 | 1.199 | 2.006 | 238 | 0 |
| Ag | 12342 | 1.297 | 3.517 | 4.218 | 5.523 | 7324 | 1444 |

## Interpretation

The TSC geometric graph remains one connected component at the archived 3.30 Å cutoff, and the principal Ag component contains 99.708% of Ag atoms at 3.05 Å. Connectivity alone does not exclude local bond distortion. The post-minimization Ti–C minimum is 1.729 Å, compared with 2.092 Å in the source; the source model has no Ti–C contacts below 2.0 Å, while the relaxed model has 1,668 such pairs including 20 below 1.8 Å. Post-minimization Ag–C and Ag–Ti contacts also contract to 1.906 and 2.027 Å. These values and the extensive displacement/reorganization are reported as a potential-induced geometry risk, not as a proven force-field error or a proof of a chemical bond. The calibration-v0 cross terms are explicitly marked unvalidated.

The minimizer's transient force spike at step 475 is localized by the neighbor table; it fell by the saved step 500. It is not itself proof of final instability, but together with the new short-contact population it warrants stopping before thermal driving under this mapping. A separate 300 K NVT screen was started on the relaxed copy and interrupted after 0.2 ps because a repeated localized interface force spike appeared; see `SCREEN_300K_AUDIT.md`. No heating/cooling, TTM, DFT, or training run was started.
