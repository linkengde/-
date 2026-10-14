# Stage E5 identifiable parameter response set

**Observed force-active columns:** `A_Ti-C, B_Ti-C, lambda1_Ti-C, lambda2_Ti-C, A_Ti-Si, B_Ti-Si`. They form a full-rank six-column matrix algebraically, but the unit-column condition number is 1104.7; singular directions are weak and several radial columns are nearly collinear. In particular, do not simultaneously free the four Ti–C radial fields or both Ti–Si amplitude fields from this dataset.

**Unexcited columns:** `A_Si-C, B_Si-C` have exactly zero response in every tested atom/component. The minimum periodic Si–C distances are all outside the Tersoff Si–C taper, whose outer limit is `R+D=3.211505 Å`:

| Structure | Minimum Si–C distance (Å) | Outer taper (Å) | Within taper |
|---|---:|---:|---|
| bulk_baseline_k10x10x4 | 3.5965 | 3.2115 | False |
| bulk_Ti4f_z_m0p02 | 3.5965 | 3.2115 | False |
| bulk_Ti4f_z_p0p02 | 3.5965 | 3.2115 | False |
| bulk_C4f_z_m0p02 | 3.5965 | 3.2115 | False |
| bulk_C4f_z_p0p02 | 3.5791 | 3.2115 | False |
| interface_derived_TSC_only_114 | 3.4982 | 3.2115 | False |

Therefore the present force labels cannot constrain `A_Si-C` or `B_Si-C`; those pair blocks need geometries that activate Si–C interactions before any fit involving them. This zero response is a structural coverage gap, not evidence that the parameters are physically irrelevant.

Among two-column subsets of the six active responses, the best-conditioned pair in this local matrix is `A_Ti-C, B_Ti-Si` (unit-column condition number 3.288). Treat this only as a candidate pair for a future, explicitly approved exploratory test: data-family independence, correct equilibrium/elastic response and validation roles are still missing. It does not establish unique parameter values or generalization.
