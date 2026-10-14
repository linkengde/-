# Parameter selection from the E5 response matrix

## Selection rules

Retain the E3/E5 tied-block mapping: one scalar field shared over the six table rows belonging to an unordered pair. Do not free individual triplet rows. Keep all Si–C fields frozen until relevant DFT-labelled geometries activate their cutoff. The numerical perturbation `±1%` is a response probe, not a fitted bound. No numeric fit interval is authorized by this design.

Candidate subsets below were chosen by enumerating the E5 six active columns, normalizing each selected column to unit norm, then minimizing the resulting submatrix condition number at its stated size. Rank and condition refer only to the current 1062 force-component rows and only to force responses.

| Candidate size | Best-conditioned tied scalar blocks | Current submatrix rank | Unit-column condition | Reading |
|---|---|---:|---:|---|
| 2 | `A_Ti-C`, `B_Ti-Si` | 2/2 | 3.288 | Best first screening subset on present geometries; not a released fit set. |
| 3 | `A_Ti-C`, `lambda2_Ti-C`, `B_Ti-Si` | 3/3 | 14.156 | More expressive; conditioning worsens and radial compensation remains. |
| 4 | `A_Ti-C`, `lambda2_Ti-C`, `A_Ti-Si`, `B_Ti-Si` | 4/4 | 46.718 | Algebraically full rank but substantially more sensitive to noise/model bias and overfit. |

Exact values and the subset search basis are recorded in `fit_candidate_conditions.csv`.

## Data requirements and risk by candidate size

### Two blocks

If a future C-approved exploratory fit is considered, start with at most `A_Ti-C` and `B_Ti-Si` while freezing every other field. Current force-only columns are comparatively separated, but the data still contain only one bulk parent family plus one interface-derived slab family; there are no independent DEV/TEST families and no stress labels for the bulk. The fit could overfit this family even with a modest condition number. It needs method-matched DFT force targets spanning bulk equilibrium and displaced/strained environments, a clean-surface training family, and fully separate DEV and TEST roots.

### Three blocks

The `A_Ti-C`, `lambda2_Ti-C`, `B_Ti-Si` subset adds a Ti–C radial decay degree of freedom. Amplitude/decay tradeoffs require a wider and physically relevant Ti–C distance distribution than the present small-displacement family. Add volume/shear responses and independent local displacements; use held-out families to test whether the extra decay parameter improves more than TRAIN residuals. Current normalized condition is 14.156, not a forecast for the expanded set.

### Four blocks

The `A_Ti-C`, `lambda2_Ti-C`, `A_Ti-Si`, `B_Ti-Si` subset adds a second Ti–Si amplitude direction. E5 already shows strong `A_Ti-Si` / `B_Ti-Si` anti-correlation. It requires multiple Ti–Si distance/coordination environments, Si displacement data, elastic response and clean surfaces. Current condition 46.718 signals markedly higher variance risk; do not use four blocks unless new-family Jacobians demonstrate robust separation and held-out performance. The current package does not establish that.

## Fit objective if a later gate passes

Use a family-balanced objective with (i) 3D atomic force vectors, first balanced across Ti/Si/C and then across source families; (ii) within-family, same-composition energy differences with explicit native/free-energy convention; (iii) displacement restoring-force response; (iv) DFT stress and bulk volume/shear response once available; (v) lattice/internal-coordinate and surface-response metrics; and (vi) regularization toward the verified source parameterization in `log(p/p0)`. Do not subtract absolute DFT and LAMMPS total energies with unrelated zeroes. Do not put slab 3D stress into a bulk stress target. Determine weights from DFT uncertainty and a frozen physical priority scheme, not atom count.

Do not specify numerical bounds or regularization strength from the ±1% probe. First independently confirm paper-source constraints, units and implementation conventions; then derive admissible bounds from the source constraints and physical checks. Fit/selection must never inspect TEST.
