# QSGS enthalpy solver: first numerical pilot

## Status

The first three-dimensional finite-volume enthalpy solver is implemented and passes unit, conservation, and one-dimensional Stefan checks. Short QSGS runs have exercised both requested scaffold limits at 5 Å and 10 Å. This is a **single-temperature continuum prototype**: heat is added directly to Ag lattice enthalpy. It is not TTM/TTM-MD, not a stable-pool production run, and not a prediction of measured Ag–Ti₃SiC₂ interface heat transfer.

## Model choices

- The geometry comes from the conservative, unthresholded QSGS phase fractions. The 5 Å and 10 Å grids preserve the input Ag fraction (0.416654).
- The enthalpy state stores energy per total cell volume. Ag sensible heat, latent heat, and liquid sensible heat are inverted to temperature and Ag liquid fraction. In mixed cells, both phase fractions enter the capacity and transport model.
- Case A is the adiabatic-scaffold limit: skeleton cells have no heat capacity or conductivity, so temperature is only meaningful in Ag-containing cells. There is zero heat flow through the Ag/skeleton geometry in this homogenized mask approximation.
- Case B is the ideal-contact limit: Ag and Ti₃SiC₂ share a temperature and have volume-weighted heat capacity and conductivity, with no added interface resistance.
- x and y are periodic; the two z faces are insulated. There is no radiation, convection, melt flow, free-surface motion, or external heat sink.
- The source retains the specified `exp(-r²/R²) exp(-d/l_skin)` shape, R=30 Å, l_skin=20 Å, and peak flux 5 eV/(ps·Å²). It is multiplied by the Ag volume fraction. The resulting power is integrated from the unrenormalized masked grid.
- The supplied thermophysical anchors remain candidates. Exact source tables/citations are not yet archived. Density is fixed at the solid reference value; melt expansion is omitted. The Ti₃SiC₂ heat capacity conversion uses 195.708 g/mol, giving 562–792 J/(kg·K), not the 656–924 J/(kg·K) in the supplied summary. The Ti₃SiC₂ model stops at 1473.15 K rather than extrapolating.

The explicit diffusion step is set from the initial 300 K maximum k/Cp bound; k/Cp decreases over the implemented candidate temperature curves. The solver also guards direct time steps against the computed stability limit. Energy exchange through internal faces is pairwise conservative.

## Checks and pilot runs

The full repository test command `python -m unittest -v` passes 42 tests; 16 directly cover thermophysical functions, QSGS projection, and this enthalpy solver. Those tests include enthalpy/temperature inversion, latent-heat plateaus, diffusion conservation, source accounting, and explicit-step guarding. The 40 µm and 20 µm 1-D Stefan grids put the interface within 0.12% and 0.08% of the similarity solution; their energy residuals are below 1.6×10⁻⁹ J/m². These checks validate the implemented numerics for those cases, not the physical accuracy of the candidate material data.

At matched 0.1 ps heat plus 0.1 ps cooling, Ag-volume-weighted mean temperature agrees closely between grids (within 0.1 K), while local maxima differ by about 14% for case A. The grid-integrated masked power differs by about 2%:

| Resolution | Boundary | Masked source power | Max T at 0.1 ps | Ag-weighted mean T at 0.1 ps | Max T after 0.1 ps cooling | Energy residual after cooling |
|---|---|---:|---:|---:|---:|---:|
| 5 Å | A | 0.8808 mW | 834.3 K | 305.86 K | 463.0 K | −3.0×10⁻³⁰ J |
| 5 Å | B | 0.8808 mW | 517.9 K | 302.46 K | 358.5 K | +1.1×10⁻³⁰ J |
| 10 Å | A | 0.8639 mW | 720.0 K | 305.77 K | 391.5 K | −6.8×10⁻³¹ J |
| 10 Å | B | 0.8639 mW | 500.0 K | 302.39 K | 352.3 K | +2.2×10⁻³¹ J |

The short grid comparison does not melt Ag. Its large peak-to-mean difference confirms that the local peak needs finer-grid convergence; 10 Å is suitable only as a coarse screening grid for this source.

A separate 10 Å screening run heated for 5 ps and then cooled for 1 ps. Case A briefly reached 1477 K and a sampled Ag liquid-volume fraction of 0.00167 at the end of heating, then returned to zero liquid fraction within the 1 ps cooling window. Case B peaked at 768 K during heating and did not melt. Energy residuals were −1.5×10⁻²⁷ J (A) and +2.1×10⁻²⁹ J (B). The tiny, one-sample case-A liquid region is **not** a stable melt pool; the boundary sensitivity is large and neither case can be declared the physical interface condition.

CSV/JSON time histories for all pilots are in this directory (`qsgs_enthalpy_*`). `qsgs_enthalpy_solver.py` can reproduce them. For example:

```bash
python phase_field_qsgs/qsgs_enthalpy_solver.py --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0.1
python phase_field_qsgs/qsgs_enthalpy_solver.py --spacing 10A --boundary B --heat-ps 5 --cool-ps 1
```

## Decision and next work

The first 3-D solver passes its implementation-level energy checks. It has not established a stable pool or a credible freezing time. Before interpreting melting quantitatively, archive and audit the cited property tables, complete spatial peak convergence on the 5 Å mesh, and extend the cooling/heat-window comparisons. The two scaffold cases must remain separate bounds because their results differ strongly.

After the enthalpy/conduction baseline is accepted, connect the electronic energy equation and Ag electron–lattice coupling as a separate two-temperature model. That integration must preserve source power and total electron-plus-lattice energy before any 100 ps run or automatic stable-pool shutoff criterion is attempted.
