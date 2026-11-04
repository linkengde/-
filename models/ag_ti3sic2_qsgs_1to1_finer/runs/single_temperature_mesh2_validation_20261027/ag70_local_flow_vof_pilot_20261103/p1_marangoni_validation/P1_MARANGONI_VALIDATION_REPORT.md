# P1 Marangoni numerical validation

## Decision scope

This report validates an isolated, dimensionless Marangoni force implementation. It does not simulate the Ag70 geometry, solve an energy equation, predict Ag melt-pool velocity, or include the Ti₃SiC₂ skeleton. The P0 automatic-curvature path was used; the known exact-curvature/wavelet-AMR branch was not called.

The installed Basilisk source mirror did not contain `src/test/marangoni.c`, and the installed build did not record its Darcs revision. Therefore this is not claimed as a code-to-code reproduction of [Basilisk's official Marangoni test](https://basilisk.fr/src/test/marangoni.c). Instead, the planar-flow test below has a direct Stokes solution, and a separate curved-interface test checks normal/tangential force closure. The exact executable and key headers are identified by SHA256 in `results/p1_solver_configuration.json`.

## Numerical setup

The transient benchmark is a flat liquid/gas interface at y=0, with equal density and viscosity in the two phases, periodic x boundaries, no-slip top and bottom walls, and `stokes = true`. The flow domain has height 2H with each fluid layer of height H=1. A constant reference surface tension σ₀=0.01 uses Basilisk `tension.h` and its automatic height-function curvature. The imposed tangential surface force per unit area is

\[
\boldsymbol{\tau}_s = \frac{d\sigma}{dT}\nabla_s T,
\qquad \mathbf f_t = \boldsymbol{\tau}_s |\nabla f|,
\]

implemented as an acceleration on the centered solver's face field. The benchmark sets dimensionless ∂T/∂x=+0.01. The `neg` and `pos` cases use equal-magnitude dimensionless coefficients ±0.5; separate `agred` and `agoxygen` cases use literature coefficient ratios described below. The constant normal surface-tension term remains supplied by `tension.h`.

For equal layer heights and viscosities, integrating the steady Stokes equation across the interface gives

\[
U_i = \frac{\tau_s H}{\mu_l+\mu_g} = \frac{\tau_s H}{2\mu}.
\]

Thus Uᵢ = ±0.0025 for the symmetric `neg`/`pos` cases. A negative dσ/dT with temperature increasing toward +x produces uₓ<0 (hot-to-cold). A positive coefficient reverses the flow.

The grid matrix uses H/Δx = 16, 32, and 64 (N=32², 64², and 128² cells). Each run advances to dimensionless t=5. The timestep check compares capillary timestep factors 0.5 and 0.25 at H/Δx=32. No exact curvature was supplied.

## Results

| H/Δx | Case | Interface velocity | Stokes value | Relative error | Mass drift |
|---:|---|---:|---:|---:|---:|
| 16 | negative coefficient | -0.002421876 | -0.002500000 | 3.1250% | -1.48×10⁻⁹% |
| 32 | negative coefficient | -0.002460949 | -0.002500000 | 1.5620% | -8.46×10⁻¹⁰% |
| 64 | negative coefficient | -0.002480488 | -0.002500000 | 0.7805% | -2.13×10⁻¹⁰% |
| 16 | positive coefficient | +0.002421876 | +0.002500000 | 3.1250% | +1.53×10⁻⁹% |
| 32 | positive coefficient | +0.002460949 | +0.002500000 | 1.5620% | +8.59×10⁻¹⁰% |
| 64 | positive coefficient | +0.002480488 | +0.002500000 | 0.7805% | +2.20×10⁻¹⁰% |

The zero-gradient control stayed exactly at zero velocity on all three grids. The velocity error falls by approximately a factor of two per grid refinement and is below 1% at H/Δx=64. This is a consistent first-order trend for the interface-force/velocity sampling used here. Maximum absolute mass drift across all runs was 1.53×10⁻⁹%, well below 0.1%.

The full transient matrix contains 16 runs: zero-gradient, equal-magnitude positive/negative slopes, the two atmosphere-specific Ag slope signs, and a timestep-halving pair. The atmosphere-specific runs give the same convergent errors while preserving their expected signs. For the reducing-gas Ag slope, the predicted dimensionless velocity is -0.00136447 and measured values are -0.00132183, -0.00134316, and -0.00135382 at H/Δx=16, 32, and 64. For pO₂=100 Pa, the predicted values are +0.0025 and measured values are +0.00242188, +0.00246095, and +0.00248049.

For the H/Δx=32 negative-coefficient case, halving the timestep changed the interface speed by 0.000369%. The largest reported pressure Poisson residual was 1.65×10⁻⁹ (the solver took one iteration); the raw residual is preserved in the CSV and is not presented as a normalized, cross-solver residual criterion.

The P0 R64 static spherical-drop parasitic velocity was 8.62×10⁻⁵ maximum and 8.11×10⁻⁶ RMS. The P1 H/Δx=64 imposed flow was 0.0024805, about 28.8 times the P0 maximum and 306 times the P0 RMS. These are different benchmark geometries, so the comparison is an order-of-magnitude signal check, not a universal error bound. At the coarsest grid the ratio to P0's R16 maximum is smaller; the fine-grid result is the relevant accepted comparison.

## Normal and tangential surface-force closure

A separate static circle has radius R=0.5 and surface tension σ(θ)=σ₀+σ₁cosθ, with σ₀=1 and σ₁=0.1. The interface is reconstructed with PLIC facets. The normal force uses Basilisk's automatic height-function curvature; the tangential force uses the projected surface gradient. The exact integrated resultants are Fₙ,x=-πσ₁ and Fₜ,x=+πσ₁, so the total resultant is zero.

| R/Δx | Normal resultant error | Tangential resultant error | Force-closure residual |
|---:|---:|---:|---:|
| 16 | 0.2248% | 0.0137% | 0.2107% |
| 32 | 0.0524% | 0.0104% | 0.0419% |
| 64 | 0.0214% | 0.0077% | 0.0137% |

Both force components converge to their analytic values and their residual resultant decreases under refinement. This verifies the decomposition and sign convention in a static curved-interface force balance. It does not validate a dynamic thermocapillary droplet migration case.

## Silver coefficient sensitivity and limits

The supplied literature reports a change in Ag surface-tension slope with atmosphere. Morohoshi et al. (2008) report -0.298 mN/(m·K) for Ar/He with 5% H₂ and 2.6 ppm H₂O, over 1219–1734 K while the oxygen partial pressure varies with temperature; their fit is σ=984-0.298(T-1234) mN/m. At pO₂=100 Pa, they report +0.546 mN/(m·K) over 1350–1520 K with fit σ=698+0.546(T-1234) mN/m. Hibiya and Ozawa (2013) likewise report negative and positive regimes (e.g. -0.245 mN/(m·K) at pO₂=10⁻¹⁰ Pa and +0.609 mN/(m·K) at 1 Pa).

The `agred` and `agoxygen` runs map the 2008 slopes to a common dimensionless reference coefficient, preserving the literature signs and ratio. They are sign/sensitivity runs only. They do not use a dimensional thermal gradient, Ag viscosity, or dimensional length and therefore do not produce an Ag velocity. The measured linear fits are not extended beyond their reported ranges, and no 2277 K surface tension is inferred. Input details are in `ag_surface_tension_reference_inputs.csv`.

## Reproduction

From this directory, with Basilisk qcc available:

```sh
python run_p1_matrix.py
python run_p1_closure.py
```

The first command writes velocity, mass, timestep, and planar stress-balance tables plus all run logs. The second writes the curved-interface normal/tangential force-closure table. Solver options, grid sizes, source hashes, and parameter normalization are recorded in `results/p1_solver_configuration.json`.

## Conclusion

The automatic-curvature VOF implementation passes the controlled planar Marangoni sign, speed, mass, timestep, and stress-balance checks at the tested resolutions. The curved static-circle normal/tangential force split also closes and converges. The official dynamic Basilisk example was not reproduced, and material-specific Ag dynamics remain unvalidated; those boundaries are reflected in the final decision file.
