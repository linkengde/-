# Finite-interface model preparation report

Date: 2026-10-28  
Base: `33a1738f1fdf8f800233e6dbb65cfb12ee2b8c2c`  
New code directory: `phase_field_qsgs/finite_interface_prototype/`  
Scope: small conduction-kernel implementation and numerical tests only. The production Case A/B solver, checkpoints, source geometry, and double-temperature project were not changed. No long transient or new morphology was run.

## Model equation and conductance

For cell `i`, constant-capacity prototype equation is

\[
C_i V_i\frac{dT_i}{dt}=\sum_jG_{ij}(T_j-T_i)+q_i'''V_i+Q_{i,\mathrm{boundary}}.
\]

For two unlike-phase cells sharing a face of area `A`, with center-to-interface distances `δᵢ/2` and `δⱼ/2`, the series conductance is

\[
G_{ij}=\frac{A}{\delta_i/(2k_i)+1/G_K+\delta_j/(2k_j)},\qquad
Q_{i\leftarrow j}=G_{ij}(T_j-T_i).
\]

The same `Q` is subtracted from cell `j`, so interfacial exchange is pairwise conservative. Same-phase faces use the bulk half-cell series resistance and have no `1/G_K` term. The local temperatures on either side of an unlike-phase face are independent; finite `G_K` permits a temperature jump. This is a phase-to-phase continuum interface law, not electron–phonon two-temperature coupling.

Limits:

- `G_K = 0`: no interphase heat exchange, while each phase continues to conduct and store heat internally.
- `G_K → ∞`: no added interface resistance; the face reduces to the usual harmonic-mean perfect-contact conductance.
- These are not identical to the historical Case A/B definitions: historical Case A also sets the scaffold's `k` and `C` to zero, and historical Case B uses a single shared temperature. The new kernel keeps finite candidate Ti₃SiC₂ capacity/conductivity in both interface limits.

## Implementation

`finite_interface_conduction.py` accepts a three-dimensional phase-label field, orthogonal cell pitches, cell conductivity/capacity fields, `G_K`, optional Dirichlet exterior faces, and optional periodic axes. It precomputes conservative face links, has an explicit timestep bound, an energy-accounted step, and a sparse steady-state solver. It is deliberately isolated from `qsgs_enthalpy_solver.py`.

The current candidate 300 K inputs loaded by the tests are:

| Phase | Candidate `k` | Candidate volumetric sensible capacity |
|---|---:|---:|
| Ag | 429 W/(m·K) | 2.47011×10⁶ J/(m³·K) |
| Ti₃SiC₂ | 37 W/(m·K) | 2.54052×10⁶ J/(m³·K) |

The capacity values use candidate density × candidate mass-specific Cp. Candidate source status and limitations are detailed in `MATERIAL_SOURCE_AUDIT.md`. The synthetic `G_K=2×10⁸ W/(m²·K)` in tests is only a numerical fixture; it is not a measured or recommended Ag/Ti₃SiC₂ interface value.

## Validation completed

Command:

```bash
python -m unittest -v models.ag_ti3sic2_qsgs_1to1_finer.phase_field_qsgs.finite_interface_prototype.test_finite_interface_conduction
```

Five tests pass:

1. A finite interface face equals the series combination of Ag half-cell resistance, `1/G_K`, and Ti₃SiC₂ half-cell resistance.
2. One-dimensional bilayers with 4, 8, and 16 cells per phase match the analytic series-resistance heat rate for zero, finite, and infinite interface conductance at each resolution.
3. An insulated transient using the candidate Ag/Ti₃SiC₂ conductivities and capacities preserves total sensible energy to the test tolerance; internal interface power sums to zero.
4. At zero interface conductance, a phase-uniform initial state remains separated and unchanged.
5. The explicit update rejects a step above its conservative row-sum stability bound.

This verifies the linear conduction kernel and its energy accounting for small tests. It does **not** validate candidate properties, a real interface conductance, the 3D QSGS model, or melting behavior.

## Remaining implementation gates

Before connecting this kernel to the production enthalpy solver:

1. Resolve or explicitly sensitivity-bound Ag/Ti₃SiC₂ Cp, `k(T)`, density/porosity, and Ag phase-change data; the current values remain candidates.
2. Find a measured/literature range for `G_K` for the actual Ag/Ti₃SiC₂ interface preparation, or report a dimensionless sensitivity study without selecting a preferred physical value.
3. Add temperature-dependent phase properties and Ag enthalpy/latent heat while keeping the unlike-phase face exchange equal and opposite. Verify the `G_K→0` and `G_K→∞` limits against independent reference solutions after coupling.
4. Test the 1D bilayer against the closed-form steady flux at multiple cell counts and confirm convergence of interface temperature jump.
5. Validate external electrode/sink boundaries and energy removed to them. Only then run a small QSGS patch, after code and properties are frozen; no high-cost run is authorized by this prototype task.

