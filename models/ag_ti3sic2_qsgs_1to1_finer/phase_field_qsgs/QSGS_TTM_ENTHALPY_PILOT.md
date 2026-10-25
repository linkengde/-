# QSGS continuum two-temperature–enthalpy pilot

## Status

The 3D QSGS continuum solver now uses a conservative semi-implicit two-temperature/enthalpy update. The earlier concern about counting the full Ag thermal conductivity in both electron and lattice equations has been corrected. The updated 10 Å, 10 fs pilot remains fully solid; it is a numerical and energy-transfer check, not a melt-pool result or a production run.

This is a continuum calculation, not atomistic TTM-MD. The Ag–Ti₃SiC₂ interface is not modeled atomistically and Cases A/B remain thermal-boundary limits rather than measured interface properties.

## Conductivity partition and equations

- The independent single-temperature enthalpy model retains the Ag total conductivity table (about 429 W/(m·K) at 300 K).
- In the two-temperature model, Ag electrons use the candidate `ke(Te)` table. Ag lattice conductivity defaults to 0 W/(m·K) so the full total Ag conductivity is not counted a second time. This is an explicit first approximation; no unsupported `k_total - ke` subtraction is used.
- Case A treats the scaffold as adiabatic. Case B adds Ti₃SiC₂ sensible heat capacity and conductivity to the shared lattice-temperature field; only the scaffold fraction conducts in this lattice equation.

The electron equation is advanced in tabulated internal energy with a conservative finite-volume operator. A backward-Euler/Newton-CG solve treats electron diffusion; the tabulated `Ue(Te)` secant slope supplies its local energy capacity. Lattice diffusion is backward Euler in enthalpy with a Newton-GMRES solve, preserving the Ag latent-heat plateau. The conductivity is frozen within each diffusion half-step, so this is semi-implicit and still needs time-step convergence checks. Local electron–lattice exchange is solved cell-by-cell by a bracketed implicit root; each cell transfers equal and opposite energy. The heat source and diffusion are applied as balanced half-steps around that exchange.

The specified source is unchanged: Ag-masked \(\exp(-r^2/R^2)\exp(-d/l_{skin})\), with \(R=30\) Å, \(l_{skin}=20\) Å and peak flux \(I_0=5\) eV/(ps·Å²). The 10 Å QSGS grid integrates to 0.863925 mW; no post-mask renormalization is applied. The Ag electron parameter table spans 300–20,000 K and extrapolation remains disabled.

By default x/y are periodic and z is insulated. `--bottom-sink` adds an ideal 300 K Dirichlet bath at the bottom z face to both the Ag electron conduction path and any conducting lattice phase; energy removed by the bath is recorded as a negative diffusion contribution in the balance. This is a boundary-condition sensitivity case, not a measured substrate interface.

## Verification

The full model test suite passes **51 tests**. The nine double-temperature tests cover 0D energy exchange, 1D electron and lattice diffusion, conductivity partition, active-cell enthalpy inversion, QSGS source/coupling accounting, sampled interval accounting, diffusion steps above the explicit stability bounds, and bottom-sink energy removal. Electron and lattice diffusion conserve energy for insulated boundaries to solver tolerance. In synthetic sink tests, both channels remove energy and the external loss closes the total-energy ledger.

For the 10 Å grid, a 10 fs heat/cool run at 0.5 fs outer steps was compared with 0.3125 fs. The final peak electron temperature changed from about 339.24 K to 339.95 K (about 0.21%); cumulative electron-to-lattice transfer changed by about 0.11% for Case A. Case B is similarly close. This is a preliminary refinement check, not proof of convergence for melting or longer transients.

## 10 Å QSGS pilot: 5 fs heat, then 5 fs source-off cooling

| Boundary | Peak \(T_e\) at heat-off | Peak \(T_l\) at heat-off | Peak \(T_e\) after 10 fs | Peak \(T_l\) after 10 fs | Energy input | Final Ag liquid fraction | Total-energy residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| A: scaffold adiabatic | 805.759 K | 300.01985 K | 339.950 K | 300.02588 K | 4.31962×10⁻¹⁸ J | 0 | −1.99×10⁻³³ J |
| B: ideal thermal contact | 805.759 K | 300.01659 K | 339.950 K | 300.02095 K | 4.31962×10⁻¹⁸ J | 0 | −1.53×10⁻³³ J |

The total source-off electron-to-lattice transfer through 10 fs is 3.5773×10⁻²⁰ J in A and 3.5773×10⁻²⁰ J in B. Net insulated-boundary diffusion is at round-off scale. Neither case approaches Ag melting; no stable-pool controller or 100 ps run has been attempted.

The 0.5 fs Case A 10 fs pilot takes about 4 s in this Python implementation on the current environment. Linear scaling would put a 100 ps, 10 Å run at roughly 11 hours per boundary case, before adding the controller or refining to 5 Å. This rough estimate is not a production performance guarantee.

## Reproduction

Run from the model directory:

```bash
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary B --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
```

The histories and machine-readable summaries are `qsgs_ttm_enthalpy_10A_BCA.csv/.json` and `qsgs_ttm_enthalpy_10A_BCB.csv/.json`. To evaluate the ideal bottom sink, add `--bottom-sink`; sink and insulated results must be reported separately.

## Remaining gates

The solver implementation and boundary accounting now pass short numerical checks. Before interpreting melt-pool size or lifetime, the next gates are (1) refine time-step accuracy over longer heating intervals, (2) test the 5 Å grid against the 10 Å result, (3) assess candidate property-table provenance and high-temperature behavior, and (4) build and validate the spatially connected liquid-pool heat-off controller. These pilots do not establish physical accuracy of the candidate Ag parameters or real Ag–Ti₃SiC₂ heat transfer.
