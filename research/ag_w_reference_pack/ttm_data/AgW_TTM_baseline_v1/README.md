# AgW_TTM_baseline_v1

This directory records the selected first exploratory parameter baseline for the Ag–W TTM work. Lin et al. (2008) supplies the temperature-dependent electronic heat capacity and electron–phonon coupling curves for both elements. W electronic conductivity uses the practical Daraszewicz et al. (2014) Drude closure. Ag electronic conductivity uses the supplied Ag-specific XTANT-3 `kappa_tot(Te)` table because no Ag `ke(Te,Tl)` closure was found in the reviewed sources.

## Parameter choices

| Material/property | Baseline choice | Data/units | Main scope note |
|---|---|---|---|
| Ag `Ce(Te)` | Lin 2008 Fig. 2(c), digitized | `Ag_TTM_Lin2008_digitized_300_20000K.csv`, J m⁻³ K⁻¹ | Figure-derived candidate; reported `gamma_exp` used through 5000 K as described in its digitization report |
| Ag `G(Te)` | Lin 2008 Fig. 2(d), digitized | same CSV, W m⁻³ K⁻¹ | Low-temperature plateau and short bridge are explicit graph-based approximations |
| Ag `ke(Te)` | XTANT-3 Ag `K_Ag.dat` `kappa_tot` | W m⁻¹ K⁻¹ | Te-only and low-temperature source values carry the source's reliability caveat |
| W `Ce(Te)` | Lin 2008 Fig. 7(c), digitized | `../W_TTM_Lin2008_digitized_300_20000K.csv`, J m⁻³ K⁻¹ | Figure-derived candidate |
| W `G(Te)` | Lin 2008 Fig. 7(d), digitized | same W CSV, W m⁻³ K⁻¹ | Below 600 K is explicitly extrapolated |
| W `ke(Te,Tl)` | Daraszewicz 2014 Drude closure | formula in `baseline_parameters.json`, W m⁻¹ K⁻¹ | Equilibrium fit uses 300–2500 K; high-Te use is extrapolative |
| W `G0` | Daraszewicz measured and ab-initio values | low-temperature comparison only | Reference anchors, not a second curve or an average |

The source files and caveats are kept with the data. The selected model does not rescale Lin's W `G(Te)` to match Daraszewicz's `G0`.

## Validation result

The low-temperature comparison finds Lin's digitized W `G(300 K)` above Daraszewicz's experimental value and outside its quoted 1σ interval. The Lin 300 K point is also in the documented extrapolation region below 600 K. This discrepancy remains visible in the baseline record.

A single-branch 0D W pilot checks the Lin `Ce/G` energy ledger with one normalized temporal Gaussian. It conserves energy to the reported numerical tolerance. This pilot does not exercise conductivity because a zero-dimensional model has no spatial gradient.

## Spatial laser source and model limits

The inherited Gaussian spatial profile and Ag-only mask remain documented in `baseline_parameters.json`. Its temporal envelope and scenario absorbed energy still need to be set. The baseline also lacks an Ag lattice-conductivity curve, Ag/W electronic interface conductance, and a validated cross-phase electronic transport rule. The 3D free-surface TTM-MD model remains gated on a compatible solver and a statically accepted Ag–W structure.

Reproduce the data extraction with `../../scripts/digitize_lin2008_ag_ttm.py` and the validation with `../../scripts/validate_agw_ttm_baseline_v1.py`.

## Primary references

- Lin, Zhigilei & Celli (2008), *Physical Review B* 77, 075133, [doi:10.1103/PhysRevB.77.075133](https://doi.org/10.1103/PhysRevB.77.075133).
- Daraszewicz et al. (2014), *Applied Physics Letters* 105, 023112, [doi:10.1063/1.4890413](https://doi.org/10.1063/1.4890413).
- Tolias (2017), tungsten thermophysical properties, [arXiv:1703.06302](https://arxiv.org/abs/1703.06302).
