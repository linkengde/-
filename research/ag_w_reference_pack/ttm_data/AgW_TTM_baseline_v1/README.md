# AgW_TTM_baseline_v1

This directory records the selected first exploratory parameter baseline for the Ag–W TTM work. Lin et al. (2008) supplies the temperature-dependent electronic heat capacity and electron–phonon coupling curves for both elements. W electronic conductivity uses the practical Daraszewicz et al. (2014) Drude closure. Ag electronic conductivity uses a Lin Drude-form `ke(Te,Tl)` closure anchored to the user-supplied 429 W m⁻¹ K⁻¹ room-temperature bulk value; the XTANT-3 K table is retained as a comparison and is not used as the room-temperature amplitude.

## Parameter choices

| Material/property | Baseline choice | Data/units | Main scope note |
|---|---|---|---|
| Ag `Ce(Te)` | Lin 2008 Fig. 2(c), digitized | `Ag_TTM_Lin2008_digitized_300_20000K.csv`, J m⁻³ K⁻¹ | Figure-derived candidate; reported `gamma_exp` used through 5000 K as described in its digitization report |
| Ag `G(Te)` | Lin 2008 Fig. 2(d), digitized | same CSV, W m⁻³ K⁻¹ | Low-temperature plateau and short bridge are explicit graph-based approximations |
| Ag `ke(Te,Tl)` | Lin 2008 Drude form with the Lin Ag `Ce(Te)` curve, room-anchored at 429 W m⁻¹ K⁻¹ | `Ag_ke_Lin_Drude_room_anchor.csv`, W m⁻¹ K⁻¹ | The A/B scattering ratio is inferred from XTANT's 300 K channel split; this is an exploratory closure, not an Ag-specific measured `ke(Te,Tl)` surface |
| W `Ce(Te)` | Lin 2008 Fig. 7(c), digitized | `../W_TTM_Lin2008_digitized_300_20000K.csv`, J m⁻³ K⁻¹ | Figure-derived candidate |
| W `G(Te)` | Lin 2008 Fig. 7(d), digitized | same W CSV, W m⁻³ K⁻¹ | Below 600 K is explicitly extrapolated |
| W `ke(Te,Tl)` | Daraszewicz 2014 Drude closure | formula in `baseline_parameters.json`, W m⁻¹ K⁻¹ | Equilibrium fit uses 300–2500 K; high-Te use is extrapolative |
| W `G0` | Daraszewicz measured and ab-initio values | low-temperature comparison only | Reference anchors, not a second curve or an average |

The raw XTANT Ag K table remains archived unchanged. Its `kappa_tot=16.56 W m⁻¹ K⁻¹` at 300 K is not used as the room-temperature Ag baseline magnitude: the total is limited by its reported electron-phonon channel (`16.66 W m⁻¹ K⁻¹`) despite an electron-electron channel of `2926.3 W m⁻¹ K⁻¹`. Instead, the Ag Drude closure is normalized to the user-provided room-temperature bulk value `429 W m⁻¹ K⁻¹`, treating Ag's total as an electronically-dominated proxy. Its A/B ratio is inferred from the XTANT 300 K channel ratio under the Lin inverse-scattering-time decomposition. This one-point calibration and high-Te behavior remain exploratory; the precise user citation for 429 W m⁻¹ K⁻¹ was not provided as a bibliographic URL.

The selected model does not rescale Lin's W `G(Te)` to match Daraszewicz's `G0`.

## Validation result

The low-temperature comparison finds Lin's digitized W `G(300 K)` above Daraszewicz's experimental value and outside its quoted 1σ interval. The Lin 300 K point is also in the documented extrapolation region below 600 K. This discrepancy remains visible in the baseline record.

A single-branch 0D W pilot checks the Lin `Ce/G` energy ledger with one normalized temporal Gaussian. It conserves energy to the reported numerical tolerance. This pilot does not exercise conductivity because a zero-dimensional model has no spatial gradient. The Ag room-anchor consistency is recorded by the conductivity table and validation JSON; no high-Te experimental conductivity validation is claimed.

## Spatial laser source and model limits

The inherited Gaussian spatial profile and Ag-only mask remain documented in `baseline_parameters.json`. Its temporal envelope and scenario absorbed energy still need to be set. The Ag–W grid prototype now uses a provisional face-based `h_e=1 GW m⁻² K⁻¹` closure between separate phase electron grids, based on Au/Pt and Au/Al method analogues. It transfers energy equally and oppositely, but it is not an Ag/W measurement and omits sub-cell interface area. An Ag lattice-conductivity curve remains open. The 3D free-surface TTM-MD model remains gated on a statically accepted Ag–W structure.

Reproduce the Ag digitization with `../../scripts/digitize_lin2008_ag_ttm.py` and regenerate the W energy pilot plus Ag conductivity grid with `../../scripts/validate_agw_ttm_baseline_v1.py`.

## Primary references

- Lin, Zhigilei & Celli (2008), *Physical Review B* 77, 075133, [doi:10.1103/PhysRevB.77.075133](https://doi.org/10.1103/PhysRevB.77.075133).
- Daraszewicz et al. (2014), *Applied Physics Letters* 105, 023112, [doi:10.1063/1.4890413](https://doi.org/10.1063/1.4890413).
- Tolias (2017), tungsten thermophysical properties, [arXiv:1703.06302](https://arxiv.org/abs/1703.06302).
