# Ag 70 vol% no-melt response: matched-state audit

## Scope and result

This is a read-only comparison of saved Case A states at 0.2 ps. It initializes the source and diffusion operators to inspect their spatial maps, but performs no transient time steps and does not change any source or material parameter.

The previously used Ag 41.6654 vol% QSGS realization has a 1251.41 K peak at 0.2 ps, just above the candidate Ag melting reference of 1235.08 K. It contains a small liquid region: Ag liquid fraction 0.00046354, 97 thresholded six-face voxels, and volume 2783.11 Å³. The Ag 69.9996 vol% geometry peaks at 748.61 K at the same time and has no liquid. Their whole-domain Ag mean temperatures are close: 311.39 and 312.35 K.

This large difference in local maximum is not explained by weaker total heating of the 70 vol% case. Its integrated source power is higher. The saved metadata confirms both cases use the same 100×100×50 grid and physical box, Case A boundaries, source center, 5 eV/(ps·Å²) peak flux, 30 Å Gaussian radius, 20 Å depth scale, exact binary-Ag source integration, and identical candidate property table/configuration/model hashes. The QSGS geometries differ: the old geometry SHA-256 begins `7e0e1509`; the Ag 70 vol% geometry SHA-256 begins `03ef47b1`.

| Matched 0.2 ps quantity | Ag 41.6654 vol% | Ag 69.9996 vol% | Change at 70 vol% |
|---|---:|---:|---:|
| Integrated source power | 0.886218 mW | 1.551359 mW | +75.05% |
| Injected energy by 0.2 ps | 1.77244×10⁻¹⁶ J | 3.10272×10⁻¹⁶ J | +75.05% |
| Case A Ag sensible capacity at 300 K | 1.47646×10⁻¹⁷ J/K | 2.48051×10⁻¹⁷ J/K | +68.01% |
| Source power / total Ag capacity | 60.02 K/ps | 62.54 K/ps | +4.20% |
| Mean Ag temperature | 311.39 K | 312.35 K | +0.95 K |
| Maximum temperature | 1251.41 K | 748.61 K | −502.80 K |
| Ag liquid volume fraction | 4.6354×10⁻⁴ | 0 | old case only |

The source and total Ag capacity therefore rise at nearly the same rate. This agrees with the nearly equal bulk Ag mean temperatures, while leaving room for a large difference in how heat concentrates locally.

## Where the source lands and how the Ag network carries heat

The integrated source has similar lateral placement in the two geometries: its power-weighted mean radius is 27.03 Å for Ag 41.7 vol% and 26.75 Å for Ag 70 vol%. It is slightly shallower in the old geometry: power-weighted mean depth is 17.04 Å versus 18.47 Å, and 71.79% versus 66.84% of source power lies within the top 20 Å.

Within the near-source region defined as radius <30 Å and depth <20 Å, the Ag voxel fraction is 38.71% in the old realization and 70.94% in the new one. Candidate Ag capacity at 300 K in this region is 5.2020×10⁻²⁰ versus 9.5323×10⁻²⁰ J/K, while source power there is 0.3963 versus 0.6622 mW. The local source-power/capacity ratio is consequently about 8.8% lower in the 70 vol% geometry (7619 versus 6947 K/ps). The maximum individual cell source power and maximum initial cell heating rate are the same in both geometries; the difference is the surrounding storage and transport network, not an altered source peak.

Case A sets both Ti₃SiC₂ conductivity and heat capacity to zero. The solver's conductivity and capacity fields are proportional to the Ag mask, and its harmonic face conductance is zero across an Ag/scaffold face. Heat can therefore travel only across shared Ag faces. A simple voxel adjacency diagnostic gives:

| Ag network proxy | Ag 41.6654 vol% | Ag 69.9996 vol% |
|---|---:|---:|
| Mean Ag face-neighbors per Ag voxel | 4.640 | 5.240 |
| Source-power-weighted mean face-neighbors | 4.417 | 5.052 |
| Ag face pairs along x/y/z | 162,653 / 162,839 / 157,788 | 308,235 / 308,255 / 300,458 |

These counts support a broader set of Ag heat-flow paths in the 70 vol% realization, but they are not an effective-conductivity calculation. The particular hottest voxel in the saved 0.2 ps fields has three Ag face-neighbors in the old geometry and two in the new one, so a single-cell neighbor count cannot explain the result. The multi-voxel path and the local temperature-dependent conductivity both matter.

At 300 K, the candidate Ag values are density 10,490 kg/m³, heat capacity 235.47 J/(kg·K), and conductivity 429 W/(m·K), giving α≈1.74×10⁻⁴ m²/s. The heuristic `sqrt(2 α t)` length is about 83 Å at 0.2 ps, comparable to the source diameter and larger than the 30 Å Gaussian radius. This makes the continuity and bottlenecks of the Ag paths especially relevant on this time scale. These are candidate values pending literature-source verification.

The snapshot diffusion operator also gives stronger instantaneous removal at the 70 vol% hottest cell: its net diffusion term is −0.735 μW against a 0.783 μW local source, compared with −0.625 μW against 0.812 μW in the old case. This is consistent with heat being carried away more effectively in the new snapshot, but it is a diagnostic of the already-heated fields, not a separate causal transient. The old hotspot is already liquid, where the candidate conductivity is 175 W/(m·K); the new 749 K hotspot remains solid, where the candidate conductivity is about 307 W/(m·K). Temperature-dependent conductivity may reinforce localization after a hotspot becomes very hot, but this comparison cannot establish how much of the original divergence it caused.

## Interpretation and limits

The strongest supported explanation is that the changed QSGS realization places a different Ag fraction and different face-connected paths around the source. In Case A, the extra Ti₃SiC₂ in the 41.7 vol% structure does not act as a heat sink at all; its zero conductivity and zero heat capacity block transport across the scaffold. The smaller Ag fraction thus also means less active heat capacity and a more restricted Ag-only conduction network. The 70 vol% geometry absorbs more total power, but has proportionally more Ag capacity and more broadly connected conducting paths, so the domain-average warming is similar and its local peak is much lower.

This is a comparison of two distinct QSGS realizations with different compositions. It does not isolate composition from morphology, prove that 70 vol% Ag is intrinsically more resistant to melting, or support an anti-welding conclusion. The old 0.2 ps liquid is also very small and cannot be called a stable melt pool. The 70 vol% 1× trajectory remains below melting through 1.0 ps (906.21 K maximum); its 2× trajectory has only progress samples through 0.38 ps, with full-field checkpoints through 0.35 ps, and must not be extrapolated to 1.0 ps.

No additional transient was run for this audit. The 2× calculation remains paused; after the cloud environment recovered, no process was present. Its 0.35 ps checkpoint was revalidated and committed as part of `16416e3742464d5991ca6de29b0e894b2de7fc91`. The comparison figure is `ag41_vs_ag70_heat_response_0p2ps.png`; the reproducible static diagnostics are `ag41_vs_ag70_static_diagnostics.json`, generated by `audit_composition_heat_response.py`.

The candidate properties remain unverified against primary literature, the scaffold/contact model is an idealized limit, and the outer domain is insulated. Before using this contrast in a paper, first source the Ag/Ti₃SiC₂ properties and constrain the physical heat-source footprint and duration. A fair micro/nano comparison should then hold composition, source, material properties, boundary assumptions, and voxel-resolution criteria fixed while varying a measured or justified morphology length scale across multiple realizations.
