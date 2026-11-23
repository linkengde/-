# Ag–W local interface reconstruction v3 audit

Date: 2026-11-23
Overall conclusion: **the Dai EFS potential produced statically relaxed local interfaces, but it did not pass the preregistered quantitative DFT screening for either orientation.** The project scope allows qualitative graduate-level exploration with this known bias, provided results are labeled exploratory and no quantitative interface-property claims are made from this comparison alone.

## Construction and run coverage

- The frozen protocol is in `V3_PROTOCOL.md`; its screening limits were written before the v3 EFS relaxations and were not changed afterward.
- The match search used the validated potential lattice constants (Ag 4.09002193 Å, W 3.15999679 Å), at most 1.00% absolute phase strain per direction, and at most 1,800 Å² cell area. The selected match is the smallest-area candidate under those limits.
- Ag(100)/W(100): Ag repeats 10×10, W 9×9, interface area 1,645.01 Å²; Ag strain −0.835%, W +0.842% in both directions; 1,400 Ag + 1,134 W atoms in seven planes per phase.
- Ag(111)/W(110): Ag repeats 4×5, W 9×9, area 1,151.39 Å²; Ag strains +0.182%/−0.835%, W −0.182%/+0.842%; 1,120 Ag + 1,134 W atoms in seven planes per phase.
- One interface was modeled per slab, with 15 Å vacuum. The central three atomic planes in each phase were fixed; all other atoms relaxed under the Dai EFS potential, without temporary springs.
- All 2×2 distinct layer-sequence terminations were scanned for Ag(100)/W(100), and all 3×2 for Ag(111)/W(110). Each termination used 3 or 4 coarse lateral translation samples and ten initial gaps (2.0–3.6 Å in 0.2 Å steps, plus 6.0 Å). Total: 360 interface minimizations.
- All 360 minimizations completed and met the 0.01 eV/Å minimizer tolerance. The largest maximum mobile force in the complete scan was 0.00152 eV/Å, below the source paper's 0.03 eV/Å force tolerance. No atoms were lost. Across all 360 relaxed points, the minimum Ag–Ag, W–W, and Ag–W distances were 2.484, 2.612, and 2.544 Å; no same-phase collapse floor was crossed. Thermodynamic closure was within 1.2×10⁻¹⁴ J/m².

## Best-energy EFS results against He et al. Table 4

| Orientation | EFS Wsep / γint (J/m²) | He et al. Wsep / γint (J/m²) | Difference (J/m²) | Strict DFT screen |
|---|---:|---:|---:|---|
| Ag(100)/W(100) | 2.685 / 1.021 | 2.80 / 0.72 | −0.115 / +0.301 | Fail: γint is 0.0009 J/m² outside the frozen ±0.30 window |
| Ag(111)/W(110) | 3.287 / 0.033 | 2.45 / 0.86 | +0.837 / −0.827 | Fail: both energy measures disagree substantially |

The selected minimum is at initial gap 2.0 Å for Ag(100)/W(100) and 2.8 Å for Ag(111)/W(110). The registries are nearly degenerate at the reported precision; this coarse translation grid does not establish that all atomic-scale registries are equivalent. Across termination choices, Ag(100)/W(100) gives Wsep 2.648–2.685 and γint 1.021–1.058 J/m². Ag(111)/W(110) gives approximately 3.287 and 0.033 J/m² for every sampled termination.

The DFT paper reports approximately 3.44% mismatch, seven layers, and 15 Å vacuum, but does not publish the complete interface cell vectors, strain partition, atomic coordinates, or DFT input. V3 instead uses a lower-strain EFS coincidence cell. The comparison therefore identifies a significant model discrepancy, especially for Ag(111)/W(110), but is not an exact matched-cell DFT reproduction.

## Use decision

- **Exploratory graduate-level atomistic work:** mechanically usable as a candidate after a valid 3D initial structure is built and passes its own density, topology, and static-force checks. Preserve the known potential limitations in reports and thesis methods.
- **Quantitative interface adhesion/energy, diffusion barriers, phase stability, and melting predictions:** not validated by this EFS/DFT comparison. The Ag(111)/W(110) discrepancy is too large for those claims.
- **Full QSGS model:** v1/v2 have not produced a valid 3D initial geometry. V3 does not solve that construction problem or authorize starting dynamics on the raw model.
- **TTM/laser heating:** still blocked by missing W electronic-property inputs, an Ag/W coupling closure, pulse definition, and a valid 3D structure. The old Ag-only 9,784.984 eV/ps preflight remains tied to the original overlapping coordinates.

The full per-gap energies and forces for all ten termination combinations are in `termination_runs/`; the aggregate machine-readable result is `local_interface_v3_report.json`, and the chosen low-strain cells are in `cell_match_report.json`. Relaxed best-gap structures are retained in `relaxed_structures/`.
