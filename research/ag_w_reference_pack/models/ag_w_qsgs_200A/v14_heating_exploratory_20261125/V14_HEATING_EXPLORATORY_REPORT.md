# V14 Ag–W heating: exploratory result

**Run status:** `AGW_TTM_SOFTWARE_LINK_PASS` for the numerical energy ledger. This is not a validated melt/quench prediction.

## Run setup

- V14 finer QSGS interpenetrating Ag/W geometry; 481,648 atoms, `p p f`, box 200 × 200 × 300 Å.
- Dai 2006 EFS candidate; fixed atoms at z=0–20 Å, 300 K Langevin heat sink at z=20–40 Å, Newton integration at z=40–300 Å, TTM map at z=40–200 Å.
- 10 Å TTM grid; timestep 0.00025 ps; 0.5 ps Gaussian pulse centered at 0.5 ps; 5 ps after the source envelope; provisional Ag/W electronic interface conductance 1×10⁸ W/(m²·K).
- **This run used `source_scale=100`, so effective I₀=500 eV/(ps·Å²), 100 times the nominal 5 eV/(ps·Å²).** The input was 449,264.3 eV. This is a high-dose exploratory run and does not represent the nominal source setting.
- The 0.5 ps pre-TTM NVT relaxation left mapped Ag/W near 196.4/195.4 K; the TTM run did not start at 300 K.

## Numerical coupling ledger

- Laser source energy: 449,264.312 eV.
- Energy removed by the 300 K bottom Langevin bath by 7.5 ps: 18,430.631 eV.
- Final combined residual: 1.736 eV (recorded tolerance 467.695 eV).
- Maximum grid-step residual: 1.573×10⁻¹⁰ eV.

The software energy ledger passes. At the end of the 5 ps source-off interval, however, the system was still very hot: mapped Newton-layer Ag/W temperatures were 2,478/2,707 K, the tracked Ag hot zone was 1,999 K, and peak electron temperatures were 7,785/6,011 K. The run therefore demonstrates source delivery and energy bookkeeping, but **not cooling to a solid state**.

## Melt-pool assessment

A bounded melt-pool size cannot be assigned from these snapshots.

At 3.0 ps, the Q6 + PTM + displacement (>2 Å) rule returned 130,742 Ag candidates (56.5% of mapped-layer Ag). The candidate fraction was nearly flat across the lateral radius: 57.7% at 20–30 Å from the beam center and 55.4% at 90–100 Å. At 7.5 ps it returned 178,407 candidates (77.1%), again nearly uniform laterally; candidate RMS displacement was 24.18 Å.

The order metrics already fail at the starting configuration: 230,307/231,493 mapped Ag atoms (99.5%) were PTM `OTHER`, and 176,275/231,493 (76.1%) had Voronoi Q6 <0.40. Thus PTM/Q6 do not distinguish a new liquid region from the disordered starting structure. In addition, the tracked hot-zone Ag displacement RMS reached 32.85 Å, with all tracked atoms displaced by more than 2 Å. These results indicate broad structural rearrangement under the 100× source, not a localized, validated liquid pool.

**Conclusion:** no defensible radius or volume for a localized melt pool is available from this run. Do not describe the 130,742 or 178,407 candidates as liquid atoms or a measured melt-pool size.

## Next run

Repeat from the same V14 300 K restart with the original nominal intensity (`source_scale=1`, effective I₀=5 eV/(ps·Å²)). Keep it as a separate run and compare spatial temperature and structural change to this high-dose record. The initial PTM/Q6 disorder must remain an explicit baseline limitation.
