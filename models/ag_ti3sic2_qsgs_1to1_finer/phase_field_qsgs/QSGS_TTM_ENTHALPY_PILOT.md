# QSGS continuum two-temperature–enthalpy pilot

## Status

The 3D QSGS continuum solver now uses a conservative semi-implicit two-temperature/enthalpy update. The earlier concern about counting the full Ag thermal conductivity in both electron and lattice equations has been corrected. The 10 Å, 0.1 ps, 1 ps and 5 ps heating screenings remain fully solid; they are heating-response checks, not melt-pool results or production runs.

This is a continuum calculation, not atomistic TTM-MD. The Ag–Ti₃SiC₂ interface is not modeled atomistically and Cases A/B remain thermal-boundary limits rather than measured interface properties.

## Conductivity partition and equations

- The independent single-temperature enthalpy model retains the Ag total conductivity table (about 429 W/(m·K) at 300 K).
- In the two-temperature model, Ag electrons use the candidate `ke(Te)` table. Ag lattice conductivity defaults to 0 W/(m·K) so the full total Ag conductivity is not counted a second time. This is an explicit first approximation; no unsupported `k_total - ke` subtraction is used.
- Case A treats the scaffold as adiabatic. Case B adds Ti₃SiC₂ sensible heat capacity and conductivity to the shared lattice-temperature field; only the scaffold fraction conducts in this lattice equation.

The electron equation is advanced in tabulated internal energy with a conservative finite-volume operator. A backward-Euler/Newton-CG solve treats electron diffusion; the tabulated `Ue(Te)` secant slope supplies its local energy capacity. Lattice diffusion is backward Euler in enthalpy with a Newton-GMRES solve, preserving the Ag latent-heat plateau. The conductivity is frozen within each diffusion half-step, so this is semi-implicit and still needs time-step convergence checks. Local electron–lattice exchange is solved cell-by-cell by a bracketed implicit root; each cell transfers equal and opposite energy. The heat source and diffusion are applied as balanced half-steps around that exchange.

The specified source is unchanged: Ag-masked \(\exp(-r^2/R^2)\exp(-d/l_{skin})\), with \(R=30\) Å, \(l_{skin}=20\) Å and peak flux \(I_0=5\) eV/(ps·Å²). The 10 Å QSGS grid integrates to 0.863925 mW; no post-mask renormalization is applied. The Ag electron parameter table spans 300–20,000 K and extrapolation remains disabled.

By default x/y are periodic and z is insulated. `--bottom-sink` adds an ideal 300 K Dirichlet bath at the bottom z face to both the Ag electron conduction path and any conducting lattice phase; energy removed by the bath is recorded as a negative diffusion contribution in the balance. This is a boundary-condition sensitivity case, not a measured substrate interface.

## Verification

The full model test suite passes **53 tests**. The nine double-temperature tests cover 0D energy exchange, 1D electron and lattice diffusion, conductivity partition, active-cell enthalpy inversion, QSGS source/coupling accounting, sampled interval accounting, diffusion steps above the explicit stability bounds, and bottom-sink energy removal. The added coarse-face audit test checks that the diagnostic distinguishes a volume-averaged link from direct source-geometry Ag contact. Electron and lattice diffusion conserve energy for insulated boundaries to solver tolerance. In synthetic sink tests, both channels remove energy and the external loss closes the total-energy ledger.

For a 0.1 ps Case A heating run, four outer time steps were compared. The energy input is 8.63925×10⁻¹⁷ J in all cases and the total-energy residual remains below 1.1×10⁻³⁰ J. Cumulative electron-to-lattice transfer varies by about 0.05%, and peak lattice temperature by 0.0278 K. Electron-peak differences decrease with refinement: 0.5 vs 0.3125 fs differ by 2.29%, 0.3125 vs 0.15625 fs by 1.66%, and 0.15625 vs 0.078125 fs by 0.78%. This is promising for the short transient, but it does not establish convergence near melting or over 10 ps.

| Outer step | Peak \(T_e\) | Peak \(T_l\) | Electron-to-lattice transfer | Max Ag liquid fraction |
|---:|---:|---:|---:|---:|
| 0.5 fs | 1088.19 K | 300.628 K | 3.62412×10⁻¹⁸ J | 0 |
| 0.3125 fs | 1063.80 K | 300.617 K | 3.62486×10⁻¹⁸ J | 0 |
| 0.15625 fs | 1046.17 K | 300.606 K | 3.62561×10⁻¹⁸ J | 0 |
| 0.078125 fs | 1038.06 K | 300.600 K | 3.62594×10⁻¹⁸ J | 0 |

## 10 Å QSGS pilot: 5 fs heat, then 5 fs source-off cooling

| Boundary | Peak \(T_e\) at heat-off | Peak \(T_l\) at heat-off | Peak \(T_e\) after 10 fs | Peak \(T_l\) after 10 fs | Energy input | Final Ag liquid fraction | Total-energy residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| A: scaffold adiabatic | 805.759 K | 300.01985 K | 339.950 K | 300.02588 K | 4.31962×10⁻¹⁸ J | 0 | −1.99×10⁻³³ J |
| B: ideal thermal contact | 805.759 K | 300.01659 K | 339.950 K | 300.02095 K | 4.31962×10⁻¹⁸ J | 0 | −1.53×10⁻³³ J |

The total source-off electron-to-lattice transfer through 10 fs is 3.5773×10⁻²⁰ J in A and 3.5773×10⁻²⁰ J in B. Net insulated-boundary diffusion is at round-off scale. Neither case approaches Ag melting; no stable-pool controller or 100 ps run has been attempted.

The 0.5 fs Case A 10 fs pilot takes about 4 s in this Python implementation on the current environment. Linear scaling would put a 100 ps, 10 Å run at roughly 11 hours per boundary case, before adding the controller or refining to 5 Å. This rough estimate is not a production performance guarantee.

The 1 ps Case A heating check at 0.3125 fs injects 8.63925×10⁻¹⁶ J; cumulative electron-to-lattice transfer is 1.58674×10⁻¹⁶ J. Peak electron and lattice temperatures are 2143.13 K and 313.177 K, respectively; the Ag liquid fraction remains zero. The total-energy residual is 4.28×10⁻²⁹ J.

The 5 ps Case A heating check at 0.3125 fs injects 4.31962×10⁻¹⁵ J; cumulative electron-to-lattice transfer is 1.70776×10⁻¹⁵ J. Peak electron and lattice temperatures are 3838.71 K and 421.243 K. All sampled Ag liquid fractions remain zero, and the final total-energy residual is 4.15×10⁻²⁸ J. The electron table ceiling is not approached. The full 5 ps run took 31 min in the current environment; a 100 ps run at this step and grid would be on the order of 10 hours per case if runtime scaled linearly. No melt-pool controller is active.

## 10 Å QSGS pilot: 10 ps heating-only screening

The 10 ps Case A run used the same 0.3125 fs outer step, no bottom heat sink, and the existing Ag-masked lateral Gaussian/depth-exponential source. It completed in 61 min 12 s. The grid-integrated source power is 0.863925 mW and total injected energy is 8.63925×10⁻¹⁵ J. Cumulative electron-to-lattice transfer is 4.45563×10⁻¹⁵ J. The final energy-balance residual is 6.24×10⁻²⁸ J (about 7.2×10⁻¹⁴ of the injected energy); the net insulated electron diffusion ledger is 3.47×10⁻²⁵ J.

| Time | Ag-weighted mean (T_e) | Peak (T_e) | Ag-weighted mean (T_l) | Peak (T_l) | Liquid Ag fraction |
|---:|---:|---:|---:|---:|---:|
| 1 ps | 1953.31 K | 2143.13 K | 310.731 K | 313.177 K | 0 |
| 5 ps | 3719.21 K | 3838.71 K | 413.918 K | 421.243 K | 0 |
| 10 ps | 4660.35 K | 4765.63 K | 590.415 K | 602.880 K | 0 |

All three 3-D snapshots have zero Ag liquid fraction. The electron and lattice fields are already domain-wide by 1 ps under the saved half-rise definition: 100% of Ag volume lies above half the peak rise from 300 K, spanning the full 148.7 Å depth and yielding a periodic lateral (2r_{RMS}) of 246.8 Å. At 10 ps, electron peak-minus-mean is 105.28 K and lattice peak-minus-mean is 12.46 K. Thus this setting produces broad, nearly uniform Ag-network heating rather than a shallow local molten pool; it does not meet the melt-pool goal. Do not call this a stable pool or increase source power based on this result.

The rapid spreading is consistent with the candidate electron properties: the local homogeneous diffusivity (k_e/C_e) is about 0.0069–0.0226 m²/s across 4765–300 K, corresponding to an uncoupled √(αt) scale of roughly 83–150 nm over 1 ps, larger than the model's roughly 30 nm lateral width. This estimate ignores phase geometry and electron–lattice coupling and is explanatory rather than a substitute for the computed fields.

The exact face-connectivity audit finds that 8.4–10.8% of positive-conductance faces at 10 Å (4.4–6.7% at 5 Å) lack direct Ag contact in the source QSGS mask, although both the coarse solver graph and source-contact graph remain a single Ag component. These local false links can affect the distribution, but the candidate Ag electron diffusivity and small domain are already sufficient to explain rapid domain-wide spreading; this run does not isolate their quantitative contribution. The 3-D data are in `qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs_fields.npz`; the cut-plane comparison is `qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs_fields_slices.png`.

The 10 Å QSGS projection contains about 92.8% mixed Ag/scaffold cells. This preserves volume fractions but can blur narrow channels and create apparent cross-face connectivity. Any first liquid region found on this grid requires a 5 Å geometry/connectivity check before it is called a connected pool.

An exact local-face audit now compares the 5 Å and 10 Å finite-volume neighbor links with direct Ag contact in the binary QSGS source geometry. The coarse solver treats two adjacent cells as conducting whenever both have positive Ag fraction. Of these positive links, 8.4–10.8% at 10 Å and 4.4–6.7% at 5 Å have no direct Ag area across that shared face in the source geometry. The solver-link graph and direct-contact graph nevertheless each have one Ag-connected component at both resolutions. Thus volume averaging adds local links but does not add a new domain-spanning Ag component here; the local heat-flow error remains unquantified. Reproduce with `python phase_field_qsgs/audit_coarse_ag_face_connectivity.py --spacing 10A` and `--spacing 5A`; machine-readable details are in the matching `qsgs_ag_face_connectivity_*.json` files.

## Reproduction

Run from the model directory:

```bash
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary B --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
```

The histories and machine-readable summaries are `qsgs_ttm_enthalpy_10A_BCA.csv/.json` and `qsgs_ttm_enthalpy_10A_BCB.csv/.json`. To evaluate the ideal bottom sink, add `--bottom-sink`; sink and insulated results must be reported separately.

The longer screening commands use `--cool-ps 0` to save a heating-only history. The 0.1 ps step-refinement, 1 ps and 5 ps histories are stored as `qsgs_ttm_enthalpy_10A_BCA_heat0p1ps_dt*.csv/.json`, `qsgs_ttm_enthalpy_10A_BCA_heat1ps_dt0p3125fs.csv/.json` and `qsgs_ttm_enthalpy_10A_BCA_heat5ps_dt0p3125fs.csv/.json`.

```bash
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.5 --sample-ps 0.01 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat0p1ps_dt0p5fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.01 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat0p1ps_dt0p3125fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.15625 --sample-ps 0.01 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat0p1ps_dt0p15625fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.1 --cool-ps 0 --dt-fs 0.078125 --sample-ps 0.01 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat0p1ps_dt0p078125fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 1 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.1 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat1ps_dt0p3125fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 5 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.25 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat5ps_dt0p3125fs.csv
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 10 --cool-ps 0 --dt-fs 0.3125 --sample-ps 0.1 --snapshot-ps 1,5,10 --output phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs.csv
python phase_field_qsgs/plot_ttm_field_snapshots.py phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs_fields.npz
```

## Remaining gates

The solver implementation and boundary accounting now pass short numerical checks. Before interpreting melt-pool size or lifetime, the next gates are (1) refine time-step accuracy over longer heating intervals, (2) test the 5 Å grid against the 10 Å result, (3) assess candidate property-table provenance and high-temperature behavior, and (4) build and validate the spatially connected liquid-pool heat-off controller. These pilots do not establish physical accuracy of the candidate Ag parameters or real Ag–Ti₃SiC₂ heat transfer.
