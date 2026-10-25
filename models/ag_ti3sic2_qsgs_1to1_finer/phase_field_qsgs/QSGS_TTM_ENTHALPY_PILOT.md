# QSGS continuum two-temperature–enthalpy pilot

## Status and scope

The finite-volume QSGS enthalpy model is now coupled to an Ag-only electron-temperature field. This is a **continuum two-temperature model**, not atomistic TTM-MD. The 0D, 1D, and short 3D numerical checks pass their energy ledgers. No Ag has melted in the completed coupled pilot, and no stable-pool controller or 100 ps run has been attempted.

The old single-temperature QSGS enthalpy pilot remains a separate method check. It directly heated the lattice and its boundary-A transient melt is not evidence from this two-temperature calculation.

## Coupled equations and implementation

The electron grid reuses `AgParameters` and `AgTTMGrid` from `ag_ttm_reference.py`. It evolves tabulated electron internal energy (U_e(T_e)), electron conduction (k_e(T_e)), source input, and (G_{e-ph}(T_e)(T_e-T_l)) exchange. The table `Ag_TTM_hybrid_figure_fit_300_20000K.csv` covers 300–20,000 K; low-temperature and high-temperature extrapolation are disabled. These are cross-paper candidate fits, not a validated Ag–Ti₃SiC₂ property set.

The QSGS Ag fraction weights electron capacity, conduction, coupling, and source exactly once. The source is converted to the grid-integrated eV/ps field from the already Ag-masked (Q_e), preserving the specified (\exp(-r^2/R^2)\exp(-d/l_{skin})\) shape, (R=30) Å, (l_{skin}=20) Å, and (I_0=5\) eV/(ps·Å²); no post-mask power renormalization is done. The 10 Å grid integral is 0.863925 mW.

The exchanged lattice energy is applied with the opposite sign to Ag electron energy. The lattice enthalpy also receives the conservative finite-volume lattice-conduction operator. Case A has Ag-only capacity/conduction and an adiabatic geometric scaffold; case B includes the Ti₃SiC₂ capacity/conduction with shared cell temperature. They remain limiting boundary assumptions rather than measured interface heat transfer. x/y are periodic and z faces are insulated; an external substrate heat sink, radiation, and convection remain absent.

Both the lattice and electron explicit solvers subcycle internally. On the 10 Å grid, the cold-state lattice conduction bounds are 0.2107 fs for case A and 0.4213 fs for case B. In addition, the reused electron solver applies its Ce/ke/G stability and source-temperature limits. This preserves stability but makes long calculations expensive.

## Numerical checks

The full repository suite passes **47 tests**. Five new tests cover zero-dimensional electron-to-lattice exchange, 1D electron diffusion, 1D lattice diffusion, sampled interval accounting, and the 3D QSGS source/coupling ledger for both boundaries. Earlier Ag enthalpy tests also cover the latent-heat plateau and a 1D Stefan benchmark. These checks validate implementation-level conservation and bounded stepping, not the physical accuracy of the candidate property tables.

The matched 10 Å QSGS pilot used 5 fs outer intervals with internal subcycling: 5 fs heat followed by 5 fs with the source off.

| Boundary | Grid source power | Energy injected | Electron-to-lattice energy | Peak (T_e) at heat-off | Peak (T_l) at heat-off | Peak (T_l) at 10 fs | Final Ag liquid fraction | Total-energy residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A: scaffold adiabatic | 0.863925 mW | 4.31962×10⁻¹⁸ J | 3.58855×10⁻²⁰ J | 727.270 K | 300.01785 K | 300.02165 K | 0 | −9.04×10⁻³¹ J |
| B: ideal contact | 0.863925 mW | 4.31962×10⁻¹⁸ J | 3.58860×10⁻²⁰ J | 727.270 K | 300.01194 K | 300.01124 K | 0 | −9.06×10⁻³¹ J |

The integrated electron and lattice diffusion terms are near zero, consistent with the periodic x/y and insulated z boundaries. At this very short interval, both cases remain solid; the difference between A and B is small. This does not contradict the much larger differences in the earlier single-temperature runs, which use a different heat-input model and must not be compared as if they were the same physics.

The A run takes 3,478 internal electron substeps during each 5 fs sampled interval. Linear scaling would put a 100 ps calculation near 70 million electron substeps at this grid resolution, before moving to 5 Å. The current Python explicit reference is therefore a verification tool, not yet a practical production solver.

Reproduce the short matched pilot from the model directory:

```bash
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.005 --cool-ps 0.005 --dt-fs 5 --sample-ps 0.005
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary B --heat-ps 0.005 --cool-ps 0.005 --dt-fs 5 --sample-ps 0.005
```

The reproducible time histories and summaries are `qsgs_ttm_enthalpy_10A_BCA.csv/.json` and `qsgs_ttm_enthalpy_10A_BCB.csv/.json`.

## Next gates

Before interpreting melt-pool lifetimes or attempting the requested 100 ps: (1) complete electron-table provenance review; (2) test timestep convergence; (3) add a separate bottom heat-sink boundary case because the current closed z boundaries only redistribute energy; and (4) replace or accelerate the explicit grid kernels for the long run. Only after those checks should a liquid-cluster persistence/volume-plateau heat-off controller be connected. No conclusion about real interface heat transfer or a stable pool follows from this pilot.
