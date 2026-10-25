# Pure-Ag TTM–MD reference pilot

This is a short numerical coupling check for scheme 3 (pure-Ag atomistic reference plus a separate QSGS phase-change model). It advances a real Ag EAM slab with LAMMPS for 100 steps at 0.5 fs while the Python TTM reference grid receives the selected source and exchanges energy with the Ag kinetic energy. A matched no-TTM NVE trajectory measures the integration-drift baseline.

It is **not** a melt-pool calculation, does not include Ti₃SiC₂ or QSGS confinement, and does not qualify the Python reference operator as a production LAMMPS TTM implementation.

## Model and source

- 4,000 FCC Ag atoms, `Ag_u3.eam`, lattice parameter 4.09 Å.
- 1 ps NVT at 300 K to settle the two free surfaces, followed by 0.05 ps NVE/TTM–MD. The measured starting Ag temperature is about 306.57 K.
- `p p f` atomic boundary conditions. The LAMMPS box has 20.45 Å of vacuum padding on both sides in z; the TTM grid covers the 40.9 Å Ag slab only, with zero-flux z boundaries.
- Source: (Q=(I_0/l_{skin})\exp(-r^2/R^2)\exp(-d/l_{skin})), with (I_0=5\ \mathrm{eV/(ps\,Å^2)}), (R=30\) Å and (l_{skin}=20\) Å. The Gaussian convention is the repository's existing `exp(-r²/R²)` definition. Source power is conservatively mapped to cells with CIC and is **not** renormalized after masking.
- Candidate Ag (C_e,G,k_e) table; low-temperature closure is explicitly enabled in this reference test. This tests the code path and energy handoff, not the physical accuracy of the candidate parameters.

## Results

| Requested grid spacing | Grid | Source power (eV/ps) | Input in 0.05 ps (eV) | Final electron mean / max (K) | Final Ag temperature (K) | Energy residual after subtracting source (eV) | Matched NVE drift (eV) | Residual minus NVE drift (eV) | Max electron substeps per MD step |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 Å | 5×5×5 | 5203.45361 | 260.172681 | 4351.111 / 4359.198 | 306.812 | 8.343×10⁻⁵ | 5.288×10⁻⁵ | 3.055×10⁻⁵ | 214 |
| 5 Å | 9×9×9 | 5203.45361 | 260.172681 | 4351.111 / 4359.563 | 306.812 | 8.319×10⁻⁵ | 5.288×10⁻⁵ | 3.031×10⁻⁵ | 697 |

For both grids, the TTM-only energy ledger closes to below (2.4\times10^{-12}) eV. The combined atom-plus-electron residual differs from the matched no-TTM NVE drift by about (3.1\times10^{-5}) eV, within the pilot's (10^{-3}) eV numerical check. The two grids give nearly identical mean electron temperature and electron-to-lattice energy transfer (about 1.373 eV). The 5 Å grid requires 33,450 electron substeps over this 100-step pilot versus 10,355 for the 10 Å grid, so finer grids carry a substantial runtime cost in the explicit Python reference solver.

The source deposits about 260 eV into the electronic system over 50 fs, while the Ag lattice remains near 307 K. This is a short nonequilibrium coupling response; it is not evidence of melting. The 50 fs duration is far below the timescale needed to assess a stable melt pool or resolidification.

## Reproduce

The prepared cloud environment uses the LAMMPS/MPICH Python environment, which has NumPy and an `mpi4py` build linked against the same MPICH runtime:

```bash
cd /workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/pure_ag_classifier
LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib OMP_NUM_THREADS=4 \
  /workspace/.venvs/lammps/bin/python pure_ag_ttm_md_pilot.py
PURE_AG_TTM_GRID_A=5 LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib OMP_NUM_THREADS=4 \
  /workspace/.venvs/lammps/bin/python pure_ag_ttm_md_pilot.py
```

The script writes one metrics JSON and matched LAMMPS logs per requested grid size. Set `PURE_AG_TTM_STEPS` to change the short MD duration and `PURE_AG_TTM_SOURCE_SCALE` only for a deliberate sensitivity run.

The checked-in 10 Å and 5 Å LAMMPS transcripts are stored as gzip-compressed copies: [10 Å TTM](log.pure_ag_ttm_md_pilot_10A.gz), [10 Å NVE reference](log.pure_ag_nve_reference_10A.gz), [5 Å TTM](log.pure_ag_ttm_md_pilot_5A.gz), and [5 Å NVE reference](log.pure_ag_nve_reference_5A.gz). Decompressing them reproduces the raw logs named in the metrics JSON.

The unsuffixed pilot log/metrics files in this directory are preserved from an earlier prototype run that lacked the 1 ps surface pre-equilibration and reported the gross energy change before subtracting source input. Use only the `_10A` and `_5A` metrics/logs above for the reported comparison.

## Interpretation and next gate

This check verifies real LAMMPS atom timesteps, the reference electron-grid source update, the electron-to-lattice energy handoff, p/p/open-grid bookkeeping, and two mesh resolutions for a dense pure-Ag slab. It does not test the composite geometry, interfacial heat resistance, phase-field model, automatic heat-off criterion, melt identification, or solidification. Continue scheme 3 by validating the continuum enthalpy/phase-field model on a simple analytic case and then map the QSGS geometry into its Ag-domain mask; keep the pure-Ag result as an atomic-scale mechanism reference only.
