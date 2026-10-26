# Electron diffusion and QSGS channel-conductance audit

## Decision

The constant-property electron diffusion operator passes its analytic checks. The remaining 5 Å/10 Å mismatch is not explained by a sign error, a nonconservative diffusion stencil, or the backward-Euler time integrator in a uniform Ag domain.

The local QSGS resistor-network check does **not** pass a geometry-transport agreement screen: both projected grids overpredict the voxel-resolved effective Ag conductance, and the 5 Å result is farther from that reference than the 10 Å result in all three sampled windows. This is evidence that exact face area alone does not preserve the effective resistance of subcell Ag paths. Do not tune Ag properties to compensate and do not extend the present melt screening on the assumption that the 5 Å transport field is converged.

No melting calculation was run in this audit.

## 1. Constant-property diffusion benchmark

Added `test_electron_diffusion_benchmarks.py`. It constructs pure-Ag periodic 1D grids with singleton y/z dimensions and an exactly linear electron energy law, `Ue = Ce (Te - 300 K)`, with constant `Ce = ke = 1` in the solver's internal units. These values are manufactured test coefficients, not Ag material data.

The test field is a sinusoidal Fourier mode. For cell spacing `dx`, mode index `m`, and `alpha = ke/Ce`, the finite-volume eigenvalue is

```text
lambda_h = 4 alpha sin²(pi m/N) / dx²
```

The one-step solution is checked against the exact backward-Euler amplification `A1 = A0/(1 + dt lambda_h)`. A second test represents a straight Ag channel by permitting x-face contacts and setting y/z contacts to zero; it checks modal decay, mean-temperature conservation, and net electron-energy conservation. The existing tabulated-property 1D test continues to exercise temperature-dependent Ag properties and heat-flow direction.

| Check | Result |
|---|---:|
| One-step modal amplitude error | < 2e-7 K |
| Mean-temperature change in insulated mode | < 2e-9 K |
| Net electron-energy change | < 1e-28 J |
| Spatial error ratios, N=16→32 and 32→64 | 3.77 and 3.94 |
| Time error ratios, dt=0.05→0.025 and 0.025→0.0125 ps | 1.995 and 1.998 |

The spatial ratios are consistent with second-order convergence of the 1D finite-volume Laplacian; the time ratios are consistent with first-order backward Euler. This validates the uniform-grid operator, not the projected QSGS transport approximation.

## 2. Local QSGS channel transport

Added `qsgs_transport_benchmark.py`. It extracts three reproducible 24×24×16-voxel windows from the repository's QSGS morphology and imposes a unit temperature difference at opposing x faces, with other faces insulated and no periodic wrap. It compares:

1. A binary source-voxel resistor network at the stored QSGS resolution (about 3.00×3.02×3.17 Å voxels).
2. Conservative 10 Å and 5 Å projections using exact binary-geometry Ag face-contact fractions.

All cases use the same manufactured intrinsic `k=1`, so the reported `k_eff/k` isolates the geometry discretization. Ag volume fraction and the imposed hot/cold boundary Ag areas agree between each voxel and projected representation. Each linear solve closes hot/cold flux to better than 1.5e-9 relative.

| Window start (source voxels) | Ag fraction | Voxel reference `k_eff/k` | 10 Å error | 5 Å error |
|---|---:|---:|---:|---:|
| (36, 60, 18) | 0.4966 | 0.1468 | +19.4% | +46.6% |
| (56, 12, 0) | 0.4959 | 0.2049 | +28.4% | +46.9% |
| (64, 76, 0) | 0.4964 | 0.1978 | +24.2% | +32.8% |

The projected result is not monotonically closer as the grid is refined. A projected cell is isothermal and can join separate Ag pieces inside that cell; face-area correction controls transport across cell faces but does not represent internal subcell path connectivity and resistance. The coarse-grained graph therefore creates effective shortcuts in these windows.

This is a numerical geometry diagnostic, not a measurement of Ag–Ti₃SiC₂ interface transport or a validated composite conductivity. The binary voxel network is the reference for the stored geometry only.

Machine-readable results are in `nested_grid_diagnostics/qsgs_local_transport_benchmark.json`. Reproduce with:

```bash
python phase_field_qsgs/qsgs_transport_benchmark.py
```

## 3. Sparse-solver profiling and bounded optimization

For one diffusion half-step on the tiled nested 10 Å QSGS grid (90×93×16 cells; 125,082 Ag-bearing cells), profiling showed approximately 0.57 s in three Jacobi-preconditioned CG calls and 0.028 s building the finite-volume operator. The dominant cost is therefore iterative sparse solution, not rebuilding the graph topology; caching just the matrix sparsity pattern would have limited impact.

Changing only the CG relative tolerance from 1e-10 to 1e-8 reduced the observed iterations from 279 to 226 (about 19%) for the same manufactured transient. The maximum Te difference was 5.0e-10 K and the change in net diffusion energy was 9.4e-30 J. A separate one-step temperature-dependent Ag diffusion/coupling comparison differed by at most 1.1e-9 K in Te and 5.7e-14 K in lattice temperature; both energy-balance residuals remained below 1.1e-33 J. The solver now uses the 1e-8 tolerance, recorded as `ELECTRON_DIFFUSION_CG_RTOL`.

Algebraic multigrid and ILU-GMRES were benchmarked but not adopted. AMG reduced iteration counts in this case but its setup and preconditioner applications increased wall time; ILU setup was also more expensive than the current solve. No new runtime dependency was added. Reproduce the tolerance comparison with `python phase_field_qsgs/benchmark_electron_diffusion_solver.py`; the output is saved as `nested_grid_diagnostics/electron_diffusion_cg_tolerance_comparison.json`.

## 4. What this resolves and what remains

- The finite-volume diffusion sign, conservation, straight-channel heat flow, constant-coefficient Fourier response, and time integration pass their targeted checks.
- The source integration and nested phase/contact-area mapping remain consistent with the earlier audits.
- The local conductance comparison identifies a remaining projected-geometry error: exact face contact area is insufficient when a thermal cell contains disconnected or tortuous subcell Ag paths.
- This audit does not validate temperature-dependent transport on QSGS at 5 Å, liquid Ag, melt-pool size, or solidification.

The next numerical step is to compare the same small windows on grids approaching the native QSGS voxel resolution, or add a conservative subcell connectivity/resistance model. Keep the material table and source fixed. Only after that transport comparison is acceptable should the matched-grid heating screen be extended.

## Reproduction and tests

```bash
cd models/ag_ti3sic2_qsgs_1to1_finer
python -m unittest discover -v
python phase_field_qsgs/qsgs_transport_benchmark.py
```

The complete unit suite passes: 72 tests.
