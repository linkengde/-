# Full-domain QSGS connected-component TTM comparison (20 fs)

## Decision

The connectivity-preserving component graph is now integrated into the dual-temperature/enthalpy solver and has passed the full 3×3-domain numerical comparison against the native Ag-voxel graph for this 20 fs screening case. The source power and Ag volume match exactly; the electron-energy, lattice-enthalpy, and cumulative electron-to-lattice transfer fields are within the project’s 5% relative-L1 screen at both time steps. Halving the time step changes those fields by less than 0.4%.

This is a transport-discretization result only. It does **not** establish a stable melt pool: the Ag liquid fraction stayed zero and the maximum lattice temperature was about 300.17 K. No 0.05 ps extension or melt calculation was run.

## Model and boundary conditions

- Geometry: repository `qsgs_morphology.npz`, a 100×100×50 binary QSGS cell tiled 3×3 in x/y, giving a 300×300×50 source-voxel domain (about 3 Å per source voxel).
- Ag mask: 1,874,943 occupied voxels (41.6654%); the complementary Ti₃SiC₂ mask is retained as geometry only for boundary case A.
- Coarse model: each 2×1×1 source block is split by actual six-face Ag connectivity. The model has 1,143,423 thermal nodes and 2,541,366 true Ag-face links. Node heat capacity and lattice enthalpy use the exact Ag voxel volume. Only real Ag shared faces conduct; the graph includes the periodic x/y supercell seams.
- Reference: one node per native Ag voxel, 1,874,943 nodes and 4,349,520 true Ag-face links.
- Boundaries: periodic x/y, insulated z. Boundary A gives the Ti₃SiC₂ mask no heat capacity or cross-phase heat transfer. This is not an atomistic interface model.
- Source: exact integral over native binary Ag voxels, then conservative summation to component nodes. The source is centered in x/y, lateral `exp(-r²/R²)` with R=30 Å, and decays from the top in depth with 20 Å absorption depth. Peak flux is 5 eV/(ps·Å²). No post-mask renormalization is used.
- Material input: candidate Ag TTM table from 300 to 20,000 K; no extrapolation. Electron heat is coupled locally to Ag lattice enthalpy, including latent heat.

## Results

| Quantity | 0.3125 fs | 0.15625 fs | Screen |
|---|---:|---:|---|
| Integrated source power, both graphs | 8.86217924228131×10⁻⁴ W | 8.86217924228131×10⁻⁴ W | exact match (relative difference 0) |
| Component electron excess-energy field vs native | 1.5726% relative L1 | 1.5767% relative L1 | pass, ≤5% |
| Component lattice-enthalpy field vs native | 1.4865% relative L1 | 1.4892% relative L1 | pass, ≤5% |
| Cumulative electron-to-lattice field vs native | 1.4865% relative L1 | 1.4892% relative L1 | pass, ≤5% |
| Component/native maximum electron temperature | 1259.91 / 1276.36 K | 1258.15 / 1271.80 K | diagnostic |
| Component/native maximum lattice temperature | 300.1727 / 300.1751 K | 300.1723 / 300.1745 K | no melt |
| Ag liquid volume fraction | 0 / 0 | 0 / 0 | no liquid |
| Largest energy residual / source energy | 6.04×10⁻¹⁴ | 1.27×10⁻¹³ | pass |
| Time-step refinement, component fields | — | 0.366% electron energy; 0.384% lattice enthalpy | pass, ≤5% |
| Time-step refinement, native fields | — | 0.364% electron energy; 0.381% lattice enthalpy | pass, ≤5% |

Peak-temperature refinement is also small: the electron peak changes by 1.76 K (0.140%) for the component graph and 4.56 K (0.357%) for the native graph. The lattice peak changes by less than 0.001 K. Liquid fraction remains zero in both runs.

A local diagnostic deserves attention alongside the passing volume-integrated metrics: the largest pointwise electron-temperature difference between the component graph and native graph is 112.8 K at 0.3125 fs and 114.3 K at 0.15625 fs, while the Ag-volume-weighted RMS difference is about 0.658 K. The 5% gate here is the declared relative-L1 criterion for energy/enthalpy fields; the coarse graph is not pointwise temperature-equivalent near every sharp hotspot. Preserve native-voxel comparison for any later local hotspot or melt-onset claim.

## Cost and reproducibility

Measured graph construction took 2.5–3.1 s for the 2×1×1 model and under 0.8 s for the native graph. Peak process memory was 2.16 GiB. With AMG preconditioning, the four full-domain runs took about 3.84 wall-clock hours in total in the recorded environment (NumPy 2.3.5, SciPy 1.17.0, PyAMG 5.3.0, one CPU thread).

A separate one-step smoke comparison on a small QSGS geometry checked AMG against the Jacobi CG preconditioner: maximum `Te` difference was 1.6×10⁻¹¹ K, maximum `Tl` difference was 5.7×10⁻¹⁴ K, and total-energy difference was 5.2×10⁻³³ J. This confirms that the faster full-domain run's optional preconditioner did not change the tested equations at meaningful precision.

Install the optional AMG dependency with `python -m pip install -r requirements-component-amg.txt`, then run from the model directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python phase_field_qsgs/qsgs_component_full_domain_match.py \
  --tile-xy 3 --duration-ps 0.02 --time-steps-fs 0.3125,0.15625 \
  --preconditioner amg
```

For a dependency-free fallback use `--preconditioner jacobi`; the default solver wrapper also retains the original projected transport backend. The run writes the scalar/field-error summary to `qsgs_component_full_domain_20fs_match.json`; progress and run metadata are in the matching `_run.log`. Full three-dimensional fields are not duplicated in Git.

Verification run for this change: `python -m unittest discover -v -s .` (80 tests), plus the full-domain two-resolution/two-time-step calculation above.

## Scope limits

The calculation uses candidate Ag properties and boundary A. Ti₃SiC₂ remains a geometric mask, so its heat capacity, its own thermal conduction, interface heat transfer, mechanical constraint, free-surface motion, and atomistic fixed/Newton/sputter/evaporation layers are not represented. A zero liquid fraction at 20 fs means only that this short numerical comparison did not melt Ag; it is not a conclusion that the longer process cannot form a melt pool. The next planned calculation is the bounded 0.05 ps screening under these explicit assumptions; it is a separate run and was not included here.
