# Ag–Ti₃SiC₂ composition audit and 2× mesh study

Date: 2026-10-27. The user confirmed the experimental target as **70:30 Ag:Ti₃SiC₂ by volume**. This directory holds an isolated mask derived from the existing QSGS growth history and its thermal-validation runs; the parent QSGS geometry and prior melt/freeze results remain unchanged.

## Composition basis

The parent thermal mask `qsgs_morphology.npz` has Ag/Ti₃SiC₂ volume fractions of `41.6654/58.3346%`. At the currently used candidate solid/reference densities (Ag `10,490 kg/m³`, Ti₃SiC₂ `4,520 kg/m³`), that mask implies a bulk-density mass ratio of `1.65763:1`.

The confirmed `70/30 vol%` target implies `5.41519:1` mass and `84.412 wt%` Ag under those same candidate densities. Therefore the earlier statement that `70/30 vol%` converts to `1:1 mass` is not supported by this density conversion. The separate atomistic candidate's approximately `1:1` mass balance belongs to its own geometry convention and is not the selected continuum composition. Exact density and porosity conventions still need a source audit before quantitative material claims.

`build_70_30_volume_geometry.py` thresholds the existing `qsgs_growth_rank.npy` realization and applies the same alternating largest six-face-connected component cleanup as `build_model.py`. It does not regenerate stochastic growth. The original mask SHA-256 remains `7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43`; it was not edited.

The independent output `qsgs_morphology_ag70vol.npz` reaches Ag/Ti₃SiC₂ volume fractions of `69.9996/30.0004%` on the `100×100×50` native grid. Both phases are a single face-connected component with periodic x/y seam connections and span the full open z extent. This checks periodic connectivity, not minimum neck thickness, statistical representativeness, or an atomistic geometry at this composition. See `composition_70vol_geometry.json` for parent/rank hashes and the selection threshold.

## Existing 2× runtime preflight (historical composition only)

The earlier 1 fs (`0.001 ps`) benchmark used the old 41.6654 vol% Ag mask. It is retained only as a code/runtime baseline and does **not** validate the newly confirmed 70/30 target.

| Metric | 1× | 2× |
|---|---:|---:|
| Grid | 100×100×50 | 200×200×100 |
| Ag volume fraction | 0.416654 | 0.416654 |
| Exact source power (W) | 8.8621792423×10⁻⁴ | 8.8621792423×10⁻⁴ |
| Stable timestep (fs) | 0.0404211 | 0.0101053 |
| Runtime for 1 fs (s) | 1.063 | 27.731 |
| Peak temperature after 1 fs (K) | 313.5727 | 313.5713 |
| Energy residual (J) | −5.78×10⁻³⁴ | 7.70×10⁻³⁴ |
| Peak RSS (MiB) | 139 | 754 |

Linear extrapolation from that short run put 0.2 ps at about 5,546 s for the old 2× geometry. This is not a cost guarantee: tracking/output and changed composition can alter throughput.

## Next validation

The matched mesh study must use the isolated 70/30 vol% mask for both resolutions, uniform 300 K initial conditions, the existing Case A boundary and Ag properties, the unchanged Gaussian-plus-depth-exponential source, and analytic source integration. Compare temperatures after block-coarsening the 2× fields and compare Ag liquid fractions, connected-cluster metrics, and integrated energy at shared times. Preserve exact 0.2 ps restart states and distinguish a short heating-window comparison from a passed mesh-convergence claim.

Current thermal property values remain first-pass candidates; literature/source-table verification remains open. Boundary A is the zero-storage, zero-conductivity scaffold limit, with insulated external faces. Neither result should be described as validated experimental heat flow or environmental cooling.

## 70:30 volume Case A heating and mesh comparison

Both resolutions were independently initialized at uniform 300 K and heated from 0 to 0.2 ps using the same Case A boundary, 5 eV/(ps·Å²) peak flux, 30 Å Gaussian radius, 20 Å absorption depth, native QSGS face network, exact source integration, and candidate thermal properties. The 2× grid is a strict voxel subdivision of the same newly selected 70:30 volume mask.

| Result at 0.2 ps | 1× | 2× |
|---|---:|---:|
| Grid | 100×100×50 | 200×200×100 |
| Exact source power (W) | 1.5513589462382562×10⁻³ | 1.5513589462382560×10⁻³ |
| Peak temperature (K) | 748.6093 | 723.1792 |
| Ag liquid volume fraction | 0 | 0 |
| Connected liquid clusters | 0 | 0 |
| Relative energy-balance error | −1.02×10⁻¹⁴ | 4.94×10⁻¹³ |
| Wall time (s) | 229.5 | 5,835.2 |
| Final checkpoint read-back | passed | passed |

Both calculations remain fully solid through 0.2 ps; this selected 70:30 vol% model does not reproduce a melt pool in the tested window. The maximum temperatures are only `0.606` and `0.586` of the candidate Ag melting temperature (`1235.08 K`). The old 41.6654 vol% Case A pool results are a separate composition and cannot be transferred to this target.

After 2× block averaging, the Ag-weighted temperature field has a mean absolute difference of `0.404 K`, RMS difference `1.415 K`, and maximum local difference `28.63 K`. Its L1 difference is `3.274%` of the 1× Ag temperature rise. The peak is 25.43 K lower at 2× (5.67% of the 1× rise), and its lateral position is closer to the heat-source center (4.39 Å at 2× versus 16.55 Å at 1×). Integrated source powers agree to relative roundoff and the injected-energy difference is about `1.0×10⁻¹²` relative.

This is useful evidence that the coarse maximum and hotspot position are resolution-sensitive, even though the domain-average temperature and the no-melt conclusion agree. A single 1×/2× pair is not a convergence proof; no melt-pool size or lifetime exists to converge in this 0.2 ps window. The prior 3×3 validation used the old 41.6654 vol% mask, so it does not establish finite-size independence for this new composition. The new report `mesh2_ag70_comparison.json` contains field, hotspot, liquid, and energy metrics; each resolution also has its timeseries and restart checkpoint under `ag70vol_melt_window/`.

The 1:1 mass volume fraction calculated from the candidate densities remains Ag `30.1133 vol%` / Ti₃SiC₂ `69.8867 vol%`. This is a distinct conditional composition. The active target remains the user's explicit `70:30 vol%` instruction, which corresponds to about `5.415:1` mass under the same candidate densities.
