# QSGS geometry input for the continuum phase-change route

`project_qsgs_mask.py` reads the saved binary QSGS morphology, maps each source voxel to a 5 Å or 10 Å thermal cell by exact axis-aligned overlap averaging, and writes Ag and Ti₃SiC₂ volume-fraction fields. It does not threshold mixed cells or infer a material interaction.

The source geometry is the repository's `qsgs_morphology.npz`: 100×100×50 voxels over 299.635×301.8575×158.6115 Å. The original voxel pitches are 2.996×3.019×3.172 Å. Its ceramic skeleton fraction is 0.583346 and complementary Ag geometry fraction is 0.416654.

## Projection results

| Requested cell size | Grid | Actual cell size (Å) | Ag volume fraction | Mixed cells | Ag volume error |
|---:|---:|---:|---:|---:|---:|
| 5 Å | 60×61×32 | 4.994×4.948×4.957 | 0.416654 | 73,912 / 117,120 | 0 Å³ |
| 10 Å | 30×31×16 | 9.988×9.737×9.913 | 0.416654 | 13,583 / 14,880 | <1×10⁻⁹ Å³ |

The 10 Å map is highly mixed at cell scale, so 5 Å is the more useful starting resolution for preserving the observed geometry. A separate 0.5-threshold, 26-neighbor diagnostic with periodic x/y and open z finds one connected skeleton component at both resolutions. The largest Ag component contains 99.9917% of occupied cells at 5 Å and 99.9831% at 10 Å; it reaches both z faces, with four and one smaller Ag fragments, respectively. This does not establish periodic winding rank. Keep the unthresholded fractions in any subsequent finite-volume calculation; do not silently replace them with thresholded masks. Details are in `qsgs_projected_connectivity_report.json`.

## Material-data limits

The first-pass Ag enthalpy model and a 1-D Stefan test now use the candidate anchors in [`thermophysical_candidates.csv`](thermophysical_candidates.csv) and [`thermophysical_candidates.json`](thermophysical_candidates.json). It uses (T_m=1235.08\) K, (L_f=11.297\) kJ/mol, a linearly interpolated solid (C_p(T)) between the supplied 300 K and melting-point anchors, a single-anchor constant liquid (C_p), and a piecewise-linear candidate (k(T)). Density is held at 10,490 kg/m³, so the model does not include melt expansion. The exact source tables and complete citations are not yet archived in the repository; treat these as the user-specified baseline candidates, not as an independently verified property database.

The 1-D finite-volume one-phase Stefan benchmark uses a hot Dirichlet boundary, an initially solid Ag domain at (T_m), and an insulated far boundary. At 0.01 s, the 40 µm and 20 µm grids put the front within 0.12% and 0.08% of the similarity solution; the enthalpy balance residual is below (1.6\times10^{-9}) J/m². This checks the implemented latent-heat plateau and conduction discretization for that benchmark only.

The Ti₃SiC₂ anchors are also first-pass values. With formula mass 195.708 g/mol, 110–155 J/(mol·K) converts to **562–792 J/(kg·K)**; this corrects the 656–924 J/(kg·K) conversion in the supplied summary. The candidate scaffold table stops at 1473.15 K and the code rejects higher temperatures rather than extrapolating.

The first QSGS 3-D finite-volume enthalpy solver is implemented in [`qsgs_enthalpy_solver.py`](qsgs_enthalpy_solver.py); the Ag electron equation and electron–lattice coupling are added in [`qsgs_ttm_enthalpy.py`](qsgs_ttm_enthalpy.py). The single-temperature checks are recorded in [`QSGS_ENTHALPY_PILOT.md`](QSGS_ENTHALPY_PILOT.md), and the coupled checks in [`QSGS_TTM_ENTHALPY_PILOT.md`](QSGS_TTM_ENTHALPY_PILOT.md). The coupled solver uses conservative semi-implicit diffusion and implicit local electron–lattice exchange. It keeps Ag total conductivity in the independent single-temperature baseline, while the two-temperature model assigns Ag conduction to electrons and defaults Ag lattice conductivity to zero. Short 0D/1D/3D tests preserve the energy ledger, but no stable melt pool has been established. These are continuum two-temperature/enthalpy calculations, not TTM-MD. Cases A and B remain limiting assumptions, not measured Ag–Ti₃SiC₂ interface properties.

## Reproduce

From the model directory:

```bash
python phase_field_qsgs/project_qsgs_mask.py
python phase_field_qsgs/analyze_projected_connectivity.py
python phase_field_qsgs/stefan_benchmark.py
python -m unittest -v test_qsgs_projection.py test_thermophysical_model.py test_qsgs_enthalpy_solver.py test_qsgs_ttm_enthalpy.py
python phase_field_qsgs/qsgs_enthalpy_solver.py --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0.1
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
```

The compressed 5 Å and 10 Å maps and their numeric report are stored beside this file as `qsgs_phase_fractions_5A.npz`, `qsgs_phase_fractions_10A.npz`, and `qsgs_projection_report.json`.
