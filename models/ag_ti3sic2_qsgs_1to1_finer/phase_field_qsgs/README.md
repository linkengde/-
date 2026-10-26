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

`audit_coarse_ag_face_connectivity.py` checks local Ag heat-conduction links introduced by volume averaging. At 10 Å, 91.28% of cells are mixed and 8.4–10.8% of coarse faces with positive Ag fractions on both sides have no direct Ag contact across that face in the saved binary QSGS geometry. At 5 Å, the corresponding values are 63.11% mixed cells and 4.4–6.7% such faces. The coarse solver graph and the direct source-contact graph each remain one Ag-connected component at both resolutions, so these extra local links do not create a new domain-spanning Ag component in this geometry. The audit integrates the piecewise-constant source geometry over each coarse face exactly; it does not quantify global heat-flow error. Machine-readable results are in `qsgs_ag_face_connectivity_10A.json` and `qsgs_ag_face_connectivity_5A.json`.

The 10 Å, Case A, no-bottom-sink 10 ps heating-only screening is archived in `qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs.csv/.json`; its 1, 5, and 10 ps three-dimensional fields are in the matching `_fields.npz`, with a cut-plane overview in `_fields_slices.png`. It stayed fully solid, while the half-rise temperature region covered all Ag volume by 1 ps. This is a domain-wide heating response, not a stable local melt pool; see `QSGS_TTM_ENTHALPY_PILOT.md` for metrics and limitations.

The follow-up 0.1 ps localization screen compares legacy volume-fraction links with exact source-QSGS Ag face-contact fractions, the original box with a 3×3 lateral periodic repeat, and a 5 Å corrected-grid check. It finds a strong small-box effect, material resolution sensitivity, and no Ag liquid; it does not establish a converged local melt pool. Results and limits are documented in [`QSGS_HEAT_LOCALIZATION_SCREEN.md`](QSGS_HEAT_LOCALIZATION_SCREEN.md), with filtered metrics in `qsgs_localization_screen_metrics.json` and the field comparison in `qsgs_heat_localization_screen_0p1ps.png`. The face-contact fields are generated by `qsgs_face_contacts.py`. The solver clips projected phase fractions below 1e-12 to remove machine-epsilon phase slivers from the active thermal domain.

The subsequent source/grid audit replaces thermal-cell-center source sampling with an analytic integral over the original binary QSGS Ag voxels. A matched 5 fs, 3×3-domain comparison gives grid-independent source power and good agreement for several bulk/local metrics, but the electron peak and q99.9 tail remain resolution-sensitive; the spatial-convergence gate is therefore not passed. A 10 Å time-step refinement is within 1.3% for the electron peak, while a 1×1 versus 3×3 repeat confirms that global temperature statistics depend strongly on lateral box size. There is still no liquid Ag. Full metrics, thresholds, run commands, and limits are recorded in [`QSGS_SOURCE_GRID_CONVERGENCE.md`](QSGS_SOURCE_GRID_CONVERGENCE.md), with the CSV/JSON histories and 3D field snapshots alongside it.

The matched 0.02 ps follow-up saves 1, 5, and 20 fs fields on both 5 Å and 10 Å grids and conservatively remaps electron energy, lattice enthalpy, and cumulative local electron-to-lattice transfer. Source and energy ledgers pass, and the common-grid hotspot diameter agrees closely, but electron high-tail and local energy/enthalpy maps still exceed the 5% spatial screen. No liquid Ag appears. The grid comparison is therefore still not converged and the runs stop at 0.02 ps; see [`QSGS_MATCHED_GRID_LOCALIZATION_0p02ps.md`](QSGS_MATCHED_GRID_LOCALIZATION_0p02ps.md) and `qsgs_localization_matched_0p02ps_grid_analysis.json`.

The subsequent strict-nesting audit creates a new 60×62×32 fine grid aligned 2:1 with the 30×31×16 coarse grid, preserving the earlier 60×61×32 projection unchanged. Phase fraction, exact source integrals, and Ag face-contact areas now coarsen to roundoff. A 20 fs source-only ablation maps to the coarse electron-energy field essentially exactly, but enabling electron diffusion yields a 10.72% common-Ag excess-energy difference; full coupling gives a 9.88% lattice-enthalpy difference and still no liquid Ag. Geometry mapping was not the main remaining error, so the 5% spatial gate remains failed and heating is not extended. See [`QSGS_NESTED_GRID_TRANSPORT_ABLATION_0p02ps.md`](QSGS_NESTED_GRID_TRANSPORT_ABLATION_0p02ps.md), the nested grid generator, and `nested_grid_diagnostics/nested_grid_20fs_mode_comparison.json`.

The follow-up [`QSGS_ELECTRON_DIFFUSION_AND_CHANNEL_AUDIT.md`](nested_grid_diagnostics/QSGS_ELECTRON_DIFFUSION_AND_CHANNEL_AUDIT.md) validates the uniform-grid electron operator against Fourier modes and tests QSGS effective conductance on three local voxel windows. The uniform operator passes; the projected 10 Å and 5 Å transport graphs overpredict the voxel-network conductance by 19–28% and 33–47%, respectively. The mismatch is not monotonic with grid refinement, showing that exact face contact fractions do not preserve the resistance of disconnected/tortuous Ag subcell paths. A bounded CG tolerance update reduces iterations by about 19% with sub-nanokelvin field differences. There is still no melt-pool result, and the local transport mismatch must be addressed before longer heating. Reproduce the local conductance and CG checks with `python phase_field_qsgs/qsgs_transport_benchmark.py` and `python phase_field_qsgs/benchmark_electron_diffusion_solver.py`.

## Reproduce

From the model directory:

```bash
python phase_field_qsgs/project_qsgs_mask.py
python phase_field_qsgs/analyze_projected_connectivity.py
python phase_field_qsgs/audit_coarse_ag_face_connectivity.py --spacing 10A
python phase_field_qsgs/audit_coarse_ag_face_connectivity.py --spacing 5A
python phase_field_qsgs/plot_ttm_field_snapshots.py phase_field_qsgs/qsgs_ttm_enthalpy_10A_BCA_heat10ps_dt0p3125fs_fields.npz
python phase_field_qsgs/stefan_benchmark.py
python -m unittest -v test_qsgs_projection.py test_thermophysical_model.py test_qsgs_enthalpy_solver.py test_qsgs_ttm_enthalpy.py
python phase_field_qsgs/qsgs_enthalpy_solver.py --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0.1
python phase_field_qsgs/qsgs_ttm_enthalpy.py --spacing 10A --boundary A --heat-ps 0.005 --cool-ps 0.005 --dt-fs 0.5 --sample-ps 0.005
```

The compressed 5 Å and 10 Å maps and their numeric report are stored beside this file as `qsgs_phase_fractions_5A.npz`, `qsgs_phase_fractions_10A.npz`, and `qsgs_projection_report.json`.
