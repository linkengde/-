# QSGS geometry input for the continuum phase-change route

`project_qsgs_mask.py` reads the saved binary QSGS morphology, maps each source voxel to a 5 Å or 10 Å thermal cell by exact axis-aligned overlap averaging, and writes Ag and Ti₃SiC₂ volume-fraction fields. It does not threshold mixed cells or infer a material interaction.

The source geometry is the repository's `qsgs_morphology.npz`: 100×100×50 voxels over 299.635×301.8575×158.6115 Å. The original voxel pitches are 2.996×3.019×3.172 Å. Its ceramic skeleton fraction is 0.583346 and complementary Ag geometry fraction is 0.416654.

## Projection results

| Requested cell size | Grid | Actual cell size (Å) | Ag volume fraction | Mixed cells | Ag volume error |
|---:|---:|---:|---:|---:|---:|
| 5 Å | 60×61×32 | 4.994×4.948×4.957 | 0.416654 | 73,912 / 117,120 | 0 Å³ |
| 10 Å | 30×31×16 | 9.988×9.737×9.913 | 0.416654 | 13,583 / 14,880 | <1×10⁻⁹ Å³ |

The 10 Å map is highly mixed at cell scale, so 5 Å is the more useful starting resolution for preserving the observed geometry. A separate 0.5-threshold, 26-neighbor diagnostic with periodic x/y and open z finds one connected skeleton component at both resolutions. The largest Ag component contains 99.9917% of occupied cells at 5 Å and 99.9831% at 10 Å; it reaches both z faces, with four and one smaller Ag fragments, respectively. This does not establish periodic winding rank. Keep the unthresholded fractions in any subsequent finite-volume calculation; do not silently replace them with thresholded masks. Details are in `qsgs_projected_connectivity_report.json`.

## Still required before a physical phase-change run

The repository does not currently contain a selected and sourced Ag lattice enthalpy model (solid/liquid heat capacity and latent heat in a consistent unit convention), nor a validated thermal transport law for the composite phase-field grid. It also does not determine heat transfer through the Ti₃SiC₂ scaffold. The next solver prototype should take these as explicit parameters and compare declared limiting boundary cases, such as zero heat flux through the scaffold versus a prescribed thermal sink. Those assumptions must remain visible in the output; they are not substitutes for an atomistic Ag–Ti₃SiC₂ interface.

No heating, melting, or solidification result is produced by this directory yet.

## Reproduce

From the model directory:

```bash
python phase_field_qsgs/project_qsgs_mask.py
python phase_field_qsgs/analyze_projected_connectivity.py
python -m unittest -v test_qsgs_projection.py
```

The compressed 5 Å and 10 Å maps and their numeric report are stored beside this file as `qsgs_phase_fractions_5A.npz`, `qsgs_phase_fractions_10A.npz`, and `qsgs_projection_report.json`.
