# Ag–Ti₃SiC₂ model geometry image

`MODEL_GEOMETRY_COMPARISON.png` compares the requested bicontinuous geometry with the structure that was actually used by the now-stopped heating run.

`BICONNECTIVE_YZ_SLICES.png` shows nine Y/Z phase slices from the 39,486-atom model. Blue is Ag and gray is Ti₃SiC₂. These cuts reveal that this model has a regular, periodically repeating wavy/gyroid-like morphology. It is geometrically bicontinuous, but its morphology does not match the irregular, labyrinth-like reference image supplied by the user. The slice pixels are assigned from the nearest atom within a ±1.2 Å X slab; they are a visualization of the atom labels, not a reconstructed continuous phase-field.

- **Left:** 39,486-atom Scheme A geometry-only bicontinuous reference. It contains 12,342 Ag atoms and 27,144 Ti₃SiC₂ atoms in a 110 × 98.085 × 100.921 Å cell with X open and Y/Z periodic. The archived geometry audit reports each phase as a connected network spanning Y/Z. The 50 Å X-vacuum is omitted from the plot. This is a geometry reference, not a force-field or dynamics pass.
- **Right:** last saved frame of the 1,084-atom layered model that was used by continuation 7 before it was stopped. It contains 484 Ag and 600 Ti₃SiC₂ atoms. It is a slab/interface model and is not bicontinuous.
- **Colors:** Ag is gold; Ti₃SiC₂ is teal.

The layered continuation was stopped by SIGTERM at 2026-10-21 09:29:14 UTC after the user requested stopping it. Its last saved checkpoint is the 9.715 ps cumulative-heating frame (absolute trajectory time 9.735 ps). The available saved metrics do not establish Ag melting. The 39,486-atom bicontinuous source model was not modified or heated by this action.

The PNGs were rendered with the reproducible `render_model_comparison.py` script from the archived OVITO-readable extxyz structures using ASE and Matplotlib.
