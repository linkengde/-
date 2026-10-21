# Ag–Ti₃SiC₂ model geometry image

`MODEL_GEOMETRY_COMPARISON.png` compares the requested bicontinuous geometry with the structure that was actually used by the now-stopped heating run.

- **Left:** 39,486-atom Scheme A geometry-only bicontinuous reference. It contains 12,342 Ag atoms and 27,144 Ti₃SiC₂ atoms in a 110 × 98.085 × 100.921 Å cell with X open and Y/Z periodic. The archived geometry audit reports each phase as a connected network spanning Y/Z. The 50 Å X-vacuum is omitted from the plot. This is a geometry reference, not a force-field or dynamics pass.
- **Right:** last saved frame of the 1,084-atom layered model that was used by continuation 7 before it was stopped. It contains 484 Ag and 600 Ti₃SiC₂ atoms. It is a slab/interface model and is not bicontinuous.
- **Colors:** Ag is gold; Ti₃SiC₂ is teal.

The layered continuation was stopped by SIGTERM at 2026-10-21 09:29:14 UTC after the user requested stopping it. Its last saved checkpoint is the 9.715 ps cumulative-heating frame (absolute trajectory time 9.735 ps). The available saved metrics do not establish Ag melting. The 39,486-atom bicontinuous source model was not modified or heated by this action.

The PNG was rendered with the reproducible `render_model_comparison.py` script from the archived OVITO-readable extxyz structures using ASE and Matplotlib.
