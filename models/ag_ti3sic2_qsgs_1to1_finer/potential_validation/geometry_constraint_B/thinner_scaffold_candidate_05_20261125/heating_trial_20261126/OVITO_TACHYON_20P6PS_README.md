# OVITO Tachyon figure: 20.6 ps

The three PNG files in `analysis_q6_ptm_msd_to20p6ps/` were rendered from the saved 20.600 ps dump using OVITO 3.17.2 and its Tachyon renderer. The full-model view uses all 276,658 atoms; the close-up keeps atoms within 72 Å of the XY center and z=112–315 Å. PySide6 is used only to arrange the two OVITO-rendered views and add titles/legend.

## Coloring and interpretation

- Black: fixed Ti/Si/C scaffold.
- Pale blue: Ag outside the two highlighted groups.
- Red: Ag with `qbar6 < 0.45` and PTM structure type `Other`. This is a local disorder screen, not standalone evidence of liquid phase.
- Orange: Ag above the initial material top, z=150 Å.

The title reports the mean TTM lattice temperature of Ag at this saved frame (about 4,172 K). The figure itself does not establish evaporation or a crater; see `HEATING_20PS_IDEAL_BOTTOM_SINK_REPORT.md` for the structural and trajectory analysis.

## Re-rendering

Run from an environment with OVITO Python 3.17.x, PySide6, a Vulkan loader, and a Vulkan ICD:

```bash
python analysis_q6_ptm_msd_to20p6ps/render_20p6ps_OVITO_Tachyon.py
```

On Debian, the software Vulkan ICD is provided by `mesa-vulkan-drivers`. In the original cloud session it was absent, so the Debian package `mesa-vulkan-drivers 25.0.7-2+deb13u1` was downloaded and extracted under `/tmp`; OVITO was launched with `VK_ICD_FILENAMES` set to its `lvp_icd.json` and `LD_LIBRARY_PATH` including the extracted library directory. No system-wide package installation was made.
