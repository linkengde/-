# step06 Ti3SiC2 skeleton figure

This folder contains the paper-style overview and local branch inset rendered from the accepted step06 model.

## Main image

- `Ti3SiC2_step06_whole_skeleton_with_local_branch.png` — 2920 × 1900 px PNG at 300 dpi.
- `Ti3SiC2_step06_whole_skeleton_with_local_branch.tif` — 300 dpi LZW TIFF.
- `Ti3SiC2_step06_whole_skeleton_with_local_branch.pdf` — PDF export.

The whole Ti3SiC2 skeleton is shown at left. The blue dashed box selects the 25 × 25 × 25 Å region shown at right. The inset is one continuous forked body: the voxel phase field has one 6-neighbor-connected component, the surface mesh has one connected component, and the extracted centerline has five terminals and two junction clusters. The orange double arrow marks the actual 3D distance between the two facing branch surfaces: **13.25 Å (1.325 nm)**.

The figure uses the step06 model's 200 × 200 × 200 Å material region. The displayed surfaces are extracted from `step06_equal_volume_Ti3SiC2_mask`; the atomic coordinates and the Ag/Ti3SiC2 model are unchanged by figure generation. The arrow tips correspond to visible points on the isosurface mesh, and their 3D separation is computed in Å.

## Rebuild

Run `build_step06_skeleton_meshes.py` in the model visualization environment, then run `render_step06_skeleton_figure.py` in OVITO 3.17.2 with a Vulkan renderer. The render uses the Mesa lavapipe ICD when hardware Vulkan is unavailable. Both scripts and their local helper modules are included here.

The exact endpoints, camera, ROI, and image dimensions are recorded in `step06_figure_report.json`.
