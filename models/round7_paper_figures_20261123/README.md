# Round7 paper figures 2 and 4

This folder contains first-pass publication figures made by post-processing the saved Round7 full-system trajectory. No MD was rerun.

## Outputs

- Figure2_round7_fiveframe_structure_evolution.png — five full-system atomic snapshots at 0.60, 0.85, 1.10, 1.35, and 1.60 ps.
- Figure2_round7_fiveframe_structure_evolution.tif — 300 dpi LZW TIFF.
- Figure2_round7_fiveframe_structure_evolution.pdf — PDF.
- Figure4_round7_ag_disorder_candidate_region.png — Q6, PTM non-FCC persistence, Ag coordination, and joint candidate positions.
- Figure4_round7_ag_disorder_candidate_region.tif — 300 dpi LZW TIFF.
- Figure4_round7_ag_disorder_candidate_region.pdf — PDF.
- round7_figure_generation_report.json — data provenance, frame counts, spatial screening counts, coordination cross-check, and limitations.
- render_round7_figures.py — reproducible renderer and post-processing script.

## Round7 source files

The source data is stored in the repository branch codex/upload-round7-fiveframe-trajectory-1p6ps-20261109, commit 878328378869ba4e5860a4081420fcdc59966c52.

- models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_fiveframe_joint_analysis_to_1p6ps_20261109/round7_5frames_full_system_0p6_to_1p6ps.dump.gz
- models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_corrected_radial_core_analysis_1p6ps_20261109/round7_1p6ps_source_core_ag.csv
- models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_corrected_radial_core_analysis_1p6ps_20261109/ptm_coordination_diagnostic.json

Example:

    python render_round7_figures.py --trajectory PATH_TO_FIVEFRAME_DUMP_GZ \
      --source-core-csv PATH_TO_SOURCE_CORE_CSV \
      --ptm-diagnostic PATH_TO_PTM_DIAGNOSTIC_JSON --out-dir .

The renderer uses OVITO 3.17.2, NumPy, Matplotlib, and Pillow. On this managed Linux environment, OVITO StandardRenderer uses the local Mesa lavapipe Vulkan driver under /workspace/.local-vulkan.

## Interpretation

Figure 2 uses the same camera, atom radius, and phase colors in all five panels. Atom types 1–3 are Ti₃SiC₂ (270,762 atoms) and type 4 is Ag (81,875 atoms). The box is 200 × 200 × 250 Å, periodic in x/y and open in z. The five stored frames are 0.25 ps apart; no independent 1.50 ps frame is present.

Figure 4 uses the corrected radial source-core selection (1,787 Ag atoms) and the stored joint Q6/PTM/displacement flags. Coordinates are from 1.60 ps; Q6 and PTM metrics summarize 0.60–1.60 ps. Ag-only coordination is recomputed at 1.60 ps using a 4 Å cutoff with periodic x/y and open z. The 874 highlighted atoms are exploratory mobile/disordered candidates, not verified liquid atoms or a validated melt pool. The Q6/PTM/displacement thresholds came from Round5 and have not been calibrated on Round7. The displacement field is a single-origin drift-corrected net displacement squared, not a time-averaged MSD.

The Round7 development MTP and Python reference TTM coupling remain exploratory. Neither figure validates the force field, melting classification, welding behavior, or interfacial heat transfer.
