# V19 CIF Ag–Si convergence bundle

This folder gathers the same frozen-geometry numerical pilot records produced across windows A and B. It contains copies; canonical calculation archives and the report remain at their original repository paths. No calculation was rerun to assemble this bundle.

## Contents

- `calculations/`: compact Gamma, 2×2×1, 3×3×1, and 4×4×1 archives, including each archive’s original `SHA256SUMS.txt`, GPAW log, progress, summary, verification, and output extxyz. Per-job launcher logs are alongside their archive folders. A owns Gamma and 3×3; B owns 2×2 and 4×4. Large `state.gpw` checkpoints are intentionally excluded.
- `inputs/` and `input_manifest.json`: exact fixed 52-atom geometry and registered method/mesh records.
- `mesh_comparison_*.json`: four pairwise analyses.
- `assessment/`: combined report, CSV, machine-readable assessment and its own SHA256 list.
- `compare_mesh_results.py` and `verify_result.py`: read-only archive verification and comparison tools. This bundle intentionally contains no DFT launch scripts.
- `bundle_manifest.json` and root `SHA256SUMS.txt`: copy provenance and checksums for this complete bundle.

Input geometry SHA256: `8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`. All four records are numerical-only pilots, not training or validation labels.

## Verify and reproduce comparisons

From this directory, run `python3 verify_result.py LABEL calculations/LABEL` for each archive. Pair comparisons can be reproduced with `python3 compare_mesh_results.py LABEL0 LABEL1`; the registered set is Gamma/2×2, 2×2/3×3, 2×2/4×4, and 3×3/4×4.

The current B handover is in `coordination/RESTART_V19_NEXT_STEPS.md`: new independent environments must claim their roles there. Do not rerun these four completed meshes. Further 5×5/6×6 calculations require the new owner-bound registrations described in that handover.
