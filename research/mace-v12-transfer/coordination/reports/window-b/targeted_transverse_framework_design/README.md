# Window B: targeted transverse and framework proposals

This is a geometry-only V16 training-acquisition design following the matched V14/V15 force-error diagnosis. It keeps Ag registry motion separate from local framework-contact motion. No DFT, MACE, MD or TTM work was run, and no labels or dataset/model files were changed.

## Deliverables

- `candidate_set_final/proposal_report.md`: design, marked-pair projection/radial changes, contact screening and A review notes.
- `candidate_set_final/proposal_manifest.json`: source and geometry hashes, persistent atom IDs, axes, displacement vectors, pair distances/projections, element-pair minima, split/input deduplication and source inventory.
- `candidate_set_final/screening_summary.json`: concise screening records.
- `candidate_set_final/*.extxyz`: 12 geometry-only proposals, two displacement modes with opposite signs for each of AgC, AgSi and AgTi.
- `candidate_set_final/generate_training_proposals.py`: reproduction script.
- `candidate_set_final/SHA256SUMS.txt`: hashes for the script, candidate geometries, manifest, summary and report.
- `SHA256SUMS.txt`: bundle-wide inventory covering this README and every candidate-set file, including its internal hash list.

The Ag arm translates all Ag atoms by ±0.15 Å transverse to the marked Ag–X vector; the framework arm moves only the selected marked contact atom by ±0.04 Å along its localized V15 force-residual direction. The selected framework IDs are C#604804, Si#366934 and Ti#479027. The generated structures preserve atom order and IDs, cell and PBC, and carry no energy or force labels.

All element-pair minimum-distance checks passed against the empirical same-composition envelope from existing splits and official geometry inputs, with a 0.05 Å allowance. This is a comparative screen, not a universal bond cutoff. The exact input-geometry screen found no exact hits, and no candidate was near a frozen V16 validation geometry. The conservative near screen did flag historical split/input neighbors: AgTi transverse minus has 17 file/frame hits, closest at 0.0565 Å all-atom RMS. A should review those hits before assigning that branch for DFT. Other near hits are listed per candidate and are not treated as proof of duplication or independence.

The geometries are proposals only. Their role and DFT timing remain for Window A to decide. Keep them out of the frozen validation set.
