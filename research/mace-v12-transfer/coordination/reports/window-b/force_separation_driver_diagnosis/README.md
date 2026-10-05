# V14/V15 force-driver diagnosis

The authoritative matched-model comparison is in [`analysis_final/diagnosis_report.md`](analysis_final/diagnosis_report.md). It uses the same three already-scored V15 AgC/AgSi/AgTi registry probes for both selected models, reproduces A's V15 scores, checks source hashes/IDs/PBC, and decomposes marked-pair projections and per-atom force errors.

## Result for A

| Interface | V14 force RMSE | V15 force RMSE | V14 separation error | V15 separation error |
|---|---:|---:|---:|---:|
| AgC | 0.14288 | 0.13528 | 0.31501 | 0.38744 |
| AgSi | 0.33589 | 0.31328 | 0.01228 | 0.03630 |
| AgTi | 0.11687 | 0.11751 | 0.30798 | 0.35013 |

Units are eV/Å. The same-geometry force vector RMSE improves slightly for AgC/AgSi and is effectively unchanged for AgTi, while separating-force absolute error rises for all three. AgSi's small scalar error hides a large transverse component. The result supports keeping V16's distance and lateral-registry data axes separate; it does not identify one causal source because V15 changed five training rows and added three residual rows at once.

**Immediate V16 training guidance:** hold the current recipe fixed to read the effect of reviewed new labels: seed 45, CPU, batch 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, same foundation model and architecture. Do not tune those values from the three scored probes. If an ablation is needed after source/split review, compare force weight 1000 vs 2000 on a fixed dataset and stronger validation set, keeping all other factors fixed; this remains a proposal for A.

## Files

- `analysis_final/common_geometry_scores.csv`: matched per-interface scores and V15−V14 deltas.
- `analysis_final/force_decomposition.json`: per-atom regions/species/components, signed Ag/X terms, transverse vectors, top IDs and neighbors.
- `analysis_final/audit_evidence.json`: model and reference hashes, split edits, training selections/settings, method provenance and six V16 input-only geometry mapping.
- `analysis_final/analyze_matched_v14_v15.py`: reproducible CPU analysis script.
- `analysis_final/SHA256SUMS.txt`: hashes for final analysis outputs; root `SHA256SUMS.txt` covers the full report directory.

The `attempt_diagnostics` directory retains the two initial report-packaging failure records. They were corrected; no reference, model, or DFT result failed verification. No V16 reference output, new DFT, training, MD, or TTM was used.
