# V9–V16 cross-version review artifacts

This directory contains B's cumulative historical review for A before the V17 recipe is frozen.

- `report.md`: evidence-linked timeline, matched metrics, error interpretation, ranked hypotheses, and actionable V17 recommendation.
- `analysis_04/common_geometry_scores.csv` and `.json`: eight frozen model versions scored on six already-scored references, with original reference roles preserved.
- `analysis_04/version_ledger.csv` and `.json`: model/data/training source paths, hashes, counts, known recipe values, and limits.
- `analyze_v9_v16_common.py`: CPU scoring and summary driver.
- `analysis_01/`–`analysis_03/`: retained failed setup/compatibility diagnostics from earlier attempts; they are not result data.
- `SHA256SUMS.txt`: hashes of the report, driver, final score files, and ledger.

The final inference used one CPU thread with float64. V15/V16 outputs and applicable V14 scores were reused only after hash checks. The report does not claim fresh validation or causal isolation.

To verify published artifacts from this directory:

```bash
sha256sum -c SHA256SUMS.txt
```

Do not execute the driver as part of a training or DFT queue. It performs model inference on the archived, already-scored references and is retained for reproducibility only.
