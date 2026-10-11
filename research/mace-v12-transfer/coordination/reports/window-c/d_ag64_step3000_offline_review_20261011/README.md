# Ag64 step3000 offline review package

This package analyzes the SHA-verified D step3000 archive using saved coordinates and logs only. It compares identical 500-step windows at step2000 and step3000, verifies exact trajectory overlap, and records thermal, RDF, Q6, coordination, CNA and MSD outputs. No MD, DFT or training is run by the analysis script.

Run with the existing analysis environment:

```bash
/workspace/.venvs/mace/bin/python analyze_step3000.py
```

Primary conclusion is in `REPORT.md`; machine-readable results are in `metrics.json`. `temperature_comparison.png` and `structure_comparison.png` plot the principal comparisons. CSV files preserve the window statistics, per-frame structural values, RDF and MSD. CNA is conventional fixed-cutoff CNA; PTM was not available in the existing workflow. `SOURCE_EVIDENCE.json` records archive, checkpoint and provenance identity.

`LIQUID_CANDIDATE_STATUS=PRELIQUID` is diagnostic only: the segment remains at 1600 K and is not a 1250 K liquid candidate. This review does not authorize further MD.
