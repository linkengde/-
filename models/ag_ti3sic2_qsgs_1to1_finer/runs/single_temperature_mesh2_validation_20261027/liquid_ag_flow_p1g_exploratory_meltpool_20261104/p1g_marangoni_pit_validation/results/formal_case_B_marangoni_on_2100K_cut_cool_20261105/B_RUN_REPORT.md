# Formal Case B — Marangoni on

- Start: cold start at 300 K; same mapped Ag70 masks and 64³ frozen grid as Case A.
- Source: fixed composite Gaussian/depth-exponential input, 8.931899 mW.
- Heat off: 0.624017255 ps when measured Ti₃SiC₂ peak reached 2117.286 K. This is a one-step overshoot of the 2100 K target and remains below the 2300 K numerical property ceiling. The source/bottom-sink difference was 6.5974 mW, so the balance stop did not trigger.
- Peak liquid Ag: 121.854858 nm³ at heat off. Largest connected cluster: 110.083468 nm³ (4035 cells), bounding size 10.263×10.565×3.320 nm.
- Maximum whole-domain speed: 62.2876 m/s at 0.781874 ps during cooling with 8.7818 nm³ still liquid; it is not a liquid-phase mean and is not physically calibrated.
- Solidification: 0.861280016 ps; Ti₃SiC₂ peak 1062.723 K, Ag peak 1057.378 K, liquid Ag volume zero; process exit 0 and `DONE_SOLIDIFIED` recorded.
- Maximum absolute Ag-volume drift: 4.7105×10⁻⁵%. Maximum `|energy residual|/enthalpy`: 2.8720×10⁻⁸.
- Corrected upper-envelope surface audit: no pit deeper than the half-cell threshold; maximum pit depth, area, and volume are zero. The maximum surface rise is one cell (0.301857 nm); source-region mean top-surface elevation change is +0.173340 nm at heat-off and +0.202756 nm after solidification.

This is an exploratory numerical result, not a validated prediction for real Ag flow. See [`P1G_FORMAL_AB_REPORT.md`](../../P1G_FORMAL_AB_REPORT.md) for the matched control, data validation, and limitations.
