# Ag64 step3000 offline structure review

**Task:** `C-D-AGL-1250K-rho097-seed01`  
**Data:** D branch commit `2ad687d7f74e035a9d2274e7d1e9a70577eaa76f`; source step2000 archive and step3000 segment archive.  
**Review:** offline only; no MD, DFT, or training was started by C. No step3000→4000 work was performed.

## Decision

`LIQUID_CANDIDATE_STATUS = PRELIQUID`

The saved structures are disordered and mobile by the chosen geometric diagnostics, consistent with a provisional hot-stage liquid-like or pre-liquid state. The simulation is still in the 1600 K hot stage, the last-window temperature is rising, and only 25 temperature samples / five saved coordinate frames cover the final 500-step window. It is not a thermally stable 1250 K liquid candidate and should not be frozen or treated as a DFT reference.

## Archive and identity checks

The step3000 archive SHA is `cd00e213b89242fa9501aacbadd5cec8c23e938a61aeebbb2eb138997deea7d2`; its sidecar, receipt, archive member manifest, and internal SHA256SUMS all agree. Task ID, step bounds 2000→3000, 64 atom IDs, cell and PBC were checked. The archived step2000→3000 trajectory contains steps1000–3000; the overlapping step1000–2000 coordinates are bitwise identical to the prior step2000 archive. Step3000 state SHA is `a3c6db34bb5828d46288c61a6065cbcf4da3bd23e99c43f97ae943d4f9d5111a`; state arrays SHA is `920c25e95455f1d490d34b0790afe5b40df3691f32e0674c6f3e075b9fda8157`.

The archive metadata records user manual authorization and says C-signature was not required. This review does not revise that provenance record; it only assesses the archived outputs.

## Temperature

Windows use the same definition as the step2000 review: compare `1500 < step <= 2000` with `2500 < step <= 3000`, each containing 25 records at 20-step spacing (10 fs between records). The final window spans 250 fs by step count; recorded samples run from step2520 to3000, a 240 fs timestamp span. Temperatures fluctuate substantially, so a short-window slope is descriptive rather than a proof of long-term heating.

| Window | Mean (K) | Sample SD (K) | Min–max (K) | OLS trend (K/100 fs) |
| --- | ---: | ---: | ---: | ---: |
| 1500 < step ≤ 2000 | 1523.94 | 154.08 | 1303.43–1789.78 | −46.28 |
| 2500 < step ≤ 3000 | 1581.97 | 169.94 | 1313.81–2069.29 | +148.69 |

The mean rose 58.03 K versus the prior matched window, and the fitted slope changed from cooling to rising. The target in both windows is 1600 K. The end of the trajectory therefore does not establish thermal stationarity; the apparent upward drift warrants continued review if a later segment is separately authorized.

## Structure: RDF, Q6, coordination and CNA

RDFs use five saved frames per comparison window, periodic minimum images, `dr=0.05 Å`, and the same finite-cell normalization. The first peak remains at 2.775 Å; its height falls from 2.743 to 2.649. The first valley moves from 3.975 Å at step2000 to 3.825 Å at step3000. Thus coordination values using each window's own valley are 12.86 and 11.81, respectively, but these are not directly identical-cutoff measurements.

For a direct comparison, coordination and Q6 were also evaluated at the fixed step2000 RDF valley, 3.975 Å:

| Window | Mean coordination | Mean local Q6 | Mean global Q6 |
| --- | ---: | ---: | ---: |
| 1500 < step ≤ 2000 | 12.856 | 0.3454 | 0.1214 |
| 2500 < step ≤ 3000 | 12.744 | 0.3378 | 0.1048 |

The fixed-cutoff changes are small and indicate no return to strong crystalline order. Fixed-cutoff conventional CNA classifies all 64 atoms as `OTHER` in each of the five frames in both windows; no FCC, HCP, or icosahedral atoms were identified. `OTHER` is not itself a liquid label. PTM was unavailable in the existing workflow and was not run. The final-frame minimum distance is 2.410 Å; minimum over the five final-window frames is 2.319 Å, with no geometry in these saved frames indicating a catastrophic close contact.

## Dynamics: MSD

MSD uses the saved trajectory at 100-step (50 fs) intervals, MIC-unwrapped increments, and center-of-mass drift removed. Referencing the common step1000 frame gives 0.659 Å² at step2000 and 1.081 Å² at step3000, an increase of 0.422 Å². The step2000→3000 segment-relative endpoint MSD is 0.716 Å². The fitted segment slope is 1.437 Å²/ps.

This is evidence of atomic rearrangement over the sampled interval, not a reliable diffusion coefficient: the trajectory is short, saved frames are 50 fs apart, and no converged diffusive regime or independent time-origin analysis is established.

## Conclusion and limitations

At step3000, the configuration is best described as **PRELIQUID**: hot-stage, structurally disordered and mobile, with rising but strongly fluctuating temperature. It has not reached the 1250 K sampling stage, so it does not meet `LIQUID_CANDIDATE`. The five-frame structural statistics and short MSD cannot establish equilibrium. CNA is geometric and PTM was not run. The outputs remain diagnostic results only; no structure was frozen as a candidate.

The requested step3000→4000 segment remains outside this review. No further MD was started.
