# Round7 1.600 ps Q6 + PTM + displacement screen: code audit

## Scope and source records

This audit uses the actual analysis script, its JSON output and per-atom CSV, the 1.6 ps run summary, and the Round5 calibration files copied under `round5_threshold_provenance/`. The relevant analysis is `extract_round7_joint_snapshot.py`; its current invocation is reflected by `round7_1p6ps_joint_screen.json`: endpoint 1.60 ps, window 0.60–1.60 ps, five frames separated by 0.25 ps. The five timestamps are 0.60, 0.85, 1.10, 1.35, and 1.60 ps.

There is no separate saved stdout/configuration log for the classifier invocation. Its effective settings are preserved in the script and JSON. The per-atom CSV retains the five-frame mean Q6, five-frame PTM non-FCC fraction, 1 ps displacement statistic, and stored joint flag for the 6,448 atoms selected by the script's source-region mask. It does not retain each frame's PTM class or RMSD, so those per-frame details cannot be reconstructed from this CSV alone.

## Parameters and provenance

| Item | Actual setting / behavior | Code or record | Provenance and validation |
|---|---|---|---|
| Frames | 0.60, 0.85, 1.10, 1.35, 1.60 ps; five snapshots, 0.25 ps apart | Script lines 28–31; JSON `window_ps`, `sample_interval_ps` | Confirmed in JSON and full-trajectory package; this is a coarse five-snapshot screen. |
| Q6 algorithm | OVITO Steinhardt per-particle `q_l`, `max_l=6`; uses `Steinhardt` component index 4 (`q6`); no neighbor-averaged `q̄6` | Lines 76–93 | Confirmed against OVITO 3.17.2 component names (`q2`…`q6`); `compute_averaged` remains its default `False`. |
| Q6 neighbors | Cutoff mode, 4.0 Å, uniform neighbor weights; `only_selected=True`, selection is particle type 4 (Ag), so only Ag centers and Ag neighbors participate | Lines 78–84 | Explicit in script; OVITO 3.17.2 API says unselected particles are treated as nonexistent. No separate neighbor-count/coordination cutoff is applied. |
| Q6 rule | Arithmetic mean of each atom's five framewise q6 values `< 0.4081069282833649` (strict inequality) | Lines 34–35, 186, 189 | Number inherited from Round5; not recalibrated on Round7. |
| PTM algorithm | OVITO Polyhedral Template Matching; `only_selected=True`, following the same Ag selection; PTM classification is calculated for Ag-only neighborhoods | Lines 76–93 | Explicit in script; modifier runtime was OVITO 3.17.2. |
| PTM RMSD | `rmsd_cutoff` is not set in this script; OVITO 3.17.2 default is 0.1 (normalized PTM RMSD, not Å). Matches above the cutoff are classified Other. `output_rmsd` is also not enabled. | Lines 85–87; OVITO 3.17.2 installed API defaults | Default was checked from `/workspace/.venvs/qsgs` runtime. The JSON does not serialize this default or per-atom RMSDs. |
| PTM categories | FCC code 1 is the sole FCC category. HCP (2), BCC (3), and Other (0) all count as non-FCC; the script does not distinguish them. Other does not bypass the other filters: an atom still needs Q6, displacement, and spatial-mask passes. | Lines 36–37, 187, 189 | The installed OVITO default has FCC/HCP/BCC enabled; other template classes are disabled by default. Any non-1 result is treated alike. |
| PTM persistence | Non-FCC in at least 80% of five frames, i.e. at least 4 of 5 frames | Lines 37, 187, 189 | Threshold inherited from Round5's solid/mobile-reference separation; not independently calibrated for Round7. |
| Displacement statistic | For each of four adjacent frame intervals, calculate coordinate increment, apply x/y minimum-image correction with 200 Å periodic lengths, subtract the mean increment of all atoms assigned initially to the Newton layer, then sum the four corrected increments. Per-atom output is the squared norm of this net vector, `|Σ Δr_corrected|²`. | Lines 96–103; called at line 188 | A single 0.60→1.60 ps start/end displacement; not time-origin-averaged MSD and not a sliding/time-averaged MSD. z is not wrapped because the boundary is open. Wrapped x/y increments are minimum-image-unwrapped interval by interval. |
| Displacement threshold | `>= 0.5171679303519052 Å²` for this 1 ps window | Lines 38–39, 189 | Numerically identical to Round5's 1 ps cutoff; script scales it linearly by `WINDOW_PS`, here 1.0. Not recalibrated on Round7. |
| Spatial prefilter | Initial `40 ≤ z < 150 Å`, type 4 Ag, and `(dx/30)² + (150-z)/20 ≤ 2`; x is minimum-image distance from 100 Å. Although `dy` is computed, it is never used in this expression. | Lines 45–61 | **Implementation defect:** this is an x–z band spanning all y, not a central radial/circular source-core region. The 6,448 count and all candidate totals below refer to this code mask, not a verified 3-D central heat-source core. |
| Joint rule | `core & (q6_mean < Q6_SCREEN) & (nonfcc >= 0.8) & (msd_i >= MSD_SCREEN_A2)` | Line 189 | Strict AND: all three indicators must pass for the same atom, in addition to the spatial prefilter. No temperature, coordination, or separate local-neighbor-count condition is used. |
| Cluster definition | Endpoint positions at 1.60 ps of joint candidates; Ag–Ag graph edges at distance ≤ 3.5 Å; connected components; no minimum cluster size | Lines 106–159, 190–196 | Explicit in script. Every singleton is counted as a component. |
| Cluster boundaries | x/y periodic with 200 Å box lengths; z is open. A constructed z-period prevents artificial z-wrap links. Dimensions unwrap x/y using the largest circular gap before bounding-box calculation. | Lines 111–139 | Explicit in script; JSON records `periodic_xy: true`. Convex-hull volume is a geometric envelope, not a thermodynamic melt volume. |

## Round5 threshold source

The numerical thresholds were copied from `round5_threshold_provenance/calibrated_thresholds.json` and `CALIBRATION_REPORT.md`. Round5 compared 864-atom pure-Ag references generated with a different, force-only development MTP: a 300 K FCC/NVT solid and a dense, mobile 2600 K candidate. The Q6 threshold 0.4081069283 is the midpoint between the solid q6 5th percentile (0.4471880204) and dense-mobile q6 95th percentile (0.3690258362). The PTM persistence threshold 0.8 lies between the measured 0% and 100% non-FCC reference fractions. The 1 ps displacement threshold 0.5171679304 Å² is the midpoint between the solid 99th percentile (0.2688601405 Å²) and dense-mobile 5th percentile (0.7654757202 Å²).

For Round7, the Q6 cutoff, PTM fraction, and 1 ps displacement cutoff were all inherited unchanged. The displacement cutoff was multiplied by `WINDOW_PS`; because this run uses a 1.0 ps window, its numerical value is unchanged. No Round7-specific solid/liquid calibration was performed. The reference potential and Round7 energy-force pilot MTP differ, so transferring the thresholds is provisional.

## Recomputed pass counts

Counts below were recomputed from the 6,448 CSV rows using the exact strict/inclusive comparisons in line 189. The stored `joint_candidate` flag matches the recomputed three-way AND for all 6,448 rows (zero mismatches).

| Screen | Number passing |
|---|---:|
| Q6 mean < 0.4081069283 | 3,078 |
| PTM non-FCC fraction ≥ 0.8 | 6,418 |
| 1 ps squared net displacement ≥ 0.5171679304 Å² | 6,304 |
| Q6 ∩ PTM | 3,067 |
| Q6 ∩ displacement | 3,000 |
| PTM ∩ displacement | 6,276 |
| All three (strict AND) | 2,989 |
| Fail at least one | 3,459 |

Failure patterns among the 3,459 that fail are: Q6 only 3,287; displacement only 78; PTM only 11; Q6 and displacement 64; Q6 and PTM 17; all three 2; PTM and displacement but Q6 passes 0. The main reason for rejection is therefore Q6, not PTM or displacement. “Fail” means only that at least one threshold was missed; it does not by itself prove that the atom is solid.

## Cluster result and interpretation

The 2,989 joint-screen atoms form 148 endpoint connected components at 3.5 Å. The largest has 1,683 atoms (56.3% of candidates); there is no minimum-size filter. The component contains 105 singleton components across the full candidate population. The largest component's periodicity-aware bounding dimensions are about 84.32 × 146.96 × 40.21 Å; its convex-hull envelope is about 201,756 Å³. These are geometrical statistics of selected candidate atoms and do not establish a spherical, continuous liquid pool. The skeleton may interrupt connectivity under this particular cutoff.

The 3,459 non-passing atoms are not all “solid”: they are candidates that missed one or more arbitrary, inherited thresholds. In particular, the Q6 rule is a low-order screen, and low coordination or missing neighbors can lower Q6; this workflow has no coordination safeguard. PTM Other can likewise reflect disorder, surfaces, or a poor local template match. Only five widely spaced frames are used for structural persistence.

## Suitability for melt-pool research

This exact workflow is suitable for exploratory visualization and for locating Ag atoms that are simultaneously low-order, persistently non-FCC, and mobile under the chosen development potential. It is not yet a validated liquid classifier or a quantitative melt-pool measurement. Round7's potential is not independently validated for liquid Ag or the Ag/Ti₃SiC₂ interface; thresholds were transferred from Round5; the original spatial prefilter omitted y; the displacement measure is a single-origin squared net displacement despite its `msd` name; and cluster connectivity depends on a provisional 3.5 Å cutoff. The original result should be called a **Q6/PTM/displacement joint candidate region**, not confirmed molten Ag. A separate reanalysis has now corrected the radial spatial mask while retaining the same five frames and thresholds: it selects 1,787 Ag, finds 874 joint candidates in 41 components, and has a largest component of 693 atoms. See `../round7_corrected_radial_core_analysis_1p6ps_20261109/`. These corrected-mask counts resolve the geometric-mask error, but Round7-specific calibration and potential validation are still outstanding.
