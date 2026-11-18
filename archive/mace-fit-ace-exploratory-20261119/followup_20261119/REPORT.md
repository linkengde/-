# Follow-up compute window: 1×6×3 DFT, ACE/MTP refits, and validation

This add-on updates the earlier compute plan without changing its frozen baseline. It uses only the MACE-derived Ag–Ti₃SiC₂ candidate data and the verified GPAW label below. It does not mix Round7 labels into this fit.

## DFT label and k-point check

The 32-atom `Ti_gap_2p3_rank_recovery_v17_01` structure converged with GPAW 26.7.0, PBE/PW 500 eV, Fermi width 0.1 eV, and a 1×6×3 mesh in 51 SCF iterations. The total energy is −207.19188074196245 eV.

Against the previous 1×4×2 result, the denser mesh changes the energy by 0.272102 eV total (8.503 meV/atom); force-component RMSE is 0.05028 eV/Å and maximum per-atom vector difference is 0.20985 eV/Å. Against the Γ-only result, the energy difference is 150.301 meV/atom and force-component RMSE is 1.604 eV/Å. The Γ-only label is not mesh-equivalent. The 1×6×3 label replaces the lower-mesh label for this frame in the follow-up candidate copy; the older label remains in the original archive for audit.

## ACE refits

Both refits meet the user’s exploratory aggregate development targets (energy RMSE ≤10 meV/atom and force-component RMSE ≤0.1 eV/Å) on the seven-frame development set:

| Ti-gap force scale | Dev energy RMSE | Dev force-component RMSE | Worst frame force RMSE | Max all-Si-site vector error (proxy) |
|---:|---:|---:|---:|---:|
| 1.0 | 4.530 meV/atom | 0.05054 eV/Å | 0.08019 eV/Å | 0.2087 eV/Å |
| 0.05 | 4.565 meV/atom | 0.05060 eV/Å | 0.07985 eV/Å | 0.2153 eV/Å |

The formal Si hotspot-center criterion cannot be scored because the center IDs have not been registered. The all-Si-site values above are proxies. The development set is small and source-correlated, so these aggregate scores do not establish transferability.

## Source-family grouped validation

Each family was held out as a complete group. These folds are diagnostic cross-validation, not a blind test. The best full-force ACE candidate gives:

| Held-out family | Frames | Energy RMSE (meV/atom) | Force-component RMSE (eV/Å) | Si-site max error proxy (eV/Å) |
|---|---:|---:|---:|---:|
| JN1967 source family | 35 | 457.454 | 1.0682 | 2.1999 |
| MACE Ag–C interface family | 10 | 1.408 | 0.1083 | 0.2211 |
| MACE Ag–Si interface family | 10 | 1.001 | 0.0964 | 1.2552 |
| MACE Ag–Ti interface family | 7 | 2.883 | 0.0874 | 0.0857 |
| v17 framework and Ti-gap family | 9 | 9.244 | 0.1404 | 0.0950 |

The JN1967 holdout fails badly, and the v17 holdout exceeds the force target. ACE is an exploratory candidate, not a source-robust or production-qualified model.

## MTP force-weight sweep

All variants used the converged 1×6×3 training label. None meets the force-component target on development data:

| Force weight | Train E RMSE (meV/atom) | Train F RMSE (eV/Å) | Dev E RMSE (meV/atom) | Dev F RMSE (eV/Å) | Dev Si-site max proxy (eV/Å) |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 4.110 | 0.1481 | 11.050 | 0.2150 | 1.411 |
| 0.10 | 6.863 | 0.1795 | 7.401 | 0.1317 | 0.668 |
| 0.20 | 68.167 | 0.5776 | 73.271 | 0.4609 | 1.210 |

The 0.10 model meets the energy target but fails the force and Si-site proxy targets. Do not use these MTP variants for production dynamics.

## LAMMPS loading checks

ACE weight 1.0 and 0.05, and all three MTP variants, loaded in the matching LAMMPS builds and completed `run 0` on representative 52–65 atom cells. PyACE and LAMMPS PACE agree for the two ACE candidates within 3.2×10⁻⁶ eV total energy and 1.1×10⁻⁶ eV/Å force-component RMSE. These checks verify model loading and implementation consistency, not physical accuracy.

## Label-design shortlist and supplied contact table

`acquisition/label_design_center_shortlist.csv` contains 81 spatially separated centers selected from the archived million-atom model-force dump. Selection is by element and z-bin, keeps centers at least 8 Å apart in the x/y-periodic, z-open cell, and uses only atoms whose model force exceeded 0.5 eV/Å. These are DFT acquisition candidates, not DFT-labeled structures; local environments still need to be extracted and reviewed before any calculation.

The separately supplied `full_model_severe_contact_edges.csv` has 12,062 edges with listed distances 1.000–2.008 Å. Comparing its atom IDs against the archived million-atom QSGS target finds all IDs in range but zero matching endpoint element types. It is therefore a different atom-ID geometry and is preserved under `unmatched_source_input/` without being added to the MACE fitting set. Its source structure is needed before those rows can be used.

I also checked the available MACE QSGS fine/finer and relaxed-backbone geometries, the Round7 V10 input and wrapped start, and three local Round7 checkpoints. None reproduces the table: the QSGS variants have zero matching endpoint types; the 347,464-atom Round7 snapshots have only 14 matching edge types out of 12,062, and none of those 14 distances agrees within 0.001 Å. The tested files and counts are in `unmatched_source_input/source_geometry_search.json`. The table's source geometry remains unidentified.

As a diagnostic of the JN1967 fold, the full-force ACE fitted on all training data fits the 21 training-frame `config_type` labels that explicitly name JN1967 at 0.155 meV/atom energy RMSE and 0.00492 eV/Å force-component RMSE. The 35-frame grouped holdout still fails across unseen configurations: the Ag–Ti subset has 643.97 meV/atom and 1.586 eV/Å RMSE; Ag–C has 618.83 meV/atom and 0.754 eV/Å; Ag–Si has 25.97 meV/atom and 0.806 eV/Å. This is a transfer failure across withheld source configurations, despite low training error on the included JN1967-named frames.

## Decision and next work

Keep the full-force ACE as the leading exploratory candidate for follow-up label acquisition. Do not run production or heating dynamics with it yet. Next, identify the source geometry for the supplied contact table, generate and DFT-label reviewed local structures from the 81-center shortlist, and add coverage for the JN1967 family. Refit and rerun source-family validation only after those labels are reviewed. No blind test set exists yet.

The zero-step million-atom force audit remains a static diagnostic of the unrelaxed target geometry; its high model forces are not DFT errors. Full target integration was not run.
