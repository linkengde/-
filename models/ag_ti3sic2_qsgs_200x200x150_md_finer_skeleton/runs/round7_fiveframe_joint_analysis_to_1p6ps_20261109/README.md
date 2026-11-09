# Round7 Ag–Ti₃SiC₂ TTM-MD continuation to 1.60 ps

This package continues the existing 352,637-atom, 200 × 200 × 250 Å Ag–Ti₃SiC₂ model from the verified 1.4005 ps restart to 1.6000 ps. It is an exploratory run using the Round7 development MTP and the Ag-only TTM reference coupling. It is not a single-temperature run and it is not a validated prediction of welding or melt-pool properties.

## Play the full-system trajectory

Open `round7_5frames_full_system_0p6_to_1p6ps.dump.gz` in OVITO and press Play. It contains exactly five whole-system frames at 0.60, 0.85, 1.10, 1.35, and 1.60 ps. OVITO validation found 352,637 atoms in every frame. No 1.50 ps frame is included.

The 1.60 ps snapshot image is `round7_1p6ps_joint_candidate_3d.png`. Orange marks Ag atoms that pass all three temporary screening criteria; dark blue is the Ti₃SiC₂ skeleton, light blue is the remaining Ag, and gray is source-core Ag that did not pass the joint screen.

## Heating continuation

- Restarted from `ttm_1p4005ps.restart` and `electrons_1p4005ps.npz`; timestep 0.1 fs; continued 1,995 steps.
- Kept the existing central Gaussian TTM source: (I_0=5.0\;\mathrm{eV\,ps^{-1}\,Å^{-2}}), (R=30\;Å), 20 Å surface skin, source surface at z=150 Å and center (100, 100) Å.
- Electronic deposition and electron–lattice exchange act on Ag only; Ti₃SiC₂ is heated by atomistic interactions. The heating operator is a deterministic Python reference coupling, not a production-validated LAMMPS TTM fix.
- At 1.60 ps the integrated source input is 6,212.286 eV, electron-to-lattice transfer is 458.950 eV, the electron-grid residual is −3.84×10⁻¹⁰ eV, and the continuation-segment combined-energy residual is −0.228 eV. The accumulated source-corrected residual is −11.959 eV.
- Final mean temperatures: electron grid 5,341.6 K; Newton layer 2,749.2 K; Newton-layer Ag 2,669.8 K; Newton-layer Ti₃SiC₂ 2,772.0 K. Maximum reported force is 103.3 eV/Å. The maximum z coordinate is 235.84 Å; one atom is in the nominal 200–250 Å evaporation layer. Atom count remains 352,637.
- Recovery snapshots were written at 1.44, 1.48, 1.52, 1.56, and 1.60 ps. The trajectory passes through 1.50 ps continuously but does not save a 1.50 ps frame.

These values demonstrate that the code advanced and kept a small numerical energy residual over the continuation segment. They do not establish that the development potential or heat-transfer model is physically accurate.

## Five-frame Q6/PTM/MSD screen

The screen uses the same five frames spanning 0.60–1.60 ps and the initial Newton-layer Ag subset selected by the legacy spatial mask. Per-atom displacement is drift-corrected and accumulated with x/y minimum-image unwrapping. The temporary rules are:

- mean Ag Q₆ < 0.408106928, using a 4.0 Å neighbor cutoff;
- PTM non-FCC fraction ≥ 0.80 over the five frames;
- drift-corrected squared displacement over the 1.00 ps window ≥ 0.517167930 Å²;
- Ag–Ag cluster connectivity cutoff 3.5 Å with periodic x/y and open z.

The numerical criteria were borrowed from Round5 pure-Ag development-potential references. Round7-specific solid/liquid calibration has not been done. Treat these as mobile, disordered Ag candidates, not confirmed liquid atoms.

The script's spatial prefilter contains 6,448 Ag atoms: 2,989 pass all three criteria and 3,459 fail at least one. A code audit found that the prefilter omits its computed y displacement, so this is an x–z band spanning y, not a true central radial source core. The screen finds 148 Ag candidate clusters. The largest contains 1,683 atoms, with a periodicity-aware bounding box of approximately 84.3 × 147.0 × 40.2 Å and a geometric convex-hull envelope of 201,756 Å³. The elongated, scaffold-interrupted shape and box dimensions do not mean a spherical or fully connected macroscopic melt pool exists. Treat these as provisional joint-screen candidates; see `ROUND7_1P6PS_Q6_PTM_MSD_AUDIT.md` for the exact rules, count breakdown, and limitations.

For comparison, the preceding 0.40–1.40 ps screen found 2,667 joint candidates, 159 clusters, and a largest cluster of 1,454 atoms (177,958 Å³ hull envelope; 83.5 × 171.3 × 33.6 Å bounding box). The rolling-window screen therefore shows more mobile/disordered candidates and a larger largest component at 1.60 ps, while the component count falls. Because the windows are shifted, candidate membership and dimensions change, and the MTP is uncalibrated, this is evidence of increased candidate connectivity only—not proof that a specific pair of pools merged.

## Files

- `round7_5frames_full_system_0p6_to_1p6ps.dump.gz` — playable full-system trajectory.
- `round7_1p6ps_joint_candidate_3d.png` — 3D snapshot of the 1.60 ps screen.
- `round7_1p6ps_joint_screen.json` — complete metrics, rules, and limitations.
- `round7_1p6ps_joint_candidates.dump` — all Ag atoms passing the joint screen in the endpoint frame.
- `round7_1p6ps_largest_candidate_cluster.dump` — endpoint atoms in the largest connected candidate component.
- `round7_1p6ps_source_core_ag.csv` — per-atom Q₆/PTM/displacement values for the Ag subset selected by the legacy spatial mask (filename retained for traceability).
- `ROUND7_1P6PS_Q6_PTM_MSD_AUDIT.md` — code-level parameter audit, Round5 threshold provenance, recomputed intersections, and known spatial-mask defect.
- `round5_threshold_provenance/` — the calibration report, JSON thresholds, reference metrics, and Round5 analyzer used to trace inherited cutoffs.
- `run_summary_1p6ps.json` — continuation summary and checkpoint records.
- `SHA256SUMS.json` — file sizes and SHA-256 checksums for the package files.
- `extract_round7_joint_snapshot.py`, `pack_fiveframe_full_trajectory.py` — analysis and trajectory packaging utilities.

The exact MTP path and source run are recorded in `run_summary_1p6ps.json`. This result branch does not include the MTP binary or restart binaries, so its files are intended for inspection and playback rather than standalone reruns.
