# Jeitschko–Nowotny 1967 Ti₃SiC₂ DFT acquisition

This folder records a reproducible PBE/GPAW calculation panel from the independently measured Ti₃SiC₂ structure of Jeitschko and Nowotny (1967), DOI [10.1007/BF00899949](https://doi.org/10.1007/BF00899949). The geometry is reconstructed from the full-text Fourier coordinates and assumes ideal full occupancy. It is a candidate source family, not a validation assignment.

## What was calculated

Twelve completed outputs passed `verify_result.py` and their per-result `SHA256SUMS.txt` checks: four Ti/C internal-coordinate displacement candidates at 10×10×2; two fixed-geometry controls at 10×10×2; five matching 6×6×2 calculations retained only as mesh controls; and one coarse Gamma diagnostic. All calculations used GPAW 26.7.0 / ASE 3.29.0, PBE, a 500 eV plane-wave cutoff, 0.10 eV Fermi smearing, and four MPI ranks. Energies (native extrapolated and finite-width free energy) are recorded separately from forces.

The four 10×10×2 displacement candidates are listed in `k10_batch_analysis.json` and `k10_batch_summary.csv`. The paired free-energy finite-difference force differs from the directly calculated baseline force by −0.000252 eV/Å for the Ti(4f) displacement and +0.000361 eV/Å for C(4f). The highest observed cgroup memory for those primitive-cell jobs was about 1.36 GiB.

## Mesh decision and limits

The same Fourier-coordinate geometry at 6×6×2 and 10×10×2 differs by 0.06198 eV/Å in RMS force-vector norm. The five 6×6×2 outputs are therefore labelled underconverged-mesh controls and must not be mixed into training. The abstract-coordinate 10×10×2 result is also a coordinate-variant control, not a training row. A separate 48-atom Gamma calculation is only a coarse-mesh diagnostic. A planned 48-atom 6×6×2 attempt was stopped after one SCF iteration while switching to a primitive-cell workflow; it was not an OOM event and produced no label. No `.gpw` checkpoint is archived.

These data do not supply the missing Ag–C, Ag–Si, or Ag–Ti interface labels that are relevant to the failed V19 development gates. The four displacement labels remain `train_acquisition_candidate_pending_audit`; the source-family/parent audit is still required before integration. No official dataset was changed and no model was fitted.

See `task_progress.json` for machine-readable status, `input_manifest.json` for the registered inputs and methods, and the JSON/CSV analysis files for exact values. `SHA256SUMS.txt` covers this task folder except itself.
