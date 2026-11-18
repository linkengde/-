# Follow-up package

This directory records the completed k-point convergence check, candidate refits, grouped validation, LAMMPS loading checks, and the next DFT acquisition shortlist. Start with `REPORT.md`.

## Contents

- `DFT/`: converged 1×6×3 frame, GPAW output, and mesh comparison. Large `.gpw` restart files were not copied.
- `candidate_data/`: training/development split with the Ti-gap label replaced by the 1×6×3 result; original frozen data remain at the parent archive root.
- `ACE/` and `MTP/`: fitted models, metrics, grouped fold models, and training logs.
- `LAMMPS/`: inputs and logs from two ACE/PACE and three MTP/MLIP `run 0` loading checks.
- `acquisition/`: 81 spatially separated high-force centers to guide future local-structure extraction and DFT labeling. These are not labeled configurations.
- `unmatched_source_input/`: user-supplied contact-edge table and a compatibility assessment. It does not match the archived million-atom target's atom-ID element map and is not part of the fit.
- `pipeline/`: scripts used in this compute window.

## Integrity and reproducibility

From the archive root run `sha256sum -c SHA256SUMS.txt`. The scripts preserve the tested invocations and paths for the current cloud image; executable installation paths may need adjustment on another machine. The candidate data and model/metric files are included for independent audit without rerunning DFT or fitting.
