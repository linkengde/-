# Ti₃SiC₂ Ti–Si–C skeleton-potential acquisition v25

This is a separate, non-official Ti–Si–C-only candidate for a host-skeleton potential. It does not reuse Ag-containing interface force labels. All new fixed-cell labels use the Fourier-coordinate Jeitschko–Nowotny 1967 structure, PBE/GPAW, 500 eV plane waves, 10×10×2 k points and 0.10 eV Fermi smearing. The fitted energy target is the finite-width free energy, consistent with the forces; native extrapolated energies are retained separately. Stress is recorded for diagnosis but is not included in the fit objective.

The existing candidate pool has five relevant records: one pristine 10×10×2 control and four paired Ti(4f)/C(4f) z-displacements. This task adds 18 configurations: paired displacements on Ti(2a), Si(2b), Ti(4f) and C(4f); four ±1% basal/c-axis strain controls; and six multi-atom random displacement snapshots. The random snapshots are split by seed into 2 training and 4 development structures. The 0.07 eV/Å force-vector RMSE working threshold is fixed before fitting; it applies per held-out skeleton configuration, not to every individual atom's maximum error.

The acquisition and fit are complete. All 18 new DFT archives passed verification. Together with five previously verified Ti–Si–C labels, the candidate dataset has 19 training and 4 development frames. No official training data, Ag cross-potential data, sealed test labels, or production model was modified.

The user-selected skeleton acceptance standard is per-development-configuration force-vector RMSE ≤0.07 eV/Å. The model was fine-tuned from MACE-MP-0b3-medium for 80 epochs with seed 45, CPU, four threads, batch size 1, learning rate 1e−4, weight decay 5e−7, energy weight 100 and force weight 1000. The fixed epoch-79 export and selected checkpoint are tensor-identical. The four development results are:

| Development configuration | Force-vector RMSE (eV/Å) | ≤0.07 standard | Additional ≤0.05 screen |
| --- | ---: | --- | --- |
| Random displacement, RMS 0.02 Å | 0.00913 | Pass | Pass |
| Random displacement, RMS 0.03 Å | 0.01240 | Pass | Pass |
| Random displacement, RMS 0.045 Å | 0.01854 | Pass | Pass |
| Random displacement, RMS 0.06 Å | 0.03295 | Pass | Pass |

Aggregate force-vector RMSE is 0.02041 eV/Å. All four configurations also pass the stricter 0.05 eV/Å per-configuration RMSE screen; this is recorded as a secondary comparison and does not replace the user-selected 0.07 acceptance standard. The largest single-atom vector error is 0.08529 eV/Å on the last configuration; neither RMSE threshold is an atomwise maximum bound. Status: `SKELETON_ACCEPTANCE_PASS`. The four development structures are not an independent blind test; this result applies only to the Ti–Si–C skeleton candidate and does not certify transferability to Ag-containing interfaces or general production use.

The crystallographic reconstruction assumes full occupancy; the source literature discusses possible nonstoichiometry. The candidate therefore inherits that structural assumption and covers only the sampled primitive-cell distortions and strains.

## Run state

`sequence_status.json` records acquisition completion. DFT jobs ran sequentially under the shared lock, four MPI ranks, a 32 GiB cgroup ceiling and a 30 GiB stop guard. Checkpoints were written every 20 SCF iterations and only a matching input-hash checkpoint could be resumed. The epoch-79 checkpoint and exported model, training logs, data, per-label archives and acceptance report are retained for audit.

## Artifacts

- `dataset/`: the 19-frame train file, 4-frame development file and dataset manifest with hashes.
- `calculations/<label>/`: verified fixed-geometry PBE energies, forces, stress diagnostics, SCF logs and per-label checksums for all 18 new jobs.
- `run/models/`: exported epoch-79 model and compiled model.
- `run/checkpoints/`: epoch-79 checkpoint and selected model; 79 intermediate epoch files were removed only after full 0–79 coverage and model-identity verification.
- `run/results/skeleton_fit_verification.json`: per-configuration 0.07 eV/Å acceptance result and model/checkpoint identity record.
- `SHA256SUMS.txt` and `ARTIFACT_SHA256SUMS.txt`: source/input and complete artifact inventories.
