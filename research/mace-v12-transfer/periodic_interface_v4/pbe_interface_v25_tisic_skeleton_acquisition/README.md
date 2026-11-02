# Ti₃SiC₂ Ti–Si–C skeleton-potential acquisition v25

This is a separate, non-official Ti–Si–C-only data acquisition for a host-skeleton potential. It must not reuse Ag-containing interface force labels. All new fixed-cell labels use the Fourier-coordinate Jeitschko–Nowotny 1967 structure, PBE/GPAW, 500 eV plane waves, 10×10×2 k points and 0.10 eV Fermi smearing. Energy targets for the skeleton fit will use the finite-width free energy so they are consistent with the forces; native extrapolated energies will also be retained. Stress is recorded for diagnosis but is not yet included in the fit recipe.

The existing candidate pool has five relevant records: one pristine 10×10×2 control and four paired Ti(4f)/C(4f) z-displacements. This task adds 18 configurations: paired displacements on Ti(2a), Si(2b), Ti(4f) and C(4f); four ±1% basal/c-axis strain controls; and six multi-atom random displacement snapshots. The random snapshots are split by seed into 2 training and 4 development structures. The 0.07 eV/Å force-vector RMSE working threshold is fixed before fitting; it applies per held-out skeleton configuration, not to every individual atom's maximum error.

The current run only prepares and acquires data. The fit remains disabled until every registered DFT archive verifies and a separate Ti–Si–C-only dataset/split passes the provenance and duplicate audit. No official training data, cross-potential data, sealed test labels, or production model is modified.

## Run state

`sequence_status.json` is the operational status record. DFT jobs run sequentially under the shared lock, four MPI ranks, a 32 GiB cgroup ceiling and a 30 GiB stop guard. Checkpoints are written every 20 SCF iterations and only a matching input-hash checkpoint may be resumed. The run does not push partial or running files. Verified results will be archived and uploaded after the batch is complete.
