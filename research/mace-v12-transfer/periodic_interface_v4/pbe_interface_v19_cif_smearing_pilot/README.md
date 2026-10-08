# Controlled Fermi-width pilot

Same frozen52atom input and5x5/6x6 meshes; width0.05eV is the only changed calculator setting versus the existing0.1eV baselines. Numerical-only; no model/data integration. Historical archives remain untouched. Each exact label binds its new owner. Root SHA256 inventory guards queue execution; method width is read by calculator and verified from full result metadata. Four MPI ranks, ScaLAPACK and one thread per rank;3GiB startup/1.5GiB runtime disk reserves. Existing/unfinished jobs stop for review.

A exact command: `python3 research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_smearing_pilot/run_queue.py window-a AgSi_COD9009647_pilot_k5x5_sigma0p05`. B substitutes `window-b AgSi_COD9009647_pilot_k6x6_sigma0p05`. Activate the verified GPAW environment first.

Compare each new label with its `baseline_label`, and compare the two new labels using `compare_mesh_results.py FIRST SECOND`. Baseline manifest is pinned. Same unchanged numerical gates; narrower width need not improve convergence and no cause is presumed. Vacuum/dipole and independent common-target labels remain downstream.
