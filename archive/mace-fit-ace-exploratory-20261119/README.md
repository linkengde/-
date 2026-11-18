# MACE-derived Ag–Ti₃SiC₂ exploratory ACE/MTP archive

Archived on 2026-11-19 from the MACE-fitting work directory. This is an additive research archive; it does not alter or merge into `main`.

## Dataset and validation split

- Source: repository MACE v19 and v25 datasets plus 16 verified v21–v24 interface labels. No Round7 labels are used.
- 71 training frames (2,823 atoms), 7 development frames (204 atoms), and no blind test set. Exact geometry hashes do not cross between train and development.
- The dataset is exploratory and source-correlated. Read `FIT_REPORT.md` and `candidate_manifest.json` before reusing it.

## Main models

- Baseline ACE: `ace_mace_cutoff6/AgCSiTi_B1616_E10000_F1_a0.001.ace.yaml` and `.yace`. Development RMSE: 5.55 meV/atom and 0.0930 eV/Å force component; some individual development frames exceed 0.1 eV/Å.
- ACE sensitivity candidate: `sensitivity_Ti_gap_weight/…weight0.05.ace.yaml` and `.yace`. Training RMSE: 0.191 meV/atom and 0.0996 eV/Å; development RMSE: 4.55 meV/atom and 0.0555 eV/Å; worst development frame force RMSE: 0.0989 eV/Å. The `Ti_gap` training frame itself remains at 0.914 eV/Å force RMSE, and its source k-point convergence is unresolved. This is not a production-qualified model.
- MTP10 baseline files and the force-weight 0.20 variant are included with logs, predictions and scores. They do not meet the force gate. A no-`Ti_gap` MTP sensitivity result is also included; it performs worse and is not a recommended fit.

Exploratory targets used in this project: energy RMSE ≤10 meV/atom; force-component RMSE ≤0.1 eV/Å; Si hotspot-center vector error ≤0.5 eV/Å. The formal hotspot-center metric cannot be scored because hotspot IDs have not been registered. The all-Si-site vector error is only a proxy.

## Million-atom static capability check

LAMMPS `run 0` loaded the repository QSGS structure (1,021,793 atoms) with the downweighted ACE model on 4 MPI ranks in 35 s. The maximum model force was 10.336 eV/Å. The force audit found 330,452 atoms above 0.5 eV/Å; the minimum nearest-neighbor distance was 2.088 Å, with no pair below 1.5 Å. These are model forces on an unrelaxed geometry, not fit errors. No integration or heating was run. The full filtered force dump is gzip-compressed beside its summary.

The million-atom geometry itself remains at `models/ag_ti3sic2_qsgs_1to1/ag_ti3sic2_qsgs.data.gz` in the repository. Its source SHA-256 is recorded in `million_atom_static_ace/benchmark_summary.json`.

## Reproduction and integrity

The scripts preserve the original analysis workflow and may reference the original working directory `/workspace/potential-comparison/mace_fit_data_20261118`. The archive includes the inputs, model files and outputs needed to audit the reported results.

From this directory, verify the archived files with:

```bash
sha256sum -c SHA256SUMS.txt
```

The next 7–8 hour compute plan is in `NEXT_8H_PLAN.md`. The separate cloud-environment setup draft is not part of this archive.
