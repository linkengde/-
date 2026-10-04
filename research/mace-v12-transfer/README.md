# Ag/Ti₃SiC₂ MACE handoff package

This repository contains the GPAW interface labels, MACE v9–v13 models and data, training/evaluation scripts, and calculation records for developing a validated Ag/C/Si/Ti potential. The Ti₃SiC₂ skeleton outline is a project constraint; the Ag fraction targets about 59.8775 wt% but may be adjusted slightly, and interface distances/model dimensions may change on documented copies.

## Current status

- The four v12 training interface labels (AgSi 2.2/2.8 Å and AgC 2.0/2.4 Å) have converged and are archived.
- The v12 dataset was built with a 4/4 energy-composition rank; training completed 80 epochs and selected the epoch-76 checkpoint.
- Three additional fixed-geometry DFT checks (AgC 2.30 Å, AgSi 2.60 Å, AgTi 2.60 Å) have also been archived.
- v13 data and training are complete: 26/2/3 split, 22 energy-labelled training frames, 80 epochs, selected epoch 72. Do not repeat v12/v13 training.
- Current model validation is **UNDETERMINED**. v13 improves related midpoint checks but lacks sufficient new-environment validation. Two additional v13 geometry DFT checks are running in window A; long MD/TTM remain later stages.

Read [`WORK_LOG_2026-10.md`](WORK_LOG_2026-10.md), including its latest continuation, for findings and evidence. The latest assessment is `periodic_interface_v4/mace_periodic_v13_interface_energy/results/v13_validation_assessment.json`.

## Next work

Use [`coordination/README.md`](coordination/README.md) and the per-window task JSON files before starting work. Window A owns the running DFT checks, two future v14 blind labels and all model training/evaluation. A separate cloud window B claims four new DFT labels and publishes them one at a time; its entry point is [`coordination/START_WINDOW_B.md`](coordination/START_WINDOW_B.md). Each window has separate result directories. The guarded v14 builder, training and comparison entry points are in `periodic_interface_v4/mace_periodic_v14_interface_energy/`; they require completed labels and a scored v13 evaluation. See [`coordination/ITERATION_PROTOCOL.md`](coordination/ITERATION_PROTOCOL.md) for the agreed screening thresholds and v14/v15/v16 data policy. Production candidate static checks, short stability tests, TTM and paired pressure simulations follow sufficient independent validation.

## Environment and artifacts

Pinned setup commands and previously verified versions are in [`environment_setup.md`](environment_setup.md) and [`environment_validation.md`](environment_validation.md). The MPI setup helper for the auxiliary window is `coordination/install_parallel_dft_environment.sh`. The MACE-MP-0b3-medium foundation model is bundled in `periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/`. Stable artifact sizes and SHA256 values are recorded in `file_manifest.json`; the manifest excludes itself and the mutable `coordination/tasks/*.json` status files. New auxiliary outputs have their own archive verification until window A incorporates them into the global manifest.
