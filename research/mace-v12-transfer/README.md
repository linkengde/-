# Ag/Ti₃SiC₂ MACE handoff package

This repository contains the GPAW interface labels, MACE v9–v12 models and data, training/evaluation scripts, and calculation records for developing a validated Ag/C/Si/Ti potential. The Ti₃SiC₂ skeleton outline is a project constraint; the Ag fraction targets about 59.8775 wt% but may be adjusted slightly, and interface distances/model dimensions may change on documented copies.

## Current status

- The four v12 training interface labels (AgSi 2.2/2.8 Å and AgC 2.0/2.4 Å) have converged and are archived.
- The v12 dataset was built with a 4/4 energy-composition rank; training completed 80 epochs and selected the epoch-76 checkpoint.
- Three additional fixed-geometry DFT checks (AgC 2.30 Å, AgSi 2.60 Å, AgTi 2.60 Å) have also been archived.
- Current model validation is **UNDETERMINED**. v12 improves external energy errors, but contact-force behavior is mixed; Ag–C force checks worsen. Do not use v12 for production candidate long MD or TTM yet.

Read [`WORK_LOG_2026-10.md`](WORK_LOG_2026-10.md) first for the chronological record, exact findings, evidence paths, and the next-window handoff. The detailed current assessment is `periodic_interface_v4/mace_periodic_v12_interface_energy/results/v12_validation_assessment.json`.

## Next work

Continue with a separate v13 dataset and model. Preserve v12 data, checkpoint, and frozen validation/test splits; add the three verified distance-scan DFT labels to the v13 training candidates; add more independent Ag–C and Ag–Si interface environments and retain unseen holdouts. Evaluate v9–v13 per structure. Only proceed to candidate static checks and short MD if independent interface force validation supports it. TTM, stable-melt preparation, paired pressure tests, and long MD remain later gated stages.

## Environment and artifacts

Pinned setup commands and previously verified versions are in [`environment_setup.md`](environment_setup.md) and [`environment_validation.md`](environment_validation.md). The MACE-MP-0b3-medium foundation model is bundled in `periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/`. File sizes and SHA256 values are recorded in `file_manifest.json` (the manifest does not hash itself).
