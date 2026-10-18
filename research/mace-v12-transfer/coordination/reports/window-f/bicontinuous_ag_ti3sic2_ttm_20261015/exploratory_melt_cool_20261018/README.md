# Exploratory Ag–Ti3SiC2 relaxation and melt–cool branch

Independent working copy for the uploaded, static-audited 39,486-atom Scheme A model. The source data hash is verified against the frozen MACE comparison source. The original file under the upload/reference area is unchanged.

This LAMMPS branch preserves the existing archived calibration-v0 mapping byte-for-byte: Ag EAM + Ti-Si-C Tersoff + three Ag–TSC Morse cross terms. MACE-MP and MACE V17 remain static comparison results and are not stacked with this classical potential.

Stage 1 runs `run 0`, then bounded FIRE minimization at fixed cell with `dmax=0.05 Å`, up to 500 iterations / 5,000 force evaluations, stopping at `Fmax<=0.5 eV/Å` if reached. This is an exploratory stability screen, not a force-field validation. The input writes pre/post force dumps, minimization frames, relaxed data and restart. No temperature, heating, cooling, TTM, or DFT step is included in this stage.

## Follow-up outcome

The bounded FIRE run reached its 500-step cap at Fmax 2.450 eV/Å. A separate 300 K NVT baseline was then started on the relaxed copy and interrupted at 0.2 ps after recurring localized interface forces above 60 eV/Å. No local heat pulse, cooling, TTM, DFT, or training was run. See `EXPLORATORY_STATUS.md`, `RELAXATION_LOCAL_GEOMETRY_AUDIT.md`, and `SCREEN_300K_AUDIT.md`. `MINIMIZATION_OVITO.extxyz` and `TRAJECTORY_300K_OVITO.extxyz` are separately readable in OVITO; the thermal segment is partial and does not establish stability or melting behavior.
