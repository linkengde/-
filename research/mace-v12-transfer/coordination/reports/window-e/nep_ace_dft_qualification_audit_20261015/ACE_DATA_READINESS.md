# ACE data readiness

## Data suitability

ACE can consume multi-element energy/force datasets through a supported fitting workflow and deploy through a compatible LAMMPS `pace` implementation. The archive has a useful starting set of confirmed PBE energy/force records and a stronger static Ag–Ti/Si/C interface slice than it has liquid or high-temperature data. The same gaps that block NEP also block ACE: no labelled Ag liquid or liquid-solid contact family, insufficient independent clean surface/Si/elastic response, incomplete provenance for 21 geometries, and only 5 explicit stress/virial records among 71 confirmed geometries.

The raw extxyz can be normalized for an ACE workflow after checking total versus free energy, units, element ordering and parent-family split. Do not pool energies from different PAW/k-point/smearing settings blindly. Keep whole families together across TRAIN/DEV/TEST.

## Compute and deployment

No PACEMAKER/ACE training tool or LAMMPS executable is installed in the current audit environment. Its 4-CPU quota and 32 GiB memory limit allow small data/software checks only; no runtime benchmark for the target approximately 700,000 atoms is available. CPU-only LAMMPS deployment may be practical but performance depends strongly on ACE basis/cutoff and neighbor count. It must be benchmarked on the target hardware before committing to that scale. Do not infer 700k feasibility from small static DFT models.

## Decision

`ACE_TRAINING_DATA_READY = FAIL`. Data qualification fails before algorithm comparison. ACE is a conditional first route for the current CPU/LAMMPS-oriented workflow once labels are harmonized and missing physical families are added. No ACE coefficients or model were generated.
