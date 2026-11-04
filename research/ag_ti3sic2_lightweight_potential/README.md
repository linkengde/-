# Ag–Ti₃SiC₂ exploration: potential and data inventory

**Status (2026-11-04):** The original Ag EAM + Ti–Si–C TSCBOP/Tersoff + Ag–framework Morse combination may be used for clearly labeled, qualitative exploratory runs. The documented full-size starting geometry has not passed a static force/relaxation screen, so it should first be equilibrated and checked for nonphysical collapse, atom loss, and loss of Ti₃SiC₂ connectivity. This is a model-specific limitation, not a requirement that every exploratory study reach production-level DFT agreement. Keep the user-selected `0.3 eV/Å` maximum per-mobile-atom force-vector target for the initial geometry relaxation; do not conflate that geometry criterion with force-error accuracy of the potential.

This record inventories data and literature relevant to a CPU-friendlier, single four-element Ag–Ti–Si–C potential. The preferred first candidate is LAMMPS's native POD model (`pair_style pod`, trained by `fitpod`); linear SNAP is the comparison candidate. Either route can replace the old Ag EAM, TSCBOP/Tersoff, and Morse stack with one model. No POD or SNAP coefficients have been fitted and no new atomistic labels were acquired in this update. The cited papers are literature leads, not verified training datasets: raw configurations with DFT energies and per-atom forces still need to be found in their supplements/repositories or generated and DFT-labeled.

## Existing repository data

The V19 train55 source contains 55 DFT-labeled frames: 4 pure Ti₃SiC₂, 51 Ag-containing, and no pure Ag frames. A separate static audit has six 32-atom FCC-Ag controls, five 48-atom Ti₃SiC₂ controls, and three 52-atom Ag–Ti₃SiC₂ interface controls. The three interface controls contain four Ag atoms each and have already served as monitored development structures; if reused to fit a new SNAP model, they cannot also be presented as independent validation.

Sources:

- [V19 train55 cold-gate inventory](../mace-v12-transfer/periodic_interface_v4/v19_train55_cold_gate1_20261104/README.md)
- [Static audit and source paths](../mace-v12-transfer/periodic_interface_v4/v19_train55_cold_gate1_20261104/results/static_audit.json)
- [Train55 dataset manifest](../mace-v12-transfer/periodic_interface_v4/v19_train55_cold_gate1_20261104/provenance/train55_dataset_manifest.json)
- [Original hybrid-potential relaxation diagnostics](../../models/ag_ti3sic2_qsgs_1to1_finer/STAGE1_MINIMIZATION_DIAGNOSTICS.md)

The original hybrid-potential full-model audit reported an initial maximum force-vector norm of `3.9878 eV/Å`. Tested minimizations did not reduce all mobile-atom forces below `1 eV/Å`; some paths also developed short Ag–C, Ag–Si, Ti–C, or C–C contacts. Thus the old potential is appropriate only for qualitative exploration with its limits stated, and the specific starting geometry needs a basic stability/equilibration check before interpreting a heat-and-cool trajectory.

## Literature leads

| Source | What it contributes | Training-data status |
| --- | --- | --- |
| Guarini et al. (2020), [“Neutron Brillouin scattering and ab initio simulation study of the collective dynamics of liquid silver”](https://doi.org/10.1103/PhysRevB.102.054210), *Physical Review B* 102, 054210 | Liquid Ag neutron measurements and AIMD-derived dynamic structure factor; the highest-priority paper to inspect for AIMD settings and supplementary/repository data. | Article is identified; raw per-frame structures, energies, and forces have not been verified. |
| Baty et al. (2021), [“Ab initio phase diagram of silver”](https://doi.org/10.1088/1361-648X/ac23fb), *Journal of Physics: Condensed Matter* | Ab initio quantum-MD melting curves and thermal equation of state, especially relevant to high-temperature/high-pressure Ag behavior. | Optional reference; raw force-labeled frames have not been verified and the pressure range may not match the present melt-pool application. |
| Kearley et al. (2014), [“Inelastic Neutron Scattering and Density Functional Theory–Molecular Dynamics Study of Si Dynamics in Ti₃SiC₂”](https://doi.org/10.1111/jace.12725), *Journal of the American Ceramic Society* 97(3), 916–922 | DFT–MD finite-temperature Si dynamics in the Ti₃SiC₂ framework. | Useful low-temperature framework reference; no raw trajectories/force labels verified. |
| Gertzen et al. (2021), [“DFT Study of MAX Phase Surfaces for Electrocatalyst Support Materials in Hydrogen Fuel Cells”](https://doi.org/10.3390/ma14010077), *Materials* 14(1), 77 | DFT comparison of (001) terminations including Ti₃SiC₂; useful for choosing slab terminations. | Surface structures/energies are not Ag-contact interface labels; no force-labeled training set verified. |
| Costa Pereira & Kenisarin (2022), [“A review of metallic materials for latent heat thermal energy storage: thermophysical properties, applications, and challenges”](https://doi.org/10.1016/j.rser.2021.111812), *Renewable and Sustainable Energy Reviews* 154, 111812 | Table 5 reports the liquid-Ag density fit range (1235–1600 K), `ρ₀ = 9294 kg/m³`, `c = 0.877 kg/(m³·K)`, and `Tref = 1234.9 K`; useful for validating predicted liquid density. | Macroscopic property compilation only; no atomic coordinates or DFT forces. The extracted `c₁` column remains ambiguous and must not be assumed zero. |

No directly relevant Ag/Ti₃SiC₂ DFT interface paper has been identified in this search. Existing repository interface labels can be reused as training for a new model, but larger and hotter Ag–surface configurations are needed for transfer to a melt pool.

## Uploaded paper review

The two user-supplied PDFs were inspected locally and are not copied into this repository. Their checksums identify the exact reviewed files:

| Supplied file | SHA256 | Finding for potential work |
| --- | --- | --- |
| `wzAbPm-1793787202933(科研通-ablesci.com).pdf` — Kearley et al., DOI `10.1111/jace.12725` | `0b4e7f6597141884f8331593a54a9b61c7e4971575c92640d5995315232037e9` | Includes a 108-atom, 3×3×1 Ti₃SiC₂ VASP PAW-PBE MD calculation at 70 K. The MD used a 300 eV cutoff, 1 fs steps, 4 ps isokinetic equilibration, 4 ps NVE equilibration, and 18 ps NVE production. The paper analyzes 18,000 production steps, but the supplied PDF contains no per-frame coordinates, energies, or forces. It is not a liquid/high-temperature or Ag-interface dataset. The trajectory would need to be requested from the authors and, even if obtained, would only add low-temperature framework configurations. |
| `materials-14-00077-v2.pdf` — Gertzen et al., DOI `10.3390/ma14010077` | `3c0dfdf8193f89a6dae9c0da56fdc55c98e7f734277f098c7baca75c2f3307b3` | Gives VASP PAW-PBE (001) surface calculations, a 0.03 eV/Å force convergence threshold, and Ti₃SiC₂ cleavage-energy pairs. The lowest reported Ti₃SiC₂ pair is Si(Ti₂)/Ti₂(C), 2.802 eV/unit cell. Its data statement says study data are available from the corresponding author on request. The supplied article has no Ag contact calculations or raw per-atom force-labelled structures. Request the relaxed slab geometries and VASP outputs if these are to seed surface training. |

The accompanying `materials-14-00077-s001.pdf` supplement supplies the following useful **bulk-only** Ti₃SiC₂ reference data:

- Table S1 initial hexagonal cell: `a = 3.0665 Å`, `c = 17.671 Å`, `V = 143.906 Å³`; fractional positions Ti(1) `(0,0,0)`, Ti(2) `(1/3,2/3,0.128)`, Si `(0,0,0.25)`, and C `(1/3,2/3,0.564)`.
- Table S2 lists a PBE bulk k-point grid of `3×3×7`; Table S3 gives a PBE bulk cutoff of `350 eV`; Table S4 lists manually fitted PBE lattice parameters `a = 3.0665 Å`, `c = 17.8651 Å`.
- Supplement checksum: `f90017ba6c5741ddf8835bdca348970bd4e1452aa50f61188098a9a71cf6e328`.

These are reproducibility/initialization data, not force labels. The supplement has no relaxed (001) slab coordinates, Ag contacts, or per-configuration energies and forces. Its settings are bulk PBE settings and should not be mixed as training labels with the repository's calculations until PAW datasets and all numerical settings are reconciled.

The Kearley article is published in the 2014 journal volume (although it was received/accepted in 2013); use 2014 in citations. The Gertzen paper's supplementary file was not among the supplied PDFs. Its article lists supplementary figures/tables for MAX unit cells, initial internal coordinates, bulk k-point grids and cutoffs, and lattice parameters; these can assist reproduction but are not a substitute for Ag-interface force labels.

Neither supplied PDF fills the priority gap for liquid Ag or hot Ag/Ti₃SiC₂ interfaces. The liquid-Ag AIMD paper by Guarini et al. remains the key source to inspect for an obtainable trajectory, while the actual interfacial configurations will likely need to be generated and DFT-labeled.

## Data needed for a SNAP fit

For each configuration, retain the cell, periodic-boundary flags, atom species and coordinates, DFT total energy, and per-atom DFT forces. Virial/stress is useful for fitting volume response. Preserve the exact exchange-correlation functional, pseudopotential identifiers/hashes, plane-wave cutoff, k-point mesh, smearing, and energy convention. Use independent train/validation/test splits by structure family, not adjacent frames from one trajectory.

Priority additions are:

1. Decorrelated liquid-Ag AIMD frames across roughly 1250–1600 K, including coordinates, energies, and forces.
2. High-temperature/distorted Ti₃SiC₂ framework configurations that retain and probe its layered bonding.
3. Ag on C-, Si-, and Ti-terminated Ti₃SiC₂ surfaces with varied registries/separations, followed by thermalized liquid-Ag contact frames.
4. Undercooled and early-solidification configurations for Ag near the framework, with an independent DFT holdout for each family.

Experimental density, lattice parameters, and surface energies are valuable validation targets but cannot replace the per-configuration DFT energy/force labels. MACE-predicted labels or trajectories are not DFT labels; a MACE trajectory can only propose configurations for subsequent DFT single points.

## Large-system fitting options

| Route | Local support | Fit/runtime characteristics | Current assessment |
| --- | --- | --- | --- |
| LAMMPS POD | `ML-POD`, `pair_style pod`, and the `fitpod` command are present in the installed LAMMPS executable. | One compact descriptor/coefficient model for all four elements; training can read extxyz DFT configurations. The pair style evaluates local environments. | First candidate to benchmark. LAMMPS documents that POD descriptor-output computes create large per-atom arrays; avoid those output computes on the production-sized structure. Benchmark `pair_style pod` itself before scaling up. |
| Linear SNAP | `ML-SNAP` and `pair_style snap` are present. | Linear coefficients and local bispectrum calculations; no neural-network inference stack. | Useful cross-check. FitSNAP 2.0.1 is in a temporary environment but currently hits Python 3.12's removed `distutils.util`; no fit has been run. |
| ACE/PACE | `pair_style pace` is not present in the installed executable. | A compact ACE model may be efficient for large systems, but requires another fitting/runtime build. | Keep as a later alternative if POD/SNAP benchmarking is insufficient. |
| MACE | Existing candidate exists but is not selected for the hundreds-of-thousands-atom production route on this CPU/memory budget. | The current serial CPU path is too costly at full size. | Do not use for the requested large-scale route. |

The runtime memory and speed of a new POD/SNAP model are not yet measured on the QSGS system. The safe scale-up sequence is a small labelled fit, held-out force/structure checks, 836–10,000-atom LAMMPS benchmarks, then gradual scale-up with measured peak memory and atom-steps/s. A few hundred thousand atoms is not cleared until those benchmarks pass.

## Software check

The local `/workspace/.venvs/lammps/bin/lmp` reports `ML-POD`, `ML-SNAP`, and `ML-UF3`, and lists `fitpod`, `pod`, `snap`, and `uf3`. No four-element POD, SNAP, or UF3 parameter file exists in the repository. FitSNAP 2.0.1 was installed in a temporary environment for a feasibility check, but import currently fails under Python 3.12 because its code imports the removed `distutils.util`; no SNAP fitting was run. This is a setup limitation, not a result about SNAP accuracy or performance.

## Machine-readable inventory

See [dataset_inventory.json](dataset_inventory.json) for the current counts, identified references, and missing label families.
