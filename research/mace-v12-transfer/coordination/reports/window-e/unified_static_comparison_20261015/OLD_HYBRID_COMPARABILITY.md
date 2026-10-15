# Old EAM + TSCBOP + Morse comparability

**COMPARABILITY = PASS for this fixed-geometry static comparison.** The source package had a complete `run 0` input, exact 135-atom LAMMPS data, pair coefficients, all potential files, environment metadata and prior output; each was SHA checked. An isolated copy was rerun as full, Ag-only, and TSC-only static calculations. The original potential files and source input were not changed. Component runs only add group deletion with `compress no`, preserving LAMMPS IDs; force dumps were mapped to the full geometry, checked by minimum-image coordinates, species, cell and PBC.

| Item | Verified value |
| --- | --- |
| LAMMPS | 2025.7.22.4.0 (LAMMPS 22 Jul 2025 - Update 4); source metadata says MPICH 5.0.2; one MPI task, OMP=1 |
| Units and pair style | `units metal`; eV, Å, eV/Å; `hybrid eam tersoff morse 8.0` |
| Type map | 1 Ti, 2 Si, 3 C, 4 Ag |
| Morse include SHA-256 | `34deced5e08174ba369290b87f5c3bfd2b284b272991adf302b2b0d28c12815f` |
| Ag EAM SHA-256 | `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807` |
| TSCBOP SHA-256 | `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` |
| Existing source structure | `0f41ad33a824beebab0cc6c2a4ded89f5508800a77172c8fd548875f6b45598d` |
| Action | `run 0`; no minimize, timestep integration, MD, or TTM |

Run-0 energies: full `-811.787971429182 eV`, Ag-only `-45.215391718684 eV`, TSC-only `-718.653197727884 eV`; composition-subtracted `Eint = -47.919381982615 eV`. The DFT total-energy increment is `-22.229525474612 eV`, so legacy calibration-v0 overbinds this geometry by `-25.690 eV`. This subtraction is meaningful because the same cell, coordinates, PBC and constituent atom IDs were used; absolute DFT and LAMMPS energies are not directly compared.

The rerun gives full-system force-vector RMSE `2.254673 eV/Å`, reproducing the archived 2.2547 eV/Å value. Effective interaction-force RMSE is `1.351321 eV/Å`; Ag normal interaction force is `-45.001092 eV/Å`, compared with DFT `-0.171011 eV/Å`. The legacy force mismatch is not solely cross-term: the TSC-only force RMSE is `2.562410 eV/Å`, while Ag-only is `0.056485 eV/Å`.

The reproduction files and raw force dumps are under `legacy_run0/`. The one initial LAMMPS invocation failed during MPI/UCX setup in the sandbox; it was rerun successfully with UCX transport selection `UCX_TLS=tcp,self`. This only selects runtime transport and does not alter force-field inputs.
