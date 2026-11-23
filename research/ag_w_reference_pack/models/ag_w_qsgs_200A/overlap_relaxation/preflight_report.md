# EFS static-minimization preflight

This is an exploratory numerical check of the raw QSGS geometry. It did not alter `../ag_w_qsgs.data.gz` and it did not run dynamics.

## Result

Two LAMMPS `eam/fs` conjugate-gradient minimizations were attempted with the original z = 0–20 Å atoms fixed in the first stage. Stage 1 used 100 iterations and stopped at its iteration cap. Its potential energy changed from −1,350,819.148 to −2,307,526.797 eV; the maximum force fell from 642.14 to 2.760 eV/Å. Stage 2 used the stage-1 coordinates and stopped on the energy criterion after 199 iterations, at −2,324,496.402 eV and maximum force 2.172 eV/Å. The force norm was 33.302 eV/Å. The fixed group in stage 2 was reselected by its current z coordinates and therefore contained 49,164 atoms, 35 more than the original fixed group; this stage is diagnostic, not a production preparation.

The closest Ag–W separation after stage 2 is still only 0.37 Å. There are 221 Ag atoms whose nearest W is closer than 1 Å. Thirty-five atoms moved above z = 200 Å during minimization. The remaining close contacts, force residual, and changes at the free surface mean this structure is **not suitable for MD initialization**.

## Decision

No 300 K run or TTM heating was launched. The geometry builder's raw QSGS data remain the canonical candidate; the EFS-minimized coordinates are rejected as an MD initial state. A valid next initialization needs a better controlled interface construction/overlap treatment and another static check. TTM also remains gated on W electronic-property data and a validated Ag/W grid closure; see `../../TTM_HEATING_READINESS.md`.

The two input files and LAMMPS logs are retained for audit. The large intermediate data files are intentionally omitted from this repository artifact; each was about 38 MB, and neither passed the MD-start gate.
