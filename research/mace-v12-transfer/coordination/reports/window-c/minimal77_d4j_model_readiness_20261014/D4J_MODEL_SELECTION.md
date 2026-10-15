# D4J candidate selection

## Selection outcome

`D4J_MODEL_READY=NO`.

All three 77-atom coordinates are suitable **geometry candidates** for later force-field screening. None has passed LAMMPS `run 0` or 300 K NVE/NVT stability. Therefore no candidate is selected as a force-field-admitted D4J model.

For a future first diagnostic, the 2.40 Å geometry is the central, inherited reference gap and is the least extrapolative starting point. It is only a **provisional geometry choice**, not a DFT-ranked, force-ranked, equilibrated, or D4J-approved model. Retain 2.16 and 2.64 Å as fixed-geometry bracketing candidates; do not infer the best contact from distance alone.

## Current gate matrix

| Candidate | Geometry file | Geometry screen | `run 0` max force | 300 K NVE/NVT | D4J selection |
|---|---|---|---|---|---|
| gap 2.16 Å | `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p16A.extxyz` | PASS, no <1.5 Å pair | BLOCKED | NOT RUN | HOLD |
| gap 2.40 Å | `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p40A.extxyz` | PASS, no <1.5 Å pair | BLOCKED | NOT RUN | HOLD |
| gap 2.64 Å | `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p64A.extxyz` | PASS, no <1.5 Å pair | BLOCKED | NOT RUN | HOLD |

## Scope limits

The common-cell model is a 0° coherent interface with the design's shared biaxial strain; it gives up the 12.51983° low-mismatch registry. It is an inexpensive local-contact geometry only. It cannot substantiate interface energetics, temperature stability, or TTM coupling. TTM heating remains unrun and outside this task.

