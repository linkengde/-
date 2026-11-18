# POD v10 full-model FIRE: first attempt failed

更新：2026-11-18

## Setup

- 结构：347,464 atoms, 200×200×300 Å simulation box.
- Potential: POD v10.
- Fixed set: ceramic atoms in 0≤z<20 Å; all remaining atoms were mobile.
- FIRE parameter: `dmax=0.1 Å`, 1000 iteration / 5000 evaluation cap, convergence force 0.3 eV/Å.
- The original input data was read-only; no MD integration was run.

## Observed trajectory in the LAMMPS log

| FIRE step | Potential energy (eV) | Fmax (eV/Å) |
|---:|---:|---:|
| 0 | −1,744,893.56 | 51.07 |
| 10 | −1,831,995.17 | 16.06 |
| 20 | −1,995,554.27 | 56.23 |
| 30 | −2,058,638.74 | 904.04 |
| 40 | −2,114,137.94 | 252.29 |
| 160 | −2,429,384.50 | 767.39 |
| 440 | −2,920,309.65 | 566.11 |
| 480 | −2,976,757.87 | 314.80 |

At step 480 LAMMPS terminated with `Domain too large for neighbor bins`. The energy fell by about 1.232×10⁶ eV (3.55 eV/atom) and force maxima repeatedly exceeded 500 eV/Å. This attempt fails the P0 stopping condition.

## What can and cannot be concluded

- The original 347,464-atom data was not overwritten.
- This attempt produced no minimized endpoint and no coordinate checkpoints. Therefore its maximum displacement, shortest contacts, Ti–C coordination, and phase connectivity cannot be measured from a final state.
- The run is insufficient to say exactly which atoms moved or whether the ceramic network broke, but the force spikes and neighbor-bin domain error show the relaxation is numerically unstable with this FIRE step size.
- No NVE/NVT or heating run is authorized by this failed P0 result. A separate retry with `dmax=0.01 Å` and coordinate checkpoints is in progress.

Raw log: [log.lammps](/workspace/potential-fit/round7_pod_snap/models/final_v10/pod_minimize_347k/log.lammps).
