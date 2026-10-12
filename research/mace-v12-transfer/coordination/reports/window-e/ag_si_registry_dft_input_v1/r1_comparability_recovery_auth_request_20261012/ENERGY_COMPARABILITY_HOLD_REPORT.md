# R1/R2 energy comparability review

**Current state: `ENERGY_COMPARABILITY=HOLD`.** The original R1 result used PPCG `niter=2` and a 160-step SCF cap, while recovered R2 used PPCG `niter=4` and a 200-step cap. These are solver-related differences, so the archived R1 and R2 energies are not accepted for comparative ranking.

The existing Phase2 result record has R1 converged at 152 iterations (`E=-188.92149225542852 eV`, `F=-189.16802546281826 eV`) and the original R2 Phase2 point failed at 160 iterations. The later C-authorized R2 recovery converged at 141 iterations (`E=-141.959943725712 eV`, `F=-142.135538402633 eV`) with PPCG `niter=4` and a 200-step cap. This request proposes only the matching R1 recovery point.

| Setting | Archived R1 | Recovered R2 | Requested R1 recovery |
|---|---:|---:|---:|
| GPAW / ASE | 26.7.0 / 3.29.0 | 26.7.0 / 3.29.0 | 26.7.0 / 3.29.0 |
| Method / cutoff | PBE/PAW / 500 eV | PBE/PAW / 500 eV | PBE/PAW / 500 eV |
| k-points | 8x8x1 | 8x8x1 | 8x8x1 |
| Fermi-Dirac width | 0.10 eV | 0.10 eV | 0.10 eV |
| Mixer | beta 0.05, nmaxold 8, weight 100 | same | same |
| PPCG niter | **2** | **4** | **4** |
| SCF max iterations | **160** | **200** | **200** |
| Convergence thresholds | 1e-5 (energy, density, eigenstates) | same | same |
| Geometry | fixed | fixed | fixed |

Requested R1 bindings: input SHA `d3f7b5542f071cbc432ade3277a22b0af5277cfd539d1eb48a0d006e7250ca6b`; structure SHA `b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890`; Phase2 lock SHA `d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78`.

If C authorizes and the new R1 converges with the bound settings, the two runs will have matching solver settings and can be reviewed for energy comparability. This request does not authorize ranking, d0 scans, or role assignment.
