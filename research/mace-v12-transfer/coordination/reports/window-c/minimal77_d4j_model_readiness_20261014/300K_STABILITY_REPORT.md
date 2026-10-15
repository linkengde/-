# 300 K short-stability report — M77 Ag/Ti₃SiC₂ candidates

## Result

Neither NVE nor NVT was run. Status is **HOLD / NOT RUN**, rather than PASS or FAIL.

The geometry files pass only a static coordinate screen. This workspace has no LAMMPS executable/library or Python LAMMPS module in its active simulation environments. In addition, the only found four-element combination is calibration-v0 EAM+Tersoff+Morse, which repository records explicitly label non-production and report as mismatching fixed-geometry DFT force responses. Running a new trajectory with that unvalidated combination would not establish a defensible D4J admission result.

The request also does not bind the short-test timestep, NVE/NVT step count, thermostat/integrator implementation, thermostat time constant, random seed, or allowed temperature/energy drift thresholds. These should be inherited from a D4J test protocol or separately frozen, not guessed in this report.

## Requested metrics

| Metric | NVE | NVT | Status |
|---|---|---|---|
| Temperature drift | — | — | Not measured |
| Maximum force | — | — | Not measured; LAMMPS `run 0` unavailable |
| Atom loss | — | — | Not measured |
| Energy change/drift | — | — | Not measured |
| Numerical stability | — | — | Not measured |

## Minimum prerequisites to resume this test

1. Identify the exact D4J-compatible four-element potential and provide its version, source/license status, type mapping, component hashes, and explicit pair definitions. Do not silently substitute calibration-v0.
2. Provide an available, versioned LAMMPS/D4J runtime and its build/package identity.
3. Freeze a short, bounded NVE/NVT protocol: timestep, step count, integrator, thermostat, damping, seed, logging cadence, and stop conditions.
4. Run `run 0` on all three coordinate candidates first; archive per-atom forces and maximum force. Stop on nonfinite forces, missing pair terms, lost atoms, or an unapproved potential identity.

No TTM heating, DFT, fitting, or potential modification was performed.

