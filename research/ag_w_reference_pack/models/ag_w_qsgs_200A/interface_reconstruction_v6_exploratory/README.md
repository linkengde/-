# Ag–W exploratory short NVE pilot

This branch tests the relaxed exploratory criteria requested for the Ag–W model. It does not modify the raw QSGS model or the v1–v5 audit files.

## Setup

- Start geometry: v5 local-normal/same-phase precondition followed by the finite Dai EFS minimization stored at `../interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz`.
- The v4 raw structure with a 0.08 Å Ag–W contact is not used.
- Boundary: `p p f`; only atoms initially at `0 ≤ z < 20 Å` are fixed. All other Ag and W atoms move.
- Dynamics: 300 K Maxwell velocities, no source and no thermostat, 100 NVE steps at 0.25 fs (0.025 ps).
- Potential: Dai 2006 Ag–W EFS in LAMMPS `eam/fs` form.

The v5 starting structure did not pass the earlier static research gate: the finite EFS minimization hit its iteration cap, left a 0.381 eV/Å maximum moving force, and displaced 1,871 atoms slightly above the nominal z=200 Å surface. Its minimum same/cross-phase distances were 2.295–2.440 Å, so this pilot tests a finite, preconditioned structure rather than the raw 0.08 Å overlap. The exploratory criterion does not require full force convergence, but all deviations remain reported.

## Scope

This is a numerical and qualitative motion screen only. A completed run, bounded atom motion, and a reasonably conserved NVE energy would show that a short trajectory can be attempted from this processed candidate. They would not validate the EFS potential, interface thermodynamics, solid-state stability, melting, or TTM heat transfer. The Ag/W morphology is based on a volume-preserving regularization of the QSGS mask; the source QSGS files remain unchanged.

## Reproduction

Run from the Ag–W model directory:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/ag-qsgs/bin/lmp -sf omp \
  -in interface_reconstruction_v6_exploratory/in.nve_pilot_0p025ps \
  -log interface_reconstruction_v6_exploratory/log.nve_pilot_0p025ps
```

The input writes initial and final atom dumps and the log. The bottom layer is a mechanical anchor; it is not a thermostatted heat sink.

## Result

The 25 fs NVE run completed for 482,601 atoms. Only the bottom 0–20 Å group (46,117 atoms) was fixed; the other 436,484 atoms moved. Atom count was unchanged, no atom moved more than 1 Å, and the largest displacement was 0.207 Å for Ag (0.130 Å for W). Total energy drift was zero at the log's 0.1 eV display resolution, so there was no numerical energy blow-up.

This passes the **level 1 program screen**. It does not produce a stable 300 K cold baseline: mobile temperature fell from 300 K to 197.5 K (W 300→160 K; Ag 300→220 K), while potential energy rose by 5,780.7 eV and kinetic energy fell by 5,780.6 eV. Maximum moving force rose from 0.381 to 3.40 eV/Å and RMS force from 0.0031 to 0.664 eV/Å. Minimum distances changed from 2.440 to 2.410 Å (W–W), 2.295 to 2.271 Å (Ag–Ag), and 2.434 to 2.405 Å (Ag–W). These values are recorded as diagnostics, not used as a paper-level pass/fail gate. The large PE-to-KE exchange shows ongoing structural relaxation, so the temperature evolution is not evidence of Ag/W heat transfer.

**Decision:** the finite, preconditioned structure can execute a short NVE trajectory under the requested bottom-only constraint, but it is not a cold equilibrium start. Stop before adding a TTM source, laser pulse, or interpreting the trajectory as melt/cool behavior. The source QSGS files remain untouched; this pilot uses the v5 volume-preserving regularized mask, which reassigned 43,054 voxels (8.41%). It therefore preserves the large-scale bicontinuous morphology but is not voxel-identical to the original mask.

Machine-readable metrics: `nve_pilot_metrics.json`. Input, outputs, source geometry, and potential hashes are listed in `SHA256SUMS.txt`.
