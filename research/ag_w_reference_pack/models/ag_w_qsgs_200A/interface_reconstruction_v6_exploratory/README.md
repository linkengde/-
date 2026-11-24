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

## Next cold-baseline check

The follow-on input `in.nvt_0p5ps_nve_0p05ps` continues from the 25 fs NVE endpoint. It runs 0.5 ps NVT at 300 K (`Tdamp=0.1 ps`, `dt=0.25 fs`), then removes the thermostat and runs 0.05 ps NVE. Only z<20 Å remains fixed. It has no laser source or TTM.

### Follow-up result

The 0.5 ps NVT and 0.05 ps NVE completed with 482,601 atoms retained. The mobile temperature began near 197.5 K and ended at 197.1 K after NVT; during NVT it ranged from about 86.5 to 197.1 K, so the 300 K target was not reached. W and Ag temperatures moved closer together, ending near 198 K and 197 K. NVT total energy increased by about 4.44 keV, with a potential-energy increase of about 4.47 keV and nearly unchanged kinetic energy. This indicates that structural rearrangement continued to absorb thermostat work.

The subsequent 0.05 ps NVE segment stayed near 195 K and conserved total energy to about 0.003 eV. Relative to the NVT start snapshot, 241 mobile atoms moved more than 1 Å and the maximum displacement was about 2.00 Å; no atoms were lost. Minimum pair distances changed only slightly during this segment, but remain short relative to the previous context screens. This verifies a short NVE integration near 197 K, not a 300 K equilibrium baseline. Therefore the cold-baseline stage remains incomplete; no TTM or laser source should be added yet. Detailed thermo, distance, and displacement data are in `nvt_followup_metrics.json`, with the LAMMPS log and compressed snapshots retained beside it.

### Corrected fixed-ID continuation

The preceding NVT input selected the bottom slab after loading the evolved coordinates, so 46,232 atoms were fixed instead of the original 46,117. The corrected continuation defines the fixed group from the original v5 EFS coordinates before loading the latest NVT positions and velocities. It restores those same 46,117 IDs. The 115 atoms that were accidentally frozen in the preceding NVT had zero velocities; their velocities are reinitialized at the current ~197 K mobile temperature, while all other positions and velocities are continued from the preceding NVT endpoint. The new fixed-ID run is limited to 1.0 ps NVT plus 0.05 ps NVE; it will determine whether the 300 K baseline criteria are met.

The first attempt at this corrected continuation omitted `pair_style/pair_coeff`. LAMMPS explicitly reported zero communication cutoff, zero pair energy, and zero forces. That entire output was invalid and has been excluded; its input and log are retained as `in.invalid_zero_interactions_1ps` and `log.invalid_zero_interactions_1ps`. The corrected input now loads Dai EFS explicitly and is the only run to use for analysis.

### Corrected fixed-ID continuation result

The corrected 1.0 ps NVT plus 0.05 ps NVE run retained all 482,601 atoms and kept the original 46,117 fixed IDs. The mobile temperature rose from 197.12 K to 293.79 K during NVT; at the endpoint Ag was 289.35 K and W was 301.08 K. During the NVE check, mobile temperature changed by only −0.084 K and total energy drifted by +0.00038 eV. There was no atom loss. Relative to the NVT start, 8,667 mobile atoms moved more than 1 Å by the NVT endpoint (maximum 2.85 Å); the NVT energy continued to rise as the thermostat heated and the structure rearranged. Thus this is an **exploratory cold-baseline candidate**, not a fully equilibrated or statically accepted structure. The NVT temperature slope over its final 0.1 ps was still +8.51 K.

The invalid no-interaction attempt is not used in these metrics. It remains separately archived to make the correction auditable.

### TTM program-link checks

The synthetic Ag–W TTM grid smoke check was rerun and passed: Ag-only source mapping matched its quadrature reference, Ag/W interface exchange was equal and opposite, and the combined electron/lattice ledger residual was about `1.48e-14 eV`. This check uses no LAMMPS atoms.

The finite-size LAMMPS/TTM link preflight uses the corrected NVE endpoint (`nve_fixedids_after_nvt_final.dump.gz`) as a restart, with the v5 minimized data file as the ID/type/box template. The exploratory runner mode requires a separate manifest that explicitly says `md_start_gate: EXPLORATORY_ONLY` and discloses the existing `STATIC_EFS_FAIL`; the ordinary runner remains gated on a static PASS. LAMMPS read the restart and completed the short 20-step 300 K setup segment, but the actual atom-to-grid map then stopped: peak Ag+W cell fill was 1.31585 on the selected 10 Å mesh. The synthetic grid check passes, but the model cannot be passed to the TTM bridge until this overfill is resolved without clipping phase fractions. No paired source/no-source LAMMPS/TTM run was completed. This is **not** a 300 K laser-melting or melt/quench result.

The intended software probe keeps `I0=5 eV/(ps Å²)`, `R=30 Å`, and `l_skin=20 Å`, with a 10 fs FWHM Gaussian at 20 fs followed by 20 fs of cooling. This pulse was not applied because the geometry mapping preflight failed. Atoms above z=200 Å are included in the Newton integration group, but only atoms in z=40–200 Å are mapped to the material TTM grid. The TTM phase fractions and atom maps are frozen in the short-run prototype.

The latest machine-readable cold baseline is `nvt_fixedids_metrics.json`; inputs, log, restart snapshots, and the exploratory state manifest are included alongside it. The TTM preflight report and failed LAMMPS setup logs are retained in `../ttm/AGW_TTM_LINK_PREFLIGHT_V6.md` and adjacent log files.

### v7 joint TTM mapping diagnosis

The next mapping audit uses the frozen corrected NVE endpoint and reconstructs LAMMPS's static Newton membership from the template data coordinates, before restart coordinates are loaded. At the frozen endpoint, 4,490 Newton-group atoms lie outside the TTM material interval 40–200 Å; 1,875 were already outside in the template and 2,619 moved outside during the NVT/NVE continuation. The mapper excludes them before CIC; the old z clipping would have deposited all 4,384 atoms above z=200 Å into the top grid layer. A synthetic check verifies this boundary filter.

A conservative joint Ag/W mapper now shares each cell's free capacity and carries phase volume, Ag source power, and atom-to-grid weights through the same transfer fractions. Synthetic conservation and TTM energy-ledger checks pass. The real 5, 10, and 20 Å maps still fail the configured 20 Å transfer limit; independent max-flow feasibility estimates show unresolved local overflow for each mesh. No actual TTM step, laser source, or source/no-source trajectory pair was run. Full details are in `../ttm/AGW_TTM_MAPPING_V7.md` and `../ttm/agw_ttm_mapping_v7_report.json`.
