# R2 Phase2 SCF recovery execution report

## Outcome

`SUCCESS`; GPAW reported `Converged in 141 steps`, and the reviewed collector returned `scf_converged=true`. The fixed-geometry R2 diagnostic single point completed in 12042.3 seconds. The result remains `REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY`; no TRAIN/DEV/TEST assignment was made.

- Total energy: `-141.959943725712 eV`
- Force-consistent free energy: `-142.135538402633 eV`
- Atomic forces: `AVAILABLE`; 20 vectors saved in `RESULTS_COLLECTED.json` (maximum force norm `1.421920285 eV/Å`)
- Stress: `AVAILABLE`; saved in the collected JSON in eV/Å³
- SCF iterations: `141`
- Runner and collector: exit status 0 / 0

## Authorization and frozen source binding

- C authorization commit: `cf6dd19fee9375a414336b008746c4881822c0b5`
- Authorization SHA-256: `a78fc33528521631f02c215c70786e49ead9ce5b11088b9351bf4631a9ab3901`
- Bound E package commit: `4d97374aefc338a714411003aa06ec35aa894611`
- Runner SHA-256: `da4b122948714b86a19ec42878f18fe24d65006238ad890383b433fd30d24a18`
- Collector SHA-256: `ca28bae1e7f2ced09b6eed161eb01ed224bfbe3b2f61fd021dc360a99dcb0bfb`
- R2 input SHA-256: `7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66`
- R2 structure SHA-256: `cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70`
- Phase2 parameter lock SHA-256: `d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78`
- Recovery config SHA-256: `c969d1cee9b8fe8014412eebbbcba4355c2a372d6baa58b5139fee2f1b252437`

The actual run used PBE/PAW, 500 eV, 8×8×1, Fermi-Dirac width 0.10 eV, Mixer beta 0.05 / nmaxold 8 / weight 100, PPCG niter 4, maxiter 200, the three locked 1×10⁻⁵ convergence thresholds, four MPI ranks and one OMP thread. Geometry, cell and PBC remained fixed. GPAW logged the exact input parameters and the 141-step convergence marker.

## Environment and preflight

Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, OpenMPI 5.0.7, NumPy 2.5.3 and SciPy 1.18.1 were used. The four Ag/Ti/Si/C PAW SHA-256 values match the authorization-bound package. The resource preflight passed with a 4-core cgroup quota, approximately 32.25 GB memory headroom, 27,939,160,064 bytes free and 2,038,282 free inodes. The preflight's resource-only JSON intentionally does not itself authorize execution; the runner's authorization gate bound the C authorization SHA above.

Two initial CLI launches exited before SCF because of GPAW CLI argument / MPI bootstrap syntax. Their stderr files are preserved. The invocation was corrected to the GPAW-supported `gpaw python --` form; the reviewed runner, configuration, geometry and scientific parameters were not changed. The actual authorized run then passed the runner gate and finished normally.

## Checkpoint record

Three periodic checkpoints were written with wavefunctions and readback-verified at completed iterations 40, 80 and 120. Their sidecars and hashes are in `checkpoints/` and `EXECUTION_MANIFEST.json`. The final SCF converged at iteration 141, so no 160-step checkpoint was due. The GPW binary payloads remain in the session workspace (each is about 1.07 GB); this GitHub report stores their verified SHA-256, byte count, and status rather than the large binaries.

## R1 comparability review

R1 and R2 both contain 20 atoms (`C4Ag7Si3Ti6`), share the same cell and PBC, and use the same PBE/PAW, cutoff, smearing, k-point mesh, convergence thresholds, GPAW/ASE versions and PAW files. The R2 recovery additionally used the C-authorized PPCG niter 4 and 200-step cap; R1's archived log does not declare an eigensolver, so the numerical solver path is not identical. Both calculations are converged, and no geometry was changed.

R1's archived 8×8×1 energy is `-188.921492255429 eV`; recovered R2 is `-141.959943725712 eV`, a difference of `+46.961548529716 eV` for the configuration. The structures share composition and cell but differ in the frozen in-plane Ag registry. This energy gap is large and needs C's review before using the pair for an energy ranking. This audit does not invalidate the R2 SCF convergence, but does not promote it as a training label or alter any role.

## Integrity

`RESULTS_COLLECTED.json` contains energy, free energy, all forces, stress, convergence state, source hashes, authorization binding and statuses. `SHA256SUMS.txt` covers every file in this report package except the checksum file itself. Checkpoint GPW payload hashes are listed separately in the manifest because those binary payloads are not committed. No R1, other registry, d0 scan, retry, or MACE training was run.
