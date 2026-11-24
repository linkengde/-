# V11 execution plan — build a usable Ag–W exploratory model

## Why this is the next move

V10 has narrowed the problem. The remaining 5 Å unassigned volume is small (1,181.24 Å³), while the largest-overflow cells listed for the best phase sit at the two z edges of the old TTM slab (z centres 42.5 and 197.5 Å). The old map begins at z=40 Å and ends at z=200 Å, so CIC contributions near both open boundaries can be forced into edge cells. V10 also excluded 7,619 atoms from the static Newton group because their evolved coordinates were outside that old slab.

This is a testable boundary-domain hypothesis, not yet a proven diagnosis. First add guard cells and include the material that is already present. Do not start by deleting atoms.

## A. Candidate A — padded physical mapping domain (first priority)

Frozen starting state:
- data: V5 `efs_minimized.data.gz`, SHA-256 `29bf0eccea5516f1e51e7454876d7f66a99fd3b4fcdb3905b2e0e27993d3bde8`
- coordinates: V9 `v9_nve2_end.dump.gz`, SHA-256 `6de2e38f0b3c01e572e90847894f6f231f18b27fea0121f06fd12c3923d6d3c2`
- original static Newton group: 388,274 atoms, as defined by the frozen V5 data and V9 manifest.

Mapping candidate:
- x/y: 200 × 200 Å periodic.
- z map: 35–205 Å (34 cells × 5 Å).
- grid: (40, 40, 34), spacing (5, 5, 5) Å.
- phases: (0,0), (0,2.5), (2.5,0), (2.5,2.5) Å.
- transfer cap: 20 Å, direct donor-to-original-free-receiver maximum flow, no cascading.
- include all original static-Newton atoms whose current z lies inside [35,205] Å. Count atoms below/above; if any are outside, stop and adjust the guard domain from measured extrema rather than clipping.
- keep Ag/W atomic volumes unchanged.
- map the full Ag and W phases together under one cell capacity.

Required outputs:
1. raw CIC fields and combined fill maps for each phase;
2. per-z-layer overfill histogram and list of top donor cells;
3. direct-flow status, unresolved volume, actual transfer edge distances;
4. Ag and W volume conservation;
5. atom-quadrature and grid-integrated source power;
6. real LAMMPS mapping-only result for the best phase.

Important: z=35–40 and z=200–205 are guard cells. They are part of the numerical mesh, not extra material. Their phase fractions must be determined only from atom mapping. They may not be filled with artificial material.

## B. Source handling on the padded mesh

The current source definition uses:
- I0 = 5 eV/(ps Å²)
- R = 30 Å
- l_skin = 20 Å
- centre = (100,100) Å
- nominal surface plane = z=200 Å
- Ag receives direct source; W receives zero direct source.

For this first exploratory check, retain z_surface=200 Å and use the existing explicit above-surface clamp for atoms with z>200 Å (depth=0). Record how many Ag atoms are clamped. This is a nominal flat-plane approximation only; it is not a moving-free-surface model. Do not set z_surface=205 Å merely because the padded mesh ends there.

If the final selected structure has significant roughness, add a separately audited local surface-height source as a later improvement. It must conserve the same atom-integrated source power.

## C. Candidate B — small lateral dilation only if Candidate A fails

Test scale factors s = 1.0025, 1.0050, 1.0075 in x/y:
- new Lx = Ly = 200s Å;
- scale x/y coordinates about (100,100) Å;
- keep z coordinates, z box, fixed IDs, and atom counts unchanged;
- recompute grid spacing sx=sy=Lx/40; use sz=5 Å and 34 z cells over 170 Å. If the production entry point currently accepts only one scalar grid size, extend it to accept a three-component spacing tuple and test the anisotropic finite-volume operator.
- scan all four x/y phases at each scale with the same 20 Å direct-flow rule.
- after a candidate maps, minimize with the actual EFS potential before cold equilibration; report the induced in-plane strain and resulting phase densities.

Do not scale z or reduce the vacuum gap. Do not alter atomic volumes to obtain a pass.

## D. Candidate C — tiny mass-balanced atom adjustment only if A and B fail

Maximum adjustment: 0.1% of the full system (482 atoms). Prefer the minimum possible number.
- identify residual donor cells from the failed production flow;
- rank atom candidates by their actual CIC volume contribution to unresolved donor cells;
- choose a mass-balanced integer set of Ag/W removals (or additions only if a genuine local vacancy is identified);
- keep Ag/W mass ratio within 0.1% of unity;
- preserve both phases' percolation and interface connectivity;
- record every changed atom ID and local environment;
- re-run minimization, all mapping phases, and topology checks.

Do not remove atoms merely because they are Ag or W rich; preserve local topology. Do not claim a pass from a partial max-flow assignment.

## E. Cold baseline rebuild

Do not inherit the V9 phase-specific velocity bias. Use the selected minimized geometry and create a fresh common Maxwell distribution for all mobile atoms:
`velocity mobile create 300.0 <seed> mom yes rot no dist gaussian`.

Use one NVT thermostat on all mobile atoms:
- timestep 0.25 fs;
- NVT 2 ps initially (8,000 steps);
- if the final 0.5 ps still has a significant potential-energy slope or phase-temperature separation, extend in 1 ps increments up to 5 ps;
- report Ag, W, all-mobile, and interface-cell temperatures;
- report both ordinary group temperatures and explicitly COM-removed diagnostics;
- do not use separate Ag and W thermostats to force equality;
- after NVT, run 0.1 ps NVE (400 steps).

Before velocity creation, remove inherited COM drift and zero fixed-bottom velocities. Keep the original 46,117 fixed IDs unless the selected geometry makes the bottom plane adjustment unavoidable.

## F. TTM gates

After A/B/C mapping and cold gates pass:
1. freeze one accepted NVE restart;
2. run mapping-only on the exact restart used for MD;
3. run a no-source TTM pair;
4. run the same restart with a low-energy Gaussian pulse;
5. start with the existing exploratory pulse envelope (FWHM 0.01 ps, centre 0.02 ps, cooling 0.02 ps) only as a code smoke test, not an experimental reproduction;
6. verify total energy, Ag/W source partition, electron–phonon handoff, and interface exchange;
7. retain the provisional 1 GW m⁻² K⁻¹ Ag/W electronic interface conductance as a labeled assumption and do not call it measured.

A successful smoke test authorizes only short exploratory TTM-MD. It does not validate quantitative melt-pool predictions or a moving solid/liquid interface.

## G. Stop conditions and reporting

Stop and report FAIL if:
- any atom is silently clipped;
- any cell exceeds total capacity after the production map;
- transfer exceeds 20 Å in the accepted production map;
- Ag/W mass ratio leaves 0.999–1.001;
- either phase loses the intended connected skeleton;
- MD loses atoms or energy runs away;
- cold-state potential energy continues a strong monotonic relaxation;
- source/no-source accounting does not close.

Report separately:
- mapping PASS/FAIL;
- geometry PASS/FAIL;
- cold-baseline PASS/FAIL;
- TTM source/no-source PASS/FAIL.

Do not collapse these into one overall PASS if one stage is still blocked.
