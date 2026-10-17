# Ti₃SiC₂ surface / DEV DFT input preparation for C review

## Scope and state

This package prepares the public E input candidates named in C's handoff commit `1b641b31e750731509a34f033eb5eeee47c029d7` (handoff file SHA-256 `7b2b2fe73306dbcb9ff764ea54a43b7a9780473de24c232952b6e31570c1c1da`). It contains four generated, unlabeled basal-surface structures: reference and one deterministic response for each termination. It does not contain DFT outputs and does not start GPAW, LAMMPS, MD, relaxation, fitting, or training.

The six previously approved bulk points remain exactly the six entries in C's separate data-use approval `52e722f30e70f4d0d7ec2a43643930d95a09dca3`. That approval does not grant a role to these new surface structures. Their TRAIN designation remains a proposal pending C review. Fitting remains unauthorized.

`PUBLIC_TASK_MANIFEST.json` records E tasks 1–4 as ready for input review and tasks 5–6 as HOLD. The only Ti₃SiC₂ CIF in the repository is COD 9009647; available alternatives either share its bulk lineage, derive from the 135-atom Ag interface whose 114-atom TSC-only descendant was already used in E5 response analysis, or have legacy MACE training history. None is silently relabelled as an independent DEV family. The sealed TEST inputs remain entirely outside this package.

## Geometry construction and source

The sole source is the tracked COD CIF `research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`, SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`. It cites Kisi et al., *Journal of Physics and Chemistry of Solids* 59 (1998), 1437–1443, DOI `10.1016/S0022-3697(98)00226-1`. The same file is the source for the archived 48-atom parent proposal (SHA-256 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`); its in-plane vectors match this construction after a 2×2 repeat.

`build_surface_candidates.py` constructs 2×2 in-plane, two-c-period slabs with one matching terminal plane at each end, then adds an 18 Å periodic-image vacuum gap. The cell is 6.115 × 6.115 × 53.247 Å (hexagonal in-plane angle 120°), PBC `[true,true,true]`, 100 atoms. Ti-terminated composition is Ti₅₂Si₁₆C₃₂; Si-terminated composition is Ti₄₈Si₂₀C₃₂. Each is a symmetric, non-stoichiometric termination slab; no passivation, atom deletion or random displacement was used. The structures are as-cleaved experimental-CIF geometries, not relaxed or claimed to be equilibrium surfaces. Surface energies are not inferred.

For each response structure, one matching top/bottom terminal-site pair is displaced symmetrically by 0.02 Å along the surface normal (bottom −z, top +z): Ti IDs 1 and 25; Si IDs 5 and 24. This is a local deterministic response within each slab family and preserves its verified mirror/inversion symmetry. Atom order and `source_atom_id` are unchanged from the corresponding reference geometry. The nearest pair distance is 2.0882 Å in both reference and perturbed inputs; no contact is below 1.5 Å.

The canonical coordinate SHA hashes element order, 8-decimal serialized Cartesian positions, cell, PBC and source atom IDs. The builder and validator both recompute this digest from the saved extxyz files.

## Candidate list

| Task | Structure | Atoms / composition | Input SHA-256 | Candidate role |
|---|---|---:|---|---|
| 1 | `TSC100_TI0001_COD9009647_2x2x2_REF` | 100 / Ti₅₂Si₁₆C₃₂ | `b1948a44b89efe17ee7c2891b984167982067feeffa774b7f6fa77d9e9158d86` | TRAIN candidate |
| 2 | `TSC100_TI0001_COD9009647_2x2x2_RESP_Z020` | 100 / Ti₅₂Si₁₆C₃₂ | `afa99ddae62e861bf641d3e7f2d908c78dfb0d466e762a25ebaeec64b0a91078` | TRAIN candidate; paired with task 1 |
| 3 | `TSC100_SI0001_COD9009647_2x2x2_REF` | 100 / Ti₄₈Si₂₀C₃₂ | `ac8b804516424215694c284d3bcc36d23da590f80390f47a0805c87cedb764c4` | TRAIN candidate |
| 4 | `TSC100_SI0001_COD9009647_2x2x2_RESP_Z020` | 100 / Ti₄₈Si₂₀C₃₂ | `f78560dfe891ea86d8a7da5843216050430f9aef8303c463231f31c19ca10da0` | TRAIN candidate; paired with task 3 |
| 5–6 | Independent DEV reference + response | — | — | HOLD: no verified independent non-TRAIN parent source |

The two surface pairs each form one parent family. Their four geometries are not four independent families. Same-parent descendants must remain together in TRAIN. No independent DEV or TEST claim is made.

## Method proposal

`PARAMETER_LOCK_PROPOSAL.json` proposes GPAW 26.7.0 / ASE 3.29.0 / gpaw-data 1.2.1; PBE/PAW; PW 500 eV; Fermi–Dirac 0.10 eV; non-spin-polarized; convergence thresholds `energy=density=eigenstates=1e-5`; maximum 160 SCF iterations; the successful P0 MSR1 settings; GPAW default Davidson; MPI 4 / OMP 1. Coordinates and cell stay fixed. No relaxation occurs.

The proposed surface mesh is standard Monkhorst–Pack `(6,6,1)` passed as a GPAW `kpts` tuple, with no explicit `gamma` flag or extra shift. `KPOINT_REVIEW.json` reproduces the archived parent `(6,6,2)` in-plane coordinates exactly and verifies the slab's single z-coordinate is Γ. The slab's in-plane cell vectors equal the parent vectors. This is a density-based mesh proposal, not a slab k-point convergence result; parent mesh convergence does not prove convergence for these slab geometries.

Each successful run should archive SCF evidence, native total energy, force-consistent free energy as a distinct field, all per-atom forces, full stress tensor, stdout/stderr/logs, runtime and hashes. The four surface candidates may contribute forces only after C approves their role. Energies and stresses are diagnostics; stress for a vacuum slab is normalized by the full 3D cell volume and is not a bulk elastic label.

## Resource and unresolved review items

The environment currently has 4 CPU cores, 32 GiB cgroup RAM, 16,659,927,040 free workspace bytes and a passing four-rank GPAW/ASE import test. An offline memory-scale estimate is approximately 6.18 GB wave-function storage plus 0.62 GB projectors per 100-atom candidate before FFT, mixer, Python and allocator overhead; no peak memory was measured. A coarse runtime scaling from the converged 48-atom 6×6×2 run is about 3.7 hours for a 41-step surface run. A four-hour cap would have little margin. One temporary checkpoint/scratch reservation is budgeted at 10 GiB; its actual size and resume mechanism are unverified. See `RESOURCE_BUDGET.md`.

Before execution, C must review the four structure hashes, surface stoichiometry/termination, response atom IDs and amplitude, k-point coordinates, method lock, resource/time budget, data roles, and a valid source for DEV tasks 5–6. A fresh C-reviewed execution manifest and runner hash must bind all inputs. Until then `execution_authorized=false`; there will be no DFT launch.
