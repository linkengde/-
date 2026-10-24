# B1 candidate-role and energy-key decision — 2026-10-24

## Decision

The 24 source-reviewed geometries are frozen into a **B1 training-candidate pool only**. None is a development or test structure. This freezes the candidate roster and lineage policy; it does not add labels, alter an existing training file, or authorize DFT, model training, or dataset integration.

The energy key for a future finite-smearing B1 thermal-label arm is frozen as `free_energy_eV` (cell total, eV), paired with the corresponding DFT forces. Preserve GPAW's native extrapolated energy in a separate field for provenance and comparison; do not train on a mixture of native energy and free energy. This is the B1 thermal-arm policy, not a retroactive change to V18/V19 or a project-wide conversion of legacy `REF_energy` data. It follows the existing [method and role rules](../unified_potential_dft_data_plan_20261010/METHOD_AND_ROLE_RULES.md): forces at finite Fermi smearing are consistent with the free-energy functional.

## Frozen candidate groups

| Group | Count | Role and lineage rule |
|---|---:|---|
| Solid Ag, COD 9013045 / Suh 1988 | 6 | Training candidates only. All supercells and displacement variants stay together. COD and Suh refer to the same experiment. |
| Ti₃SiC₂ α phase, Bai 2010 | 6 | Training candidates only. All supercells and displacement variants stay together; keep this source family grouped with any COD 9009647 Ti₃SiC₂ controls for leakage protection. |
| Liquid Ag, EAM seed 501 | 6 | Training candidates only. Keep all temperatures and density-ratio variants together. Abdullaev 2020 supports the density targets; coordinates are EAM-generated. |
| Liquid Ag, EAM seed 902 | 6 | Training candidates only. Keep all temperatures and density-ratio variants together. This shares the same EAM generation lineage as seed 501 and cannot support independent-source validation. |

The 0.97, 1.00, and 1.03 density ratios are deliberate geometry-sampling points around the Abdullaev reference, not experimental uncertainty bounds. No candidate can be reassigned to development or test. If a candidate fails DFT convergence, geometry checks, or duplicate checks, exclude it from the eventual labeled pool with a recorded reason rather than move it into validation.

The [role ledger](role_ledger.csv) enumerates every candidate path and SHA256. Its builder verifies the source manifests, current candidate bytes, label-free status, DFT-disabled flags, counts, and split restrictions. The ledger is an eligibility decision for a future training batch only; all current files remain geometry-only proposals.

## Energy and method gates

- B1 future thermal labels: use `free_energy_eV` for the energy loss and the forces computed with the same finite-smearing calculation. Retain native extrapolated energy separately; never substitute it under the free-energy key.
- Existing V19 has a complete 36-frame free-energy candidate sidecar, and its three validation frames have free energy recorded. The sidecar is not the registered training dataset. Before reuse, resolve its train/dev method differences and the V17 rank-recovery k-mesh caveat.
- Proposed common GPAW/PBE settings remain GPAW 26.7.0, ASE 3.29.0, PAW 1.2.1, 500 eV, and σ=0.10 eV. They are not yet frozen for B1. The registered Window A owner must finish the σ=0.05 and σ=0.20 V20 controls, and phase-specific k-point convergence is still needed for solid Ag, Ti₃SiC₂, and representative liquid Ag boxes.
- The latest task-card snapshot is timestamped 2026-10-24 03:28:13 UTC: σ=0.05 was at SCF step 50 and σ=0.20 was queued. This checkout has no live owner telemetry; the registered owner retains these jobs, and this decision did not duplicate them.
- Until those method gates pass, do not generate B1 labels or start model training. The current V20 task remains with its registered owner; this decision does not change that task or launch work in parallel.

## Verification result

Expected result: 24 unique geometry files; 6 Ag, 6 Ti₃SiC₂, and 12 liquid Ag; all have no DFT labels, all disable DFT launch, all have the `TRAINING_CANDIDATE_ONLY` role, and all explicitly have no development/test assignment. Liquid parent seeds 501 and 902 each contribute six files. Candidate hashes are checked against the existing manifests. See [ledger audit](ledger_audit.json).
