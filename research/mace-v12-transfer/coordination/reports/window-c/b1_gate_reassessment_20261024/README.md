# B1 source and execution gate reassessment — 2026-10-24

This read-only reassessment consolidates the new source findings and keeps the proposed B1 work behind the existing ownership and DFT authorization gates.

| Gate | Current evidence | Disposition |
|---|---|---|
| Liquid Ag density | Abdullaev et al. (2020), Eq. (2), 1234.93–1631 K; full CC BY 3.0 paper retrieved. The fit gives 9.312096308 g/cm³ at 1250 K and 8.963636308 g/cm³ at 1600 K. | Resolved for geometry proposals at 1250/1600 K. The 1050 K liquid proposals remain unsupported by this liquid fit and are excluded from the new referenced set. |
| Solid Ag structure source | COD 9013045 records fcc Ag at 293 K, `a = 4.086 Å`, Fm-3m, and cites Suh, Ohta & Waseda (1988), DOI `10.1007/BF01174717`. The supplied full paper confirms the same 293 K parameter and fcc phase. | The COD record and paper are one experimental lineage. Six geometry-only proposals are prepared; keep the family grouped. |
| Ti₃SiC₂ independent source | Bai et al. (2010) is a separate computational paper from COD 9009647. The supplied accepted manuscript's Table 1 confirms the candidate CIF parameters `a = 3.0592 Å`, `c = 17.636 Å`, `zTi2 = 0.1351`, and `zC = 0.5726`. | Numerical provenance is verified. Six derived geometries remain one α-phase source group and are not independent validation structures. |
| Geometry roles and splits | The 24-file shortlist includes six COD/Suh solid Ag geometries, six Bai Ti₃SiC₂ geometries, and twelve liquid Ag geometries whose density targets use Abdullaev 2020 while their atomic configurations remain EAM-generated. | Source/reference checks are complete, but no train/dev/test roles are frozen. Keep each parent lineage grouped; owner review and explicit batch registration remain required. |
| Free-energy key availability | A separate [V19 train36 sidecar](../v19_free_energy_completion_20261024/README.md) now has `free_energy_eV` for all 36 training-candidate frames: 30 recovered from the audited V18 archive sidecar and six reconciled to V19 owner-verified output/summary/log archives. | This resolves key availability for the V19 train candidate only. It does not add free energies to development frames, select the training target, alter the registered dataset, or authorize training. |
| DFT method and target | The earlier plan proposes GPAW 26.7.0, ASE 3.29.0, PBE/PAW 1.2.1, 500 eV, and `sigma=0.10 eV`; it calls for `free_energy` as the force-consistent target while retaining native energy separately. | Proposal only. One energy key for the complete intended training arm, converged phase-specific k-meshes, and the treatment of existing native-energy data are not frozen. Do not combine native and free energies in one target. |
| V20 controls | The latest committed Window A task snapshot is timestamped 2026-10-24 03:28:13 UTC: k6×6×1, σ=0.05 was at SCF iteration 50 and σ=0.20 was queued. Owner remains the registered Window A instance. | Snapshot is not live telemetry. Leave execution and reporting with its owner; do not duplicate or take over. |

## Geometry packages

- [Solid Ag, COD 9013045](../ag_solid_cod9013045_reference_candidates_20261024/README.md): six source-linked fcc geometries at the measured 293 K lattice parameter, in two supercell families.
- [Liquid Ag, Abdullaev 2020](../abdullaev2020_liquid_ag_density_rescaled_20261024/README.md): twelve isotropically rescaled geometry copies at 1250 and 1600 K. The earlier EAM-derived files are retained unchanged.
- [Ti₃SiC₂, Bai 2010](../independent_ti3sic2_bai2010_candidate_20261024/README.md): reconstructed independent-computation candidate and six geometry-only supercells; its Table 1 source parameters are verified.
- [Supplied-paper audit](../provided_primary_papers_audit_20261024/README.md): Bai 2010 and Suh 1988 full-text checks, including the duplicate-scan hash.

No DFT, MACE training, MD, or TTM was started in this reassessment. The source-paper checks are now complete, superseding the earlier access limitations in this report's first version. The next gates are owner review and freeze of candidate roles, selection of one complete energy convention, resolution of phase-specific DFT settings using the remaining V20 controls, and owner-bound pilot registration. The V19 train36 free-energy sidecar does not resolve missing free energies on all frames or validation data; do not mix native and free energies. The 30-label batch in the earlier plan remains a proposal, not an approved job list.

Cloud-environment revision 24 is published. The former publisher-route tunnel failures did not establish the papers' access status; the supplied full texts now resolve the relevant source checks without inferring anything about paywalls or account access.
