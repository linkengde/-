# B1 source and execution gate reassessment — 2026-10-24

This read-only reassessment consolidates the new source findings and keeps the proposed B1 work behind the existing ownership and DFT authorization gates.

| Gate | Current evidence | Disposition |
|---|---|---|
| Liquid Ag density | Abdullaev et al. (2020), Eq. (2), 1234.93–1631 K; full CC BY 3.0 paper retrieved. The fit gives 9.312096308 g/cm³ at 1250 K and 8.963636308 g/cm³ at 1600 K. | Resolved for geometry proposals at 1250/1600 K. The 1050 K liquid proposals remain unsupported by this liquid fit and are excluded from the new referenced set. |
| Solid Ag structure source | COD 9013045 records fcc Ag at 293 K, `a = 4.086 Å`, Fm-3m, and cites Suh, Ohta & Waseda (1988), DOI `10.1007/BF01174717`. The COD CIF bytes and SHA256 are recorded in the separate geometry package. | A source-linked parent family is prepared. The official Springer PDF request is currently rejected by the environment egress tunnel (HTTP 403), before a publisher response; access status is therefore unknown. B source/role review remains required. |
| Ti₃SiC₂ independent source | Bai et al. (2010) is a different computational paper from COD 9009647. The candidate CIF was reconstructed from the user-supplied audit; article metadata was checked, but the exact parameter table was not independently inspected. OpenAlex marks the paper closed, and the official ScienceDirect PDF request is currently rejected by the environment egress tunnel (HTTP 403), before a publisher response. | Computational-source candidate only; numerical provenance still needs full-text verification or a different accessible source. The existing COD-derived family remains same-source-family material. |
| Geometry roles and splits | Ag/COD 9013045, Ti₃SiC₂/Bai 2010, and liquid-Ag reference-linked geometries are all geometry-only. Candidate supercells and seeded variants must remain grouped by parent family. | No train/dev/test roles are frozen. Do not integrate any candidate into a model dataset before source/role review and explicit batch registration. |
| DFT method and target | The earlier plan proposes GPAW 26.7.0, ASE 3.29.0, PBE/PAW 1.2.1, 500 eV, and `sigma=0.10 eV`; it calls for `free_energy` as the force-consistent target while retaining native energy separately. | Proposal only. One energy key for the complete intended training arm, converged phase-specific k-meshes, and the treatment of existing native-energy data are not frozen. Do not combine native and free energies in one target. |
| V20 controls | The latest committed Window A task snapshot is timestamped 2026-10-24 03:28:13 UTC: k6×6×1, σ=0.05 was at SCF iteration 50 and σ=0.20 was queued. Owner remains the registered Window A instance. | Snapshot is not live telemetry. Leave execution and reporting with its owner; do not duplicate or take over. |

## Geometry packages

- [Solid Ag, COD 9013045](../ag_solid_cod9013045_reference_candidates_20261024/README.md): six source-linked fcc geometries at the measured 293 K lattice parameter, in two supercell families.
- [Liquid Ag, Abdullaev 2020](../abdullaev2020_liquid_ag_density_rescaled_20261024/README.md): twelve isotropically rescaled geometry copies at 1250 and 1600 K. The earlier EAM-derived files are retained unchanged.
- [Ti₃SiC₂, Bai 2010](../independent_ti3sic2_bai2010_candidate_20261024/README.md): reconstructed independent-computation candidate and six geometry-only supercells; source-table verification remains open.

No DFT, MACE training, MD, or TTM was started in this reassessment. The next useful gate is source/role review of the Ag and liquid geometry packages, plus legal verification of Bai et al.'s numerical table. After that, the data owner must freeze roles, a single energy target and phase-specific numerical settings, and obtain an owner-bound pilot registration. The 30-label batch in the earlier plan remains a proposal, not an approved job list.

The cloud-environment draft now includes `link.springer.com` and `www.sciencedirect.com` (revision 21). Those official PDF routes can be retried after the user publishes that draft. No paywall or account status is inferred from the present tunnel rejection.
