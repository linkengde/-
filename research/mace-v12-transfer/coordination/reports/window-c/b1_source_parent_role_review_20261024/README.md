# B1 source and parent-role review — 2026-10-24

## Result

The source and lineage review is complete as a proposal review. It confirms 24 current geometry-only files against their package manifests: six COD-linked fcc Ag structures, six Bai et al. Ti₃SiC₂ structures, and twelve Ag liquid geometries rescaled to the Abdullaev et al. (2020) liquid-density correlation. File hashes, formulas, atom counts, label-free status, disabled DFT launch flags, and liquid source-parent hashes pass the read-only verifier.

Eighteen geometries have an inspected source record usable for geometry proposals (six solid Ag and twelve liquid Ag). The six Bai et al. Ti₃SiC₂ geometries remain conditional because the paper's numerical structure table has not yet been checked against the reconstructed CIF. This review does not freeze training roles, assign dev/test data, authorize DFT, or integrate candidates into a dataset.

## Candidate source and grouping decisions

| Candidate family | Count | Evidence and limitation | Parent/split rule |
|---|---:|---|---|
| Solid Ag, COD 9013045 | 6 | COD records fcc Ag at 293 K with `a = 4.086 Å` and links the record to Suh, Ohta & Waseda (1988). The COD bytes and candidate manifests are audited; the cited Springer article's full text remains unread. These are one crystallographic source family, not independent polymorphs. | Keep both supercell sizes and all displacement variants in one lineage group. Candidate geometry only. |
| Liquid Ag, Abdullaev et al. (2020) density reference | 12 | The reference density is from the retrieved CC BY 3.0 primary paper at 1250 and 1600 K. The coordinates themselves are rescaled EAM-generated liquid snapshots, so the experimental source supports the density target, not independent atomic configurations. Assael et al.'s 2012 fit is a secondary cross-check only; its 1250 K estimate remains unverified against the primary paper and does not replace the Abdullaev anchor. | Six variants from each EAM seed stay together across both temperatures and all three density ratios. Do not treat the two seeds as independent validation against the shared EAM source. |
| Ti₃SiC₂, Bai et al. (2010) | 6 | Separate computational-paper lineage from COD 9009647, but the numerical table has not been checked against the reconstructed candidate CIF. Hold source acceptance pending that comparison. | Keep both supercells and their parent/displacement variants in one source group. If COD-derived Ti₃SiC₂ candidates remain, group them with this phase for leakage protection; neither source family is a dev/test holdout. |

The older 30-row `B1_ROLE_FREEZE_DRAFT.csv` is not the current source-reviewed shortlist. Its six solid-Ag rows use the EAM lattice reference and are superseded by the COD-linked geometries for source-based proposals. Its six Ti₃SiC₂ rows derive from COD 9009647, the same structural source family as existing controls, so they do not supply independent-source validation. Its 18 liquid-Ag rows use the model-derived EAM density anchor: the six 1050 K entries fall below the 1234.93 K melting point and outside the Abdullaev liquid-fit range, while the twelve 1250/1600 K entries are superseded as reference-scaled candidates by the separate twelve-copy package. The legacy files are retained unchanged for provenance.

## Role and calculation gates

- Every file reviewed here remains a geometry proposal with no DFT energy or force labels. No candidate has a train, development, or test assignment.
- The 24 candidates are not 24 independent sources. Parent groups, source lineages, and close structural relationships must remain intact in any later split. The Ti₃SiC₂ COD and Bai source families share the same α-MAX phase and must not be placed on opposite sides of a validation split.
- The V20 Window A card snapshot is dated `2026-10-24T03:28:13Z`: it records the `k6×6×1, σ=0.05 eV` slab control running at SCF iteration 50, with `σ=0.20 eV` queued after its archive verification. The card is not live telemetry; execution and task-card updates remain with its registered owner.
- The project-wide DFT method and energy target are still unfrozen. The V19 `free_energy_eV` sidecar resolves availability for the 36 training candidates only; it does not select a project-wide target or authorize training.
- The next admissible action is source-table verification for Bai et al. (2010), followed by owner review of this candidate grouping. Only after the V20 controls and the method/energy policy gates are complete should an owner-bound pilot be registered. No DFT or training was started in this review.

## Reproduction and evidence

Run the read-only audit from the repository root:

```bash
/workspace/.venvs/mace/bin/python research/mace-v12-transfer/coordination/reports/window-c/b1_source_parent_role_review_20261024/verify_role_review.py
```

The result is stored in [role_review.json](role_review.json). It verifies candidate file hashes and content against the existing manifests and records the Window A card snapshot used here. The script does not alter source candidates or task cards.

Source packages:

- [Solid Ag, COD 9013045](../ag_solid_cod9013045_reference_candidates_20261024/README.md)
- [Liquid Ag density-scaled geometry copies, Abdullaev et al. (2020)](../abdullaev2020_liquid_ag_density_rescaled_20261024/README.md)
- [Ti₃SiC₂ candidate from Bai et al. (2010)](../independent_ti3sic2_bai2010_candidate_20261024/README.md)
- [Liquid Ag source audit and Assael secondary comparison](../liquid_ag_density_source_audit_20261024/README.md)
