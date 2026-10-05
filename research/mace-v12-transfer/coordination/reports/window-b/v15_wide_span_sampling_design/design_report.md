# Window B wider-span Ag-X sampling design

Analyzed main commit: `e830a2fcb922afe4cfcbf8b21bdbae4543e2287f`. This is a geometry-only proposal for future training acquisition. No DFT, MACE inference/training, MD or TTM was run. V15 holdout labels, summaries, DFT outputs and active files were not read. Candidates have no energy/force labels and require A review and DFT labels before entering a dataset.

## Evidence and separate sampling axes

V14 force RMSE failed on AgC/AgSi/AgTi (0.154/0.118/0.147 eV/Å); AgSi/AgTi also missed separating-force screening (0.235/0.193 eV/Å absolute error). Nearby scalar distances already occur in training, so a distance-only fix is insufficient. AgC/AgTi error concentrates near the contact shell; AgSi has local/outside RMSE 0.127/0.108 and needs broader surroundings. The four undocumented force-only LCAO references are not selected as parents.

Each parent is an energy+force PW-PBE frame in V14 train. The distance axis is the marked X→Ag vector, a local pair-normal direction; these small clusters do not define a unique global planar-interface normal. All Ag atoms move rigidly by ±0.20/±0.35 Å on this axis. Registry uses the normalized cross product of that vector with z, with ±0.45/±0.90 Å translations and no added normal displacement. Its radial pair distance necessarily increases as sqrt(anchor²+registry²); that increase is reported, not hidden by a compensating normal edit. Local-shell ±0.12 Å radial framework displacements keep Ag and both marked atoms fixed, and form a separate axis.

AgC/AgTi surroundings perturb two high-error framework IDs; AgSi perturbs six, selecting outside-shell atoms where the diagnostic and IDs allow. Directions follow each neighbor’s position from the marked X atom, not an inferred force correction. Unmoved atoms and the full cell/PBC/ID order are preserved.

| Interface | parent frame | anchor (Å) | initial normal bracket (Å) | PW train distance range (Å) | perturbed framework IDs |
|---|---:|---:|---|---|---|
| AgC | 29 | 2.2600 | 1.9100–2.6100 | 2.0000–2.4000 | [604590, 604806] |
| AgSi | 30 | 2.4401 | 2.0901–2.7901 | 2.2000–2.8000 | [366931, 367113, 367117, 368832, 368823, 367115] |
| AgTi | 28 | 2.5880 | 2.2380–2.9380 | 2.5880–2.6000 | [480805, 480812] |

AgSi refinement: retain all flagged coarse registry endpoints, and additionally propose -0.30/+0.20 Å lateral shifts after their shortest Ag contacts pass the same declared screen. Add a +0.50 Å distance point (2.9401 Å) because the original +0.35 endpoint remains below its existing 2.80 Å PW training distance. The 1.85/1.95/2.15 Å and parent-relative floors limit compression; a scalar endpoint outside the old range is not forced if it fails. This separates wider exploration from explicit interface-specific narrowing.

AgSi surroundings include three local and three outside-shell framework atoms (parent distances 4.584–5.424 Å for the outside subset); the exact IDs, distances and perturbation vectors are in JSON. Parent archive outputs match each selected V14 training frame exactly in IDs/elements/cell/PBC/coordinates, and their published verification is PASS.

## Geometry and duplication audit

Screened 55 committed geometry files / 166 frames: existing train/valid/test splits, historical DFT inputs/results and acquisitions, residual V15 inputs/results, prior B geometry proposals, and only the three V15 frozen holdout **inputs**. Inventory hashes and all comparisons are in JSON. All reference frames having the same composition as a proposal have compatible persistent IDs and were compared. Missing-ID frames have different compositions; their counts are included. No geometric equivalence under rotations/relabeling is claimed. Candidate-to-candidate matches are also recorded.

Ag-C/Ag-Si/Ag-Ti/Ag-Ag review floors are 1.85/1.95/2.15/2.45 Å, strengthened by 90% of the parent species-specific shortest Ag contact. These are conservative triage choices, not universal physical limits. Framework minima may not shorten below 95% of the parent minima; existing framework neighbors ≤3.2 Å are flagged outside 0.92–1.10 of their parent bond length. The source-inherited AgC C–Ti 1.7413 Å adjacency remains a framework bond and is not judged by an Ag-contact cutoff. Geometry flags require narrowing/drop/review, not forced labeling.

Exact duplicate: same IDs/elements, cell (1e-8 Å), PBC and framework-aligned coordinate maximum ≤1e-6 Å. Near duplicate: Ag RMS ≤0.15 Å and framework RMS ≤0.15 Å, with framework maximum ≤0.25 Å. This prevents a moving small Ag sublattice being diluted by whole-system RMS. Other-PBC coordinate matches are flagged separately. Local-shell variants intentionally remain near their training parent and require review as correlated acquisition, not independent validation. Any candidate near a frozen V15 holdout is withheld.

| Candidate | axis offset (Å) | pair target / actual (Å) | Ag-contact min (Å) | exact / near existing | decision |
|---|---:|---:|---:|---|---|
| `AgC_wide_distance_m0p35` | -0.35 | 1.9100 / 1.9100 | 1.9100 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgC_wide_distance_m0p20` | -0.20 | 2.0600 / 2.0600 | 2.0600 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgC_wide_distance_p0p20` | +0.20 | 2.4600 / 2.4600 | 2.2468 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgC_wide_distance_p0p35` | +0.35 | 2.6100 / 2.6100 | 2.1882 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgC_wide_registry_m0p90` | -0.90 | 2.4326 / 2.4326 | 1.5954 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgC_wide_registry_m0p45` | -0.45 | 2.3043 / 2.3043 | 2.0442 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgC_wide_registry_p0p45` | +0.45 | 2.3043 / 2.3043 | 2.3043 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgC_wide_registry_p0p90` | +0.90 | 2.4326 / 2.4326 | 2.4326 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgC_wide_surroundings_m0p12` | -0.12 | 2.2600 / 2.2600 | 2.2600 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |
| `AgC_wide_surroundings_p0p12` | +0.12 | 2.2600 / 2.2600 | 2.2600 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |
| `AgSi_wide_distance_m0p35` | -0.35 | 2.0901 / 2.0901 | 2.0901 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgSi_wide_distance_m0p20` | -0.20 | 2.2401 / 2.2401 | 2.2401 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_distance_p0p20` | +0.20 | 2.6401 / 2.6401 | 2.3819 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_distance_p0p35` | +0.35 | 2.7901 / 2.7901 | 2.3228 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_registry_m0p90` | -0.90 | 2.6008 / 2.6008 | 1.8281 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgSi_wide_registry_m0p45` | -0.45 | 2.4812 / 2.4812 | 2.1407 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgSi_wide_registry_p0p45` | +0.45 | 2.4812 / 2.4812 | 1.9532 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgSi_wide_registry_p0p90` | +0.90 | 2.6008 / 2.6008 | 1.5037 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgSi_wide_registry_m0p30` | -0.30 | 2.4585 / 2.4585 | 2.2553 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_registry_p0p20` | +0.20 | 2.4483 / 2.4483 | 2.2031 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_distance_p0p50` | +0.50 | 2.9401 / 2.9401 | 2.2721 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgSi_wide_surroundings_m0p12` | -0.12 | 2.4401 / 2.4401 | 2.4029 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |
| `AgSi_wide_surroundings_p0p12` | +0.12 | 2.4401 / 2.4401 | 2.4029 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |
| `AgTi_wide_distance_m0p35` | -0.35 | 2.2380 / 2.2380 | 1.9489 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgTi_wide_distance_m0p20` | -0.20 | 2.3880 / 2.3880 | 2.0764 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgTi_wide_distance_p0p20` | +0.20 | 2.7880 / 2.7880 | 2.4264 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgTi_wide_distance_p0p35` | +0.35 | 2.9380 / 2.9380 | 2.4850 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgTi_wide_registry_m0p90` | -0.90 | 2.7400 / 2.7400 | 1.5243 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgTi_wide_registry_m0p45` | -0.45 | 2.6268 / 2.6268 | 1.9355 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgTi_wide_registry_p0p45` | +0.45 | 2.6268 / 2.6268 | 2.3157 | 0 / 0 | GEOMETRY_PASS_A_REVIEW |
| `AgTi_wide_registry_p0p90` | +0.90 | 2.7400 / 2.7400 | 2.0139 | 0 / 0 | FLAG_GEOMETRY_NARROW_OR_DROP |
| `AgTi_wide_surroundings_m0p12` | -0.12 | 2.5880 / 2.5880 | 2.2080 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |
| `AgTi_wide_surroundings_p0p12` | +0.12 | 2.5880 / 2.5880 | 2.2988 | 0 / 5 | REVIEW_NEAR_EXISTING_TRAINING_VARIANT |

Candidate-to-candidate screening finds zero exact duplicates and 11 near pairs. These include adjacent coarse distance endpoints and opposite small-shell perturbations. Their IDs and RMS are in JSON. A should choose a subset from each near group or document the deliberate correlated sampling rationale; the per-candidate PASS above refers to geometry and distinctness from committed reference inputs, not to independence of every proposed pair.

## Interface-specific narrowing and next selection


**AgC:** 4 proposals pass geometry/distinctness triage: `AgC_wide_distance_m0p20`, `AgC_wide_distance_p0p20`, `AgC_wide_distance_p0p35`, `AgC_wide_registry_p0p45`. Flagged geometry endpoints: `AgC_wide_distance_m0p35` (short_Ag-C), `AgC_wide_registry_m0p90` (short_Ag-Si, short_Ag-Ti), `AgC_wide_registry_m0p45` (short_Ag-Ti), `AgC_wide_registry_p0p90` (short_Ag-Si). Start review with passing ±0.20 Å distance points and the +0.45 Å registry point; retain the flagged ±0.35/±0.90 endpoints only as the documented outer exploratory bracket. Narrow any flagged compression toward the parent in ≤0.05 Å increments during a later design, not by silently replacing these proposals. If an outer labeled point later gives extreme force or SCF problems, drop/narrow it; geometry alone cannot predict those outcomes.

**AgSi:** 6 proposals pass geometry/distinctness triage: `AgSi_wide_distance_m0p20`, `AgSi_wide_distance_p0p20`, `AgSi_wide_distance_p0p35`, `AgSi_wide_registry_m0p30`, `AgSi_wide_registry_p0p20`, `AgSi_wide_distance_p0p50`. Flagged geometry endpoints: `AgSi_wide_distance_m0p35` (short_Ag-Si, short_Ag-Ti), `AgSi_wide_registry_m0p90` (short_Ag-C, short_Ag-Si), `AgSi_wide_registry_m0p45` (short_Ag-C), `AgSi_wide_registry_p0p45` (short_Ag-C), `AgSi_wide_registry_p0p90` (short_Ag-C, short_Ag-Si, short_Ag-Ti). Start review with passing ±0.20 Å distance points and the narrowed AgSi -0.30/+0.20 Å registry points; its +0.50 Å distance point explicitly extends beyond the previous upper distance; retain the flagged ±0.35/±0.90 endpoints only as the documented outer exploratory bracket. Narrow any flagged compression toward the parent in ≤0.05 Å increments during a later design, not by silently replacing these proposals. If an outer labeled point later gives extreme force or SCF problems, drop/narrow it; geometry alone cannot predict those outcomes.

**AgTi:** 4 proposals pass geometry/distinctness triage: `AgTi_wide_distance_m0p20`, `AgTi_wide_distance_p0p20`, `AgTi_wide_distance_p0p35`, `AgTi_wide_registry_p0p45`. Flagged geometry endpoints: `AgTi_wide_distance_m0p35` (short_Ag-C), `AgTi_wide_registry_m0p90` (short_Ag-Si, short_Ag-Ti), `AgTi_wide_registry_m0p45` (short_Ag-Si, short_Ag-Ti), `AgTi_wide_registry_p0p90` (short_Ag-C). Start review with passing ±0.20 Å distance points and the +0.45 Å registry point; retain the flagged ±0.35/±0.90 endpoints only as the documented outer exploratory bracket. Narrow any flagged compression toward the parent in ≤0.05 Å increments during a later design, not by silently replacing these proposals. If an outer labeled point later gives extreme force or SCF problems, drop/narrow it; geometry alone cannot predict those outcomes.

A should prioritize distinct candidates that expand distance or registry coverage, then a small controlled shell subset informed by the local/outside error pattern. Preserve all frozen V15 holdouts and their selection independence. Do not label these proposals blind tests or call V15 passed/training-ready. Shared small-cluster parent lineage limits claims about morphology, extended interfaces, hot structures and liquid Ag.

The JSON contains each parent file/frame/geometry hash, candidate SHA256, atom IDs, target/actual pair distance, translation, displacement RMS/max, species minima with atom IDs, affected neighbors/bond changes, exact/near hits and decision. All 33 proposal files are retained, including flagged endpoints, so A can review the intended coarse span. Flagged structures are not approved for DFT.
