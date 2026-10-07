# V17 geometry-only split leakage audit

## Scope and result

35 V16 training geometries, eight V17 diagnostic inputs, three development inputs and six withheld inputs: 52 structures, 405 cross-role pairs. Exact file matches: 0; normalized exact geometry matches: 0; near matches: 0; near matches involving new V17 training: 0.

## Algorithms fixed before interpretation

Cell tolerance 1e-5 A, exact maximum displacement 1e-5 A, near Ag and framework RMS each <=0.15 A. Reused split-design identity mapping is supplemented by species-preserving Hungarian assignment on squared periodic MIC distances, framework-median translation alignment and up to five iterations. No rotation or cell scaling. The assignment is a local geometric screen, not proof against every global symmetry. ID-matched displacement vectors are reported independently. Shared IDs show common lineage, not automatically leakage. File hash equality across multi-frame training files and single-frame candidates is a weak check; normalized coordinate checks are the substantive comparison.

## Interpretation and recommendation

{"structures": 52, "cross_role_pairs": 405, "exact_file_hits": 0, "exact_geometry_hits": 0, "near_geometry_hits": 0, "new_V17_training_near_hits": 0}

No exact/near cross-role hit was found under both stated algorithms. Newly added paired-training inputs do not introduce a hit under these thresholds. A may retain proposed roles subject to provenance/label review.

All candidates inherit existing small-cluster parent motifs. Absence of duplicate geometries establishes only local separation, not new morphology or thermal coverage. Same-ID maps in audit.json quantify correlated families separately from leakage. Species minima are compared against the lowest observed historical V16 distance for that element pair minus0.05 A; these flags are comparative triage, not chemical validity or a universal 1.75 A cutoff. Full species minima and flags are provided for every geometry.

Only geometry columns were parsed. Candidate hashes match the frozen role manifest. No development/withheld calculation output, label, summary, log or model score was opened; no calculation/inference/training was started. No datasets or role assignments were changed.

## Reproduce

Run audit.py with the existing GPAW Python environment; it imports the committed geometry parser without running its generator. SHA256SUMS.txt covers the report artifacts.
