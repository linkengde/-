# Independent-source Ti₃SiC₂ geometry candidate

The supplied CIF is a reconstructed derivative based on Bai et al.'s 2010 DFT study. Crossref and OpenAlex confirm the article's title, authorship, journal, and year. The article is closed access, so I could not independently inspect the paper's numerical tables; the parameter-to-paper link remains dependent on the supplied audit's transcription.

The CIF's SHA256 exactly matches the value in the supplied audit report: `cb361a88ab35948e9fc8e0e8f86b0ec7030ae0934c3bbf7c33ef68b674d3d3a5`. It is a generated derivative, not an original CIF downloaded from the paper or a crystallographic database. `source_audit_user_provided.md` preserves the source audit that accompanied it.

Parsing the supplied CIF with ASE expands it to 12 atoms (Ti₆Si₂C₄, two Ti₃SiC₂ formula units) in space group 194. The parsed volume and site-resolved nearest-neighbor distances agree with the supplied audit to its stated precision, and no atom pairs are closer than 0.5 Å.

This qualifies as a candidate with computational provenance independent of COD 9009647. It remains the same α-Ti₃SiC₂ structure type, so it is not a separate polymorph or independent experimental record. It has no training or validation role assigned.

Run `audit_candidate.py` with the workspace MACE environment's ASE installation to recheck the supplied CIF and regenerate `geometry_audit.json` without rewriting the source file.

## Supercell geometry proposals

`generate_supercell_candidates.py` creates six label-free geometries from this CIF: one parent and two seeded 0.025 Å displacement variants for each of the 72-atom 3×2×1 and 96-atom 2×2×2 supercells. The independent read-back audit in `audit_supercell_candidates.py` passes all six serialized files, including formula, atom count, finite coordinates, periodicity, source hash, and minimum-distance checks. These remain geometry-only proposals, with DFT disabled and roles unfrozen; the parent CIF's exact values still need full-text verification and B review. Reproduce them with `/workspace/.venvs/mace/bin/python generate_supercell_candidates.py`, then verify them with `/workspace/.venvs/mace/bin/python audit_supercell_candidates.py`.

One fixed-geometry, 12-atom GPAW PBE single point is recorded under `exploratory_dft/run_01/`. It converged in 32 SCF iterations using 500 eV, σ=0.10 eV, 12×12×2 k points, and four MPI ranks. Native energy and free energy are retained separately. The maximum force is 0.1351 eV/Å, so this is not a relaxed structure; one point also does not establish k-point or cutoff convergence. The calculation is diagnostic only and is excluded from training, development, and test data.
