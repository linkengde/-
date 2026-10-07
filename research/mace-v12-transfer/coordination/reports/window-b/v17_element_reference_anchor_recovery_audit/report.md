# Four unary reference anchor recovery audit

## Evidence

Ag4, Ti2, Si8 and C8 (V16 rows10–13) retain exact ordered geometry, cell/PBC, REF_energy and REF_forces through the previously documented V11–V14 lineage. Four-row ledger records original import commits, source hashes, declarations and exact missing-artifact claims. No persistent historical atom IDs exist; none were invented. Original calculator inputs, SCF logs, getter convention, GPAW/PAW version and elemental k-point meshes remain unverified. PW-PBE500/Fermi0.1 is a stored declaration, not an original-run certificate. Prior source tracing found original records absent; this focused path-history search found no independently certified exact-geometry replacement. It does not claim exhaustive Git-blob equivalence search.

## Rank and energy gauge

Rank sensitivity for all43 proposed rows is provided in CSV/JSON, including removal of each unary row and all four. For additive per-species energy offsets c, energies change by N c while coordinate forces are unaffected. Any null vector v of N leaves all composition energy equations unchanged under c -> c + t v. A rank-deficient composition matrix cannot identify all four independent offsets. Full rank alone does not certify reference method consistency. Unary solid energies are not isolated-atom chemical potentials; energy/atom here is a label quotient, not a transferability guarantee.

## Minimum safe recovery plan

First request the exact original elemental run records identified in the inherited ledger: fixed-cell input/script, k-point mesh, smearing distribution, XC/cutoff, GPAW/PAW versions, SCF stopping evidence, native/free energies and forces, with input/output hashes. Until obtained, retain explicit provisional provenance flags. If unavailable, A may approve four exact-geometry fixed-cell relabelings; preserve atom order, counts, cell and PBC, generate no synthetic historical IDs, record native extrapolated and free energies separately, and hash all inputs/outputs/logs. PW500/Gamma/Fermi0.1 is a proposed matched target only, not recovered historical settings; small bulk elemental cells require A to review k-point adequacy before adopting Gamma. Do not substitute another lattice, relaxed cell or method archive merely because the element name matches. Four relabelings close these four gaps; they do not certify the ten other partially supported interface frames. One unary row may suffice to resolve a rank-one gauge deficiency, but cannot resolve all four provenance uncertainties.

No DFT, inference, training, dataset or role edits were performed. No development-validation or withheld output content was accessed. Reproduce with audit.py in the existing GPAW Python environment. The history search is deliberately path-only and the full previous source trace is reused rather than repeated.

Rank results: [{"scenario": "full", "rank": 4}, {"scenario": "remove_Ag4", "rank": 4}, {"scenario": "remove_Ti2", "rank": 4}, {"scenario": "remove_Si8", "rank": 4}, {"scenario": "remove_C8", "rank": 4}, {"scenario": "remove_all_four", "rank": 4}]
