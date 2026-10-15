# DFT label quality audit: Ag–Ti–Si–C potential qualification

## Scope and provenance

This is a read-only audit of archived data at these fixed sources:

- E DFT handoff: `4ee3daab860a846fec4ef0c0ca4b031a36cfc87e` (archive SHA-256 `4e35ce85e6a47295e07c6e4e688925bc343fee7cdd966d39e3619b9d00d7026f`).
- E unified static comparison: `dbcf2f7d6bf966b8fc0476224d2cba52c3866de0`.
- E SevenNet/MatterSim comparison: `0a4f596a82dc2eaaa1e2b365675507187a19a22b`.
- C 77-atom three-gap geometry archive: `3c8ffa28cf6f5aaa72981d8b9f1a29389f9d8d50` (ZIP SHA-256 `f37d5c84b919ac583c9811b557c8eacdb06254d2d4095b6634c0f9ab927ba678`).
- B historical V9–V11 and V12–V19 materials: `499ddf0870235594444abad780582ea4fcc89d61`.
- E2 48-atom and 114-atom audit source: E2 inventory and its recorded per-file hashes in the E handoff/report chain.

Archive/checksum manifests for the two downloaded handoff packages were verified before inspection. Source records in `DFT_MASTER_INVENTORY.csv` include the source commit, path, raw file SHA-256, frame number, available settings, role history and label status. Original data were not changed.

## Inventory method

ASE 3.29 parsed 247 `.xyz`/`.extxyz` source files into 554 frames without parse errors. Frames were deduplicated by an order-sensitive fingerprint of species, coordinates (rounded to 1e-6 Å), cell and PBC. This produced 120 unique ordered geometries. Exact geometry duplicates copied into several MACE dataset versions are represented once; different k meshes for the same geometry are not counted as new structures. The inventory preserves coordinates, species/atom IDs when present, cell, PBC, formula, raw source hashes and closest-pair distances.

Each geometry is classified conservatively:

- `VALID_DFT_LABEL`: a source establishes converged SCF plus energy and per-atom forces.
- `INCOMPLETE_LABEL`: DFT labels or a DFT source are indicated, but a required label or reliable convergence provenance is missing.
- `FAILED_DFT`: a failed calculation with no valid final label.
- `NON_DFT_DATA`: geometry/proposal/prediction without a DFT label.

## Result

| Geometry status | Unique geometries |
|---|---:|
| `VALID_DFT_LABEL` | 71 |
| `INCOMPLETE_LABEL` | 21 |
| `FAILED_DFT` | 0 standalone geometries |
| `NON_DFT_DATA` | 28 |
| Total | 120 |

The 71 confirmed geometries contain 3,588 distinct atom-force vectors: Ag 1,027; Ti 1,190; Si 472; C 899. They have 71 total-energy records, 50 free-energy records, and 5 explicit stress/virial records. A known 135-atom 2×2×1 attempt stopped at the memory guard before SCF (`scf_started=false`); it has no usable labels and is recorded as a failed run variant of an existing geometry rather than an independent structure.

The E handoff contains converged GPAW/ASE/PBE/PAW records with fixed geometry, raw results and explicit hashes. The newer V19 source metadata commonly identifies GPAW 26.7.0, ASE 3.29.0, PBE, 500 eV, Fermi–Dirac 0.10 eV, PAW dataset 1.2.1 and tight SCF thresholds. Across the full archive, k meshes, smearing (including 0.05, 0.10 and 0.20 eV records), PAW/setup metadata and result conventions vary. Several legacy copies do not carry enough run metadata to independently prove SCF status or the energy convention. The 21 incomplete cases therefore remain excluded from a strict training set until the raw run evidence is recovered.

Total energy and force-consistent free energy are separate quantities and are not interchangeable. Absolute energies from different k meshes, smearing, PAW setup or pseudopotential conventions are not merged without normalization. Most force arrays are in eV/Å when their source metadata states that unit; missing metadata is retained as unverified rather than inferred. Only five confirmed geometries have explicit stress/virial records, and their units/sign convention is source-specific. No zero stress labels were fabricated.

The component DFT records for the 135-atom interface (full, Ag-only, and Ti₃SiC₂-only) are tied to the same frozen optimized parent geometry where their coordinates map. They are decomposition references, not independent training structures. The 77-atom 2.16/2.40/2.64 Å series has geometry only and no DFT energies or forces.

## Finding

`DFT_LABEL_INTEGRITY = INSUFFICIENT_DATA`. The core converged energy/force subset is traceable and reusable after method/role harmonization. The entire archive is not yet a uniform DFT training set because 21 geometries lack complete provenance/labels, free-energy coverage is partial, stress/virial coverage is sparse, and method settings vary. The master CSV is the source-level record for every included geometry and its raw hash.
