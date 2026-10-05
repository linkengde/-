# Window B inherited-method provenance trace

Analyzed main: `3e8577e4f30dc476960f7d5cedd5c6bc7b3649e7`. Traced exactly the 15 `unknown_inherited_method` train records from the published isolated V15 builder manifest. No dataset/builder/model/main log was changed; no DFT, MACE inference/training, MD or TTM was run. V15 blind inputs and all blind label/output/summary files were not read. Original-machine absolute paths were not accessed.

**Result: 0 fully confirmed original calculations, 14 partial recoveries, 1 unresolved.** All 15 geometries and REF energy/force arrays match their earliest committed V11 handoff rows and subsequent V12/V13/V14 copies within 1e-8 Å/eV/eVÅ⁻¹ component tolerances. This establishes serialized label lineage, not original-run provenance.

The prior audit/draft used absence of `source_method` as an unknown-method flag. Fourteen rows actually preserve **`dft_method`**, an alternate per-frame declaration. Those declarations should be exposed in A’s provenance ledger as declarations rather than replaced with invented settings. This new report corrects that interpretation while leaving the prior artifacts unchanged.

## Per-frame findings

| Draft frame | config | metadata status | declared k points | original SCF/output linkage |
|---|---|---|---|---|
| 0 | `Ti_gap_2p3` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 1 | `Ti_gap_2p6` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 2 | `Ti_gap_2p9` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 3 | `Si_gap_2p6` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 4 | `C_gap_2p6` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 5 | `Si_gap_2p4` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 6 | `C_gap_2p4` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 7 | `Si_gap_2p8` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 8 | `C_gap_2p8` | partial_recovery | [1, 4, 2] | unavailable in committed evidence |
| 9 | `AgTi_2p3871A_periodic_PW_PBE_force_only` | unresolved | unknown | unavailable in committed evidence |
| 10 | `AgTi_2p55_periodic_PW_PBE_energy_force` | partial_recovery | [1, 1, 1] | unavailable in committed evidence |
| 11 | `v10_element_reference_Ag_fcc` | partial_recovery | unknown | unavailable in committed evidence |
| 12 | `v10_element_reference_Ti_hcp` | partial_recovery | unknown | unavailable in committed evidence |
| 13 | `v10_element_reference_Si_diamond` | partial_recovery | unknown | unavailable in committed evidence |
| 14 | `v10_element_reference_C_diamond` | partial_recovery | unknown | unavailable in committed evidence |

Nine `*_gap_*` rows declare GPAW PW-PBE / 500 eV / 1x4x2 / 0.1 eV smearing. AgTi 2.55 Å declares the same basis/XC/cutoff with Gamma and Fermi smearing 0.1 eV. Four elemental solid anchors declare PW-PBE / 500 eV / Fermi smearing 0.1 eV but do not give exact k meshes. These are stored declarations, not independently verified calculator inputs. No original GPAW/ASE or PAW setup version, SCF threshold/log or verified stopping record was recovered for any row. PBC/cells are known from identical serialized geometries; absent persistent IDs for 13 rows remain absent.

The sole unresolved row is `AgTi_2p3871A_periodic_PW_PBE_force_only`: despite its old tag, it has REF energy and forces, and a policy note says matching energy was added to an older force-only geometry. Its basis/XC/cutoff/k-grid/software/convergence cannot be confirmed from a config tag or global manifest. It should be A’s first source-record recovery or matched-relabeling review target.

## Original records and Git history

The V11 training file and its manifest entered this repository at `56b2370e0b2829ca4d6f9ba17a0a559c199240b0` (Add MACE v12 handoff package). V11’s manifest/README identify a V10 energy-calibration dataset; directory/version names are not a method certificate. The earliest imported training-file SHA256 equals its current SHA256. V9/V10 original datasets/runners/logs are not preserved by this handoff’s reachable history. The JSON supplies V11–V14 import commit IDs, blob IDs and file hashes for each row.

Searched 135 historical source/data/log path names and 130 selected historical blobs over 211 reachable commits. The original external paths recorded for the nine gap summaries, V9 parent train, AgTi 2.387/2.55 outputs/summaries and four elemental outputs are absent from that path history. Where a source SHA was recorded, no matching original bytes were found elsewhere among searched blobs. No large dataset/model is copied into this publication. The history inventory records all searched paths/blob SHA256 and the exclusion scope.

The first Ti 2.3 Å source is named `summary_recovered_from_extxyz.json`. Even recovering that summary alone would not provide independent original SCF evidence if it was reconstructed from the same label file. For all nine gap rows, raw original output filenames and direct output hashes are not known; the declared V9 parent file hash does not bind a particular original GPAW run. Six other rows have more specific original-output path/hash declarations, but absent bytes still prevent identity verification.

## Energy/force linkage and conflicts

For the nine gap rows and AgTi 2.55 Å, V11 serializes both REF fields and SinglePoint/native `energy`/`forces`; the REF values equal those native aliases exactly. This is internal duplication in one data file, not an independent archived calculation. REF energies differ from stored `free_energy` (about 0.120–0.244 eV for gap cells and 0.6353 eV for AgTi 2.55). The ledger records each difference. Do not accidentally relabel with free energy while mixing the existing native-energy targets; original getter/settings conventions remain unverified. The four elemental anchors and AgTi 2.387 lack such duplicate native output fields, so only REF lineage can be checked.

Known declarations mix 1x4x2 interface cells, Gamma AgTi clusters, and unspecified dense elemental meshes. Different cell shapes can justify different sampling; the mismatch is a numerical-convergence review requirement, not proof that a method conflict caused the model error. Every row lacks original convergence/setup/version linkage. Elemental anchors are solid-reference energies, not isolated-atom chemical potentials.

## Missing evidence and exact-geometry remediation

Each JSON/CSV row lists its declared original file path and expected hash where available, source split/frame/hash, original import commit, stored PBC/cell, present IDs, parameters with evidence levels, and exact missing materials. Recover original calculator input/script, complete SCF log with tolerances/stopping evidence, software/PAW setup versions, and hash-linked original energy/force output. Compare geometry and labels to the retained frame before upgrading status to confirmed. Similar geometry, a tag, a filename or a global statement is insufficient.

If original records cannot be recovered, the remediation ledger points to the existing draft `train.extxyz` SHA256 and exact frame index plus a geometry hash; no new geometry file or relabeling is created. A may select matched PW-PBE single points preserving stored atom order/IDs/cell/PBC and recording energy convention and k-point convergence for the actual cell. For rows without original IDs, retain the stored order and document mapping rather than invent historical identifiers. Prioritize AgTi 2.387’s energy/force binding, then phase-reference numerical provenance, then gap-run convergence recovery.

Conditional counts/ranks are analysis only:

| Scenario | frames | composition rank |
|---|---:|---:|
| draft_all_30 | 30 | 4/4 |
| exclude_all_15_previous_unknown | 15 | 3/4 |
| exclude_only_unresolved_AgTi2p387 | 29 | 4/4 |
| exclude_four_element_anchors | 26 | 4/4 |

Blanket removal of all 15 previous unknown rows reduces rank to 3/4, so it is not a safe automatic fix for composition identifiability. Removing only the unresolved AgTi row would leave rank 4/4, but it does not solve the 14 partially documented references. These scenarios do not alter any dataset or authorize training; A decides exclusions/recovery/relabeling and integration.

Auxiliary packaging note: `file_manifest.json` has a literal trailing `\n` after its JSON object and fails strict JSON decoding. This trace records its raw SHA and inspects only the first object as an auxiliary inventory; the global manifest was not repaired or treated as original calculator evidence.

Artifacts: `per_frame_ledger.json`, `per_frame_ledger.csv`, `history_search_inventory.json`, this report and `SHA256SUMS.txt`. The JSON carries all evidence/input hashes; the checksum list covers current evidence files plus these compact artifacts and omits itself. Historical blobs are addressed by Git OID plus SHA256 in the inventory, without duplicating them.
