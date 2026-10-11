# Ag potential source-document supplement

Date: 2026-10-11 UTC  
Owner: F  
Purpose: record metadata and content checks for the user-supplied NIST release note, paper, and non-PDF source bundle. This is provenance documentation; it is not a new scientific execution authorization.

## Integrity and source-file identity

The uploaded `WMH2006_Ag_四个非PDF官方文件.zip` is 395,848 bytes, SHA-256 `fd3670c0cf714c60b706c6b2d485ef6134b58b9fdb75aff86ca739ecf846d51c`. ZIP integrity testing passed. Its six flat entries are `Ag.eam.alloy`, `F_ag.plt`, `fag.plt`, `pag.plt`, `DOWNLOAD_REPORT.json`, and `SHA256SUMS.txt`. All four source-data entries match the embedded SHA-256 manifest:

| File | Bytes | SHA-256 |
|---|---:|---|
| `Ag.eam.alloy` | 780,452 | `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447` |
| `F_ag.plt` | 156,069 | `a1644949665b4811629aa931ff2f2467267312490547be67be91eb5d56ec8338` |
| `fag.plt` | 156,062 | `84ee5e60eb963e149207c2b1994a9e7839a2e15ba90b93ce666a6a24ec9e8ce8` |
| `pag.plt` | 156,066 | `2f0f2979d81083f72cf0de650cff174196e7aef8a6262cf7d0302dad639aec2a` |

The `Ag.eam.alloy` hash matches the parameter identity already frozen in the C execution pack. The PLT hashes match the user's supplied download manifest. This verifies byte identity within the supplied archive; it is not an independent network re-download or a digital signature check.

## Attached PDF identity

- `Ag06_releaseNotes_1.pdf`: 43,458 bytes, SHA-256 `262d9102f5a77572792b9a1d72f6c55c01372ec798f3321285d76246446a8c53`; one page, readable and unencrypted.
- `an-embedded-atom-potential-for-the-cu-ag-system-54i0083q28.pdf`: 271,182 bytes, SHA-256 `2a1e8ec6d59e21275e671bdfd1b54b9ebcc04ab85439954a0f7aa6612a23bd57`; 17 pages, readable and unencrypted.

The paper title and DOI extracted from the document are “An embedded-atom potential for the Cu–Ag system,” P. L. Williams, Y. Mishin, and J. C. Hamilton, *Modelling and Simulation in Materials Science and Engineering* 14 (2006), 817–833, DOI `10.1088/0965-0393/14/5/002`.

The attached release note identifies `Ag.eam.alloy` as a 4 February 2009 C. A. Becker conversion from Yuri Mishin's three x,y PLT files, and describes its cubic-spline conversion and SOLD/LAMMPS comparisons. It names the 2006 Williams–Mishin–Hamilton paper as the potential's original reference and points to the NIST Interatomic Potentials Repository. It documents LAMMPS use as `units metal`, `atom_style atomic`, `pair_style eam/alloy`, and `pair_coeff * * Ag.eam.alloy Ag`.

The paper describes Ag fitting against experimental and first-principles data, notes that liquid configurations were not included in fitting, and reports transferability checks at high temperature. It gives a predicted Ag melting temperature of 1267 K and compares it with an experimental 1235 K. These are source-paper statements, not measurements made in this cloud task.

## Rights and release status

The attached files establish attribution and conversion lineage. They do not state an open-source license or an explicit right to redistribute the parameter file. The existing C `SOURCE_RELEASE.json` likewise records `open_source_license: NOT_EXPLICITLY_STATED` and `redistribution_status: UNVERIFIED_DO_NOT_REDISTRIBUTE_PARAMETER_FILE`. Project use authorization and public redistribution rights remain separate. Therefore, this supplement records hashes and citation metadata only; it does not include the PDFs or potential data.

No software was installed, no potential was loaded for calculation, no structure was generated, and no sampling or DFT was run. `C_RELEASE_CONFIG.json` was not changed; storage and warmup release remain false. A/B tasks and checkpoints were not accessed.
