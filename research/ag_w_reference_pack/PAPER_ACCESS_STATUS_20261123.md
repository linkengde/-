# Ag-W paper access check — 2026-11-23

This file records verified public access status for the priority papers and distinguishes public access links from local PDF copies.

## 1. Dai et al. (2006) — Extended Finnis–Sinclair potential for bcc and fcc metals and alloys

- DOI: https://doi.org/10.1088/0953-8984/18/19/008
- Journal: Journal of Physics: Condensed Matter 18 (2006), 4527–4542.
- Status: **Free author-uploaded full text found.**
- Author-uploaded ResearchGate page: https://www.researchgate.net/publication/231128055_Extended_Finnis-Sinclair_potential_for_bcc_and_fcc_metals_and_alloys
- The page explicitly labels the article PDF as available and the full text as public; the file was uploaded by co-author Yi Kong. Open the page and use "Download full-text PDF".
- A user-supplied copy of the 17-page paper is committed in `papers/Extended_Finnis-Sinclair_potential_for_bcc_and_fcc.pdf`. Its SHA-256 is recorded in `papers/SHA256SUMS.txt`.
- A second PDF mirror was also found, but the author-uploaded ResearchGate copy is the preferred access route.

## 2. Zhang et al. (2004) — Atomistic Modeling of Metastable Phase Selection of a Highly Immiscible Ag-W System

- DOI: https://doi.org/10.1143/JPSJ.73.2023
- Journal: Journal of the Physical Society of Japan 73(7) (2004), 2023–2027.
- Status: **No free public full-text PDF verified.**
- J-STAGE issue listing marks this article as restricted/authentication required: https://www.jstage.jst.go.jp/browse/jpsj/73/7/_contents/-char/ja?from=1
- ResearchGate record says no full text is available and offers "Request full-text PDF": https://www.researchgate.net/publication/252680752_Atomistic_Modeling_of_Metastable_Phase_Selection_of_a_Highly_Immiscible_Ag-W_System
- National Diet Library metadata/abstract: https://ndlsearch.ndl.go.jp/books/R000000004-I7004615
- DOI is provided as the reliable citation/access identifier. If needed, request a copy from the authors or use institutional/library access.

## 3. He et al. (2018) — Interfacial Stability and Electronic Properties of Ag/M (M = Ni, Cu, W, and Pd) and Cu/Cr Interfaces

- DOI: https://doi.org/10.1021/acs.jpcc.8b05920
- Journal: Journal of Physical Chemistry C 122 (2018), 17928–17935.
- Status: **Free author-hosted PDF found and verified.**
- Direct PDF: https://gr.xjtu.edu.cn/html/documents/upload/20776/jpcc2018.pdf
- The PDF opens and identifies the exact title, authors, journal citation and DOI. It contains DFT results for Ag(100)/W(100) and Ag(111)/W(110) interfaces.
- A user-supplied copy is committed in `papers/He_2018_AgW_Interface_Stability.pdf`; its SHA-256 is recorded in `papers/SHA256SUMS.txt`.

## Download note

The public access links were verified in the earlier paper-access check. In this task, the Dai and He PDFs were supplied as chat attachments and those exact files were committed to the branch. The current running environment could not fetch the ACS supplementary-information file. Its direct link and contents are documented in the Ag–W reference pack README. The Jiao manuscript was also user supplied and is covered by its hash entry.

## 4. Additional user-supplied papers and TTM data — 2026-11-24

- `papers/Vijayakumar_1988_CuW_AgW_AuW_phase_diagrams.pdf`: user-supplied scan of H. Vijayakumar, A. H. Sriramamurthy, and S. V. Nagender Naidu, “Calculated Phase Diagrams of Cu-W, Ag-W and Au-W Binary Systems,” *CALPHAD* 12(2), 177–184 (1988). The article itself says there were no experimentally determined thermodynamic/phase-diagram data for Ag-W and builds a calculated diagram with assumptions. It is not an experimental phase diagram.
- `papers/Waseda_Ohtani_1974_Liquid_noble_transition_metals_Xray.pdf`: user-supplied paper, *Physica Status Solidi (b)* 62, 535–546. It reports X-ray structure for liquid Ag and other metals, not liquid W or Ag-W alloy structure.
- `papers/Torre_2021_Nanocrystalline_W_alloy_thermodynamic_stability.pdf`: user-supplied *Materials* 14, 7179 article, DOI 10.3390/ma14237179. It reports W-Ag segregation and demixing under the study's nanocrystalline mechanical-processing/annealing conditions.
- `papers/Li_2024_Electronic_interface_thermal_conductance_metal_nanofilms.pdf`: user-supplied *International Journal of Heat and Mass Transfer* 229, 125717 paper, DOI 10.1016/j.ijheatmasstransfer.2024.125717. It measures Au/Pt and Au/Al electronic interface conductance, not Ag/W; it highlights sensitivity to interface state/oxide.
- The user-provided XTANT-3 Ag/W Ce, G, and K data were copied unchanged to `ttm_data/XTANT3_TTM_parameters_v3/` with the accompanying source README and license. `ttm_data/XTANT3_AG_W_DATA_AUDIT_20261124.md` records their contents and model limitations; hashes are stored in the TTM data folder.
- The XTANT source README cites Medvedev, Akhmetov & Milov's 2023 arXiv preprint for electronic conductivity data. Downloading that PDF returned HTTP 403 in this environment; the uploaded README and citation are retained, and no PDF copy was fabricated.
- `papers/NASA_TM_87568_Microgravity_Science_Applications_Program_Tasks_1985.pdf` is the user-uploaded NASA TM 87568 program-task compilation. Its “Measurement of the Properties of Tungsten at High Temperatures” section reports a W heat of fusion and liquid Ag heat capacity/heat of fusion; the filename does not represent a separate standalone Margrave paper.
- The uploaded `cdc_219421_DS1.pdf` has the same SHA-256 as the already archived Pankratz, *Thermodynamic Properties of Elements and Oxides*, Bulletin 672 PDF, so the duplicate was not stored twice.
