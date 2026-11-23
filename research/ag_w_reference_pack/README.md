# Ag-W research reference pack

This is a curated starter pack for the Ag:W = 1:1 mass-ratio project. It is deliberately kept separate from the legacy Ag-Ti3SiC2 research. No production MD run has been launched from these materials.

## 1. Highest-priority Ag-W potential papers

### Dai et al., 2006 — Extended Finnis-Sinclair (EFS), Ag-refractory systems
- Citation: X. D. Dai, Y. Kong, J. H. Li, B. X. Liu, *J. Phys.: Condens. Matter* 18 (2006) 4527-4542.
- DOI: https://doi.org/10.1088/0953-8984/18/19/008
- Author-uploaded full-text landing page: https://www.researchgate.net/publication/231128055_Extended_Finnis-Sinclair_potential_for_bcc_and_fcc_metals_and_alloys
- Key result: Table 6 reports explicit Ag-W cross parameters. Tables 1 and 4 report the pure W and Ag EFS parameter sets. The extracted values are recorded in `parameters/Dai_2006_EFS_AgW_parameters.json`.
- Critical caveat: the paper's EFS model predicts W melting at 4500 K versus 3695 K experimentally (~21.8% high), while Ag is 1175 K versus 1235 K. Do not use its high-temperature W response quantitatively without validation.

### Zhang et al., 2004 — direct Ag-W EAM paper
- Citation: R. F. Zhang, Y. X. Shen, H. R. Gong, L. T. Kong, B. X. Liu, *J. Phys. Soc. Jpn.* 73 (2004) 2023-2027.
- DOI: https://doi.org/10.1143/JPSJ.73.2023
- JPSJ issue page: https://www.jstage.jst.go.jp/browse/jpsj/73/7/_contents
- National Diet Library record: https://ndlsearch.ndl.go.jp/books/R000000004-I7004615
- Abstract states that an n-body Ag-W EAM was derived using ab initio cohesive energies/lattice constants and experiment, then used in MD for metastable phase selection.
- Access note: the journal record is authentication-gated; no openly downloadable parameter file was found during this collection pass.

## 2. Ag/W interface and microstructure references

### He et al., 2018 — Ag/W interface DFT
- Citation: C. He, M. Cheng, M. Zhang, W. X. Zhang, *J. Phys. Chem. C* 122 (2018) 17928-17935.
- DOI: https://doi.org/10.1021/acs.jpcc.8b05920
- Direct author-hosted PDF: https://gr.xjtu.edu.cn/html/documents/upload/20776/jpcc2018.pdf
- Publisher page: https://pubs.acs.org/doi/10.1021/acs.jpcc.8b05920
- Useful interface orientations: Ag(100)/W(100) and Ag(111)/W(110); the paper reports work of separation, interface energy and electronic properties. Use these as DFT validation targets, not as a potential by themselves.

### Jiao et al., 2018 — nanocrystalline Ag-W stability
- Open repository record: https://ira.lib.polyu.edu.hk/handle/10397/106467?mode=full
- Direct accepted-manuscript PDF: https://ira.lib.polyu.edu.hk/bitstream/10397/106467/1/Jiao_Nanocrystalline_Ag-W_Alloys.pdf
- Scope caution: this study is dilute W in nanocrystalline Ag (0.3-1.9 at.% W), not the 1:1 mass-ratio bicontinuous composite. It is useful for grain-boundary segregation mechanisms, not as a direct composition/force-field validation.

## 3. Pure-element reference potentials (downloaded text files in this folder)

- `potential_files/reference_pure_element/Ackland_1987_Ag.eam.fs`
  - NIST IPR source: https://www.ctcms.nist.gov/potentials/Download/1987--Ackland-G-J-Tichy-G-Vitek-V-Finnis-M-W--Ag/2/Ag.eam.fs
  - Original Ag EFS pure-element reference.
- `potential_files/reference_pure_element/Ackland_Thetford_1987_W.eam.fs`
  - NIST IPR source: https://www.ctcms.nist.gov/potentials/Download/1987--Ackland-G-J-Thetford-R--W/2/Ackland_Thetford_1987_W.eam.fs
  - W Finnis-Sinclair pure-element reference, adapted for LAMMPS eam/fs in 2025.

**Do not combine these two pure-element files as if they formed an Ag-W alloy potential.** They contain no fitted Ag-W cross interaction and are included only as pure-element baselines/reference implementations. They are also not the pure-element terms from Dai et al.'s 2006 EFS parameterization.

## 4. LAMMPS implementation reference
- EAM / EAM-FS documentation: https://docs.lammps.org/pair_eam.html
- LAMMPS documents `pair_style eam/fs` and the required multi-element FS table format. The Dai EFS analytic form must be tabulated/implemented consistently before a LAMMPS input can use it.

## 5. Next validation gates before production MD
1. Reconstruct the Dai-2006 Ag-W EFS functions from the tabulated parameters and verify units/sign conventions.
2. Test pure Ag and pure W lattice constants, cohesive energies, elastic constants and equations of state.
3. Test Ag and W melting/coexistence; the W melting discrepancy is a known warning.
4. Test Ag/W interface energies and work of separation for more than one orientation against He et al. (2018) DFT.
5. Test immiscibility/formation enthalpies and surface diffusion against Zhang et al. (2004) and available DFT.
6. Only then build the 1:1 mass-ratio QSGS morphology and run minimization/300 K equilibration.

## 6. Scope and provenance
- The two `.eam.fs` files are sourced from the NIST Interatomic Potentials Repository's public file listings, mirrored in the lmhale99/potentials-library repository.
- The Dai parameter JSON is a transcription of the paper's Tables 1, 4 and 6 and equations (1)-(5), for audit/reconstruction. It is not a fitted or validated potential file.
- No legacy Ag-Ti-Si-C file has been overwritten.


## 7. PDF access/download status
- The 2018 He et al. interface paper has an author-hosted full PDF link above.
- The 2018 Jiao & Schuh manuscript has an open institutional-repository PDF link above.
- The 2006 Dai et al. full text is available through the author-uploaded ResearchGate landing page above; its direct file endpoint did not resolve in this collection pass, so the PDF itself is not mirrored here.
- The 2004 Zhang et al. JPSJ paper is marked authentication-gated by the journal/ResearchGate. Use the publisher or author-request route; no openly downloadable PDF or public parameter file was verified.
- This repository branch contains the two pure-element reference `.eam.fs` text files and a machine-readable transcription of the Dai 2006 EFS parameters. It does **not** contain the publisher PDFs, and it does **not** yet contain a validated binary Ag-W alloy potential.
