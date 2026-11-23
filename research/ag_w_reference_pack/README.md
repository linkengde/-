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
- `papers/Extended_Finnis-Sinclair_potential_for_bcc_and_fcc.pdf` is the uploaded 17-page Dai et al. (2006) paper. Its title, DOI, equations (1)–(5), and Tables 1, 4, 6, and 7 were checked against the saved PDF.
- `papers/He_2018_AgW_Interface_Stability.pdf` is the uploaded He et al. (2018) interface paper. Table 4 provides DFT comparison values for Ag(100)/W(100) and Ag(111)/W(110).
- `papers/Jiao_Nanocrystalline_Ag-W_Alloys.pdf` is the uploaded 35-page accepted manuscript of Jiao and Schuh. It is an author manuscript with a CC-BY-NC-ND notice. It studies 0.3–1.9 at.% W in nanocrystalline Ag and is background on segregation, not a 1:1 composite validation.
- `papers/SHA256SUMS.txt` records hashes for the three uploaded PDFs.
- The 2004 Zhang et al. JPSJ paper is marked authentication-gated by the journal/ResearchGate. Use the publisher or author-request route; no openly downloadable PDF or public parameter file was verified.
- The He et al. paper links its supporting-information PDF at `https://pubs.acs.org/doi/suppl/10.1021/acs.jpcc.8b05920/suppl_file/jp8b05920_si_001.pdf`. Its listed contents are phonon spectra of Ag(111)/Ni(111), not Ag/W. The ACS host was not reachable under the running environment network policy, so the SI file has not been retrieved.
- The 17-page Dai article includes no separate SI citation or file. No Dai SI was supplied with the uploaded article.

## 8. Dai EFS to LAMMPS conversion and checks
- `scripts/convert_dai_efs_to_eamfs.py` tabulates the paper's equations as a LAMMPS `eam/fs` setfl potential in `potential_files/derived/Dai_2006_EFS_AgW.eam.fs`.
- Rebuild and check it from the repository root with `python3 research/ag_w_reference_pack/scripts/convert_dai_efs_to_eamfs.py --force` followed by `python3 research/ag_w_reference_pack/scripts/validate_dai_efs_eamfs.py`. The validation script accepts `--lammps /path/to/lmp` when LAMMPS is outside `/workspace/.venvs/ag-qsgs/bin/lmp`.
- `parameter_audit.md` records the primary-source equations, page numbers, parameter tables, and a unit correction. The original PDF prints A in eV Å⁻¹. A previous JSON revision incorrectly said eV¹ᐟ² Å⁻¹. The potential values and generated table did not change because of that metadata correction.
- `scripts/validate_dai_efs_eamfs.py` runs pure Ag and W bulk checks, the AgW B2 and Ag3W L12 structures from Table 7, and a small Ag(100)/W(100) BCC-site interface separation scan. Results are in `validation/dai_2006_efs_validation_report.json`.
- Pure Ag and W lattice constants and cohesive energies reproduce Tables 4 and 1 within 0.0005 eV/atom and 0.00003 Å. The B2 and L12 cohesive energies are close to the Table 7 VASP values, but the B2 result differs from Table 7's own EFS value by 0.0859 eV/atom. That inconsistency remains unresolved.
- The rigid interface scan gives a maximum work of separation of 2.923 J/m² versus the DFT value 2.80 J/m² in He et al. Table 4. The scan uses the paper's in-plane orientation relationship and BCC-site registry but does not relax the interface or reproduce a converged DFT comparison.
- Dai et al. report a W melting point of 4500 K versus 3695 K experimentally. This EFS set is not ready for production or high-temperature MD until the alloy benchmark discrepancy and remaining thermodynamic/interface tests are resolved.
- The reference pure-element Ackland files remain separate. They were not concatenated into the Ag–W potential.
