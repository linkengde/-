# Ag-W research reference pack

This is a curated starter pack for the Ag:W = 1:1 mass-ratio project. This means approximately 63.02 at.% Ag and 36.98 at.% W (not a 50:50 atomic mixture); see `composition_target.md` for the atom-count and volume-fraction conversions. It is deliberately kept separate from the legacy Ag-Ti3SiC2 research. No production MD run has been launched from these materials.

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
4. Compare constrained, zero-temperature Ag/W interface energy scans for Ag(100)/W(100) and Ag(111)/W(110) against He et al. (2018) DFT; see `validation/relaxed_interface_summary.md`. These are preliminary screenings, not a converged DFT reproduction.
5. Test immiscibility/formation enthalpies and surface diffusion against Zhang et al. (2004) and available DFT.
6. A provisional 1:1 mass-ratio QSGS geometry exists in `models/ag_w_qsgs_200A/`, but the original and v1/v2 repair candidates are not valid MD starts. V3 low-strain local slabs relaxed mechanically, while the strict DFT energy screen failed, especially for Ag(111)/W(110). The project permits qualitative graduate-level exploration with this known bias after a 3D candidate passes static geometry checks; quantitative interface energetics are not validated. Do not start 300 K MD or TTM/laser heating on the raw model.

## 6. Scope and provenance
- The two `.eam.fs` files are sourced from the NIST Interatomic Potentials Repository's public file listings, mirrored in the lmhale99/potentials-library repository.
- The Dai parameter JSON is a transcription of the paper's Tables 1, 4 and 6 and equations (1)-(5), for audit/reconstruction. It is not a fitted or validated potential file.
- No legacy Ag-Ti-Si-C file has been overwritten.


## 7. PDF access/download status
- `papers/Extended_Finnis-Sinclair_potential_for_bcc_and_fcc.pdf` is the uploaded 17-page Dai et al. (2006) paper. Its title, DOI, equations (1)–(5), and Tables 1, 4, 6, and 7 were checked against the saved PDF.
- `papers/He_2018_AgW_Interface_Stability.pdf` is the uploaded He et al. (2018) interface paper. Table 4 provides DFT comparison values for Ag(100)/W(100) and Ag(111)/W(110).
- `papers/Jiao_Nanocrystalline_Ag-W_Alloys.pdf` is the uploaded 35-page accepted manuscript of Jiao and Schuh. It is an author manuscript with a CC-BY-NC-ND notice. It studies 0.3–1.9 at.% W in nanocrystalline Ag and is background on segregation, not a 1:1 composite validation.
- `papers/SHA256SUMS.txt` records hashes for every archived PDF.
- The 2004 Zhang et al. JPSJ paper is marked authentication-gated by the journal/ResearchGate. Use the publisher or author-request route; no openly downloadable PDF or public parameter file was verified.
- The He et al. paper links its supporting-information PDF at `https://pubs.acs.org/doi/suppl/10.1021/acs.jpcc.8b05920/suppl_file/jp8b05920_si_001.pdf`. Its listed contents are phonon spectra of Ag(111)/Ni(111), not Ag/W. The ACS host was not reachable under the running environment network policy, so the SI file has not been retrieved.
- The 17-page Dai article includes no separate SI citation or file. No Dai SI was supplied with the uploaded article.

## 8. Dai EFS to LAMMPS conversion and checks
- `scripts/convert_dai_efs_to_eamfs.py` tabulates the paper's equations as a LAMMPS `eam/fs` setfl potential in `potential_files/derived/Dai_2006_EFS_AgW.eam.fs`.
- Rebuild and check it from the repository root with `python3 research/ag_w_reference_pack/scripts/convert_dai_efs_to_eamfs.py --force` followed by `python3 research/ag_w_reference_pack/scripts/validate_dai_efs_eamfs.py`. The validation script accepts `--lammps /path/to/lmp` when LAMMPS is outside `/workspace/.venvs/ag-qsgs/bin/lmp`.
- `parameter_audit.md` records the primary-source equations, page numbers, parameter tables, and a unit correction. The original PDF prints A in eV Å⁻¹. A previous JSON revision incorrectly said eV¹ᐟ² Å⁻¹. The potential values and generated table did not change because of that metadata correction.
- `scripts/validate_dai_efs_eamfs.py` runs pure Ag and W bulk checks, the AgW B2 and Ag3W L12 structures from Table 7, and a small Ag(100)/W(100) BCC-site interface separation scan. Results are in `validation/dai_2006_efs_validation_report.json`.
- `scripts/relax_and_compare_ag_w_interfaces.py` relaxes Ag(100)/W(100) and Ag(111)/W(110) interface cells with the central bulk-like layers held fixed, scans several lateral registries and initial gaps, and compares work of separation/interface energy with He et al. Table 4. Results are in `validation/relaxed_interface_report.json`.
- `scripts/audit_periodic_interface_equivalence.py` examines local coordination and per-atom energy around each of the two boundaries in the periodic interface cells. The DFT cell/strain comparison, B2 shell sum, and two-boundary findings are summarized in `validation/independent_reference_audit.md`.
- `scripts/build_ag_w_qsgs.py` generates a separate 200 Å cubic W-backbone/Ag-matrix QSGS candidate. Its dimensions, four growth controls, phase connectivity, count-based mass ratio and LAMMPS read check are documented in `models/ag_w_qsgs_200A/README.md`.
- `scripts/ag_w_ttm_source_preflight.py` maps the retained candidate heat-source profile onto actual Ag coordinates, builds bounded local Ag fractions, and reports on-state power/grid sensitivity without running MD. See `models/ag_w_qsgs_200A/ttm_source_preflight/TTM_SOURCE_PREFLIGHT.md`.
- Pure Ag and W lattice constants and cohesive energies reproduce Tables 4 and 1 within 0.0005 eV/atom and 0.00003 Å. The B2 and L12 cohesive energies are close to the Table 7 VASP values, but the B2 result differs from Table 7's own EFS value by 0.0859 eV/atom. That inconsistency remains unresolved.
- The rigid interface scan gives a maximum work of separation of 2.923 J/m² versus the DFT value 2.80 J/m² in He et al. Table 4. The scan uses the paper's in-plane orientation relationship and BCC-site registry but does not relax the interface or reproduce a converged DFT comparison.
- The relaxed interface scans minimize interfacial atoms while holding bulk-like cores fixed. The best Ag(100)/W(100) sampled registry gives Wsep = 2.894 J/m² and interface energy 1.127 J/m² versus DFT 2.80 and 0.72 J/m². Ag(111)/W(110) gives Wsep = 3.222 and interface energy 0.127 J/m² in its best sampled registry versus DFT 2.45 and 0.86 J/m². Full results and constraints are in `validation/relaxed_interface_summary.md` and `validation/relaxed_interface_report.json`; these checks are preliminary and do not clear the potential for production MD.
- The B2 energy from a direct evaluation of Dai's EFS equations agrees with the LAMMPS `eam/fs` initial energy to about 1e-12 eV/atom, while still differing by 0.0859 eV/atom from the EFS value printed in Table 7. The reported discrepancy is retained pending an explanation from the source calculation/convention.
- An independent B2 lattice-shell sum now agrees with the radial-neighbor and LAMMPS values to less than 4e-14 eV/atom. Table 7 does not document enough detail to reproduce its EFS energy convention exactly; see `validation/independent_reference_audit.md`.
- The two periodic Ag(100)/W(100) boundaries have matching local signatures. The two Ag(111)/W(110) boundaries are close but not identical; the current periodic result is their average. He et al. used one interface in a vacuum slab and reported 3.44% mismatch, while our fixed-cell strains differ. These comparisons remain preliminary.
- The new QSGS candidate has 178,435 W and 304,036 Ag atoms, with an Ag/W mass ratio of 0.999766. Its largest W and Ag components both wind through x/y and span z, but the raw perfect-lattice assignment has an Ag-W distance of 0.08 Å. An EFS static-minimization preflight still leaves a 0.37 Å pair and 2.17 eV/Å maximum force; no dynamics was started. See `models/ag_w_qsgs_200A/overlap_relaxation/preflight_report.md`.
- The v1 repair audit and v2 interface-reconstruction audit have produced no acceptable MD initial structure. V2 tested conflict-graph pruning plus EFS, full and sequential geometry projections, and same-phase lattice springs plus EFS. The pruned route deleted 49,669 atoms (10.3%) and remained unconverged after 2,700 EFS CG steps (1.969 eV/Å maximum force); the other approaches retained large short-contact or force failures. The raw geometry and failure records remain intact. See `models/ag_w_qsgs_200A/interface_reconstruction_v2/RECONSTRUCTION_V2_AUDIT.md`.
- V3 built one-interface low-strain slabs and completed 360 zero-temperature EFS relaxations across all sampled layer terminations, translations, and gaps. All met the force gate without atom loss or same-phase collapse. The best Ag(100)/W(100) result is Wsep/γint = 2.685/1.021 J/m² versus DFT 2.80/0.72; Ag(111)/W(110) gives 3.287/0.033 versus 2.45/0.86. This fails the frozen quantitative DFT screen, but is mechanically usable as a qualitative local-interface candidate. It does not fix the full QSGS geometry. See `models/ag_w_qsgs_200A/interface_reconstruction_v3/LOCAL_INTERFACE_V3_AUDIT.md`.
- The legacy Ag-only source preflight reports 9,784.984 eV/ps only for the raw coordinates. It must not be reused for a reconstructed geometry; recompute the local Ag mask and power only after static acceptance. No 300 K MD, TTM, or laser heating has been run.
- TTM remains the requested heating method. The installed `fix ttm/mod` is incompatible with the free z boundary; an Ag-only finite-volume reference solver exists elsewhere in the workspace, and the Ag/W pure-element source tables are now archived. A free-surface Ag/W solver, mixed-phase/interface closure, and pulse definition are still needed. See `TTM_HEATING_READINESS.md`.
- Dai et al. report a W melting point of 4500 K versus 3695 K experimentally. This EFS set is not ready for production or high-temperature MD until the alloy benchmark discrepancy and remaining thermodynamic/interface tests are resolved.
- The reference pure-element Ackland files remain separate. They were not concatenated into the Ag–W potential.

## 7. User-supplied TTM data, additional phase/interface papers, and full-size status

- The overlay archive `research/ag_w_reference_pack_update_20261124.zip` collects the updated papers, Ag/W TTM tables, source notices, audits, and current full-size geometry failure artifacts. Its `.sha256` sidecar checks the archive; `BUNDLE_README.md` and the in-archive manifest explain overlay paths and contents. The older `research/ag_w_reference_pack_complete_20261123.zip` remains the original dated bundle.
- The Ag and W XTANT-3 tables supplied by the user (`Ce`, `G_e-ph`, and `K`) are archived under `ttm_data/xtant3_ag_w_source/` with their upstream README, `.gitattributes`, and GPL-3.0 file. The upstream README also limits use to non-commercial peaceful purposes; preserve both notices. The row counts, temperature ranges, hashes, and cross-source comparison are in `ttm_data/XTANT3_AG_W_DATA_AUDIT_20261124.md` and `ttm_data/XTANT3_AG_W_SOURCE_AUDIT.json`.
- The standalone XTANT `Ce` and `G` tables are not interchangeable with the `Ce` columns embedded in the K tables. W XTANT `G` is substantially below the Lin 2008 digitized curve over the shared range, so these remain separate sensitivity alternatives. The XTANT K table is `Te`-only; it does not provide a validated `ke(Te,Tl)` relation or Ag/W interface conductance.
- The Lin et al. (2008), Duffy et al. (2014), Medvedev & Milov (2020), Giri et al. (2023), and Tolias (2017) papers are archived and assessed in `TTM_PARAMETER_LITERATURE_REVIEW_20261124.md`. An attempt to fetch the arXiv PDF cited by the XTANT source README returned HTTP 403; its citation and URL are retained, and the user-supplied table files remain archived.
- Additional uploaded papers now include the calculated 1988 Cu-W/Ag-W/Au-W phase diagrams, liquid Ag X-ray structure (Waseda & Ohtani, 1974), W-Ag segregation/demixing study (Torre et al., 2021), and Au/Pt and Au/Al electronic interface conductance study (Li et al., 2024). The 1988 Ag-W diagram is a model based on assumptions, not experimental phase-diagram data; the 2024 conductance values are not Ag/W measurements.
- The PDF named `OneDrive 入门.pdf` was identified by its content as the Lin et al. (2008) paper. A duplicate Duffy copy was deduplicated; an unrelated NASA microgravity program report and uploaded semiconductor tables were not added to the Ag-W model inputs.
- A v3 full-size all-atom projection candidate retained atom count and mass ratio, but 10,000-step spring-free EFS minimization (with the requested bottom anchor) still failed the practical start gate (0.650 eV/Å maximum mobile force, 2.099 Å Ag–Ag contact, 2.267 Å Ag–W contact, and 1,604 atoms above z=200 Å). No temperature dynamics or TTM run was started. See `models/ag_w_qsgs_200A/interface_reconstruction_v3/qsgs_candidate/FULL_SIZE_STATIC_AUDIT.md`.
- The Ag–W project accepts the potential for qualitative graduate-level exploration with its limitations stated. The next step remains rebuilding the 3D atomic arrangement, then repeating static checks before starting 300 K or TTM dynamics. The source archive does not make the current mechanically failed geometry an MD initial state.
