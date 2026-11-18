# Ag–Ti₃SiC₂ MTP/ACE fit report

Date: 2026-11-18

## Data basis

This fit uses the repository's MACE-fitting DFT data, not Round7 data. The local exploratory split is 71 training frames and 7 development frames (2,823 training atoms; 204 development atoms). No blind test set is available.

| Repository source | Role in this fit | Selection note |
| --- | --- | --- |
| MACE v19 | 36 train, 3 development | All 3 development frames retained. Follow-up provenance review found one V17 `Ti_gap_2p3` Gamma relabel embedded in the 36 training frames; its SCF converged, but the inherited source mesh was 1×4×2 and k-point convergence remains unresolved. |
| MACE v25 | 19 train, 4 development | Non-official Ti–Si–C skeleton candidate. The 5 reused v21 bulk records already in its training split were not added again. |
| v21–v24 interface acquisitions | 16 train | Individual result hashes, SCF and PBE/PW method records checked. The same-geometry AgSi 8×8 control was excluded. |
| v10–v18 | Not separately appended | Older sets are nested/versioned. The v18 30-frame set is an exact label subset of the current candidate. v17/v18 are marked provisional; combining them would duplicate labels. |
| v20 | Excluded | Repository README identifies these as fixed-geometry k-point/smearing controls, not training data. |

Training and development have zero exact geometry-hash intersections, and the training split has no exact duplicate geometries. Energy targets are harmonized to GPAW native extrapolated total energy; finite-width free energies remain in the staged frame metadata. Forces are GPAW/ASE forces at 0.10 eV Fermi smearing. The candidate combines older Γ-point interface labels, 6×6×1 interface-response labels, and 10×10×2 Ti₃SiC₂ bulk labels. The split is exploratory and source-correlated: the v25/v21–v24 structures share the Jeitschko–Nowotny family, and the three v19 development labels share a COD 9009647 family.

Candidate hashes:

- `train.extxyz`: `d439734db7a9162f5c1b12508da2e2ec69d285295f900bc9df01aa841c34d417`
- `dev.extxyz`: `6efada57f83be14311c2fd8a6c452de6d152ef43306e2dcca157d0002a511f48`

## Fit results

Exploratory gates: energy RMSE ≤10 meV/atom; force-component RMSE ≤0.10 eV/Å; Si hotspot-center vector error ≤0.5 eV/Å.

| Potential | Training energy RMSE | Training force-component RMSE | Development energy RMSE | Development force-component RMSE | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| MTP10, force weight 0.05 | 3.61 meV/atom | 0.200 eV/Å | 5.35 meV/atom | 0.269 eV/Å | Energy pass; force fail |
| MTP10, force weight 0.20 | 10.18 meV/atom | 0.198 eV/Å | 11.89 meV/atom | 0.254 eV/Å | Energy and force fail |
| ACE, 5 Å cutoff, 1,616 functions, α=10⁻⁴ | — | — | 2.84 meV/atom | 0.175 eV/Å | Force fail |
| **ACE, 6 Å cutoff, 1,616 functions, α=10⁻³** | **0.35 meV/atom** | **0.0932 eV/Å** | **5.55 meV/atom** | **0.0930 eV/Å** | **Aggregate energy and force gates pass** |

The selected ACE model is `ace_mace_cutoff6/AgCSiTi_B1616_E10000_F1_a0.001.ace.yaml`; the LAMMPS Ctilde export is `ace_mace_cutoff6/AgCSiTi_B1616_E10000_F1_a0.001.yace`.

The aggregate ACE development scores pass, but errors vary by frame: the worst development force-component RMSE is 0.199 eV/Å on `bulk_random_dev_seed8204_rms0p06_v25`; the worst absolute per-atom energy error is 10.86 meV/atom on `AgTi_v19_independent_dev_01`. Thus the pass applies to the stated aggregate RMSE gates. If the same thresholds are intended per frame, the current ACE candidate does not pass every frame.

The maximum vector error over all Si sites in ACE development is 0.245 eV/Å. This is a conservative site proxy below 0.5, not the formal hotspot-center score: no registered hotspot-center IDs are present in this MACE data split.

## LAMMPS availability

- ACE `.yace` export ran in LAMMPS `pace` with 2 MPI ranks. On one 52-atom development geometry, LAMMPS and Python ACE differed by 3.1×10⁻⁶ eV in energy and 6.1×10⁻⁷ eV/Å force-component RMSE.
- MTP10 ran through the MLIP LAMMPS pair style with 2 MPI ranks on the same geometry.
- The downweighted ACE `.yace` loaded and evaluated the repository's 1,021,793-atom QSGS structure in LAMMPS on 4 MPI ranks. A single `run 0` completed in 35 s; LAMMPS reported 204–214 MiB per MPI rank. This measures load plus one static energy/force evaluation, not MD steps/s.
- The million-atom structure's maximum initial ACE force was 10.34 eV/Å. This establishes large-system execution only; no integration or heating was run, and it is not a physical MD qualification.
- A follow-up `run 0` force audit found 330,452 atoms (32.3%) above 0.5 eV/Å, 203,536 (19.9%) above 1 eV/Å, and 2,143 (0.21%) above 5 eV/Å. Counts above 0.5 eV/Å by species were Ti 129,055/392,277 (32.9%), Si 41,907/130,759 (32.0%), C 77,658/261,518 (29.7%), and Ag 81,832/237,239 (34.5%). These are model force magnitudes on the unrelaxed starting geometry, not fitting residuals.
- The force-audit geometry has a 2.088 Å minimum nearest-neighbor distance and no pair below 1.5 Å. Thus the high forces are not explained by gross atomic overlaps. They occur across the material-bearing z range. The distribution warrants local reference checks and relaxation stability work before any dynamics. Details and top-force sites are in `million_atom_static_ace/force_audit_summary.json`; the filtered force dump and reproduction input are alongside it.

## Limits and next data need

The baseline ACE candidate passes the aggregate energy and force gates, while the MTP fits do not. The Si-center requirement cannot be scored exactly until hotspot atom IDs are registered. The follow-up force-weight sensitivity fit reduces the current development errors, but it underfits the rare high-force `Ti_gap` training frame; its suitability depends on whether such near-gap states are in the target simulation domain. A filtered MTP10 sensitivity fit that removed that frame became substantially worse (78.5 meV/atom and 0.495 eV/Å on development), so the frame cannot simply be removed from the MTP training pool. There is no sealed test set, and neither fit has been validated for liquid Ag, melting, impact, or high-temperature MD.

All fit files are under `/workspace/potential-comparison/mace_fit_data_20261118`. The tracked repository `/workspace/-` remains unchanged.

## Follow-up: V17 rank-recovery label sensitivity

The training split contains `Ti_gap_2p3_rank_recovery_v17_01` (32 atoms). Its exact input/output archive checks passed, and the Gamma PBE/PW 500 eV, 0.10 eV smearing calculation converged in 50 SCF iterations. The repository's V17 review records that the inherited structure declared a 1×4×2 mesh, so the Gamma relabel is a valid converged Gamma calculation but not evidence of k-point convergence. The reference forces are unusually large: 3.10 eV/Å vector RMS (1.79 eV/Å component RMS), 5.68 eV/Å maximum vector magnitude. The baseline ACE per-frame force-component RMSE on this training configuration is 0.831 eV/Å.

A sensitivity fit scaled only this frame's force sample weight to 0.05; its energy weight and the other 70 training frames were unchanged. This alternate exploratory ACE fit reports:

- Training aggregate: 0.191 meV/atom energy RMSE; 0.0996 eV/Å force-component RMSE.
- Development aggregate: 4.55 meV/atom energy RMSE; 0.0555 eV/Å force-component RMSE.
- Worst development frame: 0.0989 eV/Å force-component RMSE; maximum absolute frame energy error: 8.97 meV/atom.
- Maximum vector error over all development Si sites: 0.206 eV/Å. This remains a proxy, since formal hotspot-center IDs are absent.
- The downweighted `Ti_gap` training frame itself remains poorly fit (0.914 eV/Å force-component RMSE). Use this candidate only when accurate dynamics through that rare high-force geometry are not required.

For comparison, omitting that single frame gives 4.72 meV/atom and 0.0506 eV/Å development RMSE, with every development frame below the stated energy/force per-frame limits. Both results show that this one unresolved, high-force row controls much of the baseline error; they do not establish that its label is invalid. The maximum force-weight scale in the tested grid that also met all development per-frame limits and both train/development aggregate gates was 0.05.

The selected downweighted `.ace.yaml` was converted to `.yace` and passed a 52-atom, two-rank LAMMPS zero-step compatibility check. Compared with Python ACE, force-component RMSE was 2.8×10⁻⁷ eV/Å; the energy difference from the rounded LAMMPS thermo log was 4.2×10⁻⁶ eV. A new 1×4×2 DFT calculation was attempted but stopped after two SCF iterations due to its high CPU cost; it produced no usable label. K-point convergence therefore remains unresolved.

The same one-frame removal was tested with MTP10 (`energy-weight=1`, `force-weight=0.05`, 1,000 iterations). Its training RMSE was 58.6 meV/atom and 0.563 eV/Å; development RMSE was 78.5 meV/atom and 0.495 eV/Å. This model fails the stated gates and is retained only as a data sensitivity result. The original full-data MTP10 fit remains the more useful MTP baseline (5.35 meV/atom and 0.269 eV/Å development RMSE), but its force error still fails.
