# MatterSim-v1.0.0-5M inventory

## Frozen candidate

- Official project: [microsoft/mattersim](https://github.com/microsoft/mattersim), exact v1.0.0 source commit `afe58a03d048f46693f65423ce05165302bbbb91`.
- Package/model license: MIT. Package version `1.0.0`; checkpoint `mattersim-v1.0.0-5M.pth`, SHA-256 `e3df9fa708725e3d453140646c7d1838324b347a3d1214cf1440522146f872b5` (91,176,875 bytes). The checkpoint itself was not copied into the project report.
- Model card lists the 5M release as 4.5M parameters, trained on 6M data points. It identifies a PBE variant; exact overlap with the local structures is unverified. Model card specifically cautions that v1 is designed for bulk and surface/interface quantitative predictions are not recommended without fine-tuning. Fine-tuning was not permitted or performed.
- Runtime confirms Ag/Ti/Si/C within supported atomic numbers (max Z=94), cutoff 5 Å and float32 inference.

## Frozen DFT checks

- 48-atom Ti₃SiC₂ smoke: finite energy and forces; force-vector RMSE `0.0456193 eV/Å`; maximum force error `0.0773450 eV/Å`.
- 135-atom full interface: force-vector RMSE `0.184893 eV/Å`, P95 `0.394815`, maximum `0.419911 eV/Å` at source atom 100 (Si). Element RMSEs: Ag `0.186646`, Ti `0.122566`, Si `0.337336`, C `0.150882 eV/Å`.
- TSC-only 114 atoms: `0.205856 eV/Å`; Ti `0.125983`, Si `0.407338`, C `0.142690 eV/Å`. Ag-only 21 atoms: `0.056571 eV/Å`. All values compare finite frozen outputs to existing converged DFT labels.
- Component-subtracted interface-force RMSE is `0.184591 eV/Å`; the model predicts normal force sums Ag `+2.47918` and backbone `−2.47918 eV/Å`, opposite the DFT signs (Ag `−0.17101`, backbone `+0.16653 eV/Å`). Composition-subtracted energy is `+7.26942 eV` above DFT native-total interaction increment, or `+6.41856 eV` relative to the separately reported force-consistent free-energy increment.

## Compatibility note

The official v1.0.0 code requires its release-pinned Torch 2.2.0, NumPy<2, e3nn 0.5, PyG 2.5.3, and ASE≥3.23. The existing ASE 3.29 environment was left unchanged; inference ran in a dedicated CPU virtual environment with ASE 3.23.0, Torch 2.2.0+cpu, NumPy 1.26.4, e3nn 0.5.0, and PyG 2.5.3. This isolates release compatibility instead of altering the existing MACE/GPAW setup.
