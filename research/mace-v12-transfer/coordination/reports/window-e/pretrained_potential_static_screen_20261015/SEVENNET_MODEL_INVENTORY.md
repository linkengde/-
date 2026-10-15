# SevenNet-MF-0 inventory

## Frozen candidate

- Official project: [MDIL-SNU/SevenNet](https://github.com/MDIL-SNU/SevenNet), checked out at `209339d4eb4c1717f906c4cb715629159f00f88c`.
- Package: SevenNet `0.13.0`; code license MIT.
- Model: `SevenNet-MF-0`, official checkpoint `checkpoint_sevennet_mf_0.pth`, SHA-256 `81791329b37d445f46b531578c182c41792d98c7814222c9e5dde276402225fd` (10,322,158 bytes). The checkpoint itself was not copied into the project report. The repository carries an MIT license for the code; I found no separate checkpoint-only license statement, so this report does not assume one beyond the published repository/model distribution terms.
- Official model card identifies training on Materials Project PBE(+U) and r²SCAN crystal-relaxation trajectories. PBE modal selected for this static screen. The source states it supports these two fidelity modes.
- Runtime reports 89 element types and confirmed Ag/Ti/Si/C support; model cutoff 5 Å; inference outputs float32.

## Frozen DFT checks

- 48-atom Ti₃SiC₂ smoke: finite energy and forces; force-vector RMSE `0.0567594 eV/Å`; maximum force error `0.0943695 eV/Å`.
- 135-atom full interface: force-vector RMSE `0.262824 eV/Å`, P95 `0.541323`, maximum `0.654547 eV/Å` at source atom 76 (Si). Element RMSEs: Ag `0.228824`, Ti `0.256138`, Si `0.432287`, C `0.156909 eV/Å`.
- TSC-only 114 atoms: `0.406160 eV/Å`; Ti `0.413910`, Si `0.635487`, C `0.189725 eV/Å`. Ag-only 21 atoms: `0.069395 eV/Å`. All values are frozen predictions against the cited existing DFT labels.
- Component-subtracted interface-force RMSE is `0.418405 eV/Å`; the model predicts net normal interaction force sums Ag `+3.34694` and backbone `−3.34694 eV/Å`, whereas DFT gives Ag `−0.17101` and backbone `+0.16653 eV/Å`.

## Provenance and limitations

The Materials Project PBE(+U)/r²SCAN training provenance creates possible overlap with the Ti₃SiC₂ bulk control; exact or family overlap with local reference data was not established. The 135-atom interface is not a blind test either. Model weights stayed frozen.
