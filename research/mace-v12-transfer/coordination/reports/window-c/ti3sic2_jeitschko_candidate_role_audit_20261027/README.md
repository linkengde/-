# Jeitschko 1967 Ti₃SiC₂ candidate-label role audit

This read-only audit checks the four 10×10×2 bulk displacement outputs from the reconstructed full-text Fourier geometry. It confirms the 1967 single-crystal source family is distinct from the COD 9009647/Kisi 1998 source used by V19 dev3, reruns all four registered result verifiers and per-result checksum lists, and reads only geometry columns from the frozen 36-row training and 3-row development split. It does not inspect test/sealed data or change the approved dataset.

All four labels pass the result and checksum checks and preserve native extrapolated energy, finite-width free energy, and 0.10 eV-width forces separately. The four rows form one parent-source displacement panel; they are not four independent source families. They remain `train_acquisition_candidate_pending_audit`. No training target convention has been frozen for new data, and the existing train36 energy-key lineage review still applies.

The frozen train36 and dev3 geometries are 52-atom Ag-containing cells, while these candidates are 12-atom Ag-free Ti₃SiC₂ cells. The same-composition/cell/PBC comparator therefore finds no directly comparable frames in the official split. This is not a full supercell-equivalence audit across every historical archive. The report does not assign a train/dev role or authorize integration.

Run `audit_candidates.py` from the repository root with the GPAW project environment. It writes `candidate_role_audit.json`. `SHA256SUMS.txt` covers this report bundle, excluding itself.
