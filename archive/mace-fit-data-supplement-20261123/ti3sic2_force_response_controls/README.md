# Ti₃SiC₂ force-response DFT candidates

This folder contains four existing, converged 48-atom GPAW PBE/PW500 eV force-response calculations at 10×10×4 k points and 0.10 eV Fermi smearing. They displace one Ti4f or C4f atom along z by −0.20 or +0.20 Å relative to the same Ti₃SiC₂ parent structure. Each source `verification.json` reports PASS across its convergence, input-hash, species/order, force, energy convention, and output checks; each source SHA256SUMS file was rechecked before copying. Full logs, compact extxyz labels, summaries, verification reports, and original per-case manifests are preserved.

Native extrapolated energy and force-consistent free energy remain separately recorded; the extxyz also preserves the PBE forces. Do not collapse these energy keys without fixing the training target convention.

All four structures derive from one crystal family and are **unassigned training candidates only**. Keep them together for dataset splitting; do not use them as independent development or blind-test data. They supplement framework force response but do not add liquid Ag, Ag/Ti₃SiC₂ interface, or hot Si-contact coverage. No training file or fitted model was changed.
