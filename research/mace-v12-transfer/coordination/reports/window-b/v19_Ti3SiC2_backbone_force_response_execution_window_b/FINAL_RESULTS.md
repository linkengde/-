# Ti₃SiC₂ backbone force-response — final results

All four registered fixed-cell displacement points completed sequentially. Each archive passed the GPAW result verifier and SHA-256 inventory. Large `.gpw` checkpoints were kept local and removed only after the process exited and the compact archive passed verification.

## Method

PBE plane-wave GPAW, 500 eV cutoff, 0.10 eV Fermi–Dirac smearing, 10×10×4 k-points, PBC TTT, four MPI ranks; SCF energy/density/eigenstate thresholds 1e-5. GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1. Method-manifest SHA-256: `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`.

Native extrapolated energies (`force_consistent=False`) and free energies are kept separate. Energy differences below are relative to the unchanged 48-atom parent `Ti3SiC2_baseline_k10x10x4_v19`.

## Baseline

- Native energy: **-389.809242860 eV/cell**; free energy: **-390.061790365 eV/cell**.
- SCF: 34 iterations; elapsed 3.26 h; parent force RMS 0.056986 eV/Å; maximum force 0.087864 eV/Å.

## Per-point results

`Δforce RMS/max` is the RMS and maximum over atoms of the norm of each force-vector difference from the parent, in eV/Å.

| Point | Native E (eV/cell) | Δ native E (meV/cell; meV/atom) | Free E (eV/cell) | Δforce RMS / max (eV/Å) | SCF / time | Archive |
|---|---:|---:|---:|---:|---:|---|
| `Ti4f_id0003_z_m020A` | -389.806901983 | +2.341; +0.0488 | -390.059155012 | 0.05641 / 0.35323 | 35 / 2.90 h | PASS |
| `Ti4f_id0003_z_p020A` | -389.804590824 | +4.652; +0.0969 | -390.057433006 | 0.05525 / 0.34562 | 34 / 2.83 h | PASS |
| `C4f_id0009_z_m020A` | -389.804361797 | +4.881; +0.1017 | -390.056970752 | 0.04825 / 0.30574 | 37 / 3.18 h | PASS |
| `C4f_id0009_z_p020A` | -389.808038806 | +1.204; +0.0251 | -390.060525276 | 0.04835 / 0.30428 | 48 / 4.20 h | PASS |

## Displaced-site z forces and symmetric response

Finite-difference response is `−(Fz(+0.02 Å) − Fz(−0.02 Å)) / 0.04 Å`. This is a fixed-cell response for the selected sites, not a bulk elastic constant.

| Site | Fz at −0.02 Å (eV/Å) | Fz at +0.02 Å (eV/Å) | Symmetric response (eV/Å²) |
|---|---:|---:|---:|
| Ti 4f, local ID 3 | 0.30827444 | -0.39057549 | 17.471248 |
| C 4f, local ID 9 | 0.39360453 | -0.21641497 | 15.250488 |

## Input and archived force/energy hashes

| Label | Input SHA-256 | Extxyz SHA-256 |
|---|---|---|
| `Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19` | `5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf` | `95f3803ec53eaa5a89c6b0117446782e66f4060dbcb01588a7bb5488b24220ce` |
| `Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19` | `5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b` | `fa16c172aa4ae926b38347c230e60f8bb93dd79a2e2d66ec2481de79579d1265` |
| `Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19` | `19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f` | `ac021467637026f8728b9b25228b061f8485dda42577b5ce8ae5a8ab4b134850` |
| `Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19` | `6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72` | `3c0d7f82a748ea93e3afc1d95da44e28f79e11e3ef3959e5746c1463c7c66229` |

All archive-file hashes are included per label in `final_results.json`; each label’s archive `SHA256SUMS.txt` passed `sha256sum -c`.

## Scope and interpretation

These are fixed-cell bulk backbone response candidates sharing one parent and selected site family. They are not independent validation data and were not integrated into training. No interface DFT, model training/inference, or dev/test/holdout access was performed.
