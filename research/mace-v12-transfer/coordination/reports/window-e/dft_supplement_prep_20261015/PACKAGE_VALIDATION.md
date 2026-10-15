# Package validation — Ti₃SiC₂ P0 DFT input candidates

Validation date: 2026-10-15 (UTC workspace date).

## Integrity results

- Parent proposal SHA-256: `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`.
- Parent has 48 atoms (Ti24Si8C16), PBC `[true,true,true]`; atom order, coordinates, and cell match archived converged k10×10×4 and k6×6×4 result geometries exactly.
- k10×10×4 labelled output SHA-256: `75b083e225b8b20bbdde181c8e1847bae4b0052a854c6e9ee48885c940166415`; archived verification status: `PASS`.
- k6×6×4 labelled output SHA-256: `29e9ac254984b4a4ab6a22cbac7043c2cf9a826eb498422f1860a5798e72cee4`.
- E2 role-freeze SHA-256: `709d84b6f7205d3d437bae92dc1aaecd81278311f807b805e06c1b7d22acc27f`; TRAIN, DEV and TEST remain empty; this parent family is still diagnostic-only.
- All six candidate hashes match the manifest. Each retains atom order/IDs, composition and PBC, and has no energy, force, stress or virial labels. Transformation formula checks pass; maximum extxyz round-trip deviation from requested transformation is 4.85e-09 Å.

| Candidate | SHA-256 | Cell volume (Å³) |
|---|---|---:|
| `TSC48_BASAL_XY_SHEAR_m005` | `72be50f6e2fc3c7fcb610fadf520a7be3eadb8b14c2c7f27416da2fe7881ed7b` | 570.71030876 |
| `TSC48_BASAL_XY_SHEAR_p005` | `cd0e3d09e413d25a54dbaa38868106934b71e2d95c3c879c8346dbd0cff8f3a1` | 570.71030876 |
| `TSC48_SI_ID7_Z_m020` | `77ab3ae7dc59a41adea035003543fe19ae9c8deb8dd2e1413807ed25164e91f4` | 570.71030876 |
| `TSC48_SI_ID7_Z_p020` | `8222f1538719c859765f1ab48db483e4599d4f64ac2835f67dc5709ee860002d` | 570.71030876 |
| `TSC48_VOL_m02` | `1e8aeed059783a36e471873a3b3821369750677ae93ec5965d1d5760189ed683` | 559.29610259 |
| `TSC48_VOL_p02` | `a86ee4560676e1a9679729cd6c9721b4b4a875d41cc80862d6287a056ddb0db1` | 582.12451494 |

Volume pair: isotropic volume factors 0.98/1.02. Shear pair: `x′=x+γy`, γ=±0.005, determinant 1. Si pair: local ID 7 (Si) displaced by ±0.02 Å along Cartesian z. These are geometry-only proposals, not approved scientific settings.

## Runtime and software

- Existing private environment verified: Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, OpenMPI 5.0.7; GPAW reports MPI enabled and GPU disabled. No software was downloaded or installed because the needed DFT stack is already present.
- PAW SHA-256 values match the project baseline: `Ag.PBE.gz` `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`, `Ti.PBE.gz` `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`, `Si.PBE.gz` `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`, `C.PBE.gz` `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.
- E runtime script SHA-256: `e32d19eef429939bcb3d0aa2d365d063b6dca16b4b4052b02da226fd7247e43b`; pre-existing 4-rank PW-PBE smoke evidence SHA-256: `d06d0f479d041c1c2b731e889985a51a1669d36a098f195542615f653cf1601f` (prior evidence only; not rerun during this preparation).
- Current workspace limits: 4 CPU quota cores, 32 GiB memory, 22.43 GB free, no GPU. `/mnt/e-work` is absent and no dedicated/persistent DFT_WORK_ROOT is configured.

## Execution state

**No DFT or other scientific calculation was launched.** This package contains only six geometry candidates. Execution remains on HOLD until C freezes the parent-family role and exact parameters, issues authorization bound to the final manifest/structure SHAs, and approves a persistent DFT work root. The k-point centering, MPI rank count, wall-time cap and checkpoint policy also remain to be frozen.
