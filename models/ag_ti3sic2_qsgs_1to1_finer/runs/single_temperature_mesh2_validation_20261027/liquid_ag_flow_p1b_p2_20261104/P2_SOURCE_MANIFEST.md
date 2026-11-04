# P1b/P2 source and artifact manifest

Original user-provided PDFs are not copied into the repository; SHA256 values below identify the copies inspected in the workspace.

| Subject | Reference/source | DOI or URL | Inspected local PDF SHA256 | Audit status |
|---|---|---|---|---|
| Liquid Ag viscosity | Dubberstein et al. (2016), *Int. J. Thermophys.* 37, 100 | `10.1007/s10765-016-2104-7` | `34bb1a81fdabfd81c638dc33be45f4c8ea8c956fae22f80c818af4457c1c4963` | Table 4 and fit transcribed; NIST ThermoML page cross-check pending |
| Liquid Ag density | Abdullaev, Khairulin & Stankus (2020), *J. Phys.: Conf. Ser.* 1677, 012161 | `10.1088/1742-6596/1677/1/012161` | `f1cf9284a708b197f8233d67c06861117d004e1e3d0645b1239ee7d125cc6106` | Primary experimental paper and Eq. 2 inspected; prior density-source audit is cited |
| Ag surface tension | Morohoshi et al. (2008), *J. Japan Inst. Metals* 72, 708–713 | `10.2320/jinstmet.72.708` | `191ac3e067573b53049ad9c7c6b71e67512970fb4df061973dc4687547bba515` | Atmosphere-specific relations and validity intervals inspected |
| Ag/Ti₃SiC₂ wetting | Zhang et al. (2018), *Int. J. Miner. Metall. Mater.* 25(7), 810–816 | `10.1007/s12613-018-1629-0` | `f0da7547bb2e1dc6b8d107dd8de341245decc3edbed56c273033c22079efebcf` | Static 72° condition inspected; dynamic wall data absent |
| Dynamic Marangoni solver source | Basilisk `src/test/marangoni.c` | `https://basilisk.fr/src/test/marangoni.c` | Not acquired | Current network allowlist does not include `basilisk.fr`; no P1b run claimed |
| NIST ThermoML viscosity record | DOI-specific NIST ThermoML record | `https://trc.nist.gov/ThermoML/10.1007/s10765-016-2104-7.html` | Not acquired | Current network allowlist does not include `trc.nist.gov`; no NIST-page parse claimed |

P1 baseline: commit `66961178b383fa7c507bea6318f8672e1cb4c08e`. The environment configuration draft for restricted egress adds only `basilisk.fr` and `trc.nist.gov` to the existing custom allowlist. It is saved but requires review and publication before it can affect network access.

Geometry baseline used only for the applicability review: Ag70 QSGS mask SHA256 `03ef47b1dd534b3215cf03f42c8cf7e887dcdfcc185b33bab7453cd838c1faac`, as recorded in the prior read-only Stage 1 geometry audit. No geometry files were modified.
