# P1b/P2 source and artifact manifest

All original user-provided PDFs remain outside the repository. Their prior inspection hashes are retained below. Newly acquired machine-readable source artifacts used for calculations are archived under `results/`.

| Subject | Reference/source | DOI or URL | Source hash/status |
|---|---|---|---|
| Dynamic Marangoni test | CoMPhy Lab Basilisk snapshot `v2026-07-20`, peeled commit `9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e` | `https://github.com/comphy-lab/basilisk-C/releases/download/v2026-07-20/basilisk-linux.tar.gz` | Archive SHA256 `4686294e8f649a424ca437a7fa94bfac9984c32bea5e9f33a27a5783cc692829`; downstream-patched Basilisk/Darcs snapshot, not untouched upstream |
| Official P1b source | `src/test/marangoni.c` within the pinned snapshot | `https://basilisk.fr/src/test/marangoni.c` | SHA256 `b2d2112e012ece62876c5010e3aeca3928c12c2d18ccd4b008e74afeb1d4294f`; archived under `results/p1b_official/` |
| Liquid Ag viscosity | Dubberstein et al. (2016), *Int. J. Thermophys.* 37, 100 | `10.1007/s10765-016-2104-7` | Original NIST ThermoML JSON SHA256 `80b4ab63347b4d5abc621d52e46b003a41069fe2061eaf474225532f9294f0fc`, 17,347 bytes; 9 Ag tuples matched the uploaded extracted subset exactly |
| User-provided extracted JSON | `nist_ag_9points_extracted_not_original.json` | Same DOI-specific NIST JSON URL | SHA256 `1381e3022e3a72c5270e734dfa2067c7e53e39c368bc89a014906ffae697bb16`; explicitly a subset, not the original source file |
| User-provided fetch script | `download_nist_original.sh` | NIST DOI-specific JSON URL | SHA256 `ec582dddaa4531c27760aa66d4f10f0c6ca33d06d0d4bd74d2aa9c47f93de33a`; inspected but not executed because original JSON was already available |
| Liquid Ag density | Abdullaev, Khairulin & Stankus (2020), *J. Phys.: Conf. Ser.* 1677, 012161 | `10.1088/1742-6596/1677/1/012161` | Prior audit retained; user PDF SHA256 `f1cf9284a708b197f8233d67c06861117d004e1e3d0645b1239ee7d125cc6106` |
| Ag surface tension | Morohoshi et al. (2008), *J. Japan Inst. Metals* 72, 708–713 | `10.2320/jinstmet.72.708` | Atmosphere-specific candidates retained separately; user PDF SHA256 `191ac3e067573b53049ad9c7c6b71e67512970fb4df061973dc4687547bba515` |
| Ag/Ti₃SiC₂ wetting | Zhang et al. (2018), *Int. J. Miner. Metall. Mater.* 25(7), 810–816 | `10.1007/s12613-018-1629-0` | Static, condition-specific 72° datum only; user PDF SHA256 `f0da7547bb2e1dc6b8d107dd8de341245decc3edbed56c273033c22079efebcf` |

P1 numerical baseline before this run: commit `66961178b383fa7c507bea6318f8672e1cb4c08e`. The calculations were conducted in an isolated checkout of `codex/ag70-p1b-p2`; the original `/workspace/-` working tree was not changed.

The raw NIST file includes Au, Ag and Sn groups; the extraction selected the Ag component (`nOrgNum=2`) rather than treating all 19 property points as silver. NIST specifies viscosity in Pa·s, temperature in K, and 95% expanded uncertainty by propagation of evaluated standard uncertainties. See `results/p2_nist/nist_ag_uncertainty_summary.json` and `P2_AG_VISCOSITY_CALIBRATION.md`.
