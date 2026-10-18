# Heating contact anomaly audit

The heating run was stopped after its saved 0.09 ps frame. This is a precautionary stop for a local compressed contact and coincident force spike. It is not evidence that the material melted, and it does not by itself identify whether the cause is the 1 fs step, local interface geometry, or the provisional V17 cross-phase description.

## Evolving minimum contacts

| Heating time (ps) | Hot Ag (K) | Minimum Ag–Ag (Å) | Minimum Ag–TSC (Å) | Minimum Ti–C (Å) | Fmax (eV/Å) |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 259 | 2.648 | 2.446 | 1.869 | 2.317 |
| 0.01 | 788 | 2.592 | 2.331 | 1.914 | 2.369 |
| 0.02 | 1,069 | 2.539 | 2.225 | 1.899 | 2.394 |
| 0.03 | 1,132 | 2.496 | 2.197 | 1.851 | 2.045 |
| 0.04 | 1,271 | 2.483 | 2.166 | 1.870 | 2.240 |
| 0.05 | 1,410 | 2.457 | 2.074 | 1.841 | 2.703 |
| 0.06 | 1,499 | 2.399 | 1.986 | 1.868 | 2.592 |
| 0.07 | 1,587 | 2.381 | 1.961 | 1.845 | 3.181 |
| 0.08 | 1,624 | 2.365 | 1.971 | 1.814 | 3.668 |
| 0.09 | 1,451 | 2.243 | 1.969 | 1.825 | 6.384 |

The most compressed Ag–Ag pair at the final saved frame is IDs 1757 and 1648. Their separation decreases from 2.9065 Å at heating step 0 to 2.2434 Å at step 90. They are Ag atoms in Ag(111) layers 4 and 3, respectively. At step 90 their coordination numbers are both 12 and individual q₆ values are about 0.467, so the pair is still locally crystalline by the chosen order metric even as it becomes unusually compressed. ID 1757 carries the 6.384 eV/Å peak force; its nearest neighbor is ID 1648 at 2.2434 Å.

The minimum cross-phase contact at step 90 is C ID 468–Ag ID 1254 at 1.9689 Å. The minimum Ag–TSC distance was already 2.4463 Å at the start of the 300 K run and decreased during the early baseline. This trend means the interface deserves separate attention. It is not valid to infer a defective potential from one provisional-MACE trajectory alone, but it is also not safe to call the interface dynamically stable.

The local structural indicators do not show melting through the last frame: the hottest sampled frame is step 80 at 1,624 K, with hot-region median q₆ 0.527 versus 0.572 in the 300 K reference, median Ag coordination 12, zero hot Ag atoms below q₆=0.35, and a 0.275 Å hot-region RMS displacement. At step 90 the hot-region q₆ remains about 0.522 and there is no connected low-q₆ cluster. Thus the high force/contact is not a liquid-cluster signature.

The Ti₃SiC₂ backbone remains a single connected component at the 3.25 Å graph cutoff in every saved heating frame. At step 90 its minimum Ti–C distance is 1.825 Å and no atom has crossed the open-Z safety margin.

## Initial interpretation

The observed approach of two Ag atoms and the simultaneous shrinking cross-phase contact may reflect a high-temperature close encounter, finite-step integration, or an unvalidated interface interaction. The initial 0.09 ps record alone did not distinguish those possibilities.

## Saved-state 0.2 fs continuation

The saved 0.09 ps positions and velocities were continued for 100 steps at 0.2 fs (0.02 ps), with the same V17 model and local thermal protocol. Fmax changed from 6.384 eV/Å at the starting frame to 4.037 eV/Å at 0.11 ps. Minimum Ag–Ag distance changed from 2.243 to 2.295 Å, Ag–TSC from 1.969 to 2.013 Å, and Ti–C from 1.825 to 1.832 Å. No NaN, atom escape, or TSC component break occurred. This indicates a local force/contact recovery over that branch, but it is not a paired 1 fs control and cannot prove the original event was a timestep artifact. The TSC graph remained connected; Ag order metrics still did not show a liquid cluster. See `HEATING_DT0P2FS_REPORT.json` and `AG_MELTING_COOLING_ANALYSIS.json`.

The study then turned off the Ag heat source and ran 0.02 ps with the support-side 300 K sink active. Its short, noisy local-temperature record is not sufficient to establish cooling equilibrium. It also does not support resolidification because no melt was identified before heat-off.
