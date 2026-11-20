# v12o 10 fs no-laser contact audit

Source is the full 458,248-atom SNAP v12o cold trajectory at 10 fs (1 fs step). Mobile temperature reached 4,172.2 K without a laser; the run lost no atoms. This is a highly transient failure state, not an equilibrium liquid or a phase-change label.

## Contact counts and extrema

The independent audit found 2,922 Ti-C pairs shorter than 1.8 Å; the minimum Ti-C distance was 1.3587 Å (IDs 114696/108794). Minimum cross-phase distances were Ag-Si 1.3209 Å, Ag-C 1.5136 Å and Ag-Ti 1.9005 Å. Several high-force Ag atoms had strong model forces: e.g. Ag 455657 was 1.4348 Å from Si 339883 with SNAP force norms 79.56 and 77.40 eV/Å; Ag 412351 was 1.9616 Å from Ti 186486 with force norms 91.03 and 26.56 eV/Å.

The four candidate crops use a 5.5 Å sphere and 22 Å cubic periodic cell:

| Label | Pair | Distance | Crop atoms | SNAP force norms on pair |
|---|---|---:|---:|---:|
| ag_si_430045_250182 | Ag-Si | 1.3209 Å | 44 | 37.18, 33.08 eV/Å |
| ag_si_455657_339883 | Ag-Si | 1.4348 Å | 38 | 79.56, 62.77 eV/Å |
| ag_c_454089_336261 | Ag-C | 1.5136 Å | 40 | 27.44, 26.62 eV/Å |
| ag_ti_412351_186486 | Ag-Ti | 1.9616 Å | 37 | 91.03, 26.56 eV/Å |

They are partial-force training candidates: only the two central force rows will use DFT, while the crop boundary will use the established POD v10 force anchor. DFT labels remain pending in `labels/`; only SCF-converged labels may be used.
