# Ti₃SiC₂ 2×2 k-mesh matrix

**All four points are archived and verified.** B prepared the matrix analysis script; A ran it after the 8×8×2 archive passed verification and independently recomputed the directional contrasts. B's independent review is still requested. No DFT was launched for this analysis.

The matrix uses the same frozen 48-atom structure, source hash, PBE/PW 500 eV settings, and smearing; only k-points differ. Energy gates are 2 meV/atom for native and free energies. The force gate is 0.01 eV/Å for vector RMS.

| Directional change | Native Δ (meV/atom) | Free Δ (meV/atom) | Force-vector RMS (eV/Å) | Max atom ΔF (eV/Å) | Result |
|---|---:|---:|---:|---:|---|
| 6×6×2 → 8×8×2 (in-plane, kz=2) | −0.16335 | −0.11467 | 0.014436 | 0.022319 | Energy pass; force fail |
| 6×6×4 → 8×8×4 (in-plane, kz=4) | −0.10653 | −0.09133 | 0.010850 | 0.016824 | Energy pass; force fail |
| 6×6×2 → 6×6×4 (kz, in-plane 6×6) | +0.31834 | +0.15785 | 0.002884 | 0.004881 | All pass |
| 8×8×2 → 8×8×4 (kz, in-plane 8×8) | +0.37515 | +0.18120 | 0.001142 | 0.001504 | All pass |

The energy difference-in-differences is +0.05682 meV/atom native and +0.02334 meV/atom free. The corresponding force-vector interaction RMS is 0.003613 eV/Å; no separate threshold is defined for this interaction term.

For this frozen structure, the tested kz increment passes at both in-plane meshes, while the 6×6→8×8 in-plane increment exceeds the force budget at both kz values. Thus the force sensitivity is mainly associated with in-plane sampling in this matrix. This does **not** establish convergence beyond 8×8, physical stability, or model accuracy. A denser in-plane point, such as 10×10×2, is the next targeted check before selecting a reference mesh.

The full pairwise results, species/component breakdowns, interactions, archive hashes, and identity checks are in [`report.json`](report.json). Reproduce the matrix with:

```bash
/workspace/.venvs/gpaw-mpi/bin/python analyze_mesh_matrix.py --output /tmp/ti3sic2_mesh_matrix_reproduced.json
```

The report was generated from four verified archives; the large GPAW checkpoints were not part of the archives.
