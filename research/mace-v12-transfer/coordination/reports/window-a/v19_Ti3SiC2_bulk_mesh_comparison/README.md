# Ti3SiC2 bulk k-mesh checks

This report contains separate fixed-geometry comparisons. All controls use the same 48-atom CIF-derived Ti3SiC2 geometry and PBE/PW500 settings; each comparison changes only the stated k-point mesh components. These are numerical controls, not interface-force validation, independent validation, or a thermal-stability test.

## Earlier controls

The original 6×6×2 → 8×8×4 diagonal comparison changed both in-plane and c-axis sampling. Its energy differences met the 2 meV/atom budget, while all-atom force-vector RMS was 0.01359 eV/Å and exceeded the 0.01 eV/Å budget. It could not identify which mesh direction caused the force difference.

## Matched kz=2 refinement: 8×8×2 → 10×10×2

The 10×10×2 calculation converged in 37 SCF iterations on four MPI ranks (3.30 h). Its compact archive verification is **PASS**. The comparison script independently checks both archive SHA manifests, input/method/geometry/atom-ID matching, then recomputes the energy and atom force-vector differences.

- Native energy difference: **-0.009369 meV/atom** (absolute-value limit 2).
- Free-energy difference: **-0.009753 meV/atom** (absolute-value limit 2).
- All-atom force-vector RMS: **0.006510 eV/Å** (limit 0.01).
- Maximum per-atom force-vector difference: **0.009628 eV/Å**.
- Gate result: **PASS**.

The high-density in-plane comparison passes for this frozen geometry at kz=2. It does not establish convergence for perturbed structures, other kz values, relaxed cells, thermal configurations, or interfaces. B's matched kz=4 8×8×4 → 10×10×4 calculation remains the complementary comparison; wait for its verified archive before drawing a two-kz conclusion.

Machine-readable results are in `k8x8x2_k10x10x2_comparison.json`; reproduce them with `python3 compare_mesh_pair.py --left Ti3SiC2_baseline_k8x8x2_v19 --right Ti3SiC2_baseline_k10x10x2_v19 --output /tmp/ti3sic2_mesh_comparison.json` from the repository root. The 10×10×2 archive is under `periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/`.
