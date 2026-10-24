# B1 phase-specific k-point plan — 2026-10-24

## Outcome

I completed a read-only crosswalk between archived pure-phase k-mesh controls and the frozen B1 candidate cells. Solid Ag and Ti₃SiC₂ have useful numerical precedents, but the B1 displaced cells and alternate supercell shapes are not all the same geometry as those old controls. I prepared same-geometry mesh pairs for the B1 cells and two liquid-Ag snapshots. The plan is disabled: **no DFT was launched, no labels were added, and no dataset or model was changed.**

## Existing evidence

- Ag 32-atom cubic baseline, k6→k8: |Δfree energy| 0.7990 meV/atom; force-vector RMS 1.27×10⁻⁵ eV/Å; maximum atom difference 2.05×10⁻⁵ eV/Å. Archived k4→k6 baseline and a +0.01 Å displacement also pass the 2 meV/atom and 0.01 eV/Å RMS limits. Additional ±2% volume and ±0.5% shear k6→k8 pairs pass.
- Ti₃SiC₂ 48-atom COD-derived baseline, k8→k10 at fixed kz=2 and kz=4: free-energy differences 0.00975 and 0.01669 meV/atom; force-vector RMS 0.00651 and 0.00462 eV/Å. Both archived comparisons pass their stated budgets.
- SHA256 manifests were rechecked for all 13 cited Ag and four Ti₃SiC₂ endpoint archives.

These establish convergence evidence for those archived structures. They do not automatically certify the larger (~0.04 Å RMS) random-displacement B1 structures or the liquid boxes.

## Prepared B1 checks

| Fixed geometry | Proposed mesh pair | Why this pair |
|---|---|---|
| Ag32, seed 9073 displacement | 6×6×6 → 8×8×8 | Same cubic cell family as old 32-atom control; 8.172 Å vs 8.180 Å. Tests the actual B1 displacement amplitude. |
| Ag64, seed 9073 displacement | 5×5×5 → 7×7×7 | Different fcc supercell shape; meshes match the old Ag k6→k8 reciprocal spacing. |
| Bai Ti₃SiC₂ 72-atom, seed 9073 displacement | 5×8×2 → 7×10×2 | Reciprocal-spacing mapping from the old 8×8×2 → 10×10×2 control. |
| Bai Ti₃SiC₂ 96-atom, seed 9073 displacement | 8×8×2 → 10×10×2 | Doubled c length maps old kz=4 to new kz=2; in-plane sampling is directly matched. |
| Liquid Ag 1250 K, ρ ratio 1.00, seed 501 | Γ → 2×2×2; 3×3×3 only if needed | Tests the lower temperature end of the source-supported liquid range. |
| Liquid Ag 1600 K, ρ ratio 1.00, seed 902 | Γ → 2×2×2; 3×3×3 only if needed | Independent EAM seed at the upper intended temperature. |

Exact geometry paths and SHA256 values are in [`mesh_plan.json`](mesh_plan.json). All comparisons must keep the structure fixed and change only k points. Target the B1 force-consistent `free_energy_eV`; retain GPAW's native extrapolated energy separately.

## Pass criteria and gate

Each pair must satisfy |Δfree energy| ≤ 2 meV/atom, force-vector RMS ≤ 0.01 eV/Å, and maximum atom-force-vector difference ≤ 0.03 eV/Å. For liquid Ag, if Γ→2×2×2 misses a limit, run 3×3×3 on that same snapshot and compare again. Do not mix settings across energies and forces.

The proposed GPAW 26.7.0 / ASE 3.29.0 / PBE / PAW 1.2.1 / 500 eV / σ=0.10 eV method remains unfrozen until the registered Window A owner finishes σ=0.05 and σ=0.20 V20 controls. The latest task-card snapshot is 2026-10-24 03:28:13 UTC; it recorded σ=0.05 at iteration 50 and σ=0.20 queued. This is not live telemetry and may be stale. The main branch has no newer committed task-card update; the live owner must publish the current status and archives.

Execution order: owner completes and verifies V20 → freeze B1 method → register and run these mesh pairs → review the archived comparisons → then consider a small B1 labeling pilot. This plan does not authorize parallel DFT or change the owner's jobs.
