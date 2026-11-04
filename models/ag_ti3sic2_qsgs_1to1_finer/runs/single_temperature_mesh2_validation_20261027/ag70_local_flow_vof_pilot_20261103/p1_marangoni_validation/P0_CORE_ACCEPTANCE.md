# P0 core acceptance

**Decision: accept the automatic-curvature VOF path for the standalone static-drop verification scope only.** Keep the exact-curvature/wavelet-AMR issue as a separate open defect. P1 was restricted to automatic curvature and did not call the exact-curvature branch.

## Accepted evidence

- Automatic-curvature spherical drops completed at R/Δx = 16, 32, and 64. Laplace pressure error was 0.503%, 0.126%, and 0.031%, respectively.
- The 72° sessile-drop fit was within 0.016° at the coarsest grid and improved with refinement.
- The maximum mass drift in the accepted static-drop matrix was about 0.000281%.
- At R64, the automatic-curvature case completed its 91-step dynamic check. The measured maximum parasitic speed was 8.62×10⁻⁵ in the dimensionless P0 units, with liquid-volume-weighted RMS 8.11×10⁻⁶.

## Known limitation kept separate

At R64, prescribed exact curvature passes on the fixed geometric interface mesh but fails when combined with the wavelet-AMR topology during pressure projection. Automatic curvature on the same AMR topology remains stable. The full uniform 512³ exact-curvature check was stopped by the 30-GiB resource guard before initialization; it gives no numerical conclusion. No fix was attempted, and this auxiliary path remains prohibited for P1.

## Scope

These results validate numerical components in dimensionless static-drop tests. They do not validate a heated Ag melt pool, integrated heat-flow coupling, a Ti₃SiC₂ skeleton, or an Ag-specific flow prediction. Detailed source logs and tables remain in the neighboring `p0_static_drop_3d` directory; the P0 diagnostic commit is `968e2ea955d9f7b29aed51f8a78e28865e589bce` on `codex/ag70-vof-p0`.
