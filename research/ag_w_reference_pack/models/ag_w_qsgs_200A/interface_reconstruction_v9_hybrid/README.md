# Ag–W V9 hybrid reconstruction audit

V9 combines the complete V7-derived atom structure with V8's grid-phase search and a single direct-capacity-flow routine shared by the production mapper and diagnostics. The box, atom coordinates, species counts, and atomic-volume constants were not changed to obtain a mapping pass.

**Current decision: FAIL as a TTM-ready candidate.** The full atom structure and cold NVT/NVE runs are preserved, but no tested 5/10/20 Å grid phase can place all cell volume within the 20 Å direct-transfer limit. The final cold state also retains a large Ag/W temperature split and continuing potential-energy relaxation. No TTM step or laser source was applied.

See [AGW_V9_HYBRID_AUDIT.md](AGW_V9_HYBRID_AUDIT.md) for the stage-by-stage findings, limitations, and reproducibility commands. Input structures and V5–V8 records remain unmodified.
