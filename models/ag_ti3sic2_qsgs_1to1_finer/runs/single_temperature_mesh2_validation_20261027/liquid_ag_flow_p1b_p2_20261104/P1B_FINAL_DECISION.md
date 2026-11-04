# P1b final decision

- **Official Basilisk CLSVOF benchmark:** passed the tested reference checks. The measured migration-speed error decreased from 3.467% to 0.340% over 8/16/32 cells per radius; maximum logged mass drift was 0.0594%, below the 0.1% limit.
- **Project automatic-curvature standard-VOF comparison:** not accepted. Speed error was 4.062%, 9.437%, and 9.998%; the two finer grids did not converge toward the theoretical speed. Mass conservation alone does not resolve this discrepancy.
- **Research consequence:** do not use the project automatic-curvature route for quantitative dynamic Ag flow yet. Keep the P1 flat-interface and static-force results as their own passed tests. Diagnose the variable surface-stress discretization and repeat the dynamic benchmark before P3 or real Ag flow.
- **Scope:** the reference test prescribes a linear surface-tension gradient and does not solve heat transfer. It is not an Ag or melt-pool prediction.
