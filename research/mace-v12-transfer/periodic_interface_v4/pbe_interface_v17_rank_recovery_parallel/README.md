# Exact frame0 controlled relabel

Only V16 train frame0 Ti_gap_2p3,32 atoms; source hash and ordered force-index mapping in input_manifest.json. Input Cartesian coordinates and cell use17 significant digits and were checked exactly against source. No historical IDs or pair metadata are invented. PW-PBE500/Gamma/Fermi0.1, four MPI ranks, one OMP/BLAS thread each, native extrapolated energy plus separately recorded free energy. Source declaration was1x4x2; Gamma convergence along short periodic directions remains unresolved. This new calculation is not recovery of the original method.

run.sh refuses existing results and archives. It runs the generic calculator, verifies compact archive, then publishes completion. Checkpoints remain local under /workspace/mace_v17_rank_recovery_frame0/result; no gpw is uploaded. On any SCF/geometry/hash/MPI mismatch keep all files and stop. Never restart with another directory to evade existing-state guard.
