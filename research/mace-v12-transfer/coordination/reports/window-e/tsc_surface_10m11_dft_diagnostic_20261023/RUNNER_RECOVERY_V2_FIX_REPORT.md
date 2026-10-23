# SCF result-reader recovery v2

The first recovery rerun reached GPAW convergence and wrote `Converged in 9 steps` to the root log, but non-root MPI ranks read the shared log before rank 0 flushed it. They raised a false `no converged marker` error and MPI aborted before forces/stress were collected. The checkpoint at global step 60 remains the latest valid checkpoint; no label file was created.

The v2 runner preserves the original runner and first recovery runner. It flushes and parses the log on rank 0, then broadcasts the convergence flag, iteration count, and marker count before any rank proceeds to label collection. It does not change geometry, PBE/PAW settings, k-points, cutoff, smearing, mixer, SCF criteria, or checkpoint. A separate v2 execution lock binds the new runner SHA to the unchanged config/input/geometry/PAW hashes.

A second replay from the verified step-60 checkpoint is required because the converged wavefunctions from the aborted process were not checkpointed and the prior runner did not write forces/stress.
