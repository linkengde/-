# Ti3SiC2 48-atom k-point comparison

This record compares the converged fixed-geometry GPAW single points at standard Monkhorst-Pack 10x10x4 and 12x12x4. It does not modify the original 10x10x4 result or claim convergence for other structures.

- Input: `TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz`
- Input commit: `cc65a3e6b75be28c06fda8608b4064662d714117`
- Input SHA-256: `da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b`
- Prior 10x10x4 result: `TSC48_FINAL_RESULT_20261021/RESULT.json` from base commit `85c3f25b185c805f8baae8ec75f65e69bc7b79aa`
- 12x12x4 SCF: converged in 43 steps; GPAW log is archived here.
- Comparison: native energy difference 0.041859 meV/atom; force-vector RMS 0.00162533 eV/A. Both pass the predeclared strict gates of 2 meV/atom and 0.01 eV/A for this fixed geometry.
- Data role remains `HOLD` pending C review.

See `manifest.json` and `SHA256SUMS` for hashes and environment records. The raw 12x12x4 result JSON is preserved unchanged; SCF convergence is evidenced by the GPAW log and comparison manifest.
