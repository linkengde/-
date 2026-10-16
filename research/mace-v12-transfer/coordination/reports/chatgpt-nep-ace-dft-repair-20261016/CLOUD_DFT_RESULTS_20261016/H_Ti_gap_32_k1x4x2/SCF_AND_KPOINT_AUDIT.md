# H P0 fixed-geometry single-point audit

Source commit: `ed7d747c968b555144e9b465ff6ac041717dcc6e`. Only the authorized 32-atom Ag10Ti12Si2C8 job was run: PBE / PW 500 eV / Fermi-Dirac 0.10 eV / spin unpolarized / full 3D PBC / gamma-centered 1x4x2 / four MPI ranks. Original coordinates, cell, atomic order and input bytes are preserved. No geometry relaxation, MD, NEP/ACE training or next-grid task was started.

SCF converged: **YES**, successful GPAW energy and force calls and literal log evidence: `["  convergence={'energy': 1e-05, 'density': 1e-05, 'eigenstates': 1e-05},", '  Converge bands: occupied', 'convergence criteria:', 'Converged in 53 steps']`. SCF iterations parsed: `53`. Native energy: `-206.863236742233` eV. Force-consistent free energy: `-207.107968450995` eV. All 32x3 forces are finite; ordered values and 0-based indices are in forces_indexed.csv and DFT_P0_RESULT.json.

Wall time: 914.774 seconds. Peak sampled MPI process-tree sum RSS: 1919635456 bytes; sampled every 2 seconds and may double-count shared pages. This is not an exact per-rank or scheduler peak. Raw output is preserved outside the Git checkout; no large checkpoint is committed. The supplied runner does not write a .gpw checkpoint.

The first launcher attempt exited during CLI parsing before SCF because the `gpaw python --` script-argument separator was omitted. Its error log is retained. The corrected invocation executes the exact pinned runner; no scientific parameters were changed.

New grid centering is explicit gamma=True, confirmed in gpaw.log. The old 1x4x2 centering is **UNVERIFIED**: inspected historical scripts do not establish this particular label's original centering. Therefore no same-method historical k-point convergence claim is made. A single converged SCF is not k-point convergence; 1x6x4 and other structures remain unstarted. C/user review is needed before any next job.

All eleven supplied structures passed Git-blob, atom-count, formula, PBC, finite-coordinate and ASE-read checks; see preflight/. SHA256SUMS.txt covers result and evidence files. Hash PASS verifies archived bytes only, not model quality or resolution of historical label conflicts. No training eligibility or data role has been changed.
