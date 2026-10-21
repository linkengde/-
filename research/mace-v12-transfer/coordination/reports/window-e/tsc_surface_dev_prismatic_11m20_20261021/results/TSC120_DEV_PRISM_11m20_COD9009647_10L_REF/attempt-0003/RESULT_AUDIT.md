# Frozen-geometry reference single-point result audit

Status: **SUCCESS**. GPAW reported SCF convergence after 44 iterations. This is a frozen-geometry reference label; it is not a relaxed or equilibrium surface structure.

- Structure: `TSC120_DEV_PRISM_11m20_COD9009647_10L_REF`; 120 atoms; composition `{'C': 40, 'Si': 20, 'Ti': 60}`.
- Method: PBE/PAW, PW 500 eV, Fermi-Dirac 0.10 eV, standard Monkhorst-Pack [7, 2, 1], GPAW 26.7.0, ASE 3.29.0, 4 MPI ranks, OMP=1.
- Total energy: -939.221706735912 eV.
- Force-consistent free energy: -940.057949456019 eV.
- Maximum force: 1.71277119 eV/Angstrom. The large residual force is recorded as a feature of this unrelaxed geometry; no relaxation was performed.
- Stress (GPa): `[[0.813566079002882, -0.00021665631306847028, -2.5563382865262205e-05], [-0.00021665631306847028, -0.3373252921856701, 0.00010816908584762225], [-2.5563382865262205e-05, 0.00010816908584762225, -0.040511717693432474]]`.
- Exact input SHA-256: `fd5ff5c1cae1490bcbde9ffd08f329c7ec8cac687555f29568d725f4230cc5c3`; geometry SHA-256: `e876359dc190632d02ca2630bbd145afb4b866ddd19de65c8e1019c2e6f6c1a4`; configuration SHA-256: `23630f8ac5bd99b4b526c67b2cac40c5d7161aee5b682059623da5325d24515a`.
- Input/output check: atom count, element order, coordinates, cell and PBC match exactly; force array is 120x3 and finite; energy and stress are finite.
- Runner SHA-256: `254c3b416a2a9b35f4bec48265703eb1fe021fa8ed20ef1811329320a52b9bce`.
- Checkpoint snapshot at logical iteration 40 SHA-256 `406db04f6b32eb9deb9ddd6ff873bec3c809babd5a0a0fd908aaa935fe082fd9`; it is not the final converged state and was not used as a label. Checkpoint binary is intentionally excluded from this result archive.
- Peak sampled memory: 18063835136 bytes. Free disk at end: 14118141952 bytes.
- Role remains `DEV_GEOMETRY_FAMILY_REFERENCE_CANDIDATE`; no TRAIN/DEV/TEST assignment was made.
