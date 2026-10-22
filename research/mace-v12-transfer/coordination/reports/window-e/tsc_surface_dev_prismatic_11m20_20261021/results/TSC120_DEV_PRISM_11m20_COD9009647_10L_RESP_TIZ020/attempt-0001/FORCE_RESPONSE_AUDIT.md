# Ti-terminal paired-displacement DFT response audit

Status: **SUCCESS**. Both the frozen reference and the frozen response geometry converged (reference 44 SCF steps; response 49 steps) with the same locked GPAW/PBE/PAW configuration (GPAW 26.7.0, ASE 3.29.0, PW 500 eV, Fermi-Dirac 0.10 eV, standard Monkhorst-Pack 7×2×1, Davidson niter=2, MSR1 beta=0.05/nmaxold=10/weight=20, 4 MPI ranks, OMP=1). No relaxation or cell change was performed.

The response input differs from its reference only at Ti source atom IDs 4 and 113, with z displacements -0.020 and +0.020 Å respectively. Atom order/IDs, element list, lattice and PBC match exactly; the response output geometry exactly matches its SHA-bound input. All force and stress values are finite.

Total energy changes from -939.221706735912 eV to -939.190207387186 eV (ΔE = +0.031499348726 eV). Force-consistent free energy changes by +0.030614284240 eV. Stress difference (GPa) is `[[0.014873825762422355, -0.000433765063510476, 0.0010973260136859406], [-0.000433765063510476, 0.014559860725619767, 0.008549789272525227], [0.0010973260136859406, 0.008549789272525227, 0.03590525752954719]]`.

The force difference is `F(response) - F(reference)` mapped by `source_atom_id`. Across all atoms, component RMSE is 0.01786344 eV/Å, atom-vector RMS is 0.03094039 eV/Å, P95 atom-vector magnitude is 0.03386706 eV/Å, and the maximum is 0.21162530 eV/Å at atom ID 4.

Element-group vector RMS of the force changes: Ti 0.03984448, Si 0.01068541, C 0.02081972 eV/Å. These quantify the response to the one paired antisymmetric displacement mode; they do not independently identify derivatives for the two moved atoms and are not a fit or held-out score.

This geometry shares the COD 9009647 source bulk with existing data. It remains a geometry-family candidate only; no data role was changed and it does not resolve the independent-parent-family DEV/TEST gap.

Input SHA-256: `6f91573bcb4a0feb7e9f05329f6eca422d3e7aac723763e961b52c7b83ddf9b6`; geometry SHA-256: `817e65667a2ab51194cc69c9800b1efc4dcadfc49e7a3dbf3867a8577df54fb3`; config SHA-256: `23630f8ac5bd99b4b526c67b2cac40c5d7161aee5b682059623da5325d24515a`; runner SHA-256: `254c3b416a2a9b35f4bec48265703eb1fe021fa8ed20ef1811329320a52b9bce`.
