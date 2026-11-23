# Ag-W composition target: 1:1 mass ratio

The target for the bicontinuous Ag-W model is **equal mass of Ag and W**. It is not an equiatomic alloy.

Using standard atomic weights (M_{Ag}=107.8682\,g\,mol^{-1}) and (M_W=183.84\,g\,mol^{-1}), equal masses imply

\[
x_{Ag}=\frac{M_W}{M_{Ag}+M_W}=0.6302,\qquad
x_W=\frac{M_{Ag}}{M_{Ag}+M_W}=0.3698.
\]

Thus the target atom fractions are approximately **63.02 at.% Ag and 36.98 at.% W**. For a finite LAMMPS model, choose integer counts close to these fractions, then calculate the realized mass ratio from the actual counts and atomic weights:

\[
\frac{m_{Ag}}{m_W}=\frac{N_{Ag}M_{Ag}}{N_WM_W}.
\]

For example, 6302 Ag and 3698 W atoms give approximately 1:1 mass (with the atomic weights above). A 50:50 atom-count AgW B2 cell is only a potential benchmark and is not the target composite composition.

For an initial two-phase morphology estimate, using room-temperature bulk densities \(\rho_{Ag}=10.49\,g\,cm^{-3}\) and \(\rho_W=19.25\,g\,cm^{-3}\), equal masses correspond to approximately **64.73 vol.% Ag and 35.27 vol.% W**:

\[
\phi_{Ag}=\frac{1/\rho_{Ag}}{1/\rho_{Ag}+1/\rho_W}=0.6473,\qquad
\phi_W=1-\phi_{Ag}=0.3527.
\]

These are ideal dense-phase estimates for initializing a QSGS volume split. Porosity, interface excess volume, lattice mismatch, and any later relaxation change the realized volumes. Record atom counts, masses, and the final phase volumes for every generated structure; do not infer composition from the B2 validation cell.
