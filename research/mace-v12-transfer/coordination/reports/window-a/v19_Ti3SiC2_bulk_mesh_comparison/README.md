# Ti3SiC2 bulk mesh comparison

The verified frozen 48-atom baseline at 6×6×2 and 8×8×4 has matching geometry, IDs, PBE/PW500 settings, smearing, PAW data, and convergence thresholds. Native and free-energy differences are within 2 meV/atom, but the all-atom force-vector RMS is 0.01359 eV/Å, above the 0.01 eV/Å budget. This comparison changes both in-plane and c-axis k sampling, so it cannot identify which direction dominates.

Next, the same baseline is registered at 8×8×2 to compare with 6×6×2 while holding the c mesh fixed. This isolates the in-plane increase; it remains a numerical control, not an independent validation or thermal-stability test. The interface gap scan stays disabled until the reference method is reviewed.

Machine-readable metrics, thresholds, geometry checks, SCF metadata, and archive hashes are in `k6_k8x8x4_comparison.json`.
