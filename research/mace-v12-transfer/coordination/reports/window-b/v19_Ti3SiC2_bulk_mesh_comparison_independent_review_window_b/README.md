# Independent Ti₃SiC₂ bulk k-mesh comparison audit

This read-only audit independently recomputes A's published comparison of `Ti3SiC2_baseline_k6x6x2_v19` and `Ti3SiC2_baseline_k8x8x4_v19` from their compact archives. The comparison report and both archive SHA256 manifests pass. The two archives have 48 atoms with exactly identical positions, cell, PBC and `local_new_id`; their input/source hashes match, and their method dictionaries are identical except for k-points.

Using k8x8x4 minus k6x6x2, the native energy difference is +0.211807 meV/atom and free-energy difference is +0.066525 meV/atom. Both are within the 2 meV/atom limits. Force-vector RMS is 0.013593 eV/Å and maximum atom-vector difference is 0.021326 eV/Å, so the 0.01 eV/Å RMS force budget fails. The largest species RMS is Ti at 0.017386 eV/Å; the z-component RMS (0.010619 eV/Å) exceeds x (0.007349) and y (0.004243). All recomputed metrics agree with A's JSON report within 2e-12.

The 6×6×2 versus 8×8×4 pair changes both in-plane and z sampling, so it cannot attribute the force difference to either direction. A's registered 8×8×2 input has the same source/input hash and method except k-points as 6×6×2, while keeping kz=2. Once its calculation is archived and verified, 6×6×2 versus 8×8×2 will isolate the in-plane increment at kz=2; 8×8×2 versus 8×8×4 will isolate the z increment at in-plane 8×8. At the audit snapshot the 8×8×2 job was registered and running at SCF iteration 0, with no result available for numerical review. This three-point path does not measure a full mesh interaction term (which would also require 6×6×4) and does not establish broad convergence, relaxed stability, or interface/model accuracy.

`audit_ti3sic2_meshes.py` verifies the published report and archive manifests, identity/method invariants, recomputes energy/force/species metrics, checks the arithmetic against A's report, and records the 8×8×2 registration state. It writes `audit.json`, `audit.csv`, `verification.json`, and `SHA256SUMS.txt`. To reproduce into an empty directory from the repository root:

```bash
python3 research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_bulk_mesh_comparison_independent_review_window_b/audit_ti3sic2_meshes.py --output-dir /tmp/ti3sic2-mesh-audit-reproduction
```

The audit uses archived outputs only; it does not run DFT/MPI or modify A-owned files.
