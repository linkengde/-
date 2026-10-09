# Disabled pure-phase and short-range control pack

Prepared fourteen deterministic, label-free static controls. All launch_enabled=false, owner_instance=null; cutoff/mesh/width/energy target remain unset pending A's reference decision and phase-specific convergence. No calculator, inference, training, trajectory, thermostat, relaxation, collision input or production edit was executed.

## Sources and exact cases

Ag: conventional fcc a=4.09 A from ASE3.29 reference_states, repeated2x2x2 (32atoms). This is an explicitly declared modeling starting lattice, not an independently verified experimental lattice or an optimized value. Source/builder hashes are recorded.

Ti3SiC2: uploaded COD9009647 CIF, SHA2566abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3; symmetry-expanded cell repeated2x2x1 (48atoms, C16Si8Ti24). Internal CIF provenance is evidenced; remote checksum remains unverified. No altered occupancy or different independent parent is claimed.

Each phase has: baseline; volume ratios0.98/1.02 with isotropic length scaling ratio^(1/3); atom0 x displacement−/+0.01 A with cell unchanged; simple shear Fxy=−/+0.005 with cell and ions affinely transformed, determinant1. No combined perturbations, optimization or temperature assumptions. A baseline need not be an energy minimum; use symmetric controls to measure the response before any stability conclusion.

New local_new_id=1..N records are pack-local, explicitly not historical LAMMPS IDs. IDs/order are preserved within paired controls. No marked interface pair is invented.

ASE neighbor_list with8 A search includes periodic image neighbors. Species-dependent conservative overlap floor0.55*(covalent radius_i+radius_j) passes every case; pair minima are recorded. This is geometry triage, not chemical equilibrium or force validation. Every extxyz is read back with position/cell tolerance1e-8 A, identical species/local IDs and no energy labels/calculator. Input/source/tool hashes and complete checksum inventory are provided.

## Short-range plan, no generated collisions

Ag-Ag, Ag-Ti, Ag-Si, Ag-C proposed separation grid:0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0 A; central-difference h0.005 A, with h/2 checks. Every distance remains review-required, especially sub-Angstrom points where frozen-core PAW setups may be invalid. No grid point is an approved job or production input.

First compare neutral isolated pairs with reviewed spin/atomic-reference and vacuum>=10 A per side, including boundary/cutoff convergence. Then use a separate reviewed TRAIN-only pure/surface parent for embedded scans, fixed surroundings and species-aware minimum-distance guards for all pairs. Do not use any sealed or development/test parent as an acquisition source. Small-cell pair curves do not prove many-body interface/impact transfer.

Let rhat point i→j. Separating force S=(Fj-Fi)·rhat/2 is positive for repulsion. For symmetric pair displacement compare S with −[E(r+h)−E(r−h)]/(2h), using the explicitly force-consistent energy convention. Report native extrapolated/free energies separately. Check finite values, adjacent energy/force continuity and converged far-separation asymptote; do not assign arbitrary smoothing tolerances before reference-error review.

Multi-element MACE and a current mixed potential are distinct implementations. Before comparing require the mixed-potential source/configuration, hybrid pair mapping, element order, units/cutoffs and any ZBL splice. Those assets are not established here. ZBL assumes screened nuclear Coulomb repulsion and its atomic-number screening model: useful as a reviewed high-energy asymptotic comparison, not an equilibrium chemistry label. Explicit PAW validity and matching/splice policy are needed; do not silently assert ZBL is the full DFT reference.

## Generic pure-phase verifier changes for A

Use exact phase formula, ordered symbols/positions/cell/PBC/local IDs and input hashes. Neither a marked interface pair nor all four elements is required: pure Ag has one species, Ti3SiC2 has three. Verify finite energies/Nx3 forces, actual SCF completion/MPI ranks and exact code/PAW/cutoff/mesh/width/magnetism/Poisson/energy convention. Preserve native/free energy separately, output identity and immutable checksummed inventory. Retain owner/no-existing-run/disk guards. Stability response/derivative acceptance requires quantified reference precision and the chosen target, not merely SCF success.

Needed before launch: A fixed-target decision after6/8; reviewed phase-specific convergence; explicit jobs/owners; approved pure-phase verifier; mixed-potential configuration and short-range physical-reference review where applicable. Stage69 morphology remains a separate missing asset. These pure-phase controls do not supply independent Ag-interface validation families or authorize MD/impact.

Reproduce from /workspace/-:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_pure_phase_and_short_range_entry_design/prepare.py
```

Outputs: inputs/, input_manifest.json, collision_scan_plan.json, report.json, SHA256SUMS.txt. Generation/screening passed; physical energy/forces/stability remain uncomputed.
