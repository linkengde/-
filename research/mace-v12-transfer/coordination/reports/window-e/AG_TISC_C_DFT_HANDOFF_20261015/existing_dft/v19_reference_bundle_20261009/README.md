# DFT references for local cross-potential checks

16 completed, verified compact cases are included: six signed Ag-C/Ag-Si/Ag-Ti acquisition cases, and ten grid/width controls of ONE frozen52atom AgSi slab. Inputs, results, full logs, original inventories/verifiers, native and free energy, method, cell/PBC, atom IDs and marked contact pairs are retained. No .gpw, launch scripts, live jobs, sealed references, Stage69 trajectories or mixed-potential parameters are included. Original files remain in their canonical repository locations; bundle_manifest.json records provenance.

## Use in another window

```bash
cd /workspace/-
git pull --ff-only origin main
cd research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009
python3 verify_bundle.py
```

Read bundle_manifest.json for the exact paths/settings of every case. The inputs are label-free extxyz. The result extxyz has PW_PBE_energy_eV and the per-atom PW_PBE_forces array. atom_tables/*.csv provides index0, original lammps_id, element, pair marker, positions and force components without requiring ASE. Positions/cell are Angstrom, forces eV/Angstrom, energies eV per whole cell. Cells and PBC are explicit in the manifest.

Calculate the actual EAM+Tersoff+Morse prediction on each EXACT unchanged input, preserving cell/PBC, element order and atom-ID mapping. Use a static run0/single-point calculation; no minimization, rescaling, thermostat, cooling or trajectory extension. Map chemical symbols to the actual potential's LAMMPS types explicitly; atom IDs are not types. Keep existing Stage69 parameters unchanged. Save predicted per-atom forces with index/ID/species and total energy, plus all potential/input hashes and units. Do not substitute a MACE prediction for the existing mixed potential.

Compare all-atom force-vector RMSE as sqrt(mean(sum((F_model-F_DFT)^2,axis=1))), not per-component RMSE. Also report by species/interface region, maximum atom error, and the signed marked-pair projection dot(F_X-F_Ag,unit_Ag_to_X), using the minimum-image pair vector. Report absolute model-minus-reference pair error. Project screening gates are vector0.05eV/A and pair0.10eV/A, each interface separately; energy10meV/atom requires a justified matched reference/energy-zero convention. Raw total energies from differently referenced potentials cannot be compared blindly.

## Limits

The six signed acquisition points reuse previously scored/training parent families. They support local diagnostics, NOT untouched independent validation. The ten slab controls are numerical sensitivity data, NOT ten distinct environments and NOT training/development labels. At sigma0.20 the5/6 pair locally meets numerical budgets, but broader smearing does not establish zero-temperature convergence; width target sensitivity remains material and A's8x8/sigma0.10 follow-up is separate. Do not choose whichever reference gives the smallest model error. Freeze and disclose target settings before scoring.

GPAW native extrapolated energy (force_consistent=False) and free energy are stored separately. Forces are finite-width forces; native energy is not an exactly force-consistent energy/force pair. Preserve this distinction. Do not silently merge Gamma/sigma0.10 acquisition references with periodic sigma0.20 slab targets.

No solidification conclusion, melting temperature, impact validity, or independent cross-potential acceptance is established by this bundle. The receiving window still needs the exact mixed-potential files/type mapping and actual cooling trajectory/restart to perform those checks. Reproduce individual archive checks with the corresponding group's verify_result.py LABEL calculations/LABEL using Python with ASE/NumPy.
