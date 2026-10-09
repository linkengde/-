# GPAW 26.7 plane-wave slab and dipole semantics

This is a read-only source audit for the registered 26 Å gap-scan settings. No GPAW calculator, DFT run, or MPI job was started, and no gap-scan manifest or input was changed. Machine-readable values and installed-source SHA256s are in `evidence.json`.

## Finding

In GPAW 26.7.0's new plane-wave path, `pbc=(True, True, False)` is passed into the uniform-grid descriptor, but it does not turn the plane-wave basis or its reciprocal-space Poisson kernel into an open-z or Coulomb-truncated calculation. The grid builder sets `zerobc=[False, False, False]` for PW, and `PWDesc` constructs its domain with all three directions periodic. `PWPoissonSolver` applies the periodic reciprocal kernel and sets its G=0 component to zero.

With `poissonsolver={'dipolelayer': 'xy'}`, the PW solver treats the value as enabled and the dipole wrapper infers its axis from the sole nonperiodic grid direction. Thus, for the current z-normal slab (`pbc=[true,true,false]`), the correction is along z. The wrapper requires exactly two periodic axes and adds a sawtooth potential and dipole energy correction. The string `xy` is not parsed into an axis by this new PW wrapper; the legacy `DipoleCorrection` class does parse it. Keep the slab orientation/PBC and Poisson setting explicit in every record.

## Existing matched AgSi control

The verified `AgSi_originalvac_slab_dipole_off_v19` and `AgSi_originalvac_slab_dipole_on_v19` archives have identical geometry, cell, PBC, atom ordering/IDs, source/input hash, and all method fields except `poissonsolver`. Both are 52-atom PBE plane-wave calculations at 500 eV, 6×6×1, 0.10 eV smearing, GPAW 26.7.0 / PAW data 1.2.1, with `pbc=[true,true,false]`; their z cell length is 34.2058765 Å. Both archive verifiers report PASS and both `SHA256SUMS.txt` files validate.

Turning the correction on changes the native energy by +0.0000904919 eV/cell (+0.001740 meV/atom) and the force-consistent free energy by +0.0000881755 eV/cell (+0.001696 meV/atom). The maximum per-atom force-vector difference is 0.003321 eV/Å; RMS vector difference is 0.001211 eV/Å. These are fixed-ion SCF results, not relaxed structures: maximum residual forces are 0.9481 and 0.9487 eV/Å. The small off/on difference at this original-vacuum geometry does not establish the adequacy of a 26 Å vacuum or transfer unchanged to other gaps/geometries.

## Minimum matched controls for the 26 Å scan

Keep one electrostatic target fixed across every separation point. For the registered z-normal slab target, record `pbc=[true,true,false]` and `poissonsolver={'dipolelayer':'xy'}` explicitly, and compare native energy, free energy, full force vectors, and the interface-normal projected forces using the same convention at every point.

If the question is the effect of the dipole correction at the 26 Å geometry, the minimum contrast is an OFF/ON pair with exactly the same coordinates, cell, atom IDs/order, and all other settings. If the scan already contains its ON parent, this requires one additional matched OFF calculation after A approves/enables it. This contrast does not test vacuum-size convergence. If 26 Å is the vacuum width whose adequacy is in question, compare the same slab geometry with the same dipole setting in the 26 Å cell and one larger-z cell; change only the cell height (and preserve the slab's position relative to the vacuum convention). To separate both effects, use the 26 Å ON parent plus one same-cell OFF control and one larger-cell ON control. The registered entries remain disabled pending A's decision.

## Evidence and limits

`evidence.json` provides installed GPAW source paths, line references, SHA256s, archive hashes, and paired numerical results. The installed `gpaw/test/test_dipole.py` uses `mode='lcao'`; it is not evidence that the new PW backend was tested by that test. This report therefore relies on the installed PW source and the two verified PW archives. No new electronic-structure calculation or runtime solver test was performed.
