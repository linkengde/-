# Ti3SiC2 (10-11) fixed-geometry DFT diagnostic

This package runs the two frozen 69-atom geometries in the separately archived
`tsc_surface_geometry_dev_candidate_10m11_20261023` package. Both come from COD
9009647. They are one new facet-family diagnostic, not an independent-source
DEV or TEST family. This package assigns no TRAIN/DEV/TEST role and does not
access evaluator-held TEST data.

The reference and paired terminal-Ti response use the PBE/PW500, standard
Monkhorst-Pack 6x6x1, Fermi-Dirac 0.10 eV, non-spin-polarized GPAW settings
specified in `DFT_DIAGNOSTIC_CONFIG.json`. Coordinates, cell, and PBC remain
fixed. The runner writes energy, force-consistent free energy, atomic forces,
and stress only after a converged SCF. It saves versioned, SHA-bound GPW
checkpoints every ten SCF steps. A new checkpoint becomes current only after
the file and metadata are written; the previous valid checkpoint is removed
only after that commit. Checkpoints are removed after successful result
validation and hashing.

Run one structure at a time:

```bash
./run_one.sh TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF
./run_one.sh TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_RESP_TIZ020
```

Resume a stopped run only after the committed latest pointer, checkpoint SHA,
input SHA, geometry SHA, and config SHA validate:

```bash
./run_one.sh TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF --resume
```

The initial run is the reference only; the response task must wait until the
reference has been archived and verified. No geometry optimization,
trajectory, model scoring, or potential fitting is included.
