# Ag–W QSGS interface reconstruction v4: coarse geometry candidate

This is a fresh QSGS voxel-to-atom construction branch. It keeps the 200 × 200 × 200 Å material, the 200 × 200 × 300 Å `p p f` simulation box, approximate 1:1 Ag/W mass, and the two phases' bicontinuous topology. The cubic 20³ QSGS mask uses 10 Å voxels; the finer original geometry and earlier failed reconstruction records remain untouched.

## Mass-balanced candidate

`coarse_qsgs_massbalanced/` contains 481,638 atoms: 178,100 W and 303,538 Ag. Its Ag/W mass ratio is 1.0000059. The largest W component contains 99.50% of W, the Ag phase is one connected component, and both largest components wind through x/y and span the material in z. The geometry report records the generation parameters and box/layer checks.

This is only an **unrelaxed geometry candidate**. Its initial minimum Ag–W distance is 0.080 Å, so it fails the static MD-start gate. No minimization, 300 K equilibration, TTM, or laser simulation has been run from it. Keep the candidate as evidence of the morphology/mass-balance construction step; do not use it as an MD initial state.

`coarse_qsgs_initial/` preserves the earlier 20³ construction using the direct 35.27% W voxel target. It has a 1.00756 Ag/W mass ratio and is superseded for composition bookkeeping by the mass-balanced candidate, while also remaining unrelaxed.

## Reproduction

From the repository root, write to a new directory (the builder refuses to overwrite its four outputs):

```bash
python research/ag_w_reference_pack/scripts/build_ag_w_qsgs.py \
  --output /tmp/agw-v4-rebuild --grid 20 --seed 20261124 \
  --p-core 0.01925 --p-grow 0.145 --generations 10 \
  --target-w-fraction 0.354625
```

The 10 Å QSGS voxel fraction is an input to shape generation; the mass ratio is then obtained from the atom counts and element masses. The selected target fraction compensates for the distinct Ag and W site densities in this specific lattice projection. Always read the emitted geometry report instead of assuming the voxel fraction equals mass fraction.
