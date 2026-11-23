# Ag–W TTM source preflight

This is a static mapping of the migrated Ag-only candidate source onto the actual Ag–W QSGS atom coordinates. It does not run TTM or MD.

## Source and mask

```text
q_e = f_Ag(x,y,z) * (I0/l_skin) * exp(-r_xy^2/R^2) * exp(-(z_surface-z)/l_skin) * S(t)
```

Candidate values retained from the Ag–Ti₃SiC₂ source note: `I0=5 eV/(ps·Å²)`, `R=30 Å`, `l_skin=20 Å`, center=(100,100) Å, z_surface=200 Å.

Here `r_xy` is the x/y periodic minimum-image distance from the beam center and `d=z_surface−z` is depth into the material. The numerical powers below use `S(t)=1` and the actual Ag atom locations.

The Ag mask is mapped from actual LAMMPS atom type 2 coordinates, with an FCC atomic-volume proxy and cloud-in-cell deposition onto a 10 Å grid. Cell occupancy is locally bounded to [0,1]; it is not replaced by the global 0.6473 Ag volume fraction. The source/exchange convention follows the legacy reference: only Ag atoms in the Newton layer (40–200 Å) receive direct source. W has zero direct source in this candidate.

## On-state power, S(t)=1

| Quantity | Power (eV/ps) |
|---|---:|
| Unmasked continuum over the material box | 14136.456455 |
| Actual Ag atom quadrature, full material | 9787.700707 |
| Actual Ag atom quadrature, Newton layer | 9784.983951 |
| Newton Ag 10 Å grid sum | 9784.983951 |
| W direct source | 0 |

The Newton-layer Ag mask captures 69.2181% of the unmasked spatial source power for this actual morphology. The 10 Å grid sum agrees with atom quadrature by 2.220e-16 relative error.

## Grid sensitivity

| Nominal spacing (Å) | Grid | Raw overfilled Ag cells | Raw max Ag fill | Volume redistributed (Å³) | Newton Ag power (eV/ps) |
|---:|---:|---:|---:|---:|---:|
| 5 | 40×40×40 | 9821 | 1.123583 | 42889.584 | 9784.983951 |
| 10 | 20×20×20 | 724 | 1.008885 | 2288.783 | 9784.983951 |
| 20 | 10×10×10 | 1 | 1.003771 | 30.166 | 9784.983951 |

## Time integral and interpretation

Cumulative input energy is `E_Ag = 9784.983951 eV/ps × integral S(t) dt`. No pulse width or temporal envelope is specified, so the actual cumulative energy is **undetermined**. The coefficient per 1 ps of unit `S` integral is 9784.983951 eV; this is a normalization reference, not a selected pulse.

This assumption heats Ag electrons directly; W can gain lattice energy only through atomistic Ag–W interactions in a later MD model. Electronic transport and electron–lattice coupling in W or across Ag/W are not represented. The Ag/W property closure remains unverified, and the current QSGS geometry failed its static EFS initialization preflight. Do not interpret this report as an Ag–W laser simulation or as a production heating prescription.

Machine-readable values and the actual 10 Å mask/source arrays are in `source_preflight_report.json` and `ag_mask_source_10A.npz`. Source reference: `Migrated spatial source form/parameters from models/ag_ti3sic2_qsgs_1to1_finer/TTM_SOURCE_STATIC_PREFLIGHT.md; legacy project files remain unchanged.`

Reproduce from the repository root with `python research/ag_w_reference_pack/scripts/ag_w_ttm_source_preflight.py`. An optional `--envelope-csv` accepts `time_ps,S` columns to integrate a specified temporal envelope.
