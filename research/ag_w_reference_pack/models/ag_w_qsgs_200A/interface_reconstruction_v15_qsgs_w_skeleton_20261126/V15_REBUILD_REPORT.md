# V15 QSGS rebuild report

## Request implemented

Rebuilt the morphology as a random four-parameter QSGS W backbone, then assigned Ag to the complementary mask. The target appearance is the supplied Ag–Ti₃SiC₂ reference: black scaffold, pale-blue matrix, full-height material with vacuum above. V15 does not use the V11 Gyroid.

## Geometry result

- Seed: `20261141`
- Four QSGS controls: `p_core=0.005`, `p_grow=0.07`, `generations=12`, target W voxel fraction `0.35411`
- Mask grid: `120³`; voxel size `1.667 Å`
- Material: `200 × 200 × 200 Å`; box: `200 × 200 × 300 Å`; `p p f`
- Actual W/Ag mask volume fractions: `35.410995% / 64.589005%`
- Atoms: W `178,527`, Ag `304,151`, total `482,678`
- Ag/W mass ratio: `0.9996288` (relative error `0.0371%`)
- W largest-component share `98.20%`, Ag `99.91%`; each phase spans z and has x/y winding rank 2.

## Gate result

**PASS: QSGS complementary-mask and bicontinuous geometry checks.**

**FAIL: atomistic MD-start gate.** Minimum Ag–W contact is `0.08 Å`, caused by selecting atoms from independent ideal FCC and BCC lattices on opposite sides of the voxel mask. No static relaxation, MD, TTM, heating, melting, or cooling run is claimed for V15. The preview is therefore a geometry visualization, not evidence of a valid simulation initial state.

The next required step is interface-aware atom placement/relaxation that preserves the QSGS masks and phase connectivity, followed by force, short NVT/NVE, and only then TTM checks. Previous V5–V14 data remain untouched.
