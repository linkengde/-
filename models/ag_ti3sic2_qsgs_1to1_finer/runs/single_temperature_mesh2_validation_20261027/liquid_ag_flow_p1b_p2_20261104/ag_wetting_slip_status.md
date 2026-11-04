# Ag/Ti₃SiC₂ wetting and momentum-transfer status

## Static wetting datum

Zhang et al. (2018), “Microstructure and properties of Ag–Ti₃SiC₂ contact materials prepared by pressureless sintering,” *International Journal of Minerals, Metallurgy and Materials* 25(7) (2018), report a sessile-drop measurement: approximately 60 mg of Ag wire was placed on a polished Ti₃SiC₂ substrate, heated to 970 °C under vacuum, and the contact angle was measured as 72°. The user-provided PDF inspected for this report has SHA256 `f0da7547bb2e1dc6b8d107dd8de341245decc3edbed56c273033c22079efebcf`.

This is a static contact-angle reference for that specimen and preparation condition. It is not a universal Ag/Ti₃SiC₂ boundary condition. The paper does not provide enough information in the inspected method/results passage to reconstruct the exact surface termination, oxygen partial pressure, roughness distribution, contact-line speed, or a quantified angle uncertainty. The P0 72° test established that the numerical wall-angle implementation can represent a prescribed angle; it did not independently calibrate Ag/Ti₃SiC₂ interface physics.

## Dynamic interface quantities

| Quantity | Status | Decision |
|---|---|---|
| Wall friction coefficient λ | Not found in the inspected experimental sources | Leave unresolved; do not assign a guessed value |
| Slip length b | No directly applicable Ag/Ti₃SiC₂ measurement or qualified MD result found | Leave unresolved; do not assume a calibrated Navier-slip law |
| Temperature-dependent advancing/receding angle | Not reported by the cited static test | Do not infer from the single 72° datum |
| Interfacial momentum exchange | Not calibrated | No physical liquid-flow boundary condition yet |

An older pure-Ag EAM trajectory in the repository produced a density about 7.6% below the Abdullaev liquid-density correlation near 1600 K. That potential has not been qualified for Ag/Ti₃SiC₂ interfacial chemistry, slip, or friction. It cannot be used here as a substitute calibration.

## Consequence

The 72° value may be used as a condition-specific static wetting comparison in a future controlled case. It cannot set wall friction or slip, and it does not remove the need to characterize actual surface condition. No flow calculation should claim quantitative Ag/Ti₃SiC₂ wall transport from the data currently available.
