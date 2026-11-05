# P1H arc-pressure module design

P1H is an isolated exploratory branch based on P1G commit
`c262e2ea37f5b631da98a64b44604421150f36f4`. The P1G source, thermal
solver, geometry, and completed A/B results remain unchanged. P1H adds a
separate solver copy and pressure module; it is not a production-solver edit.

## Force definition and sign

The tracked VOF fraction is `f=1` in the composite sample and `f=0` in gas.
At the exposed upper surface, `grad(f).z < 0`. Basilisk's `iforce.h` converts
the interfacial potential to acceleration as `alpha * phi * grad(f)`, so a
positive pressure potential produces a downward force. P1H uses

```text
phi_total = sigma_effective * curvature + p_arc_effective
```

and leaves the existing tangential Marangoni force path unchanged. Pressure is
added only in VOF interfacial cells. Its mask is the product of the conditional
Ag coverage (`ag_marker / f`), Ag liquid fraction, and upper-facing normal
weight `clamp(-n_z, 0, 1)`. This avoids treating the partial VOF cell fraction
itself as a partially covered surface.

The first pressure profile is spatially uniform over the eligible upper
surface. An optional periodic Gaussian uses the minimum-image x/y distance, so
it is continuous across the periodic seams. Its radius is independent of the
thermal source radius. A Gaussian profile is only a numerical sensitivity
case unless an independent physical source for that spatial scale is supplied.

## Time schedule and units

Pressure is entered in Pa; center and radius environment values use nm; start,
ramp, and end times use ps. The pressure ramp uses a smoothstep envelope. On
the heating stop event, P1H starts matched heat-source and pressure off-ramps at
the same simulation time. If pressure is enabled, their configured ramp
durations must match. The heat input and energy ledger use the ramped source
scale, and solidification is not declared until the ramp has ended.

Restart sidecars save the pressure profile and ramp state. A restart with a
different pressure or heat-ramp configuration is rejected.

## Pressure scale

Sensitivity cases are parameterized by `Pi_p = p0 * Lc / gamma_ref`, with
`Lc=5 nm` for the local P1G melt-pool scale. `Pi_p` values are mechanism-test
levels only; no macro-scale arc-pressure value is imported or identified as a
measured nanoscale pressure. A uniform load on a flat interface is expected to
be balanced by pressure and should not be used to force a crater.

## Interface and boundary limits

The P1G VOF tracks the composite-to-gas outer surface. It does not resolve
gas-liquid interfaces in internal pores. P1H therefore applies pressure only
to the exposed outer surface and makes no claim about pressure inside the
Ti3SiC2 pore network. The inherited top boundary is closed. Gas-cap height and
open-versus-closed top-boundary sensitivity remain a prerequisite before
formal C/D heating runs.

The pressure module does not introduce evaporation, recoil, ejection, or Ag
mass loss. It does not establish a physically calibrated arc-pressure field,
contact-line law, or quantitative nanoscale flow prediction.
