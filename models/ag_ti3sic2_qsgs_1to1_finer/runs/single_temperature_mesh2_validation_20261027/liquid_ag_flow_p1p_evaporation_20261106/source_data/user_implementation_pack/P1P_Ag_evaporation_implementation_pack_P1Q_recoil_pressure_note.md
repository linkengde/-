# P1Q recoil-pressure reference note — DO NOT IMPLEMENT YET

## Why P1Q remains frozen
P1P has not yet passed multi-step real-material mass/energy conservation. Recoil pressure must not be added before the evaporation flux and interface temperature are trustworthy.

## Relevant reference
Wang, Zhang & Yan (2020), Phys. Rev. Applied 14, 064039, DOI 10.1103/PhysRevApplied.14.064039.
Full text: https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevApplied.14.064039/fulltext

The paper derives mass, momentum and energy fluxes across a Knudsen layer and couples them to the external vapor-flow structure. It explicitly distinguishes emitting molecules and back-scattered/condensing molecules. Its recoil pressure is a momentum-flux result (their Eq. 12), not simply a generic pressure multiplier.

The paper also presents the Anisimov approximation:
`P_rec = 0.54*P_b*exp[L_v*(T-T_b)/(R_u*T*T_b)]`
where the factor 0.54 follows from a sonic-vapor-velocity assumption. The paper discusses why this simplified model can be unsuitable in near-vacuum conditions. Do not apply 0.54*p_sat as a universal Ag/nanoscale law.

## Critical correction to the earlier simple formula
The expression `p_rec = j_mass*sqrt(pi*R_specific*T/2)` is at most a one-sided Maxwellian effusion momentum estimate under a specified half-range distribution. It is not a general recoil-pressure closure when back-scattered vapor, condensation, ambient-gas collisions, or a Knudsen layer are present. Do not implement it as a production formula without deriving the incoming and outgoing momentum moments consistently.

## Nanoscale decision
Wang et al. (2020) analyzes macroscopic keyhole scales (their discussion uses roughly 150–500 micrometres) and a continuum vapor-flow region outside the Knudsen layer. That gas-flow model cannot be transferred directly to a 3–10 nm Ag pool. At the current nanoscale, first estimate the Knudsen number using the actual gas mixture and local partial pressure. If Kn is large, a free-molecular/kinetic treatment or atomistic MD boundary calibration is more appropriate than continuum gas Navier–Stokes.

Holyst & Litniewski (2009), DOI 10.1063/1.3077059, is a useful MD reference: it studies liquid-film evaporation into vacuum and compares momentum flux with Hertz-Knudsen descriptions. It is not an Ag/Ti3SiC2 parameter source.

## Required P1Q decision after P1P
Choose one consistent closure:
1. A kinetic half-space momentum-flux model with explicit outgoing and incoming distributions, or
2. A validated Knudsen-layer model with specified gas composition, background pressure and collision properties, or
3. MD-calibrated interfacial momentum flux for the relevant Ag surface state.

Do not combine these closures or add a separate empirical 0.54 factor on top of a momentum-flux model; that would double-count recoil.
