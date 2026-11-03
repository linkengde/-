#!/usr/bin/env python3
"""Generate a dimensionless static-drop interFoam verification case."""

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "static_drop_case"
NX = NY = 80
R = 0.4
DOMAIN_HALF_WIDTH = 1.0
DZ = 2.0 / NX
NU = 0.01
END_TIME = 0.2
MAX_DT = 1e-4


def header(cls, obj, body, location="constant"):
    return f'''FoamFile
{{
    version 2.0;
    format ascii;
    class {cls};
    location "{location}";
    object {obj};
}}

{body}
'''


def write(rel, text):
    path = CASE / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="ascii")


def main():
    global CASE, NX, NY, DZ, NU, END_TIME, MAX_DT
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=80)
    parser.add_argument("--case-name", default=None)
    parser.add_argument("--nu", type=float, default=0.01)
    parser.add_argument("--end-time", type=float, default=0.2)
    parser.add_argument("--max-dt", type=float, default=1e-4)
    args = parser.parse_args()
    NX = NY = args.resolution
    DZ = 2.0 / NX
    NU = args.nu
    END_TIME = args.end_time
    MAX_DT = args.max_dt
    default_case = "static_drop_case" if NX == 80 else f"static_drop_case_{NX}"
    CASE = ROOT / (args.case_name or default_case)
    dx = 2 * DOMAIN_HALF_WIDTH / NX
    alpha = []
    for j in range(NY):
        y = -DOMAIN_HALF_WIDTH + (j + 0.5) * dx
        for i in range(NX):
            x = -DOMAIN_HALF_WIDTH + (i + 0.5) * dx
            alpha.append("1" if x*x + y*y <= R*R else "0")
    liquid_cells = alpha.count("1")

    write("system/blockMeshDict", header("dictionary", "blockMeshDict", f'''
scale 1;
vertices
(
    (-1 -1 {-DZ/2}) (1 -1 {-DZ/2}) (1 1 {-DZ/2}) (-1 1 {-DZ/2})
    (-1 -1 {DZ/2}) (1 -1 {DZ/2}) (1 1 {DZ/2}) (-1 1 {DZ/2})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({NX} {NY} 1) simpleGrading (1 1 1) );
edges ();
boundary
(
    walls {{ type wall; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
'''))

    write("constant/transportProperties", header("dictionary", "transportProperties", f'''
phases (liquidAg gas);
liquidAg {{ transportModel Newtonian; nu {NU}; rho 1; }}
gas      {{ transportModel Newtonian; nu {NU}; rho 1; }}
sigma 1;
'''))
    write("constant/g", header("uniformDimensionedVectorField", "g", '''
dimensions [0 1 -2 0 0 0 0];
value (0 0 0);
'''))
    write("constant/turbulenceProperties", header("dictionary", "turbulenceProperties", "simulationType laminar;"))
    write("0/alpha.liquidAg", header("volScalarField", "alpha.liquidAg", f'''
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar>
{NX*NY}
(
{" ".join(alpha)}
);
boundaryField
{{
    walls {{ type constantAlphaContactAngle; theta0 90; limit gradient; value uniform 0; }}
    frontAndBack {{ type empty; }}
}}
''', "0"))
    write("0/U", header("volVectorField", "U", '''
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField { walls { type noSlip; } frontAndBack { type empty; } }
''', "0"))
    write("0/p_rgh", header("volScalarField", "p_rgh", '''
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField { walls { type fixedFluxPressure; } frontAndBack { type empty; } }
''', "0"))

    write("system/controlDict", header("dictionary", "controlDict", f'''
application interFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {END_TIME};
deltaT 1e-5;
writeControl runTime;
writeInterval 0.02;
purgeWrite 0;
writeFormat ascii;
writePrecision 10;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable no;
adjustTimeStep yes;
maxCo 0.1;
maxAlphaCo 0.1;
maxDeltaT {MAX_DT};
''', "system"))
    write("system/fvSchemes", header("dictionary", "fvSchemes", '''
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U) Gauss upwind;
    div(phi,alpha) Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''', "system"))
    write("system/fvSolution", header("dictionary", "fvSolution", '''
solvers
{
    alpha.liquidAg { nAlphaCorr 1; nAlphaSubCycles 2; cAlpha 1; }
    "pcorr.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    p_rgh { solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0.05; }
    p_rghFinal { $p_rgh; tolerance 1e-7; relTol 0; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }
}
PIMPLE
{
    momentumPredictor no;
    nCorrectors 3;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 0;
}
''', "system"))

    print(f"case={CASE}")
    print(f"grid={NX}x{NY}x1; dx={dx}; liquid_cells={liquid_cells}")
    print(f"initial_area={liquid_cells*dx*dx:.8g}; nominal_circle_area={3.141592653589793*R*R:.8g}")
    print(f"expected_2D_cylinder_laplace_jump=sigma/R={1/R:.8g}")
    print(f"nu={NU}; end_time={END_TIME}; max_delta_t={MAX_DT}")
    print("dimensionless_benchmark_only=true")


if __name__ == "__main__":
    main()
