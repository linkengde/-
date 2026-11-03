#!/usr/bin/env python3
"""Write a dimensionless pressure-driven channel-flow verification case."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "poiseuille_case"
NX, NY = 120, 40
LENGTH, HEIGHT, THICKNESS = 12.0, 1.0, 1.0 / NY
NU = 0.01


def foam_header(cls: str, obj: str, location: str) -> str:
    return f'''FoamFile
{{
    version 2.0;
    format ascii;
    class {cls};
    location "{location}";
    object {obj};
}}
'''


def put(rel: str, body: str, cls: str, obj: str, location: str) -> None:
    path = CASE / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(foam_header(cls, obj, location) + body, encoding="ascii")


def main() -> None:
    put("system/blockMeshDict", f'''
scale 1;
vertices
(
    (0 0 {-THICKNESS/2}) ({LENGTH} 0 {-THICKNESS/2}) ({LENGTH} {HEIGHT} {-THICKNESS/2}) (0 {HEIGHT} {-THICKNESS/2})
    (0 0 {THICKNESS/2}) ({LENGTH} 0 {THICKNESS/2}) ({LENGTH} {HEIGHT} {THICKNESS/2}) (0 {HEIGHT} {THICKNESS/2})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({NX} {NY} 1) simpleGrading (1 1 1) );
edges ();
boundary
(
    inlet {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    walls {{ type wall; faces ((0 1 5 4) (3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
''', "dictionary", "blockMeshDict", "system")

    put("constant/transportProperties", f'''
transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {NU};
''', "dictionary", "transportProperties", "constant")
    put("constant/turbulenceProperties", "simulationType laminar;\n", "dictionary", "turbulenceProperties", "constant")
    put("0/U", '''
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type zeroGradient; }
    walls { type noSlip; }
    frontAndBack { type empty; }
}
''', "volVectorField", "U", "0")
    put("0/p", '''
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type fixedValue; value uniform 1; }
    outlet { type fixedValue; value uniform 0; }
    walls { type zeroGradient; }
    frontAndBack { type empty; }
}
''', "volScalarField", "p", "0")
    put("system/controlDict", '''
application simpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 5000;
deltaT 1;
writeControl timeStep;
writeInterval 250;
purgeWrite 0;
writeFormat ascii;
writePrecision 10;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable no;
''', "dictionary", "controlDict", "system")
    put("system/fvSchemes", '''
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''', "dictionary", "fvSchemes", "system")
    put("system/fvSolution", '''
solvers
{
    p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    residualControl { p 1e-8; U 1e-8; }
}
relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; }
}
''', "dictionary", "fvSolution", "system")
    print(f"case={CASE}; grid={NX}x{NY}x1; L={LENGTH}; H={HEIGHT}; nu={NU}")
    print(f"dimensionless_pressure_drop=1; expected_u_max={1 / (8 * NU * LENGTH):.9g}")


if __name__ == "__main__":
    main()
