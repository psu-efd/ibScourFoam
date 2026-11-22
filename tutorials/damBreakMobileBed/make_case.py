"""Case generator: dam break over a mobile bed, Soares-Frazao et al. (2012).

IAHR benchmark (J. Hydraul. Res. 50(4) 364-375): a 1 m breach between
two impervious blocks releases the reservoir (z0 = 0.47 m, Case 1) into a
3.6 m wide flume over an 85 mm layer of uniform sand (d50 1.61 mm, s 2.63,
porosity 0.42) on a rigid floor, held by a rigid sill at x = 9 m. Half flume
(symmetry plane y = 0) or, with --full-width, the whole flume (to test the
symmetry assumption); origin at the gate centre, z = 0 the rigid flume floor.

Mesh: multi-block blockMesh, the dam blocks cut out (body-fitted walls), the
upstream tank (9.2 m wide, treated as symmetric) as a wider block, the mesh
floor ``floor`` 2 cm below the rigid floor so every bed column has cells
below the bed. The sand is an immersed height-field bed on the whole floor
with the non-erodible floor zRigid (0, the sill 0.085), both given as
expressions in constant/immersedBoundaryDict; the bed area mesh is written by
makeFaMesh (system/finite-area/faMeshDefinition).

Usage::

    make_case.py setup  --run-dir DIR [--case 1|2] [--dx 0.05] [--dz 0.01] ...
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2606                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    {loc}object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

# Experiment (paper Sec. 2, Fig. 1, Tables 1-2)
CASES = {1: {"z0": 0.470, "z1": 0.085, "z2": 0.000},
         2: {"z0": 0.510, "z1": 0.150, "z2": 0.150}}
SAND = 0.085          # sand layer thickness [m]
SAND_X = (-1.5, 9.0)  # sand from 1.5 m upstream of the gate to the sill
RAMP = 0.20           # upstream ramp of the sand layer [m]
SILL = (9.0, 9.05)    # rigid sill, top at the sand surface
GATE_HALF = 0.5       # half width of the breach [m]
DAM = (-0.5, 0.5)     # dam blocks in x
HALF_WIDTH = 1.8      # half flume width [m]
TANK = (-12.0, -10.83)  # upstream tank in x (gate 12 m from the flume end)
TANK_HALF = 4.6       # half tank width (9.2 m, taken symmetric)
X_END = 11.0          # outlet of the model (sill at 9 m, overfall beyond)
FLOOR = -0.02         # mesh floor below the rigid floor
TOP = 0.65
GAUGES = {1: [(0.64, -0.500), (0.64, -0.165), (0.64, 0.165), (0.64, 0.500),
              (1.94, -0.990), (1.94, -0.330), (1.94, 0.330), (1.94, 0.990)],
          2: [(0.64, -0.500), (0.64, -0.165), (0.64, 0.165), (0.64, 0.500),
              (2.34, -0.990), (2.34, -0.330), (2.34, 0.330), (2.34, 0.990)]}


def foam(cls: str, obj: str, body: str, loc: str | None = None) -> str:
    """OpenFOAM file text with the standard header."""
    loc_line = f'location    "{loc}";\n    ' if loc else ""
    return HEADER.format(cls=cls, obj=obj, loc=loc_line) + body + \
        "\n// ************************************************************************* //\n"


def write(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def ncells(length: float, d: float) -> int:
    return max(1, round(length / d))


def block_mesh(a: argparse.Namespace) -> tuple[str, int]:
    """Multi-block mesh: x intervals x y bands, dam blocks and the tank band."""
    xs = [TANK[0], TANK[1], SAND_X[0], DAM[0], DAM[1], SAND_X[1], X_END]
    xd = [a.dx_tank, a.dx_res, a.dx, a.dx, a.dx, a.dx_res]
    if a.full_width:
        # y bands: tank | flume, dam | breach | flume, dam | tank
        ys = [-TANK_HALF, -HALF_WIDTH, -GATE_HALF, GATE_HALF, HALF_WIDTH, TANK_HALF]
        yd = [a.dx_tank, a.dx, a.dx, a.dx, a.dx_tank]
        tank_bands, dam_bands = (0, 4), (1, 3)
    else:
        ys = [0.0, GATE_HALF, HALF_WIDTH, TANK_HALF]
        yd = [a.dx, a.dx, a.dx_tank]
        tank_bands, dam_bands = (2,), (1,)
    nx = [ncells(xs[i + 1] - xs[i], xd[i]) for i in range(len(xd))]
    ny = [ncells(ys[j + 1] - ys[j], yd[j]) for j in range(len(yd))]

    def present(i: int, j: int) -> bool:
        if j in tank_bands:
            return i == 0                      # tank band: tank only
        return not (i == 3 and j in dam_bands)  # dam block

    # z: 2 cells below the rigid floor, uniform dz to z_fine, graded to TOP
    nzb = 2
    nzf = ncells(a.z_fine, a.dz)
    dz_top = a.dz_top
    # geometric cells from dz to dz_top over TOP - z_fine
    L = TOP - a.z_fine
    r = dz_top / a.dz
    nzt = max(2, math.ceil(L * math.log(r) / (dz_top - a.dz))) if r > 1 else ncells(L, a.dz)
    nz = nzb + nzf + nzt
    fz = [(0.0 - FLOOR) / (TOP - FLOOR), a.z_fine / (TOP - FLOOR), L / (TOP - FLOOR)]
    zgrad = (f"(({fz[0]:.6f} {nzb} 1) ({fz[1]:.6f} {nzf} 1) "
             f"({fz[2]:.6f} {nzt} {r:.4f}))")

    vid: dict[tuple[int, int, int], int] = {}
    verts: list[str] = []
    for k, z in enumerate((FLOOR, TOP)):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                vid[(i, j, k)] = len(verts)
                verts.append(f"    ({x:g} {y:g} {z:g})")

    blocks: list[str] = []
    patches: dict[str, list[tuple[int, ...]]] = {
        n: [] for n in ("floor", "atmosphere", "symmetry", "walls", "outlet")}
    for j in range(len(yd)):
        for i in range(len(xd)):
            if not present(i, j):
                continue
            v = [vid[(i, j, 0)], vid[(i + 1, j, 0)], vid[(i + 1, j + 1, 0)], vid[(i, j + 1, 0)],
                 vid[(i, j, 1)], vid[(i + 1, j, 1)], vid[(i + 1, j + 1, 1)], vid[(i, j + 1, 1)]]
            blocks.append(f"    hex ({' '.join(map(str, v))}) ({nx[i]} {ny[j]} {nz}) "
                          f"simpleGrading (1 1 {zgrad})")
            patches["floor"].append((v[0], v[3], v[2], v[1]))
            patches["atmosphere"].append((v[4], v[5], v[6], v[7]))
            # x-min / x-max neighbours
            if i == 0 or not present(i - 1, j):
                patches["walls"].append((v[0], v[4], v[7], v[3]))
            if i == len(xd) - 1:
                patches["outlet"].append((v[1], v[2], v[6], v[5]))
            elif not present(i + 1, j):
                patches["walls"].append((v[1], v[2], v[6], v[5]))
            if j == 0:
                patches["walls" if a.full_width else "symmetry"].append((v[0], v[1], v[5], v[4]))
            elif not present(i, j - 1):
                patches["walls"].append((v[0], v[1], v[5], v[4]))
            if j == len(yd) - 1 or not present(i, j + 1):
                patches["walls"].append((v[3], v[7], v[6], v[2]))

    types = {"floor": "wall", "atmosphere": "patch", "symmetry": "symmetryPlane",
             "walls": "wall", "outlet": "patch"}
    bnd = []
    for n, fl in patches.items():
        if not fl:
            continue
        faces = "\n".join(f"            ({' '.join(map(str, f))})" for f in fl)
        bnd.append(f"    {n}\n    {{\n        type {types[n]};\n        faces\n        (\n"
                   f"{faces}\n        );\n    }}")
    ncell = sum(nx[i] * ny[j] * nz for j in range(len(yd)) for i in range(len(xd)) if present(i, j))
    body = (f"// {ncell} cells; dx {a.dx} (sand, breach), {a.dx_res} (reservoir, outlet), "
            f"{a.dx_tank} (tank); dz {a.dz} to z = {a.z_fine}, graded to {dz_top}\n"
            "scale 1;\n\nvertices\n(\n" + "\n".join(verts) + "\n);\n\nblocks\n(\n" +
            "\n".join(blocks) + "\n);\n\nedges ();\n\nboundary\n(\n" + "\n".join(bnd) +
            "\n);\n\nmergePatchPairs ();\n")
    return foam("dictionary", "blockMeshDict", body), ncell


def fields_0(a: argparse.Namespace) -> dict[str, str]:
    """Initial fields (0.orig)."""
    walls = '"(floor|walls)"'
    out = {}
    out["alpha.water"] = ("volScalarField", "[0 0 0 0 0 0 0]", "0", f"""
    {walls} {{ type zeroGradient; }}
    atmosphere {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}
    outlet {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}
    symmetry {{ type symmetryPlane; }}""")
    out["U"] = ("volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", f"""
    {walls} {{ type noSlip; }}
    atmosphere {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    outlet {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
    symmetry {{ type symmetryPlane; }}""")
    out["p_rgh"] = ("volScalarField", "[1 -1 -2 0 0 0 0]", "0", f"""
    {walls} {{ type fixedFluxPressure; value uniform 0; }}
    atmosphere {{ type totalPressure; p0 uniform 0; }}
    outlet {{ type zeroGradient; }}
    symmetry {{ type symmetryPlane; }}""")
    if a.turbulence != "laminar":
        k0, w0 = 1e-6, 1.0
        out["k"] = ("volScalarField", "[0 2 -2 0 0 0 0]", f"{k0:g}", f"""
    {walls} {{ type kqRWallFunction; value uniform {k0:g}; }}
    "(atmosphere|outlet)" {{ type inletOutlet; inletValue uniform {k0:g}; value uniform {k0:g}; }}
    symmetry {{ type symmetryPlane; }}""")
        out["omega"] = ("volScalarField", "[0 0 -1 0 0 0 0]", f"{w0:g}", f"""
    {walls} {{ type omegaWallFunction; value uniform {w0:g}; }}
    "(atmosphere|outlet)" {{ type inletOutlet; inletValue uniform {w0:g}; value uniform {w0:g}; }}
    symmetry {{ type symmetryPlane; }}""")
        out["nut"] = ("volScalarField", "[0 2 -1 0 0 0 0]", "0", f"""
    {walls} {{ type nutkWallFunction; value uniform 0; }}
    "(atmosphere|outlet)" {{ type calculated; value uniform 0; }}
    symmetry {{ type symmetryPlane; }}""")
    return {n: foam(c, n, f"dimensions      {d};\n\ninternalField   uniform {v};\n\n"
                    f"boundaryField\n{{{b}\n}}\n", "0")
            for n, (c, d, v, b) in out.items()}


def setup(a: argparse.Namespace) -> None:
    run = a.run_dir
    if os.path.exists(run) and os.listdir(run):
        sys.exit(f"{run} exists and is not empty")
    c = CASES[a.case]
    bm, ncell = block_mesh(a)
    write(os.path.join(run, "system", "blockMeshDict"), bm)

    for n, t in fields_0(a).items():
        write(os.path.join(run, "0.orig", n), t)

    # water: reservoir to the gate plane x = 0 at z0; downstream z1 over the
    # sand and z2 beyond the sill (Case 1: z1 = sand top, i.e. no water).
    # The boxes start at the rigid floor or the sand surface: the buried
    # cells below the bed hold no water (the interfaceHeight gauges integrate
    # alpha over the whole column)
    def box(x0: float, x1: float, z0: float, z1: float) -> str:
        return (f"boxToCell {{ box ({x0:g} {-TANK_HALF - 1:g} {z0:g}) ({x1:g} {TANK_HALF + 1} {z1:g}); "
                "fieldValues ( volScalarFieldValue alpha.water 1 ); }")

    regions = [box(TANK[0] - 1, SAND_X[0], 0.0, c["z0"]),
               box(SAND_X[0], 0.0, SAND, c["z0"])]
    if c["z1"] > SAND + 1e-9:
        regions.append(box(0.0, SILL[1], SAND, c["z1"]))
    if c["z2"] > 1e-9:
        regions.append(box(SILL[1], X_END + 1, 0.0, c["z2"]))
    write(os.path.join(run, "system", "setFieldsDict"), foam(
        "dictionary", "setFieldsDict",
        "defaultFieldValues ( volScalarFieldValue alpha.water 0 );\n\nregions\n(\n    " +
        "\n    ".join(regions) + "\n);\n"))

    gauges = " ".join(f"({x:g} {y if a.full_width else abs(y):g} 0)" for x, y in GAUGES[a.case])
    write(os.path.join(run, "system", "controlDict"), foam("dictionary", "controlDict", f"""\
application     ibScourInterFoam;
libs            (libibScour);
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {a.end_time:g};
deltaT          1e-4;
writeControl    adjustableRunTime;
writeInterval   {a.write_interval:g};
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
adjustTimeStep  yes;
maxCo           {a.max_co:g};
maxAlphaCo      {a.max_co:g};
maxDeltaT       0.005;

functions
{{
    // water depth at the gauges of Table 2 (|y|, half flume): the integral
    // of the ghost-free volume fraction (solid cells zero) over the column,
    // so the water surface is z_b + hB
    gauges
    {{
        type            interfaceHeight;
        libs            (fieldFunctionObjects);
        alpha           alpha.waterFluid;
        locations       ({gauges});
        writeControl    adjustableRunTime;
        writeInterval   0.05;
    }}
    // water leaving the domain (for the water budget): flux of the water
    // phase through the outlet and the atmosphere [m3/s]
    waterOut
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            outlet;
        operation       sum;
        fields          (alphaPhi0.water);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   1;
    }}
    waterTop
    {{
        $waterOut;
        name            atmosphere;
    }}
}}
"""))
    write(os.path.join(run, "system", "fvSchemes"), foam("dictionary", "fvSchemes", """\
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(phiC,C)     Gauss limitedLinear01 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
wallDist        { method meshWave; }
"""))
    write(os.path.join(run, "system", "fvSolution"), foam("dictionary", "fvSolution", """\
solvers
{
    "alpha.water.*"
    {
        alphaScheme     MULES;
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;
        MULESCorr       yes;
        nLimiterIter    5;
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
    p_rgh
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-8;
        relTol          0.01;
    }
    p_rghFinal { $p_rgh; tolerance 1e-9; relTol 0; }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
    "(U|k|omega)Final" { $U; }
    "C|CFinal" { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
}
PIMPLE
{
    momentumPredictor no;
    nOuterCorrectors 1;
    nCorrectors     3;
    nNonOrthogonalCorrectors 0;
}
"""))
    write(os.path.join(run, "system", "decomposeParDict"), foam(
        "dictionary", "decomposeParDict",
        f"numberOfSubdomains {a.n_procs};\nmethod          scotch;\n"))
    write(os.path.join(run, "system", "finite-area", "faSchemes"), foam("dictionary", "faSchemes", """\
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""))
    write(os.path.join(run, "system", "finite-area", "faSolution"), foam(
        "dictionary", "faSolution", "solvers {}\n"))

    write(os.path.join(run, "constant", "g"), foam(
        "uniformDimensionedVectorField", "g",
        "dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 0 -9.81);\n"))
    write(os.path.join(run, "constant", "transportProperties"), foam(
        "dictionary", "transportProperties",
        "phases (water air);\n"
        "water { transportModel Newtonian; nu 1e-06; rho 1000; }\n"
        "air   { transportModel Newtonian; nu 1.48e-05; rho 1; }\n"
        "sigma           0.07;\n"))
    if a.turbulence == "laminar":
        turb = "simulationType  laminar;\n"
    else:
        turb = ("simulationType  RAS;\nRAS\n{\n    RASModel        kOmegaSST;\n"
                "    turbulence      on;\n    printCoeffs     on;\n}\n")
        write(os.path.join(run, "constant", "fvOptions"), foam("dictionary", "fvOptions", """\
// Larsen & Fuhrman (2018) eddy-viscosity limiter (lambda2 0.05, the paper's
// value) and the Devolder et al. (2017) buoyancy term
stabilisation
{
    type            multiphaseStabilizedTurbulence;
    active          yes;
    multiphaseStabilizedTurbulenceCoeffs
    {
        rho             rho;
        Cmu             0.09;
        C               1.51;
        lambda2         0.05;
        alpha           1.36;
    }
}
"""))
    write(os.path.join(run, "constant", "turbulenceProperties"), foam(
        "dictionary", "turbulenceProperties", turb))

    if a.suspended != "none":
        # reference level: 2 d50 for Zyserman & Fredsoe (1994), ks for van Rijn
        delta = 2 * 0.00161 if a.suspended == "zysermanFredsoe" else a.ks
        write(os.path.join(run, "constant", "sedimentProperties"), foam(
            "dictionary", "sedimentProperties", f"""\
solveSediment   yes;
diffusivity     turbulence;     // DC = nu/ScC + nut/sigmaC
ScC             1;
sigmaC          1;
d50             0.00161;
s               2.63;
settlingVelocityModel  vanRijn1984;
// reference concentration {a.suspended} at delta = {delta:g} m, Shields number
// from the total bed shear; Rouse ~3: log-law exchange from delta to the IB cell
referenceConcentrationModel  {a.suspended};
delta           {delta:g};
shieldsNumber   total;
ibBedExchange
{{
    mode        logLaw;
}}
"""))
        write(os.path.join(run, "0.orig", "C"), foam("volScalarField", "C", """\
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    "(floor|walls)" { type zeroGradient; }
    "(atmosphere|outlet)" { type inletOutlet; inletValue uniform 0; value uniform 0; }
    symmetry { type symmetryPlane; }
}
""", "0"))

    mobile = not a.rigid_bed
    gate_entry = ""
    if a.gate != "none":
        th = a.dx
        write(os.path.join(run, "constant", "triSurface", "gate.stl"),
              box_stl((-th, -GATE_HALF - 0.05 if a.full_width else -0.05, SAND),
                      (0.0, GATE_HALF + 0.05, TOP + 0.10)))
        gate_entry = (f"    // gate lifted vertically ({a.gate}, bottom clears the reservoir\n"
                      f"    // surface at t = {a.gate_time:g} s); one cell thick, upstream of x = 0\n"
                      "    gate\n    {\n        type            triSurface;\n"
                      '        file            "gate.stl";\n        closed          true;\n'
                      f"        displacement    {gate_motion(a, c['z0'])};\n    }}\n")
    write(os.path.join(run, "constant", "immersedBoundaryDict"), foam(
        "dictionary", "immersedBoundaryDict", f"""\
// Sand layer of the Soares-Frazao et al. (2012) flume as a height-field bed on
// the mesh floor (2 cm below the rigid flume floor z = 0): sand 0.085 m from
// x = -1.5 (0.20 m ramp) to the sill, rigid floor z_r = 0, sill z_r = 0.085
// at 9.00-9.05 m.
cutFaces        false;
wallSink        tangential;
fluidAlphaField true;           // alpha.waterFluid for the gauges

bed
{{
    type            heightField;
    floorPatch      floor;
    initialElevation
    {{
        type        expression;
        expression  #{{
            (pos().x() >= {SILL[0]:g} && pos().x() <= {SILL[1]:g}) ? {SAND:g}
          : (pos().x() >= {SAND_X[0]:g} && pos().x() < {SAND_X[1]:g})
            ? min({SAND:g}, {SAND:g}*(pos().x() - ({SAND_X[0]:g}))/{RAMP:g}) : 0
        #}};
    }}
    rigidElevation
    {{
        type        expression;
        expression  #{{
            (pos().x() >= {SILL[0]:g} && pos().x() <= {SILL[1]:g}) ? {SAND:g} : 0
        #}};
    }}
    openBoundaries  ();

    grainDiameter   0.00161;
    specificGravity 2.63;
    porosity        0.42;
    gravity         9.81;
    criticalShields {a.critical_shields:g};
    frictionCoefficient {math.tan(math.radians(a.repose)):.4f};   // repose {a.repose:g} deg
    bedLoad         {a.bed_load};
    bedSlope        {{ model kochFlokstra; alpha 1.5; }}
    adaptationLength {a.adaptation_length:g};

    morphStartTime  {0 if mobile else 1e9};
    morphFactor     1;
    maxBedChange    5e-4;
    updateDistance  1e-3;
    sandSlide       true;
    slideDiffusivity 5e-4;
    slideTolerance  0.02;
    slideMaxIter    500;
    wetThreshold    0.5;
    seepage         {{ type none; }}
}}

solids
{{
{gate_entry}}}

turbulence {{ kTreatment equilibrium; }}

wallModel
{{
    mode            stress;
    sampleBelow     1;
    wallLaw         {{ type roughLogLaw; ks {a.ks:g}; }}
}}
"""))
    # The bed area mesh for makeFaMesh (all edges closed: no open boundaries)
    write(os.path.join(run, "system", "finite-area", "faMeshDefinition"), foam(
        "dictionary", "faMeshDefinition",
        "polyMeshPatches ( floor );\n\nboundary\n{\n}\n\n"
        "defaultPatch\n{\n    name    bedBoundary;\n    type    patch;\n}\n"))
    setup_json = dict(vars(a))
    setup_json.update({"cells": ncell, "experiment": CASES[a.case]})
    write(os.path.join(run, "case_setup.json"), json.dumps(setup_json, indent=2, default=str))
    shutil.copy(os.path.abspath(__file__), os.path.join(run, "make_case.py"))
    print(f"{run}: Case {a.case}, {ncell} cells, {'rigid' if a.rigid_bed else 'mobile'} bed, "
          f"{a.turbulence}")


def box_stl(lo: tuple, hi: tuple) -> str:
    """ASCII STL of an axis-aligned box (outward normals)."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    p = {(i, j, k): ((x0, x1)[i], (y0, y1)[j], (z0, z1)[k])
         for i in (0, 1) for j in (0, 1) for k in (0, 1)}
    quads = [  # (outward normal, corners counter-clockwise seen from outside)
        ((-1, 0, 0), [(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)]),
        ((1, 0, 0), [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)]),
        ((0, -1, 0), [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)]),
        ((0, 1, 0), [(0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)]),
        ((0, 0, -1), [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)]),
        ((0, 0, 1), [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]),
    ]
    out = ["solid gate"]
    for n, q in quads:
        for tri in ((q[0], q[1], q[2]), (q[0], q[2], q[3])):
            out.append(f"  facet normal {n[0]} {n[1]} {n[2]}\n    outer loop")
            out += [f"      vertex {p[c][0]:g} {p[c][1]:g} {p[c][2]:g}" for c in tri]
            out.append("    endloop\n  endfacet")
    out.append("endsolid gate")
    return "\n".join(out) + "\n"


def gate_motion(a: argparse.Namespace, z0: float) -> str:
    """Displacement table of the lifted gate: the bottom clears the reservoir
    surface z0 at --gate-time (paper: 0.23 s, the time the gate leaves the
    water), uniformly accelerated (counterweight) or at constant speed, then
    rises on to 0.75 m above the sand and stops."""
    lift = z0 - SAND
    total = 0.75
    pts = []
    if a.gate == "accel":
        acc = 2 * lift / a.gate_time ** 2
        t_end = math.sqrt(2 * total / acc)
        n = 40
        pts = [(t_end * i / n, 0.5 * acc * (t_end * i / n) ** 2) for i in range(n + 1)]
    else:
        v = lift / a.gate_time
        pts = [(0.0, 0.0), (total / v, total)]
    pts.append((pts[-1][0] + 1.0, total))
    rows = " ".join(f"({t:.5f} (0 0 {d:.6f}))" for t, d in pts)
    return f"table ({rows})"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("setup")
    s.add_argument("--run-dir", required=True)
    s.add_argument("--case", type=int, choices=(1, 2), default=1)
    s.add_argument("--dx", type=float, default=0.05, help="cell size over the sand and breach [m]")
    s.add_argument("--dx-res", type=float, default=0.10, help="reservoir and outlet [m]")
    s.add_argument("--dx-tank", type=float, default=0.20, help="upstream tank [m]")
    s.add_argument("--dz", type=float, default=0.01, help="vertical cell size up to --z-fine [m]")
    s.add_argument("--z-fine", type=float, default=0.20)
    s.add_argument("--dz-top", type=float, default=0.04)
    s.add_argument("--turbulence", choices=("kOmegaSST", "laminar"), default="kOmegaSST")
    s.add_argument("--rigid-bed", action="store_true", help="bed frozen (hydrodynamics only)")
    s.add_argument("--full-width", action="store_true",
                   help="whole flume instead of the half flume with a symmetry plane")
    s.add_argument("--bed-load", choices=("meyerPeterMuller", "engelundFredsoe"),
                   default="meyerPeterMuller")
    s.add_argument("--critical-shields", type=float, default=0.047)
    s.add_argument("--suspended", choices=("none", "zysermanFredsoe", "vanRijn1984"), default="none",
                   help="suspended load: reference concentration model (log-law exchange); "
                        "van Rijn (1984) is unbounded at Shields numbers of order 1")
    s.add_argument("--repose", type=float, default=32.0)
    s.add_argument("--adaptation-length", type=float, default=0.0,
                   help="bed-load adaptation (lag) length [m], 0 = local equilibrium")
    s.add_argument("--ks", type=float, default=0.004, help="bed roughness (2.5 d50) [m]")
    s.add_argument("--gate", choices=("none", "accel", "linear"), default="none",
                   help="model the vertical gate lift (none = instantaneous removal)")
    s.add_argument("--gate-time", type=float, default=0.23,
                   help="time at which the gate bottom leaves the water [s]")
    s.add_argument("--end-time", type=float, default=20.0)
    s.add_argument("--write-interval", type=float, default=0.5)
    s.add_argument("--max-co", type=float, default=0.5)
    s.add_argument("--n-procs", type=int, default=32)
    a = ap.parse_args()
    setup(a)


if __name__ == "__main__":
    main()
