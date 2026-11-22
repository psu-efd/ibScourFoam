#!/usr/bin/env python3
"""Generate the body-fitted hump-migration test ``deformingBedHump``.

Quasi-2D channel on a body-fitted mesh that deforms with the bed
(``deformingBed`` function object, stock ``pimpleFoam`` with a
``dynamicMotionSolverFvMesh``). The bed starts with the Gaussian hump of
the immersed-boundary test ``humpMigration``:

    z_b(x) = A exp(-(x - x_h)^2 / (2 sigma^2)),
    A = 0.01 m (0.1 h), sigma = 0.05 m, x_h = 0.4 m, h = 0.1 m,

built into the mesh (the columns are stretched between the bed and the
rigid lid); deformingBed takes the initial z_b from the bed faces. U_b = 0.45 m/s gives a flat-bed Shields number
theta_0 ~ 0.1, so the whole bed is mobile. A fixed strip at the inlet
keeps the inlet bed level under the prescribed inflow profile.

Analytical celerity (quasi-steady rigid-lid flow, constant friction
factor, small amplitude), see analytic_celerity in humpMigration/check.py:

    c = MF / (1 - n) * dq_b/dtheta * 2 theta_0 / H.

Usage::

    make_case.py setup --run-dir DIR [--morph-factor 20] [--amplitude 0.01]
    blockMesh -case DIR
    make_case.py init --run-dir DIR
"""

import argparse
import logging
import os
import sys
from typing import Dict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import case_setup  # noqa: E402
import channel_profiles as prof  # noqa: E402
import foam_file  # noqa: E402

LOG = logging.getLogger("deformingBedHump.make_case")

NU = 1.0e-6            # kinematic viscosity [m^2/s]
SETUP = {
    "depth": 0.1,      # flow depth over the flat bed H [m]
    "u_bulk": 0.45,    # depth-mean velocity [m/s]
    "length": 1.2,     # channel length [m]
    "width": 5e-3,     # span (one cell, empty sides) [m]
    "dx": 5e-3,        # cell size in x [m]
    "dz": 5e-3,        # cell height [m]
    "amplitude": 0.01,  # hump height A [m]
    "sigma": 0.05,     # hump standard deviation [m]
    "x_hump": 0.4,     # initial hump centre [m]
    "d50": 2.6e-4,     # grain diameter [m]
    "s": 2.65,         # specific gravity [-]
    "theta_c0": 0.05,  # critical Shields number (flat bed) [-]
    "porosity": 0.4,   # bed porosity [-]
    "fixed_inlet": 0.05,  # fixed bed strip at the inlet [m]
    "delta_t": 0.004,  # time step [s]
    "spinup_time": 5.0,  # hydrodynamic spin-up [s]
}

FIELDS = {
    "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0.45 0 0)", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " value uniform (0 0 0);",
        "outlet": "type inletOutlet; inletValue uniform (0 0 0);"
                  " value uniform (0 0 0);",
        "top": "type slip;",
        "bed": "type movingWallVelocity; value uniform (0 0 0);",
    }),
    "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        "inlet": "type zeroGradient;",
        "outlet": "type fixedValue; value uniform 0;",
        "top": "type slip;",
        "bed": "type zeroGradient;",
    }),
    "k": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 1e-3", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " value uniform 1e-3;",
        "outlet": "type inletOutlet; inletValue uniform 1e-3;"
                  " value uniform 1e-3;",
        "top": "type slip;",
        "bed": "type kqRWallFunction; value uniform 1e-3;",
    }),
    "omega": ("volScalarField", "[0 0 -1 0 0 0 0]", "uniform 10", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " value uniform 10;",
        "outlet": "type inletOutlet; inletValue uniform 10;"
                  " value uniform 10;",
        "top": "type slip;",
        "bed": "type omegaWallFunction; blending binomial; n 2;"
               " value uniform 10;",
    }),
    "nut": ("volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        "inlet": "type calculated; value uniform 0;",
        "outlet": "type calculated; value uniform 0;",
        "top": "type calculated; value uniform 0;",
        "bed": "type nutkRoughWallFunction; Ks uniform {ks}; Cs uniform 0.5;"
               " value uniform 0;",
    }),
    "pointDisplacement": ("pointVectorField", "[0 1 0 0 0 0 0]",
                          "uniform (0 0 0)", {
        "inlet": "type slip;",
        "outlet": "type slip;",
        "top": "type fixedValue; value uniform (0 0 0);",
        "bed": "type fixedValue; value uniform (0 0 0);",
    }),
}


def bed_elevation(x: np.ndarray, setup: Dict) -> np.ndarray:
    """Return the initial bed elevation z_b(x) [m]."""
    return setup["amplitude"] * np.exp(
        -0.5 * ((x - setup["x_hump"]) / setup["sigma"]) ** 2)


# Bed conditions of the immersed-boundary wall law (--wall wallLaw); the
# wall-cell omega (and k) are fixed by the wallLawTurbulence constraint
WALL_LAW_BED = {
    "nut": """\
        type            nutWallLaw;
        wallLaw
        {{
            type        roughLogLaw;
            ks          {ks};
        }}
        value           uniform 0;""",
    "omega": """\
        type            zeroGradient;""",
}


def field_file(name: str, setup: Dict) -> str:
    """Return the body of the initial field file ``name``."""
    _, dims, internal, bcs = FIELDS[name]
    ks = "{:.6g}".format(2.5 * setup["d50"])
    lines = ["dimensions      {};".format(dims), "",
             "internalField   {};".format(internal), "", "boundaryField", "{"]
    for patch, spec in bcs.items():
        lines += ["    {}".format(patch), "    {"]
        if patch == "bed" and setup.get("wall") == "wallLaw" \
                and name in WALL_LAW_BED:
            lines += [WALL_LAW_BED[name].format(ks=ks)]
        else:
            spec = spec.replace("{ks}", ks)
            lines += ["        {}".format(e.strip() + ";")
                      for e in spec.split(";") if e.strip()]
        lines += ["    }"]
    lines += ["    frontAndBack", "    {", "        type            empty;",
              "    }", "}"]
    return "\n".join(lines)


def write_dict(run_dir: str, rel: str, obj: str, body: str) -> None:
    """Write dictionary ``rel`` (relative to the case) with header."""
    foam_file.write_foam_file(os.path.join(run_dir, rel), "dictionary", obj,
                              body)


def cmd_setup(args: argparse.Namespace) -> None:
    """Create the run directory (before meshing)."""
    setup = dict(SETUP)
    setup.update({"amplitude": args.amplitude, "sigma": args.sigma,
                  "x_hump": args.x_hump, "length": args.length,
                  "dx": args.dx, "wall": args.wall,
                  "delta_t": args.dt, "spinup_time": args.spinup,
                  "morph_factor": args.morph_factor,
                  "morph_time": args.morph_time, "nu": NU,
                  "adaptation_length": args.adaptation_length})
    nx = int(round(setup["length"] / setup["dx"]))
    nz = int(round(setup["depth"] / setup["dz"]))
    end_time = setup["spinup_time"] + setup["morph_time"] / args.morph_factor
    run = args.run_dir
    os.makedirs(run, exist_ok=True)
    case_setup.save_setup(run, setup)

    write_dict(run, "system/blockMeshDict", "blockMeshDict", """\
scale   1;

vertices
(
    (0 0 0) ({L} 0 0) ({L} {W} 0) (0 {W} 0)
    (0 0 {H}) ({L} 0 {H}) ({L} {W} {H}) (0 {W} {H})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} 1 {nz}) simpleGrading (1 1 1)
);

boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    top    {{ type patch; faces ((4 5 6 7)); }}
    bed    {{ type wall;  faces ((0 3 2 1)); }}
    frontAndBack {{ type empty; faces ((0 1 5 4) (3 7 6 2)); }}
);
""".format(L=setup["length"], W=setup["width"], H=setup["depth"], nx=nx,
           nz=nz))

    write_dict(run, "system/controlDict", "controlDict", """\
application     pimpleFoam;
libs            (ibScour);
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {end:g};
deltaT          {dt:g};
writeControl    adjustableRunTime;
writeInterval   0.5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;

functions
{{
    bed
    {{
        type            deformingBed;
        libs            (ibScour);
        bedPatch        bed;

        grainDiameter   {d50:g};
        specificGravity {s:g};
        porosity        {n:g};
        criticalShields {tc:g};
        frictionCoefficient 0.63;
        bedLoad         engelundFredsoe;
        morphFactor     {mf:g};
        morphStartTime  {t0:g};
        maxBedChange    1e-4;
        bedCourant      {bc:g};
        openBoundaries  (inlet outlet);
        fixedBed
        {{
            inlet       {fix:g};
        }}
        sandSlide       true;
        slideDiffusivity 5e-4;
        adaptationLength {la:g};
    }}
}}
""".format(la=args.adaptation_length, end=end_time, dt=setup["delta_t"], d50=setup["d50"], s=setup["s"],
           n=setup["porosity"], tc=setup["theta_c0"], mf=args.morph_factor,
           bc=args.bed_courant,
           t0=setup["spinup_time"], fix=setup["fixed_inlet"]))

    write_dict(run, "system/fvSchemes", "fvSchemes", """\
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
wallDist        { method meshWave; }
""")

    write_dict(run, "system/fvSolution", "fvSolution", """\
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }
    pFinal
    {
        $p;
        relTol          0;
    }
    pcorr
    {
        $p;
        tolerance       1e-8;
        relTol          0;
    }
    pcorrFinal
    {
        $pcorr;
    }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }
    "(U|k|omega)Final"
    {
        $U;
        relTol          0;
    }
    cellDisplacement
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-9;
        relTol          0;
    }
    cellDisplacementFinal
    {
        $cellDisplacement;
    }
}

PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors     2;
    nNonOrthogonalCorrectors 1;
    correctPhi      yes;
    moveMeshOuterCorrectors no;
    checkMeshCourantNo yes;
}
""")

    write_dict(run, "system/decomposeParDict", "decomposeParDict", """\
numberOfSubdomains 3;
method          simple;
coeffs
{
    n           (3 1 1);
}
""")

    for sub in ("bedRegion/finite-area",):
        write_dict(run, "system/{}/faSchemes".format(sub), "faSchemes", """\
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
""")
        write_dict(run, "system/{}/faSolution".format(sub), "faSolution",
                   "solvers {}\n")

    if args.wall == "wallLaw":
        write_dict(run, "constant/fvOptions", "fvOptions", """\
wallTurbulence
{{
    type            wallLawTurbulence;
    patches         (bed);
    wallLaw
    {{
        type        roughLogLaw;
        ks          {ks:g};
    }}
}}
""".format(ks=2.5 * setup["d50"]))

    write_dict(run, "constant/dynamicMeshDict", "dynamicMeshDict", """\
dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolverLibs (fvMotionSolvers);
motionSolver    displacementLaplacian;
displacementLaplacianCoeffs
{
    diffusivity     quadratic inverseDistance (bed);
}
""")
    write_dict(run, "constant/transportProperties", "transportProperties",
               "transportModel  Newtonian;\nnu              {:g};\n".format(NU))
    write_dict(run, "constant/turbulenceProperties", "turbulenceProperties",
               """\
simulationType  RAS;
RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
""")
    LOG.info("setup %s: %d x %d cells, end time %g s", run, nx, nz, end_time)


def cmd_init(args: argparse.Namespace) -> None:
    """Write the initial fields and inflow data, put the hump in the mesh."""
    run = args.run_dir
    setup = case_setup.load_setup(run)
    mesh = os.path.join(run, "constant", "polyMesh")
    points = foam_file.read_points(os.path.join(mesh, "points"))
    offsets, flat = foam_file.read_faces(os.path.join(mesh, "faces"))
    boundary = foam_file.read_boundary(os.path.join(mesh, "boundary"))

    def centres(patch: str) -> np.ndarray:
        start = int(boundary[patch]["startFace"])
        n = int(boundary[patch]["nFaces"])
        return np.array([points[flat[offsets[f]:offsets[f + 1]]].mean(axis=0)
                         for f in range(start, start + n)])

    for name in FIELDS:
        cls = FIELDS[name][0]
        foam_file.write_foam_file(os.path.join(run, "0", name), cls, name,
                                  field_file(name, setup))

    inlet = centres("inlet")
    u, k, omega = prof.turbulent_profile(inlet[:, 2], setup["depth"],
                                         setup["u_bulk"], NU)
    vel = np.zeros((len(u), 3))
    vel[:, 0] = u
    case_setup.write_boundary_data(run, "inlet", inlet,
                                   {"U": vel, "k": k, "omega": omega})

    # Initial bathymetry in the mesh: stretch the columns between the bed
    # and the rigid lid (deformingBed takes z_b from the bed-face centres)
    depth = setup["depth"]
    zb = bed_elevation(points[:, 0], setup)
    points[:, 2] = zb + points[:, 2] * (depth - zb) / depth
    foam_file.write_foam_file(os.path.join(mesh, "points"), "vectorField",
                              "points", foam_file.format_vector_list(points))
    LOG.info("init %s: bed hump in the mesh, z_b max %.4g m", run, zb.max())


def main() -> None:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("setup")
    s.add_argument("--run-dir", required=True)
    s.add_argument("--morph-factor", type=float, default=20.0)
    s.add_argument("--morph-time", type=float, default=200.0)
    s.add_argument("--bed-courant", type=float, default=0.0,
                   help="bed Courant number of the adaptive morphological "
                        "factor (0: fixed factor)")
    s.add_argument("--amplitude", type=float, default=SETUP["amplitude"])
    s.add_argument("--sigma", type=float, default=SETUP["sigma"])
    s.add_argument("--x-hump", type=float, default=SETUP["x_hump"])
    s.add_argument("--length", type=float, default=SETUP["length"])
    s.add_argument("--dx", type=float, default=SETUP["dx"])
    s.add_argument("--dt", type=float, default=SETUP["delta_t"])
    s.add_argument("--spinup", type=float, default=SETUP["spinup_time"])
    s.add_argument("--adaptation-length", type=float, default=0.0,
                   help="bed-load adaptation (saturation) length [m]; "
                        "0 = local equilibrium bed load")
    s.add_argument("--wall", choices=("wallFunction", "wallLaw"),
                   default="wallFunction",
                   help="bed treatment: OpenFOAM rough wall functions or "
                        "the immersed-boundary wall law (nutWallLaw + "
                        "wallLawTurbulence)")
    i = sub.add_parser("init")
    i.add_argument("--run-dir", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    {"setup": cmd_setup, "init": cmd_init}[args.cmd](args)


if __name__ == "__main__":
    main()
