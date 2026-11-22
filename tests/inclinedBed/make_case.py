#!/usr/bin/env python3
"""Generate the 2D inclined-bed test ``inclinedBed`` (IB vs body-fitted).

Bed shear stress on sloping beds: a 2D channel (one cell across, symmetry
sides) whose frozen bed has a trapezoidal pit,

    flat (inflow development) -> down-ramp at alpha -> flat bottom
    -> up-ramp at alpha -> flat (outflow),

of depth d below the flat bed z = 0, with a rigid lid at z = h. The same
geometry is built two ways at matched resolution (dx, dz):

* ``ib``: Cartesian box from below the pit bottom to the lid; the bed is
  the immersed-boundary height field (``heightField`` on the floor patch,
  initial z_b in ``0/finite-area/zb``), solved with ``ibScourPimpleFoam``
  and the wall model of the Roulund pier runs (stress mode, vanDriest);
* ``bf``: body-fitted reference, columns stretched between the bed and the
  lid, wall law ``nutWallLaw`` + ``wallLawTurbulence`` (vanDriest, same
  ks), stock ``pimpleFoam``; ``deformingBed`` (frozen, no mesh motion)
  computes the bed shear stress with the same convention as the sediment
  core.

Both write ``postProcessing/bed/<t>/bed.vtp`` (z_b, tauBed), which
``compare.py`` uses. The bed is frozen (morphology never starts).

Usage::

    make_case.py setup --variant ib|bf --slope-deg 30 --run-dir DIR
                       [--laminar --u-bulk 0.04 --nu 1e-5 --ks 0]
    blockMesh -case DIR
    make_case.py init --run-dir DIR
    decomposePar -case DIR                       # parallel runs
    make_case.py decompose-zb --run-dir DIR      # IB: z_b per processor

With ``--laminar`` the flow is laminar (no k, omega, nut; no wall law on
the BF bed) and the inflow is the half-Poiseuille profile of a channel
with a no-slip bed and a slip lid.
"""
from __future__ import annotations

import argparse
import logging
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import case_setup  # noqa: E402
import channel_profiles as prof  # noqa: E402
import foam_file  # noqa: E402

LOG = logging.getLogger("inclinedBed.make_case")

NU = 1.0e-6  # kinematic viscosity [m^2/s]

SETUP = {
    "depth": 0.1,         # lid height above the flat bed h [m]
    "u_bulk": 0.4,        # depth-mean inflow velocity [m/s]
    "length": 1.8,        # channel length [m]
    "dx": 2e-3,           # cell size in x [m] (IB pier runs: 2.1 mm)
    "dz": 1e-3,           # cell height [m] (IB pier runs: 1 mm)
    "pit_start": 0.6,     # x of the top of the down-ramp [m]
    "pit_depth": 0.03,    # pit depth d [m]
    "pit_bottom": 0.15,   # length of the flat pit bottom [m]
    "ks": 7e-4,           # equivalent sand roughness [m] (Roulund runs)
    "d50": 2.6e-4,        # grain diameter (sediment core inputs) [m]
    "delta_t": 0.002,     # time step [s]
    "end_time": 15.0,     # ~3.3 flow-throughs [s]
    "write_interval": 5.0,
    "below": 4,           # IB: cells below the pit bottom
}


def bed_elevation(x: np.ndarray, s: dict) -> np.ndarray:
    """Return z_b(x) [m] of the trapezoidal pit (z_b = 0 outside)."""
    t = math.tan(math.radians(s["slope_deg"]))
    d = s["pit_depth"]
    xa = s["pit_start"]
    xd = xa + 2 * d / t + s["pit_bottom"]
    return -np.clip(np.minimum((x - xa) * t, (xd - x) * t), 0.0, d)


def patch_centres(run: str, patch: str) -> np.ndarray:
    """Face centres of ``patch`` of the case mesh (point averages)."""
    mesh = os.path.join(run, "constant", "polyMesh")
    return mesh_patch_centres(mesh, patch)


def mesh_patch_centres(mesh: str, patch: str) -> np.ndarray:
    """Face centres of ``patch`` of the polyMesh directory ``mesh``."""
    points = foam_file.read_points(os.path.join(mesh, "points"))
    offsets, flat = foam_file.read_faces(os.path.join(mesh, "faces"))
    boundary = foam_file.read_boundary(os.path.join(mesh, "boundary"))
    start = int(boundary[patch]["startFace"])
    n = int(boundary[patch]["nFaces"])
    return np.array([points[flat[offsets[f]:offsets[f + 1]]].mean(axis=0)
                     for f in range(start, start + n)])


def write_dict(run: str, rel: str, obj: str, body: str) -> None:
    """Write dictionary ``rel`` (relative to the case) with a header."""
    foam_file.write_foam_file(os.path.join(run, rel), "dictionary", obj, body)


# Boundary conditions: patch -> entries; "BED" is the floor (IB: below the
# immersed bed) or the bed wall (BF)
FIELDS = {
    "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (%(u_bulk)g 0 0)", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant (0 0 0); value uniform (0 0 0)",
        "outlet": "type zeroGradient",
        "top": "type slip",
        "ib": "type zeroGradient",
        "bf": "type noSlip",
    }),
    "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        "inlet": "type zeroGradient",
        "outlet": "type fixedValue; value uniform 0",
        "top": "type zeroGradient",
        "ib": "type zeroGradient",
        "bf": "type zeroGradient",
    }),
    "k": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 1e-3", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant 0; value uniform 1e-3",
        "outlet": "type zeroGradient",
        "top": "type zeroGradient",
        "ib": "type zeroGradient",
        "bf": "type kqRWallFunction; value uniform 1e-3",
    }),
    "omega": ("volScalarField", "[0 0 -1 0 0 0 0]", "uniform 10", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant 0; value uniform 10",
        "outlet": "type zeroGradient",
        "top": "type zeroGradient",
        "ib": "type zeroGradient",
        # wall-cell omega (and k) fixed by wallLawTurbulence
        "bf": "type zeroGradient",
    }),
    "nut": ("volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        "inlet": "type calculated; value uniform 0",
        "outlet": "type calculated; value uniform 0",
        "top": "type calculated; value uniform 0",
        "ib": "type calculated; value uniform 0",
        "bf": "NUT_WALL_LAW",
    }),
}

NUT_WALL_LAW = """\
        type            nutWallLaw;
        wallLaw
        {
            type        vanDriest;
            ks          %(ks)g;
        }
        value           uniform 0;"""

SEDIMENT = """\
        grainDiameter   %(d50)g;
        specificGravity 2.63;
        porosity        0.4;
        gravity         9.81;
        criticalShields 0.05;
        frictionCoefficient 0.625;
        bedLoad         engelundFredsoe;
        openBoundaries  (inlet outlet);
        fixedBed        { inlet 0.15; outlet 0.1; }
        morphStartTime  1e9;        // frozen bed
        morphFactor     1;
        sandSlide       false;"""


TURBULENCE_FIELDS = ("k", "omega", "nut")


def laminar_profile(z: np.ndarray, h: float, u_bulk: float) -> np.ndarray:
    """Half-Poiseuille velocity (no-slip at z = 0, slip at z = h), mean u_bulk."""
    eta = np.clip(z / h, 0.0, 1.0)
    return 1.5 * u_bulk * (2.0 * eta - eta ** 2)


def field_file(name: str, s: dict) -> str:
    """Return the body of the initial field file ``name``."""
    _, dims, internal, bcs = FIELDS[name]
    floor = "bottom" if s["variant"] == "ib" else "bed"
    lines = [f"dimensions      {dims};", "",
             f"internalField   {internal % s};", "", "boundaryField", "{"]
    for patch in ("inlet", "outlet", "top", s["variant"]):
        spec = bcs[patch]
        lines += [f"    {floor if patch == s['variant'] else patch}", "    {"]
        if spec == "NUT_WALL_LAW":
            lines += [NUT_WALL_LAW % s]
        else:
            lines += [f"        {e.strip()};" for e in spec.split(";")
                      if e.strip()]
        lines += ["    }"]
    lines += ["    frontAndBack", "    {", "        type            symmetry;",
              "    }", "}"]
    return "\n".join(lines)


def cmd_setup(args: argparse.Namespace) -> None:
    """Create the run directory (before meshing)."""
    s = dict(SETUP)
    s.update({"variant": args.variant, "slope_deg": args.slope_deg,
              "pit_depth": args.pit_depth, "dx": args.dx, "dz": args.dz,
              "end_time": args.end_time, "nu": args.nu,
              "u_bulk": args.u_bulk, "ks": args.ks,
              "delta_t": args.delta_t, "write_interval": args.write_interval,
              "laminar": args.laminar})
    run = args.run_dir
    os.makedirs(run, exist_ok=True)

    nx = int(round(s["length"] / s["dx"]))
    if s["variant"] == "ib":
        z0 = -(s["pit_depth"] + s["below"] * s["dz"])
        nz = int(round((s["depth"] - z0) / s["dz"]))
        floor = "bottom"
    else:
        z0 = 0.0
        nz = int(round(s["depth"] / s["dz"]))
        floor = "bed"
    s.update({"nx": nx, "nz": nz, "z0": z0})
    case_setup.save_setup(run, s)

    write_dict(run, "system/blockMeshDict", "blockMeshDict", f"""\
scale   1;

vertices
(
    (0 0 {z0:g}) ({s['length']:g} 0 {z0:g}) ({s['length']:g} {s['dx']:g} {z0:g})
    (0 {s['dx']:g} {z0:g})
    (0 0 {s['depth']:g}) ({s['length']:g} 0 {s['depth']:g})
    ({s['length']:g} {s['dx']:g} {s['depth']:g}) (0 {s['dx']:g} {s['depth']:g})
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
    {floor} {{ type {'patch' if s['variant'] == 'ib' else 'wall'};
             faces ((0 3 2 1)); }}
    frontAndBack {{ type symmetry; faces ((0 1 5 4) (3 7 6 2)); }}
);
""")

    common = f"""\
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {s['end_time']:g};
deltaT          {s['delta_t']:g};
writeControl    adjustableRunTime;
writeInterval   {s['write_interval']:g};
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
"""
    if s["variant"] == "ib":
        write_dict(run, "system/controlDict", "controlDict",
                   'libs            ("libIBMethod.so");\n'
                   "application     ibScourPimpleFoam;\n" + common)
        write_dict(run, "constant/immersedBoundaryDict",
                   "immersedBoundaryDict", f"""\
bed
{{
    type            heightField;
    floorPatch      bottom;
    initialElevation 0;     // overridden by 0/finite-area/zb
{SEDIMENT % s}
    smoothShear     1;          // bed-shear filter (default; 0 to disable)
    maxBedChange    1e-4;
    updateDistance  5e-4;
}}

cutFaces        true;       // apertures for the height-field bed

solids
{{
}}

turbulence
{{
    kTreatment      equilibrium;
}}

wallModel
{{
    mode            stress;
    reportMaxShear  true;

    wallLaw
    {{
        type            vanDriest;
        ks              {s['ks']:g};
    }}
}}
""")
        fa_dir = "system/finite-area"
    else:
        write_dict(run, "system/controlDict", "controlDict",
                   "libs            (ibScour);\n"
                   "application     pimpleFoam;\n" + common + f"""
functions
{{
    wallShearStress
    {{
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (bed);
        writeControl    writeTime;
    }}
    bed
    {{
        type            deformingBed;
        libs            (ibScour);
        bedPatch        bed;
        moveMesh        false;
{SEDIMENT % s}
    }}
}}
""")
        if not s["laminar"]:
            write_dict(run, "constant/fvOptions", "fvOptions", f"""\
wallTurbulence
{{
    type            wallLawTurbulence;
    patches         (bed);
    wallLaw
    {{
        type        vanDriest;
        ks          {s['ks']:g};
    }}
}}
""")
        fa_dir = "system/bedRegion/finite-area"

    write_dict(run, f"{fa_dir}/faSchemes", "faSchemes", """\
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
""")
    write_dict(run, f"{fa_dir}/faSolution", "faSolution", "solvers {}\n")

    write_dict(run, "system/fvSchemes", "fvSchemes", """\
ddtSchemes      { default backward; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         Gauss upwind;
    div(phi,U)      Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
wallDist        { method meshWave; }
""")
    write_dict(run, "system/fvSolution", "fvSolution", f"""\
solvers
{{
    p
    {{
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       0;
        relTol          0.01;
    }}
    pFinal
    {{
        $p;
        tolerance       1e-6;
        relTol          0;
    }}
    "(U|k|omega)"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-5;
        relTol          0.1;
    }}
    "(U|k|omega)Final"
    {{
        $U;
        relTol          0;
    }}
}}

PIMPLE
{{
    nOuterCorrectors 1;
    nCorrectors     3;
    nNonOrthogonalCorrectors {0 if s['variant'] == 'ib' else 1};
}}
""")
    write_dict(run, "system/decomposeParDict", "decomposeParDict", """\
numberOfSubdomains 4;
method          simple;
coeffs
{
    n           (4 1 1);
}
""")
    write_dict(run, "constant/transportProperties", "transportProperties",
               f"transportModel  Newtonian;\nnu              {s['nu']:g};\n")
    write_dict(run, "constant/turbulenceProperties", "turbulenceProperties",
               "simulationType  laminar;\n" if s["laminar"] else
               "simulationType  RAS;\nRAS\n{\n    RASModel        kOmegaSST;"
               "\n    turbulence      on;\n    printCoeffs     on;\n}\n")
    LOG.info("setup %s (%s, slope %g deg): %d x %d cells", run, s["variant"],
             s["slope_deg"], nx, nz)


def write_zb(case: str, mesh: str, s: dict) -> np.ndarray:
    """Write ``<case>/0/finite-area/zb`` on the floor faces of ``mesh``.

    The faMesh of the height-field bed is built at run time from the floor
    patch (face order = patch order), so decomposePar does not decompose
    this field: ``decompose-zb`` writes it for every processor.
    """
    floor = mesh_patch_centres(mesh, "bottom")
    zb = bed_elevation(floor[:, 0], s) if len(floor) else np.zeros(0)
    foam_file.write_foam_file(
        os.path.join(case, "0", "finite-area", "zb"), "areaScalarField",
        "zb", "dimensions      [0 1 0 0 0 0 0];\n\ninternalField   "
        "nonuniform List<scalar> " + foam_file.format_scalar_list(zb)
        + ";\n\nboundaryField\n{\n    \".*\"\n    {\n        type"
        "            zeroGradient;\n    }\n    \"procBoundary.*\"\n    {"
        "\n        type            processor;\n    }\n}\n")
    return zb


def cmd_decompose_zb(args: argparse.Namespace) -> None:
    """Write the initial z_b of an IB run for every processor directory."""
    run = args.run_dir
    s = case_setup.load_setup(run)
    procs = sorted(d for d in os.listdir(run) if d.startswith("processor"))
    for p in procs:
        case = os.path.join(run, p)
        write_zb(case, os.path.join(case, "constant", "polyMesh"), s)
    LOG.info("decompose-zb %s: %d processors", run, len(procs))


def cmd_init(args: argparse.Namespace) -> None:
    """Write the initial fields and the inflow; impose the bed."""
    run = args.run_dir
    s = case_setup.load_setup(run)
    laminar = s.get("laminar", False)
    for name, (cls, *_rest) in FIELDS.items():
        if laminar and name in TURBULENCE_FIELDS:
            continue
        foam_file.write_foam_file(os.path.join(run, "0", name), cls, name,
                                  field_file(name, s))

    inlet = patch_centres(run, "inlet")
    if laminar:
        u = laminar_profile(inlet[:, 2], s["depth"], s["u_bulk"])
        inflow = {}
    else:
        u, k, omega = prof.turbulent_profile(inlet[:, 2], s["depth"],
                                             s["u_bulk"], s["nu"])
        inflow = {"k": k, "omega": omega}
    vel = np.zeros((len(u), 3))
    vel[:, 0] = u
    case_setup.write_boundary_data(run, "inlet", inlet, {"U": vel, **inflow})

    if s["variant"] == "ib":
        zb = write_zb(run, os.path.join(run, "constant", "polyMesh"), s)
    else:
        mesh = os.path.join(run, "constant", "polyMesh")
        points = foam_file.read_points(os.path.join(mesh, "points"))
        zb = bed_elevation(points[:, 0], s)
        h = s["depth"]
        points[:, 2] = zb + points[:, 2] * (h - zb) / h
        foam_file.write_foam_file(os.path.join(mesh, "points"),
                                  "vectorField", "points",
                                  foam_file.format_vector_list(points))
    LOG.info("init %s: z_b %.4g to %.4g m", run, zb.min(), zb.max())


def main() -> None:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("setup")
    st.add_argument("--run-dir", required=True)
    st.add_argument("--variant", choices=("ib", "bf"), required=True)
    st.add_argument("--slope-deg", type=float, default=30.0)
    st.add_argument("--pit-depth", type=float, default=SETUP["pit_depth"])
    st.add_argument("--dx", type=float, default=SETUP["dx"])
    st.add_argument("--dz", type=float, default=SETUP["dz"])
    st.add_argument("--end-time", type=float, default=SETUP["end_time"])
    st.add_argument("--laminar", action="store_true",
                    help="laminar flow and half-Poiseuille inflow")
    st.add_argument("--u-bulk", type=float, default=SETUP["u_bulk"])
    st.add_argument("--nu", type=float, default=NU)
    st.add_argument("--ks", type=float, default=SETUP["ks"])
    st.add_argument("--delta-t", type=float, default=SETUP["delta_t"])
    st.add_argument("--write-interval", type=float,
                    default=SETUP["write_interval"])
    it = sub.add_parser("init")
    it.add_argument("--run-dir", required=True)
    dz = sub.add_parser("decompose-zb",
                        help="after decomposePar (IB runs in parallel)")
    dz.add_argument("--run-dir", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    {"setup": cmd_setup, "init": cmd_init,
     "decompose-zb": cmd_decompose_zb}[args.cmd](args)


if __name__ == "__main__":
    main()
