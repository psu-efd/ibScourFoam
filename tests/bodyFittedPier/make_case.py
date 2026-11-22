#!/usr/bin/env python3
"""Generate a body-fitted Roulund (2005) pier case on a blockMesh O-grid.

Domain and inflow as the immersed-boundary Roulund case (pier D = 0.1 m at
the origin, x in [-0.5, 0.9] m, y in [-0.45, 0.45] m, bed at z = 0, rigid
lid at z = 0.2 m, inflow profile from ``constant/boundaryData/inlet`` of a
reference case). Four ring blocks around the pier reach a square of
half-width ``a``; eight Cartesian blocks fill the rest. With ``--r-circle
rc`` (rc > D/2) the ring is split: four circular ring blocks from the pier
to the circle of radius rc (both edges arcs, so the grid lines stay
circles where the scour hole forms) and four blocks from that circle to
the square (the circle-to-square transition outside the hole). All blocks are
extruded vertically, so the mesh columns are vertical and a bed shape is
imposed by stretching each column between the bed and the lid. The mesh is
clustered towards the pier and the bed (first cells ``--dr`` and ``--dz``).

The bed (patch ``bed``) uses the wall law of the immersed boundary
(``nutWallLaw`` + ``wallLawTurbulence``, vanDriest, ks = 0.7 mm, as the IB
runs); the pier uses the OpenFOAM wall functions of the IB runs.

Modes:

* ``frozen``: static mesh, optionally deformed to a bed surface
  (``--bed-vtp``: z_b interpolated from the face centres of a ``bed.vtp``
  of an IB run); stock pimpleFoam with wallShearStress output;
* ``live``: deforming mesh (``deformingBed``) with the sediment settings
  of the IB run sa_mf2 and an optional bed-load adaptation length
  (``--adaptation-length``). The mesh moves with displacementLayeredMotion:
  each vertical column stretches linearly between its bed point and the
  lid (cell zone ``allCells`` between the face zones ``bedZone`` and
  ``topZone``, written by ``init``). A Laplacian displacement solve with
  cell-to-point interpolation averages the bed displacement over
  neighbouring columns and collapsed the first cells on deposition mounds
  where the bed changes by more than the first-cell height between
  columns.

Usage::

    make_case.py setup --run-dir DIR --inflow REF_CASE [--mode frozen]
                       [--bed-vtp bed.vtp] [--n-theta 128] ...
    blockMesh -case DIR
    make_case.py init --run-dir DIR
"""

import argparse
import logging
import math
import os
import shutil
import sys
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import case_setup  # noqa: E402
import foam_file  # noqa: E402

LOG = logging.getLogger("bodyFittedPier.make_case")

SETUP: Dict[str, Any] = {
    "D": 0.1,          # pier diameter [m]
    "a": 0.15,         # half-width of the O-grid square [m]
    "x": (-0.5, 0.9),  # streamwise extent [m]
    "y": (-0.45, 0.45),  # spanwise extent [m]
    "depth": 0.2,      # rigid lid above the bed [m]
    "nu": 1e-6,
    "ks": 7e-4,        # bed roughness (Roulund et al. 2005, Fig. 33) [m]
    "u_bulk": 0.46,
}


def expansion(length: float, n: int, first: float) -> float:
    """Return the blockMesh expansion ratio (last/first cell) of n cells
    over ``length`` with the first cell ``first``."""
    if first * n >= length:
        return 1.0
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        if first * (r ** n - 1) / (r - 1) > length:
            hi = r
        else:
            lo = r
    r = 0.5 * (lo + hi)
    return r ** (n - 1)


def block_mesh_dict(s: Dict) -> str:
    """Return the blockMeshDict text of the O-grid."""
    R = 0.5 * s["D"]
    a = s["a"]
    xs = [s["x"][0], -a, a, s["x"][1]]
    ys = [s["y"][0], -a, a, s["y"][1]]
    zs = [0.0, s["depth"]]
    nq = s["n_theta"] // 4

    verts: List[Tuple[float, float, float]] = []
    index: Dict[Tuple, int] = {}

    def add(key: Tuple, p: Tuple[float, float, float]) -> None:
        index[key] = len(verts)
        verts.append(p)

    rc = s.get("r_circle", 0.0)
    split = rc > R
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                add(("g", i, j, k), (x, y, z))
        for q, ang in enumerate((225, 315, 45, 135)):
            t = math.radians(ang)
            add(("c", q, k), (R * math.cos(t), R * math.sin(t), z))
            if split:
                add(("o", q, k), (rc * math.cos(t), rc * math.sin(t), z))

    # Square corners (225, 315, 45, 135 deg) = grid points
    sq = [(1, 1), (2, 1), (2, 2), (1, 2)]

    def s_(q: int, k: int) -> int:
        return index[("g", sq[q][0], sq[q][1], k)]

    def c_(q: int, k: int) -> int:
        return index[("c", q, k)]

    def o_(q: int, k: int) -> int:
        return index[("o", q, k)]

    # Grading
    gz = expansion(s["depth"], s["nz"], s["dz"])
    if split:
        gr = expansion(rc - R, s["n_r"], s["dr"])
        # outer ring: continue from the last inner cell
        r_in = gr**(1.0/max(s["n_r"] - 1, 1))
        last = s["dr"]*r_in**(s["n_r"] - 1)
        gro = expansion(a*math.sqrt(2) - rc, s["n_r_outer"], last)
    else:
        gr = expansion(a - R, s["n_r"], s["dr"])
    d_sq = 2 * a / nq                      # cell size along the square
    gx_up = 1.0 / expansion(xs[1] - xs[0], s["nx_up"], d_sq)
    gx_dn = expansion(xs[3] - xs[2], s["nx_down"], d_sq)
    gy_lo = 1.0 / expansion(ys[1] - ys[0], s["ny_side"], d_sq)
    gy_hi = expansion(ys[3] - a, s["ny_side"], d_sq)

    blocks = []
    # Ring blocks (c_{q+1} c_q s_q s_{q+1}): x1 along the arc, x2 radial
    for q in range(4):
        q1 = (q + 1) % 4
        if split:
            v = [c_(q1, 0), c_(q, 0), o_(q, 0), o_(q1, 0),
                 c_(q1, 1), c_(q, 1), o_(q, 1), o_(q1, 1)]
            blocks.append(
                "hex ({}) ({} {} {}) simpleGrading (1 {:.6g} {:.6g})"
                .format(" ".join(map(str, v)), nq, s["n_r"], s["nz"],
                        gr, gz))
            v = [o_(q1, 0), o_(q, 0), s_(q, 0), s_(q1, 0),
                 o_(q1, 1), o_(q, 1), s_(q, 1), s_(q1, 1)]
            blocks.append(
                "hex ({}) ({} {} {}) simpleGrading (1 {:.6g} {:.6g})"
                .format(" ".join(map(str, v)), nq, s["n_r_outer"],
                        s["nz"], gro, gz))
            continue
        v = [c_(q1, 0), c_(q, 0), s_(q, 0), s_(q1, 0),
             c_(q1, 1), c_(q, 1), s_(q, 1), s_(q1, 1)]
        blocks.append("hex ({}) ({} {} {}) simpleGrading (1 {:.6g} {:.6g})"
                      .format(" ".join(map(str, v)), nq, s["n_r"], s["nz"],
                              gr, gz))
    # Cartesian blocks
    nxs = [s["nx_up"], nq, s["nx_down"]]
    nys = [s["ny_side"], nq, s["ny_side"]]
    gxs = [gx_up, 1, gx_dn]
    gys = [gy_lo, 1, gy_hi]
    for j in range(3):
        for i in range(3):
            if i == 1 and j == 1:
                continue
            v = [index[("g", i, j, 0)], index[("g", i + 1, j, 0)],
                 index[("g", i + 1, j + 1, 0)], index[("g", i, j + 1, 0)]]
            v += [index[("g", ii, jj, 1)] for ii, jj in
                  ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))]
            blocks.append(
                "hex ({}) ({} {} {}) simpleGrading ({:.6g} {:.6g} {:.6g})"
                .format(" ".join(map(str, v)), nxs[i], nys[j], s["nz"],
                        gxs[i], gys[j], gz))

    edges = []
    for q in range(4):
        q1 = (q + 1) % 4
        mid = math.radians((225, 315, 45, 135)[q] + 45)
        for k, z in enumerate(zs):
            edges.append("arc {} {} ({:.9g} {:.9g} {:.9g})".format(
                c_(q, k), c_(q1, k), R * math.cos(mid), R * math.sin(mid), z))
            if split:
                edges.append("arc {} {} ({:.9g} {:.9g} {:.9g})".format(
                    o_(q, k), o_(q1, k), rc * math.cos(mid),
                    rc * math.sin(mid), z))

    def face(vs: List[int]) -> str:
        return "({})".format(" ".join(map(str, vs)))

    pier = [face([c_((q + 1) % 4, 0), c_(q, 0), c_(q, 1), c_((q + 1) % 4, 1)])
            for q in range(4)]
    inlet = [face([index[("g", 0, j, 0)], index[("g", 0, j, 1)],
                   index[("g", 0, j + 1, 1)], index[("g", 0, j + 1, 0)]])
             for j in range(3)]
    outlet = [face([index[("g", 3, j, 0)], index[("g", 3, j + 1, 0)],
                    index[("g", 3, j + 1, 1)], index[("g", 3, j, 1)]])
              for j in range(3)]
    sides = [face([index[("g", i, 0, 0)], index[("g", i + 1, 0, 0)],
                   index[("g", i + 1, 0, 1)], index[("g", i, 0, 1)]])
             for i in range(3)]
    sides += [face([index[("g", i, 3, 0)], index[("g", i, 3, 1)],
                    index[("g", i + 1, 3, 1)], index[("g", i + 1, 3, 0)]])
              for i in range(3)]

    vtxt = "\n".join("    ({:.9g} {:.9g} {:.9g})".format(*p) for p in verts)
    return """\
scale   1;

vertices
(
{v}
);

blocks
(
    {b}
);

edges
(
    {e}
);

// bed, top: all block faces at z = 0 and z = depth (defaultPatch not
// used; listed by blockMesh from the hex faces below)
boundary
(
    pier   {{ type wall;  faces ({pier}); }}
    inlet  {{ type patch; faces ({inlet}); }}
    outlet {{ type patch; faces ({outlet}); }}
    frontAndBack {{ type patch; faces ({sides}); }}
    bed    {{ type wall;  faces ({bed}); }}
    top    {{ type patch; faces ({top}); }}
);
""".format(v=vtxt, b="\n    ".join(blocks), e="\n    ".join(edges),
           pier=" ".join(pier), inlet=" ".join(inlet),
           outlet=" ".join(outlet), sides=" ".join(sides),
           bed=" ".join(bed_faces(blocks, 0)),
           top=" ".join(bed_faces(blocks, 1)))


def bed_faces(blocks: List[str], level: int) -> List[str]:
    """Bottom (level 0) or top (level 1) faces of all hex blocks."""
    out = []
    for b in blocks:
        v = [int(t) for t in b[b.index("(") + 1:b.index(")")].split()]
        if level == 0:
            out.append("({} {} {} {})".format(v[0], v[3], v[2], v[1]))
        else:
            out.append("({} {} {} {})".format(v[4], v[5], v[6], v[7]))
    return out


# Boundary conditions per field: patch -> entries (";"-separated) or a raw
# block (starting with a newline)
BED_NUT = """
        type            nutWallLaw;
        wallLaw
        {
            type        vanDriest;
            ks          %(ks)g;
        }
        value           uniform 0;"""

FIELDS = {
    "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (%(u_bulk)g 0 0)", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant (0 0 0); value uniform (0 0 0)",
        "outlet": "type zeroGradient",
        "top": "type slip",
        "frontAndBack": "type slip",
        "pier": "type noSlip",
        "bed": "%(bedU)s",
    }),
    "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        "inlet": "type zeroGradient",
        "outlet": "type fixedValue; value uniform 0",
        "top": "type zeroGradient",
        "frontAndBack": "type zeroGradient",
        "pier": "type zeroGradient",
        "bed": "type zeroGradient",
    }),
    "k": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 1e-3", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant 0; value uniform 1e-3",
        "outlet": "type zeroGradient",
        "top": "type zeroGradient",
        "frontAndBack": "type zeroGradient",
        "pier": "type kqRWallFunction; value uniform 1e-3",
        "bed": "type kqRWallFunction; value uniform 1e-3",
    }),
    "omega": ("volScalarField", "[0 0 -1 0 0 0 0]", "uniform 10", {
        "inlet": "type timeVaryingMappedFixedValue; mapMethod nearest;"
                 " offset constant 0; value uniform 10",
        "outlet": "type zeroGradient",
        "top": "type zeroGradient",
        "frontAndBack": "type zeroGradient",
        "pier": "type omegaWallFunction; blending binomial; n 2;"
                " value uniform 10",
        # Wall-cell omega (and k) fixed by wallLawTurbulence
        "bed": "type zeroGradient",
    }),
    "nut": ("volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        "inlet": "type calculated; value uniform 0",
        "outlet": "type calculated; value uniform 0",
        "top": "type calculated; value uniform 0",
        "frontAndBack": "type calculated; value uniform 0",
        "pier": "type nutUWallFunction; blending stepwise; value uniform 0",
        "bed": BED_NUT,
    }),
    # Point displacement (displacementLayeredMotion): bed set by
    # deformingBed (vertical), lid fixed; the layered solver sets all
    # other points from the columns
    "pointDisplacement": ("pointVectorField", "[0 1 0 0 0 0 0]",
                          "uniform (0 0 0)", {
        "inlet": "type zeroGradient",
        "outlet": "type zeroGradient",
        "top": "type fixedValue; value uniform (0 0 0)",
        "frontAndBack": "type zeroGradient",
        "pier": "type zeroGradient",
        "bed": "type fixedValue; value uniform (0 0 0)",
    }),
}


def field_file(name: str, s: Dict) -> str:
    """Return the body of the initial field file ``name``."""
    _, dims, internal, bcs = FIELDS[name]
    lines = ["dimensions      {};".format(dims), "",
             "internalField   {};".format(internal % s), "",
             "boundaryField", "{"]
    for patch, spec in bcs.items():
        spec = spec % s
        lines += ["    {}".format(patch), "    {"]
        if spec.startswith("\n"):
            lines += [spec.strip("\n")]
        else:
            lines += ["        {};".format(e.strip())
                      for e in spec.split(";") if e.strip()]
        lines += ["    }"]
    lines += ["}"]
    return "\n".join(lines)


def write_dict(run: str, rel: str, obj: str, body: str) -> None:
    """Write dictionary ``rel`` (relative to the case) with header."""
    foam_file.write_foam_file(os.path.join(run, rel), "dictionary", obj,
                              body)


SEDIMENT = """\
        grainDiameter   0.00026;
        specificGravity 2.63;
        porosity        0.4;
        gravity         9.81;
        criticalShields %(theta_c0)g;
        frictionCoefficient 0.625;
        bedLoad         engelundFredsoe;
        bedSlope
        {
            model       kochFlokstra;
            alpha       1.5;
        }
        morphStartTime  %(morph_start)g;
        morphFactor     %(morph_factor)g;
        maxBedChange    1e-4;
        openBoundaries  (inlet outlet);
        fixedBed
        {
            inlet       0.05;
        }
        sandSlide       true;
        slideDiffusivity 5e-4;
        slideTolerance  0.02;
        slideMaxIter    500;
        adaptationLength %(adaptation_length)g;"""


def cmd_setup(args: argparse.Namespace) -> None:
    """Create the run directory (before meshing)."""
    s = dict(SETUP)
    s.update({
        "n_theta": args.n_theta, "n_r": args.n_r, "nz": args.nz,
        "dr": args.dr, "dz": args.dz, "nx_up": args.nx_up,
        "a": args.a, "r_circle": args.r_circle,
        "n_r_outer": args.n_r_outer,
        "nx_down": args.nx_down, "ny_side": args.ny_side,
        "mode": args.mode, "bed_vtp": os.path.abspath(args.bed_vtp)
        if args.bed_vtp else None,
        "inflow": os.path.abspath(args.inflow),
        "end_time": args.end_time, "theta_c0": args.theta_c0,
        "morph_factor": args.morph_factor, "morph_start": args.morph_start,
        "adaptation_length": args.adaptation_length,
    })
    s["bedU"] = ("type movingWallVelocity; value uniform (0 0 0)"
                 if args.mode == "live" else "type noSlip")
    run = args.run_dir
    os.makedirs(run, exist_ok=True)
    case_setup.save_setup(run, s)

    write_dict(run, "system/blockMeshDict", "blockMeshDict",
               block_mesh_dict(s))

    functions = """\
    wallShearStress
    {
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (bed pier);
        writeControl    writeTime;
    }
    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
    }
"""
    if args.mode == "live":
        functions += """\
    bed
    {
        type            deformingBed;
        libs            (ibScour);
        bedPatch        bed;
%s
    }
""" % (SEDIMENT % s)

    write_dict(run, "system/controlDict", "controlDict", """\
application     pimpleFoam;
libs            (ibScour);
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         %(end)g;
deltaT          1e-3;
adjustTimeStep  yes;
maxCo           0.9;
maxDeltaT       0.002;
writeControl    adjustable;
writeInterval   %(wi)g;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;

functions
{
%(fun)s}
""" % {"end": args.end_time, "wi": args.write_interval, "fun": functions})

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

    write_dict(run, "system/fvSolution", "fvSolution", """\
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-6;
        relTol          0.01;
    }
    pFinal
    {
        $p;
        relTol          0;
    }
    "pcorr.*"
    {
        $p;
        tolerance       1e-5;
        relTol          0;
    }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }
    "(U|k|omega)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors     3;
    nNonOrthogonalCorrectors 1;
    // The bed moves slowly and vertically: the mesh fluxes of the motion
    // keep continuity; correctPhi would solve pcorr (hundreds of GAMG
    // iterations) every step
    correctPhi      no;
    moveMeshOuterCorrectors no;
    pRefCell        0;
    pRefValue       0;
}
""")
    write_dict(run, "system/decomposeParDict", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod          scotch;\n"
               % args.np)

    write_dict(run, "constant/fvOptions", "fvOptions", """\
wallTurbulence
{
    type            wallLawTurbulence;
    patches         (bed);
    wallLaw
    {
        type        vanDriest;
        ks          %(ks)g;
    }
}
""" % s)
    write_dict(run, "constant/transportProperties", "transportProperties",
               "transportModel  Newtonian;\nnu              %g;\n" % s["nu"])
    write_dict(run, "constant/turbulenceProperties", "turbulenceProperties",
               "simulationType  RAS;\nRAS\n{\n    RASModel        kOmegaSST;"
               "\n    turbulence      on;\n    printCoeffs     on;\n}\n")

    if args.mode == "live":
        write_dict(run, "constant/dynamicMeshDict", "dynamicMeshDict", """\
dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolverLibs (fvMotionSolvers);
motionSolver    displacementLayeredMotion;

// Every vertical column stretches linearly between its bed point and the
// fixed lid (topological walk between the face zones): no lateral smoothing
// of the bed displacement, so the first cells cannot collapse where the bed
// changes by more than the first-cell height between neighbouring columns
displacementLayeredMotionCoeffs
{
    regions
    {
        allCells
        {
            interpolationScheme linear;
            boundaryField
            {
                bedZone
                {
                    type    follow;     // pointDisplacement of deformingBed
                }
                topZone
                {
                    type    fixedValue;
                    value   uniform (0 0 0);
                }
            }
        }
    }
}
""")
        for f in ("faSchemes", "faSolution"):
            body = ("ddtSchemes { default Euler; }\n"
                    "gradSchemes { default Gauss linear; }\n"
                    "divSchemes { default Gauss linear; }\n"
                    "laplacianSchemes { default Gauss linear corrected; }\n"
                    "interpolationSchemes { default linear; }\n"
                    "snGradSchemes { default corrected; }\n"
                    if f == "faSchemes" else "solvers {}\n")
            write_dict(run, "system/bedRegion/finite-area/" + f, f, body)

    # Inflow profile of the reference case
    src = os.path.join(args.inflow, "constant", "boundaryData", "inlet")
    dst = os.path.join(run, "constant", "boundaryData", "inlet")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

    nc = (4 * (s["n_theta"] // 4) * s["n_r"]
          + (s["nx_up"] + s["n_theta"] // 4 + s["nx_down"])
          * (2 * s["ny_side"] + s["n_theta"] // 4)
          - (s["n_theta"] // 4) ** 2) * s["nz"]
    LOG.info("setup %s (%s): about %d cells", run, args.mode, nc)


def read_bed_vtp(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return face-centre (x, y) and z_b of a bed.vtp surface."""
    import vtk  # noqa: E402
    from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    pd = reader.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(pd)
    cc.Update()
    xyz = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    return xyz[:, :2].astype(float), \
        vtk_to_numpy(pd.GetCellData().GetArray("zb")).astype(float)


def write_motion_zones(mesh_dir: str) -> None:
    """Write the zones of the layered mesh motion (before decomposePar).

    cellZone ``allCells`` (all cells) between the faceZones ``bedZone``
    (faces of patch ``bed``) and ``topZone`` (faces of patch ``top``).
    """
    boundary = foam_file.read_boundary(os.path.join(mesh_dir, "boundary"))
    with open(os.path.join(mesh_dir, "owner"), errors="ignore") as f:
        head = f.read(4000)
    n_cells = int(head.split("nCells:")[1].split()[0])

    def labels(start: int, n: int) -> str:
        return "List<label> %d\n(\n%s\n)" % (
            n, "\n".join(str(i) for i in range(start, start + n)))

    zones = []
    for zone, patch in (("bedZone", "bed"), ("topZone", "top")):
        start = int(boundary[patch]["startFace"])
        n = int(boundary[patch]["nFaces"])
        zones.append("%s\n{\n    type faceZone;\n    faceLabels %s;\n"
                     "    flipMap List<bool> %d{0};\n}\n"
                     % (zone, labels(start, n), n))
    foam_file.write_foam_file(os.path.join(mesh_dir, "faceZones"),
                              "regIOobject", "faceZones",
                              "2\n(\n%s)\n" % "".join(zones))
    foam_file.write_foam_file(
        os.path.join(mesh_dir, "cellZones"), "regIOobject", "cellZones",
        "1\n(\nallCells\n{\n    type cellZone;\n    cellLabels %s;\n}\n)\n"
        % labels(0, n_cells))
    LOG.info("motion zones: allCells %d cells, bedZone %d, topZone %d faces",
             n_cells, boundary["bed"]["nFaces"], boundary["top"]["nFaces"])


def cmd_init(args: argparse.Namespace) -> None:
    """Write the initial fields; impose the bed surface on the mesh."""
    run = args.run_dir
    s = case_setup.load_setup(run)
    names = [n for n in FIELDS if n != "pointDisplacement"
             or s["mode"] == "live"]
    for name in names:
        foam_file.write_foam_file(os.path.join(run, "0", name),
                                  FIELDS[name][0], name, field_file(name, s))
    if s["mode"] == "live":
        write_motion_zones(os.path.join(run, "constant", "polyMesh"))

    if s.get("bed_vtp"):
        from scipy.interpolate import griddata  # noqa: E402
        xy, zb = read_bed_vtp(s["bed_vtp"])
        pfile = os.path.join(run, "constant", "polyMesh", "points")
        pts = foam_file.read_points(pfile)
        zp = griddata(xy, zb, pts[:, :2], method="linear")
        miss = np.isnan(zp)
        if miss.any():
            zp[miss] = griddata(xy, zb, pts[miss, :2], method="nearest")
        h = s["depth"]
        pts[:, 2] = zp + pts[:, 2] * (h - zp) / h
        foam_file.write_foam_file(pfile, "vectorField", "points",
                                  foam_file.format_vector_list(pts))
        LOG.info("bed from %s: z_b %.4g to %.4g m (%d nearest)",
                 s["bed_vtp"], zp.min(), zp.max(), miss.sum())


def main() -> None:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd")
    st = sub.add_parser("setup")
    st.add_argument("--run-dir", required=True)
    st.add_argument("--inflow", required=True,
                    help="case with constant/boundaryData/inlet")
    st.add_argument("--mode", choices=("frozen", "live"), default="frozen")
    st.add_argument("--bed-vtp", default=None)
    st.add_argument("--n-theta", type=int, default=128)
    st.add_argument("--n-r", type=int, default=48)
    st.add_argument("--nz", type=int, default=60)
    st.add_argument("--dr", type=float, default=5e-4)
    st.add_argument("--dz", type=float, default=5e-4)
    st.add_argument("--nx-up", type=int, default=30)
    st.add_argument("--nx-down", type=int, default=50)
    st.add_argument("--ny-side", type=int, default=24)
    st.add_argument("--a", type=float, default=0.15,
                    help="half-width of the square around the ring [m]")
    st.add_argument("--r-circle", type=float, default=0.0,
                    help="radius [m] to which the ring stays circular "
                         "(0: ring from the pier to the square)")
    st.add_argument("--n-r-outer", type=int, default=8,
                    help="radial cells from the circle to the square")
    st.add_argument("--end-time", type=float, default=10.0)
    st.add_argument("--write-interval", type=float, default=1.0)
    st.add_argument("--theta-c0", type=float, default=0.05)
    st.add_argument("--morph-factor", type=float, default=2.0)
    st.add_argument("--morph-start", type=float, default=5.0)
    st.add_argument("--np", type=int, default=64)
    st.add_argument("--adaptation-length", type=float, default=2.6e-3,
                    help="bed-load adaptation (saturation) length [m] "
                         "(live mode; default 10 d50, 0 = local "
                         "equilibrium)")
    it = sub.add_parser("init")
    it.add_argument("--run-dir", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if args.cmd is None:
        ap.error("command required")
    {"setup": cmd_setup, "init": cmd_init}[args.cmd](args)


if __name__ == "__main__":
    main()
