#!/usr/bin/env python3
"""Convert a decomposed live-bed pier case to layered (column) mesh motion.

displacementLayeredMotion interpolates the point displacement linearly along
the mesh columns between the bed and the lid face zones (topological walk),
so every column stretches with its own bed point. The script

  - writes cellZones (allCells) and faceZones (bedZone, topZone) for every
    processor<N>/constant/polyMesh,
  - converts processor<N>/<time>/pointDisplacementz to the vector field
    pointDisplacement (0 0 dz) with the same boundary conditions,
  - replaces the motion solver in constant/dynamicMeshDict.

Usage:
    to_layered.py <case> <time>
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

LOG = logging.getLogger("to_layered")

HEADER = """\
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
"""

DYNAMIC_MESH = """\
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}

dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolverLibs (fvMotionSolvers);
motionSolver    displacementLayeredMotion;

// Every column of the (vertically extruded) mesh stretches linearly between
// its bed point and the fixed lid: no lateral smoothing of the bed
// displacement, so the first cells cannot collapse where the bed varies
// faster than the first-cell height between neighbouring columns (as with a
// Laplacian solve and cell-to-point interpolation)
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
                    type    follow;     // pointDisplacement set by deformingBed
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
"""


def search(pattern: str, text: str) -> str:
    """First group of ``pattern`` in ``text``; KeyError if absent."""
    m = re.search(pattern, text)
    if m is None:
        raise KeyError(pattern)
    return m.group(1)


def patch_range(boundary: str, name: str) -> tuple[int, int]:
    """Start face and size of patch ``name`` in a boundary file text."""
    m = re.search(r"\n\s+" + name + r"\s*\{([^}]*)\}", boundary)
    if m is None:
        raise KeyError(name)
    seg = m.group(1)
    start = int(search(r"startFace\s+(\d+)", seg))
    size = int(search(r"nFaces\s+(\d+)", seg))
    return start, size


def label_list(labels: range) -> str:
    """OpenFOAM ascii List<label> of a contiguous range."""
    body = "\n".join(str(i) for i in labels)
    return f"List<label> {len(labels)}\n(\n{body}\n)"


def write_zones(mesh_dir: Path) -> None:
    """Write cellZones and faceZones of one processor mesh."""
    boundary = (mesh_dir / "boundary").read_text()
    owner_head = (mesh_dir / "owner").read_text(errors="ignore")[:2000]
    n_cells = int(search(r"nCells:\s*(\d+)", owner_head))
    zones = []
    for zone, patch in (("bedZone", "bed"), ("topZone", "top")):
        start, size = patch_range(boundary, patch)
        flips = "List<bool> {0}({1})".format(
            size, " ".join("0" for _ in range(size))) if size else "List<bool> 0()"
        zones.append(f"{zone}\n{{\n    type faceZone;\n    faceLabels "
                     f"{label_list(range(start, start + size))};\n"
                     f"    flipMap {flips};\n}}\n")
    (mesh_dir / "faceZones").write_text(
        HEADER.format(cls="regIOobject", loc="constant/polyMesh",
                      obj="faceZones")
        + f"\n2\n(\n{''.join(zones)})\n")
    (mesh_dir / "cellZones").write_text(
        HEADER.format(cls="regIOobject", loc="constant/polyMesh",
                      obj="cellZones")
        + f"\n1\n(\nallCells\n{{\n    type cellZone;\n    cellLabels "
          f"{label_list(range(n_cells))};\n}}\n)\n")


def scalar_to_vector(text: str) -> str:
    """Turn a pointScalarField (z component) into a pointVectorField."""
    text = text.replace("pointScalarField", "pointVectorField")
    text = re.sub(r"object\s+pointDisplacementz;",
                  "object      pointDisplacement;", text)
    text = re.sub(r"uniform\s+([-\d.eE+]+)\s*;",
                  lambda m: f"uniform (0 0 {m.group(1)});", text)

    def lists(m: re.Match) -> str:
        values = m.group(2).split()
        body = "\n".join(f"(0 0 {v})" for v in values)
        return f"nonuniform List<vector> {m.group(1)}\n(\n{body}\n)"

    text = re.sub(r"nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", lists,
                  text)
    return text


def main() -> None:
    """Convert the case."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("time")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    procs = sorted(args.case.glob("processor*"))
    for proc in procs:
        write_zones(proc / "constant" / "polyMesh")
        src = proc / args.time / "pointDisplacementz"
        (proc / args.time / "pointDisplacement").write_text(
            scalar_to_vector(src.read_text()))
    (args.case / "constant" / "dynamicMeshDict").write_text(DYNAMIC_MESH)
    LOG.info("%s: %d processors converted at time %s", args.case, len(procs),
             args.time)


if __name__ == "__main__":
    main()
