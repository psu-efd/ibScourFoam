"""Decompose the bed level ``0/finite-area/zb`` (or another area field) of a height-field IB case.

The height-field bed lives on a finite-area mesh built at run time from the
floor patch (face order = patch face order), so ``decomposePar`` does not
decompose this field. This tool maps the global ``0/finite-area/zb`` onto
every ``processor*`` directory through ``faceProcAddressing`` and writes
``processor*/0/finite-area/zb`` with processor boundary entries. With
``--field zRigid`` the non-erodible floor elevation is decomposed the same way.

Usage::

    decompose_zb.py CASE_DIR [--time 0] [--field zb]

The floor patch is read from ``constant/immersedBoundaryDict`` (``bed``
sub-dictionary, ``floorPatch`` keyword). ASCII files only.
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foam_file

LOG = logging.getLogger("decompose_zb")


def floor_patch(case: str) -> str:
    """Name of the floor patch from the immersed boundary dictionary."""
    with open(os.path.join(case, "constant", "immersedBoundaryDict")) as src:
        text = src.read()
    m = re.search(r"floorPatch\s+(\w+)\s*;", text)
    if not m:
        raise SystemExit(f"{case}: no floorPatch entry in immersedBoundaryDict")
    return m.group(1)


def zb_body(zb: np.ndarray) -> str:
    """File body of an areaScalarField of lengths [m] with processor boundaries."""
    return (
        "dimensions      [0 1 0 0 0 0 0];\n\ninternalField   nonuniform "
        "List<scalar> " + foam_file.format_scalar_list(zb)
        + ";\n\nboundaryField\n{\n    \".*\"\n    {\n        type"
        "            zeroGradient;\n    }\n    \"procBoundary.*\"\n    {"
        "\n        type            processor;\n    }\n}\n"
    )


def decompose(case: str, time: str, field: str = "zb") -> int:
    """Write the area field for every processor directory; return the processor count."""
    patch = floor_patch(case)
    gbnd = foam_file.read_boundary(
        os.path.join(case, "constant", "polyMesh", "boundary"))
    if patch not in gbnd:
        raise SystemExit(f"{case}: floor patch {patch} not in the global mesh")
    gstart, gn = gbnd[patch]["startFace"], gbnd[patch]["nFaces"]
    zb = foam_file.read_internal_field(
        os.path.join(case, time, "finite-area", field), gn)
    if len(zb) != gn:
        raise SystemExit(f"{case}: {field} has {len(zb)} values, patch {patch} "
                         f"has {gn} faces")
    procs = sorted(d for d in os.listdir(case) if re.fullmatch(r"processor\d+", d))
    if not procs:
        raise SystemExit(f"{case}: no processor directories (run decomposePar first)")
    total = 0
    for p in procs:
        pdir = os.path.join(case, p)
        bnd = foam_file.read_boundary(
            os.path.join(pdir, "constant", "polyMesh", "boundary"))
        if patch not in bnd or bnd[patch]["nFaces"] == 0:
            local = np.zeros(0)
        else:
            addr = foam_file.read_labels(
                os.path.join(pdir, "constant", "polyMesh", "faceProcAddressing"))
            s, n = bnd[patch]["startFace"], bnd[patch]["nFaces"]
            gface = np.abs(addr[s:s + n]) - 1
            idx = gface - gstart
            if idx.min() < 0 or idx.max() >= gn:
                raise SystemExit(f"{p}: floor faces map outside the global patch")
            local = zb[idx]
        total += len(local)
        foam_file.write_foam_file(
            os.path.join(pdir, time, "finite-area", field), "areaScalarField",
            field, zb_body(local))
        LOG.info("%s: %d floor faces", p, len(local))
    if total != gn:
        raise SystemExit(f"{case}: processors hold {total} floor faces, "
                         f"global patch has {gn}")
    return len(procs)


def main() -> None:
    """Parse the arguments and decompose the bed level."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case")
    ap.add_argument("--time", default="0", help="time directory holding the field")
    ap.add_argument("--field", default="zb",
                    help="area field to decompose (zb, or zRigid for the non-erodible floor)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    n = decompose(args.case, args.time, args.field)
    LOG.info("decompose_zb %s %s: %d processors", args.case, args.field, n)


if __name__ == "__main__":
    main()
