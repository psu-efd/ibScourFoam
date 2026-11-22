#!/usr/bin/env python3
"""Reconstruct the bed level ``<time>/finite-area/zb`` of a decomposed height-field case.

``reconstructPar`` does not reconstruct the finite-area fields of the
height-field bed (its faMesh is built at run time, so there is no
``finite-area/faMesh`` to reconstruct onto). This tool gathers
``processor*/<time>/finite-area/<field>`` into the global floor-patch face
order through ``faceProcAddressing`` (the inverse of decompose_zb.py) and
writes ``<time>/finite-area/<field>`` for every time (or the given ones).

Usage::

    reconstruct_zb.py CASE_DIR [--time T ...] [--field zb]
"""
import argparse
import logging
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bed_field
import foam_file
from decompose_zb import floor_patch

LOG = logging.getLogger("reconstruct_zb")


def processor_dirs(case: str) -> list:
    """Processor directories of a case, sorted by number."""
    procs = [d for d in os.listdir(case) if re.fullmatch(r"processor\d+", d)]
    return sorted(procs, key=lambda d: int(d[9:]))


def reconstruct(case: str, time: str, field: str = "zb") -> bool:
    """Write the global area field of one time; False if a processor lacks it."""
    patch = floor_patch(case)
    gbnd = foam_file.read_boundary(os.path.join(case, "constant", "polyMesh", "boundary"))
    gstart, gn = int(gbnd[patch]["startFace"]), int(gbnd[patch]["nFaces"])
    out = np.full(gn, np.nan)
    for p in processor_dirs(case):
        pdir = os.path.join(case, p)
        bnd = foam_file.read_boundary(os.path.join(pdir, "constant", "polyMesh", "boundary"))
        if patch not in bnd or int(bnd[patch]["nFaces"]) == 0:
            continue
        path = os.path.join(pdir, time, "finite-area", field)
        if not os.path.isfile(path):
            return False
        s, n = int(bnd[patch]["startFace"]), int(bnd[patch]["nFaces"])
        addr = foam_file.read_labels(
            os.path.join(pdir, "constant", "polyMesh", "faceProcAddressing"))
        idx = np.abs(addr[s:s + n]) - 1 - gstart
        out[idx] = bed_field.read_area_field(path, n)
    if np.isnan(out).any():
        raise SystemExit(f"{case} {time}: {int(np.isnan(out).sum())} floor faces missing")
    bed_field.write_area_field(case, out, field, time)
    return True


def main() -> None:
    """Parse the arguments and reconstruct the bed level."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case")
    ap.add_argument("--time", nargs="*", default=None, help="times (default: all)")
    ap.add_argument("--field", default="zb")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    procs = processor_dirs(args.case)
    if not procs:
        raise SystemExit(f"{args.case}: no processor directories")
    times = args.time or foam_file.time_dirs(os.path.join(args.case, procs[0]))
    done = [t for t in times if reconstruct(args.case, t, args.field)]
    LOG.info("reconstruct_zb %s %s: %d of %d times", args.case, args.field,
             len(done), len(times))


if __name__ == "__main__":
    main()
