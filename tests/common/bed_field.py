"""Bed level of a height-field immersed-boundary case (``ibScourPimpleFoam``).

The bed lives on a finite-area mesh built at run time from the floor patch,
in patch-face order. This module gives the plan-view geometry of the floor
patch and reads and writes the area fields (``<time>/finite-area/zb``) of a
serial (or reconstructed) case. ASCII files only.
"""

import os
import re
from typing import Dict, List, Optional, Tuple

import foam_file
import numpy as np


def floor_faces(case: str, patch: str = "bottom") -> Dict[str, np.ndarray]:
    """Return centres ``xc`` (n, 3) and plan areas ``area`` (n,) of a patch."""
    mesh = os.path.join(case, "constant", "polyMesh")
    points = foam_file.read_points(os.path.join(mesh, "points"))
    offsets, labels = foam_file.read_faces(os.path.join(mesh, "faces"))
    bnd = foam_file.read_boundary(os.path.join(mesh, "boundary"))[patch]
    start, n = int(bnd["startFace"]), int(bnd["nFaces"])
    xc = np.zeros((n, 3))
    area = np.zeros(n)
    for i in range(n):
        p = points[labels[offsets[start + i]:offsets[start + i + 1]]]
        xc[i] = p.mean(axis=0)
        # plan area of the (planar, horizontal) face: shoelace formula
        x, y = p[:, 0], p[:, 1]
        area[i] = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return {"xc": xc, "area": area}


def write_area_field(case: str, values: np.ndarray, name: str = "zb",
                     time: str = "0") -> str:
    """Write ``<time>/finite-area/<name>`` (zero-gradient edges); return the path."""
    out = os.path.join(case, time, "finite-area")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name)
    body = "\n".join("%.10g" % v for v in values)
    with open(path, "w") as fh:
        fh.write(
            "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            "    class       areaScalarField;\n"
            "    location    \"%s/finite-area\";\n    object      %s;\n}\n\n"
            "dimensions      [0 1 0 0 0 0 0];\n\n"
            "internalField   nonuniform List<scalar>\n%d\n(\n%s\n)\n;\n\n"
            "boundaryField\n{\n    \".*\"\n    {\n        type            zeroGradient;\n"
            "    }\n}\n" % (time, name, len(values), body))
    return path


def read_area_field(path: str, n: int) -> np.ndarray:
    """Read the internal field of an ASCII area field with ``n`` faces."""
    with open(path) as fh:
        text = fh.read()
    m = re.search(r"internalField\s+uniform\s+([-+.eE\d]+)\s*;", text)
    if m:
        return np.full(n, float(m.group(1)))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;",
                  text, flags=re.S)
    if m is None:
        raise ValueError("cannot parse the internalField of %s" % path)
    vals = np.array(m.group(2).split(), dtype=float)
    if len(vals) != n:
        raise ValueError("%s: %d values, %d faces expected" % (path, len(vals), n))
    return vals


def bed_times(case: str, name: str = "zb") -> List[Tuple[float, str]]:
    """Return sorted (time, path) of every ``<time>/finite-area/<name>``."""
    out = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        p = os.path.join(case, d, "finite-area", name)
        if os.path.isfile(p):
            out.append((t, p))
    return sorted(out)


def area_to_patch(case: str, patch: str = "bottom") -> np.ndarray:
    """Floor-patch face index of every area face.

    The area mesh written by makeFaMesh lists its faces in
    ``constant/finite-area/faMesh/faceLabels``; an area mesh built by the
    solver at run time follows the patch face order.
    """
    bnd = foam_file.read_boundary(
        os.path.join(case, "constant", "polyMesh", "boundary"))[patch]
    start, n = int(bnd["startFace"]), int(bnd["nFaces"])
    path = os.path.join(case, "constant", "finite-area", "faMesh", "faceLabels")
    if not os.path.isfile(path):
        return np.arange(n)
    idx = foam_file.read_labels(path) - start
    if len(idx) != n or idx.min() < 0 or idx.max() >= n:
        raise ValueError("%s: faceLabels are not the faces of patch %s" % (path, patch))
    return idx


def bed_history(case: str, patch: str = "bottom",
                name: str = "zb") -> Tuple[Dict[str, np.ndarray], List[float], np.ndarray]:
    """Return the floor geometry, the output times and z_b (n_times, n_faces).

    z_b is returned in floor-patch face order (the order of ``geo``).
    """
    geo = floor_faces(case, patch)
    n = len(geo["area"])
    order = area_to_patch(case, patch)
    times, paths = zip(*bed_times(case, name)) if bed_times(case, name) else ((), ())
    z = np.zeros((len(paths), n))
    for k, p in enumerate(paths):
        z[k, order] = read_area_field(p, n)
    return geo, list(times), z


def grid_index(xc: np.ndarray, tol: float = 1e-9) -> Tuple[np.ndarray, np.ndarray]:
    """Map face centres of a structured floor to integer (i, j) indices."""
    xs = np.unique(np.round(xc[:, 0] / tol).astype(np.int64))
    ys = np.unique(np.round(xc[:, 1] / tol).astype(np.int64))
    i = np.searchsorted(xs, np.round(xc[:, 0] / tol).astype(np.int64))
    j = np.searchsorted(ys, np.round(xc[:, 1] / tol).astype(np.int64))
    return i, j


def log_bed_lines(log_path: str) -> List[str]:
    """Return the ``Bed: MF ...`` report lines of a solver log."""
    if not os.path.isfile(log_path):
        return []
    with open(log_path, errors="replace") as fh:
        return [line for line in fh if line.startswith("Bed: MF")]


def last_float(pattern: str, lines: List[str]) -> Optional[float]:
    """Return the first group of ``pattern`` in the last matching line."""
    for line in reversed(lines):
        m = re.search(pattern, line)
        if m:
            return float(m.group(1))
    return None


def column_index(x: np.ndarray, tol: float = 1e-4) -> Tuple[np.ndarray, np.ndarray]:
    """Group coordinates into columns: values closer than ``tol`` [m] merge.

    Returns the column coordinates (mean of the members) and the column
    index of every value.
    """
    order = np.argsort(x)
    xs = x[order]
    new = np.concatenate([[True], np.diff(xs) > tol])
    col_sorted = np.cumsum(new) - 1
    idx = np.empty(len(x), dtype=np.int64)
    idx[order] = col_sorted
    ux = np.bincount(idx, weights=x) / np.bincount(idx)
    return ux, idx
