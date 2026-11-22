"""STL helpers for bed surfaces and simple bodies (numpy only).

* :func:`height_grid_stl` triangulates a structured height field
  z_b(x, y) and writes it as an ASCII STL (the immersed bed).
* :func:`cylinder_stl` writes a closed vertical cylinder (bridge pier).
* :func:`read_stl` reads ASCII or binary STL files and merges duplicate
  vertices, so that bed elevations can be analysed per vertex.

All lengths are in metres.
"""

import os
import struct
from typing import Callable, Tuple

import numpy as np


def _write_ascii(path: str, name: str, tris: np.ndarray) -> None:
    """Write triangles (nTri, 3, 3) [m] as an ASCII STL solid ``name``."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    normals = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    mags = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(mags, 1e-300)[:, None]
    with open(path, "w") as out:
        out.write("solid {}\n".format(name))
        for nrm, tri in zip(normals, tris):
            out.write("  facet normal {:.9e} {:.9e} {:.9e}\n".format(*nrm))
            out.write("    outer loop\n")
            for vtx in tri:
                out.write(
                    "      vertex {:.9e} {:.9e} {:.9e}\n".format(*vtx)
                )
            out.write("    endloop\n  endfacet\n")
        out.write("endsolid {}\n".format(name))


def _write_binary(path: str, name: str, tris: np.ndarray) -> None:
    """Write triangles (nTri, 3, 3) [m] as a binary STL file."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    normals = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    mags = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(mags, 1e-300)[:, None]
    dtype = np.dtype([("n", "<f4", (3,)), ("v", "<f4", (3, 3)),
                      ("a", "<u2")])
    rec = np.zeros(len(tris), dtype=dtype)
    rec["n"] = normals
    rec["v"] = tris
    header = "binary STL {}".format(name).encode()[:80].ljust(80, b" ")
    with open(path, "wb") as out:
        out.write(header)
        out.write(struct.pack("<I", len(tris)))
        out.write(rec.tobytes())


def write_triangles(path: str, name: str, tris: np.ndarray,
                    binary: bool = False) -> None:
    """Write triangles (nTri, 3, 3) [m] as ASCII or binary STL."""
    if binary:
        _write_binary(path, name, tris)
    else:
        _write_ascii(path, name, tris)


def height_grid_stl(
    path: str,
    name: str,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    spacing: Tuple[float, float],
    z_func: Callable[[np.ndarray, np.ndarray], np.ndarray],
    binary: bool = False,
) -> int:
    """Write a triangulated height field z = z_func(x, y) as STL.

    Each grid quad is split into two triangles along alternating
    diagonals (a "union-jack"-free, symmetric pattern is not needed for
    these tests; alternation avoids a preferred diagonal direction).
    Normals point upwards (+z), i.e. into the fluid.

    Parameters
    ----------
    path : str
        Output file.
    name : str
        Solid name (should equal the IB object name).
    x_range, y_range : (float, float)
        Extent of the surface [m]; should exceed the mesh slightly.
    spacing : (float, float)
        Target vertex spacing in x and y [m].
    z_func : callable
        Bed elevation function z_b(x, y) [m] (vectorised).
    binary : bool
        Write binary instead of ASCII STL (smaller for large beds).

    Returns
    -------
    int
        Number of triangles written.
    """
    nx = max(1, int(round((x_range[1] - x_range[0]) / spacing[0])))
    ny = max(1, int(round((y_range[1] - y_range[0]) / spacing[1])))
    xs = np.linspace(x_range[0], x_range[1], nx + 1)
    ys = np.linspace(y_range[0], y_range[1], ny + 1)
    xg, yg = np.meshgrid(xs, ys, indexing="ij")
    zg = z_func(xg, yg)
    pts = np.stack([xg, yg, zg], axis=-1)
    tris = []
    for i in range(nx):
        for j in range(ny):
            p00, p10 = pts[i, j], pts[i + 1, j]
            p01, p11 = pts[i, j + 1], pts[i + 1, j + 1]
            if (i + j) % 2 == 0:
                tris.append((p00, p10, p11))
                tris.append((p00, p11, p01))
            else:
                tris.append((p00, p10, p01))
                tris.append((p10, p11, p01))
    tris = np.array(tris)
    write_triangles(path, name, tris, binary)
    return len(tris)


def cylinder_stl(
    path: str,
    name: str,
    centre: Tuple[float, float],
    radius: float,
    z_range: Tuple[float, float],
    n_theta: int = 64,
    n_z: int = 8,
) -> int:
    """Write a closed vertical circular cylinder as ASCII STL.

    Parameters
    ----------
    path, name : str
        Output file and solid name.
    centre : (float, float)
        Axis position (x, y) [m].
    radius : float
        Cylinder radius [m].
    z_range : (float, float)
        Bottom and top elevations [m].
    n_theta, n_z : int
        Facets around the circumference and along the height.

    Returns
    -------
    int
        Number of triangles written.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta + 1)
    zs = np.linspace(z_range[0], z_range[1], n_z + 1)
    cx, cy = centre
    ring = np.stack(
        [cx + radius * np.cos(theta), cy + radius * np.sin(theta)], axis=1
    )
    tris = []
    for k in range(n_z):
        for i in range(n_theta):
            a = (*ring[i], zs[k])
            b = (*ring[i + 1], zs[k])
            c = (*ring[i + 1], zs[k + 1])
            d = (*ring[i], zs[k + 1])
            tris.append((a, b, c))
            tris.append((a, c, d))
    for zc, flip in ((z_range[0], True), (z_range[1], False)):
        cpt = (cx, cy, zc)
        for i in range(n_theta):
            a = (*ring[i], zc)
            b = (*ring[i + 1], zc)
            tris.append((cpt, b, a) if flip else (cpt, a, b))
    tris = np.array(tris, dtype=float)
    _write_ascii(path, name, tris)
    return len(tris)


def read_stl(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read an ASCII or binary STL file and merge duplicate vertices.

    Returns
    -------
    (points, triangles)
        points: (nPoints, 3) unique vertices [m];
        triangles: (nTri, 3) vertex indices.
    """
    with open(path, "rb") as src:
        data = src.read()
    is_ascii = data[:5] == b"solid" and b"facet" in data[:400]
    if is_ascii:
        verts = [
            line.split()[1:4]
            for line in data.decode(errors="replace").splitlines()
            if line.strip().startswith("vertex")
        ]
        raw = np.array(verts, dtype=float).reshape(-1, 3, 3)
    else:
        ntri = struct.unpack("<I", data[80:84])[0]
        dtype = np.dtype(
            [("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("a", "<u2")]
        )
        rec = np.frombuffer(data[84:84 + 50 * ntri], dtype=dtype)
        raw = rec["v"].astype(float)
    flat = raw.reshape(-1, 3)
    # Merge vertices on a 1e-9 m lattice
    keys = np.round(flat / 1e-9).astype(np.int64)
    _, first, inverse = np.unique(
        keys, axis=0, return_index=True, return_inverse=True
    )
    points = flat[first]
    return points, inverse.reshape(-1, 3)


def projected_volume(points: np.ndarray, tris: np.ndarray,
                     z_ref: float = 0.0,
                     mask_xy: Callable[[np.ndarray, np.ndarray], np.ndarray]
                     = None) -> float:
    """Return the volume between the surface and z = z_ref [m^3].

    Each triangle contributes (mean z - z_ref) times its horizontal
    projected area; triangles whose centroid fails ``mask_xy`` are skipped.
    """
    tri_pts = points[tris]
    centroid = tri_pts.mean(axis=1)
    edge1 = tri_pts[:, 1, :2] - tri_pts[:, 0, :2]
    edge2 = tri_pts[:, 2, :2] - tri_pts[:, 0, :2]
    area = 0.5 * np.abs(edge1[:, 0] * edge2[:, 1] - edge1[:, 1] * edge2[:, 0])
    keep = np.ones(len(tris), dtype=bool)
    if mask_xy is not None:
        keep = mask_xy(centroid[:, 0], centroid[:, 1])
    return float(np.sum((centroid[keep, 2] - z_ref) * area[keep]))


def triangle_slopes(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Return |grad z| of every triangle (tan of the slope angle) [-]."""
    tri_pts = points[tris]
    nrm = np.cross(tri_pts[:, 1] - tri_pts[:, 0],
                   tri_pts[:, 2] - tri_pts[:, 0])
    nz = np.abs(nrm[:, 2])
    horiz = np.hypot(nrm[:, 0], nrm[:, 1])
    return horiz / np.maximum(nz, 1e-300)
