"""Minimal reader/writer for ASCII OpenFOAM files used by the test suite.

The functions here cover only what the ibScourFoam tests need:

* writing dictionaries/fields with the standard v2606 ``FoamFile`` banner,
* reading ``constant/polyMesh`` (points, faces, owner, neighbour, boundary)
  in ASCII format, including the ``faceCompactList`` layout,
* reading the internal field and patch values of vol/surface fields
  (``uniform`` and ``nonuniform List<...>`` entries),
* replacing the ``internalField`` entry of a field file.

Only numpy (>= 1.14) and the Python 3.6 standard library are required.
Binary OpenFOAM files are not supported; all test cases write ASCII.
"""

import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

# The standard v2606 file banner (80 columns wide)
BANNER = "\n".join([
    "/*--------------------------------*- C++ -*-----------------------"
    "-----------*\\",
    r"| =========                 |                                    "
    r"             |",
    r"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    r"           |",
    r"|  \\    /   O peration     | Version:  v2606                      "
    r"           |",
    r"|   \\  /    A nd           | Website:  www.openfoam.com           "
    r"           |",
    r"|    \\/     M anipulation  |                                    "
    r"             |",
    "\\*---------------------------------------------------------------"
    "------------*/",
]) + "\n"

SEPARATOR = (
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
    "* * * * //\n"
)

FOOTER = (
    "\n// *****************************************************************"
    "******** //\n"
)


# * * * * * * * * * * * * * * * * Writing * * * * * * * * * * * * * * * * * #

def foam_header(cls: str, obj: str, location: Optional[str] = None) -> str:
    """Return the banner, the ``FoamFile`` sub-dictionary and the separator.

    Parameters
    ----------
    cls : str
        OpenFOAM class name, e.g. ``dictionary`` or ``vectorField``.
    obj : str
        Object name written in the header.
    location : str, optional
        Value of the ``location`` entry (omitted when ``None``).
    """
    lines = [
        BANNER.rstrip("\n"),
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       {};".format(cls),
    ]
    if location is not None:
        lines.append('    location    "{}";'.format(location))
    lines.append("    object      {};".format(obj))
    lines.append("}")
    return "\n".join(lines) + "\n" + SEPARATOR + "\n"


def write_foam_file(
    path: str,
    cls: str,
    obj: str,
    body: str,
    location: Optional[str] = None,
) -> None:
    """Write an OpenFOAM file with header, body and footer.

    Parameters
    ----------
    path : str
        Output file; parent directories are created.
    cls, obj, location :
        See :func:`foam_header`.
    body : str
        File content between the separator and the footer.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as out:
        out.write(foam_header(cls, obj, location))
        out.write(body.rstrip("\n") + "\n")
        out.write(FOOTER)


def format_vector_list(values: np.ndarray) -> str:
    """Format an (N, 3) array as an OpenFOAM ASCII list body."""
    rows = ["({:.8g} {:.8g} {:.8g})".format(*v) for v in values]
    return "{}\n(\n{}\n)\n".format(len(rows), "\n".join(rows))


def format_scalar_list(values: np.ndarray) -> str:
    """Format a 1-D array as an OpenFOAM ASCII list body."""
    rows = ["{:.8g}".format(v) for v in values]
    return "{}\n(\n{}\n)\n".format(len(rows), "\n".join(rows))


def write_parameters(path: str, params: Dict[str, object], note: str) -> None:
    """Write a dictionary of scalar parameters (for ``#include``).

    Parameters
    ----------
    path : str
        Output file, e.g. ``system/caseParameters``.
    params : dict
        Keyword -> value (numbers or strings).
    note : str
        One-line comment written above the entries.
    """
    width = max(len(k) for k in params) + 1
    width = max(width, 16)
    lines = ["// " + note, ""]
    for key, val in params.items():
        if isinstance(val, float):
            val = "{:.10g}".format(val)
        lines.append("{}{};".format(key.ljust(width), val))
    write_foam_file(path, "dictionary", os.path.basename(path),
                    "\n".join(lines))


def replace_internal_field(
    template: str, out: str, values: np.ndarray
) -> None:
    """Copy a field file, replacing ``internalField`` by a nonuniform list.

    Parameters
    ----------
    template : str
        Source field file with an ``internalField uniform ...;`` line.
    out : str
        Destination file.
    values : numpy.ndarray
        Shape (nCells,) for scalars or (nCells, 3) for vectors.
    """
    with open(template) as src:
        text = src.read()
    if values.ndim == 1:
        entry = ("internalField   nonuniform List<scalar>\n"
                 + format_scalar_list(values))
    else:
        entry = ("internalField   nonuniform List<vector>\n"
                 + format_vector_list(values))
    entry = entry.rstrip("\n") + ";"
    new_text, count = re.subn(
        r"^internalField\s+uniform[^;]*;", lambda _m: entry, text,
        count=1, flags=re.M,
    )
    if count != 1:
        raise ValueError("no 'internalField uniform' entry in " + template)
    with open(out, "w") as dst:
        dst.write(new_text)


# * * * * * * * * * * * * * * * * Reading * * * * * * * * * * * * * * * * * #

_COMMENT_RE = re.compile(r"//.*?$|/\*.*?\*/", re.S | re.M)


def _strip_header(text: str) -> str:
    """Remove comments and the ``FoamFile`` sub-dictionary."""
    text = _COMMENT_RE.sub("", text)
    match = re.search(r"FoamFile\s*\{[^}]*\}", text)
    if match:
        text = text[match.end():]
    return text


def _read_text(path: str) -> str:
    """Return the file content without comments and header."""
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        raise IOError("compressed files are not supported: " + path)
    with open(path) as src:
        return _strip_header(src.read())


def _parse_list_at(text: str, pos: int, kind: str) -> Tuple[np.ndarray, int]:
    """Parse an ASCII list starting at ``pos`` (at the size or the brace).

    Parameters
    ----------
    text : str
        Text containing the list.
    pos : int
        Index at or before the list size.
    kind : str
        ``scalar``, ``label``, ``vector`` or ``tensor``.

    Returns
    -------
    (numpy.ndarray, int)
        The values and the index just after the closing parenthesis.
    """
    match = re.compile(r"\s*(\d+)\s*([({])").match(text, pos)
    if match is None:
        # Size-less list "( ... )"
        match = re.compile(r"\s*()\s*(\()").match(text, pos)
        if match is None:
            raise ValueError("cannot parse list near: " + text[pos:pos + 80])
    size = int(match.group(1)) if match.group(1) else -1
    start = match.end()
    if match.group(2) == "{":
        # Uniform compact list N{value}
        end = text.index("}", start)
        inner = text[start:end].replace("(", " ").replace(")", " ")
        val = np.array(inner.split(), dtype=float)
        ncomp = {"scalar": 1, "label": 1, "vector": 3, "tensor": 9}[kind]
        arr = np.tile(val, (size, 1)) if ncomp > 1 else np.full(size, val[0])
        return arr, end + 1
    if kind in ("scalar", "label"):
        end = text.index(")", start)
        arr = np.array(text[start:end].split(), dtype=float)
        if kind == "label":
            arr = arr.astype(np.int64)
        return arr, end + 1
    # Vector/tensor list: find the closing parenthesis of the outer list
    depth = 1
    idx = start
    while depth > 0:
        char = text[idx]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        idx += 1
    inner = text[start:idx - 1].replace("(", " ").replace(")", " ")
    ncomp = 3 if kind == "vector" else 9
    arr = np.array(inner.split(), dtype=float).reshape(-1, ncomp)
    if size >= 0 and arr.shape[0] != size:
        raise ValueError(
            "list size mismatch ({} vs {})".format(size, len(arr)))
    return arr, idx


def read_points(path: str) -> np.ndarray:
    """Read ``polyMesh/points`` -> (nPoints, 3) array [m]."""
    text = _read_text(path)
    arr, _ = _parse_list_at(text, 0, "vector")
    return arr


def read_labels(path: str) -> np.ndarray:
    """Read a ``labelList`` file such as ``owner`` or ``neighbour``."""
    text = _read_text(path)
    arr, _ = _parse_list_at(text, 0, "label")
    return arr


def read_faces(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read ``polyMesh/faces`` in either list or compact layout.

    Returns
    -------
    (offsets, labels)
        CSR-like layout: face i has points labels[offsets[i]:offsets[i+1]].
    """
    with open(path) as src:
        raw = src.read()
    compact = re.search(r"class\s+faceCompactList", raw) is not None
    text = _strip_header(raw)
    if compact:
        offsets, pos = _parse_list_at(text, 0, "label")
        labels, _ = _parse_list_at(text, pos, "label")
        return offsets, labels
    # Plain faceList: N ( 4(a b c d) ... )
    match = re.compile(r"\s*(\d+)\s*\(").match(text)
    body = text[match.end():]
    counts = []
    labels = []
    for size, inner in re.findall(r"(\d+)\(([^)]*)\)", body):
        counts.append(int(size))
        labels.extend(inner.split())
    offsets = np.concatenate([[0], np.cumsum(counts)])
    return offsets, np.array(labels, dtype=np.int64)


def read_boundary(path: str) -> Dict[str, Dict[str, object]]:
    """Read ``polyMesh/boundary`` -> {patch: {type, nFaces, startFace}}."""
    text = _read_text(path)
    patches = {}
    for name, body in re.findall(r"(\w+)\s*\{([^{}]*)\}", text):
        entry = {}
        for key, val in re.findall(r"(\w+)\s+([^;]+);", body):
            entry[key] = val.strip()
        if "nFaces" in entry:
            patches[name] = {
                "type": entry.get("type", ""),
                "nFaces": int(entry["nFaces"]),
                "startFace": int(entry["startFace"]),
            }
    return patches


class PolyMesh:
    """ASCII ``polyMesh`` with derived face and cell geometry.

    Attributes
    ----------
    points : (nPoints, 3) array [m]
    owner, neighbour : label arrays
    boundary : dict of patch info
    face_centres, face_areas : (nFaces, 3) arrays [m], [m^2] (area vectors)
    cell_centres : (nCells, 3) array [m] (area-weighted face-centre mean)
    """

    def __init__(self, case: str, region: str = "") -> None:
        """Load the mesh of ``case`` (optionally a sub-region)."""
        mdir = os.path.join(case, "constant", region, "polyMesh")
        self.points = read_points(os.path.join(mdir, "points"))
        self.offsets, self.labels = read_faces(os.path.join(mdir, "faces"))
        self.owner = read_labels(os.path.join(mdir, "owner"))
        self.neighbour = read_labels(os.path.join(mdir, "neighbour"))
        self.boundary = read_boundary(os.path.join(mdir, "boundary"))
        self.n_faces = len(self.owner)
        self.n_internal = len(self.neighbour)
        last = max(int(self.owner.max()),
                   int(self.neighbour.max()) if len(self.neighbour) else 0)
        self.n_cells = last + 1
        self._face_geometry()
        self._cell_centres()

    def _face_geometry(self) -> None:
        """Compute face centres and area vectors (triangle fan)."""
        sizes = np.diff(self.offsets)
        nf = self.n_faces
        centres = np.zeros((nf, 3))
        areas = np.zeros((nf, 3))
        for nv in np.unique(sizes):
            idx = np.where(sizes == nv)[0]
            pts = self.points[
                self.labels[self.offsets[idx][:, None] + np.arange(nv)]
            ]
            mean = pts.mean(axis=1)
            fc = np.zeros((len(idx), 3))
            fa = np.zeros((len(idx), 3))
            mag_sum = np.zeros(len(idx))
            for k in range(nv):
                p0 = pts[:, k]
                p1 = pts[:, (k + 1) % nv]
                tri_n = np.cross(p0 - mean, p1 - mean) * 0.5
                tri_c = (p0 + p1 + mean) / 3.0
                tri_a = np.linalg.norm(tri_n, axis=1)
                fa += tri_n
                fc += tri_c * tri_a[:, None]
                mag_sum += tri_a
            centres[idx] = fc / np.maximum(mag_sum, 1e-300)[:, None]
            areas[idx] = fa
        self.face_centres = centres
        self.face_areas = areas

    def _cell_centres(self) -> None:
        """Approximate cell centres by area-weighted face-centre averages."""
        mag_a = np.linalg.norm(self.face_areas, axis=1)
        acc = np.zeros((self.n_cells, 3))
        wsum = np.zeros(self.n_cells)
        weighted = self.face_centres * mag_a[:, None]
        for cells, n in ((self.owner, self.n_faces),
                         (self.neighbour, self.n_internal)):
            for d in range(3):
                acc[:, d] += np.bincount(
                    cells, weighted[:n, d], minlength=self.n_cells
                )
            wsum += np.bincount(cells, mag_a[:n], minlength=self.n_cells)
        self.cell_centres = acc / wsum[:, None]

    def patch_slice(self, name: str) -> slice:
        """Return the face slice of patch ``name``."""
        info = self.boundary[name]
        start = info["startFace"]
        return slice(start, start + info["nFaces"])


def _section(text: str, keyword: str) -> int:
    """Return the index just after ``keyword`` at top level."""
    match = re.search(r"(^|\s)" + keyword + r"\s", text)
    if match is None:
        raise KeyError(keyword)
    return match.end()


def _parse_value(text: str, pos: int, n: int, kind: str) -> np.ndarray:
    """Parse ``uniform X`` or ``nonuniform List<..> N(...)`` at ``pos``."""
    match = re.compile(r"\s*(uniform|nonuniform)\s*").match(text, pos)
    if match is None:
        raise ValueError("cannot parse value near: " + text[pos:pos + 60])
    pos = match.end()
    if match.group(1) == "uniform":
        if kind == "scalar":
            val = float(re.compile(r"\s*([^;\s]+)").match(text, pos).group(1))
            return np.full(n, val)
        inner = re.compile(r"\s*\(([^)]*)\)").match(text, pos).group(1)
        return np.tile(np.array(inner.split(), dtype=float), (n, 1))
    match = re.compile(r"\s*List<\w+>\s*").match(text, pos)
    pos = match.end()
    arr, _ = _parse_list_at(text, pos, kind)
    return arr


def field_kind(path: str) -> str:
    """Return ``scalar`` or ``vector`` from the field class in the header."""
    with open(path) as src:
        head = src.read(2000)
    match = re.search(r"class\s+(\w+);", head)
    cls = match.group(1) if match else ""
    return "vector" if "Vector" in cls else "scalar"


def read_internal_field(path: str, n_cells: int) -> np.ndarray:
    """Read the ``internalField`` of a vol/surface field file.

    Parameters
    ----------
    path : str
        Field file.
    n_cells : int
        Number of values expected (cells, or internal faces for
        surface fields), used to expand ``uniform`` entries.
    """
    kind = field_kind(path)
    text = _read_text(path)
    pos = _section(text, "internalField")
    return _parse_value(text, pos, n_cells, kind)


def read_patch_field(path: str, patch: str, n_faces: int) -> np.ndarray:
    """Read the ``value`` entry of one patch of a field file.

    Returns zeros if the patch has no ``value`` entry (e.g. zeroGradient
    patches of volFields, which this test suite does not need).
    """
    kind = field_kind(path)
    text = _read_text(path)
    pos = _section(text, "boundaryField")
    match = re.compile(r"\b" + re.escape(patch) + r"\s*\{").search(text, pos)
    if match is None:
        raise KeyError(patch)
    # Find the end of the patch block (skip parentheses contents)
    idx = match.end()
    depth = 1
    paren = 0
    while depth > 0:
        char = text[idx]
        if char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        elif paren == 0 and char == "{":
            depth += 1
        elif paren == 0 and char == "}":
            depth -= 1
        idx += 1
    block = text[match.end():idx - 1]
    vmatch = re.search(r"(^|\s)value\s", block)
    if vmatch is None:
        shape = (n_faces,) if kind == "scalar" else (n_faces, 3)
        return np.zeros(shape)
    return _parse_value(block, vmatch.end(), n_faces, kind)


def time_dirs(case: str) -> List[str]:
    """Return the numeric time directories of ``case`` sorted by value."""
    names = []
    for name in os.listdir(case):
        try:
            float(name)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, name)):
            names.append(name)
    return sorted(names, key=float)


# Editing of dictionary files

def set_entry(path: str, key: str, value: str) -> None:
    """Replace the value of a top-level ``key value;`` entry in a file."""
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        words = line.split()
        if words and words[0] == key:
            lines[i] = f"{key:<16}{value};"
    with open(path, "w") as out:
        out.write("\n".join(lines) + "\n")


def remove_block(path: str, key: str) -> None:
    """Remove a top-level dictionary block ``key { ... }`` from a file."""
    text = open(path).read()
    start = text.find("\n" + key)
    if start < 0:
        return
    open_brace = text.find("{", start)
    depth, i = 0, open_brace
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    with open(path, "w") as out:
        out.write(text[:start] + text[i + 1:])
