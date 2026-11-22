"""Helpers shared by the ``make_case.py`` scripts of the test cases.

* copy a template case (``system``, ``constant``, ``0.orig``) into a run
  directory,
* write ``constant/boundaryData/<patch>`` for timeVaryingMappedFixedValue,
* store/load the case set-up (JSON) used later by ``check.py``,
* write initial internal fields from cell-centre functions.
"""

import json
import logging
import os
import shutil
from typing import Callable, Dict, Tuple

import numpy as np

import foam_file

LOG = logging.getLogger(__name__)

SETUP_FILE = "case_setup.json"


def copy_template(template: str, run_dir: str,
                  parts: Tuple[str, ...] = ("system", "constant", "0.orig")
                  ) -> None:
    """Copy template sub-directories into ``run_dir`` (overwriting).

    Parameters
    ----------
    template : str
        Directory holding ``system``, ``constant`` and ``0.orig``.
    run_dir : str
        Destination case directory (created if needed).
    parts : tuple of str
        Sub-directories to copy.
    """
    os.makedirs(run_dir, exist_ok=True)
    for part in parts:
        src = os.path.join(template, part)
        dst = os.path.join(run_dir, part)
        if not os.path.isdir(src):
            continue
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    LOG.info("copied template %s -> %s", template, run_dir)


def write_boundary_data(run_dir: str, patch: str, points: np.ndarray,
                        fields: Dict[str, np.ndarray]) -> None:
    """Write ``constant/boundaryData/<patch>/{points,0/<field>}``.

    Parameters
    ----------
    run_dir : str
        Case directory.
    patch : str
        Patch name.
    points : (N, 3) array
        Sample locations [m] (here: the patch face centres).
    fields : dict
        Field name -> (N,) or (N, 3) values (SI units of the field).
    """
    base = os.path.join(run_dir, "constant", "boundaryData", patch)
    foam_file.write_foam_file(
        os.path.join(base, "points"), "vectorField", "points",
        foam_file.format_vector_list(points))
    for name, vals in fields.items():
        if vals.ndim == 1:
            body = foam_file.format_scalar_list(vals)
            cls = "scalarField"
        else:
            body = foam_file.format_vector_list(vals)
            cls = "vectorField"
        foam_file.write_foam_file(os.path.join(base, "0", name), cls, name,
                                 body)
    LOG.info("boundaryData/%s: %d points, fields %s", patch, len(points),
             sorted(fields))


def save_setup(run_dir: str, setup: Dict[str, object]) -> None:
    """Store the case set-up dictionary as JSON in ``run_dir``."""
    with open(os.path.join(run_dir, SETUP_FILE), "w") as out:
        json.dump(setup, out, indent=2, sort_keys=True)


def load_setup(run_dir: str) -> Dict[str, object]:
    """Load the case set-up dictionary stored by :func:`save_setup`."""
    with open(os.path.join(run_dir, SETUP_FILE)) as src:
        return json.load(src)


def init_fields(run_dir: str,
                func: Callable[[np.ndarray], Dict[str, np.ndarray]]) -> None:
    """Write ``0/<field>`` with internal values func(cell centres).

    ``0.orig`` is copied to ``0`` first; fields not returned by ``func``
    keep their uniform template values.

    Parameters
    ----------
    run_dir : str
        Case directory with a mesh in ``constant/polyMesh``.
    func : callable
        Maps cell centres (nCells, 3) [m] to {field: values}.
    """
    zero = os.path.join(run_dir, "0")
    if os.path.isdir(zero):
        shutil.rmtree(zero)
    shutil.copytree(os.path.join(run_dir, "0.orig"), zero)
    mesh = foam_file.PolyMesh(run_dir)
    values = func(mesh.cell_centres)
    for name, vals in values.items():
        path = os.path.join(zero, name)
        foam_file.replace_internal_field(path, path, vals)
    LOG.info("initialised %s on %d cells", sorted(values), mesh.n_cells)
