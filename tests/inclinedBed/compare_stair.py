#!/usr/bin/env python3
"""Staircase body-fitted run vs the IB and the smooth body-fitted runs.

The staircase run is the IB Cartesian mesh cut to the cells that the IB
treats as fluid (subsetMesh; exposed faces = wall patch ``bed`` with the
body-fitted wall law). On the up-ramp its bed consists of treads
(horizontal faces) and risers (vertical faces). The script tabulates per
pit segment the mean tread shear stress tau_x of the staircase, the IB
tau_x and the smooth body-fitted tau_x, and the streamwise force per unit
plan length on the staircase split into tread shear, riser shear and riser
pressure (form drag).

Usage::

    compare_stair.py --stair RUN --ib RUN --bf RUN [--time T]
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from compare import RHO, read_bed, segments

LOG = logging.getLogger("inclinedBed.compare_stair")


def stair_faces(run: Path, time: float) -> dict[str, np.ndarray]:
    """Centres, unit normals (into the fluid), areas, tau, p of the bed."""
    foam = run / "case.foam"
    foam.touch()
    r = vtk.vtkOpenFOAMReader()
    r.SetFileName(str(foam))
    r.UpdateInformation()
    r.DisableAllPatchArrays()
    r.SetPatchArrayStatus("patch/bed", 1)
    r.SetCreateCellToPoint(0)
    r.UpdateTimeStep(time)
    r.Update()
    it = r.GetOutput().NewIterator()
    it.InitTraversal()
    bed = None
    while not it.IsDoneWithTraversal():
        obj = it.GetCurrentDataObject()
        if obj.GetNumberOfCells() > 0:
            bed = obj
        it.GoToNextItem()
    norms = vtk.vtkPolyDataNormals()
    geo = vtk.vtkGeometryFilter()
    geo.SetInputData(bed)
    geo.Update()
    norms.SetInputData(geo.GetOutput())
    norms.ComputeCellNormalsOn()
    norms.ComputePointNormalsOff()
    norms.SplittingOff()
    norms.ConsistencyOff()
    norms.AutoOrientNormalsOff()
    norms.Update()
    poly = norms.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(poly)
    cc.Update()
    size = vtk.vtkCellSizeFilter()
    size.SetInputData(poly)
    size.Update()
    cd = poly.GetCellData()
    return {
        "x": vtk_to_numpy(cc.GetOutput().GetPoints().GetData()),
        "n": vtk_to_numpy(cd.GetNormals()),
        "A": vtk_to_numpy(size.GetOutput().GetCellData().GetArray("Area")),
        "tau": vtk_to_numpy(cd.GetArray("wallShearStress")),
        "p": vtk_to_numpy(cd.GetArray("p")),
    }


def main() -> None:
    """Tabulate the comparison."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--stair", type=Path, required=True)
    ap.add_argument("--ib", type=Path, required=True)
    ap.add_argument("--bf", type=Path, required=True)
    ap.add_argument("--time", type=float, default=15.0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    setup = json.loads((args.bf / "case_setup.json").read_text())
    st = stair_faces(args.stair, args.time)
    _, ib = read_bed(args.ib, None)
    _, bf = read_bed(args.bf, None)
    # wallShearStress = -(nuEff dev2 grad U) . n_boundary (outward from the
    # fluid): the force of the wall on the fluid; the shear on the bed is
    # minus that
    tread = np.abs(st["n"][:, 2]) > 0.9
    riser = np.abs(st["n"][:, 0]) > 0.9
    LOG.info("%-10s %8s %8s %8s | force per plan length [N/m per m width]:"
             " %8s %8s %8s", "segment", "stair", "IB", "BF", "treads",
             "risers", "p risers")
    for name, (x0, x1) in segments(setup).items():
        m = (st["x"][:, 0] > x0) & (st["x"][:, 0] < x1)
        mt, mr = m & tread, m & riser
        t_tread = -RHO * st["tau"][mt, 0]
        mi = (ib["x"] > x0) & (ib["x"] < x1)
        mb = (bf["x"] > x0) & (bf["x"] < x1)
        width = st["A"][mt].sum() / max((x1 - x0), 1e-12) if mt.any() else 1
        f_t = (-RHO * st["tau"][mt, 0] * st["A"][mt]).sum() / (x1 - x0)
        f_r = (-RHO * st["tau"][mr, 0] * st["A"][mr]).sum() / (x1 - x0)
        # pressure force of the fluid on the riser, x component; the
        # kinematic pressure times rho; normals point out of the fluid
        f_p = (RHO * st["p"][mr] * st["n"][mr, 0] * st["A"][mr]).sum() \
            / (x1 - x0)
        LOG.info("%-10s %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f", name,
                 t_tread.mean() if mt.any() else np.nan, ib["taux"][mi].mean(),
                 bf["taux"][mb].mean(), f_t / width, f_r / width, f_p / width)


if __name__ == "__main__":
    main()
