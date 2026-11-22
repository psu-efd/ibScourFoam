#!/usr/bin/env python3
"""Velocity profiles normal to the up-ramp: immersed boundary vs body-fitted.

At stations x0 on the up-ramp of the inclined-bed test, prints and plots the
tangential velocity U.t, the eddy viscosity and (IB) the cell marker Gamma
and wall distance ibDistance of the cells in the column through x0, against
the wall-normal distance d = (z - z_b(x)) cos(alpha) from the ramp.

Usage::

    profiles.py --ib RUN --bf RUN --time T --x 0.90 0.94 [--out fig.png]
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_case import bed_elevation  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = logging.getLogger("inclinedBed.profiles")


def cell_data(run: Path, time: float) -> dict[str, np.ndarray]:
    """Cell centres and cell fields of the reconstructed run at ``time``."""
    foam = run / "case.foam"
    foam.touch()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(foam))
    reader.UpdateInformation()
    reader.DisableAllPatchArrays()
    reader.SetPatchArrayStatus("internalMesh", 1)
    reader.SetCreateCellToPoint(0)
    reader.UpdateTimeStep(time)
    reader.Update()
    mesh = reader.GetOutput().GetBlock(0)
    centres = vtk.vtkCellCenters()
    centres.SetInputData(mesh)
    centres.Update()
    out = {"C": vtk_to_numpy(centres.GetOutput().GetPoints().GetData())}
    cd = mesh.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        out[cd.GetArrayName(i)] = vtk_to_numpy(cd.GetArray(i))
    return out


def main() -> None:
    """Print and plot the profiles."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ib", type=Path, required=True)
    ap.add_argument("--bf", type=Path, required=True)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--x", type=float, nargs="+", required=True)
    ap.add_argument("--nmax", type=int, default=8, help="cells printed")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    s = json.loads((args.bf / "case_setup.json").read_text())
    a = math.radians(s["slope_deg"])
    t_hat = np.array([math.cos(a), 0.0, math.sin(a)])    # up-ramp
    data = {"IB": cell_data(args.ib, args.time),
            "BF": cell_data(args.bf, args.time)}
    fig, ax = plt.subplots(1, len(args.x), figsize=(5 * len(args.x), 5),
                           squeeze=False)
    for j, x0 in enumerate(args.x):
        zb0 = float(bed_elevation(np.array([x0]), s)[0])
        for name, d in data.items():
            C = d["C"]
            col = np.abs(C[:, 0] - x0) < 0.5 * s["dx"]
            k = np.where(col)[0]
            dist = (C[k, 2] - bed_elevation(C[k, 0], s)) * math.cos(a)
            o = np.argsort(dist)
            k, dist = k[o], dist[o]
            ut = d["U"][k] @ t_hat
            live = dist > -0.5 * s["dz"]
            if "Gamma" in d:
                live &= d["Gamma"][k] > -0.5
            LOG.info("%s x = %.3f (z_b %.2f mm):", name, x0, zb0 * 1e3)
            for i in np.where(live)[0][:args.nmax]:
                extra = ""
                if "Gamma" in d:
                    extra = (f" Gamma {d['Gamma'][k[i]]:+.1f} ibDist "
                             f"{d['ibDistance'][k[i]] * 1e3:6.3f} mm")
                LOG.info("   d %6.3f mm  U.t %7.4f  nut/nu %7.2f%s",
                         dist[i] * 1e3, ut[i], d["nut"][k[i]] / s["nu"],
                         extra)
            ax[0, j].plot(ut[live], dist[live] * 1e3,
                          "o-" if name == "IB" else "k.-", ms=3, label=name)
        ax[0, j].set_ylim(0, 15)
        ax[0, j].set_xlabel("U . t [m/s]")
        ax[0, j].set_ylabel("wall-normal distance [mm]")
        ax[0, j].set_title(f"slope {s['slope_deg']:g} deg, x = {x0:g} m")
        ax[0, j].legend()
    if args.out:
        fig.tight_layout()
        fig.savefig(args.out, dpi=110)
        LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
