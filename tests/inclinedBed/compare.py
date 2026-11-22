#!/usr/bin/env python3
"""Compare the bed shear stress of the IB and body-fitted inclined-bed runs.

Reads ``postProcessing/bed/<t>/bed.vtp`` (face centres, z_b, tauBed) of an
immersed-boundary run and of its body-fitted reference at the same time,
plots z_b(x) and the streamwise bed shear stress tau_x(x) (the horizontal
shear stress used by the sediment core), and tabulates the mean tau_x per
segment of the pit (flat upstream, down-ramp, bottom, up-ramp, flat
downstream; two cells next to each kink excluded) with the IB/BF ratio.

Usage::

    compare.py --ib RUN --bf RUN [--time T] [--out fig.png]
"""
from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = logging.getLogger("inclinedBed.compare")
RHO = 1000.0  # [kg/m^3], kinematic -> Pa


def read_bed(run: Path, time: str | None) -> tuple[str, dict[str, np.ndarray]]:
    """Return the time and x, z_b, tau_x, |tau| of a run's bed.vtp."""
    base = run / "postProcessing" / "bed"
    times = sorted((p.name for p in base.iterdir()), key=float)
    t = time or times[-1]
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(base / t / "bed.vtp"))
    reader.Update()
    poly = reader.GetOutput()
    centres = vtk.vtkCellCenters()
    centres.SetInputData(poly)
    centres.Update()
    x = vtk_to_numpy(centres.GetOutput().GetPoints().GetData())[:, 0]
    cd = poly.GetCellData()
    zb = vtk_to_numpy(cd.GetArray("zb"))
    tau = vtk_to_numpy(cd.GetArray("tauBed"))
    o = np.argsort(x)
    return t, {"x": x[o], "zb": zb[o], "taux": tau[o, 0] * RHO,
               "tau": np.linalg.norm(tau[o], axis=1) * RHO}


def segments(setup: dict) -> dict[str, tuple[float, float]]:
    """x-ranges of the pit segments (two cells away from each kink)."""
    t = math.tan(math.radians(setup["slope_deg"]))
    d, e = setup["pit_depth"], 2 * setup["dx"]
    xa = setup["pit_start"]
    xb = xa + d / t
    xc = xb + setup["pit_bottom"]
    xd = xc + d / t
    return {"flat up": (xa - 0.2, xa - e), "down-ramp": (xa + e, xb - e),
            "bottom": (xb + e, xc - e), "up-ramp": (xc + e, xd - e),
            "flat down": (xd + e, xd + 0.2)}


def main() -> None:
    """Plot and tabulate."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ib", type=Path, required=True)
    ap.add_argument("--bf", type=Path, required=True)
    ap.add_argument("--time", default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    setup = json.loads((args.bf / "case_setup.json").read_text())
    t_ib, ib = read_bed(args.ib, args.time)
    t_bf, bf = read_bed(args.bf, args.time)
    LOG.info("slope %g deg: IB %s (t = %s), BF %s (t = %s)",
             setup["slope_deg"], args.ib, t_ib, args.bf, t_bf)
    LOG.info("%-10s %10s %10s %8s", "segment", "IB [Pa]", "BF [Pa]", "IB/BF")
    rows = {}
    for name, (x0, x1) in segments(setup).items():
        mi = (ib["x"] > x0) & (ib["x"] < x1)
        mb = (bf["x"] > x0) & (bf["x"] < x1)
        ti, tb = ib["taux"][mi].mean(), bf["taux"][mb].mean()
        rows[name] = (ti, tb, ti / tb if abs(tb) > 1e-12 else math.nan)
        LOG.info("%-10s %10.4f %10.4f %8.2f", name, *rows[name])

    if args.out:
        fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        ax[0].plot(bf["x"], bf["zb"] * 1e3, "k-", lw=1.5, label="BF bed")
        ax[0].plot(ib["x"], ib["zb"] * 1e3, "C1--", lw=1, label="IB z_b")
        ax[0].set_ylabel("z_b [mm]")
        ax[0].legend(fontsize=8)
        ax[0].set_title(f"inclined bed, slope {setup['slope_deg']:g} deg, "
                        f"dx {setup['dx'] * 1e3:g} mm, dz "
                        f"{setup['dz'] * 1e3:g} mm")
        ax[1].plot(bf["x"], bf["taux"], "k-", lw=1.5, label="body-fitted")
        ax[1].plot(ib["x"], ib["taux"], "C1-", lw=1, label="immersed boundary")
        ax[1].axhline(0, color="0.6", lw=0.5)
        for x0, x1 in segments(setup).values():
            ax[1].axvspan(x0, x1, color="0.93", zorder=0)
        ax[1].set_xlabel("x [m]")
        ax[1].set_ylabel("tau_x [Pa]")
        ax[1].legend(fontsize=8)
        ax[1].set_xlim(0.3, 1.4)
        fig.tight_layout()
        fig.savefig(args.out, dpi=120)
        LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
