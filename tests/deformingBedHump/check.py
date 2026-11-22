#!/usr/bin/env python3
"""Check the body-fitted hump-migration run against the analytic celerity.

Reads the bed surfaces ``postProcessing/bed/<t>/bed.vtp`` (z_b, tauBed per
face), and reports

* the crest position and amplitude of z_b - z_b(control) over time,
* the measured celerity (linear fit of the crest position over the
  morphological time) against c = MF/(1-n) dq_b/dtheta 2 theta_0/H, with
  theta_0 the mean Shields number of the flat bed upstream of the hump,
* the bed volume change,

and plots z_b(x) and tau_b(x) at all output times.

Usage::

    check.py RUN_DIR [--control FLAT_RUN_DIR] [--out fig.png]
"""

import argparse
import glob
import json
import logging
import math
import os
from typing import Dict, Tuple

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = logging.getLogger("deformingBedHump.check")
G = 9.81


def read_bed(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return face-centre x, z_b and tau_x of a bed surface, sorted by x."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    pd = reader.GetOutput()
    cc = vtk.vtkCellCenters()
    cc.SetInputData(pd)
    cc.Update()
    x = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())[:, 0]
    zb = vtk_to_numpy(pd.GetCellData().GetArray("zb"))
    tau = vtk_to_numpy(pd.GetCellData().GetArray("tauBed"))[:, 0]
    o = np.argsort(x)
    return x[o], zb[o], tau[o]


def series(run: str) -> Dict[float, Tuple[np.ndarray, ...]]:
    """Return {time: (x, z_b, tau_x)} for all bed surfaces of a run."""
    out = {}
    for f in glob.glob(os.path.join(run, "postProcessing", "bed", "*",
                                    "bed.vtp")):
        out[float(os.path.basename(os.path.dirname(f)))] = read_bed(f)
    return dict(sorted(out.items()))


def dqdtheta(theta: float, theta_c: float) -> float:
    """Derivative of the Engelund-Fredsoe q* with respect to theta."""
    return 18.74 * (math.sqrt(theta) - 0.7 * math.sqrt(theta_c)
                    + (theta - theta_c) / (2 * math.sqrt(theta)))


def crest(x: np.ndarray, z: np.ndarray) -> Tuple[float, float]:
    """Crest position (parabolic fit through the maximum) and height."""
    i = int(np.argmax(z))
    if 0 < i < len(z) - 1:
        a, b, c = np.polyfit(x[i - 1:i + 2], z[i - 1:i + 2], 2)
        if a < 0:
            xc = -b / (2 * a)
            return float(xc), float(np.polyval([a, b, c], xc))
    return float(x[i]), float(z[i])


def main() -> None:
    """Evaluate one run (optionally relative to a flat-bed control)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run")
    ap.add_argument("--control", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    with open(os.path.join(args.run, "case_setup.json")) as src:
        setup = json.load(src)
    mf = setup["morph_factor"]
    t0 = setup["spinup_time"]
    rgd = (setup["s"] - 1) * G * setup["d50"]
    q_scale = math.sqrt(rgd) * setup["d50"]

    data = series(args.run)
    ctrl = series(args.control) if args.control else {}

    # Flat-bed Shields number: of the control run over the hump (the flow
    # may still develop along the channel), else upstream of the hump
    t_first = min(t for t in data if t >= t0)
    if t_first in ctrl:
        x, zb, tau = ctrl[t_first]
        near = np.abs(x - setup["x_hump"]) < setup["sigma"]
    else:
        x, zb, tau = data[t_first]
        near = (x > setup["fixed_inlet"] + 0.05) & \
            (x < setup["x_hump"] - 3 * setup["sigma"])
    theta0 = float(np.mean(tau[near])) / rgd
    c_ana = mf / (1 - setup["porosity"]) * q_scale \
        * dqdtheta(theta0, setup["theta_c0"]) * 2 * theta0 / setup["depth"]
    LOG.info("theta_0 = %.4f (flat bed); analytic celerity "
             "%.4g m/s per flow second (%.4g m/s morphological)",
             theta0, c_ana, c_ana / mf)

    rows = []
    for t, (x, zb, tau) in data.items():
        z = zb - ctrl[t][1] if t in ctrl else zb
        xc, a = crest(x, z)
        dx = np.gradient(x)
        rows.append((t, xc, a, float(np.sum(zb * dx))))
    rows = np.array(rows)
    LOG.info("   t [s]   t_m [s]  x_crest [m]  A [m]   bed volume [m2]")
    for t, xc, a, v in rows:
        LOG.info("%7.2f %8.1f %10.4f %9.5f %12.6e", t,
                 max(t - t0, 0) * mf, xc, a, v)

    m = rows[:, 0] >= t0 + 1.0
    if m.sum() >= 3:
        c_num = np.polyfit(rows[m, 0], rows[m, 1], 1)[0]
        LOG.info("celerity (per flow s): numerical %.4g m/s, analytic %.4g"
                 " m/s, ratio %.3f", c_num, c_ana, c_num / c_ana)
        LOG.info("amplitude %.5f -> %.5f m (%.1f %%)", rows[m, 2][0],
                 rows[-1, 2], 100 * rows[-1, 2] / setup["amplitude"])

    if args.out:
        fig, axs = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
        cmap = plt.get_cmap("viridis")
        times = list(data)
        for k, t in enumerate(times):
            x, zb, tau = data[t]
            col = cmap(k / max(len(times) - 1, 1))
            axs[0].plot(x, 1e3 * zb, color=col, lw=1)
            axs[1].plot(x, tau / rgd, color=col, lw=1)
        axs[0].set_ylabel("z_b [mm]")
        axs[1].set_ylabel("theta_x [-]")
        axs[1].axhline(setup["theta_c0"], color="k", ls=":", lw=0.8)
        axs[1].set_xlabel("x [m]")
        axs[0].set_title("Body-fitted hump, MF %g: t = %g ... %g s"
                         % (mf, times[0], times[-1]))
        fig.tight_layout()
        fig.savefig(args.out, dpi=130)
        LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
