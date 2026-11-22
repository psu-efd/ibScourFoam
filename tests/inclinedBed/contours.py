#!/usr/bin/env python3
"""Contour plots of the flow around the pit: body-fitted vs IB runs.

Rows: one per run (e.g. body-fitted, original IB, cut-face IB); columns:
U_x, p, k, omega and nu_t/nu in the x-z plane around the pit. IB solid
cells (Gamma < -0.5) are masked; the bed z_b(x) of the case setup is
drawn on every panel. Colour ranges are common per column.

Usage::

    contours.py --time 15 --runs bf15:BF ib15z:IB ib15c2:IB-cut \
        --out fig.png [--xlim 0.55 1.05] [--zlim -0.032 0.02]
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib
import numpy as np
import matplotlib.tri as mtri

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_case import bed_elevation  # noqa: E402
from profiles import cell_data  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = logging.getLogger("inclinedBed.contours")

FIELDS = (
    ("U_x [m/s]", lambda d, nu: d["U"][:, 0], "RdBu_r", None),
    ("p [m2/s2]", lambda d, nu: d["p"], "viridis", None),
    ("k [m2/s2]", lambda d, nu: d["k"], "magma", None),
    ("log10 omega [1/s]",
     lambda d, nu: np.log10(np.maximum(d["omega"], 1e-3)), "magma", None),
    ("nu_t/nu", lambda d, nu: d["nut"] / nu, "magma", None),
)


def main() -> None:
    """Plot the contour panels."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True, help="dir:label")
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--xlim", type=float, nargs=2, default=(0.55, 1.05))
    ap.add_argument("--zlim", type=float, nargs=2, default=(-0.032, 0.02))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    runs = []
    for spec in args.runs:
        d, _, label = spec.partition(":")
        runs.append((Path(d), label or Path(d).name))
    setup = json.loads((runs[0][0] / "case_setup.json").read_text())
    nu = setup["nu"]

    data = []
    for run, label in runs:
        d = cell_data(run, args.time)
        C = d["C"]
        keep = ((C[:, 0] > args.xlim[0] - 0.01) & (C[:, 0] < args.xlim[1] + 0.01)
                & (C[:, 2] > args.zlim[0] - 0.002)
                & (C[:, 2] < args.zlim[1] + 0.002))
        if "Gamma" in d:
            keep &= d["Gamma"] > -0.5
        data.append({k: (v[keep] if hasattr(v, "__len__") and len(v) == len(C)
                         else v) for k, v in d.items()})
        LOG.info("%s: %d cells in the window", label, keep.sum())

    xb = np.linspace(*args.xlim, 400)
    zb = bed_elevation(xb, setup)
    nr, nc = len(runs), len(FIELDS)
    fig, ax = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.3 * nr + 0.6),
                           squeeze=False, sharex=True, sharey=True)
    for j, (name, fun, cmap, _) in enumerate(FIELDS):
        vals = [fun(d, nu) for d in data]
        lo = min(np.percentile(v, 1) for v in vals)
        hi = max(np.percentile(v, 99) for v in vals)
        if name.startswith("U_x"):
            m = max(abs(lo), abs(hi))
            lo, hi = -m, m
        levels = np.linspace(lo, hi, 41)
        for i, (d, (_, label)) in enumerate(zip(data, runs)):
            C = d["C"]
            tri = mtri.Triangulation(C[:, 0], C[:, 2])
            # drop triangles across the solid (long edges)
            xt, zt = C[tri.triangles, 0], C[tri.triangles, 2]
            edge = np.maximum(np.ptp(xt, axis=1), np.ptp(zt, axis=1))
            tri.set_mask(edge > 5 * setup["dx"])
            cs = ax[i, j].tricontourf(tri, np.clip(vals[i], lo, hi),
                                      levels=levels, cmap=cmap, extend="both")
            if name.startswith("U_x"):
                ax[i, j].tricontour(tri, vals[i], levels=[0], colors="k",
                                    linewidths=0.6)
            ax[i, j].plot(xb, zb, "k-", lw=1)
            ax[i, j].set_xlim(*args.xlim)
            ax[i, j].set_ylim(*args.zlim)
            if j == 0:
                ax[i, j].set_ylabel(f"{label}\nz [m]")
            if i == 0:
                ax[i, j].set_title(name)
            if i == nr - 1:
                ax[i, j].set_xlabel("x [m]")
        fig.colorbar(cs, ax=ax[:, j], orientation="horizontal",
                     fraction=0.05, pad=0.08)
    fig.suptitle(f"inclined bed {setup['slope_deg']:g} deg, t = {args.time:g} s"
                 " (black line U_x = 0)")
    fig.savefig(args.out, dpi=110, bbox_inches="tight")
    LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
