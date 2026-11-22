#!/usr/bin/env python3
"""Check one ``sandSlideRelax`` (T2) run and write ``metrics.json``.

Reads the bed level ``<time>/finite-area/zb`` written at every step and the
``Bed: MF`` report lines of the solver log. Only the bed inside the domain
(one face away from its sides) is evaluated.

Metrics:

* ``slope_max_final``: maximum slope across a bed edge (between two
  neighbouring faces) at the final time; the sand slide limits exactly this
  slope (target <= 1.02 mu_s, the slide tolerance being 0.01 mu_s);
* ``steps_to_repose``: first step with that slope <= 1.02 mu_s;
* ``volume_err_rel``: max over time of |V(t) - V(0)| / V_cone, with
  V = sum of z_b A over the evaluation area [m^3];
* ``depth_err_rel``: |d_final/d' - 1| against the axisymmetric relaxed cone
  (case_setup.json; the edge-limited slide on a Cartesian bed gives a
  slightly square cone);
* ``n_slide_not_converged``: "Sand slide not converged" warnings;
* ``log_volume_error_max``: largest Exner volume error of the bed report [m^3].

Usage::

    check.py --run-dir RUN [--out metrics.json]
"""

import argparse
import logging
import os
import sys
from typing import Dict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import bed_field  # noqa: E402
import case_setup  # noqa: E402
from metrics import Metrics, build_info, setup_logging  # noqa: E402

LOG = logging.getLogger("sandSlideRelax.check")

TOL = {
    "slope_factor": 1.02,     # max edge slope <= factor * mu_s
    "volume_err_rel": 1e-6,   # |dV| / V_cone (the slide is conservative)
    "depth_err_rel": 0.10,    # |d/d' - 1|
    "log_volume_error": 1e-12,  # m^3
}


def edge_slopes(z: np.ndarray, i: np.ndarray, j: np.ndarray, dx: float,
                keep: np.ndarray) -> np.ndarray:
    """Return |dz/dx| across all x and y edges between kept faces of a grid."""
    nx, ny = i.max() + 1, j.max() + 1
    grid = np.full((nx, ny), np.nan)
    grid[i, j] = np.where(keep, z, np.nan)
    sx = np.abs(np.diff(grid, axis=0)) / dx
    sy = np.abs(np.diff(grid, axis=1)) / dx
    s = np.concatenate([sx.ravel(), sy.ravel()])
    return s[np.isfinite(s)]


def initial_cone(xc: np.ndarray, setup: Dict) -> np.ndarray:
    """Return the initial bed [m] at the face centres ``xc``: a cone of
    height ``cone_height`` and base ``radius`` (pit: negative)."""
    r = np.hypot(xc[:, 0], xc[:, 1])
    sign = -1.0 if setup["shape"] == "pit" else 1.0
    return sign * setup["cone_height"] * np.clip(1.0 - r / setup["radius"], 0.0, None)


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    setup_logging()
    run = args.run_dir
    out = args.out or os.path.join(run, "metrics.json")
    variant = os.path.basename(os.path.normpath(run))
    met = Metrics("sandSlideRelax", variant)
    met.extra["build"] = build_info(os.path.join(HERE, "..", ".."))
    setup = case_setup.load_setup(run)
    mu_s = setup["mu_s"]

    log_path = os.path.join(run, "log.ibScourPimpleFoam")
    text = open(log_path, errors="replace").read() if os.path.isfile(log_path) else ""
    if "FOAM FATAL" in text or "\nEnd" not in text:
        met.error("solver did not finish (see log.ibScourPimpleFoam)")
    lines = bed_field.log_bed_lines(log_path)
    met.add("n_slide_not_converged", text.count("Sand slide not converged"),
            "calls", 0)
    vol_err = [abs(float(l.split("volume error ")[1].split()[0]))
               for l in lines if "volume error " in l]
    met.add("log_volume_error_max", max(vol_err) if vol_err else float("nan"),
            "m3", TOL["log_volume_error"])

    geo, times, z = bed_field.bed_history(run)
    if times and times[0] > 0.0:
        # the bed is set by the initialElevation expression and no zb is
        # written at t = 0: prepend the initial cone at the face centres
        times = [0.0] + list(times)
        z = np.vstack([initial_cone(geo["xc"], setup), z])
    if len(times) < 2:
        met.error(met.run_error or "no bed output")
        met.write(out)
        return
    xc, area = geo["xc"], geo["area"]
    lim = setup["half_width"] - setup["cell"]
    keep = (np.abs(xc[:, 0]) < lim) & (np.abs(xc[:, 1]) < lim)
    i, j = bed_field.grid_index(xc)
    vol = (z[:, keep] * area[keep]).sum(axis=1)
    ana = setup["analytic"]
    dvol = np.abs(vol - vol[0]) / ana["volume"]
    smax = np.array([edge_slopes(zt, i, j, setup["cell"], keep).max() for zt in z])
    reached = np.nonzero(smax <= TOL["slope_factor"] * mu_s)[0]
    extreme = (-z[:, keep].min(axis=1) if setup["shape"] == "pit"
               else z[:, keep].max(axis=1))

    met.add("slope_max_initial", smax[0], "-")
    met.add("slope_max_final", smax[-1], "-", TOL["slope_factor"] * mu_s)
    met.add("steps_to_repose", int(reached[0]) if len(reached) else -1, "steps",
            note="-1: not reached")
    met.add("volume_err_rel", dvol.max(), "-", TOL["volume_err_rel"])
    met.add("depth_final", extreme[-1], "m")
    met.add("depth_analytic", ana["height"], "m")
    met.add("depth_err_rel", abs(extreme[-1] / ana["height"] - 1.0), "-",
            TOL["depth_err_rel"])
    met.extra["history"] = {"time": times, "slope_max": smax.tolist(),
                            "depth": extreme.tolist(),
                            "volume_err_rel": dvol.tolist()}
    met.write(out)
    LOG.info("%s %s %s", met.case, variant, met.status())
    sys.stdout.write("%s %s %s\n" % (met.case, variant, met.status()))


if __name__ == "__main__":
    main()
