#!/usr/bin/env python3
"""Check a ``pipelineScour`` run (half-buried cylinder) -> metrics.json.

Always (from the log): exit/crash, time reached, max IB Courant number,
masked IB continuity error, max "IB blocked-face leakage flux" (printed by
builds that report it), s/step, and the order of the first-anomaly times
of the indicators in common/anomalies.py.

Hydrodynamics (last field write):

* discharge imbalance |Q_out/Q_in - 1| (tolerance 1e-3) and gross flux
  through the bed (fluid-dead faces) / Q_in (tolerance 1e-3);
* tau_inf: mean IB bed shear stress (|ibWallShearStress|, m^2/s^2)
  upstream (-0.4 m < x < -0.2 m) and its coefficient of variation
  ``tau_cv_upstream`` (uniform bed, tolerance 0.05);
* ``tau_amp_max``: max tau_b / tau_inf within 3 D of the cylinder, with
  its x position, and ``tau_lee_max`` (max downstream, D/2 < x < 3 D)
  (REPORT).

Morphology (mobile-bed runs, span-averaged bed from
``<time>/finite-area/zb``):

* ``scour_depth`` S at the pipe flanks (D/2 < |x| < 2 D, max of -z_b) and
  S/D history; ``log_volume_error_max``: largest Exner volume error of the
  bed report (tolerance 1e-12 m^3); ``volume_change_rel`` = |change of the
  bed volume over -0.4 m < x < 0.9 m| / (S_final * D) (REPORT: includes the
  bed load through the ends of the window);
* ``spike_rel`` = max second-difference residual of z_b(x) / D
  (tolerance 0.02); ``span_spread_rel`` = max across-span range / D
  (tolerance 0.02); ``slope_max_deg`` (steepest span-mean slope).

Usage::

    check.py --run-dir RUN [--variant NAME]
"""

import argparse
import logging
import math
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import case_setup  # noqa: E402
from anomalies import first_anomalies  # noqa: E402
import foam_file  # noqa: E402
import log_tools  # noqa: E402
import bed_field  # noqa: E402
from metrics import Metrics, build_info, setup_logging  # noqa: E402

# numpy >= 2.0 renamed trapz to trapezoid
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

LOG = logging.getLogger("pipelineScour.check")

TOL = {"q_imbalance": 1e-3, "leak_gross_rel": 1e-3, "tau_cv_upstream": 0.05,
       "log_volume_error": 1e-12, "spike_rel": 0.02, "span_spread_rel": 0.02}


def hydro_metrics(run: str, setup: Dict, met: Metrics) -> None:
    """Discharge, bed leakage and bed shear stress of the last write."""
    times = [t for t in foam_file.time_dirs(run) if float(t) > 0
             and os.path.exists(os.path.join(run, t, "ibWallShearStress"))]
    if not times:
        met.notes.append("no field output")
        return
    tdir = os.path.join(run, times[-1])
    met.extra["field_time"] = times[-1]
    mesh = foam_file.PolyMesh(run)
    nc = mesh.n_cells
    gamma = foam_file.read_internal_field(os.path.join(tdir, "Gamma"), nc)
    phi = foam_file.read_internal_field(os.path.join(tdir, "phi"),
                                        mesh.n_internal)

    def patch_phi(name: str) -> np.ndarray:
        return foam_file.read_patch_field(os.path.join(tdir, "phi"), name,
                                          mesh.boundary[name]["nFaces"])

    q_in = -float(np.sum(patch_phi("inlet")))
    q_out = float(np.sum(patch_phi("outlet")))
    met.add("q_imbalance", abs(q_out / q_in - 1.0), "-", TOL["q_imbalance"])
    g_o = gamma[mesh.owner[:mesh.n_internal]]
    g_n = gamma[mesh.neighbour]
    bed = ((g_o >= 0.5) & (g_n < 0)) | ((g_o < 0) & (g_n >= 0.5))
    met.add("leak_gross_rel", float(np.sum(np.abs(phi[bed]))) / q_in, "-",
            TOL["leak_gross_rel"])
    tw = foam_file.read_internal_field(
        os.path.join(tdir, "ibWallShearStress"), nc)
    ib = gamma == 0.5
    x = mesh.cell_centres[ib, 0]
    tau = np.linalg.norm(tw[ib], axis=1)
    up = (x > -0.4) & (x < -0.2)
    tau_inf = float(np.mean(tau[up]))
    met.add("tau_inf", tau_inf, "m2/s2")
    met.add("tau_cv_upstream", float(np.std(tau[up]) / tau_inf), "-",
            TOL["tau_cv_upstream"])
    near = np.abs(x) < 3 * setup["diameter"]
    lee = (x > 0.5 * setup["diameter"]) & (x < 3 * setup["diameter"])
    if np.any(lee):
        met.add("tau_lee_max", float(tau[lee].max() / tau_inf), "-")
    if np.any(near):
        i = int(np.argmax(tau[near]))
        met.add("tau_amp_max", float(tau[near][i] / tau_inf), "-")
        met.add("tau_amp_max_x_over_D", float(x[near][i]
                                               / setup["diameter"]), "-")
    order = np.argsort(x)
    met.extra["tau_profile"] = {"x": x[order][::4].tolist(),
                                "tau": tau[order][::4].tolist()}


def span_profiles(run: str, y_lo: float, y_hi: float
                  ) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
    """Return times [s], x [m], span-mean z [m] and span range [m] (time, x)."""
    geo, times, z = bed_field.bed_history(run)
    xc = geo["xc"]
    keep = (xc[:, 1] > y_lo) & (xc[:, 1] < y_hi)
    ux, inv = bed_field.column_index(xc[keep, 0])
    zk = z[:, keep]
    mean = np.zeros((len(times), len(ux)))
    rng = np.zeros((len(times), len(ux)))
    for col in range(len(ux)):
        v = zk[:, inv == col]
        mean[:, col] = v.mean(axis=1)
        rng[:, col] = v.max(axis=1) - v.min(axis=1)
    return times, ux, mean, rng


def morph_metrics(run: str, setup: Dict, met: Metrics) -> None:
    """Scour depth history, volume and oscillations of the bed level."""
    d = setup["diameter"]
    if setup["dim"] == "3d":
        y_lo, y_hi = -0.25 * setup["cyl_length"], 0.25 * setup["cyl_length"]
    else:
        y_lo, y_hi = -1e-9, setup["width"] + 1e-9
    times, x, z, zr = span_profiles(run, y_lo, y_hi)
    if len(times) < 2:
        met.notes.append("no bed output")
        return
    win = (x > -0.4) & (x < 0.9)
    # pipe flanks; |x| < D/2 lies inside the (solid) cylinder
    under = (np.abs(x) > 0.5 * d) & (np.abs(x) < 2.0 * d)
    scour = (-z[:, under].min(axis=1)).tolist()
    volume = [float(_trapz(zt[win], x[win])) for zt in z]
    met.extra["scour_history"] = {"time": times, "scour": scour,
                                  "volume": volume}
    resid = z[-1, 1:-1] - 0.5 * (z[-1, :-2] + z[-1, 2:])
    slope = np.abs(np.diff(z[-1]) / np.diff(x))
    s_fin = scour[-1]
    met.add("scour_depth_final", s_fin, "m")
    met.add("scour_over_D_final", s_fin / d, "-")
    dvol = abs(volume[-1] - volume[0])
    met.add("volume_change_rel", dvol / max(s_fin * d, 1e-12), "-",
            note="bed volume change over -0.4..0.9 m / (S D); includes the "
            "bed load through the window ends")
    lines = bed_field.log_bed_lines(os.path.join(run, "log.ibScourPimpleFoam"))
    err = [abs(float(l.split("volume error ")[1].split()[0]))
           for l in lines if "volume error " in l]
    met.add("log_volume_error_max", max(err) if err else float("nan"), "m3",
            TOL["log_volume_error"], note="Exner volume balance per step")
    met.add("spike_rel", float(np.max(np.abs(resid[win[1:-1]])) / d), "-",
            TOL["spike_rel"])
    # a 2D bed must stay uniform across the span; in 3D it is not
    met.add("span_spread_rel", float(np.max(zr[:, win]) / d), "-",
            TOL["span_spread_rel"] if setup["dim"] == "2d" else None)
    met.add("slope_max_deg", math.degrees(math.atan(float(slope.max()))), "deg")


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--variant", default=None)
    args = parser.parse_args()
    setup_logging()
    run = args.run_dir
    variant = args.variant or os.path.basename(os.path.normpath(run))
    met = Metrics("pipelineScour", variant)
    met.extra["build"] = build_info(os.path.join(HERE, "..", ".."))
    setup = case_setup.load_setup(run)
    logs = sorted((os.path.join(run, f) for f in os.listdir(run)
                   if f.startswith("log.ibScourPimpleFoam")), key=os.path.getmtime)
    steps: List[Dict] = []
    fatal = None
    for path in logs:
        log = log_tools.parse_log(path)
        steps += log["steps"]
        fatal = fatal or log["fatal"]
    summ = log_tools.summarise({"steps": steps, "fatal": fatal,
                                "finished": not fatal})
    met.extra["log"] = summ
    met.add("t_end", summ["t_end"], "s")
    met.add("max_co_ib", summ["max_co_ib"], "-")
    met.add("max_cont_ib_local", summ["max_cont_ib_local"], "-")
    met.add("max_leak_sum_abs_phi", summ["max_leak_sum_abs_phi"], "m3/s",
            note="'IB blocked-face leakage flux' printed by the solver")
    met.add("s_per_step", summ["s_per_step"], "s")
    if fatal:
        met.error(fatal[:400])
    first = first_anomalies(steps, 0.5 * setup["diameter"] / 16)
    met.add("anomaly_order", [n for t, n in sorted(
        (t, n) for n, t in first.items() if t is not None)], "-")
    try:
        hydro_metrics(run, setup, met)
        if setup["mobile"]:
            morph_metrics(run, setup, met)
    except (OSError, KeyError, ValueError, RuntimeError) as err:
        LOG.error("check failed: %s", err)
        met.notes.append("check failed: {}".format(err))
    met.write(os.path.join(run, "metrics.json"))
    print("{} {} {}".format(met.case, variant, met.status()))


if __name__ == "__main__":
    main()
