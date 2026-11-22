#!/usr/bin/env python3
"""Check a ``humpMigration`` (T3) run and write ``metrics.json``.

Inputs: the bed level ``<time>/finite-area/zb`` at every output time of the
mobile-bed phase (after ``spinup_time``) and the kinematic wall shear
stress ``ibWallShearStress`` of the IB cells at the end of the spin-up (to
obtain the flat-bed Shields number). With ``--ref-dir`` the span-mean bed of
a flat-bed control run is subtracted at each output time, so the metrics
describe the hump alone.

Metrics:

* ``theta0``: flat-bed Shields number upstream of the hump [-];
* ``celerity`` [m/s] of the hump centroid x_c(t) = int x z+ dx / int z+ dx
  (linear fit in time) and ``celerity_ratio`` to the analytical value
  (see ``analytic_celerity``; REPORT: the estimate leaves out the slope correction
  of the threshold and the phase lag of the shear);
* ``log_volume_error_max``: largest Exner volume error of the bed report
  [m^3] (sediment conservation of the bed update);
* ``volume_change_rel``: max |V(t) - V(0)| / V_hump over the interior
  0.15 m < x < 1.05 m [-] (REPORT: the hump lowers the shear back to the
  inlet, so more of the equilibrium feed deposits than in the control);
* ``amplitude_ratio_final``: max z_b / A at the end [-];
* ``spike_rel``: max |z - (z[i-1] + z[i+1])/2| / A of the span-averaged
  profile (grid-scale oscillations) [-];
* ``span_spread_rel``: max of (max_y z - min_y z) / A over the hump, the
  columns where the bed (less the control) is above ``HUMP_LEVEL`` A at
  that time (a quasi-2D hump should stay uniform across the span) [-];
* ``span_spread_wake_rel``: the same over the rest of the window, mainly
  the scour trough behind the front (REPORT: the flow separates at the
  steep front and the wake varies across the narrow slip-sided span, so
  its spread depends on the start-up and peaks transiently at 0.04-0.05 A);
* ``flat_noise_rel``: std of z_b / A on the flat bed 0.15 m < x < 0.25 m
  after removing a linear trend (grid-scale noise; the raw std, which
  includes the deposition ramp near the inlet, is ``flat_std_rel``).

Usage::

    check.py --run-dir RUN [--ref-dir FLAT_RUN] [--out metrics.json]
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
sys.path.insert(0, HERE)

import bed_field  # noqa: E402
import case_setup  # noqa: E402
import foam_file  # noqa: E402
from metrics import Metrics, build_info, setup_logging  # noqa: E402

# numpy >= 2.0 renamed trapz to trapezoid
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

LOG = logging.getLogger("humpMigration.check")

TOL = {
    "log_volume_error": 1e-12,
    "spike_rel": 0.05,
    "span_spread_rel": 0.05,
    "flat_noise_rel": 0.02,
}
X_LO, X_HI = 0.15, 1.05   # interior evaluation window [m]
HUMP_LEVEL = 0.1          # hump columns: z_b - z_control > HUMP_LEVEL A


def bedload_derivative(theta: float, setup: Dict) -> float:
    """Return dq_b/dtheta [m^2/s] of the Roulund (2005) bed-load law."""
    scale = math.sqrt((setup["s"] - 1.0) * 9.81 * setup["d50"] ** 3)
    tc = setup["theta_c0"]
    return 18.74 * scale * (math.sqrt(theta) - 0.7 * math.sqrt(tc)
                            + (theta - tc) / (2.0 * math.sqrt(theta)))


def analytic_celerity(theta0: float, setup: Dict) -> float:
    """Return the hump celerity [m/s of flow time] incl. the morphological factor."""
    return (setup["morph_factor"] / (1.0 - setup["porosity"])
            * bedload_derivative(theta0, setup) * 2.0 * theta0
            / setup["depth"])


def profiles(run: str, setup: Dict) -> Tuple[np.ndarray, List[float],
                                            np.ndarray, np.ndarray]:
    """Return x [m], morphological-phase times [s], span-mean and span-range z [m].

    The first profile is the initial bed (t_m = 0); the others are the
    outputs after the spin-up.
    """
    geo, times, z = bed_field.bed_history(run)
    ux, inv = bed_field.column_index(geo["xc"][:, 0])
    t0 = setup["spinup_time"]
    # the initial bed: a zb file at t = 0, else the last output of the
    # spin-up (the bed is fixed until morphStartTime)
    before = [k for k, t in enumerate(times) if t <= t0 + 1e-9]
    after = [k for k, t in enumerate(times) if t > t0 + 1e-9]
    sel = ([before[-1]] if before else []) + after
    mean = np.zeros((len(sel), len(ux)))
    spread = np.zeros((len(sel), len(ux)))
    for row, k in enumerate(sel):
        for col in range(len(ux)):
            v = z[k, inv == col]
            mean[row, col] = v.mean()
            spread[row, col] = v.max() - v.min()
    tm = [max(times[k] - t0, 0.0) for k in sel]
    return ux, tm, mean, spread


def flat_bed_theta(run: str, setup: Dict) -> float:
    """Return the mean Shields number on the flat bed after the spin-up."""
    tname = "{:g}".format(setup["spinup_time"])
    mesh = foam_file.PolyMesh(run)
    tdir = os.path.join(run, tname)
    tw = foam_file.read_internal_field(os.path.join(tdir, "ibWallShearStress"),
                                      mesh.n_cells)
    gamma = foam_file.read_internal_field(os.path.join(tdir, "Gamma"),
                                         mesh.n_cells)
    x = mesh.cell_centres[:, 0]
    sel = (gamma == 0.5) & (x > 0.1) & (x < 0.25)
    tau = np.mean(np.linalg.norm(tw[sel], axis=1))
    return float(tau / ((setup["s"] - 1.0) * 9.81 * setup["d50"]))


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--ref-dir", default=None,
                        help="flat-bed control run, subtracted from z_b")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    setup_logging()
    run = args.run_dir
    variant = os.path.basename(os.path.normpath(run))
    met = Metrics("humpMigration", variant)
    met.extra["build"] = build_info(os.path.join(HERE, "..", ".."))
    setup = case_setup.load_setup(run)
    out = args.out or os.path.join(run, "metrics.json")
    log_path = os.path.join(run, "log.ibScourPimpleFoam")
    text = ""
    if os.path.isfile(log_path):
        with open(log_path, errors="replace") as fh:
            text = fh.read()
    if "FOAM FATAL" in text or "\nEnd" not in text:
        met.error("solver did not finish (see log.ibScourPimpleFoam)")

    amp = setup["amplitude"]
    v_hump = amp * setup["sigma"] * np.sqrt(2 * np.pi)  # per unit width
    x, times, z, spread = profiles(run, setup)
    if args.ref_dir:
        xr, tr, zr, _ = profiles(args.ref_dir, setup)
        ref = {round(t, 6): zr[k] for k, t in enumerate(tr)}
        for k, t in enumerate(times):
            if t > 0 and round(t, 6) in ref:
                z[k] = z[k] - ref[round(t, 6)]
            elif t > 0:
                LOG.warning("no control bed at t = %g", t)
        met.notes.append("bed change of the flat control run subtracted")
    if len(times) < 3:
        met.error(met.run_error or "fewer than 3 bed snapshots")
        met.write(out)
        return
    win = (x > X_LO) & (x < X_HI)
    flat = (x > 0.15) & (x < 0.25)
    cents, vols = [], []
    for zt in z:
        zp = np.clip(zt[win], 0.0, None)
        cents.append(float(np.sum(x[win] * zp) / np.sum(zp)))
        vols.append(float(_trapz(zt[win], x[win])))
    resid = z[-1, 1:-1] - 0.5 * (z[-1, :-2] + z[-1, 2:])
    met.extra["history"] = {"t_morph": times, "centroid": cents, "volume": vols}
    vol_err = max(abs(v - vols[0]) for v in vols) / v_hump
    celerity = float(np.polyfit(times[1:], cents[1:], 1)[0])
    lines = bed_field.log_bed_lines(log_path)
    err = [abs(float(line.split("volume error ")[1].split()[0]))
           for line in lines if "volume error " in line]
    met.add("log_volume_error_max", max(err) if err else float("nan"), "m3",
            TOL["log_volume_error"])
    met.add("volume_change_rel", vol_err, "-",
            note="includes the extra deposition of the inlet feed")
    met.add("celerity", celerity, "m/s", note="flow-time units (x MF)")
    try:
        theta0 = flat_bed_theta(run, setup)
        c_an = analytic_celerity(theta0, setup)
        met.add("theta0", theta0, "-")
        met.add("celerity_analytic", c_an, "m/s")
        met.add("celerity_ratio", celerity / c_an, "-",
                note="small-amplitude quasi-steady estimate")
    except (OSError, KeyError, ValueError) as err:
        LOG.warning("no spin-up wall shear stress: %s", err)
    met.add("amplitude_ratio_final", float(np.max(z[-1, win]) / amp), "-")
    met.add("spike_rel", float(np.max(np.abs(resid[win[1:-1]])) / amp), "-",
            TOL["spike_rel"])
    hump = (z > HUMP_LEVEL * amp) & win
    wake = win & ~hump
    met.add("span_spread_rel", float(np.max(spread[hump])) / amp, "-",
            TOL["span_spread_rel"], note="hump columns only")
    met.add("span_spread_wake_rel", float(np.max(spread[wake])) / amp, "-",
            note="outside the hump (wake trough)")
    zf = z[-1, flat]
    trend = np.polyval(np.polyfit(x[flat], zf, 1), x[flat])
    met.add("flat_std_rel", float(np.std(zf) / amp), "-")
    met.add("flat_noise_rel", float(np.std(zf - trend) / amp), "-",
            TOL["flat_noise_rel"])
    met.write(out)
    sys.stdout.write("{} {} {}\n".format(met.case, variant, met.status()))


if __name__ == "__main__":
    main()
