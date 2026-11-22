#!/usr/bin/env python3
"""Check one ``ibChannelFlat`` (T1) run and write ``metrics.json``.

Metrics (latest time directory unless ``--time`` is given):

* **Discharge balance**: Q_in (inlet), Q_out through fluid cells and
  through dead (below-bed) cells of the outlet [m^3/s];
  ``q_imbalance`` = |Q_out/Q_in - 1|.
* **Bed leakage**: flux through faces between fluid (Gamma >= 0.5) and
  dead (Gamma < 0) cells, net and gross (sum |phi|), relative to Q_in.
* **Velocity profile**: RMS of (u - u_ref)/U_b over live cells in the
  evaluation window, against the analytical laminar profile and/or the
  body-fitted (BF) reference run (interpolated in eta = (z - z_b)/h).
* **Bed shear stress** (kinematic, tau/rho [m^2/s^2]) per bed column:
  RAS: |ibWallShearStress| in IB cells (new IB library), or
  |shearVelocity|^2 (the legacy wall-model output);
  laminar: nu du/dz at the bed from a quadratic through u(z_b) = 0 and the
  first two fluid cells with centres more than 0.1 dz above the bed.
  Mean, coefficient of variation (std/mean, linear trend in (x, y)
  removed) and ratio to the reference
  (laminar analytic 3 nu U_b/h; BF ``wallShearStress`` on patch ``bed``).

Usage::

    check.py --run-dir RUN [--ref-dir BF_RUN] [--out metrics.json]

Tolerances are the targets for the new IB library (documented in
tests/README.md); the legacy solver is expected to fail some of them.
"""

import argparse
import logging
import math
import os
import sys
from typing import Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))

import case_setup  # noqa: E402
import channel_profiles as prof  # noqa: E402
import foam_file  # noqa: E402
import log_tools  # noqa: E402
from metrics import Metrics, build_info, setup_logging  # noqa: E402

LOG = logging.getLogger("ibChannelFlat.check")

TOL = {
    "q_imbalance": 1e-3,        # |Q_out/Q_in - 1|
    "leak_gross_rel": 1e-3,     # sum|phi_bed| / Q_in
    "u_rms_err": 0.05,          # RMS (u - u_ref)/U_b
    "tau_ratio_err": 0.10,      # |tau_mean/tau_ref - 1| (laminar)
    "tau_ratio_err_ras": 0.15,  # same, RAS
    "tau_cv": 0.05,             # std/mean of tau_b over the bed
}


class Result:
    """Fields of one run at one time (cell-centre arrays, SI units)."""

    def __init__(self, run_dir: str, time: Optional[str]) -> None:
        """Load mesh, set-up and fields of ``run_dir``."""
        self.run_dir = run_dir
        self.setup = case_setup.load_setup(run_dir)
        self.mesh = foam_file.PolyMesh(run_dir)
        times = foam_file.time_dirs(run_dir)
        times = [t for t in times if float(t) > 0]
        if not times:
            raise RuntimeError("no time directories in " + run_dir)
        self.time = time if time else times[-1]
        tdir = os.path.join(run_dir, self.time)
        nc = self.mesh.n_cells
        self.U = foam_file.read_internal_field(os.path.join(tdir, "U"), nc)
        self.phi_int = foam_file.read_internal_field(
            os.path.join(tdir, "phi"), self.mesh.n_internal)
        gpath = os.path.join(tdir, "Gamma")
        self.gamma = (foam_file.read_internal_field(gpath, nc)
                      if os.path.exists(gpath) else np.ones(nc))
        self.tdir = tdir
        cc = self.mesh.cell_centres
        tilt = math.tan(math.radians(self.setup["tilt_deg"]))
        self.zb = (self.setup["zb_centre"]
                   + (cc[:, 1] - 0.5 * self.setup["width"]) * tilt)
        self.depth = self.setup["z1"] - self.zb
        self.eta = (cc[:, 2] - self.zb) / self.depth

    def patch_phi(self, patch: str) -> np.ndarray:
        """Return phi [m^3/s] on the faces of ``patch``."""
        n = self.mesh.boundary[patch]["nFaces"]
        return foam_file.read_patch_field(os.path.join(self.tdir, "phi"),
                                         patch, n)

    def window(self) -> np.ndarray:
        """Return a mask of cells in the evaluation window (x range)."""
        x = self.mesh.cell_centres[:, 0] / self.setup["length"]
        lo, hi = (0.5, 0.9) if self.setup["mode"] == "laminar" else (0.6, 0.9)
        return (x > lo) & (x < hi)


def discharge(res: Result, met: Metrics) -> float:
    """Add discharge and bed-leakage metrics; return Q_in [m^3/s]."""
    mesh = res.mesh
    q_in = -float(np.sum(res.patch_phi("inlet")))
    phi_out = res.patch_phi("outlet")
    out_cells = mesh.owner[mesh.patch_slice("outlet")]
    dead = res.gamma[out_cells] < 0
    q_out_fluid = float(np.sum(phi_out[~dead]))
    q_out_dead = float(np.sum(phi_out[dead]))
    met.add("Q_in", q_in, "m3/s")
    met.add("Q_out_fluid", q_out_fluid, "m3/s")
    met.add("Q_out_dead", q_out_dead, "m3/s",
            note="flux leaving through dead (below-bed) outlet cells")
    met.add("q_imbalance", abs((q_out_fluid + q_out_dead) / q_in - 1.0), "-",
            TOL["q_imbalance"])
    if not np.any(res.gamma < 0):
        return q_in
    own = mesh.owner[:mesh.n_internal]
    nei = mesh.neighbour
    g_o, g_n = res.gamma[own], res.gamma[nei]
    fluid_dead = (g_o >= 0.5) & (g_n < 0)
    dead_fluid = (g_o < 0) & (g_n >= 0.5)
    # Flux oriented from the fluid side into the dead region
    into_dead = np.concatenate([res.phi_int[fluid_dead],
                                -res.phi_int[dead_fluid]])
    ib_ib = (g_o == 0.5) & (g_n == 0.5)
    met.add("leak_net_rel", float(np.sum(into_dead)) / q_in, "-",
            note="net flux fluid->dead across the bed / Q_in")
    met.add("leak_gross_rel", float(np.sum(np.abs(into_dead))) / q_in, "-",
            TOL["leak_gross_rel"])
    met.add("ib_ib_gross_rel", float(np.sum(np.abs(res.phi_int[ib_ib])))
            / q_in, "-", note="sum|phi| on IB-IB faces (masked in phiHbyA)")
    return q_in


def profile_error(res: Result, ref: Optional[Result],
                  met: Metrics) -> None:
    """Add RMS velocity-profile errors vs analytic and BF reference."""
    ub = res.setup["u_bulk"]
    sel = res.window() & (res.gamma == 1.0) & (res.eta > 0) & (res.eta < 1)
    u = res.U[sel, 0]
    eta = res.eta[sel]
    if res.setup["mode"] == "laminar":
        u_an = prof.laminar_profile(eta, np.ones_like(eta), ub)
        err = np.sqrt(np.mean(((u - u_an) / ub) ** 2))
        met.add("u_rms_err_analytic", err, "-", TOL["u_rms_err"])
    if ref is None:
        return
    rsel = ref.window() & (ref.eta > 0) & (ref.eta < 1)
    order = np.argsort(ref.eta[rsel])
    # Mean BF profile u(eta) over the window (columns are equivalent)
    r_eta = ref.eta[rsel][order]
    r_u = ref.U[rsel, 0][order]
    bins = np.unique(np.round(r_eta, 6))
    r_mean = np.array([r_u[np.round(r_eta, 6) == b].mean() for b in bins])
    u_ref = np.interp(eta, bins, r_mean)
    err = np.sqrt(np.mean(((u - u_ref) / ub) ** 2))
    met.add("u_rms_err_bf", err, "-", TOL["u_rms_err"])
    ib_sel = res.window() & (res.gamma == 0.5)
    if np.any(ib_sel):
        u_ib_ref = np.interp(res.eta[ib_sel], bins, r_mean)
        met.add("u_ib_cells_rel_err_mean",
                float(np.mean((res.U[ib_sel, 0] - u_ib_ref) / ub)), "-",
                note="mean (u_IB - u_BF)/U_b in IB cells")


def column_tau(res: Result) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-column (x, y, tau_b [m^2/s^2]) of the IB run."""
    cc = res.mesh.cell_centres
    win = res.window()
    if res.setup["mode"] == "ras":
        sel = win & (res.gamma == 0.5)
        # New IB library (ibScourPimpleFoam): kinematic wall shear stress
        path = os.path.join(res.tdir, "ibWallShearStress")
        if os.path.exists(path):
            tw = foam_file.read_internal_field(path, res.mesh.n_cells)
            tau = np.linalg.norm(tw[sel], axis=1)
            return cc[sel, 0], cc[sel, 1], tau
        # Legacy solver: friction velocity u_tau
        path = os.path.join(res.tdir, "shearVelocity")
        us = foam_file.read_internal_field(path, res.mesh.n_cells)
        tau = np.sum(us[sel] ** 2, axis=1)
        return cc[sel, 0], cc[sel, 1], tau
    nu = res.setup["nu"]
    # Fluid cells clearly above the bed (an IB cell centre can lie on it)
    sel = win & (res.gamma >= 0.5) & \
        (cc[:, 2] - res.zb > 0.1 * res.setup["dz"])
    keys = np.round(cc[sel, :2] / 1e-7).astype(np.int64)
    xs, ys, taus = [], [], []
    idx = np.where(sel)[0]
    for key in np.unique(keys, axis=0):
        col = idx[np.all(keys == key, axis=1)]
        col = col[np.argsort(cc[col, 2])][:2]
        y1, y2 = cc[col, 2] - res.zb[col]
        u1, u2 = res.U[col, 0]
        slope = (u1 * y2 ** 2 - u2 * y1 ** 2) / (y1 * y2 * (y2 - y1))
        xs.append(cc[col[0], 0])
        ys.append(cc[col[0], 1])
        taus.append(nu * slope)
    return np.array(xs), np.array(ys), np.array(taus)


def tau_stats(x: np.ndarray, y: np.ndarray, tau: np.ndarray,
              detrend: bool) -> Tuple[float, float]:
    """Return (mean, coefficient of variation) of tau, optionally detrended
    by a least-squares plane in (x, y)."""
    mean = float(np.mean(tau))
    resid = tau - mean
    if detrend and len(tau) > 3:
        amat = np.column_stack([np.ones_like(x), x, y])
        coef, _, _, _ = np.linalg.lstsq(amat, tau, rcond=None)
        resid = tau - amat.dot(coef)
    return mean, float(np.std(resid) / abs(mean)) if mean else float("nan")


def shear_stress(res: Result, ref: Optional[Result], met: Metrics) -> None:
    """Add bed-shear-stress metrics."""
    x, y, tau = column_tau(res)
    tilted = abs(res.setup["tilt_deg"]) > 0
    mean, cv = tau_stats(x, y, tau, detrend=True)
    met.add("tau_mean", mean, "m2/s2")
    met.add("tau_cv", cv, "-", TOL["tau_cv"],
            note="std/mean after removing a plane fit in (x, y)")
    met.add("tau_min_max", [float(tau.min()), float(tau.max())], "m2/s2")
    met.extra["tau_columns"] = {"x": x.tolist(), "y": y.tolist(),
                                "tau": tau.tolist()}
    ras = res.setup["mode"] == "ras"
    tol = TOL["tau_ratio_err_ras"] if ras else TOL["tau_ratio_err"]
    if not ras:
        depth = res.setup["z1"] - res.setup["zb_centre"]
        tau_an = 3.0 * res.setup["nu"] * res.setup["u_bulk"] / depth
        met.add("tau_ratio_err_analytic", abs(mean / tau_an - 1.0), "-", tol)
    if ref is None:
        return
    wss = foam_file.read_patch_field(os.path.join(ref.tdir, "wallShearStress"),
                                    "bed", ref.mesh.boundary["bed"]["nFaces"])
    fc = ref.mesh.face_centres[ref.mesh.patch_slice("bed")]
    lo, hi = (0.5, 0.9) if not ras else (0.6, 0.9)
    sel = (fc[:, 0] > lo * ref.setup["length"]) & \
        (fc[:, 0] < hi * ref.setup["length"])
    tau_bf = np.linalg.norm(wss[sel], axis=1)
    bf_mean, bf_cv = tau_stats(fc[sel, 0], fc[sel, 1], tau_bf, True)
    met.add("tau_mean_bf", bf_mean, "m2/s2")
    met.add("tau_cv_bf", bf_cv, "-")
    met.add("tau_ratio_err_bf", abs(mean / bf_mean - 1.0), "-", tol)
    if tilted:
        met.notes.append("tilted bed: CVs are computed after detrending")


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--ref-dir", default=None,
                        help="body-fitted reference run directory")
    parser.add_argument("--time", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--variant", default=None,
                        help="variant name (default: run directory name)")
    args = parser.parse_args()
    setup_logging()
    variant = (args.variant
               or os.path.basename(os.path.normpath(args.run_dir)))
    met = Metrics("ibChannelFlat", variant)
    met.extra["build"] = build_info(os.path.join(HERE, "..", ".."))
    out = args.out or os.path.join(args.run_dir, "metrics.json")
    logs = [f for f in os.listdir(args.run_dir) if f.startswith("log.")
            and f.endswith(("ibScourPimpleFoam", "pimpleFoam"))]
    if logs:
        summ = log_tools.summarise(log_tools.parse_log(
            os.path.join(args.run_dir, logs[0])))
        met.extra["log"] = summ
        for key in ("max_co", "max_co_ib", "s_per_step", "wall_s"):
            met.add(key, summ[key])
        if summ["fatal"]:
            met.error(summ["fatal"])
    try:
        res = Result(args.run_dir, args.time)
        ref = Result(args.ref_dir, None) if args.ref_dir else None
        met.extra["time"] = res.time
        discharge(res, met)
        if not res.setup["body_fitted"]:
            profile_error(res, ref, met)
            shear_stress(res, ref, met)
    except (OSError, RuntimeError, KeyError, ValueError) as err:
        LOG.error("check failed: %s", err)
        met.error("check failed: {}".format(err))
    met.write(out)
    print("{} {} {}".format(met.case, variant, met.status()))


if __name__ == "__main__":
    main()
