#!/usr/bin/env python3
"""Laminar inclined-bed pit: IB vs body-fitted comparison (L3 gate).

Compares the immersed-boundary run (``ibScourPimpleFoam``) with the
body-fitted reference (``pimpleFoam``) of the same pit, using in each run the
differential operators of its own solver:

* wall-normal pressure gradient: IB ``ibGradP`` (written by the solver; the
  gradient of the momentum equation, blocked faces at the cell value), BF
  ``grad(p)`` (``postProcess -func 'grad(p)'``, the Gauss gradient that
  pimpleFoam uses). ``grad(p)`` of an IB run is never used: in the IB cells it
  takes the decoupled solid-cell pressure on the blocked faces;
* bed shear stress: ``tauBed`` of ``postProcessing/bed/<t>/bed.vtp`` (IB wall
  model and BF ``deformingBed``, same convention);
* velocity and pressure: cell values interpolated along the column to a
  fixed normal distance from the bed (no differential operator);
* bed pressure: the wall value of each solver's own condition, zero normal
  gradient at the bed (BF zeroGradient patch, IB ib.grad), i.e. the pressure
  of the first non-solid cell of the column. It drives the flow in the
  sediment (pore pressure) in coupled applications.

Pressure in both solvers is kinematic, p/rho [m2/s2]; the bed pressure is
plotted in Pa with --rho (default 1000 kg/m3), the other pressure
quantities stay kinematic (q = U_b^2/2 in m2/s2).

Reports, per bed region and write time, the signed tangential shear, the
IB/BF shear ratio, the normal pressure gradient, the near-bed pressure
difference relative to q = U_b^2/2, a sawtooth index of the shear (energy at
the two-cell scale) and the change between the last two writes (steadiness).

Usage::

    compare_laminar.py --ib RUN --bf RUN [--times 150 175 200] [--rho 1000]
                       [--out fig.png]
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_case import bed_elevation  # noqa: E402
from profiles import cell_data  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOG = logging.getLogger("inclinedBed.compare_laminar")

GAMMA_IB = 0.5      # Gamma marker of the IB (cut) cells
GAMMA_FLUID = 1.0   # Gamma marker of the fluid cells
D_SAMPLE = 2.0e-3   # normal distance of the near-bed samples [m]
D_NORMAL = 0.5e-3   # normal distance of the wall-normal velocity sample [m]


def bed_slope(x: np.ndarray, s: dict) -> np.ndarray:
    """Return dz_b/dx of the trapezoidal pit (central difference)."""
    h = 1.0e-5
    return (bed_elevation(x + h, s) - bed_elevation(x - h, s)) / (2.0 * h)


def bed_frame(x: np.ndarray, s: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return the unit bed tangent (downstream) and normal (into the flow)."""
    m = bed_slope(x, s)
    w = np.sqrt(1.0 + m ** 2)
    zero = np.zeros_like(m)
    t = np.stack([1.0 / w, zero, m / w], axis=1)
    n = np.stack([-m / w, zero, 1.0 / w], axis=1)
    return t, n


def regions(s: dict) -> list[tuple[str, float, float]]:
    """Return the named bed regions (x ranges) of the pit, corners excluded."""
    xa = s["pit_start"]
    run = s["pit_depth"] / np.tan(np.radians(s["slope_deg"]))
    xb, xc = xa + run, xa + run + s["pit_bottom"]
    xd = xc + run
    e = 4.0 * s["dx"]
    return [("flatUp", xa - 0.15, xa - e), ("downRamp", xa + e, xb - e),
            ("pitFloor", xb + e, xc - e), ("upRamp", xc + e, xd - e),
            ("flatDown", xd + e, xd + 0.15)]


def bed_shear(run: Path, time: str, s: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return x and the signed tangential bed shear of ``bed.vtp`` at ``time``."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(run / "postProcessing" / "bed" / time / "bed.vtp"))
    reader.Update()
    out = reader.GetOutput()
    arr = out.GetCellData().GetArray("tauBed")
    if arr is not None:
        cc = vtk.vtkCellCenters()
        cc.SetInputData(out)
        cc.Update()
        x = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())[:, 0]
    else:
        arr = out.GetPointData().GetArray("tauBed")
        x = vtk_to_numpy(out.GetPoints().GetData())[:, 0]
    tau = vtk_to_numpy(arr)
    t, _ = bed_frame(x, s)
    tt = np.einsum("ij,ij->i", tau, t) if tau.ndim == 2 else tau
    # Front and back points of the 2D bed share x: average the duplicates
    xu, inv = np.unique(np.round(x, 7), return_inverse=True)
    return xu, np.bincount(inv, tt) / np.bincount(inv)


def column_samples(d: dict, s: dict, xs: np.ndarray, grad: str,
                   ib: bool) -> dict[str, np.ndarray]:
    """Near-bed samples of each mesh column at the cell-centre x of ``xs``.

    Returns the normal pressure gradient of the first cell above the bed
    (IB cell for an IB run) and of the next cell (first full fluid cell),
    U.t and p interpolated to the normal distance ``D_SAMPLE``, the
    wall-normal velocity U.n at ``D_NORMAL`` (zero at the wall) and the bed
    pressure (first non-solid cell, zero normal gradient at the wall).
    """
    C = d["C"]
    t_all, n_all = bed_frame(xs, s)
    out = {k: np.full(len(xs), np.nan)
           for k in ("dpdn0", "dpdn1", "d0", "d1", "ut", "p", "un", "pbed")}
    col_of = np.round(C[:, 0] / s["dx"] - 0.5).astype(int)
    for i, x0 in enumerate(xs):
        col = np.where(col_of == int(round(x0 / s["dx"] - 0.5)))[0]
        zb = float(bed_elevation(np.array([x0]), s)[0])
        if ib:
            gam = d["Gamma"][col]
            col = col[(np.isclose(gam, GAMMA_IB)) | (np.isclose(gam, GAMMA_FLUID))]
        col = col[np.argsort(C[col, 2])]
        dn = (C[col, 2] - zb) * n_all[i, 2]
        if len(col) < 3:
            continue
        g = d[grad][col] @ n_all[i]
        out["dpdn0"][i], out["dpdn1"][i] = g[0], g[1]
        out["d0"][i], out["d1"][i] = dn[0], dn[1]
        ut = d["U"][col] @ t_all[i]
        out["ut"][i] = np.interp(D_SAMPLE, dn, ut)
        out["p"][i] = np.interp(D_SAMPLE, dn, d["p"][col])
        out["un"][i] = np.interp(D_NORMAL, dn, d["U"][col] @ n_all[i])
        out["pbed"][i] = d["p"][col[0]]
    return out


def sawtooth_index(tau: np.ndarray) -> float:
    """Two-cell-scale content of a shear profile: rms(second difference)/rms."""
    if len(tau) < 5:
        return float("nan")
    d2 = tau[2:] - 2.0 * tau[1:-1] + tau[:-2]
    return float(np.sqrt(np.mean(d2 ** 2)) / max(np.sqrt(np.mean(tau ** 2)), 1e-30))


def oscillation(g: np.ndarray) -> float:
    """Rms of a profile minus its 5-point running mean (two-cell content)."""
    g = g[np.isfinite(g)]
    if len(g) < 7:
        return float("nan")
    r = g - np.convolve(g, np.ones(5) / 5.0, mode="same")
    return float(np.sqrt(np.mean(r[2:-2] ** 2)))


def written_times(run: Path) -> list[str]:
    """Return the write times of ``postProcessing/bed`` (sorted)."""
    ts = [p.name for p in (run / "postProcessing" / "bed").iterdir()]
    return sorted(ts, key=float)


def main() -> None:
    """Compare the two runs and write the figure and the summary."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ib", type=Path, required=True)
    ap.add_argument("--bf", type=Path, required=True)
    ap.add_argument("--times", nargs="+", default=None,
                    help="write times for the shear history (default: all)")
    ap.add_argument("--rho", type=float, default=1000.0,
                    help="density for the bed pressure in Pa [kg/m3]")
    ap.add_argument("--out", type=Path, default=Path("compare_laminar.png"))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    s = json.loads((args.bf / "case_setup.json").read_text())
    q = 0.5 * s["u_bulk"] ** 2
    times = args.times or [t for t in written_times(args.ib)
                           if t in written_times(args.bf)]
    t_last = times[-1]
    regs = regions(s)

    # --- Bed shear: history and IB/BF per region
    shear = {}
    for t in times:
        shear[t] = (bed_shear(args.ib, t, s), bed_shear(args.bf, t, s))
    LOG.info("q = U_b^2/2 = %.3e Pa; times %s", q, " ".join(times))
    LOG.info("%-9s %6s %11s %11s %7s %8s %8s", "region", "t", "mean|tau|IB",
             "mean|tau|BF", "IB/BF", "sawIB", "sawBF")
    for name, lo, hi in regs:
        for t in times:
            (xi, ti), (xf, tf) = shear[t]
            mi, mf = (xi > lo) & (xi < hi), (xf > lo) & (xf < hi)
            ai, af = np.abs(ti[mi]).mean(), np.abs(tf[mf]).mean()
            LOG.info("%-9s %6s %11.3e %11.3e %7.3f %8.3f %8.3f", name, t, ai, af,
                     ai / af, sawtooth_index(ti[mi]), sawtooth_index(tf[mf]))
    if len(times) > 1:
        for k, lab in ((0, "IB"), (1, "BF")):
            (_, ta), (_, tb) = shear[times[-2]][k], shear[t_last][k]
            LOG.info("steadiness %s: max|tau(%s)-tau(%s)|/max|tau| = %.4f", lab,
                     t_last, times[-2], np.abs(tb - ta).max() / np.abs(tb).max())

    # --- Cell fields at the last time: dp/dn, near-bed U.t and p
    di = cell_data(args.ib, float(t_last))
    df = cell_data(args.bf, float(t_last))
    if "ibGradP" not in di:
        sys.exit(f"{args.ib}/{t_last}: no ibGradP (solver build without it); "
                 "grad(p) of an IB run is not used, see the module docstring")
    if "grad(p)" not in df:
        sys.exit(f"{args.bf}/{t_last}: run postProcess -func 'grad(p)' first")
    xs = s["dx"] * (np.arange(int(0.45 / s["dx"]), int(1.15 / s["dx"])) + 0.5)
    si = column_samples(di, s, xs, "ibGradP", ib=True)
    sf = column_samples(df, s, xs, "grad(p)", ib=False)
    LOG.info("\nt = %s, d = %.1f mm: %-9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s",
             t_last, 1e3 * D_SAMPLE, "region", "|dpdn|IB", "|dpdn|BF",
             "|dUt|/Ub", "|dp|/q", "chkIB", "|dUn|/Ub", "osc0IB", "osc1IB",
             "osc1BF", "|dpb|/q")
    for name, lo, hi in regs:
        m = (xs > lo) & (xs < hi)
        g0 = si["dpdn0"][m]
        # checkerboard: fraction of sign changes of the detrended IB dp/dn
        r = g0 - np.convolve(g0, np.ones(5) / 5.0, mode="same")
        chk = np.mean(np.signbit(r[2:-3]) != np.signbit(r[3:-2]))
        LOG.info("%33s %9.4f %9.4f %9.4f %9.4f %9.2f %9.4f %9.5f %9.5f %9.5f %9.4f",
                 name, np.nanmax(np.abs(g0)),
                 np.nanmax(np.abs(sf["dpdn0"][m])),
                 np.nanmax(np.abs(si["ut"][m] - sf["ut"][m])) / s["u_bulk"],
                 np.nanmax(np.abs(si["p"][m] - sf["p"][m])) / q, chk,
                 np.nanmax(np.abs(si["un"][m] - sf["un"][m])) / s["u_bulk"],
                 oscillation(g0), oscillation(si["dpdn1"][m]),
                 oscillation(sf["dpdn1"][m]),
                 np.nanmax(np.abs(si["pbed"][m] - sf["pbed"][m])) / q)

    # --- Figure
    fig, ax = plt.subplots(5, 1, figsize=(10, 16), sharex=True)
    cmap = plt.get_cmap("viridis")
    for k, t in enumerate(times):
        c = cmap(k / max(len(times) - 1, 1))
        (xi, ti), (xf, tf) = shear[t]
        lw = 1.4 if t == t_last else 0.7
        ax[0].plot(xi, ti, "-", color=c, lw=lw, label=f"IB t={t}")
        ax[0].plot(xf, tf, "--", color=c, lw=lw, label=f"BF t={t}")
    ax[0].axhline(0.0, color="k", lw=0.4)
    ax[0].set_ylabel(r"$\tau_b \cdot t$ [m$^2$/s$^2$] (kinematic)")
    ax[0].legend(ncol=len(times), fontsize=7)
    ax[0].set_title("(a) signed bed shear, IB solid / BF dashed")
    ax[1].plot(xs, si["dpdn0"], "C0", lw=1.0, label="IB cell (ibGradP)")
    ax[1].plot(xs, si["dpdn1"], "C2", lw=1.0, label="IB first fluid cell")
    ax[1].plot(xs, sf["dpdn0"], "C3--", lw=1.0, label="BF first cell")
    ax[1].plot(xs, sf["dpdn1"], "C1--", lw=1.0, label="BF second cell")
    ax[1].set_ylabel(r"$\nabla p \cdot n$ [m/s$^2$] (kinematic)")
    ax[1].legend(fontsize=8, ncol=2)
    ax[1].set_title(f"(b) wall-normal pressure gradient, t = {t_last}")
    ax[2].plot(xs, si["ut"], "C0", lw=1.0, label="IB")
    ax[2].plot(xs, sf["ut"], "C3--", lw=1.0, label="BF")
    ax[2].set_ylabel(r"$U \cdot t$ [m/s]")
    ax[2].legend(fontsize=8)
    ax[2].set_title(f"(c) tangential velocity {1e3 * D_SAMPLE:.0f} mm from the bed")
    rho = args.rho
    ax[3].plot(xs, rho * si["pbed"], "C0", lw=1.0, label="IB")
    ax[3].plot(xs, rho * sf["pbed"], "C3--", lw=1.0, label="BF")
    ax[3].set_ylabel(r"$\rho\, p_{bed}$ [Pa]")
    ax[3].legend(fontsize=8)
    ax[3].set_title(rf"(d) bed pressure ($\rho$ = {rho:g} kg/m$^3$; outlet p = 0)")
    ax[4].plot(xs, (si["pbed"] - sf["pbed"]) / q, "C0", lw=1.0, label="at the bed")
    ax[4].plot(xs, (si["p"] - sf["p"]) / q, "C2", lw=0.8,
               label=f"{1e3 * D_SAMPLE:.0f} mm from the bed")
    ax[4].set_ylabel(r"$(p_{IB} - p_{BF})/q$")
    ax[4].set_xlabel("x [m]")
    ax[4].legend(fontsize=8)
    ax[4].set_title(r"(e) pressure difference, $q = U_b^2/2$")
    for a in ax:
        for _, lo, hi in regs[1:4]:
            a.axvspan(lo, hi, color="0.92", zorder=0)
    ax[0].set_xlim(xs[0], xs[-1])
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
