#!/usr/bin/env python3
"""Verify the bed-load adaptation length on the body-fitted hump.

The bed core relaxes the bed-load magnitude to its equilibrium value along
the transport direction,

    L (e . grad) q = q_eq - q,

discretised with upwind inflow edges. On the one-cell-wide hump channel
(uniform cells dx, transport direction e_i = sign(tau_x,i)) the edge
coefficient is F = (e_i + e_(i+1))/2, and the upwind system is

    q_i (1 + a sum_in |F|) - a sum_in |F| q_N = q_eq,i,   a = L/dx,

which in a region of flow in +x is the recursion q_i (1 + a) - a q_(i-1)
= q_eq,i; there is no coupling across a separation point (F = 0) and
q = q_eq at the open inlet. The script recomputes q_eq from the written
bed shear stress and bed elevation (Engelund-Fredsoe with the slope-corrected
critical Shields number of Roulund et al. 2005, as in bedSediment), and
compares the written (adapted) flux with

  - the discrete recursion (tests the implementation: agreement to round-off
    of the written precision), and
  - the continuous solution q(x) = int q_eq(x') exp(-(x - x')/L) dx'/L
    (tests the model; first-order upwind differs by O(dx/L)).

Usage:
    check_adaptation.py <run> [--time T] [--out fig.png]
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

LOG = logging.getLogger("check_adaptation")


def read_bed(path: Path) -> dict[str, np.ndarray]:
    """Read a bed.vtp surface; return x of the face centres and fields."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    poly = reader.GetOutput()
    centres = vtk.vtkCellCenters()
    centres.SetInputData(poly)
    centres.Update()
    x = vtk_to_numpy(centres.GetOutput().GetPoints().GetData())[:, 0]
    cd = poly.GetCellData()
    out = {"x": x}
    for name in ("zb", "tauBed", "qBed"):
        out[name] = vtk_to_numpy(cd.GetArray(name))
    order = np.argsort(x)
    return {k: v[order] for k, v in out.items()}


def bed_dict(run: Path) -> dict[str, float]:
    """Sediment parameters of the deformingBed entry of system/controlDict."""
    text = (run / "system" / "controlDict").read_text()
    keys = {"grainDiameter": "d50", "specificGravity": "s",
            "criticalShields": "theta_c0", "frictionCoefficient": "mu_s",
            "adaptationLength": "L"}
    par: dict[str, float] = {"L": 0.0}
    for line in text.splitlines():
        words = line.replace(";", " ").split()
        if len(words) >= 2 and words[0] in keys:
            par[keys[words[0]]] = float(words[1])
    return par


def equilibrium_flux(bed: dict[str, np.ndarray], par: dict[str, float],
                     g: float = 9.81) -> np.ndarray:
    """Equilibrium bed load q_eq (signed, x) from tau_b and z_b."""
    x, z = bed["x"], bed["zb"]
    tau = bed["tauBed"][:, 0]
    sg = (par["s"] - 1) * g * par["d50"]
    theta = np.abs(tau) / sg
    # fac::grad (Gauss linear) on a uniform 1D strip: central differences,
    # zero gradient beyond the ends
    zp = np.concatenate(([z[0]], z, [z[-1]]))
    dzdx = (zp[2:] - zp[:-2]) / (x[1] - x[0]) / 2
    beta = np.minimum(np.arctan(np.abs(dzdx)), np.arctan(par["mu_s"]))
    t_hat = np.sign(tau)
    cos_phi = np.where(np.abs(dzdx) > 0, -t_hat * np.sign(dzdx), 0.0)
    theta_c = np.maximum(
        par["theta_c0"] * (np.cos(beta) - cos_phi * np.sin(beta) / par["mu_s"]),
        0.0)
    q_star = np.where(
        theta > theta_c,
        18.74 * (theta - theta_c) * (np.sqrt(theta) - 0.7 * np.sqrt(theta_c)),
        0.0)
    return q_star * np.sqrt(sg) * par["d50"] * t_hat


def discrete(q_eq: np.ndarray, e: np.ndarray, dx: float,
             L: float) -> np.ndarray:
    """Upwind discrete adapted magnitude, as in bedSediment (signed)."""
    a = L / dx
    m = np.abs(q_eq)
    n = len(m)
    F = 0.5 * (e[:-1] + e[1:])      # edge i | i+1, positive towards +x
    sum_in = np.zeros(n)
    sum_in[1:] += np.maximum(F, 0)
    sum_in[:-1] += np.maximum(-F, 0)
    b = m.copy()
    if e[0] > 0:                    # open inlet: equilibrium inflow
        sum_in[0] += 1
        b[0] += a * m[0]
    if e[-1] < 0:                   # open outlet, reversed flow
        sum_in[-1] += 1
        b[-1] += a * m[-1]
    q = m.copy()
    for _ in range(100 * n):
        rhs = b.copy()
        rhs[1:] += a * np.maximum(F, 0) * q[:-1]
        rhs[:-1] += a * np.maximum(-F, 0) * q[1:]
        new = rhs / (1 + a * sum_in)
        if np.max(np.abs(new - q)) < 1e-14 * max(m.max(), 1e-30):
            q = new
            break
        q = new
    return q * e


def adapted(q_eq: np.ndarray, dx: float, L: float) -> tuple[np.ndarray,
                                                            np.ndarray]:
    """Discrete (upwind) and continuous adapted flux for flow in +x."""
    a = L / dx
    q_d = np.empty_like(q_eq)
    prev = q_eq[0]                  # open inlet: equilibrium inflow
    for i, qe in enumerate(q_eq):
        q_d[i] = (qe + a * prev) / (1 + a)
        prev = q_d[i]
    # Continuous: exponential kernel, piecewise-constant q_eq per cell
    q_c = np.empty_like(q_eq)
    acc = q_eq[0]
    decay = np.exp(-dx / L)
    for i, qe in enumerate(q_eq):
        acc = acc * decay + qe * (1 - decay)
        q_c[i] = acc
    return q_d, q_c


def main() -> None:
    """Compare the written adapted flux with the discrete and exact ones."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--time", default=None, help="bed output time (latest)")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    base = args.run / "postProcessing" / "bed"
    times = sorted((p.name for p in base.iterdir()), key=float)
    t = args.time or times[-1]
    bed = read_bed(base / t / "bed.vtp")
    par = bed_dict(args.run)
    x = bed["x"]
    dx = float(np.mean(np.diff(x)))
    q = bed["qBed"][:, 0]
    q_eq = equilibrium_flux(bed, par)

    e = np.sign(bed["tauBed"][:, 0])
    rev = e < 0
    if par["L"] > 0:
        q_d = discrete(q_eq, e, dx, par["L"])
        _, q_c = adapted(q_eq, dx, par["L"])
    else:
        q_d = q_c = q_eq
    scale = np.max(np.abs(q_eq))
    # Interior only: the end faces see the one-sided gradient of the bed;
    # the continuous (+x) solution only upstream of the first reversal
    inner = np.zeros(len(x), bool)
    inner[2:-2] = True
    upstream = inner & (np.arange(len(x)) < (np.argmax(rev) if rev.any()
                                             else len(x)))
    LOG.info("faces with reversed shear: %d", int(rev.sum()))
    err_d = np.max(np.abs(q - q_d)[inner]) / scale
    err_c = np.max(np.abs(q - q_c)[upstream]) / scale
    err_eq = np.max(np.abs(q - q_eq)[inner]) / scale
    worst = int(np.argmax(np.where(inner, np.abs(q - q_d), 0)))
    LOG.info("  worst discrete mismatch at x = %.3f m (tau_x %.2e, "
             "z_b %.2f mm)", x[worst], bed["tauBed"][worst, 0],
             bed["zb"][worst] * 1e3)
    LOG.info("%s t = %s: L = %g m (L/dx = %.2f), max q_eq %.3e m2/s",
             args.run, t, par["L"], par["L"] / dx, scale)
    LOG.info("  max |q - q_discrete|/max q_eq   = %.2e", err_d)
    LOG.info("  max |q - q_continuous|/max q_eq = %.2e (upstream of the "
             "first reversal)", err_c)
    LOG.info("  max |q - q_eq|/max q_eq         = %.2e", err_eq)
    imax_eq, imax = int(np.argmax(q_eq)), int(np.argmax(q))
    LOG.info("  peak shift x(q max) - x(q_eq max) = %.4f m (L = %g m)",
             x[imax] - x[imax_eq], par["L"])

    if args.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(x, bed["zb"] * 1e3, "k-")
        ax[0].set_ylabel("z_b [mm]")
        ax[0].set_title(f"{args.run.name}, t = {t} s, L = {par['L']:g} m")
        ax[1].plot(x, q_eq, "C0-", lw=1, label="q_eq (from tau_b, z_b)")
        ax[1].plot(x, q, "C3-", lw=2, label="q written (adapted)")
        ax[1].plot(x, q_d, "k:", lw=1.5, label="discrete recursion")
        ax[1].plot(x, q_c, "C2--", lw=1, label="continuous solution")
        ax[1].set_xlabel("x [m]")
        ax[1].set_ylabel("q_x [m2/s]")
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.out, dpi=120)
        LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
