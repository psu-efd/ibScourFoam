#!/usr/bin/env python3
"""Compare the bed shear stress along the bed of pipeline runs.

For each run, the kinematic bed shear stress |tau_b| [m^2/s^2] of the IB
cells (Gamma == 0.5), or of the faces of the wall patch ``bed`` of a
body-fitted run, is averaged over the span and over the time directories
in [t0, t1], then normalised by its upstream value tau_inf (mean over
x < x_inf). The profiles tau_b/tau_inf versus x/D are drawn in one figure.

Field per solver:

* new library (ibScourPimpleFoam): ``ibWallShearStress``;
* legacy (ibScourFoam): ``shearVelocity`` (tau_b = u_tau^2);
* body-fitted (e.g. pimpleFoam with the ``wallShearStress`` function
  object): patch ``bed`` of ``wallShearStress`` (detected from the mesh
  having a patch ``bed`` and no ``Gamma`` field).

The log reports per run tau_inf, the peaks of tau_b/tau_inf upstream
(x < -D/2) and downstream (x > D/2) of the pipe, the upstream separation
and downstream reattachment points (sign changes of the streamwise
component, sign fixed by the upstream flow) and, with ``--reference``,
the RMS difference of tau_b/tau_inf against the reference run over
``--rms-range``.

Requires numpy and matplotlib.

Usage::

    plot_bed_shear.py --run NAME=DIR [--run NAME=DIR ...] --diameter 0.05 \\
        --t0 1.5 --t1 3 [--reference NAME] --out bed_shear.png
"""

import argparse
import logging
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foam_file  # noqa: E402

logger = logging.getLogger(__name__)

# Categorical slots 1 and 2 of the reference palette, text inks, surface
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e4e3df"


def _column_mean(x: np.ndarray, values: np.ndarray
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """Average ``values`` over the span: bin by x (rounded to 0.1 um)."""
    ux, inv = np.unique(np.round(x, 7), return_inverse=True)
    count = np.bincount(inv)
    if values.ndim == 1:
        return ux, np.bincount(inv, weights=values) / count
    cols = [np.bincount(inv, weights=values[:, d]) / count
            for d in range(values.shape[1])]
    return ux, np.column_stack(cols)


def _ib_shear(case: str, tdir: str, n: int) -> Optional[np.ndarray]:
    """Return the IB-cell shear vector [m^2/s^2] (NaN elsewhere) or None."""
    gpath = os.path.join(tdir, "Gamma")
    if not os.path.exists(gpath):
        return None
    gamma = foam_file.read_internal_field(gpath, n)
    if os.path.exists(os.path.join(tdir, "ibWallShearStress")):
        tw = foam_file.read_internal_field(
            os.path.join(tdir, "ibWallShearStress"), n)
    elif os.path.exists(os.path.join(tdir, "shearVelocity")):
        us = foam_file.read_internal_field(
            os.path.join(tdir, "shearVelocity"), n)
        # tau = u_tau^2 along the direction of u_tau
        tw = us * np.linalg.norm(us, axis=1)[:, None]
    else:
        return None
    return np.where((gamma == 0.5)[:, None], tw, np.nan)


def bed_shear(case: str, t0: float, t1: float, patch: str = "bed"
              ) -> Tuple[np.ndarray, np.ndarray, int]:
    """Return x [m], time- and span-averaged shear vector, n_times.

    The shear vector is the kinematic bed shear stress [m^2/s^2] of the
    IB cells, or of the faces of the wall patch ``patch`` of a
    body-fitted run (``wallShearStress`` field).

    Args:
        case: Run directory (serial or reconstructed).
        t0: Start of the averaging window [s].
        t1: End of the averaging window [s].
        patch: Bed patch of a body-fitted run.
    """
    mesh = foam_file.PolyMesh(case)
    n = mesh.n_cells
    body_fitted = patch in mesh.boundary
    if body_fitted:
        faces = mesh.patch_slice(patch)
        x = mesh.face_centres[faces, 0]
        n_faces = mesh.boundary[patch]["nFaces"]
    else:
        x = mesh.cell_centres[:, 0]

    times = [t for t in foam_file.time_dirs(case) if t0 <= float(t) <= t1]
    acc = None
    used = 0
    for t in times:
        tdir = os.path.join(case, t)
        if body_fitted:
            path = os.path.join(tdir, "wallShearStress")
            if not os.path.exists(path):
                continue
            tw = foam_file.read_patch_field(path, patch, n_faces)
        else:
            tw = _ib_shear(case, tdir, n)
            if tw is None:
                continue
        acc = tw if acc is None else acc + tw
        used += 1

    if used == 0:
        raise RuntimeError(f"no usable time directories in {case}")

    tau_mean = acc / used
    sel = np.all(np.isfinite(tau_mean), axis=1)
    ux, tau_col = _column_mean(x[sel], tau_mean[sel])
    logger.info("%s: %s, %d bed columns", case,
                "body-fitted patch " + patch if body_fitted else "IB cells",
                len(ux))
    return ux, tau_col, used


def sign_changes(xd: np.ndarray, f: np.ndarray) -> List[float]:
    """Return the x/D where f changes sign (linear interpolation)."""
    out = []
    for i in np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]:
        out.append(xd[i] - f[i] * (xd[i + 1] - xd[i]) / (f[i + 1] - f[i]))
    return out


def local_maxima(xd: np.ndarray, ratio: np.ndarray, lo: float, hi: float,
                 floor: float = 0.05) -> List[Tuple[float, float]]:
    """Return (x/D, value) of the local maxima of ratio in lo < x/D < hi."""
    i = np.arange(1, len(xd) - 1)
    peak = (ratio[i] > ratio[i - 1]) & (ratio[i] >= ratio[i + 1])
    peak &= (xd[i] > lo) & (xd[i] < hi) & (ratio[i] > floor)
    return [(xd[k], ratio[k]) for k in i[peak]]


def profile_stats(xd: np.ndarray, ratio: np.ndarray, tau_x: np.ndarray,
                  xlim: Tuple[float, float]) -> Dict[str, object]:
    """Peaks of tau_b/tau_inf and separation/reattachment points.

    Args:
        xd: x/D of the bed columns.
        ratio: tau_b/tau_inf.
        tau_x: Streamwise shear, positive for downstream-directed shear.
        xlim: x/D window of the peak search.

    Returns:
        Local maxima upstream (xlim[0] < x/D < -0.5) and downstream
        (0.5 < x/D < xlim[1]) of the pipe, all sign changes of tau_x on
        either side, the most upstream one before the pipe (separation)
        and the last one behind it (reattachment); NaN where none exists.
    """
    up, down = xd < -0.5, xd > 0.5
    up_z = sign_changes(xd[up], tau_x[up])
    down_z = sign_changes(xd[down], tau_x[down])
    return {
        "up_peaks": local_maxima(xd, ratio, xlim[0], -0.5),
        "down_peaks": local_maxima(xd, ratio, 0.5, xlim[1]),
        "zeros": up_z + down_z,
        "separation_x": min(up_z) if up_z else np.nan,
        "reattachment_x": max(down_z) if down_z else np.nan,
    }


def _fmt_peaks(peaks: List[Tuple[float, float]]) -> str:
    """Format [(x/D, value), ...] for the log."""
    return ", ".join(f"{v:.3f} at {x:.2f}" for x, v in peaks) or "none"


def main() -> None:
    """Parse arguments, compute the profiles and write the figure."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", action="append", required=True,
                        help="NAME=DIR, repeatable")
    parser.add_argument("--diameter", type=float, required=True,
                        help="pipe diameter D [m]")
    parser.add_argument("--t0", type=float, default=0.0)
    parser.add_argument("--t1", type=float, default=1e30)
    parser.add_argument("--x-inf", type=float, default=-0.2,
                        help="upstream reference region x < x_inf [m]")
    parser.add_argument("--xlim", type=float, nargs=2, default=(-4, 8),
                        help="x/D range of the plot")
    parser.add_argument("--reference", default=None,
                        help="NAME of the run the others are compared to")
    parser.add_argument("--rms-range", type=float, nargs=2, default=(-4, 8),
                        help="x/D range of the RMS difference")
    parser.add_argument("--out", default="bed_shear.png")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 3.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    profiles: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for i, spec in enumerate(args.run):
        name, case = spec.split("=", 1)
        x, tau_vec, used = bed_shear(case, args.t0, args.t1)
        tau = np.linalg.norm(tau_vec, axis=1)
        upstream = x < args.x_inf
        tau_inf = np.mean(tau[upstream])
        # Streamwise shear, positive where it acts like the upstream one
        tau_x = tau_vec[:, 0] * np.sign(np.mean(tau_vec[upstream, 0]))
        xd = x / args.diameter
        ratio = tau / tau_inf
        st = profile_stats(xd, ratio, tau_x, args.xlim)
        logger.info("%s: %d times, tau_inf = %.4e m2/s2 (u_tau = %.4f m/s)",
                    name, used, tau_inf, np.sqrt(tau_inf))
        logger.info("  local maxima of tau_b/tau_inf (value at x/D): "
                    "upstream %s; downstream %s",
                    _fmt_peaks(st["up_peaks"]), _fmt_peaks(st["down_peaks"]))
        logger.info("  tau_x sign changes at x/D = %s; separation %.2f, "
                    "reattachment %.2f",
                    ", ".join(f"{z:.2f}" for z in st["zeros"]) or "none",
                    st["separation_x"], st["reattachment_x"])
        # Break the line across the pipe (no bed under the cylinder)
        gap = np.where(np.diff(xd) > 0.5)[0] + 1
        ax.plot(np.insert(xd, gap, np.nan), np.insert(ratio, gap, np.nan),
                color=SERIES[i], linewidth=2, label=name)
        profiles[name] = (xd, ratio)

    if args.reference is not None:
        xr, rr = profiles[args.reference]
        lo, hi = args.rms_range
        sel = (xr > lo) & (xr < hi)
        for name, (xd, ratio) in profiles.items():
            if name == args.reference:
                continue
            # Compare on the reference columns, outside the pipe footprint
            other = np.interp(xr[sel], xd, ratio, left=np.nan,
                              right=np.nan)
            ok = np.isfinite(other) & (np.abs(xr[sel]) > 0.5)
            rms = np.sqrt(np.mean((other[ok] - rr[sel][ok]) ** 2))
            logger.info("RMS(%s - %s) of tau_b/tau_inf over %g < x/D < %g:"
                        " %.3f", name, args.reference, lo, hi, rms)

    ax.axhline(1.0, color=INK_2, linewidth=0.8, linestyle=(0, (2, 3)))
    ax.set_xlim(*args.xlim)
    ax.set_xlabel("x / D", color=INK_2)
    ax.set_ylabel(r"$\tau_b / \tau_\infty$", color=INK_2)
    ax.set_title(
        f"Bed shear stress, time-averaged over t = {args.t0:g}-{args.t1:g} s",
        color=INK, fontsize=11, loc="left")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.tick_params(colors=INK_2, labelsize=9)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(args.out, facecolor=SURFACE)
    logger.info("wrote %s", args.out)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()
