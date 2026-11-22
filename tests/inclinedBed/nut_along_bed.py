"""Along the bed: wall shear and eddy viscosity of the first cell, IB vs BF.

Plots tau_x of both runs, the first-cell nu_t/nu and the SST ratio
nu_t*omega/k (1 = no strain limiter, < 1 = limiter active) against x.
"""
import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "../../../ibScourFoam/tests/inclinedBed")
from profiles import cell_data
import compare as cmp  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(message)s")

LOG = logging.getLogger("nutAlongBed")
NU = 1e-6


def parse_args() -> argparse.Namespace:
    """Command-line arguments."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ib", type=Path, required=True)
    ap.add_argument("--bf", type=Path, required=True)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--xlim", type=float, nargs=2, default=(0.56, 0.66))
    return ap.parse_args()
ARGS = parse_args()
ib_run, bf_run, t, out = ARGS.ib, ARGS.bf, ARGS.time, ARGS.out

def first_cells(d: dict, mask_fn) -> tuple[np.ndarray, dict]:
    """Lowest fluid cell of each column (x rounded to 1 mm)."""
    C = d["C"]
    sel = mask_fn(d)
    idx = np.where(sel)[0]
    cols: dict[float, int] = {}
    for c in idx:
        key = round(C[c, 0], 4)
        if key not in cols or C[c, 2] < C[cols[key], 2]:
            cols[key] = c
    keys = np.array(sorted(cols))
    sel_c = np.array([cols[k] for k in keys])
    return keys, {k: d[k][sel_c] for k in ("nut", "omega", "k", "U", "C")}

ibd = cell_data(ib_run, t)
bfd = cell_data(bf_run, t)
x_ib, fc_ib = first_cells(ibd, lambda d: d["Gamma"] == 0.5)
x_bf, fc_bf = first_cells(bfd, lambda d: d["k"] > 1e-12)
_, bed_ib = cmp.read_bed(ib_run, None)
_, bed_bf = cmp.read_bed(bf_run, None)

fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
ax[0].plot(bed_bf["x"], bed_bf["taux"], "k-", label="body-fitted")
ax[0].plot(bed_ib["x"], bed_ib["taux"], "C1-", label="immersed boundary")
ax[0].set_ylabel("tau_x [Pa]"); ax[0].set_xlim(*ARGS.xlim); ax[0].set_ylim(-0.05, 0.9); ax[0].legend()
ax[0].axhline(0, color="0.6", lw=0.5)

ax[1].plot(x_bf, fc_bf["nut"] / NU, "k-", label="body-fitted, first cell")
ax[1].plot(x_ib, fc_ib["nut"] / NU, "C1-", label="IB wall cell (y+ ~ 11)")
ax[1].set_ylabel("nu_t/nu (first cell)"); ax[1].legend()

ax[2].plot(x_bf, fc_bf["nut"] * fc_bf["omega"] / fc_bf["k"], "k-")
ax[2].plot(x_ib, fc_ib["nut"] * fc_ib["omega"] / fc_ib["k"], "C1-")
ax[2].set_ylabel("nu_t*omega/k"); ax[2].set_xlabel("x [m]")
ax[2].axhline(1, color="0.6", lw=0.5)
fig.suptitle(f"15 deg inclined bed, t = {t:g} s: crest region")
fig.tight_layout(); fig.savefig(out, dpi=130)
for name, x, f in (("BF", x_bf, fc_bf), ("IB", x_ib, fc_ib)):
    for x0 in (0.56, 0.58, 0.59, 0.60, 0.61):
        i = np.argmin(abs(x - x0))
        LOG.info("%s x=%.3f nut/nu=%.2f nut*om/k=%.2f", name, x[i], f["nut"][i]/NU,
                 f["nut"][i]*f["omega"][i]/f["k"][i])
LOG.info("wrote %s", out)
