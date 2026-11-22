#!/usr/bin/env python3
"""Check one bedInput run: the written bed against the analytic level.

z_b = 0.01 sin(pi x/0.1) cos(pi y/0.1) [m] at the floor-face centres, at the
last write time (and, for ``rigid``, zRigid = z_b - 0.02 m). Exact inputs
(expression, rigid, faMesh, parallel) must agree to 1e-10 m; interpolated
ones (points, surface) within the linear-interpolation error of their data
(max 5e-4 m, rms 1e-4 m, i.e. 5 % and 1 % of the amplitude).

Usage: check.py --run-dir RUN --variant NAME
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import bed_field

TOL = {"exact": (1e-10, 1e-10), "interp": (5e-4, 1e-4)}  # (max, rms) [m]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--variant", required=True)
    a = ap.parse_args()
    log = os.path.join(a.run_dir, "log.ibScourPimpleFoam")
    text = ""
    if os.path.isfile(log):
        with open(log, errors="replace") as fh:
            text = fh.read()
    if "FOAM FATAL" in text or "\nEnd" not in text:
        print(f"bedInput {a.variant}: ERROR (solver did not finish, see {log})")
        return
    geo, times, z = bed_field.bed_history(a.run_dir)
    if not times:
        print(f"bedInput {a.variant}: ERROR (no zb written)")
        return
    x, y = geo["xc"][:, 0], geo["xc"][:, 1]
    exact = 0.01 * np.sin(math.pi * x / 0.1) * np.cos(math.pi * y / 0.1)
    kind = "interp" if a.variant in ("points", "surface") else "exact"
    tmax, trms = TOL[kind]
    msgs, ok = [], True
    fields = [("zb", z[-1], exact)]
    if a.variant == "rigid":
        _, tr, zr = bed_field.bed_history(a.run_dir, name="zRigid")
        fields.append(("zRigid", zr[-1] if len(tr) else np.full_like(exact, np.nan),
                       exact - 0.02))
    for name, val, ref in fields:
        err = val - ref
        emax, erms = float(np.max(np.abs(err))), float(np.sqrt(np.mean(err ** 2)))
        good = emax <= tmax and erms <= trms
        ok &= good
        msgs.append(f"{name} at t = {times[-1]:g}: max error {emax:.2e} m, rms {erms:.2e} m")
    if a.variant in ("faMesh", "parallel"):
        read = "finite-area mesh read from constant" in text
        ok &= read
        msgs.append("area mesh read from constant" if read else "area mesh NOT read from constant")
    if a.variant == "parallel":
        msgs.append(f"{len(times)} reconstructed bed times")
        ok &= len(times) >= 2
    print(f"bedInput {a.variant}: " + "; ".join(msgs) + (": PASS" if ok else ": FAIL"))


if __name__ == "__main__":
    main()
