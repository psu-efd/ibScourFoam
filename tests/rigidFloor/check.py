#!/usr/bin/env python3
"""Checks of the rigidFloor test (see Allrun).

Reads z_b and z_r at every written time (serial or reconstructed case) and the
bed lines of the log. The initial bed is the rule of the expressions in
constant/immersedBoundaryDict: sand 10 mm thick for 0.10 < x < 0.32 m with the
rigid sill z_r = 0.01 for 0.30 < x < 0.32 m.

Usage: check.py CASE
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import bed_field  # noqa: E402


def main(case: str) -> int:
    geo, times, z = bed_field.bed_history(case)
    _, times_r, zr_all = bed_field.bed_history(case, name="zRigid")
    if not times or not times_r:
        print("no zb/zRigid written: FAIL")
        return 1
    x = geo["xc"][:, 0]
    area = geo["area"]
    zr = zr_all[-1]
    z0 = np.where((x > 0.10) & (x < 0.32), 0.01, 0.0)
    zr0 = np.where((x > 0.30) & (x < 0.32), 0.01, 0.0)
    vol0 = float(np.sum((z0 - zr0) * area))
    worst_below = float(np.max(zr - z))
    sill = zr > 0
    sill_drop = float(np.max(zr[sill] - z[:, sill]))
    vol1 = float(np.sum((z[-1] - zr) * area))
    stripped = int(np.sum((z0 > zr0) & (z[-1] - zr < 1e-9)))
    with open(os.path.join(case, "log.ibScourPimpleFoam")) as fh:
        lines = [line for line in fh if line.startswith("Bed: MF")]
    err = max(abs(float(re.search(r"volume error ([-\d.eE+]+)", line).group(1)))
              for line in lines)
    clamp = float(re.search(r"rigid-floor clamp [-\d.eE+]+ m3 \(cumulative ([-\d.eE+]+)\)",
                            lines[-1]).group(1))
    rel = abs(vol1 / vol0 - 1)
    zr_err = float(np.max(np.abs(zr - zr0)))
    ok = (worst_below <= 1e-12 and sill_drop <= 1e-12 and rel < 1e-12
          and abs(clamp) < 1e-15 and err < 1e-12 and stripped > 0 and zr_err < 1e-14)
    print(f"times {times[0]:g}..{times[-1]:g}: z_r - z_b max {worst_below:.1e}, sill drop "
          f"{sill_drop:.1e}, z_r vs input {zr_err:.1e}, sand volume {vol0:.6e} -> {vol1:.6e} "
          f"(rel. {rel:.1e}), clamp {clamp:.1e}, max Exner volume error {err:.1e}, "
          f"patch faces stripped to the floor {stripped}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
