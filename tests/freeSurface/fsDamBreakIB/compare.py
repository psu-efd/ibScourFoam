#!/usr/bin/env python3
"""fsDamBreakIB: front position, column height and the alpha field of the immersed-floor run against the body-fitted run."""
import re
import sys

import numpy as np

NX, DX = 117, 0.005


def alpha(case, t):
    with open(f"{case}/{t}/alpha.water") as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", s)
    body = s[m.end():]
    return np.array(body[:body.index("\n)")].split(), float).reshape(-1, NX)


worst_front = worst_height = worst_l1 = 0.0
for t in ("0.1", "0.2", "0.3", "0.4", "0.5"):
    try:
        aIB = alpha("ib", t)[2:]          # drop the two solid rows below the floor
        aBF = alpha("bodyFitted", t)
    except (OSError, AttributeError):
        print(f"t = {t}: fields missing")
        sys.exit(1)
    res = {}
    for lab, a in (("ib", aIB), ("bodyFitted", aBF)):
        wetRow = np.where(a[0] > 0.5)[0]
        front = (wetRow.max() + 1)*DX if wetRow.size else 0.0
        height = a[:, 0].sum()*DX
        res[lab] = (front, height)
    l1 = np.abs(aIB - aBF).sum()*DX*DX/(aBF.sum()*DX*DX)
    dFront = abs(res["ib"][0] - res["bodyFitted"][0])/0.585
    dHeight = abs(res["ib"][1] - res["bodyFitted"][1])/0.292
    worst_front, worst_height, worst_l1 = max(worst_front, dFront), max(worst_height, dHeight), max(worst_l1, l1)
    print(f"t = {t} s: front ib {res['ib'][0]:.3f} bf {res['bodyFitted'][0]:.3f} m ({dFront*100:.2f} % of the box), "
          f"column height ib {res['ib'][1]:.4f} bf {res['bodyFitted'][1]:.4f} m ({dHeight*100:.2f} %), "
          f"alpha L1 difference {l1*100:.2f} % of the water area")
ok = worst_front < 0.02 and worst_height < 0.02 and worst_l1 < 0.05
print("PASS" if ok else "FAIL")
