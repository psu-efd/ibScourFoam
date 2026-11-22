#!/usr/bin/env python3
"""fsChannelIB: depth, bed shear and discharge of the uniform flow at the last written time."""
import re
import sys

import numpy as np

NX, NZ, DX, DZ = 150, 42, 0.02, 0.005
Q, S, G, N_MANNING = 0.05, 8e-4, 9.81, 0.0122
H_NORMAL = (Q*N_MANNING/np.sqrt(S))**0.6

with open("log.ibScourInterFoam") as f:
    log = f.read()
if "End\n" not in log:
    print("FAIL: the run did not finish")
    print(log[-1500:])
    sys.exit(1)
tEnd = re.findall(r"^Time = ([\d.eE+-]+)$", log, re.MULTILINE)[-1]


def field(name, vec=False):
    with open(f"{tEnd}/{name}") as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*\n?(\d+)\s*\n?\(", s)
    body = s[m.end():]
    body = body[:body.index("\n)")]
    if vec:
        return np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body), float).reshape(NZ, NX, 3)
    return np.array(body.split(), float).reshape(NZ, NX)


a = field("alpha.water"); U = field("U", True); tau = field("ibWallShearStress", True)
rows = slice(2, NZ)                       # fluid rows above the immersed bed
x = (np.arange(NX) + 0.5)*DX
h = a[rows].sum(axis=0)*DZ                # depth per column
q = (a[rows]*U[rows, :, 0]).sum(axis=0)*DZ
tauIB = np.linalg.norm(tau[2], axis=1)    # the IB row is the first fluid row
mid = (x > 0.75) & (x < 2.25)
hMid, hSlope = h[mid].mean(), (h[(x > 2.0) & (x < 2.5)].mean() - h[(x > 0.5) & (x < 1.0)].mean())/h[mid].mean()
tauMid, tauRef = tauIB[mid].mean(), G*h[mid].mean()*S
qErr = np.abs(q[mid] - Q).max()/Q
print(f"t = {tEnd} s: depth mid-reach {hMid:.4f} m vs normal depth {H_NORMAL:.4f} ({(hMid/H_NORMAL-1)*100:+.1f} %), "
      f"depth change over the middle half {hSlope*100:+.2f} %; inlet {h[:5].mean():.4f}, outlet {h[-5:].mean():.4f} m")
print(f"bed shear (IB row) mid-reach {tauMid:.3e} vs g h S {tauRef:.3e} m2/s2 ({(tauMid/tauRef-1)*100:+.1f} %); "
      f"discharge error max {qErr*100:.2f} %")
ok = abs(hMid/H_NORMAL - 1) < 0.03 and abs(hSlope) < 0.02 and abs(tauMid/tauRef - 1) < 0.10 and qErr < 0.02
print("PASS" if ok else "FAIL")
