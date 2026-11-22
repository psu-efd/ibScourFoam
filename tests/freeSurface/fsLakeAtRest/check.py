#!/usr/bin/env python3
"""fsLakeAtRest: spurious velocity history, water volume and the hydrostatic pressure at the end."""
import re
import sys

import numpy as np

with open("log.ibScourInterFoam") as f:
    log = f.read()
if "End\n" not in log:
    print("FAIL: the run did not finish")
    print(log[-2000:])
    sys.exit(1)
u = [float(m) for m in re.findall(r"waterSpeed: max \|U\| in water = ([-\d.eE+]+)", log)]
vol = [float(m) for m in re.findall(r"Phase-1 volume fraction = ([-\d.eE+]+)", log)]
tEnd = re.findall(r"^Time = ([\d.eE+-]+)$", log, re.MULTILINE)[-1]


def field(name):
    with open(f"{tEnd}/{name}") as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", s)
    body = s[m.end():]
    return np.array(body[:body.index("\n)")].split(), float)


# hydrostatic pressure in the water: p = rho g (0.3 - z) with p = 0 at the surface
p = field("p"); a = field("alpha.water"); dist = field("ibDistance")
NX, DX = 100, 0.01
z = (np.arange(len(p))//NX + 0.5)*DX
wet = (a > 0.99) & (dist > 0) & (z < 0.28)
pErr = np.abs(p[wet] - 1000*9.81*(0.3 - z[wet])).max()/(1000*9.81*0.3)
print(f"max |U| in water: first {u[0]:.3e}, max {max(u):.3e}, last {u[-1]:.3e} m/s")
print(f"water volume fraction: first {vol[0]:.8f}, last {vol[-1]:.8f}, drift {vol[-1]-vol[0]:.2e}")
print(f"hydrostatic pressure error in the water at t = {tEnd}: {pErr:.2e} of rho g h")
ok = u[-1] < 1e-6 and max(u) < 1e-5 and abs(vol[-1] - vol[0]) < 1e-8 and pErr < 2e-3
print("PASS" if ok else "FAIL")
