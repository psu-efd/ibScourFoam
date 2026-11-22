#!/usr/bin/env python3
"""fsDamBreakMobile: sediment confined to the water, water volume under the moving bed, dry bed untouched."""
import re
import sys

import numpy as np

NX, NZ, DX = 200, 90, 0.01

with open("log.ibScourInterFoam") as f:
    log = f.read()
if "End\n" not in log:
    print("FAIL: the run did not finish")
    print(log[-1500:])
    sys.exit(1)
times = sorted({float(t) for t in re.findall(r"^Time = ([\d.eE+-]+)$", log, re.MULTILINE)})
writes = ["0.1", "0.2", "0.5", "0.7"]


def field(t, name):
    with open(f"{t}/{name}") as f:
        s = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", s)
    body = s[m.end():]
    return np.array(body[:body.index("\n)")].split(), float).reshape(NZ, NX)


lost = sum(float(v) for v in re.findall(r"buried cells moved [-\d.eE+]+ m3 of water up, lost ([-\d.eE+]+) m3", log))
moved = sum(float(v) for v in re.findall(r"buried cells moved ([-\d.eE+]+) m3 of water up", log))
dry = [int(v) for v in re.findall(r"Bed: (\d+) dry IB cells", log)]
vol = [float(v) for v in re.findall(r"Phase-1 volume fraction = ([-\d.eE+]+)", log)]
worstAir = 0.0
worstBound = 0.0
water = []
V0 = 0.5*0.4*DX            # initial water column times the slab width 0.01 m
for t in writes:
    a = field(t, "alpha.water"); c = field(t, "C"); d = field(t, "ibDistance")
    fluid = d > 0
    air = fluid & (a < 1e-3)
    mostlyAir = fluid & (a < 0.5)
    inWater = (a*c)[fluid].sum()*DX*DX*0.01
    cAir = c[air].max() if air.any() else 0.0
    cWater = c[fluid & (a >= 0.5)].max()
    bound = (c[mostlyAir].max() if mostlyAir.any() else 0.0)/max(cWater, 1e-30)
    worstAir = max(worstAir, cAir); worstBound = max(worstBound, bound)
    water.append(a[fluid].sum()*DX*DX*0.01)
    print(f"t = {t} s: sediment {inWater:.3e} m3 (alpha C summed), C max in water cells {cWater:.2e}, "
          f"in cells with alpha < 0.5 {bound*cWater:.2e} (ratio {bound:.2f}), in air cells (alpha < 1e-3) {cAir:.1e}; water {water[-1]:.6e} m3")
# mass balance: the sediment in suspension equals the integrated net entrainment
dts = [float(v) for v in re.findall(r"^deltaT = ([\d.eE+-]+)", log, re.MULTILINE)]
ed = [(float(e), float(dd)) for e, dd in re.findall(r"entrainment ([\d.eE+-]+) deposition ([\d.eE+-]+) m3/s", log)]
net = sum((e - dd)*dt for (e, dd), dt in zip(ed, dts))
sus = [float(v) for v in re.findall(r"Suspended load: C min/max [-\d.eE+]+ / [-\d.eE+]+, volume ([-\d.eE+]+) m3", log)]
balance = abs(sus[-1] - net)/max(net, 1e-30)
print(f"sediment mass balance: in suspension at the end {sus[-1]:.3e} m3, integrated net entrainment {net:.3e} m3 (mismatch {balance:.2e})")
drift = (water[-1] - water[0])/V0
print(f"water volume: initial box {V0:.6e} m3, first write {water[0]:.6e}, last {water[-1]:.6e} (drift {drift:.2e}); "
      f"moved up under the bed {moved:.3e} m3, lost {lost:.3e} m3")
print(f"dry IB cells at the first bed step {dry[0] if dry else 'n/a'}, at the last {dry[-1] if dry else 'n/a'}")
# worstBound: the concentration of the half-air cells may equal the water
# maximum to round-off (1.00007 with the symmetric IB velocity
# reconstruction, 2026-11-01); it must not exceed it
ok = worstAir == 0.0 and worstBound <= 1.0 + 1e-3 and balance < 0.10 and abs(drift) < 1e-3 and abs(water[0]/V0 - 1) < 1e-3 and lost < 1e-6 and (dry and dry[0] > 100)
print("PASS" if ok else "FAIL")
