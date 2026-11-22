#!/usr/bin/env python3
"""fsWaveIB: wave height and surface elevation at four gauges, immersed bed against the body-fitted run."""
import sys

import numpy as np

T_WAVE, H_WAVE, T0 = 1.2, 0.02, 9.0       # period, target height, start of the comparison window


def gauges(case):
    d = np.loadtxt(f"{case}/postProcessing/gauges/0/height.dat", comments="#")
    t = d[:, 0]
    eta = d[:, 2::2]                        # hL: interface height above the gauge location (still water level)
    return t, eta


try:
    tI, eI = gauges("ib")
    tB, eB = gauges("bodyFitted")
except OSError as e:
    print("FAIL:", e)
    sys.exit(1)
tEnd = min(tI[-1], tB[-1])
if tEnd < T0 + 4*T_WAVE:
    print(f"FAIL: runs end at {tEnd:.2f} s, window needs {T0 + 4*T_WAVE:.1f} s")
    sys.exit(1)
tc = np.arange(T0, tEnd, 0.002)
worstH = worstEta = 0.0
for g, x in enumerate((2, 4, 6, 8)):
    ei = np.interp(tc, tI, eI[:, g]); eb = np.interp(tc, tB, eB[:, g])
    Hi, Hb = ei.max() - ei.min(), eb.max() - eb.min()
    dH = abs(Hi - Hb)/Hb
    dEta = np.abs(ei - eb).mean()/Hb
    worstH, worstEta = max(worstH, dH), max(worstEta, dEta)
    print(f"gauge x = {x} m: wave height ib {Hi*1e3:.2f} mm, body-fitted {Hb*1e3:.2f} mm ({dH*100:.2f} %); "
          f"mean |eta_ib - eta_bf| {dEta*100:.2f} % of H over {T0}-{tEnd:.1f} s")
print(f"body-fitted wave height at the first gauge {(eB[tB > T0, 0].max() - eB[tB > T0, 0].min())*1e3:.2f} mm (target {H_WAVE*1e3:.0f} mm)")
ok = worstH < 0.02 and worstEta < 0.05
print("PASS" if ok else "FAIL")
