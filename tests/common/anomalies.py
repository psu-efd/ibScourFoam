"""First-anomaly detection on parsed solver logs (divergence diagnosis).

For each indicator the first time step at which it leaves its normal
range is returned; the order of these times points to the mechanism of a
divergence:

* ``co``: IB Courant number > 2,
* ``cont``: masked IB continuity error (sum local) > 1e-4,
* ``p_iters``: p iterations > 3x the median of the first 20 steps,
* ``dz``: |bed change per update| (``dzPoint``) > 0.5 ``dz_cell``,
* ``tau``: max dual-mesh wall shear stress > 5x its median of the first
  20 steps,
* ``nan``: non-finite residual.
"""

from typing import Any, Dict, List, Optional

import numpy as np


def first_anomalies(steps: List[Dict[str, Any]],
                    dz_cell: float) -> Dict[str, Optional[float]]:
    """Return {indicator: first time [s] it is anomalous, or None}."""
    head = steps[:20]

    def median(key: str, sub: Optional[int] = None) -> float:
        vals = []
        for s in head:
            if key in s:
                vals.append(s[key][sub] if sub is not None else s[key])
        return float(np.median(vals)) if vals else float("nan")

    p_ref = median("p_iters")
    tau_ref = median("sed_wallShearStressMag", 1)
    tests = {
        "co": lambda s: s.get("co_ib_max", 0.0) > 2.0,
        "cont": lambda s: s.get("cont_ib_local", 0.0) > 1e-4,
        "p_iters": lambda s: s.get("p_iters", 0) > 3.0 * p_ref,
        "dz": lambda s: "dz_point" in s and max(
            abs(s["dz_point"][0]), abs(s["dz_point"][1])) > 0.5 * dz_cell,
        "tau": lambda s: "sed_wallShearStressMag" in s and
        s["sed_wallShearStressMag"][1] > 5.0 * tau_ref,
        "nan": lambda s: "nan_residual" in s,
    }
    first = {}  # type: Dict[str, Optional[float]]
    for name, test in tests.items():
        first[name] = None
        for st in steps:
            try:
                bad = test(st)
            except (TypeError, IndexError):
                bad = False
            if bad:
                first[name] = st["time"]
                break
    return first
