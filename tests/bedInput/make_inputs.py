#!/usr/bin/env python3
"""Write the input data of the bedInput test (run once; the files are kept).

* inputs/boundaryData/bottom/points, .../0/zbSurvey: 1200 scattered survey
  points (fixed seed) in the domain plus a ring just outside it, with the
  analytic bed level z_b = 0.01 sin(pi x/0.1) cos(pi y/0.1) [m];
* inputs/bed.stl: the same bed as a triangulated surface (4 mm grid).

Usage: make_inputs.py
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import stl_tools


def bed(x, y):
    """Analytic bed level [m]."""
    return 0.01 * np.sin(math.pi * x / 0.1) * np.cos(math.pi * y / 0.1)


def header(cls, obj):
    return ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            f"    class       {cls};\n    object      {obj};\n}}\n\n")


rng = np.random.default_rng(1)
xy = rng.uniform(-0.1, 0.1, size=(1200, 2))
ring = np.linspace(-0.11, 0.11, 23)
frame = np.array([(s, t) for s in ring for t in (-0.11, 0.11)]
                 + [(t, s) for s in ring[1:-1] for t in (-0.11, 0.11)])
xy = np.vstack([xy, frame])
z = bed(xy[:, 0], xy[:, 1])
d = os.path.join(HERE, "inputs", "boundaryData", "bottom")
os.makedirs(os.path.join(d, "0"), exist_ok=True)
with open(os.path.join(d, "points"), "w") as f:
    f.write(header("vectorField", "points") + f"{len(xy)}\n(\n")
    f.write("\n".join(f"({x:.8g} {y:.8g} 0)" for x, y in xy) + "\n)\n")
with open(os.path.join(d, "0", "zbSurvey"), "w") as f:
    f.write(header("scalarField", "zbSurvey") + f"{len(z)}\n(\n")
    f.write("\n".join(f"{v:.8g}" for v in z) + "\n)\n")
n = stl_tools.height_grid_stl(os.path.join(HERE, "inputs", "bed.stl"), "bed",
                              (-0.11, 0.11), (-0.11, 0.11), (0.004, 0.004),
                              bed, binary=True)
print(f"{len(xy)} survey points, {n} STL triangles")
