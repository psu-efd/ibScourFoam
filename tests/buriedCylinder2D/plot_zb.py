"""Plot the bed profiles z_b(x) of the buried-cylinder test at the written
times, with the cylinder, to inspect the activation and pass-through.

Usage: uv run --with matplotlib python plot_zb.py [case_dir] [out.png]
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_zb(path):
    """Return the internal field of an ascii finite-area scalar file."""
    text = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;", text, re.S)
    if m:
        return np.array([float(v) for v in m.group(1).split()])
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)", text)
    return np.array([float(m.group(1))]*80)


def main():
    case = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out = sys.argv[2] if len(sys.argv) > 2 else "zb_profiles.png"
    times = sorted(
        (float(d.name), d) for d in case.iterdir()
        if re.fullmatch(r"[\d.]+", d.name) and (d / "finite-area" / "zb").exists()
    )
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for t, d in times:
        zb = read_zb(d / "finite-area" / "zb")
        x = (np.arange(len(zb)) + 0.5)*0.4/len(zb)
        ax.plot(x, zb*1e3, label=f"t = {t:g} s")
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(0.2 + 0.02*np.cos(th), 20*np.sin(th), "k-", lw=1.5, label="cylinder")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z_b [mm]")
    ax.set_ylim(-85, 25)
    ax.legend(fontsize=7, ncol=3)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == "__main__":
    main()
