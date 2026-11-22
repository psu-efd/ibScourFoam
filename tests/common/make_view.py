#!/usr/bin/env python3
"""Create a ParaView-friendly view of an ibScourFoam run directory.

The legacy solver writes the bed (``<time>/<object>_dual``) at
``extrudeMeshOutputTimeStep`` and the flow fields at ``writeInterval``. When
the bed is written more often, many time directories contain only the bed
region, and ParaView's OpenFOAM reader (which skips time 0 by default)
shows an almost empty field list.

This script builds ``<run>/view/`` with symbolic links to ``constant``,
``system`` and only those time directories that contain flow fields, plus a
``view.foam`` file to open in ParaView. The run directory itself is not
modified, so the view can be rebuilt at any time while a run is going.

Usage::

    python3 make_view.py <runDir> [<runDir> ...] [--field U]
"""

import argparse
import logging
import os
import shutil
from pathlib import Path
from typing import List

logger = logging.getLogger(__name__)


def is_time_dir(path: Path) -> bool:
    """Return True if ``path`` is a directory whose name is a number."""
    if not path.is_dir():
        return False
    try:
        float(path.name)
    except ValueError:
        return False
    return True


def field_times(run_dir: Path, field: str) -> List[Path]:
    """Return time directories of ``run_dir`` that contain ``field``.

    Args:
        run_dir: OpenFOAM case directory (serial run).
        field: Name of a flow field that marks a full write, e.g. ``U``.

    Returns:
        Time directories sorted by time value.
    """
    times = [
        p for p in run_dir.iterdir()
        if is_time_dir(p) and (p / field).exists()
    ]
    return sorted(times, key=lambda p: float(p.name))


def make_view(run_dir: Path, field: str) -> Path:
    """Build (or rebuild) ``run_dir/view`` and return its path."""
    view = run_dir / "view"
    if view.exists():
        shutil.rmtree(view)
    view.mkdir()

    # Mesh and settings are shared with the run
    for name in ("constant", "system"):
        os.symlink(Path("..") / name, view / name)

    times = field_times(run_dir, field)
    for t in times:
        os.symlink(Path("..") / t.name, view / t.name)

    (view / "view.foam").touch()
    logger.info(
        "%s: %d time directories with '%s' (%s)",
        view, len(times), field, ", ".join(t.name for t in times),
    )
    return view


def main() -> None:
    """Parse the command line and build a view for each run directory."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+", type=Path, help="run directories")
    parser.add_argument(
        "--field", default="U",
        help="field that marks a full write (default: U)",
    )
    args = parser.parse_args()

    for run in args.runs:
        if not (run / "system").is_dir():
            logger.warning("%s is not an OpenFOAM case; skipped", run)
            continue
        make_view(run.resolve(), args.field)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()
