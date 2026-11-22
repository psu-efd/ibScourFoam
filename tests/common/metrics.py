"""Collect test metrics and write the ``metrics.json`` of a test case.

Each metric has a value, a unit, an optional tolerance and a status:

* ``PASS``/``FAIL`` when a tolerance is given (``value <= tol`` or
  ``value >= tol`` depending on ``op``),
* ``REPORT`` when the metric is informative only.

The overall status is ``FAIL`` if any metric fails, ``ERROR`` if the run
itself failed (set with :meth:`Metrics.error`), else ``PASS`` (or
``REPORT`` when nothing was checked).
"""

import json
import logging
import math
import os
import subprocess
import time
from typing import Any, Dict, Optional

LOG = logging.getLogger(__name__)


def _clean(value: Any) -> Any:
    """Convert numpy scalars and non-finite floats for JSON output."""
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, dict):
        return {k: _clean(v) for k, v in value.items()}
    return value


def build_info(repo_dir: str) -> Dict[str, str]:
    """Return git and library information identifying the tested build."""
    info = {}
    try:
        info["git_head"] = subprocess.check_output(
            ["git", "-C", repo_dir, "rev-parse", "--short", "HEAD"],
            universal_newlines=True).strip()
        info["git_branch"] = subprocess.check_output(
            ["git", "-C", repo_dir, "rev-parse", "--abbrev-ref", "HEAD"],
            universal_newlines=True).strip()
        dirty = subprocess.check_output(
            ["git", "-C", repo_dir, "status", "--short", "--", "src",
             "solvers"], universal_newlines=True).strip()
        info["src_dirty_files"] = str(len(dirty.splitlines()))
    except (OSError, subprocess.CalledProcessError) as err:
        LOG.warning("git information unavailable: %s", err)
    libdir = os.environ.get("FOAM_USER_LIBBIN", "")
    lib = os.path.join(libdir, "libIBMethod.so")
    if os.path.exists(lib):
        info["libIBMethod_mtime"] = time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(lib)))
    return info


class Metrics:
    """Container for the metrics of one test (case + variant)."""

    def __init__(self, case: str, variant: str) -> None:
        """Create an empty record for ``case``/``variant``."""
        self.case = case
        self.variant = variant
        self.items = {}  # type: Dict[str, Dict[str, Any]]
        self.notes = []  # type: list
        self.run_error = None  # type: Optional[str]
        self.extra = {}  # type: Dict[str, Any]

    def add(self, name: str, value: Any, unit: str = "-",
            tol: Optional[float] = None, op: str = "<=",
            note: str = "") -> None:
        """Add a metric; with ``tol`` it is checked as ``value op tol``."""
        status = "REPORT"
        if tol is not None:
            try:
                val = float(value)
                ok = val <= tol if op == "<=" else val >= tol
                status = "PASS" if (ok and math.isfinite(val)) else "FAIL"
            except (TypeError, ValueError):
                status = "FAIL"
        self.items[name] = {"value": _clean(value), "unit": unit,
                            "tolerance": tol, "op": op if tol is not None
                            else None, "status": status, "note": note}

    def error(self, message: str) -> None:
        """Mark the run as failed (crash, divergence, missing output)."""
        self.run_error = message

    def status(self) -> str:
        """Return the overall status string."""
        if self.run_error:
            return "ERROR"
        stats = [m["status"] for m in self.items.values()]
        if "FAIL" in stats:
            return "FAIL"
        return "PASS" if "PASS" in stats else "REPORT"

    def write(self, path: str) -> None:
        """Write the JSON record to ``path``."""
        record = {
            "case": self.case,
            "variant": self.variant,
            "status": self.status(),
            "run_error": self.run_error,
            "metrics": self.items,
            "notes": self.notes,
            "extra": _clean(self.extra),
        }
        with open(path, "w") as out:
            json.dump(record, out, indent=2, sort_keys=False)
        LOG.info("wrote %s (%s)", path, record["status"])


def setup_logging(verbose: bool = False) -> None:
    """Configure the root logger for the test scripts."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
