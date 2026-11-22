"""Parse ibScourFoam / pimpleFoam logs into per-step records.

Extracted per time step (when present in the log):

* ``time`` [s], Courant numbers (stock and "ImmersedBoundary" variants),
* continuity errors (stock and the masked "IB" variant),
* p-solver iterations (sum over correctors) and final residuals,
* sediment information printed by ``sediment_dual`` (dIta, ita, Shields
  number, wall shear stress on the dual mesh) and the sand-slide mass
  change,
* timings printed with ``DebugSwitches { immersedBoundaryFvMesh 1; }``.

Also detected: FOAM FATAL errors, NaN/Inf in residuals, sand-slide
non-convergence warnings, and the last time reached.
"""

import math
import re
from typing import Any, Dict, List

# Markers of a crash (not the "FPE trapping enabled" banner line)
FATAL_TAGS = ("FOAM FATAL", "sigFpe::sigHandler", "sigSegv::sigHandler",
              "Floating point exception (core", "Segmentation fault",
              "mpirun noticed", "MPI_ABORT")

RE_TIME = re.compile(r"^Time = (\S+)")
RE_CO = re.compile(
    r"^(ImmersedBoundary )?Courant Number mean: (\S+) max: (\S+)")
RE_CONT = re.compile(
    r"^(IB )?time step continuity errors : sum local = (\S+), "
    r"global = (\S+), cumulative = (\S+)"
)
RE_SOLVE = re.compile(
    r"Solving for (\w+), Initial residual = (\S+), "
    r"Final residual = (\S+), No Iterations (\d+)"
)
RE_SED = re.compile(
    r"^\s*(dIta\(m\)|ita\(m\)|shieldsNumber|shieldsNumberC|wallShearStressMag)"
    r"\s+(\S+)\s+(\S+)\s+(\S+)"
)
RE_SLIDE = re.compile(r"Mass change after sandSlide: (\S+)%, (\S+)")
RE_TIMING = re.compile(
    r"^(Total_Time|IB_Update_Time|IB_CFD_Time|Post_Time) = (\S+) s"
)
RE_EXEC = re.compile(r"^ExecutionTime = (\S+) s\s+ClockTime = (\S+) s")
RE_DZ = re.compile(r"^dzPoint (\S+) (\S+) (\S+)")
RE_BOUND = re.compile(r"^bounding (\w+), min: (\S+) max: (\S+)")
RE_SCALE = re.compile(r"^Scaling down (\w+) flux with factor = (\S+)")
RE_LEAK = re.compile(r"^IB blocked-face leakage flux : sum \|phi\| = (\S+)")


def _num(text: str) -> float:
    """Convert a log token to float; returns NaN on failure."""
    try:
        return float(text)
    except ValueError:
        return float("nan")


def parse_log(path: str) -> Dict[str, Any]:
    """Parse a solver log.

    Parameters
    ----------
    path : str
        Log file (``log.ibScourFoam`` or ``log.pimpleFoam``).

    Returns
    -------
    dict
        ``steps`` (list of per-step dicts), ``fatal`` (str or ``None``),
        ``finished`` (bool, "End" reached), ``n_slide_warnings`` (int),
        ``startup_clock_s`` (float, ClockTime at the first step).
    """
    steps: List[Dict[str, Any]] = []
    cur = None  # type: Any
    fatal = None
    finished = False
    n_slide_warn = 0
    in_fatal = False
    fatal_lines: List[str] = []
    with open(path, errors="replace") as log:
        for raw in log:
            line = raw.rstrip("\n")
            if in_fatal:
                if line.strip():
                    fatal_lines.append(line.strip())
                if len(fatal_lines) > 12:
                    in_fatal = False
                continue
            match = RE_TIME.match(line)
            if match:
                cur = {"time": _num(match.group(1)), "p_iters": 0,
                       "solves": {}, "bounding": []}
                steps.append(cur)
                continue
            if any(tag in line for tag in FATAL_TAGS):
                in_fatal = True
                fatal_lines = [line.strip()]
                continue
            if "Sand slide not converged" in line:
                n_slide_warn += 1
            if line.strip() == "End":
                finished = True
            if cur is None:
                continue
            _parse_step_line(line, cur)
    if fatal_lines:
        fatal = " | ".join(fatal_lines)
    startup = None
    for st in steps:
        if "clock_s" in st:
            startup = st["clock_s"]
            break
    return {
        "steps": steps,
        "fatal": fatal,
        "finished": finished,
        "n_slide_warnings": n_slide_warn,
        "first_step_clock_s": startup,
    }


def _parse_step_line(line: str, cur: Dict[str, Any]) -> None:
    """Update the current step record ``cur`` from one log line."""
    match = RE_CO.match(line)
    if match:
        key = "co_ib" if match.group(1) else "co"
        cur[key + "_mean"] = _num(match.group(2))
        cur[key + "_max"] = _num(match.group(3))
        return
    match = RE_CONT.match(line)
    if match:
        key = "cont_ib" if match.group(1) else "cont"
        cur[key + "_local"] = _num(match.group(2))
        cur[key + "_global"] = _num(match.group(3))
        cur[key + "_cum"] = _num(match.group(4))
        return
    match = RE_SOLVE.search(line)
    if match:
        field = match.group(1)
        init = _num(match.group(2))
        cur["solves"].setdefault(field, []).append(init)
        if field == "p":
            cur["p_iters"] += int(match.group(4))
        if not math.isfinite(init):
            cur["nan_residual"] = field
        return
    match = RE_SED.match(line)
    if match:
        key = match.group(1).replace("(m)", "")
        cur["sed_" + key] = [_num(match.group(i)) for i in (2, 3, 4)]
        return
    match = RE_SLIDE.search(line)
    if match:
        cur["slide_mass_change_pct"] = _num(match.group(1))
        cur["slide_mass_change"] = _num(match.group(2))
        return
    match = RE_TIMING.match(line)
    if match:
        cur[match.group(1)] = _num(match.group(2))
        return
    match = RE_DZ.match(line)
    if match:
        cur["dz_point"] = [_num(match.group(i)) for i in (1, 2, 3)]
        return
    match = RE_BOUND.match(line)
    if match:
        cur["bounding"].append(match.group(1))
        return
    match = RE_SCALE.match(line)
    if match:
        fac = _num(match.group(2))
        cur["adjust_phi_min_factor"] = min(
            fac, cur.get("adjust_phi_min_factor", fac))
        return
    match = RE_LEAK.match(line)
    if match:
        cur["leak_sum_abs_phi"] = _num(match.group(1))
        return
    match = RE_EXEC.match(line)
    if match:
        cur["cpu_s"] = _num(match.group(1))
        cur["clock_s"] = _num(match.group(2))


def summarise(log: Dict[str, Any]) -> Dict[str, Any]:
    """Return scalar summary metrics of a parsed log.

    Keys: ``n_steps``, ``t_end`` [s], ``max_co``, ``max_co_ib``,
    ``max_cont_local``, ``max_cont_ib_local``, ``mean_p_iters``,
    ``wall_s`` [s] (ClockTime at the last step), ``s_per_step`` [s]
    (mean clock time per step after the first), ``fatal``, ``finished``,
    ``first_nan_time`` [s], ``min_adjust_phi_factor`` (smallest
    global IB flux-rescaling factor of immersedBoundaryAdjustPhi) and
    ``max_leak_sum_abs_phi`` [m^3/s] (when the solver prints them).
    """
    steps = log["steps"]

    def _max(key: str) -> float:
        vals = [s[key] for s in steps if key in s]
        return max(vals) if vals else float("nan")

    def _min(key: str) -> float:
        vals = [s[key] for s in steps if key in s]
        return min(vals) if vals else float("nan")

    clocks = [s["clock_s"] for s in steps if "clock_s" in s]
    s_per_step = float("nan")
    if len(clocks) > 2:
        s_per_step = (clocks[-1] - clocks[0]) / (len(clocks) - 1)
    first_nan = None
    for st in steps:
        if "nan_residual" in st:
            first_nan = st["time"]
            break
    piters = [s["p_iters"] for s in steps if s["p_iters"]]
    return {
        "n_steps": len([s for s in steps if "clock_s" in s]),
        "t_end": steps[-1]["time"] if steps else float("nan"),
        "max_co": _max("co_max"),
        "max_co_ib": _max("co_ib_max"),
        "max_cont_local": _max("cont_local"),
        "max_cont_ib_local": _max("cont_ib_local"),
        "mean_p_iters": sum(piters) / len(piters) if piters else float("nan"),
        "wall_s": clocks[-1] if clocks else float("nan"),
        "s_per_step": s_per_step,
        "fatal": log["fatal"],
        "finished": log["finished"],
        "first_nan_time": first_nan,
        "min_adjust_phi_factor": _min("adjust_phi_min_factor"),
        "max_leak_sum_abs_phi": _max("leak_sum_abs_phi"),
    }
