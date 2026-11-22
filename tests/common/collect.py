#!/usr/bin/env python3
"""Collect ``metrics.json`` files of the test suite into one summary.

Usage::

    collect.py RUN_ROOT [--cases c1 c2 ...] [--json summary.json]

Prints one Markdown table per case (variant x metric, value and status)
and writes a combined JSON file.
"""

import argparse
import glob
import json
import logging
import os
from typing import Any, Dict, List

LOG = logging.getLogger("collect")


def fmt(value: Any) -> str:
    """Format a metric value compactly for a table cell."""
    if isinstance(value, float):
        return "{:.4g}".format(value)
    if isinstance(value, list):
        return "[" + ", ".join(fmt(v) for v in value[:4]) + \
            (", ..." if len(value) > 4 else "") + "]"
    return str(value)


def table(records: List[Dict[str, Any]]) -> str:
    """Return a Markdown table for the records of one case."""
    names = []  # type: List[str]
    for rec in records:
        for name in rec["metrics"]:
            if name not in names:
                names.append(name)
    head = "| variant | status | " + " | ".join(names) + " |"
    sep = "|" + "---|" * (len(names) + 2)
    rows = [head, sep]
    for rec in records:
        cells = []
        for name in names:
            item = rec["metrics"].get(name)
            if item is None:
                cells.append("")
                continue
            mark = {"PASS": "", "FAIL": " ✗", "REPORT": ""}[item["status"]]
            cells.append(fmt(item["value"]) + mark)
        rows.append("| {} | {} | {} |".format(rec["variant"], rec["status"],
                                              " | ".join(cells)))
    return "\n".join(rows)


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_root")
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--json", default=None)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    cases = args.cases or sorted(
        d for d in os.listdir(args.run_root)
        if os.path.isdir(os.path.join(args.run_root, d)))
    files = []  # type: List[str]
    for case in cases:
        files += sorted(glob.glob(os.path.join(args.run_root, case, "**",
                                               "metrics.json"),
                                  recursive=True))
    by_case = {}  # type: Dict[str, List[Dict[str, Any]]]
    for path in files:
        with open(path) as src:
            rec = json.load(src)
        rec["path"] = path
        rec.get("extra", {}).pop("tau_columns", None)
        if args.cases and rec["case"] not in args.cases:
            continue
        by_case.setdefault(rec["case"], []).append(rec)
    for case in sorted(by_case):
        print("\n### {}\n".format(case))
        print(table(by_case[case]))
        for rec in by_case[case]:
            if rec.get("run_error"):
                print("\n* {}: ERROR {}".format(rec["variant"],
                                                rec["run_error"][:300]))
    if args.json:
        with open(args.json, "w") as out:
            json.dump(by_case, out, indent=2)
        LOG.info("wrote %s", args.json)


if __name__ == "__main__":
    main()
