"""Locations of the stored data, figures and printed results, and the quick mode.

Every script reads and writes through these helpers, so it can be run from any directory.

Quick mode (environment variable ``HG_QUICK=1``, or ``python run_all.py compute --quick``) runs the
compute scripts on small grids and writes to ``data_quick/`` and ``results_quick/``; inputs are read
from ``data_quick/`` when present and from ``data/`` otherwise. It checks that a script runs end to
end in minutes; its numbers are not those of the paper.
"""

import os
import sys
import pickle
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
QUICK = os.environ.get("HG_QUICK", "0") == "1"
DATA = ROOT / ("data_quick" if QUICK else "data")
STORED = ROOT / "data"
FIGURES = ROOT / "figures"
RESULTS = ROOT / ("results_quick" if QUICK else "results")


def pick(full, quick):
    """The production value, or the reduced one in quick mode."""
    return quick if QUICK else full


def load(name):
    """Load data/<name>.pkl (in quick mode, the quick version if it exists)."""
    path = DATA / f"{name}.pkl"
    if not path.exists():
        path = STORED / f"{name}.pkl"
    with open(path, "rb") as fh:
        return pickle.load(fh)


def exists(name):
    return (DATA / f"{name}.pkl").exists()


def save(obj, name):
    """Write data/<name>.pkl (data_quick/ in quick mode)."""
    DATA.mkdir(exist_ok=True)
    with open(DATA / f"{name}.pkl", "wb") as fh:
        pickle.dump(obj, fh)


def log(*args):
    """Progress messages go to stderr, so that stdout holds only the results of a script."""
    print(*args, file=sys.stderr, flush=True)


def main(name, compute, report):
    """Entry point of a compute script: compute, save data/<name>.pkl and print the report;
    with ``--report``, print the report from the stored data only."""
    if "--report" in sys.argv[1:]:
        report(load(name))
    else:
        out = compute()
        save(out, name)
        report(out)
