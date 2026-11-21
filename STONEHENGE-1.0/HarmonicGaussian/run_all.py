#!/usr/bin/env python3
"""Reproduce the data, figures and printed results of the paper.

    python run_all.py figures              # Figures 1-10 from the stored data (about 1 min)
    python run_all.py results              # all printed numbers from the stored data and analytic checks
    python run_all.py compute [names]      # recompute data (all, or the named scripts, e.g. c22 c23)
    python run_all.py compute --quick      # every compute script on small grids, into data_quick/
    python run_all.py all                  # compute, then results, then figures (several hours)
    python run_all.py list                 # list the scripts
    python run_all.py clean-cache          # remove the Numba caches (needed after editing model.py)

Figures go to figures/, data to data/, printed results to results/<script>.txt. Progress messages
of the compute scripts are shown on the terminal.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
COMPUTE = sorted(p.stem for p in (SCRIPTS / "compute").glob("c[0-9][0-9]_*.py"))
FIGURES = sorted(p.stem for p in (SCRIPTS / "figures").glob("fig[0-9][0-9]_*.py"))
CHECKS = sorted(p.stem for p in (SCRIPTS / "checks").glob("*.py"))


def run(kind, name, args=(), quick=False):
    """Run scripts/<kind>/<name>.py as a module from the root; return its stdout."""
    env = dict(
        os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""), HG_QUICK="1" if quick else "0"
    )
    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-m", f"scripts.{kind}.{name}", *args], cwd=ROOT, env=env, stdout=subprocess.PIPE, text=True
    )
    status = "ok" if r.returncode == 0 else f"FAILED (exit {r.returncode})"
    print(f"  {kind}/{name}: {status}, {time.time() - t0:.0f} s", flush=True)
    if r.returncode != 0:
        sys.exit(1)
    return r.stdout


def write_result(name, text, quick=False):
    out = ROOT / ("results_quick" if quick else "results")
    out.mkdir(exist_ok=True)
    (out / f"{name}.txt").write_text(text)


def figure01():
    """Figure 1 is a TikZ drawing; it needs pdflatex with circuitikz."""
    if shutil.which("pdflatex") is None:
        print("  figures/fig01_device: skipped (pdflatex not found)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(SCRIPTS / "figures" / "fig01_device.tex", tmp)
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "fig01_device.tex"], cwd=tmp, capture_output=True)
        if r.returncode == 0:
            (ROOT / "figures").mkdir(exist_ok=True)
            shutil.copy(Path(tmp) / "fig01_device.pdf", ROOT / "figures" / "fig01_device.pdf")
        print(f"  figures/fig01_device: {'ok' if r.returncode == 0 else 'FAILED'}")


def select(names, pool):
    if not names:
        return pool
    chosen = [p for p in pool if any(p == n or p.startswith(n + "_") for n in names)]
    unknown = [n for n in names if not any(p == n or p.startswith(n + "_") for p in pool)]
    if unknown:
        sys.exit(f"unknown script(s): {unknown}; see python run_all.py list")
    return chosen


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return
    cmd, rest = argv[0], argv[1:]
    quick = "--quick" in rest
    names = [a for a in rest if not a.startswith("--")]
    if cmd == "list":
        print("compute:", *COMPUTE, sep="\n  ")
        print("figures:", "fig01_device (TikZ)", *FIGURES, sep="\n  ")
        print("checks:", *CHECKS, sep="\n  ")
    elif cmd == "figures":
        figure01()
        for n in select(names, FIGURES):
            run("figures", n)
    elif cmd == "compute":
        for n in select(names, COMPUTE):
            write_result(n, run("compute", n, quick=quick), quick)
    elif cmd == "results":
        for n in COMPUTE:
            write_result(n, run("compute", n, ["--report"]))
        for n in CHECKS:
            write_result(n, run("checks", n))
    elif cmd == "all":
        main(["compute"] + rest)
        if not quick:
            main(["results"])
            main(["figures"])
    elif cmd == "clean-cache":
        for d in ROOT.rglob("__pycache__"):
            shutil.rmtree(d)
        print("removed __pycache__ folders")
    else:
        sys.exit(f"unknown command {cmd!r}; see python run_all.py --help")


if __name__ == "__main__":
    main(sys.argv[1:])
