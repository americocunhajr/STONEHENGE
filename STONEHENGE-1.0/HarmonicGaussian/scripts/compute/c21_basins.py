"""Basins of the resonant orbit at Omega = 1 (phase 0, v = 0) on a 401 x 321 grid of initial states,
below (f = 0.21) and above (f = 0.225) the crisis, and a Poincare section of 2 x 10^5 points of the
chaotic interwell attractor at f = 0.21 (Fig. 9a, b).

Output: data/basins.pkl (about 5 MB). Runtime: a few minutes on two cores.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot, A_WELL
from harmonic_gaussian.ensembles import run_capture, strobe
from harmonic_gaussian.storage import main, pick, log


def compute():
    c = ctot(0.01, 0.05)
    xs = np.linspace(-3, 3, pick(401, 41))
    ys = np.linspace(-3, 3, pick(321, 33))
    X, Y = np.meshgrid(xs, ys)
    out = dict(xs=xs, ys=ys)
    for f, nmax in [(0.21, 600), (0.225, 4000)]:
        o = run_capture(
            X.ravel(), Y.ravel(), np.zeros(X.size), params("gauss", f=f), 1.0, pick(nmax, 300), 256, 0.5 * f / c
        )
        out[("basin", f)] = o[:, 0].reshape(X.shape)
        out[("nmax", f)] = nmax
        log("basin", f, "done")
    out[("strobe", 0.21)] = strobe(A_WELL, 0.0, 0.0, params("gauss", f=0.21), 1.0, 500, pick(200000, 2000), 256, 0.0)
    return out


def report(out):
    for f in (0.21, 0.225):
        T = out[("basin", f)]
        print(
            f"basin f={f}: {T.shape[1]} x {T.shape[0]} initial states, captured {np.mean(T >= 0):.4f}, latest capture {int(T.max())} periods"
        )
    print(f"Poincare section at f=0.21: {len(out[('strobe', 0.21)])} points")


if __name__ == "__main__":
    main("basins", compute, report)
