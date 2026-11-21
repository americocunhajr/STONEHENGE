"""Capture from rest at Omega = 1 for the point-dipole magnet pair (Sect. 8.3): fraction captured
within 8000 periods, within 1000, median capture time, for three damping levels.

Output: data/magnet_capture.pkl.
"""

import numpy as np
from scipy.optimize import brentq
from harmonic_gaussian.model import params, dU, ctot, MAGNET_S
from harmonic_gaussian.ensembles import run_capture, rest_starts, interwell
from harmonic_gaussian.storage import main, pick, log


def compute():
    a = brentq(lambda x: dU(x, params("magnet", s=MAGNET_S)), 0.2, 1.5)
    x0, y0, phi = rest_starts(pick(64, 8), -a)
    out = {}
    for xi in pick([0.01, 0.02, 0.03], [0.01]):
        cc = ctot(xi, 0.05)
        rows = []
        for f in pick([0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32], [0.22, 0.30]):
            p = params("magnet", xi=xi, lam=0.05, f=f, s=MAGNET_S)
            o = run_capture(x0, y0, phi, p, 1.0, pick(8000, 1000), 256, max(0.5 * f / cc, 2.0))
            T = o[:, 0]  # the resonant orbit must exceed twice the gap
            unc = T < 0
            rows.append(
                dict(
                    f=f,
                    cap=float(np.mean(T >= 0)),
                    cap1000=float(np.mean((T >= 0) & (T < 1000))),
                    medT=float(np.median(T[T >= 0])) if (T >= 0).any() else -1,
                    unc_inter=float(np.mean(unc & interwell(o[:, 2], o[:, 3]))),
                    unc_per=float(np.mean(unc & (o[:, 1] < 1e-3))),
                )
            )
            log("xi", xi, rows[-1])
        out[xi] = rows
    return dict(well=a, cap=out)


def report(D):
    print(f"magnet pair, well at {D['well']:.4f}; capture from rest at Omega = 1")
    for xi, rows in D["cap"].items():
        for r in rows:
            print(
                f"  xi={xi} f={r['f']:.2f}: captured {r['cap']:.3f} (within 1000: {r['cap1000']:.3f}), median T {r['medT']:.1f}, "
                f"uncaptured interwell {r['unc_inter']:.3f}, uncaptured periodic {r['unc_per']:.3f}"
            )


if __name__ == "__main__":
    main("magnet_capture", compute, report)
