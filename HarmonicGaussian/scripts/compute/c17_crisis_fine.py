"""Survival of the chaotic interwell motion just above the capture threshold (Omega = 1, baseline):
first-capture times of 256 starts up to 40000 periods at ten forcing levels from 0.216 to 0.23
(Fig. 9c, d); and a stroboscopic sample of the coexisting attractor at f = 0.24.

Output: data/crisis_fine.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot, A_WELL
from harmonic_gaussian.ensembles import run_capture, strobe, rest_starts
from harmonic_gaussian.storage import main, pick, log


def compute():
    c = ctot(0.01, 0.05)
    x0, y0, phi = rest_starts(pick(128, 8), -A_WELL)
    res = []
    for f in pick([0.216, 0.2165, 0.217, 0.2175, 0.218, 0.219, 0.22, 0.222, 0.225, 0.23], [0.22, 0.23]):
        o = run_capture(x0, y0, phi, params("gauss", f=f), 1.0, pick(40000, 4000), 256, 0.5 * f / c)
        res.append(dict(f=f, T=o[:, 0], std=o[:, 1]))
        log("f", f, "done")
    S = strobe(
        -A_WELL, 0.0, 0.0, params("gauss", f=0.24), 1.0, pick(20000, 2000), 60, 256, 7 * 2 * np.pi / 64
    )  # start 7 of the 128-start set
    return dict(res=res, other=S)


def report(D):
    for r in D["res"]:
        T = r["T"]
        cap = T >= 0
        print(
            f"f={r['f']}: captured {cap.mean():.3f}, mean T {np.mean(T[cap]) if cap.any() else -1:.0f}, median {np.median(T[cap]) if cap.any() else -1:.0f}"
        )
    print("f=0.24 other attractor, stroboscopic x:", np.round(D["other"][:12, 0], 4))


if __name__ == "__main__":
    main("crisis_fine", compute, report)
