"""First-capture times from rest at Omega = 1 (baseline Gaussian), up to 20000 periods, 128 starts,
twelve forcing levels from 0.12 to 0.30 (Sect. 8.3, Supplementary Sect. S2.2).

Output: data/crisis.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot, A_WELL
from harmonic_gaussian.ensembles import run_capture, rest_starts, interwell
from harmonic_gaussian.storage import main, pick, log


def compute():
    c = ctot(0.01, 0.05)
    x0, y0, phi = rest_starts(pick(64, 8), -A_WELL)
    res = []
    for f in pick([0.12, 0.14, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.24, 0.27, 0.30], [0.21, 0.24]):
        o = run_capture(x0, y0, phi, params("gauss", f=f), 1.0, pick(20000, 2000), 256, 0.5 * f / c)
        res.append(dict(f=f, T=o[:, 0], strobe_std=o[:, 1], inter=interwell(o[:, 2], o[:, 3])))
        log("f", f, "done")
    return res


def report(res):
    for r in res:
        T = r["T"]
        cap = T >= 0
        print(
            f"f={r['f']}: captured {cap.mean():.3f} (by 1000: {np.mean(cap & (T < 1000)):.3f}), median T {np.median(T[cap]) if cap.any() else -1}, "
            f"uncaptured interwell {np.mean(~cap & r['inter']):.3f}, periodic {np.mean(~cap & (r['strobe_std'] < 1e-3)):.3f}"
        )


if __name__ == "__main__":
    main("crisis", compute, report)
