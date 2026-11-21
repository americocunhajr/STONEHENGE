"""Upper fold and lower turning point of the Gaussian interwell branch near the cusp
(f = 1.2 to 1.70), for the fit of the hysteresis width in the check cusp_fit.

Output: data/cusp.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params
from harmonic_gaussian.orbits import march, continuation, analyse, refine_fold
from harmonic_gaussian.branches import solve_case
from harmonic_gaussian.storage import main, pick, log


def compute():
    res = []
    for f in pick([1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.62, 1.64, 1.66, 1.67, 1.68, 1.69, 1.70], [1.5, 1.6]):
        o = solve_case(f=f, xi=0.01, lam=0.05, Wmax=2.5, down=False)
        if not o.get("ok") or not np.isfinite(o.get("W_fold", np.nan)):
            res.append(dict(f=f, ok=False))
            log(f, "no fold")
            continue
        p = params(
            "gauss", xi=0.01, lam=0.05, f=f, W=o["W_peak"]
        )  # restart from the peak, continue past both turning points
        z = march(np.array([0.0, o["A_peak"] * o["W_peak"], 0.0]), p, 512, 600)
        up = continuation(z, p, N=512, direction=+1, Wmax=2.5, Wmin=0.3, ds=0.01, ds_max=0.05, maxpts=20000, Astop=0.9)
        W, A, Pe, mu, un, folds = analyse(up)
        res.append(
            dict(f=f, ok=True, folds=[refine_fold(up, i)[0] for i in folds[:2]], Afolds=[A[i] for i in folds[:2]])
        )
        log(f, res[-1]["folds"])
    return res


def report(res):
    print("Turning points of the interwell branch near the cusp (upper fold, lower turning point)")
    for q in res:
        print(f"  f={q['f']}:", np.round(q["folds"], 5) if q.get("ok") else "no fold")


if __name__ == "__main__":
    main("cusp", compute, report)
