"""Low end of the stable interwell band (pitchfork in Fig. 8a): bifurcation type and orbit
symmetry, followed downwards with small continuation steps; Gaussian barrier, symmetric and
offset by delta = 0.1.

Output: data/lowend.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot
from harmonic_gaussian.orbits import march, continuation, analyse, multiplier_type
from harmonic_gaussian.storage import main, pick, log


def compute():
    res = []
    for delta, fs in [(0.0, np.round(np.geomspace(0.05, 1.2, 16), 4)), (0.1, np.round(np.geomspace(0.05, 1.8, 14), 4))]:
        for f in pick(fs, fs[[6, 12]]):
            p = params("gauss", xi=0.01, lam=0.05, f=f, W=1.0, delta=delta)
            A0 = f / ctot(0.01, 0.05)
            z = march(np.array([delta, A0, 0.0]), p, 512, 400)
            dn = continuation(z, p, N=512, direction=-1, Wmax=2.0, Wmin=0.25, ds=0.02, ds_max=0.05, maxpts=6000)
            W, A, Pe, mu, un, folds = analyse(dn)
            k = np.where(un)[0]
            if not len(k):
                res.append(dict(delta=delta, f=f, type="none"))
                continue
            i = k[0]
            res.append(
                dict(
                    delta=delta,
                    f=f,
                    W=0.5 * (W[i] + W[i - 1]),
                    type=multiplier_type(dn[i]["mu"]),
                    asym_before=dn[i - 1]["xmax"] + dn[i - 1]["xmin"],
                    asym_after=dn[i]["xmax"] + dn[i]["xmin"],
                    fold_before=bool(folds and folds[0] < i),
                )
            )
            log(res[-1])
    return res


def report(res):
    print("Low end of the stable interwell band: (delta, f, Omega, multiplier, x_max + x_min before/after)")
    for q in res:
        if q["type"] == "none":
            print(f"  delta={q['delta']} f={q['f']}: stable down to Omega = 0.25")
            continue
        print(
            f"  delta={q['delta']} f={q['f']:<7} Omega={q['W']:.4f} {q['type']:3s} asymmetry {q['asym_before']:+.2e} -> {q['asym_after']:+.2e}"
        )


if __name__ == "__main__":
    main("lowend", compute, report)
