"""Point-dipole magnet pair with tip rotation (K = 0.3, eta = 0.046): area, free-oscillation
backbone by high-precision quadrature against the area law, and peak and fold of the forced
interwell branch (Sects. 4 and 5; Figs. 4b and 6).

Output: data/magnet.pkl.
"""

import numpy as np
import mpmath as mp
from scipy.integrate import quad
from harmonic_gaussian.model import MAGNET_S, magnet_barrier, ctot, beff
from harmonic_gaussian.theory import period_mp, epsilon
from harmonic_gaussian.branches import solve_case
from harmonic_gaussian.storage import main, pick, log

K, ETA = MAGNET_S[2], MAGNET_S[3]


def compute():
    M0 = quad(lambda u: magnet_barrier(u, K, ETA), -200, 200, limit=2000)[0]
    out = dict(K=K, eta=ETA, M0=M0)
    mp.mp.dps = 30
    V = lambda u: K * ((2 - u * u) * mp.cos(ETA * u) - 3 * u * mp.sin(ETA * u)) / (1 + u * u) ** mp.mpf(2.5)
    bb = []
    for A in pick([2, 3, 4, 5, 6, 8, 10, 14, 20], [4, 10]):
        w = 2 * mp.pi / period_mp(V, mp.mpf(A))
        bb.append((A, float(1 - w * w), 2 * M0 / (np.pi * A**3)))
    out["backbone"] = bb
    c = ctot(0.01, 0.05)
    be = beff(0.05)
    res = []
    for f in pick([0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.45], [0.15]):
        o = solve_case(shape="magnet", s=MAGNET_S, f=f, xi=0.01, lam=0.05, Wmax=2.0)
        o["eps"] = epsilon(f, c, be, M0)
        res.append(o)
        log("f", f, "done")
    out["forced"] = res
    return out


def report(out):
    print(f"magnet pair K={out['K']} eta={out['eta']}: area M0 = {out['M0']:.6f}")
    print("backbone (A, 1-w^2, area law):", [(a, f"{x:.3e}", f"{y:.3e}") for a, x, y in out["backbone"]])
    for o in out["forced"]:
        print(
            "f",
            o["f"],
            "ok",
            o.get("ok"),
            "W_peak",
            o.get("W_peak"),
            "W_fold",
            o.get("W_fold"),
            "A_peak",
            o.get("A_peak"),
            "eps",
            o["eps"],
        )


if __name__ == "__main__":
    main("magnet", compute, report)
