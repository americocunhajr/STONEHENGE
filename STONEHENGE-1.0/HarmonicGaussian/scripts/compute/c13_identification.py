"""Identification of the barrier area M0 = -int x [F(x) - beta x] dx from a sampled static force
(Sect. 3, Fig. 3b): relative RMS error against the window half-width X, for noise levels of 0, 1
and 3% of the peak barrier force and stiffness errors of 0, 0.2 and 1%, with beta known or fitted
on |x| >= 3. Monte Carlo with 4000 replicates per case and a fixed seed.

Output: data/identification.pkl.
"""

import math
import numpy as np
from harmonic_gaussian.storage import main, pick

M0 = np.sqrt(np.pi)
FMAX = np.sqrt(2 / np.e)  # area and peak barrier force for alpha = gamma = 1
DX = 0.05
F = lambda x: x - 2 * x * np.exp(-x * x)


def compute():
    rng = np.random.default_rng(20261120)
    R = pick(4000, 200)

    def est(X, sig, dbeta=0.0, fit=False, X0=3.0):
        x = np.arange(-X, X + 1e-12, DX)
        w = np.full(x.size, DX)
        w[0] = w[-1] = DX / 2
        Fm = F(x)[None, :] + sig * FMAX * rng.standard_normal((R, x.size))
        if fit:
            far = np.abs(x) >= X0
            b = (Fm[:, far] @ x[far]) / (x[far] @ x[far])
        else:
            b = np.full(R, 1.0 + dbeta)
        return -(Fm * x[None, :] - b[:, None] * x[None, :] ** 2) @ w

    stats = lambda M: (np.mean(M) / M0 - 1, np.std(M) / M0, np.sqrt(np.mean((M / M0 - 1) ** 2)))
    out = {}
    Xs = np.array([1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8])
    for sig in [0.0, 0.01, 0.03]:
        for db in [0.0, 0.002, 0.01]:
            out[("known", sig, db)] = np.array([stats(est(X, sig, db)) for X in Xs])
        out[("fit", sig)] = np.array([(np.nan,) * 3 if X < 4 else stats(est(X, sig, fit=True)) for X in Xs])
    out["Xs"] = Xs
    return out


def report(out):
    np.set_printoptions(precision=4, suppress=True)
    for k, v in out.items():
        if k != "Xs":
            print(k, "rel. RMS error:", v[:, 2])
    Xs = out["Xs"]
    print("X:", Xs)
    print(
        "analytic: truncation 2XV(X)+2int_X^inf V ->",
        [f"{(2*X*np.exp(-X*X)+np.sqrt(np.pi)*math.erfc(X))/M0:.2e}" for X in Xs],
    )
    print("stiffness term 2X^3 dbeta/3/M0 for dbeta=0.002:", [f"{2*X**3*0.002/3/M0:.3f}" for X in Xs])


if __name__ == "__main__":
    main("identification", compute, report)
