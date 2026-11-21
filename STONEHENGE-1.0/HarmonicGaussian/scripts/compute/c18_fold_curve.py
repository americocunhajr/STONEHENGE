"""Fold (saddle-node) curve of the interwell periodic orbits in the (Omega, f) plane, continued in
two parameters from f = 1.4 through the cusp (Sect. 5.4; Fig. 10a). Unknowns (z, Omega, f);
equations: periodicity and det(M - I) = 0.

Output: data/fold_curve.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params
from harmonic_gaussian.orbits import march, continuation, analyse, fold_curve, cusp_of
from harmonic_gaussian.storage import main, pick, log


def compute():
    f0 = 1.40
    p = params("gauss", xi=0.01, lam=0.05, f=f0, W=1.2)
    z = march(np.array([0.0, 1.5, 0.0]), p, 512, 600)
    rec = continuation(z, p, N=512, direction=+1, Wmax=2.5, Wmin=0.3, ds=0.01, ds_max=0.05, maxpts=20000, Astop=0.9)
    W, A, Pe, mu, un, folds = analyse(rec)
    log("turning points at f=1.4:", [W[i] for i in folds])
    i = folds[0]
    C = fold_curve(np.r_[rec[i]["z"], rec[i]["W"], f0], p, N=1024, maxpts=pick(400, 12), f_stop=f0 - 0.05)
    fc, Wc = cusp_of(C)
    return dict(curve=C, fc=fc, Wc=Wc)


def report(D):
    C = D["curve"]
    print(f"fold curve: {len(C)} points, f from {C[:, 4].min():.4f} to {C[:, 4].max():.5f}")
    print("cusp: f_c=%.5f Omega_c=%.5f (max on curve %.5f)" % (D["fc"], D["Wc"], C[:, 4].max()))


if __name__ == "__main__":
    main("fold_curve", compute, report)
