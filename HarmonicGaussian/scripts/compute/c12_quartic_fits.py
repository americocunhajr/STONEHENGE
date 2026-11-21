"""Quartic models fitted by least squares to the Gaussian restoring force over |x| <= X, and the
frequency of maximum power they predict, including fits whose range equals the stroke they
predict (Sect. 7, crosses in Fig. 7a).

Output: data/quartic_fits.pkl.
"""

import numpy as np
from harmonic_gaussian.model import BASE, ctot
from harmonic_gaussian.storage import main

CHI, KAP = BASE["chi"], BASE["kap"]
F = lambda x: x - 2 * x * np.exp(-x * x)


def fit(X, n=4001):
    """(k1, k3) with F ~ -k1 x + k3 x^3 on [-X, X]."""
    x = np.linspace(-X, X, n)
    k, *_ = np.linalg.lstsq(np.c_[-x, x**3], F(x), rcond=None)
    return k


def jump(k1, k3, f, xi=0.01, lam=0.05):
    """Maximum-power (jump) frequency of the quartic from single-harmonic balance at quadrature."""
    W = 1.0
    for _ in range(200):
        ke = CHI * KAP * W * W / (lam**2 + W * W)
        c = ctot(xi, lam, W)
        W = np.sqrt(0.5 * ((ke - k1) + np.sqrt((ke - k1) ** 2 + 3 * k3 * f * f / c**2)))
    return W


def compute():
    res = []
    for X in [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 14.0, 28.0]:
        k1, k3 = fit(X)
        res.append(dict(X=X, k1=k1, k3=k3, J=[jump(k1, k3, f) if k3 > 0 else np.nan for f in (0.15, 0.3, 0.6)]))
    sc = []
    for f in [0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.9, 1.2]:
        X = f / ctot(0.01, 0.05)
        for _ in range(100):
            k1, k3 = fit(X)
            W = jump(k1, k3, f)
            Xn = f / (ctot(0.01, 0.05, W) * W)
            if abs(Xn - X) < 1e-8:
                break
            X = 0.5 * (X + Xn)
        sc.append(dict(f=f, X=X, k1=k1, k3=k3, W=W))
    return dict(fixed=res, sc=sc)


def report(Fi):
    print("  X     k1       k3      wells?   a      | jump at f=0.15, 0.3, 0.6")
    for q in Fi["fixed"]:
        k1, k3 = q["k1"], q["k3"]
        wells = k1 > 0
        print(
            f"{q['X']:5.1f} {k1:+.4f} {k3:+.5f}  {str(wells):5}  {np.sqrt(k1/k3) if wells and k3 > 0 else np.nan:.3f} | "
            + " ".join(f"{j:.4f}" for j in q["J"])
        )
    print("self-consistent fits (X = stroke at the quartic jump):")
    for q in Fi["sc"]:
        print(f"  f={q['f']}: X={q['X']:.2f} k1={q['k1']:+.4f} k3={q['k3']:.2e} jump {q['W']:.4f}")


if __name__ == "__main__":
    main("quartic_fits", compute, report)
