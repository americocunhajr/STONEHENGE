"""Similarity parameter eps = 2 beta_eff M0 c^2/(pi f^3) (Sect. 5.3): peak gap, fold distance,
half-power bandwidth and Floquet modulus at the peak, from data/branches_eps.pkl.

Output: data/similarity.pkl (Fig. 6). Prints the ratios quoted in Sect. 5.3.
"""

import numpy as np
from harmonic_gaussian.model import beff, ctot, Wlin
from harmonic_gaussian.theory import epsilon
from harmonic_gaussian.branches import DUFF
from harmonic_gaussian.storage import main, load


def halfpower(br):
    """Half-power band of the stable branch around its maximum: (Omega_lo, Omega_hi, Omega_max, Pe_max)."""
    W, P, mu = br[:, 0], br[:, 2], br[:, 3]
    st = mu < 1 + 1e-6
    ip = int(np.argmax(np.where(st, P, -1)))
    Pm = P[ip]
    lo = ip
    while lo > 0 and st[lo - 1] and P[lo - 1] >= Pm / 2:
        lo -= 1
    hi = ip
    while hi < len(W) - 1 and st[hi + 1] and P[hi + 1] >= Pm / 2:
        hi += 1
    return W[lo], W[hi], W[ip], Pm


def compute():
    R = load("branches_eps")["res"]
    rows = []
    duff = []
    for key, lst in R.items():
        if key[0] != "gauss":
            continue
        _, a, g, xi, lam = key
        M0 = g * np.sqrt(np.pi / a)
        be = beff(lam)
        c = ctot(xi, lam)
        Wl = Wlin(xi, lam)
        for o in lst:
            if not o.get("ok") or o["br"] is None:
                continue
            f = o["f"]
            eps = epsilon(f, c, be, M0)
            lo, hi, Wp, Pm = halfpower(o["br"])
            rows.append(
                dict(
                    key=key,
                    a=a,
                    g=g,
                    xi=xi,
                    lam=lam,
                    f=f,
                    c=c,
                    fc=f / c,
                    eps=eps,
                    M0=M0,
                    gap=(Wl - o["W_peak"]) / (c / 2),
                    fold=(o["W_fold"] - np.sqrt(be)) / (c / 2) if np.isfinite(o.get("W_fold", np.nan)) else np.nan,
                    A_peak=o["A_peak"],
                    bw=(hi - lo) / c,
                    bw_lo=lo,
                    bw_hi=hi,
                    W_fold=o.get("W_fold", np.nan),
                    W_peak=o["W_peak"],
                    mu=o["peak_maxmu"],
                    mu_lin=np.exp(-np.pi * ctot(xi, lam, o["W_peak"]) / o["W_peak"]),
                )
            )
    for key, lst in R.items():
        if key[0] != "duff":
            continue
        for o in lst:
            if not o.get("ok") or o["br"] is None:
                continue
            lo, hi, Wp, Pm = halfpower(o["br"])
            duff.append(
                dict(name=key[1], f=o["f"], fc=o["f"] / ctot(0.01, 0.05), bw=(hi - lo) / ctot(0.01, 0.05), Wp=Wp, Pm=Pm)
            )
    return dict(rows=rows, duff=duff)


def report(E):
    rows, duff = E["rows"], E["duff"]
    r = [q for q in rows if q["eps"] < 0.1]
    g = np.array([q["gap"] / q["eps"] for q in r])
    g2 = np.array([q["gap"] / (q["eps"] * (1 + 3 / (4 * q["a"] * q["A_peak"] ** 2))) for q in r])
    print(
        f"peak gap / eps for eps<0.1: {g.min():.3f}-{g.max():.3f}  (n={len(g)}); with M2 correction {g2.min():.3f}-{g2.max():.3f}"
    )
    for lo_, hi_ in [(0.1, 1), (1, 10)]:
        r = [q for q in rows if lo_ <= q["eps"] < hi_]
        if r:
            print(
                f'  eps in [{lo_},{hi_}): gap/eps {min(q["gap"]/q["eps"] for q in r):.3f}-{max(q["gap"]/q["eps"] for q in r):.3f}'
            )
    r = [q for q in rows if np.isfinite(q["fold"]) and q["eps"] < 0.02]
    fr = np.array([q["fold"] / (2 / (3 * np.sqrt(3 * q["eps"]))) for q in r])
    print(f"fold/(2/(3 sqrt(3 eps))) for eps<0.02: {fr.min():.3f}-{fr.max():.3f} (n={len(fr)})")
    r = [q for q in rows if q["eps"] < 0.3]
    bw = np.array([q["bw"] for q in r])
    print(f"Gaussian bandwidth/c for eps<0.3: {bw.min():.3f}-{bw.max():.3f} (n={len(bw)})")
    for q in sorted([q for q in rows if q["key"] == ("gauss", 1, 1, 0.01, 0.05)], key=lambda q: q["f"]):
        print(
            f"   base f={q['f']:.4f} eps={q['eps']:.3g} bw/c={q['bw']:.3f} gap/eps={q['gap']/q['eps']:.3f} mu={q['mu']:.4f} lin {q['mu_lin']:.4f}"
        )
    for name in DUFF:
        d = [q for q in duff if q["name"] == name]
        print(name, "Duffing bandwidth/c:", " ".join(f"{q['fc']:.0f}:{q['bw']:.1f}" for q in d))
    mm = np.array([abs(q["mu"] - q["mu_lin"]) for q in rows if q["eps"] < 0.1])
    print("max |mu_peak - exp(-pi c/W)| for eps<0.1:", mm.max())
    for key in sorted({q["key"] for q in rows if q["key"][1:3] != (1, 1)}):
        lst = [q for q in rows if q["key"] == key]
        print(
            key,
            "M0=%.3f" % lst[0]["M0"],
            "n=",
            len(lst),
            "eps range %.2g-%.2g" % (min(q["eps"] for q in lst), max(q["eps"] for q in lst)),
        )


if __name__ == "__main__":
    main("similarity", compute, report)
