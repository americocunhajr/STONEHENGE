"""Interwell branches for the similarity analysis (Sect. 5.3, Fig. 6): Gaussian barrier at four
damping levels, at the near-optimal load lambda = 1, with four other barrier sizes (alpha, gamma),
and the four quartic models. Each record keeps the whole branch (Omega, A, Pe, max|mu|) in 'br'.
Peaks of orbits with A_peak sqrt(alpha) > 60 are re-maximized with 2048 steps per period.

Output: data/branches_eps.pkl (about 3 MB).
"""

import numpy as np
from harmonic_gaussian.branches import solve_case, solve_duffing, refine_peak_direct, DUFF
from harmonic_gaussian.storage import main, pick, log


def keep(o):
    o = dict(o)
    br = o.pop("branch", None)
    o["br"] = np.array(br) if br is not None else None
    return o


def compute():
    fs = pick(np.round(np.geomspace(0.05, 1.2, 16), 4), np.array([0.15, 0.6]))
    fs2 = pick(np.round(np.geomspace(0.03, 1.5, 18), 4), np.array([0.15]))
    out = {}
    for xi in [0.005, 0.01, 0.02, 0.04]:
        out[("gauss", 1, 1, xi, 0.05)] = [keep(solve_case(f=f, xi=xi, lam=0.05, Wmax=2.0, down=True)) for f in fs]
        log("xi", xi, "done")
    out[("gauss", 1, 1, 0.01, 1.0)] = [keep(solve_case(f=f, xi=0.01, lam=1.0, Wmax=2.0, down=True)) for f in fs]
    log("lambda=1 done")
    for a, g in [(1.0, 0.75), (1.0, 2.0), (0.5, 2.0), (2.0, 0.5)]:
        out[("gauss", a, g, 0.01, 0.05)] = [
            keep(solve_case(s=(a, 1, g, 0), f=f, xi=0.01, lam=0.05, Wmax=2.5, down=True)) for f in fs2
        ]
        log("alpha, gamma", a, g, "done")
    for name, (k1, k3) in DUFF.items():
        out[("duff", name)] = [keep(solve_duffing(k1, k3, f, Wmax=16, Wmin=0.2, down=True)) for f in fs]
        log("quartic", name, "done")
    patched = []
    for key, lst in out.items():
        if key[0] != "gauss":
            continue
        _, a, g, xi, lam = key
        for o in lst:
            if o.get("ok") and o["A_peak"] * np.sqrt(a) > 60:
                (W1, P1), (W2, P2) = refine_peak_direct(o, xi=xi, lam=lam, s=(a, 1, g, 0))
                patched.append((key, o["f"], o["W_peak"], W1, W2))
                o["W_peak"], o["Pe_peak"] = W2, P2
                o["peak_note"] = "N=2048 direct"
    return dict(fs=fs, fs2=fs2, res=out, patched=patched)


def report(D):
    for key, lst in D["res"].items():
        ok = [o for o in lst if o.get("ok")]
        print(key, f"{len(ok)} of {len(lst)} branches found")
    for key, f, W0, W1, W2 in D.get("patched", []):
        print(f"peak re-maximized: {key} f={f}: continuation {W0:.8f}, N=1024 {W1:.8f}, N=2048 {W2:.8f}")


if __name__ == "__main__":
    main("branches_eps", compute, report)
