"""Peak of power and fold of the interwell branch against forcing.

Gaussian barrier at four damping levels and the quartic models (nominal and three calibrations),
16 forcing levels from 0.05 to 1.2. Peaks of very large orbits (A_peak > 60) are re-maximized with
2048 steps per period. Used by Figs. 7, 8 and 10 and by scripts c09 and checks.
Output: data/scaling.pkl.
"""

import numpy as np
from harmonic_gaussian.branches import solve_case, solve_duffing, refine_peak_direct, DUFF
from harmonic_gaussian.storage import main, pick, log

XIS = [0.005, 0.01, 0.02, 0.04]


def compute():
    fs = pick(np.round(np.geomspace(0.05, 1.2, 16), 4), np.array([0.15, 0.6]))
    out = dict(fs=fs, gauss={}, duff={}, calib={})
    for xi in XIS:
        out["gauss"][xi] = [solve_case(f=f, xi=xi) for f in fs]
        out["duff"][xi] = [solve_duffing(*DUFF["nominal"], f, xi=xi) for f in fs]
        log("xi", xi, "done")
    for name in ["(a,h)", "(a,k_w)", "Taylor"]:
        out["calib"][name] = [solve_duffing(*DUFF[name], f, xi=0.01) for f in fs]
        log(name, "done")
    for v in list(out["gauss"].values()) + list(out["duff"].values()) + list(out["calib"].values()):
        for o in v:
            o.pop("branch", None)
    out["patched"] = []
    for xi in XIS:
        for o in out["gauss"][xi]:
            if o.get("ok") and o["A_peak"] > 60:
                (W1, P1), (W2, P2) = refine_peak_direct(o, xi=xi)
                out["patched"].append((xi, o["f"], o["W_peak"], W1, W2))
                o["W_peak"], o["Pe_peak"] = W2, P2
                o["peak_note"] = "N=2048 direct"
    return out


def report(S):
    print("Peak and fold of the Gaussian interwell branch (lambda = 0.05)")
    for xi in XIS:
        print(f"xi={xi}")
        for o in S["gauss"][xi]:
            if o.get("ok"):
                print(
                    f"  f={o['f']:<7} W_peak={o['W_peak']:.6f} Pe_peak={o['Pe_peak']:.4e} A_peak={o['A_peak']:8.3f} W_fold={o['W_fold']:.5f} W_low={o.get('W_low', np.nan):.4f} ({o.get('low_type', '')})"
                )
            else:
                print(f"  f={o['f']:<7} no interwell orbit")
    for xi, f, W0, W1, W2 in S.get("patched", []):
        print(f"peak re-maximized: xi={xi} f={f}: continuation {W0:.8f}, N=1024 {W1:.8f}, N=2048 {W2:.8f}")
    print("Jump (fold) frequency of the quartic models, xi = 0.01")
    for name, lst in [("nominal", S["duff"][0.01])] + list(S["calib"].items()):
        print(f"  {name:8s}", " ".join(f"{o['W_fold']:.3f}" if o.get("ok") else "nan" for o in lst))


if __name__ == "__main__":
    main("scaling", compute, report)
