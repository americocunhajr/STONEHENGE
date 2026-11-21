"""Peak and fold of the Gaussian interwell branch with beam hardening k3 x^4/4 (Sect. 6).

Ten hardening levels from 1e-5 to 1e-2 and ten forcing levels from 0.08 to 0.9 (xi = 0.01).
Used by c09 (capture at the peak with hardening) and the check branch_numbers.
Output: data/hardening.pkl.
"""

import numpy as np
from harmonic_gaussian.branches import solve_case
from harmonic_gaussian.storage import main, pick, log


def compute():
    k3s = pick(np.geomspace(1e-5, 1e-2, 10), np.array([1e-4, 1e-3]))
    fs = pick(np.round(np.geomspace(0.08, 0.9, 10), 4), np.array([0.2, 0.6]))
    res = []
    for k3 in k3s:
        for f in fs:
            o = solve_case(f=f, k3h=k3, Wmax=6.0, down=False)
            o.pop("branch", None)
            o["k3h"] = k3
            res.append(o)
        log("k3", k3, "done")
    return dict(k3s=k3s, fs=fs, res=res)


def report(H):
    print("Peak and fold with hardening k3 (xi = 0.01)")
    for k3 in H["k3s"]:
        rows = [o for o in H["res"] if abs(o["k3h"] - k3) < 1e-12]
        print(
            f"k3={k3:.2e}:",
            " ".join(
                f"f={o['f']}: {o['W_peak']:.4f}/{o.get('W_fold', np.nan):.4f}" if o.get("ok") else f"f={o['f']}: -"
                for o in rows
            ),
        )


if __name__ == "__main__":
    main("hardening", compute, report)
