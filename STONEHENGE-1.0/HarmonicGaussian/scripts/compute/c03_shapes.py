"""Peak and fold for three barrier shapes with the same area and top curvature (Gaussian,
Lorentzian, dipole-type), and for the Gaussian offset by delta = 0.1 (Sect. 5.4).

Output: data/shapes.pkl. Used by the check branch_numbers.
"""

import numpy as np
from harmonic_gaussian.branches import solve_case, SHAPES
from harmonic_gaussian.storage import main, pick, log


def compute():
    fs = pick(np.round(np.geomspace(0.05, 1.8, 14), 4), np.array([0.15, 0.6]))
    out = dict(fs=fs)
    for name in ["gauss", "lorentz", "dipole"]:
        out[name] = [solve_case(shape=name, s=SHAPES[name], f=f, Wmax=2.5, down=False) for f in fs]
        log(name, "done")
    out["gauss_d0.1"] = [solve_case(f=f, delta=0.1, Wmax=2.5, down=True) for f in fs]
    log("offset done")
    for k, v in out.items():
        if k != "fs":
            for o in v:
                o.pop("branch", None)
    return out


def report(SH):
    print("Peak / fold frequency for barriers of equal area (xi = 0.01)")
    for i, f in enumerate(SH["fs"]):
        print(
            f"  f={f:<7}",
            "  ".join(
                f"{k}: {SH[k][i].get('W_peak', np.nan):.5f}/{SH[k][i].get('W_fold', np.nan):.4f}"
                for k in ["gauss", "lorentz", "dipole", "gauss_d0.1"]
            ),
        )


if __name__ == "__main__":
    main("shapes", compute, report)
