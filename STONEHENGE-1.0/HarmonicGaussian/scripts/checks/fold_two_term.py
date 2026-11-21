"""Fold law with the leading and the two-term large-amplitude stiffness, against continuation
(Sect. 5.2, Supplementary Sect. S1.4)."""

import numpy as np
from harmonic_gaussian.model import beff, ctot
from harmonic_gaussian.theory import M0_GAUSS as M0
from harmonic_gaussian.storage import load

S = load("scaling")
M2 = M0 / 2
be = beff(0.05)
fold1 = lambda f: (2 / 3) * np.sqrt(np.pi / (6 * M0)) * f**1.5


def fold2(f):
    A = np.sqrt(((6 * M0 / np.pi) + np.sqrt((6 * M0 / np.pi) ** 2 + 60 * f * M2 / np.pi)) / (2 * f))
    return -(2 / np.pi) * (M0 / A**3 + 1.5 * M2 / A**5) + f / A


for xi in [0.005, 0.01, 0.02, 0.04]:
    c = ctot(xi, 0.05)
    r1 = []
    r2 = []
    for o in S["gauss"][xi]:
        if not o.get("ok") or not np.isfinite(o["W_fold"]) or o["f"] / c < 13:
            continue
        d = o["W_fold"] ** 2 - be
        r1.append(d / fold1(o["f"]))
        r2.append(d / fold2(o["f"]))
    print(
        f"xi={xi}: f/c>=13  leading-order ratio {min(r1):.3f}-{max(r1):.3f};  two-term ratio {min(r2):.3f}-{max(r2):.3f}"
    )
