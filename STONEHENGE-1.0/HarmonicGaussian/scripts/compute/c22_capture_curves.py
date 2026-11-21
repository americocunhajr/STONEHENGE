"""Fraction captured from rest at Omega = 1 against forcing, in steps of 0.0025 from 0.15, for
three damping levels; 256 starts, horizon 4000 periods (Figs. 8e and 10a).

Output: data/capture_curves.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot, A_WELL
from harmonic_gaussian.ensembles import run_capture, rest_starts
from harmonic_gaussian.storage import main, pick, log


def compute():
    x0, y0, phi = rest_starts(pick(128, 8), -A_WELL)
    out = {}
    for xi in [0.005, 0.01, 0.02]:
        c = ctot(xi, 0.05)
        rows = []
        for f in pick(np.round(np.arange(0.15, 0.3001, 0.0025), 4), np.round(np.arange(0.18, 0.3001, 0.02), 4)):
            T = run_capture(
                x0, y0, phi, params("gauss", xi=xi, f=f), 1.0, pick(4000, 1000), 256, max(0.5 * f / c, 2.0)
            )[:, 0]
            rows.append((f, float(np.mean(T >= 0)), float(np.mean((T >= 0) & (T < 1000)))))
            log(xi, rows[-1])
            if len(rows) > 3 and all(r[1] == 1.0 for r in rows[-3:]):
                break
        out[xi] = dict(c=c, rows=rows)
    return out


def report(H):
    for xi, d in H.items():
        print(f'xi={xi} (c={d["c"]:.5f}): (f, captured within 4000, within 1000)')
        print("  " + "  ".join(f"{f:.4f}:{a:.3f}/{b:.3f}" for f, a, b in d["rows"]))


if __name__ == "__main__":
    main("capture_curves", compute, report)
