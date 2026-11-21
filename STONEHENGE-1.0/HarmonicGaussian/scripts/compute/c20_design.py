"""Design example (Sect. 9.4): point-dipole magnet pair, strong coupling chi kap = 0.1, matched
load lambda = 1. Capture from rest at Omega_lin and the resonant orbit (stroke, power, stability)
at six forcing levels, for two damping levels.

Output: data/design.pkl.
"""

import numpy as np
from scipy.optimize import brentq
from harmonic_gaussian.model import params, dU, ctot, beff, Wlin, MAGNET_S
from harmonic_gaussian.orbits import march, flow, newton_shoot, Pe_from
from harmonic_gaussian.ensembles import run_capture, rest_starts
from harmonic_gaussian.storage import main, pick, log

CHI, KAP, LAM = 0.2, 0.5, 1.0


def compute():
    out = {}
    a = brentq(lambda x: dU(x, params("magnet", s=MAGNET_S)), 0.2, 1.5)
    x0, y0, phi = rest_starts(pick(64, 8), -a)
    for xi in [0.0056, 0.01]:
        c = ctot(xi, LAM, chi=CHI, kap=KAP)
        be = beff(LAM, chi=CHI, kap=KAP)
        Wl = Wlin(xi, LAM, chi=CHI, kap=KAP, bounds=(0.8, 1.3))
        rows = []
        for f in pick([0.16, 0.20, 0.24, 0.28, 0.32, 0.36], [0.24]):
            p = params("magnet", xi=xi, chi=CHI, kap=KAP, lam=LAM, f=f, s=MAGNET_S)
            p[5] = Wl
            T = run_capture(x0, y0, phi, p, Wl, pick(4000, 500), 256, 0.5 * f / c)[:, 0]
            z = march(np.array([0.0, f / c * Wl, 0.0]), p, 1024, 300)
            z, ok, _ = newton_shoot(z, p, 1024)  # resonant orbit at Omega_lin
            ze, M, xmax, xmin, v2 = flow(z, p, 1024)
            rows.append(
                dict(
                    f=f,
                    cap=float(np.mean(T >= 0)),
                    A=max(xmax, -xmin),
                    Pe=Pe_from(v2, p),
                    mu=float(np.max(np.abs(np.linalg.eigvals(M)))),
                    ok=ok,
                )
            )
            log("xi", xi, "f", f, "done")
        out[xi] = dict(c=c, be=be, Wl=Wl, well=a, rows=rows)
    return out


def report(out):
    for xi, d in out.items():
        for r in d["rows"]:
            print(
                f"xi={xi} c={d['c']:.4f} Wlin={d['Wl']:.5f} f={r['f']}: capture {r['cap']:.3f}  A={r['A']:.3f} Pe={r['Pe']:.4f} maxmu={r['mu']:.4f} ok={r['ok']}"
            )


if __name__ == "__main__":
    main("design", compute, report)
