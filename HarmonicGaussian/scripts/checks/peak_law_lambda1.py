"""Peak law at the near-optimal load lambda = 1 against the baseline lambda = 0.05 (Sect. 9.3):
peak of power by Newton shooting and direct maximization over Omega, against Eq. (16)."""

import numpy as np
from scipy.optimize import minimize_scalar
from harmonic_gaussian.model import params, beff, ctot, Wlin, linear_state
from harmonic_gaussian.orbits import newton_shoot, flow, Pe_from
from harmonic_gaussian.theory import peak_law

XI = 0.01


def Pe(lam, f, W, N=512):
    p = params("gauss", xi=XI, lam=lam, f=f, W=W)
    z, ok, _ = newton_shoot(linear_state(f, W, XI, lam), p, N)
    return Pe_from(flow(z, p, N)[4], p) if ok else np.nan


for lam in [1.0, 0.05]:
    be = beff(lam)
    c = ctot(XI, lam)
    Wl = Wlin(XI, lam, bounds=(0.9, 1.1))
    print(f"lambda={lam}: sqrt(beff)={np.sqrt(be):.5f} c={c:.5f} Omega_lin={Wl:.5f}")
    for f in [0.15, 0.3, 0.6]:
        r = minimize_scalar(
            lambda W: -Pe(lam, f, W), bounds=(Wl - 0.02, Wl + 0.005), method="bounded", options={"xatol": 1e-9}
        )
        print(f"   f={f}: peak {r.x:.5f} law {peak_law(f, Wl, c, be):.5f}  Pe {-r.fun:.4e}")
