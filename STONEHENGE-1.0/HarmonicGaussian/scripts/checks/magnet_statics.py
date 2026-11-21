"""Static properties of the point-dipole magnet pair (Sect. 2.3): stiffness at the top, well
position a, barrier height h, well stiffness k_w, Q = h/(k_w a^2) and area M0, without and with
tip rotation (eta = 0.046), for four magnet strengths K."""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from harmonic_gaussian.model import magnet_barrier

for e in (0.0, 0.046):
    for K in (0.1, 0.15, 0.2, 0.3):
        U = lambda x: x * x / 2 + magnet_barrier(x, K, e)
        Up = lambda x, h=1e-6: (U(x + h) - U(x - h)) / (2 * h)
        try:
            a = brentq(Up, 0.05, 3)
        except Exception:
            a = np.nan
        h = U(0) - U(a) if np.isfinite(a) else np.nan
        kw = (U(a + 1e-4) - 2 * U(a) + U(a - 1e-4)) / 1e-8 if np.isfinite(a) else np.nan
        M0 = quad(lambda u: magnet_barrier(u, K, e), -60, 60, limit=400)[0]
        print(
            f"eta={e} K={K}: U''(0)={1-K*2*(6+3*e+e*e):+.3f} a={a:.3f} h={h:.4f} k_w={kw:.3f} Q={h/(kw*a*a) if np.isfinite(a) else np.nan:.3f} M0={M0:.4f}"
        )
print("Gaussian baseline: a=0.833 h=0.153 M0=1.772")
