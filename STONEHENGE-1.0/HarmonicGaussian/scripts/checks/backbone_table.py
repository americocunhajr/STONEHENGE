"""Supplementary Table S1 (accuracy of the backbone approximations) and static numbers of the
baseline Gaussian barrier quoted in Sects. 2-5 and 9."""

import numpy as np
import mpmath as mp
from scipy.optimize import brentq
from harmonic_gaussian.model import beff, ctot
from harmonic_gaussian.theory import omega_inter, keff_hb_gauss, Q_gauss, M0_GAUSS as M0

U = lambda x: np.exp(-x * x) + x * x / 2
M2 = M0 / 2
print("A   exact   HB   errHB  errLead  errTwo   (relative errors in deficit omega^2-beta)")
for A in [1.6, 2, 3, 4, 6, 8]:
    w = omega_inter(U, A)
    d = w * w - 1
    dHB = keff_hb_gauss(A) - 1
    dL = -2 * M0 / (np.pi * A**3)
    dT = -(2 / (np.pi * A**3)) * (M0 + 1.5 * M2 / A**2)
    print(
        f"{A:<4}{w:.4f} {np.sqrt(keff_hb_gauss(A)):.4f} {(dHB-d)/d*100:+.0f}% {(dL-d)/d*100:+.0f}% {(dT-d)/d*100:+.0f}%"
    )
print("separatrix turning point", brentq(lambda A: U(A) - 1.0, 0.5, 5))
L2 = np.log(2)
print("a", np.sqrt(L2), "k_w", 2 * L2, "h", (1 - L2) / 2, "Q", Q_gauss(2.0))
A = np.linspace(1.27, 6, 20000)
g = A**2 * np.gradient(keff_hb_gauss(A), A)
i = np.argmax(g)
print("Gauss fmax", g[i], "at A", A[i])
print("sqrt k_w", np.sqrt(2 * L2))
s = brentq(lambda s: np.exp(s * s) / (1 - 2 * s * s) - 2, 1e-6, 0.7)
print("bias boundary Delta at r=2", 2 * s**3 / (1 - 2 * s * s))
be = beff(0.05)
c = ctot(0.01, 0.05)
sb = np.sqrt(be)
fl = (2 * M0 / np.pi) ** (1 / 3) * (c * sb) ** (2 / 3)
print("f_l", fl, "f_l/c", fl / c)
k3s = (4 / 3) * (c * sb) ** (5 / 3) * (np.pi / (2 * M0)) ** (2 / 3)
print("k3* ", k3s, " c/k3*", c / k3s)
print(
    "A* at f=0.6",
    np.sqrt(6 * M0 / (np.pi * 0.6)),
    "fold law f=0.6",
    np.sqrt(be + (2 / 3) * np.sqrt(np.pi / (6 * M0)) * 0.6**1.5),
)
print("lambda=1: sqrt beff", np.sqrt(beff(1.0)), "c", ctot(0.01, 1.0))
Qf = lambda y: (mp.e**y - 1 - y) / (4 * y**2)  # conditioning of the inverse problem from Q
for r in [1.2, 2.0]:
    y = mp.log(r)
    dyq = Qf(y) / mp.diff(Qf, y)
    print(f"r={r}: amplification into alpha,beta {float(abs(dyq/y)):.1f}, into gamma {float(abs((1-2/y)*dyq)):.1f}")
print(
    "(a,h): k1",
    4 * 0.5 * (1 - L2) / L2,
    "k3",
    4 * 0.5 * (1 - L2) / L2**2,
    "k_w",
    2 * 4 * 0.5 * (1 - L2) / L2,
    " (a,k_w): h",
    L2**2 / 4,
)
