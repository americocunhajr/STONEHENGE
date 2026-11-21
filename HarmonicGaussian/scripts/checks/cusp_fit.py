"""Hysteresis width near the cusp (Sect. 5.4): fit Omega_1^2 - Omega_2^2 = C (f_c - f)^p to the
turning points of data/cusp.pkl, and the harmonic-balance normal-form prediction (exponent 3/2)."""

import numpy as np
from scipy.optimize import curve_fit
from harmonic_gaussian.theory import keff_hb_gauss
from harmonic_gaussian.storage import load

C = [q for q in load("cusp") if q.get("ok") and len(q["folds"]) == 2]
f = np.array([q["f"] for q in C])
W1 = np.array([q["folds"][0] for q in C])
W2 = np.array([q["folds"][1] for q in C])
d = W1**2 - W2**2
m = f >= 1.45
fun = lambda f, Cc, fc, p: Cc * np.clip(fc - f, 1e-12, None) ** p
(Cc, fc, p), cov = curve_fit(fun, f[m], d[m], p0=[0.1, 1.7, 1.5])
err = np.sqrt(np.diag(cov))
print(f"continuation: f_c={fc:.4f}+-{err[1]:.4f}  exponent p={p:.3f}+-{err[2]:.3f}  C={Cc:.4f}")
A = np.linspace(1.27, 4, 40001)
h = A**2 * np.gradient(keff_hb_gauss(A), A)
i = np.argmax(h)
Ac = A[i]
hpp = np.gradient(np.gradient(h, A), A)[i]
CHB = 4 * np.sqrt(2) / (3 * Ac**2 * np.sqrt(abs(hpp)))
print(f"HB: f_max={h[i]:.4f} A_c={Ac:.3f} h''={hpp:.3f} prefactor {CHB:.4f} (exponent 3/2)")
(Cc2, fc2), _ = curve_fit(lambda f, Cc, fc: Cc * np.clip(fc - f, 1e-12, None) ** 1.5, f[m], d[m], p0=[0.1, 1.7])
print(f"fit with p=3/2 fixed: C={Cc2:.4f} f_c={fc2:.4f}")
for q in C:
    print(q["f"], np.round(q["folds"], 5))
