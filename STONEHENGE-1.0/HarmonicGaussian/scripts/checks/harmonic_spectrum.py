"""Harmonics of the barrier force (Sect. 4, Eq. 14): b_k ~ -(-1)^((k-1)/2) (2k/(pi A^2)) FT[U_m](k/A)
against exact Fourier coefficients (FFT) for the Gaussian and the Lorentzian, and the crossover of
the Gaussian spectrum from algebraic to Gaussian decay at k ~ 2 A sqrt(alpha)."""

import numpy as np
from harmonic_gaussian.theory import bk_exact, M0_GAUSS as M0

N = 2**16
th = 2 * np.pi * np.arange(N) / N
ks = np.array([1, 5, 11, 21, 31, 41, 61, 81])
for name, Vp, Vh in [
    ("Lorentz", lambda x: -2 * x / (1 + x * x) ** 2, lambda q: np.pi * np.exp(-abs(q))),
    ("Gauss", lambda x: -2 * x * np.exp(-x * x), lambda q: np.sqrt(np.pi) * np.exp(-q * q / 4)),
]:
    for A in (10.0, 14.0, 20.0):
        c = np.fft.rfft(Vp(A * np.cos(th))) / N * 2
        ex = c[ks].real
        pr = -((-1) ** ((ks - 1) // 2)) * 2 * ks / (np.pi * A * A) * Vh(ks / A)
        print(name, A, " ".join(f"{k}:{e:.3e}/{p:.3e}({abs(p/e-1)*100:.1f}%)" for k, e, p in zip(ks, ex, pr)))
print("Gaussian crossover: |b_k| exact, algebraic 2kM0/(pi A^2), algebraic x exp(-k^2/(4A^2))")
for A in (7.0, 14.0):
    print(f"A={A}")
    for k in (1, 5, 11, 21, 31, 41, 61, 81):
        alg = 2 * k * M0 / (np.pi * A * A)
        print(f"   {k:3d}  {abs(bk_exact(k, A)):.3e}  {alg:.3e}  {alg*np.exp(-k*k/(4*A*A)):.3e}")
