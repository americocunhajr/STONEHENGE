"""Next-order term of the exact interwell backbone of x^2/2 + exp(-x^2) (Sect. 4):
beta - w^2 = 2M0/(pi A^3) + (3/pi)(M2 + N2/beta) A^-5 + O(A^-6), checked with 50-digit quadrature
after removing the known A^-6 term 3 M0^2/(pi^2 A^6)."""

import mpmath as mp
from harmonic_gaussian.theory import period_mp

mp.mp.dps = 50
M0 = mp.sqrt(mp.pi)
M2 = M0 / 2
N2 = mp.sqrt(mp.pi / 2)
V = lambda x: mp.e ** (-x * x)
pred = -(3 / mp.pi) * (M2 + N2)
rows = []
for A in [4, 6, 8, 10, 12]:
    A = mp.mpf(A)
    w = 2 * mp.pi / period_mp(V, A)
    r = (w**2 - 1 + 2 * M0 / (mp.pi * A**3) - 3 * M0**2 / (mp.pi**2 * A**6)) * A**5
    rows.append((float(A), float(r)))
    print(
        f"A={float(A):4.1f}  scaled residual {mp.nstr(r, 6)}   prediction {mp.nstr(pred, 6)}   (HB-only coefficient -3M2/pi = {mp.nstr(-(3/mp.pi)*M2, 5)})"
    )
(A1, r1), (A2, r2) = rows[-2], rows[-1]
print("Richardson extrapolation (error ~ A^-2):", (r2 * A2**2 - r1 * A1**2) / (A2**2 - A1**2))
