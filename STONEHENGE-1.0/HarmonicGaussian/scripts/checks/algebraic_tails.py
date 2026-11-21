"""Area law for barriers with algebraic tails (Sect. 4, Supplementary Sect. S1.2): residual of the
exact backbone after the leading term 2M0/(pi A^3), its decay exponent, and the coefficient
(4 Lambda/pi) G_p of the next term, with 40-digit quadrature."""

import mpmath as mp
from harmonic_gaussian.theory import period_mp

mp.mp.dps = 40


def residual(V, M0, A):
    w = 2 * mp.pi / period_mp(V, mp.mpf(A))
    return (1 - w * w) - 2 * M0 / (mp.pi * mp.mpf(A) ** 3)


cases = [
    ("Lorentzian  p=2", lambda x: 1 / (1 + x * x), mp.pi),
    ("p=3 type     ", lambda x: (1 + x * x) ** mp.mpf(-1.5), mp.mpf(2)),
    ("dipole pair, tail p=3 (negative)", lambda x: (2 - x * x) / (1 + x * x) ** mp.mpf(2.5), mp.mpf(2)),
    ("Gaussian     ", lambda x: mp.e ** (-x * x), mp.sqrt(mp.pi)),
]
for name, V, M0 in cases:
    As = [10, 20, 40, 80, 160]
    res = [residual(V, M0, A) for A in As]
    sl = [float(mp.log(abs(res[i + 1] / res[i])) / mp.log(2)) for i in range(len(As) - 1)]
    print(name, " resid:", " ".join(f"{float(r):+.2e}" for r in res), " slopes:", " ".join(f"{s:.2f}" for s in sl))
print("dipole-pair area:", mp.quad(lambda x: (2 - x * x) / (x * x + 1) ** mp.mpf(2.5), [-mp.inf, 0, mp.inf]))
print("--- intermediate tails, predicted slope -(p+2) ---")
for p in (1.5, 2.5):
    V = lambda x, p=p: (1 + x * x) ** mp.mpf(-p / 2)
    M0 = mp.quad(V, [-mp.inf, 0, mp.inf])
    As = [20, 40, 80, 160, 320]
    res = [residual(V, M0, A) for A in As]
    print(p, [f"{float(mp.log(abs(res[i+1]/res[i]))/mp.log(2)):.2f}" for i in range(len(As) - 1)])
print("--- p=3: resid*A^5/log A ---")
V = lambda x: (1 + x * x) ** mp.mpf(-1.5)
for A in (40, 160, 640, 2560):
    r = residual(V, mp.mpf(2), A)
    print(A, float(r * mp.mpf(A) ** 5), float(r * mp.mpf(A) ** 5 / mp.log(A)))
print("--- next-term coefficient: predicted (4/pi) G_p, G_p = -sqrt(pi) Gamma((1-p)/2)/Gamma(-p/2) ---")
for p in (1.5, 2.5):
    V = lambda x, p=p: (1 + x * x) ** mp.mpf(-p / 2)
    M0 = mp.quad(V, [-mp.inf, 0, mp.inf])
    Gp = -mp.sqrt(mp.pi) * mp.gamma((1 - p) / 2) / mp.gamma(-p / 2)
    for A in (80, 320, 1280):
        print(p, A, "coef num", float(residual(V, M0, A) * mp.mpf(A) ** (p + 2)), " pred", float(4 * Gp / mp.pi))
