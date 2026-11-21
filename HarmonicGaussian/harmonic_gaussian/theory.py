"""Closed-form and quadrature results used by several scripts: free-oscillation frequencies,
single-harmonic stiffness, harmonics of the barrier force, the Q invariant, and the similarity
parameter with the peak law.
"""

import numpy as np
import mpmath as mp
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ive

M0_GAUSS = np.sqrt(np.pi)  # area of the baseline barrier exp(-x^2)


def area_gauss(alpha=1.0, gamma=1.0):
    """Area M0 = gamma sqrt(pi/alpha) of the Gaussian barrier."""
    return gamma * np.sqrt(np.pi / alpha)


# ---------------------------------------------------------------- free oscillations
def omega_inter(U, A):
    """Angular frequency of the free interwell orbit with turning point A (scipy quadrature)."""
    E = U(A)
    f = lambda th: A * np.cos(th) / np.sqrt(max(2 * (E - U(A * np.sin(th))), 1e-300))
    return 2 * np.pi / (4 * quad(f, 0, np.pi / 2, limit=800, epsabs=1e-14, epsrel=1e-13)[0])


def omega_intra(U, xw, xs, E):
    """Angular frequency and outer turning point of the free intrawell orbit of energy E
    in the well at xw (barrier top at xs)."""
    x1 = brentq(lambda x: U(x) - E, xs + 1e-12, xw)
    x2 = brentq(lambda x: U(x) - E, xw, 20)
    c = (x1 + x2) / 2
    h = (x2 - x1) / 2
    f = lambda th: h * np.cos(th) / np.sqrt(max(2 * (E - U(c + h * np.sin(th))), 1e-300))
    return 2 * np.pi / (2 * quad(f, -np.pi / 2, np.pi / 2, limit=800)[0]), x2


def period_mp(V, A):
    """Period of x'' + x + V'(x) = 0 with turning point A, in mpmath arithmetic.

    With x = A sin(th), T = 4 int_0^{pi/2} dth / sqrt(1 + 2[V(A) - V(A sin th)]/(A^2 cos^2 th)),
    which has no cancellation at the turning point."""
    VA = V(A)
    lim = 1 / mp.sqrt(1 + mp.diff(V, A) / A)

    def g(th):
        c = mp.cos(th)
        if abs(c) < mp.mpf("1e-12"):
            return lim
        return 1 / mp.sqrt(1 + 2 * (VA - V(A * mp.sin(th))) / (A * A * c * c))

    return 4 * mp.quad(g, [0, mp.pi / 8, mp.pi / 4, 3 * mp.pi / 8, mp.pi / 2])


def keff_hb_gauss(A):
    """Single-harmonic stiffness of x^2/2 + exp(-x^2): 1 - 2 e^{-z}[I0(z) - I1(z)], z = A^2/2."""
    return 1 - 2 * (ive(0, A * A / 2) - ive(1, A * A / 2))


def Q_gauss(r):
    """Q = h/(k_w a^2) of the Gaussian barrier as a function of the stiffness ratio r = 2 alpha gamma/beta."""
    return (r - 1 - np.log(r)) / (4 * np.log(r) ** 2)


# ---------------------------------------------------------------- harmonics of the barrier force
def bk_exact(k, A, alpha=1.0, gamma=1.0):
    """Fourier coefficient of order k of -U_m'(A cos th) for the Gaussian barrier (Bessel form)."""
    m = (k - 1) // 2
    z = alpha * A * A / 2
    return -2 * alpha * gamma * A * (-1) ** m * (ive(m, z) - ive(m + 1, z))


def bk_asym(k, A, M0=M0_GAUSS):
    """Large-amplitude, low-order limit of bk: -(-1)^((k-1)/2) (2/pi) k M0 / A^2."""
    return -((-1) ** ((k - 1) // 2)) * (2 / np.pi) * k * M0 / A**2


# ---------------------------------------------------------------- similarity parameter and peak law
def epsilon(f, c, be, M0=M0_GAUSS):
    """Similarity parameter eps = 2 beta_eff M0 c^2 / (pi f^3), Eq. (18)."""
    return 2 * be * M0 * c * c / (np.pi * f**3)


def peak_law(f, Wl, c, be, M0=M0_GAUSS):
    """Frequency of maximum power from the area law, Eq. (16): sqrt(Wlin^2 - (2M0/pi)(c sqrt(beta_eff)/f)^3)."""
    return np.sqrt(Wl**2 - (2 * M0 / np.pi) * (c * np.sqrt(be) / f) ** 3)
