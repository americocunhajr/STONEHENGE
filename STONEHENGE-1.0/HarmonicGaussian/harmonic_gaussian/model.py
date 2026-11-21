"""Model of the piezoelectric harvester with a localized magnetic barrier.

Nondimensional equations (Sect. 2 of the paper)::

    x'' + 2 xi x' + U'(x) - chi v = f cos(Omega t)
    v'  + lam v + kap x'          = 0

with U(x) = beta x^2/2 + U_m(x) + k3h x^4/4. The barrier U_m is one of

    'gauss'    gamma exp(-alpha (x - delta)^2)                     s = (alpha, beta, gamma, -)
    'lorentz'  gamma / (1 + u^2),            u = (x - delta)/w     s = (w, beta, gamma, -)
    'dipole'   gamma / (1 + u^2)^(3/2)                             s = (w, beta, gamma, -)
    'magnet'   point-dipole pair with tip rotation, Eq. (7):
               K [(2 - u^2) cos(eta u) - 3 u sin(eta u)] / (1 + u^2)^(5/2)   s = (w, beta, K, eta)
    'duffing'  quartic potential, U' = -k1 x + k3 x^3                         s = (k1, k3, -, -)

All parameters travel in one float vector ``p`` (layout below) so that the compiled kernels of
:mod:`orbits` and :mod:`ensembles` take a single array argument.
"""

import numpy as np
from numba import njit
from scipy.optimize import brentq, minimize_scalar

# ---------------------------------------------------------------- parameter vector
#: positions in the parameter vector
I_XI, I_CHI, I_KAP, I_LAM, I_F, I_W, I_SHAPE, I_S1, I_S2, I_S3, I_S4, I_K3H, I_DEL = range(13)
SHAPES_ID = {"gauss": 0, "lorentz": 1, "dipole": 2, "duffing": 3, "magnet": 4}

#: baseline coupling, circuit and damping (Sect. 2.5)
BASE = dict(xi=0.01, chi=0.05, kap=0.5, lam=0.05)
#: point-dipole magnet pair used throughout the paper (K = 0.3, eta = 0.046)
MAGNET_S = (1.0, 1.0, 0.3, 0.046)
#: well position of the baseline Gaussian barrier (alpha = beta = gamma = 1): a = sqrt(ln 2)
A_WELL = np.sqrt(np.log(2))


def params(
    shape="gauss", xi=0.01, chi=0.05, kap=0.5, lam=0.05, f=0.083, W=1.0, s=(1.0, 1.0, 1.0, 0.0), k3h=0.0, delta=0.0
):
    """Parameter vector for the compiled kernels (default: baseline Gaussian barrier)."""
    p = np.zeros(13)
    p[0:6] = xi, chi, kap, lam, f, W
    p[6] = SHAPES_ID[shape]
    p[7:11] = s
    p[11] = k3h
    p[12] = delta
    return p


# ---------------------------------------------------------------- restoring force
@njit(cache=True)
def dU(x, p):
    """Restoring force U'(x)."""
    sid = int(p[I_SHAPE])
    k3h = p[I_K3H]
    if sid == 3:  # Duffing: U' = -k1 x + k3 x^3
        return -p[I_S1] * x + p[I_S2] * x**3
    d = x - p[I_DEL]
    if sid == 0:  # Gaussian: s = (alpha, beta, gamma)
        a, b, g = p[I_S1], p[I_S2], p[I_S3]
        return b * x - 2 * a * g * d * np.exp(-a * d * d) + k3h * x**3
    w, b, g = p[I_S1], p[I_S2], p[I_S3]
    u = d / w
    if sid == 1:  # Lorentzian g/(1+u^2)
        return b * x - 2 * g * u / (w * (1 + u * u) ** 2) + k3h * x**3
    if sid == 4:  # point-dipole pair with tip rotation, s = (w, beta, K, eta)
        e = p[I_S4]
        ce = np.cos(e * u)
        se = np.sin(e * u)
        g1 = (
            5 * u * (3 * u * se + (u * u - 2) * ce)
            + (u * u + 1) * (-3 * e * u * ce + e * (u * u - 2) * se - 2 * u * ce - 3 * se)
        ) / (u * u + 1) ** 3.5
        return b * x + g * g1 / w + k3h * x**3
    return b * x - 3 * g * u / (w * (1 + u * u) ** 2.5) + k3h * x**3  # dipole-type g/(1+u^2)^1.5


@njit(cache=True)
def d2U(x, p):
    """Stiffness U''(x)."""
    sid = int(p[I_SHAPE])
    k3h = p[I_K3H]
    if sid == 3:
        return -p[I_S1] + 3 * p[I_S2] * x * x
    d = x - p[I_DEL]
    if sid == 0:
        a, b, g = p[I_S1], p[I_S2], p[I_S3]
        return b - 2 * a * g * np.exp(-a * d * d) * (1 - 2 * a * d * d) + 3 * k3h * x * x
    w, b, g = p[I_S1], p[I_S2], p[I_S3]
    u = d / w
    if sid == 1:
        return b - 2 * g * (1 - 3 * u * u) / (w * w * (1 + u * u) ** 3) + 3 * k3h * x * x
    if sid == 4:
        e = p[I_S4]
        ce = np.cos(e * u)
        se = np.sin(e * u)
        q = u * u + 1
        g2 = (
            10 * u * q * (3 * e * u * ce - e * (u * u - 2) * se + 2 * u * ce + 3 * se)
            + (5 - 30 * u * u) * (3 * u * se + (u * u - 2) * ce)
            + q * q * (3 * e * e * u * se + e * e * (u * u - 2) * ce + 4 * e * u * se - 6 * e * ce - 2 * ce)
        ) / q**4.5
        return b + g * g2 / (w * w) + 3 * k3h * x * x
    return b - 3 * g * (1 - 4 * u * u) / (w * w * (1 + u * u) ** 3.5) + 3 * k3h * x * x


def magnet_barrier(u, K=0.3, eta=0.046):
    """Point-dipole barrier U_m(u) of Eq. (7) (numpy, vectorized)."""
    return K * ((2 - u * u) * np.cos(eta * u) - 3 * u * np.sin(eta * u)) / (1 + u * u) ** 2.5


def wells(p, lo=0.3, hi=3.0):
    """Left and right well positions (zeros of U') of a symmetric or offset barrier."""
    d = p[I_DEL]
    return (brentq(lambda x: dU(x, p), -hi, -lo + d), brentq(lambda x: dU(x, p), lo + d, hi))


# ---------------------------------------------------------------- linear beam and circuit
def ctot(xi, lam, W=1.0, chi=0.05, kap=0.5):
    """Total damping c = 2 xi + chi kap lam / (lam^2 + W^2), Eq. (13)."""
    return 2 * xi + chi * kap * lam / (lam**2 + W**2)


def beff(lam, beta=1.0, chi=0.05, kap=0.5):
    """Electrically stiffened beam stiffness beta_eff = beta + chi kap beta / (lam^2 + beta)."""
    return beta + chi * kap * beta / (lam**2 + beta)


def dyn_stiffness(W, xi=0.01, lam=0.05, beta=1.0, chi=0.05, kap=0.5):
    """Complex dynamic stiffness of the linear harvester, beta - W^2 + 2i xi W + i chi kap W/(lam + iW)."""
    return beta - W**2 + 2j * xi * W + chi * kap * 1j * W / (lam + 1j * W)


def P_lin(W, f, xi=0.01, lam=0.05, chi=0.05, kap=0.5, beta=1.0):
    """Mean electrical power of the linear harvester (no magnets)."""
    X = f / np.abs(dyn_stiffness(W, xi, lam, beta, chi, kap))
    return chi * lam * kap * W**2 * X**2 / (2 * (lam**2 + W**2))


def Wlin(xi, lam, beta=1.0, chi=0.05, kap=0.5, bounds=None):
    """Frequency of maximum power of the linear harvester, Omega_lin."""
    mP = lambda W: -(W**2) / abs(dyn_stiffness(W, xi, lam, beta, chi, kap)) ** 2 / (lam**2 + W**2)
    if bounds is None:
        bounds = (0.8 * np.sqrt(beta), 1.2 * np.sqrt(beta))
    return minimize_scalar(mP, bounds=bounds, method="bounded", options={"xatol": 1e-12}).x


def linear_state(f, W, xi=0.01, lam=0.05, chi=0.05, kap=0.5, beta=1.0):
    """State (x, x', v) at phase 0 of the linear periodic response; an initial guess for shooting."""
    X = f / dyn_stiffness(W, xi, lam, beta, chi, kap)
    V = -kap * 1j * W * X / (lam + 1j * W)
    return np.array([X.real, (1j * W * X).real, V.real])


def power_factor(p):
    """Factor chi lam / kap that turns the mean of v^2 into the mean electrical power."""
    return p[I_CHI] * p[I_LAM] / p[I_KAP]
