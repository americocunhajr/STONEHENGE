"""Harmonic content of the barrier force and of computed periodic orbits at the Gaussian peak,
against the closed form b_k / H(k Omega); and of the quartic (a, k_w) model below its jump
(Sect. 7, Fig. 7b).

Output: data/harmonics.pkl.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from harmonic_gaussian.model import BASE, params, dU, dyn_stiffness, linear_state, ctot
from harmonic_gaussian.orbits import newton_shoot
from harmonic_gaussian.theory import bk_exact, bk_asym
from harmonic_gaussian.branches import DUFF
from harmonic_gaussian.storage import main

CHI, KAP = BASE["chi"], BASE["kap"]


def orbit_harmonics(shape, s, f, W, kmax=41, xi=0.01, lam=0.05, z0=None):
    """Amplitudes |X_k| of the displacement harmonics of the periodic orbit (DOP853 over one period)."""
    p = params(shape, xi=xi, lam=lam, f=f, W=W, s=s)
    if z0 is None:
        z0 = linear_state(f, W, xi, lam)
    z, ok, _ = newton_shoot(z0, p, 1024)
    rhs = lambda t, y: [y[1], -2 * xi * y[1] - dU(y[0], p) + CHI * y[2] + f * np.cos(W * t), -lam * y[2] - KAP * y[1]]
    n = 1 << 13
    T = 2 * np.pi / W
    t = np.linspace(0, T, n, endpoint=False)
    sol = solve_ivp(rhs, (0, T), z, t_eval=t, method="DOP853", rtol=1e-12, atol=1e-13)
    return ok, (2 * np.abs(np.fft.rfft(sol.y[0])) / n)[: kmax + 1]


def compute():
    res = {}
    lines = []
    th = np.linspace(0, 2 * np.pi, 1 << 14, endpoint=False)  # closed form against FFT of U_m'(A cos th)
    for A in [3.0, 10.0]:
        Fk = 2 * np.fft.rfft(-2 * A * np.cos(th) * np.exp(-((A * np.cos(th)) ** 2))).real / len(th)
        lines.append(
            f"A {A} closed form vs FFT (k=1,3,5,11): "
            + str([(round(bk_exact(k, A), 6), round(Fk[k], 6)) for k in (1, 3, 5, 11)])
        )
    H = lambda w: dyn_stiffness(w, 0.01, 0.05)
    for f, W in [(0.15, 1.0106), (0.3, 1.0121), (0.6, 1.0123)]:
        ok, X = orbit_harmonics("gauss", (1, 1, 1, 0), f, W)
        A = X[1]
        ks = np.arange(3, 42, 2)
        pred = np.array([abs(bk_exact(k, A) / H(k * W)) for k in ks])
        pa = np.array([abs(bk_asym(k, A) / H(k * W)) for k in ks])
        res[("gauss", f)] = dict(W=W, A=A, ks=ks, X=X[ks], pred=pred, pred_asym=pa)
        lines.append(
            f"Gauss f={f} ok={ok} A1={A:.3f}  X3/X1={X[3]/A:.2e} pred {pred[0]/A:.2e}; X11/X1={X[11]/A:.2e} pred {pred[4]/A:.2e}; X21/X1={X[21]/A:.2e} pred {pred[9]/A:.2e}"
        )
    k1, k3 = DUFF["(a,k_w)"]
    f = 0.3
    W = 3.0  # quartic (a, k_w) just below its jump
    ke = CHI * KAP * W**2 / (0.05**2 + W**2)
    cW = ctot(0.01, 0.05, W) * W
    Gd = lambda A: (-k1 + 0.75 * k3 * A * A + ke - W * W) - np.sqrt(max(f * f / (A * A) - cW * cW, 0))
    Ab = np.sqrt((W * W - ke + k1) / (0.75 * k3))
    A = brentq(Gd, Ab, Ab * 1.2)
    phi = np.arctan2(cW, -k1 + 0.75 * k3 * A * A + ke - W * W)
    Xc = A * np.exp(-1j * phi)
    Vc = -1j * KAP * W * Xc / (0.05 + 1j * W)
    lines.append(f"Duffing HB amplitude {A}")
    ok, X = orbit_harmonics(
        "duffing", (k1, k3, 0, 0), f, W, z0=np.array([A * np.cos(phi), A * W * np.sin(phi), Vc.real])
    )
    ks = np.arange(3, 42, 2)
    res[("duff", f)] = dict(W=W, A=X[1], ks=ks, X=X[ks])
    lines.append(
        f"Duffing (a,k_w) f=0.3 W=3.0 ok={ok} A1={X[1]:.3f} X3/X1={X[3]/X[1]:.2e} X5/X1={X[5]/X[1]:.2e} X7/X1={X[7]/X[1]:.2e} X11/X1={X[11]/X[1]:.2e}"
    )
    res["summary"] = lines
    return res


def report(res):
    print("\n".join(res["summary"]))


if __name__ == "__main__":
    main("harmonics", compute, report)
