"""Melnikov threshold for transverse intersection of the manifolds of the hilltop saddle, with
weak damping and coupling (Sect. 8.3; orange curve in Fig. 8a).

M(t0) = int xd_h [f cos W(t + t0) - 2 xi xd_h + chi v_h] dt has simple zeros iff
f |S(W)| > 2 xi I0 + E_e, where xd_h is the velocity on the homoclinic loop, S its Fourier
transform, I0 = int xd_h^2 dt, and E_e the energy taken by the circuit driven along the loop.
Output: data/melnikov.pkl (thresholds f_M(Omega) for the symmetric and the offset barrier).
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from harmonic_gaussian.model import BASE, ctot
from harmonic_gaussian.storage import main

CHI, KAP = BASE["chi"], BASE["kap"]


def loops(alpha=1.0, beta=1.0, gamma=1.0, delta=0.0, T=16.0, n=40001):
    """Homoclinic loops of the hilltop saddle on either side, as functions of time."""
    Up = lambda x: beta * x - 2 * alpha * gamma * (x - delta) * np.exp(-alpha * (x - delta) ** 2)
    U = lambda x: beta * x * x / 2 + gamma * np.exp(-alpha * (x - delta) ** 2)
    xs = brentq(Up, delta - 0.5, delta + 0.5) if delta != 0 else 0.0
    out = {}
    for side in (+1, -1):
        xt = brentq(lambda x: U(x) - U(xs), xs + side * 1e-3, xs + side * 5.0)  # turning point of the loop
        t = np.linspace(0, T, n)
        sol = solve_ivp(
            lambda t, y: [y[1], -Up(y[0])], (0, T), [xt, 0.0], t_eval=t, method="DOP853", rtol=1e-12, atol=1e-14
        )
        x = sol.y[0]
        v = sol.y[1]
        tt = np.r_[-t[::-1], t[1:]]
        xh = np.r_[x[::-1], x[1:]]
        vh = np.r_[-v[::-1], v[1:]]  # even about the turning point
        I0q = 2 * quad(lambda s: np.sqrt(max(2 * (U(xs) - U(s)), 0)), min(xs, xt), max(xs, xt), limit=400)[0]
        out[side] = dict(t=tt, x=xh, xd=vh, I0=np.trapezoid(vh**2, tt), I0q=I0q, xt=xt)
    return xs, out


def threshold(L, xi, lam, Ws):
    """Dissipation 2 xi I0 + E_e, its electrical part E_e, and |S(W)| along a loop."""
    t, xd = L["t"], L["xd"]
    dt = t[1] - t[0]
    v = np.zeros_like(t)
    e = np.exp(-lam * dt)
    for i in range(1, len(t)):
        v[i] = v[i - 1] * e - KAP * 0.5 * dt * (xd[i] + xd[i - 1] * e)  # v' = -lam v - kap xd_h
    Ee = -CHI * np.trapezoid(xd * v, t)
    D = 2 * xi * L["I0"] + Ee
    S = np.array([abs(np.trapezoid(xd * np.exp(1j * W * t), t)) for W in Ws])
    return D, Ee, S


def compute():
    Ws = np.linspace(0.2, 1.6, 281)
    res = {}
    lines = []
    xs, L = loops()
    lines.append(f"symmetric: I0 = {L[1]['I0']}  quadrature {L[1]['I0q']}  turning point {L[1]['xt']}")
    for xi in [0.005, 0.01, 0.02]:
        D, Ee, S = threshold(L[1], xi, 0.05, Ws)
        I0 = L[1]["I0"]
        cM = D / I0
        c = ctot(xi, 0.05)
        R = I0 / S
        fM = cM * R
        i = np.argmin(fM)
        res[("sym", xi)] = dict(Ws=Ws, fM=fM, R=R, cM=cM, c=c)
        lines.append(
            f"xi={xi}: c={c:.5f} c_M={cM:.5f} (electrical part {Ee/I0:.5f})  min f_M={fM[i]:.4f} at W={Ws[i]:.3f}  f_M/c={fM[i]/c:.2f}  R_min={R[i]:.3f};  at W=1: f_M/c={fM[np.argmin(abs(Ws-1))]/c:.2f}"
        )
    xs, Lb = loops(delta=0.1)
    for side in (+1, -1):
        D, Ee, S = threshold(Lb[side], 0.01, 0.05, Ws)
        I0 = Lb[side]["I0"]
        fM = D / S
        i = np.argmin(fM)
        c = ctot(0.01, 0.05)
        res[("bias", side)] = dict(Ws=Ws, fM=fM, R=I0 / S, cM=D / I0, c=c)
        lines.append(
            f'delta=0.1 loop {"right" if side > 0 else "left"}: saddle {xs:.4f} turning {Lb[side]["xt"]:.4f} I0={I0:.4f} min f_M/c={fM[i]/c:.2f} at W={Ws[i]:.3f}; at W=1: {fM[np.argmin(abs(Ws-1))]/c:.2f}'
        )
    res["summary"] = lines
    return res


def report(res):
    print("\n".join(res["summary"]))


if __name__ == "__main__":
    main("melnikov", compute, report)
