"""Periodic orbits: period map with variational equations, Newton shooting, pseudo-arclength
continuation in the forcing frequency, Floquet multipliers, and continuation of the fold curve
in the (Omega, f) plane.

The period map integrates one forcing period in the rescaled time tau = Omega t with a fixed-step
RK4 scheme (N steps per period) together with the 3 x 3 variational equations, so that one call
returns the end state and the monodromy matrix.
"""

import sys
import numpy as np
from numba import njit
from .model import dU, d2U, I_W, I_F, I_CHI, I_LAM, I_KAP


@njit(cache=True)
def _rhs(s, tau, p, out):
    # s = (x, y, v, Phi row-major); derivative with respect to tau = Omega t
    W = p[I_W]
    x, y, v = s[0], s[1], s[2]
    xi, chi, kap, lam, f = p[0], p[1], p[2], p[3], p[4]
    out[0] = y / W
    out[1] = (-2 * xi * y - dU(x, p) + chi * v + f * np.cos(tau)) / W
    out[2] = (-lam * v - kap * y) / W
    k = d2U(x, p)
    # J = [[0,1,0],[-k,-2xi,chi],[0,-kap,-lam]] / W ; dPhi = J Phi
    for j in range(3):
        p0 = s[3 + j]
        p1 = s[6 + j]
        p2 = s[9 + j]
        out[3 + j] = p1 / W
        out[6 + j] = (-k * p0 - 2 * xi * p1 + chi * p2) / W
        out[9 + j] = (-kap * p1 - lam * p2) / W


@njit(cache=True)
def flow(z, p, N):
    """Integrate one forcing period from state z = (x, x', v).

    Returns the end state, the monodromy matrix, max x, min x and the time average of v^2."""
    s = np.zeros(12)
    s[0:3] = z
    s[3] = 1.0
    s[7] = 1.0
    s[11] = 1.0
    h = 2 * np.pi / N
    k1 = np.zeros(12)
    k2 = np.zeros(12)
    k3 = np.zeros(12)
    k4 = np.zeros(12)
    tmp = np.zeros(12)
    xmax = z[0]
    xmin = z[0]
    v2 = 0.0
    for i in range(N):
        t = i * h
        _rhs(s, t, p, k1)
        for q in range(12):
            tmp[q] = s[q] + 0.5 * h * k1[q]
        _rhs(tmp, t + 0.5 * h, p, k2)
        for q in range(12):
            tmp[q] = s[q] + 0.5 * h * k2[q]
        _rhs(tmp, t + 0.5 * h, p, k3)
        for q in range(12):
            tmp[q] = s[q] + h * k3[q]
        _rhs(tmp, t + h, p, k4)
        v_old = s[2]
        for q in range(12):
            s[q] += h / 6 * (k1[q] + 2 * k2[q] + 2 * k3[q] + k4[q])
        v2 += 0.5 * (v_old * v_old + s[2] * s[2]) * h
        if s[0] > xmax:
            xmax = s[0]
        if s[0] < xmin:
            xmin = s[0]
    M = s[3:12].reshape(3, 3).copy()
    return s[0:3].copy(), M, xmax, xmin, v2 / (2 * np.pi)


@njit(cache=True)
def march(z, p, N, nper):
    """Integrate nper forcing periods (transient) and return the final state."""
    for _ in range(nper):
        z, M, a, b, c = flow(z, p, N)
    return z


def Pe_from(v2, p):
    """Mean electrical power from the time average of v^2."""
    return p[I_CHI] * p[I_LAM] / p[I_KAP] * v2


def newton_shoot(z, p, N, tol=1e-11, maxit=30):
    """Newton iteration for a fixed point of the period map. Returns (z, converged, iterations)."""
    for it in range(maxit):
        ze, M, *_ = flow(z, p, N)
        G = ze - z
        if np.max(np.abs(G)) < tol:
            return z, True, it
        dz = np.linalg.solve(M - np.eye(3), -G)
        z = z + dz
    return z, False, maxit


def periodic_orbit(z0, p, N=512):
    """Converged periodic orbit through z0: dict with state, amplitude, power and Floquet multipliers."""
    z, ok, _ = newton_shoot(z0, p, N)
    ze, M, xmax, xmin, v2 = flow(z, p, N)
    return dict(z=z, ok=ok, A=max(xmax, -xmin), xmax=xmax, xmin=xmin, Pe=Pe_from(v2, p), mu=np.linalg.eigvals(M))


def _GW(z, p, N, eps=1e-6):
    pp = p.copy()
    pp[I_W] += eps
    e1, *_ = flow(z, pp, N)
    pm = p.copy()
    pm[I_W] -= eps
    e2, *_ = flow(z, pm, N)
    return (e1 - e2) / (2 * eps)


def continuation(
    z0,
    p,
    N=512,
    ds=0.05,
    ds_min=1e-5,
    ds_max=0.5,
    direction=+1,
    Wmin=0.3,
    Wmax=3.0,
    Astop=None,
    maxpts=2000,
    sW=10.0,
    verbose=False,
    stop_after_fold=None,
):
    """Pseudo-arclength continuation of a periodic orbit in Omega.

    The unknowns are u = (z, sW * Omega); the scaling sW balances the state and the frequency.
    Stops when Omega leaves [Wmin, Wmax], when the amplitude drops below Astop, or after the given
    number of turning points. Returns a list of dicts (W, z, mu, xmax, xmin, A, Pe)."""
    p = p.copy()
    z, ok, _ = newton_shoot(z0, p, N)
    if not ok:
        raise RuntimeError("initial Newton failed")
    rec = []

    def record(z, p):
        ze, M, xmax, xmin, v2 = flow(z, p, N)
        rec.append(
            dict(
                W=p[I_W],
                z=z.copy(),
                mu=np.linalg.eigvals(M),
                xmax=xmax,
                xmin=xmin,
                A=max(xmax, -xmin),
                Pe=Pe_from(v2, p),
            )
        )

    def jac(z, p):
        ze, M, *_ = flow(z, p, N)
        return np.hstack([M - np.eye(3), _GW(z, p, N)[:, None] / sW])

    t = np.linalg.svd(jac(z, p))[2][-1]  # initial tangent: null space of [M - I, G_W/sW]
    if np.sign(t[3]) != direction:
        t = -t
    record(z, p)
    u = np.r_[z, sW * p[I_W]]
    for k in range(maxpts):
        success = False
        while not success:
            up = u + ds * t
            v = up.copy()
            for it in range(12):
                p[I_W] = v[3] / sW
                ze, M, *_ = flow(v[:3], p, N)
                G = ze - v[:3]
                Jz = np.hstack([M - np.eye(3), _GW(v[:3], p, N)[:, None] / sW])
                R = np.r_[G, t @ (v - up)]
                if np.max(np.abs(R)) < 1e-10:
                    break
                v = v - np.linalg.solve(np.vstack([Jz, t]), R)
            if np.max(np.abs(R)) < 1e-9 and np.all(np.isfinite(v)):
                success = True
            else:
                ds *= 0.5
                if ds < ds_min:
                    if verbose:
                        print("step too small, stop at W", u[3] / sW, file=sys.stderr)
                    return rec
        tn = np.linalg.svd(Jz)[2][-1]
        if tn @ t < 0:
            tn = -tn
        t = tn
        u = v
        p[I_W] = u[3] / sW
        record(u[:3], p)
        if it <= 3:
            ds = min(ds * 1.4, ds_max)
        W_ = u[3] / sW
        if W_ < Wmin or W_ > Wmax:
            break
        if Astop is not None and rec[-1]["A"] < Astop:
            break
        if stop_after_fold is not None and len(rec) > 3:
            Ws = [r["W"] for r in rec]
            if any((Ws[i] - Ws[i - 1]) * (Ws[i + 1] - Ws[i]) < 0 for i in range(1, len(Ws) - 1)):
                stop_after_fold -= 1
                if stop_after_fold <= 0:
                    break
    return rec


def analyse(rec):
    """Arrays (W, A, Pe, |mu| sorted, unstable mask, indices of turning points) of a continuation record."""
    W = np.array([r["W"] for r in rec])
    A = np.array([r["A"] for r in rec])
    Pe = np.array([r["Pe"] for r in rec])
    mu = np.array([np.sort(np.abs(r["mu"]))[::-1] for r in rec])
    unstable = mu[:, 0] > 1 + 1e-7
    folds = [i for i in range(1, len(W) - 1) if (W[i] - W[i - 1]) * (W[i + 1] - W[i]) < 0]
    return W, A, Pe, mu, unstable, folds


def multiplier_type(mu):
    """Type of the dominant Floquet multiplier at a loss of stability: 'PD', 'NS' or '+1'."""
    big = mu[np.argmax(np.abs(mu))]
    return "PD" if (big.real < 0 and abs(big.imag) < 1e-9) else ("NS" if abs(big.imag) > 1e-9 else "+1")


def refine_peak(rec, p, N=512, key="Pe"):
    """Golden-section refinement of the maximum of Pe (or A) near the best point of a branch.
    Returns (Omega_peak, maximum, amplitude)."""
    val = np.array([r[key] for r in rec])
    i = int(np.argmax(val))
    i = min(max(i, 1), len(rec) - 2)
    a, b = sorted([rec[i - 1]["W"], rec[i + 1]["W"]])
    z0 = rec[i]["z"]
    pp = p.copy()

    def obj(Wv):
        nonlocal z0
        pp[I_W] = Wv
        z, ok, _ = newton_shoot(z0, pp, N)
        if ok:
            z0 = z
        ze, M, xmax, xmin, v2 = flow(z, pp, N)
        return Pe_from(v2, pp) if key == "Pe" else max(xmax, -xmin), max(xmax, -xmin)

    g = (np.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = obj(c)[0], obj(d)[0]
    for _ in range(40):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = obj(c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = obj(d)[0]
        if b - a < 1e-7:
            break
    Wp = 0.5 * (a + b)
    val, A = obj(Wp)
    return Wp, val, A


def refine_fold(rec, i):
    """Turning point near index i from a parabola W(A) through the neighbouring points. Returns (W_fold, A_fold)."""
    j0, j1 = max(i - 2, 0), min(i + 3, len(rec))
    A = np.array([rec[j]["A"] for j in range(j0, j1)])
    W = np.array([rec[j]["W"] for j in range(j0, j1)])
    c = np.polyfit(A, W, 2)
    Af = -c[1] / (2 * c[0])
    Wf = np.polyval(c, Af)
    return Wf, Af


# ---------------------------------------------------------------- fold curve in (Omega, f)
def _fold_residual(u, p, N):
    """Periodicity (3 equations) and det(M - I) = 0 for u = (z0, z1, z2, Omega, f)."""
    q = p.copy()
    q[I_W] = u[3]
    q[I_F] = u[4]
    ze, M, *_ = flow(u[:3], q, N)
    return np.r_[ze - u[:3], np.linalg.det(M - np.eye(3))]


def _fold_jacobian(u, p, N, h=1e-7):
    h0 = _fold_residual(u, p, N)
    Jm = np.zeros((4, 5))
    for k in range(5):
        e = np.zeros(5)
        e[k] = h * max(1, abs(u[k]))
        Jm[:, k] = (_fold_residual(u + e, p, N) - _fold_residual(u - e, p, N)) / (2 * e[k])
    return h0, Jm


def _null_tangent(Jm, tprev=None):
    t = np.linalg.svd(Jm)[2][-1]
    if tprev is not None and t @ tprev < 0:
        t = -t
    return t


def _fold_correct(u, p, N, t=None, upred=None, fixed_f=False):
    for it in range(25):
        h0, Jm = _fold_jacobian(u, p, N)
        if fixed_f:
            du = np.linalg.solve(Jm[:, :4], -h0)
            u = u + np.r_[du, 0]
        else:
            u = u - np.linalg.solve(np.vstack([Jm, t]), np.r_[h0, t @ (u - upred)])
        if np.max(np.abs(h0)) < 1e-10 and it > 0:
            return u, True
    return u, np.max(np.abs(h0)) < 1e-8


def fold_curve(u0, p, N=1024, ds=0.02, ds_max=0.05, maxpts=400, f_stop=None, verbose=True):
    """Two-parameter pseudo-arclength continuation of the fold (saddle-node) of periodic orbits.

    u0 = (z, Omega, f) near a fold; the curve is followed towards larger f, through the cusp,
    until f drops below f_stop. The Jacobian is formed by central differences.
    Returns the array of points (z0, z1, z2, Omega, f) along the curve."""
    u, ok = _fold_correct(np.asarray(u0, float), p, N, fixed_f=True)
    if verbose:
        print("fold at f=%.3f: Omega=%.7f converged=%s" % (u[4], u[3], ok), file=sys.stderr, flush=True)
    curve = [u.copy()]
    _, Jm = _fold_jacobian(u, p, N)
    t = _null_tangent(Jm)
    if t[4] < 0:
        t = -t
    for k in range(maxpts):
        upred = u + ds * t
        un, ok = _fold_correct(upred.copy(), p, N, t, upred)
        if not ok:
            ds *= 0.5
            if ds < 1e-5:
                break
            continue
        _, Jm = _fold_jacobian(un, p, N)
        tn = _null_tangent(Jm, t)
        u, t = un, tn
        curve.append(u.copy())
        ds = min(ds * 1.3, ds_max)
        if f_stop is not None and u[4] < f_stop and len(curve) > 10:
            break
        if verbose and k % 10 == 0:
            print(k, "Omega=%.5f f=%.5f" % (u[3], u[4]), file=sys.stderr, flush=True)
    return np.array(curve)


def cusp_of(curve):
    """Cusp = maximum of f along the fold curve, from a parabola in the arclength index. Returns (f_c, Omega_c)."""
    i = int(np.argmax(curve[:, 4]))
    j = slice(max(i - 2, 0), i + 3)
    s = np.arange(len(curve))[j]
    cf = np.polyfit(s, curve[j, 4], 2)
    sc = -cf[1] / (2 * cf[0])
    return np.polyval(cf, sc), np.interp(sc, s, curve[j, 3])
