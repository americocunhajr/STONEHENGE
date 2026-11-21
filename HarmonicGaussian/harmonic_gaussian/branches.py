"""Interwell resonance branches: harmonic-balance initial guesses and the branch solver that
returns the peak of power, the fold (jump) frequency and the low end of the stable band.

``solve_case`` is the workhorse behind the scaling, shape, hardening and similarity data:
it finds the large interwell orbit, continues it upwards past the fold and downwards to the
first loss of stability, and refines the peak of power and the fold.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from .model import params, dU, ctot, linear_state
from .orbits import flow, march, newton_shoot, continuation, analyse, refine_peak, refine_fold, multiplier_type, Pe_from

LN2 = np.log(2)
_wL = (1 / np.sqrt(np.pi)) ** (1 / 3)
_wD = (3 * np.sqrt(np.pi) / 4) ** (1 / 3)
#: barrier shapes with the same area M0 = sqrt(pi) and the same curvature at the top (Sect. 5.4)
SHAPES = {"gauss": (1, 1, 1, 0), "lorentz": (_wL, 1, _wL**2, 0), "dipole": (_wD, 1, 2 * _wD**2 / 3, 0)}
#: quartic models (k1, k3) of Sect. 7: nominal, and calibrated to (a, h), to (a, k_w) and by Taylor expansion
DUFF = {
    "nominal": (0.5, 0.5),
    "(a,h)": (4 * 0.5 * (1 - LN2) / LN2, 4 * 0.5 * (1 - LN2) / LN2**2),
    "(a,k_w)": (LN2, 1.0),
    "Taylor": (1.0, 2.0),
}


# ---------------------------------------------------------------- harmonic balance
def keff_num(p, A):
    """Single-harmonic effective stiffness (1/(pi A)) int U'(A cos th) cos th dth."""
    I = quad(lambda th: dU(A * np.cos(th), p) * np.cos(th), 0, 2 * np.pi, limit=400)[0]
    return I / (np.pi * A)


def hb_state(p, W0, Amin=1.3, Amax=400):
    """Largest single-harmonic interwell solution at frequency W0 (with circuit terms).
    Returns (state (x, x', v) at phase 0, amplitude), or None."""
    xi, chi, kap, lam, f = p[0], p[1], p[2], p[3], p[4]
    ke = chi * kap * W0**2 / (lam**2 + W0**2)
    c = 2 * xi + chi * kap * lam / (lam**2 + W0**2)
    G = lambda A: (keff_num(p, A) + ke - W0**2) ** 2 + (c * W0) ** 2 - (f / A) ** 2
    As = np.geomspace(Amin, Amax, 400)
    g = np.array([G(a) for a in As])
    idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    if len(idx) == 0:
        return None
    A = brentq(G, As[idx[-1]], As[idx[-1] + 1])
    phi = np.arctan2(c * W0, keff_num(p, A) + ke - W0**2)
    X = A * np.exp(-1j * phi)
    V = -1j * kap * W0 * X / (lam + 1j * W0)
    return np.array([A * np.cos(phi), A * W0 * np.sin(phi), V.real]), A


# ---------------------------------------------------------------- branch solver
def run_case(
    shape="gauss",
    s=(1, 1, 1, 0),
    f=0.15,
    xi=0.01,
    lam=0.05,
    delta=0.0,
    W0=1.0,
    A0=None,
    Wmax=1.6,
    Wmin=0.3,
    thr=1.4,
    N=512,
    down=True,
    z_init=None,
    k3h=0.0,
):
    """Continue the interwell branch through (W0, z_init or the linear amplitude f/(c W0)).

    Returns a dict with W_fold, A_fold, W_peak, Pe_peak, A_peak, the low end of the stable band
    (W_low, its multiplier type), the largest Floquet modulus at the peak, and the branch itself
    as a list of (W, A, Pe, max|mu|). ``ok`` is False when no large interwell orbit is found."""
    p = params(shape, xi=xi, lam=lam, f=f, W=W0, s=s, delta=delta, k3h=k3h)
    c = ctot(xi, lam, W0)
    if A0 is None:
        A0 = f / (c * W0)
    if z_init is None:
        z = march(np.array([0.0, A0 * W0, 0.0]), p, N, 400)
    else:
        z, _, _ = newton_shoot(z_init, p, N)
    ze, M, xmax, xmin, v2 = flow(z, p, N)
    if max(xmax, -xmin) < thr or min(xmax, -xmin) < 0.2:  # not a large interwell orbit
        return dict(f=f, xi=xi, lam=lam, ok=False)
    up = continuation(
        z, p, N=N, direction=+1, Wmax=Wmax, Wmin=Wmin, ds=0.05, ds_max=0.2, maxpts=4000, stop_after_fold=10
    )
    W, A, Pe, mu, un, folds = analyse(up)
    out = dict(f=f, xi=xi, lam=lam, delta=delta, ok=True)
    if folds:
        Wf, Af = refine_fold(up, folds[0])
        out.update(W_fold=Wf, A_fold=Af, Pe_fold=up[folds[0]]["Pe"])
    else:
        out.update(W_fold=np.nan, A_fold=np.nan, Pe_fold=np.nan)
    pk = up[: folds[0] + 1] if folds else up
    if down:
        dn = continuation(z, p, N=N, direction=-1, Wmax=Wmax, Wmin=Wmin, ds=0.2, ds_max=0.8, maxpts=800)
        Wd, Ad, Ped, mud, und, fd = analyse(dn)
        k = np.where(und)[0]  # first instability going down
        if len(k):
            i = k[0]
            out.update(
                W_low=0.5 * (Wd[i] + Wd[i - 1]),
                A_low=Ad[i],
                low_type=multiplier_type(dn[i]["mu"]),
                low_fold=bool(fd and fd[0] <= i),
            )
        else:
            out.update(W_low=np.nan, A_low=np.nan, low_type="none", low_fold=False)
        allrec = dn[::-1] + pk
    else:
        allrec = pk
    Wp, Pp, Ap = refine_peak(allrec, p, N)
    out.update(W_peak=Wp, Pe_peak=Pp, A_peak=Ap)
    pp = p.copy()
    pp[5] = Wp  # Floquet modulus of the orbit at the peak
    zz, ok, _ = newton_shoot(allrec[int(np.argmax([r["Pe"] for r in allrec]))]["z"], pp, N)
    ze, M, *_ = flow(zz, pp, N)
    out["peak_maxmu"] = np.max(np.abs(np.linalg.eigvals(M)))
    out["branch"] = [(r["W"], r["A"], r["Pe"], float(np.max(np.abs(r["mu"])))) for r in allrec]
    return out


def solve_case(shape="gauss", s=(1, 1, 1, 0), f=0.15, xi=0.01, lam=0.05, delta=0.0, k3h=0.0, Wmax=2.0, down=True):
    """``run_case`` from the linear amplitude at Omega = 1, with harmonic-balance restarts if that fails."""
    try:
        o = run_case(shape=shape, s=s, f=f, xi=xi, lam=lam, delta=delta, k3h=k3h, Wmax=Wmax, down=down)
        if o.get("ok"):
            return o
    except Exception:
        pass
    p = params(shape, xi=xi, lam=lam, f=f, W=1.0, s=s, delta=delta, k3h=k3h)
    for W0 in [1.0, 0.99, 0.98, 0.95, 1.005, 0.9]:
        p[5] = W0
        r = hb_state(p, W0, Amin=1.3 if shape != "duffing" else 0.3, Amax=300)
        if r is None:
            continue
        z, ok, _ = newton_shoot(r[0], p, 512)
        if not ok:
            continue
        try:
            o = run_case(
                shape=shape,
                s=s,
                f=f,
                xi=xi,
                lam=lam,
                delta=delta,
                k3h=k3h,
                W0=W0,
                z_init=z,
                Wmax=Wmax,
                thr=0.0,
                down=down,
            )
            if o.get("ok"):
                return o
        except Exception:
            continue
    return dict(f=f, xi=xi, lam=lam, ok=False)


def solve_duffing(k1, k3, f, xi=0.01, lam=0.05, Wmax=16, Wmin=0.3, down=False):
    """Large-orbit branch of the quartic model U' = -k1 x + k3 x^3, started from harmonic balance
    below its jump frequency."""
    c = ctot(xi, lam)
    WJ = max((3 * k3 / 4) ** 0.25 * (f / c) ** 0.5, 1.0)
    p = params("duffing", xi=xi, lam=lam, f=f, W=1.0, s=(k1, k3, 0, 0))
    for fr in [0.6, 0.75, 0.5, 0.9]:
        W0 = max(fr * WJ, 0.9)
        p[5] = W0
        r = hb_state(p, W0, Amin=0.3, Amax=300)
        if r is None:
            continue
        z, ok, _ = newton_shoot(r[0], p, 512)
        if not ok:
            continue
        try:
            o = run_case(
                shape="duffing",
                s=(k1, k3, 0, 0),
                f=f,
                xi=xi,
                lam=lam,
                W0=W0,
                z_init=z,
                Wmax=Wmax,
                Wmin=Wmin,
                thr=0.0,
                down=down,
            )
            if o.get("ok"):
                return o
        except Exception:
            continue
    return dict(f=f, xi=xi, ok=False)


def pack(rec):
    """Continuation record as an array of rows (W, A, Pe, max|mu|, xmax, xmin)."""
    return np.array([(r["W"], r["A"], r["Pe"], float(np.max(np.abs(r["mu"]))), r["xmax"], r["xmin"]) for r in rec])


def refine_peak_direct(o, xi=0.01, lam=0.05, s=(1, 1, 1, 0), Ns=(1024, 2048), halfwidth=2e-4):
    """Peak of power of a very large orbit, re-maximized over Omega by Newton shooting from the
    linear response with N steps per period (continuation with 512 steps is not accurate enough
    when A_peak sqrt(alpha) > 60). Returns [(Omega_peak, Pe_peak)] for each N in Ns."""
    res = []
    for N in Ns:

        def mP(W):
            p = params("gauss", xi=xi, lam=lam, f=o["f"], W=W, s=s)
            z, ok, _ = newton_shoot(linear_state(o["f"], W, xi, lam), p, N)
            return -Pe_from(flow(z, p, N)[4], p)

        r = minimize_scalar(
            mP, bounds=(o["W_peak"] - halfwidth, o["W_peak"] + halfwidth), method="bounded", options={"xatol": 1e-10}
        )
        res.append((r.x, -r.fun))
    return res
