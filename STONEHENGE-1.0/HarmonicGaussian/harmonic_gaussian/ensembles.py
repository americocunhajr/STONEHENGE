"""Ensembles of trajectories integrated in parallel (Numba ``prange``), used for everything that
concerns reaching the resonance from rest (Sect. 8): fixed-frequency runs, slow frequency sweeps,
first-capture times, stroboscopic (Poincare) sections, and maps of final states with early stopping.

All kernels use the same fixed-step RK4 scheme on (x, x', v) with forcing f cos(Omega t + phi);
a start is a point (x0, y0) with v = 0 and a forcing phase phi.
"""

import numpy as np
from numba import njit, prange
from .model import dU


@njit(cache=True)
def _acc(x, y, v, p, F):
    return -2 * p[0] * y - dU(x, p) + p[1] * v + F


@njit(cache=True)
def _step(x, y, v, p, F0, Fm, F1, h):
    """One RK4 step; F0, Fm, F1 are the forcing at the start, middle and end of the step."""
    lam, kap = p[3], p[2]
    k1x = y
    k1y = _acc(x, y, v, p, F0)
    k1v = -lam * v - kap * y
    x2 = x + 0.5 * h * k1x
    y2 = y + 0.5 * h * k1y
    v2 = v + 0.5 * h * k1v
    k2x = y2
    k2y = _acc(x2, y2, v2, p, Fm)
    k2v = -lam * v2 - kap * y2
    x3 = x + 0.5 * h * k2x
    y3 = y + 0.5 * h * k2y
    v3 = v + 0.5 * h * k2v
    k3x = y3
    k3y = _acc(x3, y3, v3, p, Fm)
    k3v = -lam * v3 - kap * y3
    x4 = x + h * k3x
    y4 = y + h * k3y
    v4 = v + h * k3v
    k4x = y4
    k4y = _acc(x4, y4, v4, p, F1)
    k4v = -lam * v4 - kap * y4
    return (
        x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x),
        y + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y),
        v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v),
    )


@njit(cache=True)
def _period(x, y, v, p, W, f, ph, h, spp, i0):
    """One forcing period of spp steps; the time of step i is (i0 + i) h.

    Returns the end state and, over the period: max |x|, max x, min x, the trapezoidal integral
    of v^2 and the plain sum of v^2 at the step ends."""
    amax = 0.0
    xmx = -1e9
    xmn = 1e9
    strap = 0.0
    sv2 = 0.0
    for i in range(spp):
        t = (i0 + i) * h
        vo = v
        x, y, v = _step(
            x, y, v, p, f * np.cos(W * t + ph), f * np.cos(W * (t + 0.5 * h) + ph), f * np.cos(W * (t + h) + ph), h
        )
        if abs(x) > amax:
            amax = abs(x)
        if x > xmx:
            xmx = x
        if x < xmn:
            xmn = x
        strap += 0.5 * (vo * vo + v * v) * h
        sv2 += v * v
    return x, y, v, amax, xmx, xmn, strap, sv2


@njit(parallel=True, cache=True)
def run_fixed(x0, y0, phi, p, W, nper, nwin, spp):
    """Integrate nper forcing periods at frequency W; statistics over the last nwin periods.

    Returns per start: max x, min x, mean power, standard deviation of the stroboscopic x."""
    n = x0.shape[0]
    out = np.zeros((n, 4))
    h = 2 * np.pi / W / spp
    f = p[4]
    cpe = p[1] * p[3] / p[2]
    for k in prange(n):
        x = x0[k]
        y = y0[k]
        v = 0.0
        ph = phi[k]
        xmax = -1e9
        xmin = 1e9
        acc = 0.0
        s1 = 0.0
        s2 = 0.0
        cnt = 0
        for j in range(nper):
            x, y, v, am, mx, mn, st, sv = _period(x, y, v, p, W, f, ph, h, spp, j * spp)
            if j >= nper - nwin:
                if mx > xmax:
                    xmax = mx
                if mn < xmin:
                    xmin = mn
                acc += st
                s1 += x
                s2 += x * x
                cnt += 1
        out[k, 0] = xmax
        out[k, 1] = xmin
        out[k, 2] = cpe * acc / (nwin * spp * h)
        m = s1 / cnt
        out[k, 3] = np.sqrt(max(s2 / cnt - m * m, 0.0))
    return out


@njit(parallel=True, cache=True)
def run_capture(x0, y0, phi, p, W, nmax, spp, Athr):
    """Integrate up to nmax periods; a start is captured when max|x| over a period exceeds Athr
    for 20 consecutive periods. Returns per start: capture period (first of the 20) or -1,
    standard deviation of the stroboscopic x over the last 40 periods, max x and min x over them."""
    n = x0.shape[0]
    out = np.zeros((n, 4))
    h = 2 * np.pi / W / spp
    f = p[4]
    for k in prange(n):
        x = x0[k]
        y = y0[k]
        v = 0.0
        ph = phi[k]
        cap = -1
        run = 0
        s1 = 0.0
        s2 = 0.0
        cnt = 0
        xmx = -1e9
        xmn = 1e9
        for j in range(nmax):
            x, y, v, amax, mx, mn, st, sv = _period(x, y, v, p, W, f, ph, h, spp, 0)
            if j >= nmax - 40:
                if mx > xmx:
                    xmx = mx
                if mn < xmn:
                    xmn = mn
                s1 += x
                s2 += x * x
                cnt += 1
            if amax > Athr:
                run += 1
                if run >= 20:
                    cap = j - 19
                    break
            else:
                run = 0
        out[k, 0] = cap
        if cnt > 0:
            m = s1 / cnt
            out[k, 1] = np.sqrt(max(s2 / cnt - m * m, 0.0))
        out[k, 2] = xmx
        out[k, 3] = xmn
    return out


@njit(parallel=True, cache=True)
def run_map(x0, y0, phi, W, F, p, nmax, spp, tol):
    """Final state of each start, with early stopping (maps of Fig. 8a-d).

    Frequencies W and forcing levels F are given per start. A start stops once its stroboscopic
    state (x, x') repeats with period k in {1, 2, 3} to within tol for 10 consecutive periods; power
    and extremes are then measured over k more periods. Otherwise it runs nmax periods and is
    measured over the last 40. Returns per start: [periodic (1/0), period k, mean v^2, max x,
    min x, periods used]."""
    n = x0.shape[0]
    out = np.zeros((n, 6))
    for k in prange(n):
        Wk = W[k]
        fk = F[k]
        ph = phi[k]
        h = 2 * np.pi / Wk / spp
        x = x0[k]
        y = y0[k]
        v = 0.0
        sx = np.zeros(4)
        sy = np.zeros(4)
        cnt = 0
        kper = 0
        done = False
        j = 0
        s1 = 0.0
        xmx = -1e9
        xmn = 1e9
        nw = 0
        while j < nmax and not done:
            x, y, v, am, mx, mn, st, sv = _period(x, y, v, p, Wk, fk, ph, h, spp, 0)
            if j >= nmax - 40:
                s1 += sv
                nw += spp
                if mx > xmx:
                    xmx = mx
                if mn < xmn:
                    xmn = mn
            j += 1
            for q in range(3, 0, -1):
                sx[q] = sx[q - 1]
                sy[q] = sy[q - 1]
            sx[0] = x
            sy[0] = y
            if j > 4:
                found = 0
                for q in range(1, 4):
                    if abs(sx[0] - sx[q]) < tol and abs(sy[0] - sy[q]) < tol:
                        found = q
                        break
                if found > 0 and found == kper:
                    cnt += 1
                elif found > 0:
                    kper = found
                    cnt = 1
                else:
                    kper = 0
                    cnt = 0
                if cnt >= 10:  # periodic: measure over kper periods
                    s1 = 0.0
                    nw = 0
                    xmx = -1e9
                    xmn = 1e9
                    for jj in range(kper):
                        x, y, v, am, mx, mn, st, sv = _period(x, y, v, p, Wk, fk, ph, h, spp, 0)
                        s1 += sv
                        nw += spp
                        if mx > xmx:
                            xmx = mx
                        if mn < xmn:
                            xmn = mn
                    done = True
        out[k, 0] = 1.0 if done else 0.0
        out[k, 1] = kper if done else 0
        out[k, 2] = s1 / max(nw, 1)
        out[k, 3] = xmx
        out[k, 4] = xmn
        out[k, 5] = j
    return out


@njit(parallel=True, cache=True)
def run_ramp(x0, phi, p, W0, W1, Tr, Thold, h, Twin):
    """Start at rest at x0; the frequency rises linearly from W0 to W1 over Tr and is then held for
    Thold. The forcing phase is theta = int W dt + phi. Returns max x, min x and mean power over the
    final Twin."""
    n = x0.shape[0]
    out = np.zeros((n, 3))
    f = p[4]
    cpe = p[1] * p[3] / p[2]
    N = int((Tr + Thold) / h)
    Nw = int(Twin / h)
    for k in prange(n):
        x = x0[k]
        y = 0.0
        v = 0.0
        th = phi[k]
        xmax = -1e9
        xmin = 1e9
        acc = 0.0
        for i in range(N):
            t = i * h
            Wa = W0 + (W1 - W0) * min(t / Tr, 1.0)
            Wm = W0 + (W1 - W0) * min((t + 0.5 * h) / Tr, 1.0)
            Wb = W0 + (W1 - W0) * min((t + h) / Tr, 1.0)
            thm = th + 0.5 * h * 0.5 * (Wa + Wm)
            thb = th + h / 6 * (Wa + 4 * Wm + Wb)
            vo = v
            x, y, v = _step(x, y, v, p, f * np.cos(th), f * np.cos(thm), f * np.cos(thb), h)
            th = thb
            if i >= N - Nw:
                if x > xmax:
                    xmax = x
                if x < xmin:
                    xmin = x
                acc += 0.5 * (vo * vo + v * v) * h
        out[k, 0] = xmax
        out[k, 1] = xmin
        out[k, 2] = cpe * acc / (Nw * h)
    return out


@njit(cache=True)
def strobe(x, y, v, p, W, ntr, nrec, spp, ph):
    """Stroboscopic samples (x, x') at forcing phase ph, after ntr transient periods."""
    h = 2 * np.pi / W / spp
    f = p[4]
    out = np.zeros((nrec, 2))
    for j in range(ntr + nrec):
        x, y, v, am, mx, mn, st, sv = _period(x, y, v, p, W, f, ph, h, spp, 0)
        if j >= ntr:
            out[j - ntr, 0] = x
            out[j - ntr, 1] = y
    return out


def rest_starts(nph, left, right=None):
    """Starts at rest at the bottom of each well with nph uniformly spaced forcing phases per well.

    Returns (x0, y0, phi) with the nph left-well starts first; right defaults to -left."""
    right = -left if right is None else right
    phi = np.tile(np.linspace(0, 2 * np.pi, nph, endpoint=False), 2)
    x0 = np.r_[np.full(nph, left), np.full(nph, right)]
    return x0, np.zeros(2 * nph), phi


def interwell(xmax, xmin, gap=0.3):
    """True for motion that visits both wells (max x > gap and min x < -gap)."""
    return (xmax > gap) & (xmin < -gap)
