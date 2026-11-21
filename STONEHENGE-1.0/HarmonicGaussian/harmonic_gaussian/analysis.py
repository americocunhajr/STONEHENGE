"""Post-processing shared by a figure and the script that prints its numbers: the lifetime of the
chaotic transient near the boundary crisis (Fig. 9), capture thresholds and band-averaged power
ratios (Fig. 10), and the decomposition of the map of final states (Fig. 8)."""

import numpy as np
from scipy.optimize import curve_fit
from .model import P_lin


# ---------------------------------------------------------------- boundary crisis (Fig. 9)
def crisis_lifetimes(res, horizon=40000, f_censor=0.218, min_captured=10):
    """Mean capture time tau at each forcing level, as (f, tau, record).

    Near the crisis (f < f_censor) the uncaptured starts are chaotic transients censored at the
    horizon, and tau is the maximum-likelihood estimate of an exponential law with censoring;
    further up, the uncaptured starts sit on a coexisting period-3 orbit and tau is the plain mean."""
    out = []
    for r in sorted(res, key=lambda r: r["f"]):
        T = r["T"]
        capd = T >= 0
        if capd.sum() < min_captured:
            continue
        tau = (T[capd].sum() + (~capd).sum() * horizon) / capd.sum() if r["f"] < f_censor else T[capd].mean()
        out.append((r["f"], tau, r))
    return out


def fit_crisis(fv, tv, fmax_fit=0.26):
    """Fit tau = C (f - f_cr)^(-gamma). Returns parameters (ln C, f_cr, gamma) and their standard errors."""
    fun = lambda f, lC, fc, g: lC - g * np.log(np.clip(f - fc, 1e-9, None))
    m = fv < fmax_fit
    p, cov = curve_fit(
        fun, fv[m], np.log(tv[m]), p0=(0, 0.2171, 1.5), bounds=([-30, 0.2165, 0.1], [30, min(fv) - 1e-6, 10])
    )
    return p, np.sqrt(np.diag(cov))


# ---------------------------------------------------------------- capture thresholds and band power (Fig. 10)
def capture_thresholds(H):
    """For each damping level of the capture curves: (c, last f below 50% capture, first f at or above)."""
    pts = []
    for xi, d in H.items():
        rows = d["rows"]
        i = next(k for k, r in enumerate(rows) if r[1] >= 0.5)
        pts.append((d["c"], rows[i - 1][0], rows[i][0]))
    return sorted(pts)


def band_ratio(Ws, PE, f, W0, W1, xi=0.01):
    """Expected power from rest integrated over [W0, W1], divided by that of the linear harvester."""
    m = (Ws >= W0 - 1e-9) & (Ws <= W1 + 1e-9)
    G = np.trapezoid(PE[m], Ws[m])
    Wl = np.linspace(W0, W1, 40001)
    L = np.trapezoid(P_lin(Wl, f, xi=xi), Wl)
    return G / L


# ---------------------------------------------------------------- map of final states (Fig. 8)
def map_fields(X, xi=0.01):
    """Fields of the map of final states with all rows computed: frequencies, forcing levels,
    probabilities (period-1 interwell, period-2/3 interwell, chaotic interwell, intrawell),
    expected power, mean amplitude, median settling time, and the gain over the linear harvester."""
    Ws, fs, K = X["Ws"], X["fs"], X["K"]
    done = np.isfinite(K[:, :, 0]).all(axis=1)
    fs, K = fs[done], K[done]
    P1, PS, PC, PI, PE, AMP, TS = [K[:, :, k] for k in range(7)]
    GAIN = PE / P_lin(Ws[None, :], fs[:, None], xi=xi)
    return dict(Ws=Ws, fs=fs, P1=P1, PS=PS, PC=PC, PI=PI, PE=PE, AMP=AMP, TS=TS, GAIN=GAIN)
