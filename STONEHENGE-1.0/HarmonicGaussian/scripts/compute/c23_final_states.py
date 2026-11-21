"""Map of final states from rest (Fig. 8a-d): 181 frequencies (0.4 to 1.3, step 0.005) x 121 forcing
levels (0.01 to 1, log-spaced), 64 starts per cell (32 phases in each well), up to 1000 periods with
early stopping once a start is periodic (period 1, 2 or 3).

K[i, j, :] = probabilities of period-1 interwell, period-2/3 interwell, chaotic interwell and
intrawell motion; mean power; mean amplitude; median settling time of the captured starts.
Output: data/final_states.pkl. Runtime: about 50 min on two cores.
"""

import warnings
import numpy as np
from harmonic_gaussian.model import params, A_WELL
from harmonic_gaussian.ensembles import run_map, rest_starts, interwell
from harmonic_gaussian.storage import main, pick, log


def compute():
    Ws = pick(np.round(np.arange(0.40, 1.30001, 0.005), 4), np.round(np.arange(0.40, 1.30001, 0.1), 4))
    fs = pick(np.geomspace(0.01, 1.0, 121), np.geomspace(0.05, 0.5, 4))
    nph = pick(32, 4)
    p = params("gauss")
    xs0, _, ph0 = rest_starts(nph, -A_WELL)
    ns = 2 * nph
    W = np.repeat(Ws, ns)
    phi = np.tile(ph0, len(Ws))
    x0 = np.tile(xs0, len(Ws))
    y0 = np.zeros(x0.size)
    K = np.full((len(fs), len(Ws), 7), np.nan)
    for i, f in enumerate(fs):
        o = run_map(x0, y0, phi, W, np.full(x0.size, f), p, 1000, 128, 1e-5).reshape(len(Ws), ns, 6)
        conv = o[:, :, 0] == 1
        per = o[:, :, 1]
        inter = interwell(o[:, :, 3], o[:, :, 4])
        hi = conv & inter & (per == 1)
        K[i, :, 0] = np.mean(hi, axis=1)
        K[i, :, 1] = np.mean(conv & inter & (per > 1), axis=1)
        K[i, :, 2] = np.mean(~conv & inter, axis=1)
        K[i, :, 3] = np.mean(~inter, axis=1)
        K[i, :, 4] = np.mean(o[:, :, 2], axis=1) * p[1] * p[3] / p[2]
        K[i, :, 5] = np.mean(np.maximum(o[:, :, 3], -o[:, :, 4]), axis=1)
        with warnings.catch_warnings():  # cells without captures give all-NaN medians
            warnings.simplefilter("ignore", RuntimeWarning)
            K[i, :, 6] = np.where(hi.sum(1) >= 8, np.nanmedian(np.where(hi, o[:, :, 5], np.nan), axis=1), np.nan)
        log(f"row {i+1}/{len(fs)} f={f:.4f}")
    return dict(Ws=Ws, fs=fs, K=K, nstarts=ns, horizon=1000, tol=1e-5)


def report(X):
    print(
        f"map of final states: {len(X['Ws'])} frequencies x {len(X['fs'])} forcing levels, {X['nstarts']} starts per cell, "
        f"horizon {X['horizon']} periods; summary in results/final_states_summary.txt"
    )


if __name__ == "__main__":
    main("final_states", compute, report)
