"""Expected power from rest (64 starts: uniform phases, both wells) on a frequency grid refined
near resonance, at six forcing levels (Sect. 9.2; Fig. 10b, c). Time step at most 0.025,
1000 periods, power averaged over the last 40.

Output: data/expected_power.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, A_WELL
from harmonic_gaussian.ensembles import run_fixed, rest_starts, interwell
from harmonic_gaussian.storage import main, pick, log


def compute():
    Ws = pick(
        np.unique(
            np.round(
                np.r_[np.arange(0.40, 0.95, 0.02), np.arange(0.95, 1.10001, 0.0025), np.arange(1.12, 1.401, 0.02)], 5
            )
        ),
        np.array([0.6, 0.9, 1.0, 1.2]),
    )
    fs = pick([0.10, 0.15, 0.20, 0.25, 0.30, 0.40], [0.15, 0.30])
    x0, y0, phi = rest_starts(pick(32, 8), -A_WELL)
    PE = np.zeros((len(fs), len(Ws)))
    AMP = np.zeros_like(PE)
    CAP = np.zeros_like(PE)
    for i, f in enumerate(fs):
        p = params("gauss", f=f)
        for j, W in enumerate(Ws):
            o = run_fixed(x0, y0, phi, p, W, pick(1000, 200), 40, int(np.ceil(2 * np.pi / W / 0.025)))
            PE[i, j] = o[:, 2].mean()
            AMP[i, j] = np.maximum(o[:, 0], -o[:, 1]).mean()
            CAP[i, j] = np.mean(interwell(o[:, 0], o[:, 1]) & (o[:, 3] < 1e-3))
        log("f", f, "done")
    return dict(Ws=Ws, fs=np.array(fs), PE=PE, AMP=AMP, CAP=CAP)


def report(P):
    print(f"{len(P['Ws'])} frequencies from {P['Ws'][0]} to {P['Ws'][-1]}")
    for i, f in enumerate(P["fs"]):
        j = int(np.argmax(P["PE"][i]))
        print(
            f"f={f}: max expected power {P['PE'][i, j]:.4e} at Omega={P['Ws'][j]:.4f}; fraction on periodic interwell orbits there {P['CAP'][i, j]:.2f}"
        )


if __name__ == "__main__":
    main("expected_power", compute, report)
