"""(i) Fold at f = 0.15, xi = 0.01 with maximum continuation step 0.2 (production) against 0.1.
(ii) Existence of the Gaussian fold just below and above the cusp (Sect. 5.4)."""

import numpy as np
from harmonic_gaussian.model import params, ctot
from harmonic_gaussian.orbits import march, continuation, analyse, refine_fold
from harmonic_gaussian.branches import solve_case

f = 0.15
p = params("gauss", xi=0.01, lam=0.05, f=f, W=1.0)
z = march(np.array([0.0, f / ctot(0.01, 0.05), 0.0]), p, 512, 400)
for dsm in [0.2, 0.1]:
    up = continuation(
        z, p, N=512, direction=+1, Wmax=2.0, Wmin=0.3, ds=0.05, ds_max=dsm, maxpts=8000, stop_after_fold=10
    )
    W, A, Pe, mu, un, folds = analyse(up)
    print(f"f=0.15 fold with ds_max={dsm}: {refine_fold(up, folds[0])[0]:.7f}")
for f in [1.55, 1.6, 1.65, 1.70, 1.75]:
    o = solve_case(f=f, Wmax=2.5, down=False)
    print(f'f={f}: ok={o.get("ok")} fold={o.get("W_fold", np.nan):.4f} peak={o.get("W_peak", np.nan):.5f}')
