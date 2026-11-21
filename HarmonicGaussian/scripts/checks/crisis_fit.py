"""Boundary crisis at Omega = 1 (Sect. 8.3, Fig. 9d): mean capture time against forcing and the fit
tau ~ (f - f_cr)^(-gamma); basins of Fig. 9a, b."""

import numpy as np
from harmonic_gaussian.analysis import crisis_lifetimes, fit_crisis
from harmonic_gaussian.storage import load

out = crisis_lifetimes(load("crisis_fine")["res"])
fv = np.array([o[0] for o in out])
tv = np.array([o[1] for o in out])
p, pe = fit_crisis(fv, tv)
print("fit: f_cr=%.5f +- %.5f, gamma=%.3f +- %.3f" % (p[1], pe[1], p[2], pe[2]))
for f, tau, r in out:
    print(f, round(tau, 1), "captured", np.mean(r["T"] >= 0))
B = load("basins")
for f in (0.21, 0.225):
    print(f"basin f={f}: captured {np.mean(B[('basin', f)] >= 0):.4f}")
