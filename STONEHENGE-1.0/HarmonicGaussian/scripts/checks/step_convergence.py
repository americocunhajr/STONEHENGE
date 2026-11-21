"""Time-step convergence of capture fractions and power near the accessibility thresholds
(Supplementary Sect. S2.2): fixed-frequency ensembles with 128, 256 and 512 steps per period, and
sweeps with time steps 0.04, 0.02 and 0.01."""

import numpy as np
from harmonic_gaussian.model import params, A_WELL
from harmonic_gaussian.ensembles import run_fixed, run_ramp, rest_starts, interwell

x0, y0, phi = rest_starts(32, -A_WELL)
for f, W in [(0.1587, 0.90), (0.2, 0.98), (0.2, 1.0), (0.252, 1.0), (0.126, 0.86)]:
    p = params("gauss", f=f)
    row = []
    for spp in [128, 256, 512]:
        o = run_fixed(x0, y0, phi, p, W, 1000, 40, spp)
        row.append(((interwell(o[:, 0], o[:, 1]) & (o[:, 3] < 1e-3)).mean(), o[:, 2].mean()))
    print(
        f"fixed f={f} W={W}: "
        + "  ".join(f"spp={s}: cap={c:.3f} Pe={P:.3e}" for s, (c, P) in zip([128, 256, 512], row))
    )
for f, Tr in [(0.15, 2000), (0.15, 8000)]:
    p = params("gauss", f=f)
    row = []
    for h in [0.04, 0.02, 0.01]:
        o = run_ramp(x0, phi, p, 0.4, 1.0, float(Tr), 2000.0, h, 251.0)
        row.append(np.mean(np.maximum(o[:, 0], -o[:, 1]) > 3))
    print(f"ramp f={f} Tr={Tr}: " + "  ".join(f"h={h}: {c:.3f}" for h, c in zip([0.04, 0.02, 0.01], row)))
