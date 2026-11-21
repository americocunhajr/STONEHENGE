"""Figure 9: capture from rest as a boundary crisis (baseline Gaussian, Omega = 1).
(a) below the threshold: the chaotic attractor (Poincare section) over the basin of the resonant
orbit; (b) above it: capture time over the plane of initial states (data/basins.pkl);
(c) survival of the chaotic transient against t/tau; (d) mean capture time against the distance
to the crisis, with the power-law fit (data/crisis_fine.pkl)."""

import numpy as np
from matplotlib.colors import LogNorm
from harmonic_gaussian.analysis import crisis_lifetimes, fit_crisis
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, TXT, RAMP

setup()
B = load("basins")
xs, ys = B["xs"], B["ys"]
keys = sorted([k[1] for k in B if isinstance(k, tuple) and k[0] == "basin"])
fb, fa = keys[0], keys[-1]
fig, axs = plt.subplots(2, 2, figsize=(W2, 5.4), constrained_layout=True)
a, b, c, d = axs.ravel()

# (a) below threshold: chaotic attractor over the basin of the resonant orbit
a.contourf(xs, ys, (B[("basin", fb)] >= 0).astype(float), levels=[0.5, 1.5], colors=[RAMP[0]], alpha=0.55)
S = B[("strobe", fb)]
a.plot(S[:, 0], S[:, 1], ",", color=TXT, alpha=0.35, rasterized=True)
a.set_xlim(xs[0], xs[-1])
a.set_ylim(ys[0], ys[-1])
a.set_xlabel("$x$")
a.set_ylabel("$\\dot x$")
a.set_title(f"$f={fb}$", fontsize=8)
panel(a, "(a)")

# (b) above threshold: capture time
Ta = B[("basin", fa)].astype(float)
Ta[Ta < 0] = np.nan
im = b.pcolormesh(
    xs, ys, Ta + 1, norm=LogNorm(vmin=1, vmax=np.nanmax(Ta) + 1), cmap="Blues_r", shading="auto", rasterized=True
)
cb = fig.colorbar(im, ax=b, shrink=0.85, pad=0.02)
cb.set_label("capture time (periods)", fontsize=7)
cb.ax.tick_params(labelsize=6)
b.set_xlabel("$x$")
b.set_ylabel("$\\dot x$")
b.set_xticks([-2, 0, 2])
a.set_xticks([-2, 0, 2])
b.set_title(f"$f={fa}$", fontsize=8)
panel(b, "(b)")

# (c) survival collapse and (d) lifetime against distance to the crisis
out = crisis_lifetimes(load("crisis_fine")["res"])
fv = np.array([o[0] for o in out])
tv = np.array([o[1] for o in out])
p, pe = fit_crisis(fv, tv)
for (f, tau, r), col in zip(out[:6], RAMP + ["#0b0b0b", "#52514e"]):
    T = np.sort(r["T"][r["T"] >= 0])
    n = len(r["T"]) if f < 0.218 else len(T)
    St = 1 - np.arange(1, len(T) + 1) / n
    c.semilogy(T / tau, np.clip(St, 1e-3, 1), "-", color=col, lw=1.0, label=f"{f:g}")
tt = np.linspace(0, 5, 50)
c.semilogy(tt, np.exp(-tt), "--", color=TXT, lw=0.8)
c.set_xlim(0, 5)
c.set_ylim(5e-3, 1.1)
c.set_xlabel("$t/\\tau$")
c.set_ylabel("not yet captured")
c.legend(fontsize=5.6, title="$f$", title_fontsize=6, loc="lower left", handlelength=1.2, labelspacing=0.2)
panel(c, "(c)")
grid(c)
d.loglog(fv - p[1], tv, "o", color=BLUE, ms=3.5)
xx = np.geomspace(min(fv - p[1]) * 0.8, max(fv - p[1]) * 1.2, 50)
d.loglog(xx, np.exp(p[0]) * xx ** (-p[2]), "--", color=TXT, lw=0.8)
d.set_xlabel("$f-f_{\\rm cr}$")
d.set_ylabel("mean capture time $\\tau$ (periods)")
panel(d, "(d)")
grid(d)
d.text(0.05, 0.08, f"slope $-{p[2]:.2f}$", transform=d.transAxes, fontsize=7)
save(fig, "fig09_mechanism", dpi=400)
