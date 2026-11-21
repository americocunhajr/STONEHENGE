"""Figure 6: one similarity parameter controls the resonance. (a) Peak gap and (b) fold distance
against eps for all Gaussian cases and the magnet pair; (c) half-power bandwidth against f/c for the
Gaussian and the quartic models (data/similarity.pkl, data/magnet.pkl)."""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from harmonic_gaussian.model import Wlin, ctot, beff
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, ORANGE, YELLOW, TXT, TXT2, RAMP, ORAMP, MK

setup()
E = load("similarity")
rows, duff = E["rows"], E["duff"]
OTH = {(1.0, 0.75): "v", (1.0, 2.0): ">", (0.5, 2.0): "<", (2.0, 0.5): "p"}


def style(q):
    if q["a"] == 1 and q["g"] == 1 and q["lam"] == 0.05:
        i = [0.005, 0.01, 0.02, 0.04].index(q["xi"])
        return dict(color=RAMP[i], marker=MK[q["xi"]], ms=3.6, mfc=RAMP[i], ls="none")
    if q["lam"] == 1.0:
        return dict(color=TXT2, marker="+", ms=5, ls="none", mew=0.9)
    return dict(color=TXT2, marker=OTH[(q["a"], q["g"])], ms=3.8, mfc="none", ls="none", mew=0.8)


fig, axs = plt.subplots(1, 3, figsize=(W2, 2.45), constrained_layout=True)
a, b, c = axs
for q in rows:
    if q["gap"] > 0:
        a.plot(q["eps"], q["gap"], **style(q))
    if np.isfinite(q["fold"]) and q["fold"] > 0:
        b.plot(q["eps"], q["fold"], **style(q))
Wl = Wlin(0.01, 0.05)
cc = ctot(0.01, 0.05)
be = beff(0.05)
for o in load("magnet")["forced"]:
    if not o.get("ok"):
        continue
    gp = (Wl - o["W_peak"]) / (cc / 2)
    fd = (o["W_fold"] - np.sqrt(be)) / (cc / 2)
    if gp > 0:
        a.plot(o["eps"], gp, "*", color=YELLOW, ms=6, mec=TXT, mew=0.4)
    if np.isfinite(fd) and fd > 0:
        b.plot(o["eps"], fd, "*", color=YELLOW, ms=6, mec=TXT, mew=0.4)
e = np.geomspace(3e-5, 10, 100)
a.plot(e, e, "--", color=TXT, lw=0.9)
a.set_xscale("log")
a.set_yscale("log")
a.set_ylim(2e-7, 40)
a.set_xlabel("similarity parameter $\\varepsilon$")
a.set_ylabel("peak gap $(\\Omega_{\\rm lin}-\\Omega_{\\rm peak})/(c/2)$")
a.text(2e-4, 3e-2, "slope 1", fontsize=7.5, color=TXT)
panel(a, "(a)")
grid(a)
b.plot(e, 2 / (3 * np.sqrt(3 * e)), "--", color=TXT, lw=0.9)
b.set_xscale("log")
b.set_yscale("log")
b.text(1.2e-2, 18, "$\\frac{2}{3\\sqrt{3\\varepsilon}}$", fontsize=9, color=TXT)
b.set_xlabel("similarity parameter $\\varepsilon$")
b.set_ylabel("fold distance $(\\Omega_J-\\sqrt{\\beta_{\\rm eff}})/(c/2)$")
panel(b, "(b)")
grid(b)
hand = [
    Line2D([], [], color=RAMP[i], marker=MK[x], ls="none", ms=3.6, label=f"$\\xi={x}$")
    for i, x in enumerate([0.005, 0.01, 0.02, 0.04])
]
hand += [Line2D([], [], color=TXT2, marker="+", ls="none", ms=5, label="$\\lambda=1$")]
hand += [
    Line2D(
        [], [], color=TXT2, marker=m, mfc="none", ls="none", ms=3.8, label=f"$(\\alpha,\\gamma)=({k[0]:g},{k[1]:g})$"
    )
    for k, m in OTH.items()
]
hand += [Line2D([], [], color=YELLOW, marker="*", mec=TXT, mew=0.4, ls="none", ms=6, label="magnet pair")]
a.legend(
    handles=hand,
    loc="lower right",
    ncol=2,
    columnspacing=0.5,
    fontsize=5.6,
    handletextpad=0.1,
    labelspacing=0.25,
    borderaxespad=0.2,
)
for q in rows:
    c.plot(q["fc"], q["bw"], **{**style(q), "ms": 3.0})
for i, name in enumerate(["nominal", "(a,k_w)", "(a,h)", "Taylor"]):
    d = sorted([q for q in duff if q["name"] == name], key=lambda q: q["fc"])
    c.plot([q["fc"] for q in d], [q["bw"] for q in d], "-", color=ORAMP[i], lw=1.2)
c.axhline(1, color=TXT, lw=0.8, ls="--")
c.set_xscale("log")
c.set_yscale("log")
c.set_xlim(1.5, 120)
c.set_ylim(0.4, 200)
c.text(2, 140, "quartic models", color=ORANGE, fontsize=7.5)
c.text(9, 0.55, "Gaussian: $\\Delta\\Omega_{1/2}\\approx c$", color=BLUE, fontsize=7.5)
c.set_xlabel("forcing over damping, $f/c$")
c.set_ylabel("half-power bandwidth $\\Delta\\Omega_{1/2}/c$")
panel(c, "(c)")
grid(c)
for ax in axs:
    ax.xaxis.set_minor_formatter(NullFormatter())
save(fig, "fig06_scaling")
