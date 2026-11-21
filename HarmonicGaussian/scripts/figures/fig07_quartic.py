"""Figure 7: what quartic models predict. (a) Frequency of maximum power against forcing for the
Gaussian and three calibrated quartics, and for quartics fitted over the predicted stroke;
(b) displacement harmonics of the Gaussian orbits against the closed form, and of a quartic orbit
(data/scaling.pkl, data/quartic_fits.pkl, data/harmonics.pkl)."""

import numpy as np
from matplotlib.ticker import NullFormatter
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, ORANGE, RAMP

setup()
S = load("scaling")
Fi = load("quartic_fits")
Fp = load("harmonics")
OR = ["#f09a78", "#eb6834", "#81391d"]
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.5), constrained_layout=True)
a, b = axs
G = [o for o in S["gauss"][0.01] if o.get("ok")]
a.plot([o["f"] for o in G], [o["W_peak"] for o in G], "-", color=BLUE, lw=1.5, label="Gaussian peak")
a.plot(
    [o["f"] for o in G if np.isfinite(o["W_fold"])],
    [o["W_fold"] for o in G if np.isfinite(o["W_fold"])],
    "--",
    color=BLUE,
    lw=1.2,
    label="Gaussian fold",
)
for i, (name, lab) in enumerate(
    [("(a,h)", "quartic, $(a,h)$"), ("(a,k_w)", "quartic, $(a,k_w)$"), ("Taylor", "quartic, Taylor")]
):
    d = [o for o in S["calib"][name] if o.get("ok")]
    a.plot([o["f"] for o in d], [o["W_fold"] for o in d], "-", color=OR[i], lw=1.2, label=lab)
sc = Fi["sc"]
a.plot(
    [q["f"] for q in sc],
    [q["W"] for q in sc],
    "x",
    color=ORANGE,
    ms=4.5,
    mew=1.0,
    label="quartic fitted over amplitude",
)
a.set_xscale("log")
a.set_yscale("log")
a.set_xticks([0.05, 0.1, 0.2, 0.5, 1])
a.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1"])
a.set_yticks([1, 2, 4, 8])
a.set_yticklabels(["1", "2", "4", "8"])
a.yaxis.set_minor_formatter(NullFormatter())
a.xaxis.set_minor_formatter(NullFormatter())
a.legend(loc="upper left", fontsize=6.6, labelspacing=0.25)
a.set_xlabel("forcing $f$")
a.set_ylabel("frequency of maximum power")
panel(a, "(a)")
grid(a)
for i, f in enumerate([0.15, 0.3, 0.6]):
    q = Fp[("gauss", f)]
    A = q["A"]
    b.semilogy(q["ks"], q["X"] / A, "o", color=RAMP[i + 1], ms=3.2, ls="none", label=f"Gaussian, $f={f}$ ($A={A:.0f}$)")
    b.semilogy(q["ks"], q["pred"] / A, "-", color=RAMP[i + 1], lw=0.9)
q = Fp[("duff", 0.3)]
m = q["X"] / q["A"] > 1e-10
b.semilogy(q["ks"][m], (q["X"] / q["A"])[m], "s-", color=ORANGE, ms=3.2, lw=0.9, label="quartic $(a,k_w)$, $f=0.3$")
b.set_ylim(1e-12, 1e-1)
b.set_xlim(1, 42)
b.legend(loc="upper right", fontsize=6.6, labelspacing=0.25)
b.set_xlabel("harmonic order $k$")
b.set_ylabel("displacement harmonic $|X_k|/|X_1|$")
panel(b, "(b)")
grid(b)
save(fig, "fig07_quartic")
