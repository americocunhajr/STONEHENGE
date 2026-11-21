"""Figure 3: (a) the static invariant Q = h/(k_w a^2) for three barrier shapes against the quartic
value 1/8; (b) error of the barrier area identified from a sampled static force (data/identification.pkl)."""

import numpy as np
from scipy.optimize import brentq
from harmonic_gaussian.theory import Q_gauss
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, ORANGE, AQUA, YELLOW, TXT, TXT2

setup()
LN2 = np.log(2)


def Q_shape(phi, dphi, curv, r):
    """Q of x^2/2 + g phi(x) with g chosen so that the stiffness ratio is r."""
    g = r / curv
    dU = lambda x: x + g * dphi(x)
    a_ = brentq(dU, 1e-6, 60)
    h = g * phi(0) - (a_**2 / 2 + g * phi(a_))
    eps = 1e-5
    kw = (dU(a_ + eps) - dU(a_ - eps)) / (2 * eps)
    return h / (kw * a_**2)


fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.8, 2.4), constrained_layout=True)
rs = np.geomspace(1.02, 30, 90)
QL = [Q_shape(lambda x: 1 / (1 + x * x), lambda x: -2 * x / (1 + x * x) ** 2, 2.0, r) for r in rs]
QD = [Q_shape(lambda x: (1 + x * x) ** -1.5, lambda x: -3 * x * (1 + x * x) ** -2.5, 3.0, r) for r in rs]
b = axs[0]
b.axhline(1 / 8, color=ORANGE, lw=1.2)
b.text(3.2, 0.132, "every quartic: $Q=1/8$", color=TXT, fontsize=7, va="bottom")
b.plot(rs, Q_gauss(rs), color=BLUE)
b.plot(rs, QD, color=AQUA)
b.plot(rs, QL, color=YELLOW)
b.text(9.0, 0.235, "Gaussian", color=TXT, fontsize=7, ha="left")
b.text(4.6, 0.37, "dipole-type and\nLorentzian", color=TXT, fontsize=7, ha="right")
b.plot([2], [(1 - LN2) / (4 * LN2**2)], "o", color=BLUE, ms=4, mec="white", mew=0.6)
b.plot([3.684], [0.220], "*", color=TXT, ms=7, mec="white", mew=0.4)
b.annotate(
    "magnet pair",
    (3.684, 0.220),
    xytext=(1.12, 0.262),
    fontsize=7,
    color=TXT,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=TXT2),
)
b.set_xscale("log")
b.set_xlim(1, 30)
b.set_ylim(0.1, 0.8)
b.set_xticks([1, 2, 5, 10, 30])
b.set_xticklabels(["1", "2", "5", "10", "30"])
b.set_xlabel("stiffness ratio $r=2\\alpha\\gamma/\\beta$")
b.set_ylabel("$Q=h/(k_w a^2)$")
panel(b, "(a)")
grid(b)

I = load("identification")
Xs = I["Xs"]
d = axs[1]
d.plot(Xs, I[("known", 0.01, 0.0)][:, 2], "-o", color=BLUE, ms=3, label="$\\beta$ exact")
d.plot(Xs, I[("known", 0.01, 0.002)][:, 2], "--s", color=BLUE, ms=3, label="$\\delta\\beta=0.2\\%$")
d.plot(Xs, I[("known", 0.01, 0.01)][:, 2], ":^", color=BLUE, ms=3, label="$\\delta\\beta=1\\%$")
m = np.isfinite(I[("fit", 0.01)][:, 2])
d.plot(Xs[m], I[("fit", 0.01)][m, 2], "-D", color=AQUA, ms=3, label="$\\beta$ fitted, $|x|\\geq3$")
d.plot(Xs[m], I[("fit", 0.03)][m, 2], "--D", color=AQUA, ms=3, mfc="white", label="same, 3% noise")
d.set_yscale("log")
d.set_ylim(1e-3, 80)
d.set_xlim(1.3, 8.2)
d.set_xlabel("window half-width $X$")
d.set_ylabel("relative RMS error of $M_0$")
d.legend(loc="upper center", fontsize=6.0, handlelength=2.0, borderaxespad=0.1, ncol=2, columnspacing=0.8)
panel(d, "(b)")
grid(d)
save(fig, "fig03_statics")
