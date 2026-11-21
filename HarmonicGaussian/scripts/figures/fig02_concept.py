"""Figure 2: how the harmonic-Gaussian potential is built, and the orbits a harvester can follow.
Computed from U(x) = x^2/2 + exp(-x^2) (alpha = beta = gamma = 1, delta = 0)."""

import numpy as np
from matplotlib.patches import FancyArrowPatch
from harmonic_gaussian.style import plt, setup, panel, save, W2, BLUE, ORANGE, AQUA, TXT, TXT2, GRID, RAMP

setup()
LN2 = np.log(2)
U = lambda x: x**2 / 2 + np.exp(-(x**2))
Um = lambda x: np.exp(-(x**2))
a = np.sqrt(LN2)
Ua = U(a)
U0 = U(0.0)
k1, k3 = LN2, 1.0  # quartic matched to (a, k_w)
UD = lambda x: -k1 * x**2 / 2 + k3 * x**4 / 4 - (-k1 * a**2 / 2 + k3 * a**4 / 4) + Ua


def darrow(ax, p, q, color=TXT2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>", mutation_scale=6, lw=0.6, color=color, shrinkA=0, shrinkB=0))


fig, axs = plt.subplots(1, 2, figsize=(W2, 2.45), constrained_layout=True)

# (a) anatomy of the potential
ax = axs[0]
x = np.linspace(-3.6, 3.6, 1200)
for s in (-1, 1):
    ax.axvspan(s * 2.3, s * 3.6, color=GRID, alpha=0.55, lw=0)
    ax.text(s * 2.95, 0.12, "linear\nfar field", ha="center", va="bottom", fontsize=7, color=TXT2)
ax.fill_between(x, 0, Um(x), color=AQUA, alpha=0.22, lw=0)
(lb,) = ax.plot(x, Um(x), color=AQUA, lw=1.0, label="magnetic barrier $U_m$")
(lp,) = ax.plot(x, x**2 / 2, ":", color=TXT2, lw=1.0, label="beam, $\\beta x^2/2$")
(lq,) = ax.plot(x, UD(x), color=ORANGE, lw=1.1, label="quartic, same $a$ and $k_w$")
(lu,) = ax.plot(x, U(x), color=BLUE, lw=1.6, label="$U=\\beta x^2/2+U_m$")
ax.legend(
    handles=[lu, lp, lb, lq], loc="upper center", fontsize=6.6, handlelength=1.6, borderaxespad=0.2, labelspacing=0.25
)
ax.text(0.0, 0.32, "area $M_0$", ha="center", fontsize=7, color=TXT)
ax.plot([-1.75, 0], [U0, U0], lw=0.5, color=TXT2, ls="--")  # barrier height h and well position a
ax.plot([-1.75, -a], [Ua, Ua], lw=0.5, color=TXT2, ls="--")
darrow(ax, (-1.6, Ua), (-1.6, U0))
ax.text(-1.67, (Ua + U0) / 2, "$h$", ha="right", va="center", fontsize=7.5, color=TXT)
ax.plot([a, a], [-0.05, Ua], lw=0.5, color=TXT2, ls="--")
ax.plot([-a, -a], [-0.05, Ua], lw=0.5, color=TXT2, ls="--")
ax.set_xticks([-3, -a, 0, a, 3])
ax.set_xticklabels(["$-3$", "$-a$", "0", "$a$", "3"])
ax.set_xlim(-3.6, 3.6)
ax.set_ylim(0, 3.5)
ax.set_xlabel("displacement $x$")
ax.set_ylabel("potential energy")
panel(ax, "(a)")

# (b) phase plane of the free oscillator
ax = axs[1]


def orbit(E, branch):
    xs = np.linspace(-4, 4, 40001)
    v2 = 2 * (E - U(xs))
    m = v2 >= 0
    if branch == "right":
        m &= xs > 0
    if branch == "left":
        m &= xs < 0
    xx = xs[m]
    vv = np.sqrt(v2[m])
    return np.r_[xx, xx[::-1], xx[:1]], np.r_[vv, -vv[::-1], vv[:1]]


for E, col in [(Ua + 0.02, RAMP[0]), (Ua + 0.09, RAMP[1])]:
    for br in ("left", "right"):
        X, V = orbit(E, br)
        ax.plot(X, V, color=col, lw=1.0)
Xs, Vs = orbit(U0, "all")
ax.plot(Xs, Vs, "--", color=TXT2, lw=0.9)
A = 3.0
X, V = orbit(U(A), "all")
ax.plot(X, V, color=BLUE, lw=1.7)
ax.plot([-a, a], [0, 0], "o", color=TXT2, ms=2.5)
ax.plot([0], [0], "x", color=TXT2, ms=4, mew=0.8)
darrow(ax, (0, -3.15), (A, -3.15), BLUE)
ax.text(-0.08, -3.15, "amplitude $A$", ha="right", va="center", fontsize=7, color=TXT)
ax.annotate(
    "intrawell orbits:\nsmall amplitude",
    (a + 0.35, 0.28),
    xytext=(1.55, 1.25),
    fontsize=7,
    color=TXT,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=TXT2),
)
ax.annotate(
    "separatrix",
    (-1.1, 0.62),
    xytext=(-2.5, 1.55),
    fontsize=7,
    color=TXT2,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=TXT2),
)
ax.text(0, 2.95, "interwell (high-energy) orbit", ha="center", fontsize=7, color=TXT)
ax.set_xlim(-3.6, 3.6)
ax.set_ylim(-3.45, 3.45)
ax.set_xlabel("displacement $x$")
ax.set_ylabel("velocity $\\dot x$")
panel(ax, "(b)")

save(fig, "fig02_concept")
