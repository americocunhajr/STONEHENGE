"""Figure 4: (a) free-oscillation frequency against amplitude for the Gaussian barrier and the
matched quartic (exact, harmonic balance, area law); (b) frequency deficit over the area law for
four barriers of equal area, including the point-dipole magnet pair. Computed by quadrature."""

import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import brentq
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from harmonic_gaussian.model import magnet_barrier
from harmonic_gaussian.theory import omega_inter, omega_intra, keff_hb_gauss, M0_GAUSS as M0
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, ORANGE, AQUA, YELLOW, TXT, TXT2

setup()
warnings.simplefilter("ignore", IntegrationWarning)  # integrable endpoint singularity of orbits close to the separatrix
LN2 = np.log(2)
wL = (1 / np.sqrt(np.pi)) ** (1 / 3)
gL = wL**2
wD = (3 * np.sqrt(np.pi) / 4) ** (1 / 3)
gD = 2 * wD**2 / 3
UG = lambda x: np.exp(-x * x) + x * x / 2
UL = lambda x: gL / (1 + (x / wL) ** 2) + x * x / 2
UDp = lambda x: gD * (1 + (x / wD) ** 2) ** -1.5 + x * x / 2
k1, k3 = LN2, 1.0
UD = lambda x: -k1 * x * x / 2 + k3 * x**4 / 4
fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.82, 2.55), constrained_layout=True)

# (a) backbones
a = axs[0]
xw = np.sqrt(LN2)
Uw = UG(xw)
As = brentq(lambda x: UG(x) - 1.0, 0.5, 5)
Es = Uw + (1 - Uw) * np.r_[np.geomspace(1e-4, 0.5, 40), 1 - np.geomspace(0.5, 1e-6, 50)[1:]]
intra = [omega_intra(UG, xw, 0.0, E) for E in Es]
a.plot([x for w, x in intra], [w for w, x in intra], color=BLUE)
Ai = np.r_[As + np.geomspace(1e-6, 0.3, 40), np.linspace(As + 0.31, 7, 120)]
a.plot(Ai, [omega_inter(UG, A) for A in Ai], color=BLUE)
Ah = np.linspace(1.5, 7, 200)
a.plot(Ah, np.sqrt(keff_hb_gauss(Ah)), "--", color=BLUE, lw=0.9)
Aa = np.linspace(1.75, 7, 200)
a.plot(Aa, np.sqrt(1 - 2 * M0 / (np.pi * Aa**3)), ":", color=TXT2, lw=1.1)
xwd = np.sqrt(k1 / k3)
Uwd = UD(xwd)
Asd = np.sqrt(2 * k1 / k3)
Esd = Uwd + (0 - Uwd) * np.r_[np.geomspace(1e-4, 0.5, 40), 1 - np.geomspace(0.5, 1e-6, 50)[1:]]
intrad = [omega_intra(UD, xwd, 0.0, E) for E in Esd]
a.plot([x for w, x in intrad], [w for w, x in intrad], color=ORANGE)
Aid = np.r_[Asd + np.geomspace(1e-6, 0.3, 40), np.linspace(Asd + 0.31, 4.0, 100)]
a.plot(Aid, [omega_inter(UD, A) for A in Aid], color=ORANGE)
a.axhline(1.0, color=TXT2, lw=0.6, ls=(0, (1, 2)))
a.text(6.95, 1.03, "$\\sqrt{\\beta}$", color=TXT2, fontsize=7.5, ha="right", va="bottom")
a.text(2.65, 1.45, "Duffing", color=TXT, fontsize=7.5)
a.text(4.4, 0.86, "Gaussian", color=TXT, fontsize=7.5)
a.text(0.62, 1.24, "intrawell", color=TXT2, fontsize=7)
a.text(1.45, 0.33, "interwell", color=TXT2, fontsize=7)
a.annotate(
    "separatrix",
    (As, 0.0),
    xytext=(1.9, 0.12),
    fontsize=7,
    color=TXT2,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=TXT2),
)
a.set_xlim(0.5, 7)
a.set_ylim(0, 2.6)
a.set_xlabel("amplitude $A$ (turning point)")
a.set_ylabel("free-oscillation frequency $\\omega$")
panel(a, "(a)")
grid(a)
a.legend(
    handles=[
        Line2D([], [], color=TXT2, lw=1.4, label="exact"),
        Line2D([], [], color=TXT2, lw=0.9, ls="--", label="harmonic balance"),
        Line2D([], [], color=TXT2, lw=1.1, ls=":", label="area law"),
    ],
    loc="upper left",
    fontsize=7,
    handlelength=2.2,
)

# (b) deficit relative to the area law
b = axs[1]
A2 = np.geomspace(2.0, 30, 40)
UM = lambda x: magnet_barrier(x) + x * x / 2
M0M = quad(magnet_barrier, -200, 200, limit=2000)[0]
for U, M, col, lab, ls in [
    (UG, M0, BLUE, "Gaussian", "-"),
    (UDp, M0, AQUA, "dipole-type", "-"),
    (UL, M0, YELLOW, "Lorentzian", "-"),
    (UM, M0M, TXT, "magnet pair", "--"),
]:
    d = np.array([1 - omega_inter(U, A) ** 2 for A in A2])
    b.plot(A2, d / (2 * M / (np.pi * A2**3)), ls, color=col, lw=1.3, label=lab)
b.axhline(1, color=TXT2, lw=0.7, ls=":")
b.legend(loc="upper right", fontsize=7)
b.set_xscale("log")
b.set_xlim(2, 30)
b.set_ylim(0.75, 1.6)
b.set_xticks([2, 5, 10, 20, 30])
b.set_xticklabels(["2", "5", "10", "20", "30"])
b.set_xlabel("amplitude $A$")
b.set_ylabel("deficit over area law, $(\\beta-\\omega^2)\\,\\pi A^3/2M_0$")
panel(b, "(b)")
grid(b)
b.xaxis.set_minor_formatter(NullFormatter())
save(fig, "fig04_backbones")
