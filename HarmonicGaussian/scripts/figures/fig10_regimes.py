"""Figure 10: (a) operating regimes in the forcing-damping plane (pinning, cusp, capture threshold,
hardening and stroke limits); (b) expected power from rest against the linear harvester at two
forcing levels; (c) band-averaged power relative to the linear harvester
(data/fold_curve.pkl, data/capture_curves.pkl, data/expected_power.pkl)."""

import numpy as np
import matplotlib.patheffects as pe
from harmonic_gaussian.model import beff, P_lin
from harmonic_gaussian.theory import M0_GAUSS as M0
from harmonic_gaussian.analysis import capture_thresholds, band_ratio
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, BLUE, ORANGE, AQUA, TXT, TXT2, RAMP

setup()
sb = np.sqrt(beff(0.05))
fmax = load("fold_curve")["fc"]
pts = capture_thresholds(load("capture_curves"))  # (c, f just below 50% capture, f at or above)
cs = np.geomspace(0.004, 0.4, 300)
fl = (2 * M0 / np.pi) ** (1 / 3) * (cs * sb) ** (2 / 3)  # eps = 1
k3 = 1e-4
fh = np.sqrt(4 / (3 * k3)) * (cs * sb) ** 1.5  # hardening ceiling
Amax = 10
fs_ = Amax * cs * sb  # stroke limit
cP = np.array([q[0] for q in pts])
fP = np.array([0.5 * (q[1] + q[2]) for q in pts])
fthr = np.interp(np.log(cs), np.log(cP), fP)
fig, axs = plt.subplots(1, 3, figsize=(W2, 2.45), constrained_layout=True)
a, b, c = axs
lo = np.maximum(fl, fthr)
hi = np.minimum(np.minimum(fh, fs_), fmax)
a.fill_betweenx(cs, np.maximum(fl, 0), np.minimum(fmax, 10), where=fl < fmax, color=BLUE, alpha=0.10, lw=0)
a.fill_betweenx(cs, lo, hi, where=lo < hi, color=BLUE, alpha=0.35, lw=0)
a.plot(fl, cs, "-", color=BLUE, lw=1.2)
a.plot(3 ** (1 / 3) * fl, cs, ":", color=BLUE, lw=1.0)
a.axvline(fmax, color=TXT, lw=1.0, ls="--")
a.axvline(6 * M0 / np.pi, color=TXT2, lw=0.7, ls=":")
a.plot(fh, cs, "-", color=ORANGE, lw=1.1)
a.plot(fs_, cs, "-", color=AQUA, lw=1.1)
for cc, f1, f2 in pts:
    a.plot([f1, f2], [cc, cc], "-", color=TXT, lw=1.6, solid_capstyle="butt")
    a.plot([0.5 * (f1 + f2)], [cc], "o", color=TXT, ms=3)
a.plot(fthr, cs, "-", color=TXT, lw=0.6, alpha=0.6)
a.set_xscale("log")
a.set_yscale("log")
a.set_xlim(0.02, 5)
a.set_ylim(0.005, 0.3)
a.set_xlabel("forcing $f$")
a.set_ylabel("total damping $c$")
HL = [pe.withStroke(linewidth=2.2, foreground="white")]
a.text(fmax * 1.06, 0.0062, "cusp", color=TXT, fontsize=7, rotation=90, path_effects=HL)
a.text(0.0215, 0.0105, "$\\varepsilon=1$", color=BLUE, fontsize=7, path_effects=HL)
a.text(1.85, 0.045, "$k_3=10^{-4}$", color=ORANGE, fontsize=6.5, path_effects=HL)
a.text(2.05, 0.125, "$A_p=10$", color=AQUA, fontsize=6.5, path_effects=HL)
a.text(0.088, 0.0066, "capture", color=TXT, fontsize=7, path_effects=HL)
panel(a, "(a)")
grid(a)
P = load("expected_power")
Ws, fsP, PE = P["Ws"], P["fs"], P["PE"]
Wf = np.linspace(0.4, 1.4, 20001)
for f, col in [(0.15, RAMP[1]), (0.3, RAMP[3])]:
    i = int(np.argmin(abs(fsP - f)))
    b.plot(Wf, P_lin(Wf, f), "-", color=TXT2, lw=0.8)
    b.plot(Ws, PE[i], "-", color=col, lw=1.2, marker="o", ms=1.8)
b.set_yscale("log")
b.set_xlim(0.4, 1.4)
b.set_ylim(1e-5, 2)
b.text(0.42, 0.6, "$f=0.3$", color=RAMP[3], fontsize=7)
b.text(0.42, 1.5e-4, "$f=0.15$", color=RAMP[1], fontsize=7)
b.annotate(
    "linear harvester",
    xy=(1.2, 3e-4),
    xytext=(1.08, 0.04),
    fontsize=7,
    color=TXT2,
    arrowprops=dict(arrowstyle="-", lw=0.5, color=TXT2),
)
b.set_xlabel("forcing frequency $\\Omega$")
b.set_ylabel("expected power $\\mathbb{E}[\\overline{P}_{\\!e}\\,|\\,\\Omega]$")
panel(b, "(b)")
grid(b)
for (W0, W1), ls, lab in [
    ((0.4, 1.4), "-", "$0.4\\leq\\Omega\\leq1.4$"),
    ((0.9, 1.1), "--", "$0.9\\leq\\Omega\\leq1.1$"),
    ((0.4, 0.9), ":", "$0.4\\leq\\Omega\\leq0.9$"),
]:
    c.plot(
        fsP, [band_ratio(Ws, PE[i], f, W0, W1) for i, f in enumerate(fsP)], ls, color=BLUE, marker="o", ms=3, label=lab
    )
c.axhline(1, color=TXT2, lw=0.7)
c.set_yscale("log")
c.set_xlabel("forcing $f$")
c.set_ylabel("band power, magnets / linear")
c.legend(fontsize=6.3, loc="lower right", handlelength=3.2)
panel(c, "(c)")
grid(c)
save(fig, "fig10_regimes")
