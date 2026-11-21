"""Figure 8: reaching the resonance from rest (baseline Gaussian, xi = 0.01, lambda = 0.05).
(a) capture probability, (b) composition of final states, (c) expected power relative to the linear
harvester, (d) settling time of captured starts, on the 181 x 121 map of data/final_states.pkl;
(e) capture at Omega = 1 against forcing for three damping levels (data/capture_curves.pkl);
(f) slow frequency sweeps (data/capture_at_peak.pkl). Lines from continuation: data/scaling.pkl,
data/lowend.pkl; Melnikov threshold: data/melnikov.pkl."""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm, ListedColormap, BoundaryNorm, to_rgb
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import NullFormatter, FixedLocator
from matplotlib.lines import Line2D
from harmonic_gaussian.analysis import map_fields
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, ORANGE, TXT, TXT2, RAMP, HALO

setup()
M = map_fields(load("final_states"))
Ws, fs = M["Ws"], M["fs"]
P1, PS, PC, PI, GAIN, TS = (M[k] for k in ["P1", "PS", "PC", "PI", "GAIN", "TS"])
S = load("scaling")
L = {round(q["f"], 4): q for q in load("lowend") if q["delta"] == 0.0}
MEL = load("melnikov")[("sym", 0.01)]
R = load("capture_at_peak")
C = load("capture_curves")
fe = np.sqrt(fs[1:] * fs[:-1])
fe = np.r_[fs[0] ** 2 / fe[0], fe, fs[-1] ** 2 / fe[-1]]  # cell edges
dW = Ws[1] - Ws[0]
We = np.r_[Ws - dW / 2, Ws[-1] + dW / 2]
G = [o for o in S["gauss"][0.01] if o.get("ok")]
fg = np.array([o["f"] for o in G])
Wpk = np.array([o["W_peak"] for o in G])
Wfo = np.array([o["W_fold"] for o in G])
pf = [(L[round(o["f"], 4)]["W"], o["f"]) for o in G if round(o["f"], 4) in L and L[round(o["f"], 4)]["type"] == "+1"]


def frame(ax, ylabel=True, xlabel=True):
    ax.set_yscale("log")
    ax.set_xlim(We[0], We[-1])
    ax.set_ylim(fe[0], fe[-1])
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.03, 0.1, 0.3, 1.0]))
    ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3", "1"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    for s in ["top", "right"]:
        ax.spines[s].set_visible(True)
    if xlabel:
        ax.set_xlabel("forcing frequency $\\Omega$")
    if ylabel:
        ax.set_ylabel("forcing $f$")
    else:
        ax.set_yticklabels([])


def branches(ax, lw=1.0, col=TXT, pitch=True):
    ax.plot(Wpk, fg, "-", color=col, lw=lw, path_effects=HALO)
    ax.plot(Wfo, fg, "--", color=col, lw=lw * 0.95, path_effects=HALO)
    if pitch:
        ax.plot([q[0] for q in pf], [q[1] for q in pf], ":", color=col, lw=lw * 1.2, path_effects=HALO)


fig = plt.figure(figsize=(W2, 7.9), layout="constrained")
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.82])
a = fig.add_subplot(gs[0, 0])
b = fig.add_subplot(gs[0, 1])
c = fig.add_subplot(gs[1, 0])
d = fig.add_subplot(gs[1, 1])
e = fig.add_subplot(gs[2, 0])
g = fig.add_subplot(gs[2, 1])

# (a) capture probability
cmap = LinearSegmentedColormap.from_list(
    "seq", ["#ffffff", "#dbe9fb", "#a9cbf4", "#6aa2ea", "#3176d4", "#1c5cab", "#0d366b"]
)
m = a.pcolormesh(We, fe, P1, cmap=cmap, vmin=0, vmax=1, shading="flat", rasterized=True)
a.contour(Ws, fs, P1, levels=[0.5], colors=[TXT], linewidths=0.45)
branches(a)
a.plot(MEL["Ws"], MEL["fM"], "-", color=ORANGE, lw=1.2, path_effects=HALO)
a.text(1.135, 0.52, "fold", color=TXT, fontsize=7, path_effects=HALO, ha="left", clip_on=True)
a.text(0.99, 0.62, "peak", color=TXT, fontsize=7, path_effects=HALO, ha="right", clip_on=True)
a.text(0.455, 0.1, "pitchfork", color=TXT, fontsize=7, path_effects=HALO, ha="left", clip_on=True)
a.text(0.47, 0.0205, "Melnikov", color=ORANGE, fontsize=7, path_effects=HALO, clip_on=True)
frame(a, xlabel=False)
cb = fig.colorbar(m, ax=a, pad=0.015, fraction=0.05, aspect=28)
cb.set_label("capture probability", fontsize=7.5)
cb.ax.tick_params(labelsize=7)
panel(a, "(a)")

# (b) composition of final states, colours mixed in proportion to the probabilities
cols = {"res": "#1c5cab", "sub": "#eda100", "cha": "#c4314b", "intra": "#d4d2cc"}
RGB = np.zeros(P1.shape + (3,))
for p_, k in [(P1, "res"), (PS, "sub"), (PC, "cha"), (PI, "intra")]:
    RGB += p_[:, :, None] * np.array(to_rgb(cols[k]))
RGB = np.clip(RGB / np.maximum(P1 + PS + PC + PI, 1e-12)[:, :, None], 0, 1)
b.imshow(
    RGB,
    origin="lower",
    aspect="auto",
    extent=(0, 1, 0, 1),
    transform=b.transAxes,
    interpolation="nearest",
    rasterized=True,
)
branches(b, lw=0.8, pitch=False)
frame(b, ylabel=False, xlabel=False)
smb = ScalarMappable(
    norm=BoundaryNorm([0, 1, 2, 3, 4], 4), cmap=ListedColormap([cols["intra"], cols["cha"], cols["sub"], cols["res"]])
)
cbb = fig.colorbar(smb, ax=b, pad=0.015, fraction=0.05, aspect=28)
cbb.set_ticks([0.5, 1.5, 2.5, 3.5])
cbb.set_ticklabels(["intrawell", "chaotic\ninterwell", "period 2, 3\ninterwell", "period 1\ninterwell"])
cbb.ax.tick_params(labelsize=6.8, length=0)
cbb.set_label("final state", fontsize=7.5)
panel(b, "(b)")

# (c) expected power from rest relative to the linear harvester
dcm = LinearSegmentedColormap.from_list("div", ["#b0273a", "#e4808a", "#f3f2ef", "#7fb0ec", "#1c5cab"])
mc = c.pcolormesh(We, fe, GAIN, cmap=dcm, norm=LogNorm(vmin=1 / 30, vmax=30), shading="flat", rasterized=True)
c.contour(Ws, fs, GAIN, levels=[1.0], colors=[TXT2], linewidths=0.6)
branches(c, lw=0.8, pitch=False)
frame(c)
cbc = fig.colorbar(mc, ax=c, pad=0.015, fraction=0.05, aspect=28)
cbc.set_label("$\\mathbb{E}[\\overline{P}_{\\!e}]\\,/\\,\\overline{P}_{\\!\\mathrm{lin}}$", fontsize=7.5)
cbc.set_ticks([1 / 30, 0.1, 0.3, 1, 3, 10, 30])
cbc.set_ticklabels(["1/30", "0.1", "0.3", "1", "3", "10", "30"])
cbc.ax.tick_params(labelsize=7)
cbc.ax.yaxis.set_minor_formatter(NullFormatter())
panel(c, "(c)")

# (d) settling time of the captured starts (median, in forcing periods)
tcm = LinearSegmentedColormap.from_list("teal", ["#e3f4ee", "#9fd9c3", "#4fb894", "#1b8a64", "#0c5a40", "#06311f"])
md = d.pcolormesh(
    We, fe, np.ma.masked_invalid(TS), cmap=tcm, norm=LogNorm(vmin=20, vmax=1000), shading="flat", rasterized=True
)
d.set_facecolor("#f3f2ef")
branches(d, lw=0.8, pitch=False)
frame(d, ylabel=False)
cbd = fig.colorbar(md, ax=d, pad=0.015, fraction=0.05, aspect=28)
cbd.set_label("median settling time (periods)", fontsize=7.5)
cbd.ax.tick_params(labelsize=7)
panel(d, "(d)")

# (e) capture at Omega=1 against forcing, three damping levels
cl = {0.005: RAMP[0], 0.01: RAMP[1], 0.02: RAMP[3]}
mks = {0.005: "o", 0.01: "s", 0.02: "^"}
for xi_ in sorted(C):
    r = np.array(C[xi_]["rows"])
    e.plot(r[:, 0], r[:, 1], "-", marker=mks[xi_], ms=2.6, lw=1.1, color=cl[xi_], label=f"$\\xi={xi_}$")
e.axvline(0.21715, color=TXT2, lw=0.7, ls=":")
e.text(0.2185, 0.06, "$f_\\mathrm{cr}$", color=TXT2, fontsize=7)
e.set_xlim(0.15, 0.28)
e.set_ylim(-0.03, 1.05)
e.set_xlabel("forcing $f$  ($\\Omega=1$)")
e.set_ylabel("capture probability")
for xi_, xl, yl in [(0.005, 0.1945, 0.72), (0.01, 0.2160, 0.45), (0.02, 0.2515, 0.45)]:
    e.text(
        xl,
        yl,
        f"$\\xi={xi_}$",
        color=cl[xi_] if xi_ != 0.005 else "#5d93d6",
        fontsize=7,
        ha="right",
        va="center",
        path_effects=HALO,
    )
panel(e, "(e)")
grid(e)

# (f) frequency sweeps ending at the peak frequency
mk = {0.12: "o", 0.15: "s", 0.2: "^"}
colr = {0.12: RAMP[0], 0.15: RAMP[1], 0.2: RAMP[3]}
for delta, ls, fill in [(0.0, "-", True), (0.1, "--", False)]:
    for f in [0.12, 0.15, 0.2]:
        q = [r for r in R["ramps"] if r["delta"] == delta and r["f"] == f]
        g.plot(
            [r["Tr"] for r in q],
            [r["cap"] for r in q],
            ls,
            marker=mk[f],
            color=colr[f],
            mfc=colr[f] if fill else "white",
            ms=3.6,
            lw=1.1,
        )
g.set_xscale("log")
g.set_ylim(-0.03, 1.05)
g.set_xlabel("sweep duration $T_r$")
g.set_ylabel("fraction on resonant branch")
h = [Line2D([], [], color=colr[f], marker=mk[f], ls="-", ms=3.6, label=f"$f={f}$") for f in mk]
h += [
    Line2D([], [], color=TXT2, marker="o", mfc=TXT2, ls="-", ms=3.4, label="$\\delta=0$"),
    Line2D([], [], color=TXT2, marker="o", mfc="white", ls="--", ms=3.4, label="$\\delta=0.1$"),
]
g.legend(
    handles=h,
    loc="upper left",
    bbox_to_anchor=(0.0, 0.97),
    ncol=2,
    columnspacing=0.8,
    fontsize=6.6,
    labelspacing=0.25,
    frameon=False,
)
panel(g, "(f)")
grid(g)

save(fig, "fig08_access", dpi=450)
