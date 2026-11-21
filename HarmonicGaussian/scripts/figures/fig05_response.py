"""Figure 5: interwell branches (power against frequency) of the Gaussian barrier and of the quartic
matched to (a, k_w), at four forcing levels; dashed where unstable (data/branches.pkl)."""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from harmonic_gaussian.model import beff
from harmonic_gaussian.storage import load
from harmonic_gaussian.style import plt, setup, panel, grid, save, W2, TXT, TXT2, RAMP, ORAMP

setup()
d = load("branches")
FS = [0.083, 0.15, 0.3, 0.6]


def draw(ax, arr, col, lw, mask):
    """Branch in runs of equal stability; points outside the mask (not interwell) are skipped."""
    W, Pe, mu = arr[:, 0], arr[:, 2], arr[:, 3]
    st = mu < 1 + 1e-6
    i = 0
    while i < len(W) - 1:
        j = i
        while j + 1 < len(W) and st[j + 1] == st[i] and mask[j + 1] == mask[i]:
            j += 1
        if mask[i]:
            ax.plot(
                W[i : j + 2],
                Pe[i : j + 2],
                ls="-" if st[i] else (0, (2.2, 1.6)),
                color=col,
                lw=lw if st[i] else lw * 0.8,
                solid_capstyle="round",
            )
        i = j + 1


fig, axs = plt.subplots(1, 2, figsize=(W2, 2.7), constrained_layout=True)
for ax, key, R, xl in [(axs[0], "gauss", RAMP, (0.3, 1.25)), (axs[1], "duff", ORAMP, (0.3, 5.6))]:
    for k, f in enumerate(FS):
        hb = d[(key, f)]["high"]
        draw(ax, hb, R[k], 1.4, (hb[:, 4] > 0.05) & (hb[:, 5] < -0.05))
    ax.set_yscale("log")
    ax.set_xlim(*xl)
    ax.set_ylim(1e-4, 2.0)
    ax.set_xlabel("forcing frequency $\\Omega$")
    grid(ax)
    ax.legend(
        handles=[Line2D([], [], color=R[k], lw=1.4, label=f"$f={f}$") for k, f in enumerate(FS)],
        loc="upper left",
        fontsize=7,
        frameon=False,
    )
sb = np.sqrt(beff(0.05))
axs[0].axvline(sb, color=TXT2, lw=0.6, ls=":")
axs[0].text(sb + 0.01, 1.3e-4, "$\\sqrt{\\beta_{\\rm eff}}$", color=TXT2, fontsize=7.5)
axs[0].set_ylabel("electrical power $\\overline{P}_e$")
axs[0].set_title("Gaussian barrier", loc="left", color=TXT, fontsize=8.5)
axs[1].set_title("quartic matched to $(a,k_w)$", loc="left", color=TXT, fontsize=8.5)
panel(axs[0], "(a)")
panel(axs[1], "(b)")
axs[1].add_artist(
    Legend(
        axs[1],
        [Line2D([], [], color=TXT2, lw=1.4), Line2D([], [], color=TXT2, lw=1.1, ls=(0, (2.2, 1.6)))],
        ["stable", "unstable"],
        loc="lower right",
        fontsize=7,
        frameon=False,
    )
)
save(fig, "fig05_response")
