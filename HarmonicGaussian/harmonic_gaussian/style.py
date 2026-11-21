"""Common style of the figures: Computer Modern text, a blue ordinal ramp with orange, aqua and
yellow accents, light grids, and the journal's single and double column widths."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from .storage import FIGURES

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TXT, TXT2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # ordinal blue (light -> dark)
ORAMP = ["#f09a78", "#eb6834", "#b04e27", "#81391d"]  # ordinal orange (light -> dark)
MK = {0.005: "o", 0.01: "s", 0.02: "^", 0.04: "D"}  # marker per damping level
W1, W2 = 3.31, 6.85  # single and double column widths (in)
HALO = [pe.withStroke(linewidth=2.4, foreground="white")]


def setup():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["cmr10"],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "axes.unicode_minus": False,
            "font.size": 8,
            "axes.labelsize": 8.5,
            "axes.titlesize": 8.5,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.edgecolor": TXT2,
            "axes.labelcolor": TXT,
            "xtick.color": TXT2,
            "ytick.color": TXT2,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.minor.size": 1.6,
            "ytick.minor.size": 1.6,
            "lines.linewidth": 1.4,
            "lines.markersize": 4,
            "figure.facecolor": SURF,
            "axes.facecolor": SURF,
            "savefig.facecolor": SURF,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )


def panel(ax, label):
    """Panel label (a), (b), ... above the top left corner."""
    ax.text(
        -0.02, 1.02, label, transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5, fontweight="bold", color=TXT
    )


def grid(ax):
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def save(fig, name, dpi=None):
    """Write figures/<name>.pdf."""
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{name}.pdf", **({} if dpi is None else dict(dpi=dpi)))
    plt.close(fig)
