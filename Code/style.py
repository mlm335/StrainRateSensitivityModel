"""style.py -- one place for the plot style and for saving figures."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
RATE_COLORS = ["#1b4f72", "#148f77", "#d4ac0d", "#ca6f1e", "#922b21"]
DATA_COLORS = ["#c0392b", "#2471a3", "#239b56", "#7d3c98", "#d68910", "#6e2c00", "#17202a"]
MARKERS = ["o", "s", "^", "D", "v", "P", ">", "<"]

plt.rcParams.update({"font.size": 11, "font.family": "serif", "mathtext.fontset": "dejavuserif",
                     "axes.linewidth": 0.9, "lines.linewidth": 1.8, "legend.frameon": False,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
                     "ytick.right": True, "figure.figsize": (5.2, 4.0), "savefig.dpi": 300,
                     "savefig.bbox": "tight"})


def new_figure():
    return plt.subplots()


def save(fig, name):
    """Save one figure as PNG and PDF in Code/figures and close it."""
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote figures/" + name + ".png/.pdf")


def rate_label(e):
    """'6.4e-04' -> '$6.4\\times10^{-4}$ s$^{-1}$'"""
    import numpy as np
    k = int(np.floor(np.log10(e)))
    m = e / 10.0**k
    head = f"{m:.1f}".rstrip("0").rstrip(".")
    return (fr"${head}\times10^{{{k}}}$ s$^{{-1}}$" if head != "1" else fr"$10^{{{k}}}$ s$^{{-1}}$")
