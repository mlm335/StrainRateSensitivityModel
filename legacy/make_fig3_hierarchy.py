"""Fig. 3: nested model hierarchy (W level -> Zr level -> Eurofer level) and Lambda."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from material_cards import *
from material_cards import _diff

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.8, "lines.linewidth": 1.4,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False, "legend.fontsize": 6.3, "savefig.dpi": 400,
                     "font.family": "serif", "mathtext.fontset": "dejavuserif"})
LV = ["0.55", "#2471a3", "#c0392b"]


def lab(ax, t, right=False):
    ax.text(0.97 if right else 0.03, 0.95, t, transform=ax.transAxes, va="top",
            ha="right" if right else "left", fontweight="bold")


fig, axs = plt.subplots(1, 4, figsize=(7.0, 2.05))
# (a) Zr: level I -> level II
TZ = np.linspace(400, 950, 111); e = 3.3e-5
for (kw, name, c, ls) in [(dict(aging=False, climb=False), "I", LV[0], "--"),
                          (dict(climb=False), "I + aging", LV[1], ":"),
                          ({}, "II (+ recovery)", LV[1], "-")]:
    axs[0].plot(TZ, make_Zr(**kw).solve(e / ZR_M, TZ)["tau"] / 1e6, color=c, ls=ls, label=name)
T, y = ZR_DATA["crss"][e]
axs[0].scatter(T, y, s=12, facecolors="none", edgecolors="k", lw=0.8, zorder=5)
axs[0].set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", ylim=(0, 90), xlim=(400, 950))
axs[0].legend(loc="lower left"); lab(axs[0], "(a) Zr", True)
# (b) Eurofer tensile: I -> II -> III
TE = np.linspace(300, 900, 121); e = 3e-4
for (kw, name, c, ls) in [(dict(aging=False, precipitates=False, recovery=False), "I", LV[0], "--"),
                          (dict(precipitates=False), "II", LV[1], "-"),
                          ({}, "III (+ particles)", LV[2], "-")]:
    axs[1].plot(TE, make_Eurofer(**kw).solve(e / EU_M, TE)["tau"] / 1e6, color=c, ls=ls, label=name)
y = EU_DATA["sigma"][:, 1] / 3
axs[1].scatter(EU_DATA["T"], y, s=12, facecolors="none", edgecolors="k", lw=0.8, zorder=5)
axs[1].set(xlabel=r"$T$ (K)", ylabel=r"$\sigma/3$, $\tau$ (MPa)", ylim=(0, 270), xlim=(300, 900))
axs[1].legend(loc="center left", bbox_to_anchor=(0, 0.42)); lab(axs[1], "(b) Eurofer-97", True)
# (c) Eurofer creep at 823 K: with / without the climb route
sig = np.linspace(120, 300, 60) * 1e6
for (kw, name, c, ls) in [(dict(climb=False), "III, no climb", LV[2], "--"), ({}, "III", LV[2], "-")]:
    cm = CreepModel(make_Eurofer(**kw), None)
    r = cm.rate(sig[None, :] * EU_M, np.array([823.0]))
    axs[2].plot(sig / 1e6, r["total"][0] * EU_M, color=c, ls=ls, label=name)
k = EU_CREEP[:, 0] == 823
axs[2].scatter(EU_CREEP[k, 1], EU_CREEP[k, 2], s=12, facecolors="none", edgecolors="k", lw=0.8, zorder=5)
axs[2].set(xscale="log", yscale="log", xlabel=r"$\sigma$ (MPa)", ylabel=r"$\dot\epsilon_{\min}$ (s$^{-1}$)",
           xlim=(120, 300), ylim=(1e-11, 1e-5))
axs[2].set_xticks([150, 200, 300]); axs[2].set_xticklabels(["150", "200", "300"])
axs[2].minorticks_off() if False else None
axs[2].legend(loc="lower right"); lab(axs[2], "(c) 823 K")
# (d) Lambda
for y0, y1, c in [(-1e6, -1, "#d6e6f5"), (-1, 1, "#eef5ee"), (1, 1e6, "#f8e1e1")]:
    axs[3].axhspan(y0, y1, color=c, zorder=0)
TW = np.linspace(200, 900, 141); TZ = np.linspace(300, 1000, 141)
for (name, mdl, M, T, er), c in zip([("W", make_W(), W_M, TW, 5e-3), ("Zr", make_Zr(), ZR_M, TZ, 1e-4),
                                     ("Eurofer-97", make_Eurofer(), EU_M, TE, 3e-4)], ["k", LV[1], LV[2]]):
    axs[3].plot(T, mdl.srs(er / M, T)["Lambda"], color=c, label=name)
axs[3].set_yscale("symlog", linthresh=0.1)
axs[3].set(xlabel=r"$T$ (K)", ylabel=r"$\Lambda=AB$", ylim=(-30, 1e4), xlim=(200, 1000))
axs[3].text(990, -9, r"$m<0$", ha="right", fontsize=6.5)
axs[3].text(990, 2.2, r"$B$-controlled", ha="right", va="bottom", fontsize=6.5)
axs[3].text(990, 0.22, r"$A$-controlled", ha="right", fontsize=6.5)
axs[3].text(880, 0.03, "W", ha="right", fontsize=7)
axs[3].text(800, 2500, "Zr", ha="right", color=LV[1], fontsize=7)
axs[3].text(575, 25, "Eurofer", ha="left", color=LV[2], fontsize=7)
lab(axs[3], "(d)")
plt.tight_layout(pad=0.3, w_pad=0.5)
fig.savefig(f"{OUT}/fig3_hierarchy.pdf"); fig.savefig(f"{OUT}/fig3_hierarchy.png")
