"""Publication figures for the PRL letter.

Usage (from any directory):
    python make_paper_figs.py

Needs numpy and matplotlib, plus unified_plasticity.py, material_cards.py and
make_fig3_hierarchy.py in the same folder. Writes paper/figures/fig2_srs,
fig3_hierarchy and fig4_maps_creep (.pdf and .png). Takes ~1-2 minutes.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from material_cards import *

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "paper", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.8, "lines.linewidth": 1.4,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False, "legend.fontsize": 6.5, "savefig.dpi": 400,
                     "font.family": "serif", "mathtext.fontset": "dejavuserif"})
MK = ["o", "s", "^", "D", "v", "P", ">", "<"]
CL = ["#c0392b", "#2471a3", "#239b56", "#7d3c98", "#d68910", "#6e2c00"]
RC = ["#1b4f72", "#148f77", "#d4ac0d"]


def pts(ax, sets, s=12):
    for k, (n, (x, y)) in enumerate(sets.items()):
        ax.scatter(x, y, s=s, marker=MK[k % 8], facecolors="none", edgecolors=CL[k % 6], linewidths=0.8,
                   label=n, zorder=5)


def lab(ax, t):
    ax.text(0.03, 0.95, t, transform=ax.transAxes, va="top", ha="left", fontweight="bold")


# ---------------------------------------------------------------- Fig. 2
W, Zr, EU = make_W(), make_Zr(), make_Eurofer()
TW = np.linspace(200, 900, 141); TZ = np.linspace(300, 1000, 141); TE = np.linspace(300, 900, 121)
fig, axs = plt.subplots(3, 3, figsize=(7.0, 6.3))
# W (bulk ATW, axial flow stress)
for c, e in zip(RC, W_RATES):
    r = W.solve(e / W_M, TW)
    axs[0, 0].plot(TW, r["tau"] / W_M / 1e6, color=c, label=fr"${e/10**np.floor(np.log10(e)):.1f}\times10^{{{int(np.floor(np.log10(e)))}}}$ s$^{{-1}}$")
    k = np.isclose(W_DATA["flow"][:, 1], e)
    axs[0, 0].scatter(W_DATA["flow"][k, 0], W_DATA["flow"][k, 2], s=16, facecolors="none", edgecolors=c, lw=0.9, zorder=5)
axs[0, 0].plot(TW, r["tau_a"] / W_M / 1e6, "k-.", lw=1, label=r"$\tau_a/M_s$")
axs[0, 0].set(ylabel=r"$\sigma$ (MPa)", xlim=(300, 900), ylim=(400, 1700))
axs[0, 0].legend(loc="upper right", fontsize=5.5)
g = W.global_m_V(W_RATES, TW, W_M)
axs[0, 1].plot(TW, g["m"], "k"); pts(axs[0, 1], W_DATA["m"])
axs[0, 1].set(ylabel=r"$m$", xlim=(200, 850), ylim=(0, 0.06))
axs[0, 2].plot(TW, g["Vstar_b3"], "k"); pts(axs[0, 2], W_DATA["V"])
axs[0, 2].set(ylabel=r"$V^\ast$ ($b^3$)", xlim=(200, 850), ylim=(0, 250))
axs[0, 2].legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=5.5)
# Zr
for c, e in zip(["#2471a3", "#d68910"], [3.3e-3, 3.3e-5]):
    r = Zr.solve(e / ZR_M, TZ)
    axs[1, 0].plot(TZ, r["tau"] / 1e6, color=c, label=fr"{e:.1e} s$^{{-1}}$")
    T, y = ZR_DATA["crss"][e]
    axs[1, 0].scatter(T, y, s=12, facecolors="none", edgecolors=c, lw=0.8)
axs[1, 0].set(ylabel=r"$\tau$ (MPa)", xlim=(300, 1000), ylim=(0, 150)); axs[1, 0].legend(loc="upper right")
g = Zr.global_m_V(ZR_RATES, TZ, ZR_M)
axs[1, 1].plot(TZ, g["m"], "k"); pts(axs[1, 1], ZR_DATA["m"]); axs[1, 1].axhline(0, color="0.5", lw=0.5)
axs[1, 1].set(ylabel=r"$m$", ylim=(-0.02, 0.1), xlim=(280, 1000))
axs[1, 1].legend(loc="lower right", fontsize=4.8, ncol=2, bbox_to_anchor=(1.0, 0.0))
axs[1, 2].plot(TZ, g["Vstar_b3"], "k"); pts(axs[1, 2], ZR_DATA["V"])
axs[1, 2].set(ylabel=r"$V^\ast$ ($b^3$)", ylim=(0, 300), xlim=(50, 1000))
axs[1, 2].legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=5.5)
# Eurofer
for j, (c, e) in enumerate(zip(RC, EU_RATES)):
    r = EU.solve(e / EU_M, TE)
    axs[2, 0].plot(TE, r["tau"] / 1e6, color=c, label=fr"$3\times10^{{{int(np.log10(e/3))}}}$ s$^{{-1}}$")
    y = EU_DATA["sigma"][:, j] / 3; k = np.isfinite(y)
    axs[2, 0].scatter(EU_DATA["T"][k], y[k], s=12, marker=MK[j], facecolors="none", edgecolors=c, lw=0.8)
axs[2, 0].set(ylabel=r"$\sigma/3$, $\tau$ (MPa)", xlim=(300, 900)); axs[2, 0].legend(loc="lower left")
axs[2, 0].scatter(EU_MM[:, 0], EU_MM[:, 1] / 3, s=16, marker="v", facecolors="none", edgecolors="#239b56",
                  lw=0.8, label="Materna-Morris")
axs[2, 0].legend(loc="lower left", fontsize=5.5)
g = EU.global_m_V(EU_RATES, TE, EU_M)
dsa = EU.srs(3e-4 / EU_M, TE)["Lambda"] < -0.5          # aging-dominated window
for ax_ in (axs[2, 1], axs[2, 2]):
    ax_.fill_between(TE, 0, 1, where=dsa, transform=ax_.get_xaxis_transform(), color="0.9", lw=0, zorder=0)
axs[2, 1].plot(TE, g["m"], "k"); axs[2, 1].axhline(0, color="0.5", lw=0.5)
axs[2, 1].scatter(EU_SRS_M[:, 0], EU_SRS_M[:, 1], s=14, facecolors="none", edgecolors=CL[0], lw=0.8,
                  label="Vanaja et al.")
axs[2, 1].scatter(SHAH_SRS[:, 0], SHAH_SRS[:, 1], s=14, marker="v", facecolors="none", edgecolors=CL[4], lw=0.8,
                  label="Shah et al. (IN-RAFM)")
axs[2, 1].set(ylabel=r"$m$", xlim=(300, 900), ylim=(-0.01, 0.12))
axs[2, 1].legend(loc="upper left", bbox_to_anchor=(0, 0.93), fontsize=5.5)
axs[2, 2].plot(TE, g["Vstar_b3"], "k")
axs[2, 2].scatter(EU_SRS_V[:, 0], EU_SRS_V[:, 1], s=14, facecolors="none", edgecolors=CL[0], lw=0.8)
axs[2, 2].scatter(SHAH_SRS[:, 0], SHAH_SRS[:, 2], s=14, marker="v", facecolors="none", edgecolors=CL[4], lw=0.8)
axs[2, 2].set(ylabel=r"$V^\ast$ ($b^3$)", ylim=(0, 350), xlim=(300, 900))
for i, name in enumerate(["W", "Zr", "Eurofer97"]):
    for j in range(3):
        axs[i, j].set_xlabel(r"$T$ (K)" if i == 2 else "")
        lab(axs[i, j], f"({'abcdefghi'[3*i+j]})" + (f" {name}" if j == 0 else ""))
plt.tight_layout(pad=0.3, h_pad=0.4, w_pad=0.6)
fig.savefig(f"{OUT}/fig2_srs.pdf"); fig.savefig(f"{OUT}/fig2_srs.png")
plt.close(fig)

# ---------------------------------------------------------------- Fig. 3 : hierarchy (separate script)
import runpy
runpy.run_path(os.path.join(HERE, "make_fig3_hierarchy.py"))
print("wrote fig2_srs, fig3_hierarchy", flush=True)

# ---------------------------------------------------------------- Fig. 4 : maps + parity + creep
FIELD_COL = ["#e6e6e6", "#d3e3f3", "#f4c9a8", "#cfe8c6"]
cmap = ListedColormap(FIELD_COL); norm = BoundaryNorm([-0.5, 0.5, 1.5, 3.5, 4.5], cmap.N)
SPEC = {"W": (0.05, 0.70), "Zr": (0.14, 1136 / ZR_MAT.Tm), "Eurofer97": (0.16, 1100 / EU_Tm)}
fig = plt.figure(figsize=(7.0, 4.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05])
allp = tensile_points()
for k, name in enumerate(SPEC):
    ax = fig.add_subplot(gs[0, k]); d = MATERIALS[name]; Tm = d["mat"].Tm
    T = np.linspace(SPEC[name][0], SPEC[name][1], 80) * Tm
    tm = np.logspace(-6, -2, 100)
    mp = d["creep"](diffusion=True).mechanism_map(T, tm)
    F = mp["field"].T.astype(float); F = np.where(F == 3, 2, F); F = np.where(F == 2, 2, F)
    ax.pcolormesh(T / Tm, tm, np.ma.masked_less(np.where(mp["field"].T == 3, 2, np.where(mp["field"].T == 4, 4, mp["field"].T)), 0),
                  cmap=cmap, norm=norm, shading="auto", rasterized=True)
    lo, hi = mp["dsa_band"]
    ax.fill_between(T / Tm, lo / 1.15, hi * 1.15, where=np.isfinite(lo), facecolor="#1f4e9a", alpha=0.45, lw=0)
    lg = np.log10(np.clip(mp["rate"].T, 1e-300, None))
    cs = ax.contour(T / Tm, tm, lg, levels=np.arange(-12, 3, 2), colors="k", linewidths=0.5)
    ax.clabel(cs, fmt=lambda v: fr"$10^{{{int(v)}}}$", fontsize=5)
    P = [p for p in allp if p[0] == name]
    Tp = np.array([p[1] for p in P]); tp = np.array([p[3] for p in P]) * 1e6
    ax.scatter(Tp / Tm, tp / d["mat"].mu(Tp), s=8, c="k", zorder=5, label="tensile")
    if name == "Eurofer97":
        Cc = np.vstack([EU_CREEP, EU_CREEP_YU])
        ax.scatter(Cc[:, 0] / Tm, Cc[:, 1] * 1e6 * EU_M / d["mat"].mu(Cc[:, 0]), s=8, marker="^",
                   facecolors="none", edgecolors="#c0392b", lw=0.7, zorder=5, label="creep")
        ax.legend(loc="lower left", fontsize=5.5)
    ax.set(yscale="log", xlim=(T[0] / Tm, T[-1] / Tm), ylim=(1e-6, 1e-2), xlabel=r"$T/T_m$")
    if k == 0:
        ax.set_ylabel(r"$\tau/\mu$")
        ax.legend(handles=[Patch(color=FIELD_COL[0], label="kink pair"), Patch(color=FIELD_COL[1], label="obstacle"),
                           Patch(color="#1f4e9a", alpha=0.45, label="DSA ($m<0$)"),
                           Patch(color=FIELD_COL[2], label="climb/recovery"), Patch(color=FIELD_COL[3], label="diffusional")],
                  loc="lower left", fontsize=5)
    lab(ax, f"({'abc'[k]}) {name}")
# parity (stress)
ax = fig.add_subplot(gs[1, 0])
cols = {"W": "k", "Zr": "#2471a3", "Eurofer97": "#239b56"}
for name in MATERIALS:
    P = [p for p in allp if p[0] == name]
    T = np.array([p[1] for p in P]); e = np.array([p[2] for p in P]); t = np.array([p[3] for p in P])
    tm_ = MATERIALS[name]["make"]().solve(e / MATERIALS[name]["M"], T, details=False)["tau"] / 1e6
    ax.scatter(t, tm_, s=9, facecolors="none", edgecolors=cols[name], lw=0.7, label=name)
ax.plot([20, 700], [20, 700], "k-", lw=0.6)
ax.fill_between([20, 700], [18, 630], [22, 770], color="0.88", zorder=0)
ax.set(xscale="log", yscale="log", xlim=(25, 700), ylim=(25, 700), xlabel=r"measured $\tau$ (MPa)",
       ylabel=r"model $\tau$ (MPa)")
ax.legend(loc="lower right"); lab(ax, "(d)")
# creep Norton
ax = fig.add_subplot(gs[1, 1])
cm = creep_Eurofer()
tc = {723: "#1b4f72", 773: "#148f77", 823: "#d4ac0d", 873: "#ca6f1e", 923: "#922b21"}
s_grid = np.logspace(np.log10(40), np.log10(420), 80)
for T, c in tc.items():
    r = cm.rate(s_grid[None, :] * 1e6 * EU_M, np.array([float(T)]))["total"][0] * EU_M
    ax.loglog(s_grid, r, color=c, lw=1.1, label=f"{T - 273:.0f} °C")
    k = EU_CREEP[:, 0] == T
    ax.scatter(EU_CREEP[k, 1], EU_CREEP[k, 2], s=9, color=c, edgecolors="k", lw=0.3, zorder=5)
ax.scatter(EU_CREEP_YU[:, 1], EU_CREEP_YU[:, 2], s=14, marker="*", color=tc[823], edgecolors="k", lw=0.3, zorder=6)
ax.set(xlim=(40, 420), ylim=(1e-11, 1e-4), xlabel=r"$\sigma$ (MPa)", ylabel=r"$\dot\epsilon_{\min}$ (s$^{-1}$)")
ax.legend(loc="upper left", fontsize=5.5); lab(ax, "(e)")
# creep parity
ax = fig.add_subplot(gs[1, 2])
Cc = np.vstack([EU_CREEP, EU_CREEP_YU]); pred = np.empty(len(Cc))
for T in np.unique(Cc[:, 0]):
    k = Cc[:, 0] == T
    pred[k] = cm.rate(Cc[k, 1][None, :] * 1e6 * EU_M, np.array([T]))["total"][0] * EU_M
for T, c in tc.items():
    k = Cc[:, 0] == T
    ax.scatter(Cc[k, 2], pred[k], s=9, color=c, edgecolors="k", lw=0.3)
ax.plot([1e-11, 1e-4], [1e-11, 1e-4], "k-", lw=0.6)
ax.fill_between([1e-11, 1e-4], [1e-12, 1e-5], [1e-10, 1e-3], color="0.88", zorder=0)
ax.set(xscale="log", yscale="log", xlim=(1e-11, 1e-4), ylim=(1e-11, 1e-4),
       xlabel=r"measured $\dot\epsilon_{\min}$ (s$^{-1}$)", ylabel=r"model $\dot\epsilon_{\min}$ (s$^{-1}$)")
lab(ax, "(f)")
plt.tight_layout(pad=0.3, h_pad=0.5, w_pad=0.5)
fig.savefig(f"{OUT}/fig4_maps_creep.pdf"); fig.savefig(f"{OUT}/fig4_maps_creep.png")
plt.close(fig)

le = np.log10(np.clip(pred, 1e-20, None)) - np.log10(Cc[:, 2])
print("creep: rms log10 error = %.2f, within one decade: %d/%d" % (np.sqrt(np.mean(le**2)), np.sum(np.abs(le) < 1), len(le)))
