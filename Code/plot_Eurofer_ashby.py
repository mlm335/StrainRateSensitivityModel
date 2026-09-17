"""
plot_Eurofer_ashby.py -- Eurofer-97: deformation-mechanism (Ashby) map.

Run:  python plot_W_ashby.py
Diffusional (Nabarro-Herring + Coble) flow is included here, and only here.
Axes: tau/mu against T/Tm, with contours of log10(gdot).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import physics as ph
import parameters as pm
from style import new_figure, save

P = pm.EU
T = np.linspace(0.16, 0.60, 90) * P["Tm"]
tau_mu = np.logspace(-6, -2, 90)
mp = ph.mechanism_map(T, tau_mu, P, diffusion=True)

FIELDS = {0: ("kink-pair glide", "#e6e6e6"), 1: ("obstacle / athermal", "#d3e3f3"),
          3: ("climb, recovery", "#f4c9a8"), 4: ("diffusional flow", "#cfe8c6")}
cmap = ListedColormap([FIELDS[k][1] for k in (0, 1, 3, 4)])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 3.5, 4.5], cmap.N)

fig, ax = new_figure()
ax.pcolormesh(T / P["Tm"], tau_mu, np.ma.masked_less(mp["field"].T, 0), cmap=cmap, norm=norm,
              shading="auto", rasterized=True)
lg = np.log10(np.clip(mp["rate"].T, 1e-300, None))
cs = ax.contour(T / P["Tm"], tau_mu, lg, levels=np.arange(-12, 3, 2), colors="k", linewidths=0.6)
ax.clabel(cs, fmt=lambda v: fr"$10^{{{int(v)}}}$", fontsize=7)

# tensile and creep data (converted to shear)
for j in range(3):
    ok = np.isfinite(pm.EU_SIGMA[:, j])
    ax.scatter(pm.EU_T[ok] / P["Tm"], pm.EU_SIGMA[ok, j] * 1e6 * P["M_s"] / ph.mu(pm.EU_T[ok], P),
               s=20, c="k", zorder=5)
CR = np.vstack([pm.EU_CREEP, pm.EU_CREEP_YU])
ax.scatter(CR[:, 0] / P["Tm"], CR[:, 1] * 1e6 * P["M_s"] / ph.mu(CR[:, 0], P), s=20, marker="^",
           facecolors="none", edgecolors="#c0392b", linewidths=0.8, zorder=5)
lo, hi = mp["dsa_band"]
ax.fill_between(T / P["Tm"], lo / ph.mu(T, P), hi / ph.mu(T, P), where=np.isfinite(lo),
                facecolor="#1f4e9a", alpha=0.45, lw=0)
ax.set(yscale="log", xlabel=r"$T/T_m$", ylabel=r"$\tau/\mu$", xlim=(T[0] / P["Tm"], T[-1] / P["Tm"]),
       ylim=(1e-6, 1e-2), title="Eurofer-97")
handles = [Patch(color=c, label=n) for n, c in FIELDS.values()]
ax.legend(handles=handles + [plt.Line2D([], [], ls="", marker="o", color="k", label="tensile data"),
                             plt.Line2D([], [], ls="", marker="^", mfc="none", color="#c0392b", label="creep data")],
          loc="lower left", fontsize=7.5)
save(fig, "Eurofer_ashby_map")
