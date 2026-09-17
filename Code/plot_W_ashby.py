"""
plot_W_ashby.py -- tungsten: deformation-mechanism (Ashby) map.

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

P = pm.W
T = np.linspace(0.05, 0.70, 90) * P["Tm"]
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

# tensile data (converted to shear)
ax.scatter(pm.W_FLOW[:, 0] / P["Tm"], pm.W_FLOW[:, 2] * 1e6 * P["M_s"] / ph.mu(pm.W_FLOW[:, 0], P),
           s=22, c="k", zorder=5, label="tensile data")
ax.set(yscale="log", xlabel=r"$T/T_m$", ylabel=r"$\tau/\mu$", xlim=(T[0] / P["Tm"], T[-1] / P["Tm"]),
       ylim=(1e-6, 1e-2), title="W")
handles = [Patch(color=c, label=n) for n, c in FIELDS.values()]
ax.legend(handles=handles + [plt.Line2D([], [], ls="", marker="o", color="k", label="tensile data")],
          loc="lower left", fontsize=7.5)
save(fig, "W_ashby_map")
