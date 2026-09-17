"""
plot_Eurofer_creep.py -- Eurofer-97: minimum creep rate against the data.

Run:  python plot_Eurofer_creep.py
Dislocation creep only (no diffusional flow).  Data: Fernandez et al. (2005),
53 points, 450-650 C, plus the 550 C / 300 MPa point of Yu et al. (2005).
Axial data are converted to shear:  tau = M_s sigma,  gdot = eps_dot / M_s.
"""
import numpy as np
import physics as ph
import parameters as pm
from style import new_figure, save, RATE_COLORS

P = pm.EU
M = P["M_s"]
DATA = np.vstack([pm.EU_CREEP, pm.EU_CREEP_YU])
T_LIST = np.unique(DATA[:, 0])
COL = dict(zip(T_LIST, RATE_COLORS))
tau_grid = np.logspace(np.log10(15.0), np.log10(140.0), 80) * 1e6

# model at the measured points
pred = np.empty(len(DATA))
for T in T_LIST:
    k = DATA[:, 0] == T
    pred[k] = ph.creep_rate(DATA[k, 1][None, :] * 1e6 * M, np.array([T]), P)["total"][0]
err = np.log10(pred) - np.log10(DATA[:, 2] / M)

# ---------------------------------------------------------------- 1. creep curves
fig, ax = new_figure()
for T in T_LIST:
    g = ph.creep_rate(tau_grid[None, :], np.array([T]), P)["total"][0]
    ax.loglog(tau_grid / 1e6, g, color=COL[T], label=f"{T-273:.0f} °C")
    k = DATA[:, 0] == T
    ax.scatter(DATA[k, 1] * M, DATA[k, 2] / M, s=26, color=COL[T], edgecolors="k", linewidths=0.4, zorder=5)
ax.set(xlabel=r"$\tau$ (MPa)", ylabel=r"$\dot\gamma_{\min}$ (s$^{-1}$)", xlim=(15, 140), ylim=(1e-11, 1e-4))
ax.legend(title="lines: model, symbols: data", title_fontsize=8, fontsize=8.5)
save(fig, "Eurofer_creep_curves")

# ---------------------------------------------------------------- 2. parity
fig, ax = new_figure()
for T in T_LIST:
    k = DATA[:, 0] == T
    ax.loglog(DATA[k, 2] / M, pred[k], "o", color=COL[T], mec="k", mew=0.4, ms=5, label=f"{T-273:.0f} °C")
lim = [1e-11, 1e-4]
ax.plot(lim, lim, "k-", lw=0.9)
ax.fill_between(lim, [x / 10 for x in lim], [x * 10 for x in lim], color="0.9", zorder=0)
ax.set(xlabel=r"measured $\dot\gamma_{\min}$ (s$^{-1}$)", ylabel=r"model $\dot\gamma_{\min}$ (s$^{-1}$)",
       xlim=lim, ylim=lim,
       title=f"rms {np.sqrt(np.mean(err**2)):.2f} decades, {np.sum(np.abs(err) < 1)}/{len(err)} within one decade")
ax.legend(fontsize=8)
save(fig, "Eurofer_creep_parity")

# ---------------------------------------------------------------- 3. stress exponent
fig, ax = new_figure()
n_model, n_data = [], []
for T in T_LIST:
    k = DATA[:, 0] == T
    g = ph.creep_rate(DATA[k, 1][None, :] * 1e6 * M, np.array([T]), P)["total"][0]
    n_model.append(np.polyfit(np.log(DATA[k, 1]), np.log(g), 1)[0])
    n_data.append(np.polyfit(np.log(DATA[k, 1]), np.log(DATA[k, 2]), 1)[0])
ax.plot(T_LIST - 273, n_data, "o-", color="k", label="data")
ax.plot(T_LIST - 273, n_model, "s--", color=RATE_COLORS[0], label="model")
ax.set(xlabel=r"$T$ (°C)", ylabel=r"$n = \partial\ln\dot\gamma/\partial\ln\tau$")
ax.legend(fontsize=9)
save(fig, "Eurofer_creep_exponent")

# ---------------------------------------------------------------- 4. role of local climb
fig, ax = new_figure()
T = 823.0
no_climb = dict(P, D0_climb=0.0)
ax.loglog(tau_grid / 1e6, ph.creep_rate(tau_grid[None, :], np.array([T]), P)["total"][0],
          color=RATE_COLORS[2], label="detachment + local climb")
ax.loglog(tau_grid / 1e6, ph.creep_rate(tau_grid[None, :], np.array([T]), no_climb)["total"][0],
          color=RATE_COLORS[2], ls="--", label="detachment only")
k = DATA[:, 0] == T
ax.scatter(DATA[k, 1] * M, DATA[k, 2] / M, s=26, color=RATE_COLORS[2], edgecolors="k", linewidths=0.4, zorder=5)
ax.set(xlabel=r"$\tau$ (MPa)", ylabel=r"$\dot\gamma_{\min}$ (s$^{-1}$)", xlim=(40, 120), ylim=(1e-11, 1e-4),
       title="550 °C")
ax.legend(fontsize=9)
save(fig, "Eurofer_creep_climb_route")

print(f"creep: rms {np.sqrt(np.mean(err**2)):.2f} decades, "
      f"{np.sum(np.abs(err) < 1)}/{len(err)} within one decade")
