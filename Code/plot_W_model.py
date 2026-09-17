"""
plot_W_model.py -- tungsten: flow stress, rate sensitivity, activation volume.

Run:  python plot_W_model.py
Each plot is a separate figure in Code/figures/.  All stresses are shear (tau).
"""
import numpy as np
import physics as ph
import parameters as pm
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS, MARKERS

P = pm.W
T = np.linspace(300.0, 900.0, 121)

# ---------------------------------------------------------------- 1. flow stress
fig, ax = new_figure()
for c, e in zip(RATE_COLORS, P["rates"]):
    r = ph.flow_stress(np.full(T.shape, e / P["M_s"]), T, P)
    ax.plot(T, r["tau"] / 1e6, color=c, label=fr"$\dot\epsilon$ = " + rate_label(e))
    k = np.isclose(pm.W_FLOW[:, 1], e)
    ax.scatter(pm.W_FLOW[k, 0], pm.W_FLOW[k, 2] * P["M_s"], s=34, facecolors="none",
               edgecolors=c, linewidths=1.2, zorder=5)
ax.plot(T, r["tau_a"] / 1e6, "k-.", lw=1.2, label=r"$\tau_a$")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 900))
ax.legend(title="lines: model, symbols: Miranda 2026", title_fontsize=8, fontsize=8)
save(fig, "W_flow_stress")

# ---------------------------------------------------------------- 2. stress parts
fig, ax = new_figure()
r = ph.flow_stress(np.full(T.shape, P["rates"][1] / P["M_s"]), T, P)
ax.plot(T, r["tau"] / 1e6, "k", label=r"$\tau$")
ax.plot(T, r["tau_a"] / 1e6, color=DATA_COLORS[1], label=r"$\tau_a$ (forest)")
ax.plot(T, r["tau_star"] / 1e6, color=DATA_COLORS[0], label=r"$\tau^\star$ (kink pairs)")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 900))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(P["rates"][1]), title_fontsize=8, fontsize=9)
save(fig, "W_stress_parts")

# ---------------------------------------------------------------- 3. m(T)
g = ph.global_m_V(P["rates"], T, P)
fig, ax = new_figure()
ax.plot(T, g["m"], "k", label="model")
for i, (name, (Td, y)) in enumerate(pm.W_M_DATA.items()):
    ax.scatter(Td, y, s=28, marker=MARKERS[i % 8], facecolors="none",
               edgecolors=DATA_COLORS[i % 7], linewidths=1.1, label=name)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$m = \partial\ln\tau/\partial\ln\dot\gamma$",
       xlim=(300, 900), ylim=(0, 0.06))
ax.legend(fontsize=7.5)
save(fig, "W_srs_m")

# ---------------------------------------------------------------- 4. V*(T)
fig, ax = new_figure()
ax.plot(T, g["Vstar_b3"], "k", label="model")
for i, (name, (Td, y)) in enumerate(pm.W_V_DATA.items()):
    ax.scatter(Td, y, s=28, marker=MARKERS[i % 8], facecolors="none",
               edgecolors=DATA_COLORS[i % 7], linewidths=1.1, label=name)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$V^\star$ ($b^3$)", xlim=(300, 900), ylim=(0, 250))
ax.legend(fontsize=7.5)
save(fig, "W_activation_volume")

# ---------------------------------------------------------------- 5. Lambda(T)
fig, ax = new_figure()
r = ph.srs(P["rates"][1] / P["M_s"] * np.ones_like(T), T, P)
ax.plot(T, r["Lambda"], "k")
ax.axhline(0, color="0.6", lw=0.8)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\Lambda = AB$", xlim=(300, 900), ylim=(-1.5, 1.5),
       title=r"W: no rate-dependent resistance, $\Lambda \approx 0$")
save(fig, "W_lambda")
