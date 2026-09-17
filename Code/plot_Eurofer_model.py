"""
plot_Eurofer_model.py -- Eurofer-97: flow stress, m, V*, DSA and stress parts.

Run:  python plot_Eurofer_model.py
Eurofer adds boundary storage and MX / M23C6 particles (thermal detachment or
local climb) to the Zr-level model.  Axial data are converted to shear with M_s.
"""
import numpy as np
import physics as ph
import parameters as pm
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS, MARKERS

P = pm.EU
T = np.linspace(300.0, 900.0, 121)
M = P["M_s"]

# ---------------------------------------------------------------- 1. flow stress
fig, ax = new_figure()
for j, (c, e) in enumerate(zip(RATE_COLORS, P["rates"])):
    r = ph.flow_stress(np.full(T.shape, e / M), T, P)
    ax.plot(T, r["tau"] / 1e6, color=c, label=r"$\dot\epsilon$ = " + rate_label(e))
    y = pm.EU_SIGMA[:, j] * M
    ax.scatter(pm.EU_T, y, s=32, marker=MARKERS[j], facecolors="none", edgecolors=c,
               linewidths=1.1, zorder=5)
ax.scatter(pm.EU_MM[:, 0], pm.EU_MM[:, 1] * M, s=34, marker="v", facecolors="none",
           edgecolors="k", linewidths=1.0, label="Materna-Morris (rate n/a)")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 900))
ax.legend(title="lines: model, symbols: Vanaja et al.", title_fontsize=8, fontsize=8)
save(fig, "Eurofer_flow_stress")

# ---------------------------------------------------------------- 2. stress parts
fig, ax = new_figure()
r = ph.flow_stress(np.full(T.shape, P["rates"][1] / M), T, P)
ax.plot(T, r["tau"] / 1e6, "k", label=r"$\tau$ (total)")
ax.plot(T, r["tau_a"] / 1e6, color=DATA_COLORS[1], label=r"$\tau_a$ (forest + boundaries + aging)")
ax.plot(T, r["tau_p"] / 1e6, color=DATA_COLORS[2], label=r"$\tau_p$ (MX, M$_{23}$C$_6$)")
ax.plot(T, r["tau_star"] / 1e6, color=DATA_COLORS[0], label=r"$\tau^\star$ (kink pairs)")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 900))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(P["rates"][1]), title_fontsize=8, fontsize=8.5)
save(fig, "Eurofer_stress_parts")

# ---------------------------------------------------------------- 3. build-up
fig, ax = new_figure()
e = P["rates"][1]
level_I = dict(P, C_sol=0.0, K_cl=0.0, tau_ppt=0.0, Lambda_b=0.0)
level_II = dict(P, tau_ppt=0.0, Lambda_b=0.0)
level_II_b = dict(P, tau_ppt=0.0)
for Pi, lab, c, ls in [(level_I, "I: kink pairs + forest", "0.55", "--"),
                       (level_II, "II: + aging + climb recovery", DATA_COLORS[1], ":"),
                       (level_II_b, "+ boundaries", DATA_COLORS[3], "-."),
                       (P, "III: + particles (full)", DATA_COLORS[0], "-")]:
    ax.plot(T, ph.flow_stress(np.full(T.shape, e / M), T, Pi)["tau"] / 1e6, color=c, ls=ls, label=lab)
ax.scatter(pm.EU_T, pm.EU_SIGMA[:, 1] * M, s=32, facecolors="none", edgecolors="k", linewidths=1.1, zorder=5)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 900))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(e), title_fontsize=8, fontsize=8)
save(fig, "Eurofer_model_build_up")

# ---------------------------------------------------------------- 4. m(T)
g = ph.global_m_V(P["rates"], T, P)
fig, ax = new_figure()
ax.plot(T, g["m"], "k", label="model")
ax.axhline(0, color="0.6", lw=0.8)
ax.scatter(pm.EU_M_DATA[:, 0], pm.EU_M_DATA[:, 1], s=30, facecolors="none",
           edgecolors=DATA_COLORS[0], linewidths=1.2, label="Vanaja et al.")
ax.scatter(pm.SHAH[:, 0], pm.SHAH[:, 1], s=30, marker="v", facecolors="none",
           edgecolors=DATA_COLORS[4], linewidths=1.2, label="Shah et al. (IN-RAFM)")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$m$", xlim=(300, 900), ylim=(-0.01, 0.12))
ax.legend(fontsize=8.5)
save(fig, "Eurofer_srs_m")

# ---------------------------------------------------------------- 5. V*(T)
fig, ax = new_figure()
ax.plot(T, g["Vstar_b3"], "k", label="model")
ax.scatter(pm.EU_V_DATA[:, 0], pm.EU_V_DATA[:, 1], s=30, facecolors="none",
           edgecolors=DATA_COLORS[0], linewidths=1.2, label="Vanaja et al.")
ax.scatter(pm.SHAH[:, 0], pm.SHAH[:, 2], s=30, marker="v", facecolors="none",
           edgecolors=DATA_COLORS[4], linewidths=1.2, label="Shah et al. (IN-RAFM)")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$V^\star$ ($b^3$)", xlim=(300, 900), ylim=(0, 350))
ax.legend(fontsize=8.5)
save(fig, "Eurofer_activation_volume")

# ---------------------------------------------------------------- 6. Lambda and aging
fig, ax = new_figure()
r = ph.srs(np.full(T.shape, P["rates"][1] / M), T, P)
ax.plot(T, r["Lambda"], "k")
ax.set_yscale("symlog", linthresh=0.1)
ax.axhspan(-1e4, -1, color="#d6e6f5"); ax.axhspan(-1, 1, color="#eef5ee"); ax.axhspan(1, 1e4, color="#f8e1e1")
ax.text(0.98, 0.06, r"$m<0$ (DSA)", transform=ax.transAxes, ha="right", fontsize=8)
ax.text(0.98, 0.93, r"resistance ($B$) controlled", transform=ax.transAxes, ha="right", va="top", fontsize=8)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\Lambda = AB$", xlim=(300, 900), ylim=(-30, 1e4))
save(fig, "Eurofer_lambda")

fig, ax = new_figure()
for c, e in zip(RATE_COLORS, P["rates"]):
    S = ph.state(np.full(T.shape, e / M), T, P)
    ax.plot(T, S["C_C0"], color=c, label=r"$\dot\epsilon$ = " + rate_label(e))
ax.set(xlabel=r"$T$ (K)", ylabel=r"$C/C_0$ (C, N atmospheres)", xlim=(300, 900))
ax.legend(fontsize=9)
save(fig, "Eurofer_aging")
