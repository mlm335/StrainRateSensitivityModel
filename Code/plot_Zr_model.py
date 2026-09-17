"""
plot_Zr_model.py -- zirconium (prismatic slip): flow stress, m, V*, DSA diagnostics.

Run:  python plot_Zr_model.py
Zr adds oxygen aging and climb-controlled recovery to the tungsten-level model.
"""
import numpy as np
import physics as ph
import parameters as pm
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS, MARKERS

P = pm.ZR
T = np.linspace(300.0, 1000.0, 141)

# ---------------------------------------------------------------- 1. flow stress
fig, ax = new_figure()
for c, e in zip([RATE_COLORS[0], RATE_COLORS[2]], sorted(pm.ZR_CRSS)):
    r = ph.flow_stress(np.full(T.shape, e / P["M_s"]), T, P)
    ax.plot(T, r["tau"] / 1e6, color=c, label=r"$\dot\epsilon$ = " + rate_label(e))
    Td, y = pm.ZR_CRSS[e]
    ax.scatter(Td, y, s=34, facecolors="none", edgecolors=c, linewidths=1.2, zorder=5)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 1000), ylim=(0, 150))
ax.legend(title="lines: model, symbols: Derep et al.", title_fontsize=8, fontsize=9)
save(fig, "Zr_flow_stress")

# ---------------------------------------------------------------- 2. what each term does
fig, ax = new_figure()
e = 3.3e-5
full = ph.flow_stress(np.full(T.shape, e / P["M_s"]), T, P)
no_climb = dict(P, K_cl=0.0)
base = dict(P, C_sol=0.0, K_cl=0.0)
ax.plot(T, ph.flow_stress(np.full(T.shape, e / P["M_s"]), T, base)["tau"] / 1e6,
        color="0.55", ls="--", label="kink pairs + forest only")
ax.plot(T, ph.flow_stress(np.full(T.shape, e / P["M_s"]), T, no_climb)["tau"] / 1e6,
        color=DATA_COLORS[1], ls=":", label="+ oxygen aging")
ax.plot(T, full["tau"] / 1e6, color=DATA_COLORS[0], label="+ climb recovery (full)")
Td, y = pm.ZR_CRSS[3.3e-5]
ax.scatter(Td, y, s=34, facecolors="none", edgecolors="k", linewidths=1.2, zorder=5)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)", xlim=(300, 1000), ylim=(0, 150))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(e), title_fontsize=8, fontsize=9)
save(fig, "Zr_model_build_up")

# ---------------------------------------------------------------- 3. m(T)
g = ph.global_m_V(P["rates"], T, P)
fig, ax = new_figure()
ax.plot(T, g["m"], "k", label="model")
ax.axhline(0, color="0.6", lw=0.8)
for i, (name, (Td, y)) in enumerate(pm.ZR_M_DATA.items()):
    ax.scatter(Td, y, s=28, marker=MARKERS[i % 8], facecolors="none",
               edgecolors=DATA_COLORS[i % 7], linewidths=1.1, label=name)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$m$", xlim=(300, 1000), ylim=(-0.02, 0.10))
ax.legend(fontsize=7)
save(fig, "Zr_srs_m")

# ---------------------------------------------------------------- 4. V*(T)
fig, ax = new_figure()
ax.plot(T, g["Vstar_b3"], "k", label="model")
for i, (name, (Td, y)) in enumerate(pm.ZR_V_DATA.items()):
    ax.scatter(Td, y, s=28, marker=MARKERS[i % 8], facecolors="none",
               edgecolors=DATA_COLORS[i % 7], linewidths=1.1, label=name)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$V^\star$ ($b^3$)", xlim=(300, 1000), ylim=(0, 300))
ax.legend(fontsize=7.5)
save(fig, "Zr_activation_volume")

# ---------------------------------------------------------------- 5. Lambda(T)
fig, ax = new_figure()
r = ph.srs(np.full(T.shape, 1e-4 / P["M_s"]), T, P)
ax.plot(T, r["Lambda"], "k")
ax.set_yscale("symlog", linthresh=0.1)
ax.axhspan(-1e4, -1, color="#d6e6f5"); ax.axhspan(-1, 1, color="#eef5ee"); ax.axhspan(1, 1e4, color="#f8e1e1")
ax.text(0.98, 0.06, r"$m<0$ (DSA)", transform=ax.transAxes, ha="right", fontsize=8)
ax.text(0.98, 0.93, r"resistance ($B$) controlled", transform=ax.transAxes, ha="right", va="top", fontsize=8)
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\Lambda = AB$", xlim=(300, 1000), ylim=(-30, 1e4))
save(fig, "Zr_lambda")

# ---------------------------------------------------------------- 6. aging state
fig, ax = new_figure()
for c, e in zip(RATE_COLORS, P["rates"]):
    S = ph.state(np.full(T.shape, e / P["M_s"]), T, P)
    ax.plot(T, S["C_C0"], color=c, label=r"$\dot\epsilon$ = " + rate_label(e))
ax.set(xlabel=r"$T$ (K)", ylabel=r"$C/C_0$ (oxygen atmospheres)", xlim=(300, 1000))
ax.legend(fontsize=9)
save(fig, "Zr_aging")
