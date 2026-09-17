"""
plot_Eurofer_ss.py -- Eurofer-97: stress-strain curves.

Run:  python plot_Eurofer_ss.py
Each plot is a separate figure in Code/figures/.  All stresses are shear (tau).

SS_physics.py integrates the same density equation along strain instead of taking
its steady-state root, so each curve ends on the flow stress the model plots show.
Eurofer-97 carries most of its strength in the boundary term and in the particles,
neither of which changes much with strain, so the curves are nearly flat -- and a
low initial mobile density produces a yield drop, since tau* falls as rho_m grows.
"""
import numpy as np
import physics as ph
import parameters as pm
import SS_physics as ss
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS

P = pm.EU

# =====================================================================
# NEW INPUTS (not needed by the steady-state model)
# =====================================================================
RHO_0 = 1e14          # initial dislocation density (m^-2); tempered martensite
                      # already starts dense.  Steady state here is 2-3e14 m^-2.
EPS_SAT = 0.15        # axial strain over which hardening saturates -> sets the
                      # hardening scale S (k1, k2 and K all divided by S, which
                      # leaves the steady state, and every existing result, alone)
EPS_MAX = 0.25        # axial strain the curves run to (Eurofer necks early)
N_PTS = 301           # points per curve

T_LIST = np.array([300.0, 573.0, 723.0, 873.0])       # for the temperature sweep
RATE_LIST = P["rates"]                                 # axial rates, from parameters.py
T_RATE = 573.0                                         # T for the rate sweep
T_PARTS, E_PARTS = 573.0, P["rates"][1]                # condition for the decomposition
RHO_0_LIST = [1e12, 1e13, 1e14, 1e15]                  # yield-drop sweep
T_DROP = 300.0                                         # T for the yield-drop sweep
                                                       # (tau* only matters when cold)

# adiabatic heating: dT/dgamma = eta tau / (rho_mass c_p)
E_FAST = 1.0e2        # axial rate at which to show adiabatic heating (s^-1)
T_FAST = 300.0
ADIABATIC = dict(eta=0.9, rho_mass=7750.0, c_p=490.0)  # Eurofer-97: kg/m^3, J/kg/K
# =====================================================================

S = ss.hardening_scale(P, T_LIST, EPS_SAT)
GMAX = EPS_MAX / P["M_s"]
print(f"hardening scale S = {S:.1f}  (rho_0 = {RHO_0:.1e} m^-2)")

# ---------------------------------------------------------------- 1. curves vs T
fig, ax = new_figure()
r = ss.stress_strain(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, T in enumerate(T_LIST):
    ax.plot(r["gamma"][i], r["tau"][i] / 1e6, color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "Eurofer_ss_temperature")

# ---------------------------------------------------------------- 2. curves vs rate
fig, ax = new_figure()
g = np.array([e / P["M_s"] for e in RATE_LIST])
r = ss.stress_strain(g, np.full(g.shape, T_RATE), P, RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, e in enumerate(RATE_LIST):
    ax.plot(r["gamma"][i], r["tau"][i] / 1e6, color=RATE_COLORS[i], label=rate_label(e))
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T$ = {T_RATE:.0f} K", title_fontsize=8, fontsize=9)
save(fig, "Eurofer_ss_rate")

# ---------------------------------------------------------------- 3. stress parts
fig, ax = new_figure()
r = ss.stress_strain(np.array([E_PARTS / P["M_s"]]), np.array([T_PARTS]), P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
ax.plot(r["gamma"][0], r["tau"][0] / 1e6, "k", label=r"$\tau$ (total)")
ax.plot(r["gamma"][0], r["tau_a"][0] / 1e6, color=DATA_COLORS[1],
        label=r"$\tau_a$ (forest + boundaries + aging)")
ax.plot(r["gamma"][0], r["tau_p"][0] / 1e6, color=DATA_COLORS[2],
        label=r"$\tau_p$ (MX, M$_{23}$C$_6$)")
ax.plot(r["gamma"][0], r["tau_star"][0] / 1e6, color=DATA_COLORS[0],
        label=r"$\tau^\star$ (kink pairs)")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T$ = {T_PARTS:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_PARTS),
          title_fontsize=8, fontsize=8.5)
save(fig, "Eurofer_ss_stress_parts")

# ---------------------------------------------------------------- 4. density and aging
fig, ax = new_figure()
r = ss.stress_strain(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, T in enumerate(T_LIST):
    ax.semilogy(r["gamma"][i], r["rho"][i], color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.axhline(RHO_0, color="0.6", lw=0.8, ls=":")
ax.text(0.02, RHO_0 * 1.1, r"$\rho_0$", fontsize=9, color="0.4")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\rho$ (m$^{-2}$)", xlim=(0, GMAX))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "Eurofer_ss_density")

fig, ax = new_figure()
for i, T in enumerate(T_LIST):
    ax.plot(r["gamma"][i], r["C_C0"][i], color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.set(xlabel=r"$\gamma$", ylabel=r"$C/C_0$ (C, N atmospheres)", xlim=(0, GMAX))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "Eurofer_ss_aging")

# ---------------------------------------------------------------- 5. yield drop
# A low initial mobile density forces a high tau* to carry the imposed rate; as
# rho_m multiplies, tau* falls faster than tau_a rises and the curve drops.
fig, ax = new_figure()
for i, r0 in enumerate(RHO_0_LIST):
    q = ss.stress_strain(np.array([E_PARTS / P["M_s"]]), np.array([T_DROP]), P,
                         r0, S=S, eps_max=EPS_MAX, n=N_PTS)
    ax.plot(q["gamma"][0], q["tau"][0] / 1e6, color=RATE_COLORS[i],
            label=r"$\rho_0$ = " + f"{r0:.0e}".replace("e+", r"$\times10^{") + r"}$ m$^{-2}$")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T$ = {T_DROP:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_PARTS),
          title_fontsize=8, fontsize=8.5)
save(fig, "Eurofer_ss_yield_drop")

# ---------------------------------------------------------------- 6. adiabatic
fig, ax = new_figure()
gf = np.array([E_FAST / P["M_s"]])
Tf = np.array([T_FAST])
iso = ss.stress_strain(gf, Tf, P, RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
adi = ss.stress_strain(gf, Tf, P, RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS, adiabatic=ADIABATIC)
ax.plot(iso["gamma"][0], iso["tau"][0] / 1e6, color=DATA_COLORS[1], label="isothermal")
ax.plot(adi["gamma"][0], adi["tau"][0] / 1e6, color=DATA_COLORS[0], label="adiabatic")
ax2 = ax.twinx()
ax2.plot(adi["gamma"][0], adi["T"][0] - T_FAST, color="0.55", lw=1.1, ls="--")
ax2.set_ylabel(r"$\Delta T$ (K)", color="0.4")
ax2.tick_params(axis="y", colors="0.4")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T_0$ = {T_FAST:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_FAST),
          title_fontsize=8, fontsize=9, loc="lower right")
save(fig, "Eurofer_ss_adiabatic")

# ---------------------------------------------------------------- verification
d = ss.check_saturation(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P, RHO_0, S=S)
print(f"  saturation vs physics.flow_stress: max rel. difference {d:.2e}")
