"""
plot_W_ss.py -- tungsten: stress-strain curves from the same flow rule.

Run:  python plot_W_ss.py
Each plot is a separate figure in Code/figures/.  All stresses are shear (tau).

Nothing in physics.py or parameters.py is used differently here: SS_physics.py
integrates the same Kocks-Mecking density equation along strain instead of taking
its steady-state root, so the end of every curve is exactly the flow stress the
model plots already show.  Two inputs are needed that steady state does not fix,
and they are set just below.
"""
import numpy as np
import physics as ph
import parameters as pm
import SS_physics as ss
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS

P = pm.W

# =====================================================================
# NEW INPUTS (not needed by the steady-state model)
# =====================================================================
RHO_0 = 1e13          # initial dislocation density (m^-2); stress-relieved W.
                      # Steady state at these conditions is 5-9e13 m^-2.
EPS_SAT = 0.25        # axial strain over which hardening saturates.  Fixes the
                      # hardening scale S: k1, k2 and K are all divided by S, so
                      # the steady state -- every result already in the folder --
                      # is unchanged and only the approach to it is slowed.
EPS_MAX = 0.5         # axial strain the curves run to
N_PTS = 301           # points per curve

T_LIST = np.array([300.0, 523.0, 673.0, 823.0])      # for the temperature sweep
RATE_LIST = P["rates"]                                # axial rates, from parameters.py
T_RATE = 523.0                                        # T for the rate sweep
T_PARTS, E_PARTS = 673.0, P["rates"][1]               # condition for the decomposition

# adiabatic heating: dT/dgamma = eta tau / (rho_mass c_p)
E_FAST = 1.0e3        # axial rate at which to show adiabatic heating (s^-1)
T_FAST = 300.0
ADIABATIC = dict(eta=0.9, rho_mass=19250.0, c_p=134.0)   # W: kg/m^3, J/kg/K
# =====================================================================

S = ss.hardening_scale(P, T_LIST, EPS_SAT)
print(f"hardening scale S = {S:.1f}  (rho_0 = {RHO_0:.1e} m^-2)")

# ---------------------------------------------------------------- 1. curves vs T
fig, ax = new_figure()
r = ss.stress_strain(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, T in enumerate(T_LIST):
    ax.plot(r["gamma"][i], r["tau"][i] / 1e6, color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, EPS_MAX / P["M_s"]), ylim=(0, None))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "W_ss_temperature")

# ---------------------------------------------------------------- 2. curves vs rate
fig, ax = new_figure()
g = np.array([e / P["M_s"] for e in RATE_LIST])
r = ss.stress_strain(g, np.full(g.shape, T_RATE), P, RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, e in enumerate(RATE_LIST):
    ax.plot(r["gamma"][i], r["tau"][i] / 1e6, color=RATE_COLORS[i], label=rate_label(e))
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, EPS_MAX / P["M_s"]), ylim=(0, None))
ax.legend(title=f"$T$ = {T_RATE:.0f} K", title_fontsize=8, fontsize=9)
save(fig, "W_ss_rate")

# ---------------------------------------------------------------- 3. stress parts
fig, ax = new_figure()
r = ss.stress_strain(np.array([E_PARTS / P["M_s"]]), np.array([T_PARTS]), P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
ax.plot(r["gamma"][0], r["tau"][0] / 1e6, "k", label=r"$\tau$")
ax.plot(r["gamma"][0], r["tau_a"][0] / 1e6, color=DATA_COLORS[1], label=r"$\tau_a$ (forest)")
ax.plot(r["gamma"][0], r["tau_star"][0] / 1e6, color=DATA_COLORS[0], label=r"$\tau^\star$ (kink pairs)")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, EPS_MAX / P["M_s"]), ylim=(0, None))
ax.legend(title=f"$T$ = {T_PARTS:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_PARTS),
          title_fontsize=8, fontsize=9)
save(fig, "W_ss_stress_parts")

# ---------------------------------------------------------------- 4. density
fig, ax = new_figure()
r = ss.stress_strain(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, T in enumerate(T_LIST):
    ax.semilogy(r["gamma"][i], r["rho"][i], color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.axhline(RHO_0, color="0.6", lw=0.8, ls=":")
ax.text(0.02, RHO_0 * 1.15, r"$\rho_0$", fontsize=9, color="0.4")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\rho$ (m$^{-2}$)", xlim=(0, EPS_MAX / P["M_s"]))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "W_ss_density")

# ---------------------------------------------------------------- 5. adiabatic
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
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, EPS_MAX / P["M_s"]), ylim=(0, None))
ax.legend(title=f"$T_0$ = {T_FAST:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_FAST),
          title_fontsize=8, fontsize=9, loc="lower right")
save(fig, "W_ss_adiabatic")

# ---------------------------------------------------------------- 6. Voce check
# With Lambda_b = 0, K_cl = 0 and no aging, the same equation has the closed-form
# Voce solution.  The integrator must reproduce it.
fig, ax = new_figure()
Q = dict(P, K_cl=0.0, Lambda_b=0.0, C_sol=0.0)
r = ss.stress_strain(np.array([RATE_LIST[1] / P["M_s"]]), np.array([T_PARTS]), Q,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
va = ss.voce_analytic(r["gamma"][0], RATE_LIST[1] / P["M_s"], T_PARTS, Q, RHO_0, S=S)
ax.plot(r["gamma"][0], r["tau_a"][0] / 1e6, "k", label=r"integrated $\tau_a$")
ax.plot(r["gamma"][0], va / 1e6, color=DATA_COLORS[0], ls="--", label="analytic Voce")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau_a$ (MPa)", xlim=(0, EPS_MAX / P["M_s"]), ylim=(0, None))
ax.legend(title=r"$\Lambda_b=0$, $K=0$", title_fontsize=8, fontsize=9)
save(fig, "W_ss_voce_check")
print(f"  Voce check: max rel. error {np.max(np.abs(r['tau_a'][0] - va) / va):.2e}")

# ---------------------------------------------------------------- verification
d = ss.check_saturation(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P, RHO_0, S=S)
print(f"  saturation vs physics.flow_stress: max rel. difference {d:.2e}")
