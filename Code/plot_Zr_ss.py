"""
plot_Zr_ss.py -- zirconium (prismatic slip): stress-strain curves.

Run:  python plot_Zr_ss.py
Each plot is a separate figure in Code/figures/.  All stresses are shear (tau).

SS_physics.py integrates the same density equation along strain instead of taking
its steady-state root, so each curve ends on the flow stress the model plots show.
Zr adds oxygen aging and climb-controlled recovery to the tungsten-level model, so
the curves also show the atmospheres building up as the material deforms.
"""
import numpy as np
import physics as ph
import parameters as pm
import SS_physics as ss
from style import new_figure, save, rate_label, RATE_COLORS, DATA_COLORS

P = pm.ZR

# =====================================================================
# NEW INPUTS (not needed by the steady-state model)
# =====================================================================
RHO_0 = 1e13          # initial dislocation density (m^-2); annealed Zr.
                      # Steady state at these conditions is 1-2e14 m^-2.
EPS_SAT = 0.25        # axial strain over which hardening saturates -> sets the
                      # hardening scale S (k1, k2 and K all divided by S, which
                      # leaves the steady state, and every existing result, alone)
EPS_MAX = 0.4         # axial strain the curves run to
N_PTS = 301           # points per curve

T_LIST = np.array([300.0, 500.0, 650.0, 800.0])       # for the temperature sweep
RATE_LIST = P["rates"]                                 # axial rates, from parameters.py
T_RATE = 500.0                                         # T for the rate sweep
T_PARTS, E_PARTS = 500.0, P["rates"][1]                # condition for the decomposition

# rate-jump test: deform at the base rate, jump up, then back down
T_JUMP = [500.0, 665.0]        # below the aging window, and inside it
JUMP_GAMMA = [0.6, 1.2]        # shear strains at which the rate changes
JUMP_RATES = [P["rates"][1], P["rates"][2], P["rates"][1]]   # axial rates, in order
EPS_MAX_JUMP = 0.8             # axial strain for the jump test
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
save(fig, "Zr_ss_temperature")

# ---------------------------------------------------------------- 2. curves vs rate
fig, ax = new_figure()
g = np.array([e / P["M_s"] for e in RATE_LIST])
r = ss.stress_strain(g, np.full(g.shape, T_RATE), P, RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, e in enumerate(RATE_LIST):
    ax.plot(r["gamma"][i], r["tau"][i] / 1e6, color=RATE_COLORS[i], label=rate_label(e))
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T$ = {T_RATE:.0f} K", title_fontsize=8, fontsize=9)
save(fig, "Zr_ss_rate")

# ---------------------------------------------------------------- 3. stress parts
fig, ax = new_figure()
r = ss.stress_strain(np.array([E_PARTS / P["M_s"]]), np.array([T_PARTS]), P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
ax.plot(r["gamma"][0], r["tau"][0] / 1e6, "k", label=r"$\tau$")
ax.plot(r["gamma"][0], r["tau_a"][0] / 1e6, color=DATA_COLORS[1], label=r"$\tau_a$ (forest + aging)")
ax.plot(r["gamma"][0], r["tau_star"][0] / 1e6, color=DATA_COLORS[0], label=r"$\tau^\star$ (kink pairs)")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau$ (MPa)", xlim=(0, GMAX), ylim=(0, None))
ax.legend(title=f"$T$ = {T_PARTS:.0f} K, " + r"$\dot\epsilon$ = " + rate_label(E_PARTS),
          title_fontsize=8, fontsize=9)
save(fig, "Zr_ss_stress_parts")

# ---------------------------------------------------------------- 4. density
fig, ax = new_figure()
r = ss.stress_strain(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P,
                     RHO_0, S=S, eps_max=EPS_MAX, n=N_PTS)
for i, T in enumerate(T_LIST):
    ax.semilogy(r["gamma"][i], r["rho"][i], color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.axhline(RHO_0, color="0.6", lw=0.8, ls=":")
ax.text(0.02, RHO_0 * 1.15, r"$\rho_0$", fontsize=9, color="0.4")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\rho$ (m$^{-2}$)", xlim=(0, GMAX))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "Zr_ss_density")

# ---------------------------------------------------------------- 5. aging along the curve
# The waiting time t_a = rho_m b lambda_f / gdot grows as the forest builds up, so
# the oxygen atmospheres thicken as the sample is strained.
fig, ax = new_figure()
for i, T in enumerate(T_LIST):
    ax.plot(r["gamma"][i], r["C_C0"][i], color=RATE_COLORS[i], label=f"{T:.0f} K")
ax.set(xlabel=r"$\gamma$", ylabel=r"$C/C_0$ (oxygen atmospheres)", xlim=(0, GMAX))
ax.legend(title=r"$\dot\epsilon$ = " + rate_label(RATE_LIST[1]), title_fontsize=8, fontsize=9)
save(fig, "Zr_ss_aging")

# ---------------------------------------------------------------- 6. rate-jump test
# The clearest signature of m < 0.  Outside the aging window a rate increase
# raises the stress; inside it, slower flow has had more time to build atmospheres,
# so the same increase LOWERS the stress.  Normalised by the stress just before
# the first jump, because the two temperatures are far apart in absolute stress.
fig, ax = new_figure()
jump = ss.rate_steps(JUMP_GAMMA, [e / P["M_s"] for e in JUMP_RATES])
for i, T in enumerate(T_JUMP):
    q = ss.stress_strain(jump, np.array([T]), P, RHO_0, S=S, eps_max=EPS_MAX_JUMP, n=801)
    gq, tq = q["gamma"][0], q["tau"][0]
    ref = tq[np.searchsorted(gq, JUMP_GAMMA[0]) - 1]
    ax.plot(gq, tq / ref, color=DATA_COLORS[i],
            label=f"{T:.0f} K" + ("  (aging window)" if i else ""))
    print("  rate jump at %.0f K: %+.2f%% on the increase, %+.2f%% on the decrease"
          % (T, 100 * (tq[np.searchsorted(gq, JUMP_GAMMA[0]) + 1] / ref - 1),
             100 * (tq[np.searchsorted(gq, JUMP_GAMMA[1]) + 1]
                    / tq[np.searchsorted(gq, JUMP_GAMMA[1]) - 1] - 1)))
for gj in JUMP_GAMMA:
    ax.axvline(gj, color="0.75", lw=0.8, ls=":")
ax.set(xlabel=r"$\gamma$", ylabel=r"$\tau/\tau_{\rm ref}$", xlim=(0, EPS_MAX_JUMP / P["M_s"]))
ax.legend(title=r"$\dot\epsilon$: " + " $\\to$ ".join(rate_label(e) for e in JUMP_RATES),
          title_fontsize=7.5, fontsize=9)
save(fig, "Zr_ss_rate_jump")

# ---------------------------------------------------------------- verification
d = ss.check_saturation(np.full(T_LIST.shape, RATE_LIST[1] / P["M_s"]), T_LIST, P, RHO_0, S=S)
print(f"  saturation vs physics.flow_stress: max rel. difference {d:.2e}")
