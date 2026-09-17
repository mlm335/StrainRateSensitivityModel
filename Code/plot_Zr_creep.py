"""
plot_Zr_creep.py -- zirconium: creep predictions (no diffusional flow).

Run:  python plot_W_creep.py
Stress control: the shear rate is read off tau(gdot) at the imposed shear stress.
Zr has no creep data in this set, so these curves are predictions from the tensile fit
with lattice self-diffusion for the recovery term.
"""
import numpy as np
import physics as ph
import parameters as pm
from style import new_figure, save, RATE_COLORS

P = pm.ZR
T_LIST = np.array([700.0, 800.0, 900.0, 1000.0])
tau = np.logspace(np.log10(1.0), np.log10(120.0), 80) * 1e6      # shear stress (Pa)

# ---------------------------------------------------------------- 1. creep curves
fig, ax = new_figure()
for c, T in zip(RATE_COLORS, T_LIST):
    g = ph.creep_rate(tau[None, :], np.array([T]), P, diffusion=False)["total"][0]
    ax.loglog(tau / 1e6, g, color=c, label=f"{T:.0f} K")
ax.set(xlabel=r"$\tau$ (MPa)", ylabel=r"$\dot\gamma$ (s$^{-1}$)", ylim=(1e-12, 1e-2))
ax.legend(title="Zr (prediction)", title_fontsize=9, fontsize=9)
save(fig, "Zr_creep_curves")

# ---------------------------------------------------------------- 2. stress exponent
fig, ax = new_figure()
for c, T in zip(RATE_COLORS, T_LIST):
    g = ph.creep_rate(tau[None, :], np.array([T]), P, diffusion=False)["total"][0]
    ok = np.isfinite(g) & (g > 1e-14)
    n = np.gradient(np.log(g[ok]), np.log(tau[ok]))
    ax.semilogx(tau[ok] / 1e6, n, color=c, label=f"{T:.0f} K")
ax.axhline(3, color="0.6", lw=0.9, ls="--")
ax.text(0.98, 0.06, r"natural creep law $n=3$", transform=ax.transAxes, ha="right", fontsize=8, color="0.35")
ax.set(xlabel=r"$\tau$ (MPa)", ylabel=r"$n = \partial\ln\dot\gamma/\partial\ln\tau$", ylim=(0, 25))
ax.legend(fontsize=9)
save(fig, "Zr_creep_exponent")

# ---------------------------------------------------------------- 3. iso-rate stresses
fig, ax = new_figure()
T = np.linspace(600.0, 1100.0, 60)
for c, gdot in zip(RATE_COLORS, [1e-10, 1e-8, 1e-6, 1e-4]):
    r = ph.flow_stress(np.full(T.shape, gdot), T, P)
    ax.semilogy(T, r["tau"] / 1e6, color=c, label=fr"$\dot\gamma = 10^{{{int(np.log10(gdot))}}}$ s$^{{-1}}$")
ax.set(xlabel=r"$T$ (K)", ylabel=r"$\tau$ (MPa)")
ax.legend(fontsize=9)
save(fig, "Zr_creep_isorate")
