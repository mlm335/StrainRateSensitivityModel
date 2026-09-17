"""
SS_physics.py -- stress-strain curves from the same flow rule.

No new physics.  physics.py solves the Kocks-Mecking density equation for its
steady state, which is the saturation stress at large strain.  This module
integrates the identical equation along strain instead,

    dx/dgamma = [ k1 (x + Lambda_b) - k2 x^2 - K D(T) x^4 / gdot ] / (2 x S),
    x = sqrt(rho),   x(0) = sqrt(rho_0),

and evaluates the same three stress terms at every point of the curve.

Two new inputs, both set at the top of the plotting scripts:

  rho_0   initial dislocation density (m^-2).

  S       hardening scale.  The tensile, creep and map results calibrate only the
          *steady state* of the density equation, i.e. the root of its right-hand
          side.  Multiplying k1, k2 and K by the same factor 1/S leaves that root
          -- and therefore every result already in the folder -- untouched, while
          making the approach to it S times slower.  One number is therefore free
          and is fixed here by the strain over which hardening saturates.  Use
          hardening_scale() to get S from a chosen saturation strain.

Optional adiabatic heating adds  dT/dgamma = eta tau / (rho_mass c_p).

Everything is in SHEAR (tau, gamma) and SI units, as in physics.py; axial values
follow from sigma = tau / M_s and eps = M_s gamma.

Functions
---------
dx_dgamma         right-hand side of the density equation, in x = sqrt(rho)
hardening_scale   S from a chosen saturation strain
state_from_rho    rho_f, rho_m, lambda_f, C/C0 and tau_a for a given rho
tau_star_of       kink-pair stress for a given rho_m
stress_from_rho   all stress terms for a given rho
stress_strain     the curve: integrate rho along gamma, evaluate the stresses
voce_analytic     exact Voce solution of the same equation (Lambda_b = K = 0)
check_saturation  compare the end of a long curve with physics.flow_stress
"""
import numpy as np
import physics as ph

KB, EV, KB_EV = ph.KB, ph.EV, ph.KB_EV


# ---------------------------------------------------------------- density ODE
def _coeffs(gdot, T, P):
    """k2(T) and the climb coefficient K D(T)/gdot of the density equation."""
    T = np.asarray(T, float)
    k2 = P["k2_0"] * np.exp(P["q2_eV"] / (KB_EV * T))
    c3 = (P["K_cl"] * ph._diffusivity(T, P["D0_rec"], P["QD_rec_eV"]) / gdot
          if P["K_cl"] else np.zeros_like(T))
    return k2, c3


def dx_dgamma(x, gdot, T, P, S=1.0):
    """d sqrt(rho) / d gamma.  Same right-hand side that rho_steady sets to zero."""
    k2, c3 = _coeffs(gdot, T, P)
    x = np.maximum(x, 1e-3)
    return (P["k1"] * (x + P["Lambda_b"]) - k2 * x**2 - c3 * x**4) / (2.0 * x * S)


def rate_steps(edges, rates):
    """Piecewise-constant shear rate for a rate-jump test.

    rate_steps([g1, g2], [r0, r1, r2]) gives r0 for gamma < g1, r1 between g1 and
    g2, and r2 after g2.  Pass the result as `gdot` to stress_strain.
    """
    edges = np.atleast_1d(np.asarray(edges, float))
    rates = np.atleast_1d(np.asarray(rates, float))

    def f(g):
        return rates[min(int(np.searchsorted(edges, g, side="right")), len(rates) - 1)]
    f.first = float(rates[0])
    return f


def hardening_scale(P, T, eps_sat):
    """S that makes hardening saturate (to 98%) at axial strain eps_sat.

    From the Voce limit x = x_sat + (x_0-x_sat) exp(-k2 gamma / 2S): 98% is
    reached at gamma = 7.8 S / k2, and gamma = eps/M_s.
    """
    k2, _ = _coeffs(1.0, np.mean(np.atleast_1d(T)), P)
    return float(k2 * eps_sat / (7.8 * P["M_s"]))


def _relax_rate(gdot, T, P, S):
    """Fastest relaxation rate of the ODE, for choosing a stable step."""
    x = np.sqrt(ph.rho_steady(gdot, T, P))
    e = 1e-4 * x
    d = (dx_dgamma(x + e, gdot, T, P, S) - dx_dgamma(x - e, gdot, T, P, S)) / (2 * e)
    return float(np.max(np.abs(d)))


# ---------------------------------------------------------------- stresses
def state_from_rho(rho, gdot, T, P):
    """Densities, solute enrichment and tau_a for a density that is *given*.

    Identical to physics.state except that rho comes from the integration rather
    than from the steady-state root.
    """
    rho = np.asarray(rho, float)
    rho_f = P["beta"] * rho
    rho_m = P["f_m"] * rho_f
    lam_f = 1.0 / np.sqrt(rho_f)
    if P["C_sol"] > 0:
        t_a = rho_m * P["b"] * lam_f / gdot
        y = (t_a * P["nu_sol"] * np.exp(-P["Q_sol_eV"] / (KB_EV * T)))**P["alpha_sol"]
        C_C0 = 1.0 + P["C_sol"] * (1.0 - np.exp(-y))
    else:
        C_C0 = np.ones_like(rho)
    tau_a = P["alpha"] * ph.mu(T, P) * P["b"] * np.sqrt(C_C0) * (np.sqrt(rho_f) + P["Lambda_b"])
    return dict(rho=rho, rho_f=rho_f, rho_m=rho_m, lam_f=lam_f, C_C0=C_C0, tau_a=tau_a)


def tau_star_of(gdot, T, P, rho_m, iters=90):
    """Kink-pair stress that carries gdot with the given mobile density."""
    shape = np.broadcast(np.asarray(gdot, float), np.asarray(T, float), rho_m).shape
    lo = np.full(shape, np.log(1e-9))
    hi = np.full(shape, np.log(5e10))
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        up = rho_m * P["b"] * ph.glide_velocity(np.exp(mid), T, P) < gdot
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return np.exp(0.5 * (lo + hi))


def stress_from_rho(rho, gdot, T, P):
    """Every stress term at a point of the curve.  Same three terms as flow_stress."""
    S = state_from_rho(rho, gdot, T, P)
    S["tau_star"] = tau_star_of(gdot, T, P, S["rho_m"])
    S["tau_p"] = ph.particle_stress(gdot, T, P, S["rho_m"])
    S["tau"] = S["tau_a"] + S["tau_star"] + S["tau_p"]
    return S


# ---------------------------------------------------------------- the curve
def _march_isothermal(x0, gam, gd, T, P, S):
    """RK4 on x alone; the stress is not needed inside the loop."""
    x = np.empty((len(x0), len(gam)))
    x[:, 0] = x0
    for i in range(len(gam) - 1):
        h = gam[i + 1] - gam[i]
        xi, g0, gh, g1 = x[:, i], gd(gam[i]), gd(gam[i] + 0.5 * h), gd(gam[i + 1])
        k1 = dx_dgamma(xi, g0, T, P, S)
        k2 = dx_dgamma(xi + 0.5 * h * k1, gh, T, P, S)
        k3 = dx_dgamma(xi + 0.5 * h * k2, gh, T, P, S)
        k4 = dx_dgamma(xi + h * k3, g1, T, P, S)
        x[:, i + 1] = np.maximum(xi + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4), 1e-3)
    return x, np.repeat(np.asarray(T, float)[:, None], len(gam), axis=1)


def _march_adiabatic(x0, T0, gam, gd, P, S, eta, rho_mass, c_p):
    """Heun on (x, T); the plastic work tau dgamma heats the sample as it deforms."""
    c = eta / (rho_mass * c_p)
    x = np.empty((len(x0), len(gam)))
    T = np.empty_like(x)
    x[:, 0], T[:, 0] = x0, T0
    for i in range(len(gam) - 1):
        h = gam[i + 1] - gam[i]
        xi, Ti, g0, g1 = x[:, i], T[:, i], gd(gam[i]), gd(gam[i + 1])
        k1x = dx_dgamma(xi, g0, Ti, P, S)
        k1T = c * stress_from_rho(xi**2, g0, Ti, P)["tau"]
        xp = np.maximum(xi + h * k1x, 1e-3)
        Tp = Ti + h * k1T
        k2x = dx_dgamma(xp, g1, Tp, P, S)
        k2T = c * stress_from_rho(xp**2, g1, Tp, P)["tau"]
        x[:, i + 1] = np.maximum(xi + 0.5 * h * (k1x + k2x), 1e-3)
        T[:, i + 1] = Ti + 0.5 * h * (k1T + k2T)
    return x, T


def stress_strain(gdot, T, P, rho0, S=1.0, eps_max=0.5, n=301, adiabatic=None,
                  cfl=0.15, n_int_max=400000):
    """Stress-strain curve(s) at imposed SHEAR rate(s) gdot and temperature(s) T.

    gdot and T are scalars or 1-D arrays of the same length: one curve each.
    gdot may instead be a callable gdot(gamma) -- see rate_steps -- for a rate
    jump partway through the test.
    rho0  initial dislocation density (m^-2);  S  hardening scale (see module doc);
    eps_max  AXIAL strain the curve runs to;  n  points on the returned curve;
    adiabatic  None or dict(eta=, rho_mass=, c_p=).

    The equation is stiff, so it is integrated on an internal grid fine enough for
    the fastest relaxation rate (cfl) and then sampled onto the returned grid.

    Returns arrays of shape (n_curves, n): gamma, eps, tau, sigma, tau_a,
    tau_star, tau_p, rho, C_C0, T, gdot.
    """
    gd = gdot if callable(gdot) else (lambda g, v=np.atleast_1d(np.asarray(gdot, float)): v)
    g_ref = np.atleast_1d(np.asarray(gd(0.0), float))
    g_ref, T0 = np.broadcast_arrays(g_ref, np.atleast_1d(np.asarray(T, float)))
    g_ref, T0 = np.array(g_ref, float), np.array(T0, float)
    gam = np.linspace(0.0, eps_max / P["M_s"], n)
    x0 = np.full(T0.shape, np.sqrt(float(rho0)))

    lam = _relax_rate(g_ref, T0, P, S)                      # stable internal step
    n_int = int(min(max(n, gam[-1] * lam / cfl + 2), n_int_max))
    gi = np.linspace(0.0, gam[-1], n_int)

    if adiabatic is None:
        xi, Ti = _march_isothermal(x0, gi, gd, T0, P, S)
    else:
        xi, Ti = _march_adiabatic(x0, T0, gi, gd, P, S, **adiabatic)

    nC = len(T0)
    x = np.array([np.interp(gam, gi, xi[k]) for k in range(nC)])
    TT = np.array([np.interp(gam, gi, Ti[k]) for k in range(nC)])
    GG = np.broadcast_to(np.array([np.broadcast_to(gd(g), (nC,)) for g in gam]).T,
                         (nC, n)).copy()

    out = stress_from_rho(x**2, GG, TT, P)
    out.update(gamma=np.repeat(gam[None, :], nC, axis=0), T=TT, gdot=GG, n_int=n_int)
    out["eps"] = out["gamma"] * P["M_s"]
    out["sigma"] = out["tau"] / P["M_s"]
    return out


# ---------------------------------------------------------------- checks
def voce_analytic(gam, gdot, T, P, rho0, S=1.0):
    """Exact solution of the same equation when Lambda_b = 0, K_cl = 0, C_sol = 0.

    x = x_sat + (x_0 - x_sat) exp(-k2 gamma / 2S),  x_sat = k1 / k2,
    so tau_a is a Voce law.  Used only to verify the integrator.
    """
    k2, _ = _coeffs(gdot, T, P)
    x_sat = P["k1"] / k2
    x = x_sat + (np.sqrt(rho0) - x_sat) * np.exp(-k2 * np.asarray(gam, float) / (2.0 * S))
    return P["alpha"] * ph.mu(T, P) * P["b"] * np.sqrt(P["beta"]) * x


def check_saturation(gdot, T, P, rho0, S=1.0, eps_max=None, n=1201):
    """Largest relative difference between the end of a long curve and flow_stress.

    The consistency test between SS_physics and physics: the integrated curve must
    saturate on the steady-state solution the rest of the code uses.
    """
    if eps_max is None:
        eps_max = 40.0 * S / hardening_scale(P, T, 1.0)      # ~40 saturation strains
    r = stress_strain(gdot, T, P, rho0, S=S, eps_max=eps_max, n=n)
    ss = ph.flow_stress(np.atleast_1d(np.asarray(gdot, float)),
                        np.atleast_1d(np.asarray(T, float)), P)
    return float(np.max(np.abs(r["tau"][:, -1] - ss["tau"]) / ss["tau"]))
