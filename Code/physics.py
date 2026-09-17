"""
physics.py -- every equation of the unified kinetic flow rule, as plain functions.

Everything is in SHEAR (tau) space and SI units.  A material is a plain dict of
parameters (see parameters.py); terms that a material does not have are switched
off by setting their parameter to zero (C_sol = 0: no aging, tau_ppt = 0: no
particles, K_cl = 0: no climb recovery, Lambda_b = 0: no boundaries).

Flow rule for an imposed shear rate gdot and temperature T:

    gdot  = rho_m b v_g(tau_star)                       lattice (kink pairs)
    gdot  = rho_m b lambda_p [nu_th(tau_p)+nu_cl(tau_p)] particles (if present)
    tau   = tau_a(gdot,T) + tau_star + tau_p

with the state (rho, C/C0) fixed by (gdot, T) and not by tau.

Functions
---------
mu                shear modulus
glide_velocity    kink-pair velocity v_g(tau_star, T)
rho_steady        steady-state dislocation density (Kocks-Mecking + climb)
state             rho, rho_f, rho_m, C/C0, tau_a for given (gdot, T)
particle_stress   tau_p from the detachment || climb kinetics
flow_stress       tau and all its parts
srs               local m, V*, A, B, Lambda
global_m_V        m and V* fitted over several rates, as in experiment
diffusional_rate  Nabarro-Herring + Coble rate (maps only)
creep_rate        gdot at an imposed tau (stress control)
mechanism_map     rate contours and controlling-mechanism field
"""
import numpy as np

KB = 1.380649e-23          # J/K
EV = 1.602176634e-19       # J
KB_EV = KB / EV            # eV/K


# ---------------------------------------------------------------- basics
def mu(T, P):
    """Shear modulus (Pa).  Linear softening with a floor if mu_slope > 0."""
    T = np.asarray(T, float)
    m = P["mu0"] * (1.0 - P.get("mu_slope", 0.0) * (T - 300.0) / (P["Tm"] - 300.0))
    return np.maximum(m, P.get("mu_floor", 0.0) * P["mu0"]) if P.get("mu_slope", 0.0) else m * np.ones_like(T)


def _diffusivity(T, D0, Q_eV):
    return D0 * np.exp(-Q_eV / (KB_EV * np.asarray(T, float)))


# ---------------------------------------------------------------- travel
def glide_velocity(tau_star, T, P):
    """Kink-pair glide velocity (m/s), Po et al. (2016) form."""
    tau_star = np.asarray(tau_star, float)
    T = np.asarray(T, float)
    x = np.clip(tau_star / P["tau_P"], 0.0, 0.999999)
    dg = (1.0 - x**P["p"])**P["q"] - T / P["T0"]
    dg = np.maximum(dg, 0.0)
    B = P["B_kink"]
    if P.get("B_free"):                      # Zr: kink drag -> free-flight drag
        s = 1.0 / (1.0 + np.exp(-(0.05 - dg) / 0.05))
        B = P["B_kink"] * (1.0 - s) + P["B_free"] * s
    return tau_star * P["b"] / B * np.exp(-P["dH0_eV"] * dg / (2.0 * KB_EV * T))


# ---------------------------------------------------------------- state
def rho_steady(gdot, T, P):
    """Steady state of  drho/dgamma = k1(sqrt(rho)+Lambda_b) - k2 rho - K D rho^2/gdot."""
    gdot, T = np.broadcast_arrays(np.asarray(gdot, float), np.asarray(T, float))
    k2 = P["k2_0"] * np.exp(P["q2_eV"] / (KB_EV * T))
    c3 = P["K_cl"] * _diffusivity(T, P["D0_rec"], P["QD_rec_eV"]) / gdot if P["K_cl"] else np.zeros_like(T)
    f = lambda x: P["k1"] * (x + P["Lambda_b"]) - k2 * x**2 - c3 * x**4
    lo = np.full(T.shape, 1e-3)              # bisection in sqrt(rho)
    hi = np.full(T.shape, 1e10)
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        up = f(mid) > 0
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return np.sqrt(lo * hi)**2


def state(gdot, T, P):
    """Microstructural state for an imposed (gdot, T): densities, aging, tau_a."""
    gdot, T = np.broadcast_arrays(np.asarray(gdot, float), np.asarray(T, float))
    rho = rho_steady(gdot, T, P)
    rho_f = P["beta"] * rho
    rho_m = P["f_m"] * rho_f
    lam_f = 1.0 / np.sqrt(rho_f)
    if P["C_sol"] > 0:                        # solute atmospheres on arrested segments
        t_a = rho_m * P["b"] * lam_f / gdot
        x = (t_a * P["nu_sol"] * np.exp(-P["Q_sol_eV"] / (KB_EV * T)))**P["alpha_sol"]
        C_C0 = 1.0 + P["C_sol"] * (1.0 - np.exp(-x))
    else:
        C_C0 = np.ones_like(T)
    tau_a = P["alpha"] * mu(T, P) * P["b"] * np.sqrt(C_C0) * (np.sqrt(rho_f) + P["Lambda_b"])
    return dict(gdot=gdot, T=T, rho=rho, rho_f=rho_f, rho_m=rho_m, lam_f=lam_f,
                C_C0=C_C0, tau_a=tau_a)


# ---------------------------------------------------------------- particles
def particle_rate(tau_p, T, P, rho_m):
    """Shear rate carried by segments escaping particles at a stress tau_p."""
    tau_p = np.asarray(tau_p, float)
    T = np.asarray(T, float)
    kT = KB * T
    r = mu(T, P) / P["mu0"]
    tau_hat = P["tau_ppt"] * r                         # Orowan strength
    F = P["F_det_eV"] * EV * r                         # detachment energy
    y = np.clip(tau_p / tau_hat, 0.0, 1.0)
    nu_th = P["nu0_det"] * np.exp(-F * (1.0 - y**P["p_det"])**P["q_det"] / kT)
    u = np.minimum(tau_p * P["Omega"] / kT, 30.0)
    D_c = _diffusivity(T, P["D0_climb"], P["QD_climb_eV"])
    nu_cl = (2 * np.pi * D_c / (P["b"] * P["lnR_rc"])) * np.expm1(u) / \
            (P["h_ppt"] * np.maximum(1.0 - y, 1e-3))
    lam_p = 1.0 / np.sqrt(P["inv_lambda2_ppt"])
    return rho_m * P["b"] * lam_p * (nu_th + nu_cl)


def particle_stress(gdot, T, P, rho_m):
    """tau_p that makes the particle kinetics carry the imposed rate."""
    if P["tau_ppt"] <= 0:
        return np.zeros_like(np.asarray(T, float))
    lo = np.full(np.asarray(T, float).shape, np.log(1e-6))
    hi = np.full(np.asarray(T, float).shape, np.log(5e10))
    below = particle_rate(np.exp(lo), T, P, rho_m) >= gdot   # already fast enough
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        up = particle_rate(np.exp(mid), T, P, rho_m) < gdot
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return np.where(below, 0.0, np.exp(0.5 * (lo + hi)))


# ---------------------------------------------------------------- flow stress
def flow_stress(gdot, T, P):
    """Shear stress (Pa) and its parts for an imposed shear rate."""
    S = state(gdot, T, P)
    lo = np.full(S["T"].shape, np.log(1e-9))
    hi = np.full(S["T"].shape, np.log(5e10))
    for _ in range(90):                                  # bisection on tau_star
        mid = 0.5 * (lo + hi)
        up = S["rho_m"] * P["b"] * glide_velocity(np.exp(mid), S["T"], P) < S["gdot"]
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    tau_star = np.exp(0.5 * (lo + hi))
    tau_p = particle_stress(S["gdot"], S["T"], P, S["rho_m"])
    S["tau_star"], S["tau_p"] = tau_star, tau_p
    S["tau"] = S["tau_a"] + tau_star + tau_p
    return S


# ---------------------------------------------------------------- sensitivity
def srs(gdot, T, P, h=0.01):
    """Local rate sensitivity: m, V* (axial convention), A, B, Lambda = A B."""
    gdot = np.asarray(gdot, float)
    r = flow_stress(gdot, T, P)
    tp = flow_stress(gdot * np.exp(h), T, P)["tau"]
    tm = flow_stress(gdot * np.exp(-h), T, P)["tau"]
    dtau = (tp - tm) / (2 * h)                            # dtau/dln gdot
    ts = r["tau_star"]
    e = 1e-4                                              # A = dln(rho_m b v_g)/dtau*
    lnp = np.log(glide_velocity(ts * np.exp(e), r["T"], P))
    lnm = np.log(glide_velocity(ts * np.exp(-e), r["T"], P))
    A = (lnp - lnm) / (ts * 2 * np.sinh(e))
    B = dtau - 1.0 / A
    r.update(A=A, B=B, Lambda=A * B, dtau_dlngdot=dtau, m=dtau / r["tau"],
             Vstar=P["M_s"] * KB * r["T"] / dtau)
    r["Vstar_b3"] = r["Vstar"] / P["b"]**3
    return r


def global_m_V(rates, T, P):
    """m and V* fitted over several imposed AXIAL rates, as done in experiment.

    m  = dln(tau)/dln(gdot),  V* = kB T / (dsigma/dln eps_dot) with sigma = tau/M_s.
    """
    rates = np.asarray(rates, float)
    T = np.asarray(T, float)
    TT, EE = np.meshgrid(T, rates)
    tau = flow_stress(EE / P["M_s"], TT, P)["tau"]
    L = np.log(rates) - np.log(rates).mean()
    lt = np.log(tau)
    m = (L[:, None] * (lt - lt.mean(0))).sum(0) / (L**2).sum()
    slope = (L[:, None] * (tau - tau.mean(0))).sum(0) / (L**2).sum()     # dtau/dln gdot
    V = np.where(slope > 0, P["M_s"] * KB * T / slope, np.nan)
    return dict(T=T, tau=tau, m=m, Vstar=V, Vstar_b3=V / P["b"]**3)


# ---------------------------------------------------------------- creep
def diffusional_rate(tau, T, P):
    """Nabarro-Herring + Coble shear rate (used for the mechanism maps only)."""
    T = np.asarray(T, float)
    D_v = _diffusivity(T, P["D0_rec"], P["QD_rec_eV"])
    D_b = _diffusivity(T, P["Db0"], P["Qb_eV"])
    D_eff = D_v + np.pi * 2 * P["b"] / P["d_grain"] * D_b
    return 42.0 * P["Omega"] * np.asarray(tau, float) * D_eff / (KB * T * P["d_grain"]**2)


def rate_table(T, P, lg_min=-14.0, lg_max=6.0, n=241):
    """tau(gdot) on a log rate grid for each T (used to invert to stress control)."""
    T = np.atleast_1d(np.asarray(T, float))
    lg = np.linspace(lg_min, lg_max, n)
    GG, TT = np.meshgrid(10.0**lg, T)
    return lg, flow_stress(GG, TT, P)["tau"]


def creep_rate(tau, T, P, diffusion=False, table=None):
    """Dislocation shear rate at an imposed shear stress (plus diffusional flow).

    Returns dict(total, disl, diff).  Where tau(gdot) is non-monotonic (DSA) the
    lowest-rate branch that reaches the stress is used.
    """
    T = np.atleast_1d(np.asarray(T, float))
    lg, tab = table if table is not None else rate_table(T, P)
    tau = np.atleast_2d(np.asarray(tau, float))
    if tau.shape[0] != len(T):
        tau = np.broadcast_to(tau, (len(T), tau.shape[-1]))
    disl = np.empty(tau.shape)
    for i in range(len(T)):
        env = np.maximum.accumulate(np.nan_to_num(tab[i], nan=np.inf))
        env = env + 1e-9 * np.arange(len(env))            # strictly increasing
        with np.errstate(divide="ignore"):
            disl[i] = 10.0**np.interp(np.log(tau[i]), np.log(np.maximum(env, 1e-30)), lg,
                                      left=-np.inf, right=np.inf)
    diff = diffusional_rate(tau, T[:, None], P) if diffusion else np.zeros_like(disl)
    return dict(total=disl + diff, disl=disl, diff=diff)


# ---------------------------------------------------------------- maps
def mechanism_map(T, tau_over_mu, P, diffusion=True):
    """Iso-rate contours and controlling mechanism on a (T, tau/mu) grid.

    field: 0 kink-pair glide (tau*/tau > 0.5)     1 obstacle / athermal plateau
           3 climb or recovery controlled          4 diffusional flow
          -1 outside the tabulated rate range.  The m < 0 (DSA) band is returned
    separately as (tau_lo, tau_hi) per temperature: it is unstable at fixed stress.
    """
    T = np.atleast_1d(np.asarray(T, float))
    tau = tau_over_mu[None, :] * mu(T, P)[:, None]
    lg, tab = rate_table(T, P)
    parts = creep_rate(tau, T, P, diffusion=diffusion, table=(lg, tab))
    GG, TT = np.meshgrid(10.0**lg, T)
    d = srs(GG, TT, P)                                    # diagnostics on the rate grid
    field = np.zeros(tau.shape, int)
    lg_d = np.log10(np.clip(parts["disl"], 1e-300, None))
    for i in range(len(T)):
        j = np.clip(np.searchsorted(lg, lg_d[i]), 0, len(lg) - 1)
        frac = d["tau_star"][i, j] / d["tau"][i, j]
        lam = d["Lambda"][i, j]
        f = np.where(frac > 0.5, 0, 1)
        field[i] = np.where((frac <= 0.5) & (lam > 1), 3, f)
    field = np.where(parts["diff"] > parts["disl"], 4, field)
    field = np.where(np.isinf(parts["disl"]) | (lg_d > lg[-1] - 1.0), -1, field)
    dsa = d["Lambda"] < -1
    lo = np.where(dsa.any(1), np.nanmin(np.where(dsa, d["tau"], np.inf), 1), np.nan)
    hi = np.where(dsa.any(1), np.nanmax(np.where(dsa, d["tau"], -np.inf), 1), np.nan)
    return dict(T=T, tau_over_mu=tau_over_mu, tau=tau, rate=parts["total"], field=field,
                dsa_band=(lo, hi))
