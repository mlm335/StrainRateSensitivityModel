"""
unified_plasticity.py
=====================

One kinetic framework for rate- and temperature-dependent plasticity and
creep of BCC metals and alloys (and HCP prismatic slip), generalising
bcc_unified_model.py.

Flow rule (time per unit glide distance, Rules S and P):

    1/gdot = 1/gdot_glide(tau*) + 1/gdot_wait(tau*)

    gdot_glide = rho_m b v_g(tau*, T)                       travel
    gdot_wait  = rho_m b / sum_j [ n_j / nu_j(tau*) ]       waiting
    nu_j       = sum_k nu_jk(tau*)                          ways past obstacle j
    n_j        = encounters of obstacle j per unit glide distance
                 point obstacles : n_j = lambda_j^-2 * lambda_pt

    tau = tau_a(state) + tau*

tau_a holds *athermal* resistances (Orowan-type thresholds, long-range
stresses): Taylor forest stress, boundary/lath spacing stress.  These
superpose linearly with the lattice friction (Kocks: dense weak + sparse
strong obstacles).  Thermally penetrable obstacles and dislocation sources
enter as waiting steps.

The state (rho, C/C0, d, ...) may depend on (gdot, T, gamma) but NOT on tau*,
so gdot(tau*) is monotonic and the stress for an imposed rate is unique.
Rate dependence of the state produces the "rate-sensitive resistance" B:

    1/m = tau A / (1 + A B),   A = dlnG/dtau*|state,
    B   = -(dlnG/dln gdot|tau*)/A,   Lambda = A B,   V* = kT A/(1+Lambda)

Rate-dependent state in this module:
  * solute aging (DSA)            -> B < 0   (Lambda < -1  <=>  m < 0)
  * precipitate bypass by climb   -> B > 0   (threshold-type creep)
  * climb-controlled recovery     -> B > 0   (creep-like softening, m -> 1/3)

Only numpy is required.  SI units.
"""
from __future__ import annotations

import numpy as np

KB = 1.380649e-23
EV = 1.602176634e-19
NA = 6.02214076e23
KB_EV = KB / EV


def eV(x):
    return np.asarray(x, dtype=float) * EV


def kJmol(x):
    return np.asarray(x, dtype=float) * 1e3 / NA


# =============================================================================
# helpers
# =============================================================================
def _scale(S, kind, target):
    return S["scale"].get((kind, target), 1.0)


def _mul_scale(S, kind, target, factor):
    S["scale"][(kind, target)] = _scale(S, kind, target) * factor


def _logsumexp(arrs):
    arrs = np.broadcast_arrays(*arrs)
    stack = np.stack(arrs, axis=0)
    mx = np.max(stack, axis=0)
    mx_safe = np.where(np.isfinite(mx), mx, 0.0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out = mx_safe + np.log(np.sum(np.exp(stack - mx_safe), axis=0))
    return np.where(np.isfinite(mx), out, mx)


def _full(S, v):
    return np.broadcast_to(np.asarray(v, dtype=float), S["T"].shape).astype(float)


# =============================================================================
# Material
# =============================================================================
class Material:
    """mu_model: 'constant' | 'linear' | callable(T)."""

    def __init__(self, name, b, mu0, Tm, Omega=None, density=None, cp=None,
                 taylor_quinney=0.95, mu_model="constant", mu_slope=0.55,
                 mu_floor=0.15):
        self.name, self.b, self.mu0, self.Tm = name, b, mu0, Tm
        self.Omega = Omega if Omega is not None else b**3
        self.density, self.cp, self.taylor_quinney = density, cp, taylor_quinney
        self.mu_model, self.mu_slope, self.mu_floor = mu_model, mu_slope, mu_floor

    def mu(self, T):
        T = np.asarray(T, dtype=float)
        if callable(self.mu_model):
            return self.mu_model(T)
        if self.mu_model == "constant":
            return self.mu0 * np.ones_like(T)
        return self.mu0 * np.maximum(self.mu_floor, 1.0 - self.mu_slope * T / self.Tm)


class Diffusivity:
    """D_eff = D_L0 exp(-Q_L/kT) + f_p * D_p0 exp(-Q_p/kT)
    f_p = pi a_c^2 rho_f (Hart) if f_pipe is None, else f_pipe."""

    def __init__(self, D_L0, Q_L, D_p0=0.0, Q_p=None, a_c=None, f_pipe=None):
        self.D_L0, self.Q_L, self.D_p0 = D_L0, Q_L, D_p0
        self.Q_p = Q_p if Q_p is not None else Q_L
        self.a_c, self.f_pipe = a_c, f_pipe

    def __call__(self, S, mat, rho=None):
        kT = KB * S["T"]
        D = self.D_L0 * np.exp(-self.Q_L / kT)
        if self.D_p0:
            if self.f_pipe is None:
                a_c = self.a_c if self.a_c is not None else 2.0 * mat.b
                r = rho if rho is not None else S.get("rho_f", 0.0)
                fp = np.pi * a_c**2 * r
            else:
                fp = self.f_pipe
            D = D + fp * self.D_p0 * np.exp(-self.Q_p / kT)
        return D


# =============================================================================
# Microstructure (state) providers -- must not depend on tau
# =============================================================================
class Microstructure:
    def update(self, S, mat):
        raise NotImplementedError


class FixedMicrostructure(Microstructure):
    def __init__(self, **values):
        self.values = values

    def update(self, S, mat):
        for k, v in self.values.items():
            S[k] = _full(S, v)


class KocksMecking(Microstructure):
    """Forest density with boundary storage and climb-controlled recovery.

        d rho/d gamma = k1 (sqrt(rho) + Lambda) - k2(T, gdot) rho
                        - (K_cl D_eff(T) / gdot) rho^2

      k2      = k20 exp(q2[eV]/kT) [ (gdot_ref/gdot)^(kT/A_r) ]   (cross-slip type)
      Lambda  = sum_i c_i / d_i  (boundaries, laths, particles; 1/m)
      last term: diffusion (climb) controlled recovery, dRho/dt = -K_cl D rho^2.
                 It makes rho, hence tau_a, rate dependent (B > 0) and gives the
                 steady-state creep law tau ~ (gdot/D)^(1/3) (n = 3) when dominant.

    saturate=True: steady state (root of RHS = 0; unique, RHS concave in sqrt(rho)).
    saturate=False: integrate from rho0 to gamma (RK4).
    rho_f = beta * rho,  rho_m = f_m * rho_f.
    """

    def __init__(self, k1, k20, q2=0.0, A_r=None, gdot_ref=1e7, Lambda=0.0,
                 K_cl=0.0, diffusivity=None, rho0=1e12, saturate=True,
                 beta=1.0, f_m=1.0, n_steps=200):
        self.k1, self.k20, self.q2, self.A_r, self.gdot_ref = k1, k20, q2, A_r, gdot_ref
        self.Lambda, self.K_cl, self.diff = Lambda, K_cl, diffusivity
        self.rho0, self.saturate, self.beta, self.f_m = rho0, saturate, beta, f_m
        self.n_steps = n_steps

    def k2(self, T, gdot):
        kT_eV = KB_EV * T
        k2 = self.k20 * np.exp(self.q2 / kT_eV)          # q2 in eV
        if self.A_r is not None:
            k2 = k2 * np.power(self.gdot_ref / gdot, KB * T / self.A_r)
        return k2

    def _rhs_x(self, x, k2, c3, lam):
        # f(x) with x = sqrt(rho): d rho/d gamma
        return self.k1 * (x + lam) - k2 * x**2 - c3 * x**4

    def update(self, S, mat):
        T, g = S["T"], S["gdot"]
        k2 = self.k2(T, g)
        lam = _full(S, self.Lambda)
        if self.K_cl and self.diff is not None:
            # pipe weighting uses a first guess of rho (saturated, no climb)
            rho_guess = ((self.k1 + np.sqrt(self.k1**2 + 4 * self.k1 * k2 * lam)) / (2 * k2))**2
            c3 = self.K_cl * self.diff(S, mat, rho=rho_guess) / g
        else:
            c3 = np.zeros_like(T)
        if self.saturate:
            lo = np.full(T.shape, 1e-3)
            hi = np.full(T.shape, 1e10)
            for _ in range(80):
                mid = np.sqrt(lo * hi)
                pos = self._rhs_x(mid, k2, c3, lam) > 0
                lo = np.where(pos, mid, lo)
                hi = np.where(pos, hi, mid)
            x = np.sqrt(lo * hi)
            rho = x**2
        else:
            rho = _full(S, self.rho0)
            gam = S["gamma"]
            h = gam / self.n_steps
            f = lambda r: self._rhs_x(np.sqrt(r), k2, c3, lam)
            for _ in range(self.n_steps):
                a1 = f(rho); a2 = f(np.maximum(rho + 0.5 * h * a1, 1.0))
                a3 = f(np.maximum(rho + 0.5 * h * a2, 1.0)); a4 = f(np.maximum(rho + h * a3, 1.0))
                rho = np.maximum(rho + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4), 1.0)
        S["rho"] = rho
        S["rho_f"] = self.beta * rho
        S["rho_m"] = self.f_m * S["rho_f"]
        S["Lambda_b"] = lam


class IrradiationLoops(Microstructure):
    def __init__(self, N0=1e22, d0=5e-9, gamma_c=None):
        self.N0, self.d0, self.gamma_c = N0, d0, gamma_c

    def update(self, S, mat):
        N = _full(S, self.N0)
        if self.gamma_c is not None:
            N = N * np.exp(-S["gamma"] / self.gamma_c)
        S["N_loop"] = N
        S["d_loop"] = _full(S, self.d0)


# =============================================================================
# Athermal (threshold) stresses
# =============================================================================
class Athermal:
    name = "athermal"

    def tau(self, S, mat):
        raise NotImplementedError


class ConstantAthermal(Athermal):
    def __init__(self, tau0, scale_with_mu=False, name="const"):
        self.tau0, self.scale_with_mu, self.name = tau0, scale_with_mu, name

    def tau(self, S, mat):
        t = _full(S, self.tau0)
        return t * mat.mu(S["T"]) / mat.mu0 if self.scale_with_mu else t


class TaylorAthermal(Athermal):
    """tau_a = alpha mu(T) b [ sqrt(rho_f) + Lambda_b ]
    multiplied by S scale ('tau_hat', name) -> solute aging can strengthen it
    (junction/forest pinning by atmospheres)."""

    def __init__(self, alpha, include_boundaries=True, name="forest"):
        self.alpha, self.inc_b, self.name = alpha, include_boundaries, name

    def tau(self, S, mat):
        inv_len = np.sqrt(S["rho_f"])
        if self.inc_b and "Lambda_b" in S:
            inv_len = inv_len + S["Lambda_b"]
        return self.alpha * mat.mu(S["T"]) * mat.b * inv_len * _scale(S, "tau_hat", self.name)


class PowerLawAthermal(Athermal):
    def __init__(self, tau0, n1, name="powerlaw"):
        self.tau0, self.n1, self.name = tau0, n1, name

    def tau(self, S, mat):
        return self.tau0 * np.power(S["gamma"], self.n1)


# =============================================================================
# Travel (glide) laws -> ln(rho_m b v)
# =============================================================================
class GlideLaw:
    name = "glide"

    def ln_rate(self, tau_s, S, mat):
        raise NotImplementedError


class KinkPairGlide(GlideLaw):
    """Po et al. (2016)-type screw mobility (long-segment kink-pair regime)

        v = tau* b / B * exp(-dH0 max[(1-(tau*/tau_P)^p)^q - T/T0, 0] / 2kT)

    B = B_kink, or (drag_switch=True) the smooth kink->free-drag interpolation
    of the original law:  B = B_kink (1-s) + B_free s,  s = s(dg).
    Setting dH0 = 0 gives pure drag (FCC).
    """

    def __init__(self, dH0, tau_P, p, q, B_kink, T0=None, B_free=None,
                 drag_switch=False, B_T=0.0, name="glide"):
        self.dH0, self.tau_P, self.p, self.q = dH0, tau_P, p, q
        self.B_kink, self.T0, self.B_free = B_kink, T0, B_free
        self.drag_switch, self.B_T, self.name = drag_switch, B_T, name

    def ln_rate(self, tau_s, S, mat):
        T = S["T"]
        dH = self.dH0 * _scale(S, "F", self.name)
        tP = self.tau_P * _scale(S, "tau_hat", self.name)
        x = np.clip(tau_s / tP, 0.0, 0.999999)
        dg = (1.0 - x**self.p) ** self.q
        if self.T0 is not None:
            dg = dg - T / self.T0
        dg = np.maximum(dg, 0.0)
        B = self.B_kink + self.B_T * T
        if self.drag_switch and self.B_free is not None:
            sg = 0.5 * 2.0 / (1.0 + np.exp(2.0 * (-0.5 * (0.05 - dg) / 0.05)))
            B = B * (1 - sg) + self.B_free * sg
        with np.errstate(divide="ignore"):
            ln_v = np.log(tau_s * mat.b / B) - dH * dg / (2.0 * KB * T)
        return np.log(S["rho_m"] * mat.b) + ln_v


class KocksGlide(GlideLaw):
    def __init__(self, gdot0, F0, tau_hat0, p=2 / 3, q=2.0, name="glide"):
        self.gdot0, self.F0, self.tau_hat0, self.p, self.q, self.name = gdot0, F0, tau_hat0, p, q, name

    def ln_rate(self, tau_s, S, mat):
        F = self.F0 * _scale(S, "F", self.name)
        th = self.tau_hat0 * _scale(S, "tau_hat", self.name)
        x = np.clip(tau_s / th, 0.0, 1.0)
        return np.log(self.gdot0) - F * (1.0 - x**self.p) ** self.q / (KB * S["T"])


# =============================================================================
# Obstacles / sources and the ways past them
# =============================================================================
class Route:
    name = "route"

    def ln_nu(self, tau_s, S, mat, obs):
        raise NotImplementedError


class ThermalRoute(Route):
    """nu0 exp{-(F/kT)[1-(tau*/tau_hat)^p]^q}; barrier vanishes at tau_hat."""

    def __init__(self, nu0=1e11, p=2 / 3, q=1.5, name="thermal"):
        self.nu0, self.p, self.q, self.name = nu0, p, q, name

    def ln_nu(self, tau_s, S, mat, obs):
        F, th = obs.F(S, mat), obs.tau_hat(S, mat)
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.nan_to_num(np.clip(tau_s / th, 0.0, 1.0), nan=1.0)
        return np.log(self.nu0) - F * (1.0 - x**self.p) ** self.q / (KB * S["T"])


class ClimbRoute(Route):
    """Glide-assisted climb: nu = v_c/h,
    v_c = 2 pi D_eff /(b ln(R/rc)) [exp(beta tau* Omega/kT) - 1]."""

    def __init__(self, diffusivity, beta=1.0, ln_R_rc=5.0, name="climb",
                 stress_assisted_height=False, h_min=1e-3):
        self.diff, self.beta, self.ln_R_rc, self.name = diffusivity, beta, ln_R_rc, name
        self.sah, self.h_min = stress_assisted_height, h_min

    def ln_nu(self, tau_s, S, mat, obs):
        kT = KB * S["T"]
        u = self.beta * tau_s * mat.Omega / kT
        with np.errstate(divide="ignore", over="ignore"):
            ln_drive = np.where(u < 30, np.log(np.expm1(np.minimum(u, 30))), u)
            h = obs.height(S, mat)
            if self.sah:
                # local climb: the height still to be climbed shrinks as the
                # line is pushed against the obstacle (Brown-Ham type)
                x = np.clip(np.nan_to_num(tau_s / obs.tau_hat(S, mat), nan=1.0), 0.0, 1.0)
                h = h * np.maximum(1.0 - x, self.h_min)
            return (np.log(2 * np.pi * self.diff(S, mat) / (mat.b * self.ln_R_rc))
                    + ln_drive - np.log(h))


class Obstacle:
    kind = "point"

    def __init__(self, name, routes, mixed_spacing=True):
        self.name, self.routes, self.mixed_spacing = name, list(routes), mixed_spacing

    # geometry
    def inv_lambda2(self, S, mat):
        raise NotImplementedError

    def encounters(self, S, mat):
        """number met per unit glide distance"""
        return self.inv_lambda2(S, mat) * S["lambda_pt"]

    def height(self, S, mat):
        return 10 * mat.b * np.ones_like(S["T"])

    # strength
    def _tau_hat(self, S, mat):
        raise NotImplementedError

    def _F(self, S, mat):
        raise NotImplementedError

    def spacing_inv2(self, S, mat):
        if self.mixed_spacing and "inv_lambda2" in S:
            return S["inv_lambda2"]
        return self.inv_lambda2(S, mat)

    def tau_hat(self, S, mat):
        return self._tau_hat(S, mat) * _scale(S, "tau_hat", self.name)

    def F(self, S, mat):
        return self._F(S, mat) * _scale(S, "F", self.name)


class ForestObstacle(Obstacle):
    """Thermally penetrable forest junctions (e.g. FCC cutting)."""

    def __init__(self, routes, alpha=0.3, g_F=0.5, n_h=10.0, name="forest_cut",
                 mixed_spacing=True):
        super().__init__(name, routes, mixed_spacing)
        self.alpha, self.g_F, self.n_h = alpha, g_F, n_h

    def inv_lambda2(self, S, mat):
        return S["rho_f"]

    def _tau_hat(self, S, mat):
        return self.alpha * mat.mu(S["T"]) * mat.b * np.sqrt(self.spacing_inv2(S, mat))

    def _F(self, S, mat):
        return self.g_F * mat.mu(S["T"]) * mat.b**3

    def height(self, S, mat):
        return self.n_h * mat.b * np.ones_like(S["T"])


class LoopObstacle(Obstacle):
    def __init__(self, routes, alpha=0.3, g_F=0.3, name="loops", mixed_spacing=True):
        super().__init__(name, routes, mixed_spacing)
        self.alpha, self.g_F = alpha, g_F

    def inv_lambda2(self, S, mat):
        return S["N_loop"] * S["d_loop"]

    def _tau_hat(self, S, mat):
        return self.alpha * mat.mu(S["T"]) * mat.b * np.sqrt(self.spacing_inv2(S, mat))

    def _F(self, S, mat):
        return self.g_F * mat.mu(S["T"]) * mat.b**3

    def height(self, S, mat):
        return S["d_loop"]


class PrecipitateObstacle(Obstacle):
    """Non-shearable particles / solute clusters with a given strength.
    tau_hat = tau_hat0 * mu(T)/mu0 (Orowan-type), F = g_F mu b^3,
    1/lambda^2 = inv_lambda2 (given), climb height h."""

    def __init__(self, routes, tau_hat0, inv_lambda2, h, g_F=5.0, name="precipitates"):
        super().__init__(name, routes, mixed_spacing=False)
        self.tau_hat0, self.il2, self.h, self.g_F = tau_hat0, inv_lambda2, h, g_F

    def inv_lambda2(self, S, mat):
        return _full(S, self.il2)

    def _tau_hat(self, S, mat):
        return self.tau_hat0 * mat.mu(S["T"]) / mat.mu0

    def _F(self, S, mat):
        return self.g_F * mat.mu(S["T"]) * mat.b**3

    def height(self, S, mat):
        return _full(S, self.h)


class ResistanceGroup:
    """Obstacles that act SIMULTANEOUSLY with lattice friction along the line
    (sparse, strong pinning points: particles, clusters).  Their stress adds
    to the stress of the main (sequential) process -- Kocks' linear
    superposition of dense-weak and sparse-strong obstacles:

        tau = tau_a + tau*_main(gdot) + sum_g tau_g(gdot)

    tau_g solves  rho_m b / sum_j [n_j / nu_j(tau_g)] = gdot   (Rules S, P
    inside the group).  An unbypassable group gives a rate-independent
    (athermal) tau_g; a climb route makes it rate sensitive at high T
    (threshold-type creep)."""

    def __init__(self, name, obstacles):
        self.name, self.obstacles = name, list(obstacles)

    def ln_rate(self, tau_g, S, mat):
        saved = {k: S.get(k) for k in ("inv_lambda2", "lambda_pt")}
        inv = sum(o.inv_lambda2(S, mat) for o in self.obstacles)
        S["inv_lambda2"], S["lambda_pt"] = inv, 1.0 / np.sqrt(inv)
        terms = []
        for o in self.obstacles:
            ln_nu = _logsumexp([r.ln_nu(tau_g, S, mat, o) for r in o.routes])
            with np.errstate(divide="ignore"):
                terms.append(np.log(o.encounters(S, mat)) - ln_nu)
        for k, v in saved.items():
            if v is None:
                S.pop(k, None)
            else:
                S[k] = v
        return np.log(S["rho_m"] * mat.b) - _logsumexp(terms)

    def solve(self, S, mat, n_iter=80):
        ln_t = np.log(S["gdot"])
        lo = np.full(S["T"].shape, np.log(1e-6))
        hi = np.full(S["T"].shape, np.log(5e10))
        f_lo = self.ln_rate(np.exp(lo), S, mat) - ln_t
        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            up = self.ln_rate(np.exp(mid), S, mat) - ln_t < 0
            lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
        t = np.exp(0.5 * (lo + hi))
        return np.where(f_lo >= 0, 0.0, t)


# =============================================================================
# Solute aging
# =============================================================================
class SoluteAging:
    """C/C0 = 1 + sum_i C_i [1 - exp(-x_i)],  scale sqrt(C/C0) on targets.

    mode='waiting': x_i = (t_a nu_i exp(-Q_i/kT))^alpha_i,
                    t_a = rho_m b lambda / gdot  (lambda: point-obstacle spacing,
                    or 1/sqrt(rho_f) if there are no point obstacles)
    mode='cheng'  : Cheng et al. (2001) Eq. 33-34 (species key 'Omega').
    Targets are names of glide laws, obstacles or TaylorAthermal terms.
    exponent: scale = (C/C0)^exponent (0.5 = Cheng).
    """

    def __init__(self, species, targets=("forest",), mode="waiting", m=0.8,
                 m0=0.25, gdot00=1.0, exponent=0.5, kinds=("F", "tau_hat")):
        self.species, self.targets, self.mode = species, tuple(targets), mode
        self.m, self.m0, self.gdot00, self.exponent, self.kinds = m, m0, gdot00, exponent, kinds

    def C_over_C0(self, S, mat):
        T, g = S["T"], S["gdot"]
        kT = KB * T
        cc = np.ones_like(T)
        for sp in self.species:
            if self.mode == "cheng":
                Om = sp["Omega"] * np.power(g / self.gdot00, self.m * np.exp(-g / self.gdot00) - self.m0)
                x = (Om / g * np.exp(-sp["Q"] / kT)) ** sp["alpha"]
            else:
                lam = S["lambda"]
                t_a = S["rho_m"] * mat.b * lam / g
                x = (t_a * sp["nu"] * np.exp(-sp["Q"] / kT)) ** sp["alpha"]
            cc = cc + sp["C"] * (1.0 - np.exp(-x))
        return cc

    def update(self, S, mat):
        cc = self.C_over_C0(S, mat)
        S["C_C0"] = cc
        s = cc ** self.exponent
        for t in self.targets:
            for k in self.kinds:
                _mul_scale(S, k, t, s)


# =============================================================================
# Model
# =============================================================================
class DislocationModel:
    def __init__(self, material, glide, obstacles=(), athermal=(),
                 microstructure=(), aging=None, groups=(), name=""):
        self.groups = list(groups)
        self.mat, self.glide = material, glide
        self.obstacles, self.athermal = list(obstacles), list(athermal)
        self.microstructure, self.aging, self.name = list(microstructure), aging, name

    def build_state(self, gdot, T, gamma=0.0):
        gdot, T, gamma = np.broadcast_arrays(np.asarray(gdot, float), np.asarray(T, float),
                                             np.asarray(gamma, float))
        S = {"gdot": gdot, "T": T, "gamma": gamma, "scale": {}}
        for ms in self.microstructure:
            ms.update(S, self.mat)
        pts = [o for o in self.obstacles if o.kind == "point"]
        if pts:
            inv = sum(o.inv_lambda2(S, self.mat) for o in pts)
            S["inv_lambda2"] = inv
            S["lambda_pt"] = 1.0 / np.sqrt(inv)
            S["lambda"] = S["lambda_pt"]
        elif "rho_f" in S:
            S["lambda"] = 1.0 / np.sqrt(S["rho_f"])
        if self.aging is not None:
            self.aging.update(S, self.mat)
        tau_a = np.zeros_like(T)
        for a in self.athermal:
            tau_a = tau_a + a.tau(S, self.mat)
        S["tau_a"] = tau_a
        return S

    def ln_rates(self, tau_s, S):
        out = {"glide": self.glide.ln_rate(tau_s, S, self.mat)}
        terms = []
        for o in self.obstacles:
            routes = {r.name: r.ln_nu(tau_s, S, self.mat, o) for r in o.routes}
            ln_nu = _logsumexp(list(routes.values()))
            for rn, v in routes.items():
                out[f"frac_{o.name}_{rn}"] = v - ln_nu
            with np.errstate(divide="ignore"):
                t = np.log(o.encounters(S, self.mat)) - ln_nu
            out[f"ln_t_{o.name}"] = t
            terms.append(t)
        if terms:
            ln_tw = _logsumexp(terms)                       # time per unit length
            out["wait"] = np.log(S["rho_m"] * self.mat.b) - ln_tw
        else:
            out["wait"] = np.full_like(out["glide"], np.inf)
        out["total"] = -_logsumexp([-out["glide"], -out["wait"]])
        return out

    def solve(self, gdot, T, gamma=0.0, tau_lo=1e-9, tau_hi=5e10, n_iter=90,
              details=True, S=None):
        S = S if S is not None else self.build_state(gdot, T, gamma)
        ln_t = np.log(S["gdot"])
        lo = np.full(S["T"].shape, np.log(tau_lo))
        hi = np.full(S["T"].shape, np.log(tau_hi))
        f_lo = self.ln_rates(np.exp(lo), S)["total"] - ln_t
        f_hi = self.ln_rates(np.exp(hi), S)["total"] - ln_t
        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            up = self.ln_rates(np.exp(mid), S)["total"] - ln_t < 0
            lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
        ts = np.exp(0.5 * (lo + hi))
        ts = np.where(f_lo >= 0, 0.0, ts)
        ts = np.where(f_hi < 0, np.nan, ts)
        res = {k: np.array(S[k], dtype=float) for k in
               ("gdot", "T", "gamma", "tau_a", "rho_f", "rho_m", "lambda", "C_C0")
               if k in S}
        res["tau_star"] = ts
        tau_groups = np.zeros_like(ts)
        for g in self.groups:
            tg = g.solve(S, self.mat)
            res["tau_" + g.name] = tg
            tau_groups = tau_groups + tg
        res["tau_groups"] = tau_groups
        res["tau"] = S["tau_a"] + ts + tau_groups
        if details:
            r = self.ln_rates(np.where(ts > 0, ts, tau_lo), S)
            res["w_glide"] = np.exp(r["total"] - r["glide"])
            res["w_wait"] = 1 - res["w_glide"]
            for k, v in r.items():
                if k.startswith("frac_"):
                    res[k] = np.exp(v)
        return res

    # ---------------- strain-rate sensitivity ------------------------------
    def ln_rate_intra(self, tau_s, S):
        """Intrinsic mobility process: glide + point obstacles (series)."""
        g = self.glide.ln_rate(tau_s, S, self.mat)
        terms = []
        for o in self.obstacles:
            if o.kind != "point":
                continue
            ln_nu = _logsumexp([r.ln_nu(tau_s, S, self.mat, o) for r in o.routes])
            terms.append(np.log(o.encounters(S, self.mat)) - ln_nu)
        if not terms:
            return g
        w = np.log(S["rho_m"] * self.mat.b) - _logsumexp(terms)
        return -_logsumexp([-g, -w])

    def srs(self, gdot, T, gamma=0.0, h=0.01):
        """Local strain-rate sensitivity and its A-B decomposition.

        A      = d ln(gdot_intra)/d tau*   (intrinsic mobility sensitivity:
                 glide + intragranular obstacles in series, evaluated at the
                 stress tau*_intra this process needs on its own at gdot)
        B      = d tau/d ln gdot - 1/A     (everything else: sources/boundaries,
                 superposed groups, rate-dependent state such as aging or
                 climb recovery)
        Lambda = A B ;  1/m = tau A/(1+Lambda) ;  V* = kT/(m tau) = kT A/(1+Lambda)
        """
        gdot = np.asarray(gdot, float)
        S = self.build_state(gdot, T, gamma)
        r = self.solve(gdot, T, gamma, S=S)
        tp = self.solve(gdot * np.exp(h), T, gamma, details=False)["tau"]
        tm = self.solve(gdot * np.exp(-h), T, gamma, details=False)["tau"]
        dtau = (tp - tm) / (2 * h)
        # intrinsic process solved on its own at the imposed rate (its share of tau)
        ln_t = np.log(S["gdot"])
        lo = np.full(S["T"].shape, np.log(1e-9)); hi = np.full(S["T"].shape, np.log(5e10))
        for _ in range(90):
            mid = 0.5 * (lo + hi)
            up = self.ln_rate_intra(np.exp(mid), S) - ln_t < 0
            lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
        ts = np.exp(0.5 * (lo + hi))
        e = 1e-4
        A = (self.ln_rate_intra(ts * np.exp(e), S) - self.ln_rate_intra(ts * np.exp(-e), S)) / (ts * 2 * np.sinh(e))
        r["tau_star_intra"] = ts
        B = dtau - 1.0 / A
        Lam = A * B
        m = dtau / r["tau"]
        with np.errstate(divide="ignore", invalid="ignore"):
            V = KB * r["T"] / dtau
        r.update(A=A, B=B, Lambda=Lam, m=m, dtau_dlngdot=dtau, Vstar=V,
                 Vstar_b3=V / self.mat.b**3)
        return r

    def global_m_V(self, axial_rates, T, M, gamma=0.0):
        """Fit over several rates, as done experimentally:
        m = dln(tau)/dln(eps_dot), V* = kT / (d sigma / dln eps_dot), sigma = tau/M."""
        axial_rates = np.asarray(axial_rates, float)
        T = np.asarray(T, float)
        TT, EE = np.meshgrid(T, axial_rates)
        tau = self.solve(EE / M if np.ndim(M) == 0 else EE * M, TT, gamma, details=False)["tau"]
        L = np.log(axial_rates)
        Lc = L - L.mean()
        lt = np.log(tau)
        m = (Lc[:, None] * (lt - lt.mean(0))).sum(0) / (Lc**2).sum()
        sig = tau / M
        slope = (Lc[:, None] * (sig - sig.mean(0))).sum(0) / (Lc**2).sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            V = np.where(slope > 0, KB * T / slope, np.nan)
        return dict(T=T, tau=tau, m=m, Vstar=V, Vstar_b3=V / self.mat.b**3)

    def stress_strain(self, gdot, T0, gamma_max, n=200, adiabatic="auto",
                      gdot_adiabatic=1.0, gamma_min=1e-4):
        T0 = np.atleast_1d(np.asarray(T0, float))
        gam = np.linspace(gamma_min, gamma_max, n)
        if adiabatic == "auto":
            adiabatic = gdot >= gdot_adiabatic
        T = T0.copy(); taus, Ts = [], []
        for i, g in enumerate(gam):
            r = self.solve(gdot, T, g, details=False)
            taus.append(r["tau"]); Ts.append(T.copy())
            if adiabatic and i < n - 1:
                T = T + self.mat.taylor_quinney * np.nan_to_num(r["tau"]) * (gam[i + 1] - g) / (
                    self.mat.density * self.mat.cp)
        return gam, np.array(taus).T, np.array(Ts).T


# =============================================================================
# Creep: diffusional flow (independent carrier) and stress-controlled solution
# =============================================================================
class DiffusionalCreep:
    """Nabarro-Herring + Coble flow (Frost & Ashby form), an independent strain
    carrier whose rate ADDS to the dislocation rate:

        gdot_diff = A_d * Omega * tau * D_eff / (k T d^2),
        D_eff     = D_v + (pi delta / d) D_b
    """

    def __init__(self, d, D_v, D_b=None, delta=None, A_d=42.0):
        self.d, self.D_v, self.D_b, self.delta, self.A_d = d, D_v, D_b, delta, A_d

    def D_eff(self, T, mat):
        S = {"T": np.asarray(T, float)}
        D = self.D_v(S, mat)
        if self.D_b is not None:
            delta = self.delta if self.delta is not None else 2 * mat.b
            D = D + np.pi * delta / self.d * self.D_b(S, mat)
        return D

    def rate(self, tau, T, mat):
        return self.A_d * mat.Omega * tau * self.D_eff(T, mat) / (KB * np.asarray(T) * self.d**2)


class CreepModel:
    """Dislocation model (solved in rate control) + optional diffusional flow,
    evaluated in STRESS control, as for creep tests and deformation maps.

    For each T the dislocation stress tau(gdot) is tabulated on a log grid and
    inverted.  Where DSA makes tau(gdot) non-monotonic the lowest-rate branch
    that reaches the stress is used (monotone envelope).
    """

    def __init__(self, disl, diffusional=None, lg_min=-14.0, lg_max=6.0, n_rates=241):
        self.disl, self.diff = disl, diffusional
        self.lg = np.linspace(lg_min, lg_max, n_rates)
        self.mat = disl.mat

    def table(self, T):
        T = np.atleast_1d(np.asarray(T, float))
        GG, TT = np.meshgrid(10.0**self.lg, T)                 # (nT, nrate)
        r = self.disl.solve(GG, TT, details=False)
        return GG, TT, r["tau"]

    def disl_rate(self, tau, T, table=None):
        """gdot_disl at shear stress tau (array broadcast with T rows)."""
        T = np.atleast_1d(np.asarray(T, float))
        GG, TT, tab = table if table is not None else self.table(T)
        tau = np.atleast_2d(np.asarray(tau, float))
        if tau.shape[0] != len(T):
            tau = np.broadcast_to(tau, (len(T), tau.shape[-1]))
        out = np.empty(tau.shape)
        for i in range(len(T)):
            t = np.nan_to_num(tab[i], nan=np.inf)
            env = np.maximum.accumulate(t) + 1e-9 * np.arange(len(t))
            with np.errstate(divide="ignore"):
                lg = np.interp(np.log(tau[i]), np.log(np.maximum(env, 1e-30)), self.lg,
                               left=-np.inf, right=np.inf)
            out[i] = 10.0**lg
        return out

    def rate(self, tau, T):
        """Total shear rate and its parts at stress tau and temperature T."""
        T = np.atleast_1d(np.asarray(T, float))
        tab = self.table(T)
        gd = self.disl_rate(tau, T, tab)
        tau2 = np.broadcast_to(np.atleast_2d(tau), gd.shape)
        gf = self.diff.rate(tau2, T[:, None], self.mat) if self.diff is not None else np.zeros_like(gd)
        return dict(total=gd + gf, disl=gd, diff=gf)

    def mechanism_map(self, T, tau_over_mu, n_diag=None):
        """Rate contours and controlling-process field on a (T, tau/mu) grid.

        field codes: 0 kink-pair glide (tau*/tau > 0.5)
                     1 obstacle-limited plateau (|Lambda| < 1, athermal dominated)
                     2 dynamic strain aging (Lambda < -1, m < 0); in stress control
                       this branch is unstable, so it is returned separately as
                       'dsa_band' (stress interval per T) rather than as a field
                     3 climb / recovery controlled creep (Lambda > 1)
                     4 diffusional flow (gdot_diff > gdot_disl)
        """
        T = np.asarray(T, float)
        mu = self.mat.mu(T)[:, None]
        tau = tau_over_mu[None, :] * mu
        parts = self.rate(tau, T)
        # diagnostics on the rate grid, mapped to stress through the table
        GG, TT = np.meshgrid(10.0**self.lg, T)
        d = self.disl.srs(GG, TT)
        field = np.zeros(tau.shape, int)
        lg_d = np.log10(np.clip(parts["disl"], 1e-300, None))
        for i in range(len(T)):
            j = np.clip(np.searchsorted(self.lg, lg_d[i]), 0, len(self.lg) - 1)
            frac_th = d["tau_star"][i, j] / d["tau"][i, j]
            lam = d["Lambda"][i, j]
            f = np.where(frac_th > 0.5, 0, 1)
            f = np.where((frac_th <= 0.5) & (lam < -1), 2, f)
            f = np.where((frac_th <= 0.5) & (lam > 1), 3, f)
            field[i] = f
        field = np.where(parts["diff"] > parts["disl"], 4, field)
        field = np.where(np.isinf(parts["disl"]) | (lg_d > self.lg[-1] - 1.0), -1, field)  # beyond table
        # DSA band: stresses reached on the m<0 branch (unstable in stress control)
        dsa = d["Lambda"] < -1
        band_lo = np.where(dsa.any(1), np.nanmin(np.where(dsa, d["tau"], np.inf), 1), np.nan)
        band_hi = np.where(dsa.any(1), np.nanmax(np.where(dsa, d["tau"], -np.inf), 1), np.nan)
        rate = np.where(np.isinf(parts["disl"]), np.nan, parts["total"])
        return dict(T=T, tau_over_mu=tau_over_mu, tau=tau, rate=rate,
                    disl=parts["disl"], diff=parts["diff"], field=field,
                    dsa_band=(band_lo / mu[:, 0], band_hi / mu[:, 0]))
