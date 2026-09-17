"""
material_cards.py
=================
Experimental data and calibrated model instances for W (BCC), Zr (HCP,
prismatic slip, kink-pair controlled) and Eurofer-97 (BCC ferritic-
martensitic), all built from the same components in unified_plasticity.py.

Conventions follow the original PRL notebooks:
  * shear rate gdot = eps_dot / M_s  (M_s = 0.333 W/Eurofer, 0.45 Zr)
  * activation volume from sigma = tau / M_s
  * CRSS comparison: W data as given (MPa), Zr and Eurofer axial data / 3
"""
import numpy as np
from unified_plasticity import *

# =============================================================================
# W  (BCC)  -- kink-pair glide + athermal forest (Kocks-Mecking)
#   Parameters = the original BCC notebook (the framework reproduces that
#   model exactly: no waiting steps, no rate-dependent state -> Lambda = 0)
# =============================================================================
W_MAT = Material("W", b=0.2722e-9, mu0=161.0e9, Tm=3695.0, Omega=1.585e-29)
W_P = dict(dH0_eV=2.52, tau_P=0.70e9, p=0.6, q=1.4, T0_frac=0.9,
           B_kink=8.3e-5 * 25.0 / (2.0 * (2.0 * np.sqrt(2.0) / 3.0)),
           alpha=0.549, k1=5.5e9, k20=891.0, q2_eV=-0.011,
           # climb-controlled recovery (lattice self-diffusion; not calibrated to creep yet)
           K_cl=1e9, D0_rec=4e-5, QD_rec_eV=6.0,
           # diffusional creep
           d_grain=1e-4, Db0=4e-5, Qb_eV=3.9)
W_M = 0.333
W_RATES = np.array([6.4e-4, 6.4e-3, 6.4e-2])


def make_W(P=W_P, recovery=True):
    return DislocationModel(
        W_MAT,
        glide=KinkPairGlide(eV(P["dH0_eV"]), P["tau_P"], P["p"], P["q"], P["B_kink"],
                            T0=P["T0_frac"] * W_MAT.Tm),
        athermal=[TaylorAthermal(P["alpha"], name="forest")],
        microstructure=[KocksMecking(P["k1"], P["k20"], q2=P["q2_eV"],
                                     K_cl=P["K_cl"] if recovery else 0.0,
                                     diffusivity=Diffusivity(P["D0_rec"], eV(P["QD_rec_eV"])))],
        name="W")


W_DATA = dict(
    # Bulk ITER-grade tungsten (ATW), axial flow stress at eps_p = 0.02 (MPa): Miranda et al.,
    # Int. J. Refract. Met. Hard Mater. 139 (2026) 107811 (digitized from Abaalkhail et al. 2026, Fig. 3a)
    flow=np.array([(523, 6.4e-4, 668.), (523, 6.4e-3, 743.), (673, 6.4e-4, 555.), (673, 6.4e-3, 569.),
                   (673, 6.4e-2, 615.), (823, 6.4e-4, 515.), (823, 6.4e-3, 527.), (823, 6.4e-2, 539.)]),
    m={"ATW Miranda 2026": (np.array([523, 673, 823]), np.array([0.0400, 0.0170, 0.0090])),
       "Micropillar W": (np.array([300, 500, 700]), np.array([0.035, 0.046, 0.025])),
       "Nanoindentation": (np.array([300, 380, 480, 550, 670]), np.array([0.03, 0.035, 0.045, 0.05, 0.028])),
       "Indentation literature": (np.array([25, 100, 200, 300, 400, 500, 600, 700, 800]) + 273.15,
                                  np.array([0.023, 0.03, 0.037, 0.038, 0.03, 0.02, 0.018, 0.011, 0.01])),
       "IGW": (np.array([523.9, 672.1, 820.3]), np.array([0.0488, 0.0140, 0.0067])),
       "W plate": (np.array([475.9, 573.5, 672.3]), np.array([0.0350, 0.0211, 0.0149]))},
    V={"ATW Miranda 2026": (np.array([673, 823]), np.array([45.7, 117.2])),
       "Micropillar W": (np.array([300, 500, 700]), np.array([11, 20, 130])),
       "Nanoindentation": (np.array([300, 380, 480, 550, 670]), np.array([6, 8, 10, 15, 60])),
       "Indentation literature": (np.array([25, 100, 200, 300, 400, 500, 600, 700, 800]) + 273.15,
                                  np.array([8, 8.5, 10, 16, 29, 52, 73, 108, 125])),
       "IGW": (np.array([522.2, 672.7, 822.2]), np.array([12.9, 65.3, 179.5])),
       "W plate": (np.array([473.1, 573.2, 672.8]), np.array([11.8, 28.1, 55.4])),
       "Arc-melted W": (np.array([314.4, 549.5, 679.8]), np.array([9.8, 19.2, 84.8]))},
)

# =============================================================================
# Zr  (HCP prismatic) -- kink-pair glide + athermal forest
#   + DSA by O atmospheres on the forest (waiting-time aging)
#   + climb-controlled dynamic recovery of the forest (glide-assisted climb)
# =============================================================================
ZR_MAT = Material("Zr", b=0.3233e-9, mu0=33.0e9, Tm=2128.0, Omega=2.33e-29)
ZR_P = dict(dH0_eV=3.75, tau_P=300e6, p=0.86, q=1.69, T0_frac=0.7,
            B_kink=50e-6 * 25.0 / (2.0 * (np.sqrt(8.0 / 3.0) / 2.0)), B_free=50e-6,
            alpha=0.36, k1=5.5e9, k20=500.0, q2_eV=-0.005,
            # O aging (effective Q = migration - binding), nu, kinetics exponent
            C_O=0.30, Q_O_eV=1.80, nu_O=1e13, alpha_O=2 / 3,
            # climb recovery: d rho/dt = -K_cl D rho^2, D = D0 exp(-Q/kT)
            K_cl=1e9, D0=1e-4, QD_eV=3.2,
            # diffusional creep
            d_grain=2e-5, Db0=1e-4, Qb_eV=2.1)
ZR_M = 0.45
ZR_RATES = np.array([1e-5, 1e-4, 1e-3])


def make_Zr(P=ZR_P, aging=True, climb=True):
    D = Diffusivity(P["D0"], eV(P["QD_eV"]))
    return DislocationModel(
        ZR_MAT,
        glide=KinkPairGlide(eV(P["dH0_eV"]), P["tau_P"], P["p"], P["q"], P["B_kink"],
                            T0=P["T0_frac"] * ZR_MAT.Tm, B_free=P["B_free"], drag_switch=True),
        athermal=[TaylorAthermal(P["alpha"], name="forest")],
        microstructure=[KocksMecking(P["k1"], P["k20"], q2=P["q2_eV"],
                                     K_cl=P["K_cl"] if climb else 0.0, diffusivity=D)],
        aging=(SoluteAging([dict(C=P["C_O"], Q=eV(P["Q_O_eV"]), nu=P["nu_O"], alpha=P["alpha_O"])],
                           targets=("forest",), kinds=("tau_hat",)) if aging else None),
        name="Zr")


ZR_DATA = dict(
    crss={3.3e-3: (np.array([470, 575, 660, 720, 790, 870]), np.array([285, 195, 155, 150, 151, 103]) / 3),
          3.3e-5: (np.array([470, 570, 605, 675, 705, 785]), np.array([223, 168, 160, 165, 162, 106]) / 3)},
    m={"Cantilever prism slip": (np.array([20, 150, 300]) + 273, np.array([0.0155, 0.0475, 0.0378])),
       "Cantilever basal slip": (np.array([20, 150, 300]) + 273, np.array([0.0185, 0.0480, 0.0288])),
       "Lee 1972 slow": (np.array([293, 383, 463, 563, 693]), np.array([0.024, 0.038, 0.035, 0.020, 0.087])),
       "Lee 1972 fast": (np.array([563, 673]), np.array([0.020, 0.041])),
       "Lee 2007 1.33e-4": (np.array([300, 500, 590, 630, 680, 730]),
                            np.array([0.007, 0.045, 0.015, 0.011, 0.046, 0.080])),
       "Lee 2007 3.33e-3": (np.array([300, 540, 580, 630, 720]), np.array([0.013, 0.050, 0.031, 0.007, 0.020]))},
    V={"Cantilever prism slip": (np.array([20, 150, 300]) + 273, np.array([26.5, 18.0, 40.5])),
       "Cantilever basal slip": (np.array([20, 150, 300]) + 273, np.array([19.5, 15.8, 52.5])),
       "Derep (a)": (np.array([80, 190, 240, 270, 330, 590, 790]), np.array([16, 29, 42, 47, 47, 168, 130])),
       "Derep (b)": (np.array([100, 150, 200, 250, 300, 400, 500, 550, 580, 600, 620, 635, 780, 820, 860, 900]),
                     np.array([18, 24, 32, 45, 50, 42, 60, 90, 120, 165, 220, 270, 140, 125, 115, 108]))},
)

# =============================================================================
# Eurofer-97  (BCC ferritic-martensitic)
#   kink-pair glide + athermal forest/lath (KM with boundary storage)
#   + DSA by C/N on the forest
#   + precipitates (MX, M23C6) acting simultaneously (superposed group):
#     thermally activated detachment OR glide-assisted local climb
# =============================================================================
EU_Tm = 1811.0


def _mu_eurofer(T):
    return np.maximum(80e9 * (1.0 - 0.35 * (np.asarray(T) - 300.0) / (EU_Tm - 300.0)), 0.45 * 80e9)


A_FE = 0.287e-9                       # lattice parameter of alpha-Fe
EU_MAT = Material("Eurofer97", b=A_FE * np.sqrt(3) / 2, mu0=80e9, Tm=EU_Tm, mu_model=_mu_eurofer,
                  Omega=A_FE**3 / 2, density=7760.0, cp=460.0)
# boundaries (lath, block, PAG) and particles (MX, M23C6), 1/m
EU_LAMBDA = 1 / 5e-7 + 0.3 / 3.1e-6 + 0.1 / 2.1e-5 + 0.3 / 1.5e-6 + 0.2 / 4.56e-7
# particle areal density N_A = sum N_v D (MX: 2.4e19 m^-3, 18.6 nm; M23C6: 5.7e19 m^-3, 67.6 nm)
EU_NA_PPT = 2.4e19 * 18.6e-9 + 5.7e19 * 67.6e-9
EU_P = dict(# kink-pair mobility (Abaalkhail et al. 2026 / Mahler et al. 2021)
            dH0_eV=2.17, tau_P=380e6, p=0.6, q=1.95, T0_frac=0.845,
            B_kink=6.6e-5 * 25.0 / (2.0 * (2.0 * np.sqrt(2.0) / 3.0)),
            alpha=0.1283, k1=7.8e10, k20=6750.0, q2_eV=-0.007, beta=0.1,
            # DSA (C, N): aging of forest junctions (multiplies the Taylor stress)
            C_CN=6.000, Q_CN_eV=1.625, nu_CN=1e13, alpha_CN=0.288,
            # particles (MX, M23C6): detachment || local climb (D0 lumped with h, ln R/rc)
            tau_ppt=254.0e6, F_det_eV=3.556, p_det=1.000, q_det=1.000,
            inv_lambda2_ppt=EU_NA_PPT, h_ppt=1e-9,
            D0_climb=10**-0.9661, QD_climb_eV=3.226,
            # climb-controlled recovery of the forest
            K_cl=10**5.9296, D0_rec=2.0e-4, QD_rec_eV=2.772,
            # diffusional flow (used only for the deformation-mechanism maps)
            d_grain=2.1e-5, Db0=2.0e-4, Qb_eV=1.8)
EU_M = 0.333
EU_RATES = np.array([3e-5, 3e-4, 3e-3])


def make_Eurofer(P=EU_P, aging=True, precipitates=True, climb=True, recovery=True):
    groups = []
    if precipitates:
        routes = [ThermalRoute(1e11, P.get("p_det", 1.0), P.get("q_det", 1.0), name="detachment")]
        if climb:
            routes.append(ClimbRoute(Diffusivity(P["D0_climb"], eV(P["QD_climb_eV"])), name="climb",
                                     stress_assisted_height=True))
        ppt = PrecipitateObstacle(routes, tau_hat0=P["tau_ppt"], inv_lambda2=P["inv_lambda2_ppt"],
                                  h=P["h_ppt"], g_F=eV(P["F_det_eV"]) / (EU_MAT.mu0 * EU_MAT.b**3),
                                  name="MX_M23C6")
        groups = [ResistanceGroup("precipitates", [ppt])]
    athermal = [TaylorAthermal(P["alpha"], name="forest")]
    return DislocationModel(
        EU_MAT,
        glide=KinkPairGlide(eV(P["dH0_eV"]), P["tau_P"], P["p"], P["q"], P["B_kink"],
                            T0=P["T0_frac"] * EU_Tm),
        athermal=athermal,
        microstructure=[KocksMecking(P["k1"], P["k20"], q2=P["q2_eV"], Lambda=EU_LAMBDA,
                                     beta=P["beta"], f_m=1.0,
                                     K_cl=P["K_cl"] if recovery else 0.0,
                                     diffusivity=Diffusivity(P["D0_rec"], eV(P["QD_rec_eV"])))],
        aging=(SoluteAging([dict(C=P["C_CN"], Q=eV(P["Q_CN_eV"]), nu=P["nu_CN"], alpha=P["alpha_CN"])],
                           targets=("forest",), kinds=("tau_hat",)) if aging else None),
        groups=groups, name="Eurofer97")


EU_DATA = dict(
    T=np.array([300, 373, 423, 473, 523, 573, 623, 673, 723, 773, 823, 873]),
    rates=EU_RATES,
    sigma=np.array([[655, 720, 750], [620, 665, 670], [np.nan, 660, np.nan], [565, 585, 585],
                    [560, 600, 580], [520, 570, 535], [515, 540, 530], [505, 515, 500],
                    [450, 465, 480], [400, 440, 470], [350, 400, 440], [235, 305, 365]]),
)


# Minimum creep rates of Eurofer-97 (constant load, plate and bar, air/vacuum),
# digitized from Fernandez et al., Fusion Eng. Des. 75-79 (2005) 1003, Fig. 1
# (T [K], sigma [MPa], eps_dot_min [1/s]); plus the steady-state rate quoted by
# Yu, Nita & Baluc, Fusion Eng. Des. 75-79 (2005) 1037 (550 C, 300 MPa).
EU_CREEP = np.array([
    (923, 49.7, 3.35e-10), (923, 59.8, 1.41e-09), (923, 59.6, 2.54e-09), (923, 69.7, 3.99e-09),
    (923, 69.3, 1.56e-08), (923, 79.8, 7.02e-09), (923, 79.8, 1.48e-08), (923, 89.5, 1.98e-08),
    (923, 100.1, 6.61e-08), (923, 100.1, 5.06e-08), (923, 109.6, 1.16e-07),
    (873, 100.1, 8.75e-10), (873, 100.1, 6.21e-10), (873, 100.1, 4.97e-10), (873, 109.6, 1.49e-09),
    (873, 119.4, 3.90e-09), (873, 120.5, 6.80e-09), (873, 129.5, 6.87e-09), (873, 139.9, 1.60e-08),
    (873, 139.9, 1.98e-08), (873, 149.3, 3.56e-08), (873, 159.9, 6.61e-08),
    (823, 159.4, 6.35e-10), (823, 169.8, 8.75e-10), (823, 179.4, 2.54e-09), (823, 180.3, 2.83e-09),
    (823, 189.2, 4.93e-09), (823, 189.2, 3.90e-09), (823, 189.9, 3.51e-09), (823, 199.9, 1.26e-08),
    (823, 221.6, 5.94e-08), (823, 229.7, 1.28e-07), (823, 240.2, 2.38e-07), (823, 260.8, 6.71e-07),
    (773, 199.2, 7.74e-11), (773, 240.2, 9.43e-10), (773, 241.9, 7.37e-10), (773, 249.0, 5.29e-10),
    (773, 250.7, 5.21e-09), (773, 250.7, 3.82e-09), (773, 250.7, 3.43e-09), (773, 260.8, 7.98e-09),
    (773, 267.6, 7.56e-09), (773, 289.1, 6.40e-08), (773, 300.8, 1.23e-07), (773, 302.9, 8.19e-08),
    (723, 299.3, 3.53e-10), (723, 320.5, 2.49e-09), (723, 330.0, 3.51e-09), (723, 340.4, 1.05e-08),
    (723, 362.7, 2.67e-08), (723, 361.5, 3.20e-08), (723, 369.0, 3.79e-08),
])
EU_CREEP_YU = np.array([(823, 300.0, 1.0e-5)])

# Strain-rate sensitivity and activation volume reported by Vanaja et al., J. Nucl. Mater.
# 424 (2012) 116 (digitized from Abaalkhail et al. 2026, Fig. 6): (T [K], m) and (T [K], V*/b^3)
EU_SRS_M = np.array([(300, .0296), (373, .0154), (473, .0057), (523, .0000), (573, .0062), (623, .0103),
                     (673, .0045), (723, .0210), (773, .0409), (823, .0623), (873, .0960)], float)
EU_SRS_V = np.array([(300, 12.5), (373, 33.1), (473, 129.7), (573, 157.5), (623, 101.4), (673, 277.8),
                     (723, 65.3), (773, 39.2), (823, 29.2), (873, 27.0)], float)
# Shah et al., Metall. Mater. Trans. A 49 (2018) 2644 -- IN-RAFM steel (comparison only): T, m, V*/b^3
SHAH_SRS = np.array([(773., 0.0179, 263.), (823., 0.0391, 155.), (873., 0.0600, 135.)])
# Materna-Morris et al., J. Nucl. Mater. 442 (2013) S62 -- unirradiated flow stress (MPa), rate not stated
EU_MM = np.array([(300., 635.), (573., 541.), (723., 467.)])

# =============================================================================
# Creep models (stress controlled). Diffusional flow is added only on request
# (diffusion=True), which is used for the deformation-mechanism maps.
# =============================================================================
def _diff(P, D0, Q, mat):
    return DiffusionalCreep(P["d_grain"], Diffusivity(D0, eV(Q)),
                            D_b=Diffusivity(P["Db0"], eV(P["Qb_eV"])))


def creep_W(P=W_P, diffusion=False):
    d = _diff(P, P["D0_rec"], P["QD_rec_eV"], W_MAT) if diffusion else None
    return CreepModel(make_W(P), d)


def creep_Zr(P=ZR_P, diffusion=False):
    d = _diff(P, P["D0"], P["QD_eV"], ZR_MAT) if diffusion else None
    return CreepModel(make_Zr(P), d)


def creep_Eurofer(P=EU_P, diffusion=False):
    d = _diff(P, P["D0_rec"], P["QD_rec_eV"], EU_MAT) if diffusion else None
    return CreepModel(make_Eurofer(P), d)


MATERIALS = {
    "W": dict(make=make_W, creep=creep_W, mat=W_MAT, M=W_M, rates=W_RATES),
    "Zr": dict(make=make_Zr, creep=creep_Zr, mat=ZR_MAT, M=ZR_M, rates=ZR_RATES),
    "Eurofer97": dict(make=make_Eurofer, creep=creep_Eurofer, mat=EU_MAT, M=EU_M, rates=EU_RATES),
}


def tensile_points():
    """All flow-stress data as (material, T, eps_dot, tau_MPa) rows (shear)."""
    rows = []
    for T, e, sig in W_DATA["flow"]:
        rows.append(("W", float(T), float(e), float(sig * W_M)))   # bulk ATW flow stress -> shear
    for e, (T, y) in ZR_DATA["crss"].items():
        rows += [("Zr", float(a), e, float(b)) for a, b in zip(T, y)]
    for i, T in enumerate(EU_DATA["T"]):
        for j, e in enumerate(EU_RATES):
            if np.isfinite(EU_DATA["sigma"][i, j]):
                rows.append(("Eurofer97", float(T), float(e), float(EU_DATA["sigma"][i, j] / 3)))
    return rows
