"""
parameters.py -- model parameters and experimental data for W, Zr and Eurofer-97.

Each material is one flat dict.  Terms a material does not have are set to zero:
    C_sol   = 0   no solute aging (W)
    tau_ppt = 0   no particles     (W, Zr)
    Lambda_b= 0   no boundary storage (W, Zr)

Conventions
-----------
Everything in the model is in shear:  gdot = eps_dot / M_s  and  tau = M_s sigma.
Experimental tables below are stored as measured (axial sigma in MPa for W and
Eurofer-97, already-resolved tau in MPa for Zr) and converted where they are used.
"""
import numpy as np

# =====================================================================
# Tungsten: kink-pair glide against a Kocks-Mecking forest
# =====================================================================
W = dict(
    name="W",
    # lattice and elasticity
    b=0.2722e-9, 
    mu0=161.0e9, 
    mu_slope=0.0, 
    Tm=3695.0, 
    Omega=1.585e-29, 
    M_s=0.333,
    # kink-pair mobility
    dH0_eV=2.52, 
    tau_P=0.70e9, 
    p=0.6, 
    q=1.4, 
    T0=0.9 * 3695.0,
    B_kink=8.3e-5 * 25.0 / (2.0 * (2.0 * np.sqrt(2.0) / 3.0)), 
    B_free=0.0,
    # forest: Taylor stress and Kocks-Mecking kinetics
    alpha=0.549, 
    k1=5.5e9, 
    k2_0=891.0, 
    q2_eV=-0.011, 
    beta=1.0, 
    f_m=1.0, 
    Lambda_b=0.0,
    # climb-controlled recovery (lattice self-diffusion; not calibrated to creep)
    K_cl=1e9, 
    D0_rec=4e-5, 
    QD_rec_eV=6.0,
    # solute aging: none
    C_sol=0.0, 
    Q_sol_eV=0.0, 
    nu_sol=1e13, 
    alpha_sol=1 / 3,
    # particles: none
    tau_ppt=0.0, 
    F_det_eV=0.0, 
    p_det=1.0, 
    q_det=1.0, 
    nu0_det=1e11,
    inv_lambda2_ppt=1.0, 
    h_ppt=1e-9, 
    lnR_rc=5.0, 
    D0_climb=0.0, 
    QD_climb_eV=0.0,
    # diffusional flow (mechanism map only)
    d_grain=1e-4, 
    Db0=4e-5, 
    Qb_eV=3.9,
    rates=np.array([6.4e-4, 6.4e-3, 6.4e-2]),
)

# Bulk ITER-grade W (ATW), axial flow stress at eps_p = 0.02: Miranda et al.,
# Int. J. Refract. Met. Hard Mater. 139 (2026) 107811.  (T [K], eps_dot [1/s], sigma [MPa])
W_FLOW = np.array([(523, 6.4e-4, 668.), (523, 6.4e-3, 743.), (673, 6.4e-4, 555.),
                   (673, 6.4e-3, 569.), (673, 6.4e-2, 615.), (823, 6.4e-4, 515.),
                   (823, 6.4e-3, 527.), (823, 6.4e-2, 539.)])
# strain-rate sensitivity m(T) and activation volume V*/b^3 from the literature
W_M_DATA = {
    "ATW Miranda 2026": (np.array([523, 673, 823]), np.array([0.0400, 0.0170, 0.0090])),
    "Nanoindentation": (np.array([300, 380, 480, 550, 670]), np.array([0.030, 0.035, 0.045, 0.050, 0.028])),
    "Indentation literature": (np.array([25, 100, 200, 300, 400, 500, 600, 700, 800]) + 273.15,
                               np.array([0.023, 0.030, 0.037, 0.038, 0.030, 0.020, 0.018, 0.011, 0.010])),
    "IGW": (np.array([523.9, 672.1, 820.3]), np.array([0.0488, 0.0140, 0.0067])),
    "W plate": (np.array([475.9, 573.5, 672.3]), np.array([0.0350, 0.0211, 0.0149])),
}
W_V_DATA = {
    "ATW Miranda 2026": (np.array([673, 823]), np.array([45.7, 117.2])),
    "Nanoindentation": (np.array([300, 380, 480, 550, 670]), np.array([6, 8, 10, 15, 60])),
    "Indentation literature": (np.array([25, 100, 200, 300, 400, 500, 600, 700, 800]) + 273.15,
                               np.array([8, 8.5, 10, 16, 29, 52, 73, 108, 125])),
    "IGW": (np.array([522.2, 672.7, 822.2]), np.array([12.9, 65.3, 179.5])),
    "W plate": (np.array([473.1, 573.2, 672.8]), np.array([11.8, 28.1, 55.4])),
    "Arc-melted W": (np.array([314.4, 549.5, 679.8]), np.array([9.8, 19.2, 84.8])),
}

# =====================================================================
# Zirconium (prismatic slip): W level + oxygen aging + climb recovery
# =====================================================================
ZR = dict(
    name="Zr",
    b=0.3233e-9, 
    mu0=33.0e9, 
    mu_slope=0.0, 
    Tm=2128.0, 
    Omega=2.33e-29, 
    M_s=0.45,
    dH0_eV=3.75, 
    tau_P=300e6, 
    p=0.86, 
    q=1.69, 
    T0=0.7 * 2128.0,
    B_kink=50e-6 * 25.0 / (2.0 * (np.sqrt(8.0 / 3.0) / 2.0)), 
    B_free=50e-6,
    alpha=0.36, 
    k1=5.5e9, 
    k2_0=500.0, 
    q2_eV=-0.005, 
    beta=1.0, 
    f_m=1.0, 
    Lambda_b=0.0,
    K_cl=1e9, 
    D0_rec=1e-4, 
    QD_rec_eV=3.2,
    # oxygen atmospheres (effective Q = migration - binding)
    C_sol=0.30, 
    Q_sol_eV=1.80, 
    nu_sol=1e13, 
    alpha_sol=2 / 3,
    tau_ppt=0.0, 
    F_det_eV=0.0, 
    p_det=1.0, 
    q_det=1.0, 
    nu0_det=1e11,
    inv_lambda2_ppt=1.0, 
    h_ppt=1e-9, 
    lnR_rc=5.0, 
    D0_climb=0.0, 
    QD_climb_eV=0.0,
    d_grain=2e-5, 
    Db0=1e-4, 
    Qb_eV=2.1,
    rates=np.array([1e-5, 1e-4, 1e-3]),
)

# critical resolved shear stress (MPa) at two axial rates: Derep et al.
ZR_CRSS = {3.3e-3: (np.array([470, 575, 660, 720, 790, 870.]), np.array([285, 195, 155, 150, 151, 103.]) / 3),
           3.3e-5: (np.array([470, 570, 605, 675, 705, 785.]), np.array([223, 168, 160, 165, 162, 106.]) / 3)}
ZR_M_DATA = {
    "Cantilever prism slip": (np.array([20, 150, 300]) + 273, np.array([0.0155, 0.0475, 0.0378])),
    "Cantilever basal slip": (np.array([20, 150, 300]) + 273, np.array([0.0185, 0.0480, 0.0288])),
    "Lee 1972 slow": (np.array([293, 383, 463, 563, 693]), np.array([0.024, 0.038, 0.035, 0.020, 0.087])),
    "Lee 1972 fast": (np.array([563, 673]), np.array([0.020, 0.041])),
    "Lee 2007 1.33e-4": (np.array([300, 500, 590, 630, 680, 730]),
                         np.array([0.007, 0.045, 0.015, 0.011, 0.046, 0.080])),
    "Lee 2007 3.33e-3": (np.array([300, 540, 580, 630, 720]), np.array([0.013, 0.050, 0.031, 0.007, 0.020])),
}
ZR_V_DATA = {
    "Cantilever prism slip": (np.array([20, 150, 300]) + 273, np.array([26.5, 18.0, 40.5])),
    "Cantilever basal slip": (np.array([20, 150, 300]) + 273, np.array([19.5, 15.8, 52.5])),
    "Derep (a)": (np.array([80, 190, 240, 270, 330, 590, 790]), np.array([16, 29, 42, 47, 47, 168, 130])),
    "Derep (b)": (np.array([100, 150, 200, 250, 300, 400, 500, 550, 580, 600, 620, 635, 780, 820, 860, 900]),
                  np.array([18, 24, 32, 45, 50, 42, 60, 90, 120, 165, 220, 270, 140, 125, 115, 108])),
}

# =====================================================================
# Eurofer-97: Zr level + boundary storage + particles (detachment || climb)
# =====================================================================
_A_FE = 0.287e-9
# boundaries: lath 0.5 um, block 3.1 um, PAG 21 um, MX 1.5 um, M23C6 0.456 um
EU_LAMBDA_B = 1 / 5e-7 + 0.3 / 3.1e-6 + 0.1 / 2.1e-5 + 0.3 / 1.5e-6 + 0.2 / 4.56e-7
# particle areal density: MX (2.4e19 m^-3, 18.6 nm) + M23C6 (5.7e19 m^-3, 67.6 nm)
EU_NA_PPT = 2.4e19 * 18.6e-9 + 5.7e19 * 67.6e-9

EU = dict(
    name="Eurofer97",
    b=_A_FE * np.sqrt(3) / 2, mu0=80e9, mu_slope=0.35, mu_floor=0.45, Tm=1811.0,
    Omega=_A_FE**3 / 2, M_s=0.333,
    # kink-pair mobility (Abaalkhail et al. 2026 / Mahler et al. 2021)
    dH0_eV=2.17, tau_P=380e6, p=0.6, q=1.95, T0=0.845 * 1811.0,
    B_kink=6.6e-5 * 25.0 / (2.0 * (2.0 * np.sqrt(2.0) / 3.0)), B_free=0.0,
    alpha=0.1283, k1=7.8e10, k2_0=6750.0, q2_eV=-0.007, beta=0.1, f_m=1.0,
    Lambda_b=EU_LAMBDA_B,
    K_cl=10**5.9296, D0_rec=2.0e-4, QD_rec_eV=2.772,
    # C, N atmospheres on forest junctions (DSA)
    C_sol=6.000, Q_sol_eV=1.625, nu_sol=1e13, alpha_sol=0.288,
    # MX and M23C6 particles: thermal detachment || glide-assisted local climb
    tau_ppt=254.0e6, F_det_eV=3.556, p_det=1.0, q_det=1.0, nu0_det=1e11,
    inv_lambda2_ppt=EU_NA_PPT, h_ppt=1e-9, lnR_rc=5.0,
    D0_climb=10**-0.9661, QD_climb_eV=3.226,
    d_grain=2.1e-5, Db0=2.0e-4, Qb_eV=1.8,
    rates=np.array([3e-5, 3e-4, 3e-3]),
)

# axial flow stress (MPa) at eps_p = 0.02, Vanaja et al., J. Nucl. Mater. 424 (2012) 116
EU_T = np.array([300, 373, 423, 473, 523, 573, 623, 673, 723, 773, 823, 873.])
EU_SIGMA = np.array([[655, 720, 750], [620, 665, 670], [np.nan, 660, np.nan], [565, 585, 585],
                     [560, 600, 580], [520, 570, 535], [515, 540, 530], [505, 515, 500],
                     [450, 465, 480], [400, 440, 470], [350, 400, 440], [235, 305, 365.]])
# m and V*/b^3 reported by Vanaja et al. (digitized from Abaalkhail et al. 2026, Fig. 6)
EU_M_DATA = np.array([(300, .0296), (373, .0154), (473, .0057), (523, .0000), (573, .0062), (623, .0103),
                      (673, .0045), (723, .0210), (773, .0409), (823, .0623), (873, .0960)])
EU_V_DATA = np.array([(300, 12.5), (373, 33.1), (473, 129.7), (573, 157.5), (623, 101.4), (673, 277.8),
                      (723, 65.3), (773, 39.2), (823, 29.2), (873, 27.0)])
# Shah et al., Metall. Mater. Trans. A 49 (2018) 2644 -- IN-RAFM steel (comparison only)
SHAH = np.array([(773., 0.0179, 263.), (823., 0.0391, 155.), (873., 0.0600, 135.)])
# Materna-Morris et al., J. Nucl. Mater. 442 (2013) S62 -- unirradiated (rate not stated)
EU_MM = np.array([(300., 635.), (573., 541.), (723., 467.)])

# minimum creep rates (T [K], sigma [MPa], eps_dot_min [1/s]):
# Fernandez et al., Fusion Eng. Des. 75-79 (2005) 1003, Fig. 1 (digitized)
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
# Yu, Nita and Baluc, Fusion Eng. Des. 75-79 (2005) 1037: 550 C, 300 MPa
EU_CREEP_YU = np.array([(823., 300.0, 1.0e-5)])
