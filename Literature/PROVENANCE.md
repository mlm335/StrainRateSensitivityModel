# Where every term comes from

Equation numbers are those of `paper/main.tex`. "[in folder]" means the PDF is here
under its `year - author - title` name; "[to fetch]" means it is cited but not yet
downloaded.

---

## Eq. (1) Orowan relation  `gdot = rho_m b v_bar`

| component | source |
|---|---|
| the relation itself | Orowan, *Problems of plastic gliding*, Proc. Phys. Soc. 52, 8 (1940) [in folder] |
| dislocation-density picture of slip it rests on | Taylor, Proc. R. Soc. A 145, 362 (1934) [in folder] |

## Eqs. (2)-(3) Travel time plus waiting time

| component | source |
|---|---|
| `v_bar = lambda/(t_r + t_w)`, obstacle-by-obstacle waiting | Kocks, Argon & Ashby, *Thermodynamics and Kinetics of Slip*, Prog. Mater. Sci. 19 (1975) [to fetch] |
| alternative escape routes at one obstacle add their rates (`nu_j = sum_k nu_jk`) | Kocks, Argon & Ashby (1975) [to fetch]; used in this form by Rosler & Arzt (1990) [in folder] |
| assembling both into one `1/gdot` balance | this work; no single prior paper writes it this way |

## Eq. (4) Stress superposition  `tau = tau_a + tau* + sum_g tau_g`

| component | source |
|---|---|
| adding an athermal and a thermal stress | Seeger, *Handbuch der Physik* VII/2 (1958) [to fetch] |
| superposition of dense-weak and sparse-strong obstacle populations | Kocks, Argon & Ashby (1975) [to fetch] |
| same split used in a BCC constitutive law | Cheng, Nemat-Nasser & Guo, Mech. Mater. 33, 603 (2001) [in folder] |
| adding internal-stress terms in a state-variable law | Estrin & McCormick, Acta Metall. 39, 2977 (1991) [to fetch] |

## Eq. (5) Kink-pair glide velocity `v_g`

| component | source |
|---|---|
| the complete law as used here | Po, Cui, Rivera, Cereceda, Swinburne, Marian & Ghoniem, *A phenomenological dislocation mobility law for bcc metals*, Acta Mater. 119, 123 (2016) [in folder] |
| `tau* b / B_d` phonon/kink drag prefactor | Hirth & Lothe, *Theory of Dislocations*, 2nd ed. [to fetch] |
| `(p,q)` barrier profile `dH0[1-(tau*/tau_P)^p]^q` | Kocks, Argon & Ashby (1975) [to fetch] |
| factor 1/2 on `dH0` for long screw segments (rate = sqrt(nucleation x migration)) | Dorn & Rajnak, *Nucleation of kink pairs and the Peierls mechanism*, 1964 [in folder] |
| kink-pair nucleation as the BCC rate-controlling step | Dorn & Rajnak (1964) [in folder]; Hirth & Lothe [to fetch] |
| `T/T_0` entropic term that removes the barrier at the athermal transition | Po et al. (2016) [in folder] |
| independent check of `tau_P` and of double-kink nucleation in BCC | Ghafarollahi & Curtin, Acta Mater. 196, 635 (2020) [in folder] |
| `(p,q)` values for a BCC metal fitted to flow-stress data | Cheng, Nemat-Nasser & Guo (2001) [in folder] |

## Eq. (6) Athermal stress `tau_a = alpha mu b sqrt(C/C0) (sqrt(rho_f) + Lambda_b)`

| component | source |
|---|---|
| Taylor forest stress `alpha mu b sqrt(rho)` | Taylor (1934) [in folder] |
| `alpha` values and the modern form of the Taylor law | Kocks & Mecking, Prog. Mater. Sci. 48, 171 (2003) [in folder] |
| boundary contribution `Lambda_b` (inverse boundary spacing) | Estrin & Mecking, *A unified phenomenological description of work hardening and creep*, Acta Metall. 32, 57 (1984) [in folder]; Estrin (1996) [to fetch] |
| boundary / sub-boundary dislocation storage as a separate density | Nes, Prog. Mater. Sci. 41, 129 (1998) [in folder] |
| `sqrt(C/C0)` solute strengthening of the junctions | Cheng, Nemat-Nasser & Guo (2001) [in folder]; Mulford & Kocks, Acta Metall. 27, 1125 (1979) [in folder] |
| lath/block/PAG boundaries as the strengthening hierarchy in RAFM steel | Mahler, Po, Cui, Ghoniem & Aktaa, Nucl. Mater. Energy 26, 100814 (2021) [in folder] |

## Eq. (7) Density kinetics `drho/dgamma = k1(sqrt(rho)+Lambda_b) - k2 rho - K D rho^2 / gdot`

| component | source |
|---|---|
| `k1 sqrt(rho)` athermal storage, `k2 rho` dynamic recovery | Kocks, *Laws for work-hardening and low-temperature creep*, J. Eng. Mater. Technol. 98, 76 (1976) [in folder] |
| the same, with the temperature and rate dependence of `k2` | Mecking & Kocks, *Kinetics of flow and strain-hardening*, Acta Metall. 29, 1865 (1981) [in folder] |
| `k2 = k2_0 exp(q2/kT)`, and the whole one-parameter framework | Kocks & Mecking, Prog. Mater. Sci. 48, 171 (2003) [in folder] |
| `k1 Lambda_b` storage at boundaries | Estrin & Mecking (1984) [in folder]; Estrin (1996) [to fetch] |
| `-K D(T) rho^2` static (climb-controlled) recovery: second-order annihilation of networks and dipoles | Sandstrom & Lagneborg, *A model for hot working occurring by recrystallization*, Acta Metall. 23, 387 (1975) [in folder] |
| the same term in a work-hardening framework | Nes, Prog. Mater. Sci. 41, 129 (1998) [in folder] |
| later use of `-K D rho^2` for creep of ferritic steels | Sandstrom & Hallgren, J. Nucl. Mater. 422, 51 (2012) [to fetch] |
| mechanistic (DDD) justification that static recovery is climb-controlled and second order | Kohnert & Capolungo, npj Comput. Mater. 8, 92 (2022) [in folder] |
| resulting network creep law `tau_a ~ (gdot/D)^(1/3)` | Blum, Eisenlohr & Breutinger, *Understanding creep*, Metall. Mater. Trans. A 33, 291 (2002) [in folder]; Poirier, *Creep of Crystals* (1985) [to fetch] |

**This answers the question you asked earlier.** The `K D rho^2 / gdot` term is *not*
dynamic recovery. Dynamic recovery is `k2 rho`: it is driven by the deformation itself
(cross slip and glide annihilation) and acts per unit strain. The `K D rho^2` term is
*static* recovery: diffusion-controlled climb annihilation that runs in wall-clock time
whether or not the crystal deforms, which is why dividing by `gdot` is what puts it on a
strain basis. Sandstrom & Lagneborg (1975) is the paper that introduces it in this form;
Kohnert & Capolungo (2022) is the modern mechanistic confirmation.

## Eq. (8) Solute aging `C/C0 = 1 + sum_i C_i {1 - exp[-(t_a nu_i e^{-Q_i/kT})^alpha_i]}`

| component | source |
|---|---|
| `t_a = rho_m b lambda_f / gdot`, the waiting time at forest junctions | McCormick, *Theory of flow localisation due to dynamic strain ageing*, Acta Metall. 36, 3061 (1988) [in folder] |
| the same waiting time inside a density-based PLC model | Kubin & Estrin, *Evolution of dislocation densities and the critical conditions for the PLC effect*, Acta Metall. Mater. 38, 697 (1990) [in folder] |
| `t^(2/3)` pipe-diffusion solute accumulation | Cottrell & Bilby, *Dislocation theory of yielding and strain ageing of iron*, Proc. Phys. Soc. A 62, 49 (1949) [in folder] |
| saturating form `1 - exp[-(...)^alpha]` that fixes Cottrell-Bilby at long times | Louat, *On the theory of the Portevin-Le Chatelier effect*, Scr. Metall. 15, 1167 (1981) [in folder] |
| `sqrt(C/C0)` strengthening, and DSA as extra forest strength | Mulford & Kocks, *New observations on the mechanisms of dynamic strain aging and of jerky flow*, Acta Metall. 27, 1125 (1979) [in folder] |
| aging term in a BCC constitutive model of this type | Cheng, Nemat-Nasser & Guo (2001) [in folder] |
| DSA data for Zr that the Zr aging parameters are checked against | Lee et al., J. Alloys Compd. 428, 99 (2007) [in folder]; Derep et al., Acta Metall. 28, 607 (1980) [in folder] |

## Eq. (9) Thermal detachment from a particle `nu_th`

| component | source |
|---|---|
| `(p,q)` activation form with `tau_p/tau_hat` | Kocks, Argon & Ashby (1975) [to fetch] |
| detachment from dispersoids as the rate-controlling event, and the detachment parameter | Rosler & Arzt, *A new model-based creep equation for dispersion strengthened materials*, Acta Metall. Mater. 38, 671 (1990) [in folder] |
| threshold stress `tau_hat` for a dispersion, and the Orowan limit | Arzt & Ashby, *Threshold stresses in materials containing dispersed particles*, Scr. Metall. 16, 1285 (1982) [in folder] |

## Eq. (10) Local climb over a particle `nu_cl = v_c / h_eff`

| component | source |
|---|---|
| climb velocity `v_c = 2 pi D / (b ln(R/r_c)) [exp(tau Omega/kT) - 1]` | Hirth & Lothe, *Theory of Dislocations* [to fetch] |
| local (as opposed to general) climb over particles | Brown & Ham, in *Strengthening Methods in Crystals* (1971) [to fetch] |
| residual height `h_eff = h(1 - tau_p/tau_hat)` shrinking as the line is pressed on | Arzt & Ashby (1982) [in folder]; Rosler & Arzt (1990) [in folder] |
| climb as the diffusion-controlled bypass in creep | Blum, Eisenlohr & Breutinger (2002) [in folder] |

## Eq. (11) Particle stress `gdot = rho_m b lambda_p (nu_th + nu_cl)`

| component | source |
|---|---|
| parallel escape routes adding their rates at one obstacle | Kocks, Argon & Ashby (1975) [to fetch] |
| the same competition (detachment vs climb) setting a threshold creep law | Rosler & Arzt (1990) [in folder]; Arzt & Ashby (1982) [in folder] |
| particle spacing `lambda_p` from the measured populations | Mahler et al. (2021) [in folder] |

## Diffusional flow (maps only)

| component | source |
|---|---|
| lattice (Nabarro-Herring) flow | Nabarro, Bristol Conf. on Strength of Solids (1948) [to fetch]; Herring, *Diffusional viscosity of a polycrystalline solid*, J. Appl. Phys. 21, 437 (1950) [in folder] |
| boundary (Coble) flow | Coble, *A model for boundary diffusion controlled creep*, J. Appl. Phys. 34, 1679 (1963) [in folder] |
| the combined rate expression and its coefficients | Frost & Ashby, *Deformation-Mechanism Maps* (1982) [to fetch] |

## Eq. (12) Rate-sensitivity decomposition `1/m = tau A/(1+Lambda)`, `Lambda = AB`

| component | source |
|---|---|
| the decomposition itself | this work (derived from Eq. 4; no prior source) |
| definitions of `m` and `V*` and their interpretation as activation parameters | Kocks, Argon & Ashby (1975) [to fetch] |
| `m < 0` as the precondition for serrated flow | McCormick (1988) [in folder]; Louat (1981) [in folder] |
| DSA producing negative rate sensitivity without a new mechanism | Mulford & Kocks (1979) [in folder] |
| critical conditions for the PLC effect from density evolution | Kubin & Estrin (1990) [in folder] |

## Deformation-mechanism maps

Frost & Ashby, *Deformation-Mechanism Maps* (1982) [to fetch] -- the axes, the field
construction, and the diffusional-flow coefficients. The field labels here come from the
rule's own diagnostics rather than from separate rate equations, which is the difference
from Frost & Ashby.

---

# Data, by figure

| data | source |
|---|---|
| W bulk flow stress vs T at three rates | Miranda et al., Int. J. Refract. Met. Hard Mater. 139, 107811 (2026) [in folder] |
| W `m`, `V*` -- cold-rolled sheet | Wei & Kecskes, Mater. Sci. Eng. A 491, 62 (2008) [in folder] |
| W `m`, `V*` -- high-temperature nanoindentation | Kappacher et al., Mater. Des. 189, 108499 (2020) [in folder] |
| W -- K-doped and pure cold-rolled sheet, BDTT | Lied et al., J. Nucl. Mater. 544, 152664 (2021) [in folder] |
| W -- kink-pair theory for comparison | Ghafarollahi & Curtin, Acta Mater. 196, 635 (2020) [in folder] |
| W -- remaining micropillar / IGW / arc-melted sets | Sasaki et al. (2015), Schmalbach (2022 thesis) [to fetch] |
| Zr prismatic CRSS at two rates, activation parameters | Derep et al., Acta Metall. 28, 607 (1980) [in folder] |
| Zr prismatic slip strength, micro-cantilever | Gong et al., Acta Mater. 96, 249 (2015) [in folder] |
| Zr DSA, `m(T)` | Lee et al., J. Alloys Compd. 428, 99 (2007) [in folder]; Lee (1972) [to fetch] |
| Eurofer-97 tensile at three rates, 300-873 K | Vanaja et al., J. Nucl. Mater. 424, 116 (2012) [in folder] |
| Eurofer-97 tensile, irradiated and unirradiated | Materna-Morris et al., J. Nucl. Mater. 442, S62 (2013) [in folder] |
| Eurofer-97, further mechanical data | Luzginova et al., J. Nucl. Mater. 409, 153 (2011) [in folder] |
| Eurofer-97 microstructure and particle populations | Mahler et al., Nucl. Mater. Energy 26, 100814 (2021) [in folder] |
| Eurofer-97 minimum creep rate, 450-650 C | Fernandez et al., Fusion Eng. Des. 75-79, 1003 (2005) [in folder] |
| Eurofer-97 thermal creep, 550 C / 300 MPa point | Yu, Nita & Baluc, Fusion Eng. Des. 75-79, 1037 (2005) [in folder] |
| IN-RAFM steel, shown for comparison only | Shah, Sunil & Sarkar, Metall. Mater. Trans. A 49, 2644 (2018) [in folder] |

---

# Still to fetch

**Load-bearing (nothing in the folder substitutes for these):**

* Kocks, Argon & Ashby, *Thermodynamics and Kinetics of Slip*, Prog. Mater. Sci. 19
  (1975) -- the `(p,q)` barrier profile in Eqs. (5) and (9), the dense-weak /
  sparse-strong superposition rule behind Eq. (4), the travel-plus-waiting
  bookkeeping of Eqs. (2)-(3), and the definitions of `m` and `V*`.  It is cited
  more times in `main.tex` than any other reference.
* Hirth & Lothe, *Theory of Dislocations*, 2nd ed. -- the drag prefactor in Eq. (5)
  and the climb velocity in Eq. (10).
* Frost & Ashby, *Deformation-Mechanism Maps* (1982) -- the diffusional-flow
  coefficients and the map construction.  Free PDF posted by Dartmouth; check
  `engineering.dartmouth.edu`.

**Secondary (a paper in the folder now covers the same point):**

* Seeger, *Handbuch der Physik* VII/2 (1958) -- stress superposition; Kocks, Argon &
  Ashby and Cheng et al. (2001, in folder) cover it.
* Estrin (1996), *Unified Constitutive Laws of Plastic Deformation*, ch. 2 -- boundary
  storage; Estrin & Mecking (1984, in folder) is the primary reference.  (This is the
  chapter whose Fig. 5 shows the Al 1100 stress-strain curves.)
* Estrin & McCormick, Acta Metall. Mater. 39, 2977 (1991) -- transient DSA.
* Brown & Ham (1971) -- local climb; Arzt & Ashby (1982, in folder) covers it.
* Nabarro (1948) -- lattice diffusional flow; Herring (1950, in folder) covers it.
* Poirier, *Creep of Crystals* (1985) -- network creep; Blum et al. (2002, in folder)
  covers it.
* Sandstrom (1977); Nes (1995); Sandstrom & Hallgren (2012); Roters, Raabe &
  Gottstein (2000); Kuhlmann, Masing & Raffelsieper (1949) -- all recovery-term
  supporting citations; Sandstrom & Lagneborg (1975), Nes (1998) and Kohnert &
  Capolungo (2022) are in the folder.

**Data sets still missing:**

* Sasaki et al., J. Nucl. Mater. 461, 357 (2015) -- K-doped W tensile.
* Schmalbach, PhD thesis, Univ. of Minnesota (2022) -- W rate effects.
* Lee (1972) -- Zr rate sensitivity; journal and volume still unknown.
* The remaining W `m` / `V*` sets behind `W_M_DATA` and `W_V_DATA` in
  `parameters.py` (micropillar, indentation literature, IGW, W plate, arc-melted W).
  These are the `\todo{Wdata}` citation in `main.tex`; the Abaalkhail manuscript
  cites several of them and is the quickest place to resolve them.
