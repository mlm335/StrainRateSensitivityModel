# Unified kinetic flow rule for BCC metals -- code

Two library files and nine plotting scripts.  Everything is in shear (tau) space
and SI units; experimental axial data are converted with the Schmid factor M_s.

    physics.py      every equation of the model
    parameters.py   parameters and experimental data for W, Zr and Eurofer-97
    style.py        plot style and figure saving

    plot_W_model.py         plot_Zr_model.py         plot_Eurofer_model.py
    plot_W_creep.py         plot_Zr_creep.py         plot_Eurofer_creep.py
    plot_W_ashby.py         plot_Zr_ashby.py         plot_Eurofer_ashby.py

Run any script on its own, e.g.

    python plot_Eurofer_model.py

Each script writes one file per plot (PNG and PDF) into `figures/`.


## The model in one place

    state (set by the imposed rate and temperature, never by the stress)
        k1(sqrt(rho)+Lambda_b) - k2 rho - K D(T) rho^2/gdot = 0      -> rho
        C/C0 = 1 + C_sol [1 - exp(-(t_a nu exp(-Q/kT))^alpha)]       -> aging
    stresses
        tau_a     = alpha mu b sqrt(C/C0) (sqrt(rho_f) + Lambda_b)
        tau_star  from   gdot = rho_m b v_g(tau_star)                (kink pairs)
        tau_p     from   gdot = rho_m b lambda_p [nu_th + nu_cl]     (particles)
        tau       = tau_a + tau_star + tau_p

Switch a term off by zeroing its parameter: `C_sol = 0` (no aging), `tau_ppt = 0`
(no particles), `K_cl = 0` (no climb recovery), `Lambda_b = 0` (no boundaries).
W uses the first block only, Zr adds aging and climb recovery, Eurofer-97 adds
boundaries and particles.

## Conventions

* `gdot = eps_dot / M_s` and `tau = M_s sigma`, with M_s = 0.333 (W, Eurofer-97)
  and 0.45 (Zr).
* `m = dln(tau)/dln(gdot)` and `V* = kB T/(dsigma/dln eps_dot)`, the experimental
  definition, so V* can be compared directly with measured values.
* Diffusional (Nabarro-Herring + Coble) flow is used **only** in the Ashby maps;
  the creep scripts are dislocation creep alone.

## Fit quality with the parameters in `parameters.py`

* W: bulk flow stress (Miranda 2026, 8 points) to about 2%.
* Zr: Derep CRSS at two rates to within 7%.
* Eurofer-97: flow stress 4.5% mean error; creep rms 0.53 decades, 52/54 points
  within one decade.



## The files

- **`parameters.py`** — three dicts (`W`, `ZR`, `EU`) plus the data arrays. Nothing is computed here.
- **`physics.py`** — all equations as functions, each taking `(…, P)` where `P` is one of those dicts.
- **`style.py`** — rcParams, colors, and `save(fig, name)`.
- **`plot_*.py`** — nine scripts that call the physics and draw one figure per plot.

## Rate control (everything in the model plots)

`flow_stress(gdot, T, P)` is the entry point, and the order inside it matters because each step uses only the steps before it.

1. **`state(gdot, T, P)`** — builds the microstructure. This never sees the stress.
   - **`rho_steady`**: computes `k2 = k2_0·exp(q2/kT)` and `c3 = K_cl·D(T)/gdot`, then bisects 80 times on x = √ρ to find the root of `k1(x + Λ_b) − k2·x² − c3·x⁴ = 0`. The bisection is geometric (`mid = sqrt(lo·hi)`) because ρ spans many decades. Concavity guarantees one positive root.
   - Splits it: `rho_f = beta·rho`, `rho_m = f_m·rho_f`, `lam_f = 1/√rho_f`.
   - **Aging**, only if `C_sol > 0`: `t_a = rho_m·b·lam_f/gdot`, then `C/C0 = 1 + C_sol[1 − exp(−(t_a·ν·e^{−Q/kT})^α)]`.
   - **`tau_a = alpha·mu(T)·b·√(C/C0)·(√rho_f + Λ_b)`**, a direct formula.
2. **τ\*** — bisects 90 times in ln τ\* on `rho_m·b·glide_velocity(τ*) < gdot`. `glide_velocity` is the kink-pair law; it clips τ\*/τ_P below 1, applies `max(barrier, 0)` and, for Zr only (`B_free > 0`), blends the kink drag into the free-flight drag as the barrier vanishes. Strictly increasing in τ\*, so the root is unique.
3. **τ_p** — returns zero immediately if `tau_ppt == 0` (W and Zr). Otherwise `particle_rate` sums the detachment and climb rates and `particle_stress` bisects 80 times in ln τ_p. If the particles are already faster than the imposed rate at essentially zero stress, it returns 0.
4. **Sum**: `tau = tau_a + tau_star + tau_p`, returned in a dict with all the parts and the state.

Three independent one-dimensional root finds, each on a monotonic function, no iteration between them. That is the whole solve.

**Derived quantities:**
- **`srs`** calls `flow_stress` at γ̇·e^±0.01 for dτ/dlnγ̇, gets A by differentiating ln v_g analytically at the solved τ\*, then B = dτ/dlnγ̇ − 1/A and Λ = AB.
- **`global_m_V`** does what the experiments do: solves at three rates on a T grid and least-squares fits the slopes. V\* uses the axial convention, V\* = M_s·k_B·T/(dτ/dlnγ̇), so it compares directly with measured values.

## Stress control (creep and Ashby plots)

Here γ̇ is the unknown and the state depends on it, so the problem is implicit. The code doesn't iterate; it tabulates and inverts.

1. **`rate_table(T, P)`** — runs `flow_stress` on 241 rates from 10⁻¹⁴ to 10⁶ s⁻¹ at each T, giving τ(γ̇).
2. **`creep_rate(tau, T, P, diffusion)`** — for each T it takes a running maximum of that table (`np.maximum.accumulate`, plus a tiny ramp to keep it strictly increasing), then interpolates in ln τ to get log γ̇. The running maximum is what picks the lowest-rate branch where DSA makes τ(γ̇) non-monotonic. Stresses outside the table return ±inf. `diffusion=True` adds `diffusional_rate`; the creep scripts leave it off.

## Mechanism maps

`mechanism_map` builds a (T, τ/μ) grid, calls `rate_table` once, inverts it with `creep_rate` to get the rate everywhere, then calls `srs` on the same rate grid for diagnostics. For each grid point it looks up the nearest tabulated rate and labels the field: kink-pair if τ\*/τ > 0.5, climb or recovery if Λ > 1, otherwise obstacle-controlled; diffusional wherever that rate wins; −1 outside the table. The m < 0 region comes back separately as `dsa_band`, since it is unstable at fixed stress.

## The plot scripts

Each one sets `P = pm.W` (or `ZR`, `EU`) and a temperature or stress grid, then calls the physics and draws figures one at a time, each ending in `save(fig, name)`. The model scripts use `flow_stress`, `global_m_V`, `srs` and `state`; the creep scripts use `creep_rate`; the Ashby scripts use `mechanism_map`.

**Switching physics on and off** is done by copying the dict, which is how the build-up plots work:

```python
base = dict(P, C_sol=0.0, K_cl=0.0)      # kink pairs + forest only
no_climb = dict(P, D0_climb=0.0)          # particles, detachment route only
```

`C_sol = 0` skips the aging block, `tau_ppt = 0` short-circuits the particle solve, `K_cl = 0` drops the climb term from the density equation, and `Λ_b = 0` removes boundary storage. That's how W runs the same code as Eurofer without touching a line of `physics.py`.
