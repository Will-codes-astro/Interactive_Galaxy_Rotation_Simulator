# Physics reference

This file lists the equations the code uses, what each symbol means, its units,
and the function that implements it. Read it before changing any calculation.
Data columns and file formats are covered in [DATA_DICTIONARY.md](DATA_DICTIONARY.md).

Where it says "page JS", the code is inside `HTML_TEMPLATE` in
`generate_interactive_html.py`. To find a function there, search for `function <name>`.

## Units and constants

The code works in astronomical units: **kpc, km/s, M☉, Myr**. It only
converts to SI for MOND.

| Symbol | Value | Units | Where |
|---|---|---|---|
| G | 4.30091 × 10⁻⁶ | kpc (km/s)² M☉⁻¹ | `G_KPC_KMS2_PER_MSUN` (generator), `G` (page JS) |
| H₀ | 0.07 (= 70 km/s/Mpc) | km/s/kpc | `H0_KMS_PER_KPC`, `H0` |
| ρ_crit = 3H₀²/(8πG) | ≈ 136 | M☉/kpc³ | computed inline (`_nfw_velocity`, `nfwVelocityModel`) |
| 1 km/s | 1.0227121650537077 × 10⁻³ | kpc/Myr | `KM_S_TO_KPC_MYR` |
| 1 kpc/Myr | 977.8 | km/s | `KPC_PER_MYR_TO_KMS` |
| 1 kpc | 3.08567758 × 10¹⁹ | m | `M_PER_KPC` |
| a₀ (MOND) | 1.2 × 10⁻¹⁰ (default) | m/s² | `MOND_A0_DEFAULT` (Milgrom 1983; value as in Famaey & McGaugh 2012) |

The analysis scripts use the same G-free conversions. `simulate_baryonic_orbits.py`
rounds km/s→kpc/Myr to 0.00102271.

## Observations versus models

| Quantity | Kind | Meaning |
|---|---|---|
| V_obs ± σ | **Observation** (SPARC) | measured circular speed and its 1σ error |
| V_gas, V_disk, V_bulge | **Derived from observations** (SPARC mass models) | speed each component would produce alone; V_disk and V_bulge are for M/L = 1 |
| Υ_disk = 0.5, Υ_bulge = 0.7 | **Assumption** (3.6 μm averages; Lelli, McGaugh & Schombert 2016; McGaugh & Schombert 2014) | mass-to-light ratios. They are fixed, not fitted, except in the MOND best-fit search |
| V_bar | **Model** (from the above) | total visible-matter speed |
| V_halo, V_MOND, V_model | **Model prediction** | dark-matter or modified-gravity speed |

Every curve or column whose name has no `obs` in it is a model.

## Baryonic (visible-matter) speed

V_bar² = V_gas·|V_gas| + Υ_disk·V_disk² + Υ_bulge·V_bulge²

- **Keep the sign of V_gas.** A negative SPARC V_gas is real: at small radii the
  gas pulls outwards. Never use V_gas². The physics checks scan all Python for this.
- If the sum is negative, V_bar is set to 0 (`max(0, …)`).
- Implemented in: `fit_galaxy_halo` (generator); `baryonsFor` / `applyBaryonState`
  (page JS); `analyse_galaxy` (`analysis/baryonic_baseline.py`);
  `simulate_baryonic_orbits.py`; `_load_sparc_profile` (Blender).
- Between SPARC radii the page interpolates linearly (`baryonicVelocity`).
  Inside the first radius V rises linearly from 0 (solid-body assumption).
  Beyond the last radius V_bar ∝ r^(−1/2) (Keplerian: all visible mass enclosed).
  The offline orbit script uses PCHIP interpolation with the same two edge rules.

## NFW dark-matter halo (Navarro, Frenk & White 1996)

Parameters: M₂₀₀ (M☉), the mass inside r₂₀₀, and the scale radius r_s (kpc).

- r₂₀₀ = [3M₂₀₀ / (4π · 200 ρ_crit)]^(1/3), the radius where the mean density is 200 ρ_crit
- c = r₂₀₀ / r_s (concentration)
- f(x) = ln(1 + x) − x/(1 + x)
- M(<r) = M₂₀₀ f(r/r_s) / f(c) for r < r₂₀₀; **M(<r) = M₂₀₀ for r ≥ r₂₀₀**.
  This is a truncation assumption: the halo is treated as ending at r₂₀₀.
- V_halo² = G M(<r) / r
- V_model² = V_bar² + V_halo² (the contributions add in quadrature)

Implemented in: `_nfw_shape`, `_nfw_velocity` (generator) and `nfwVelocityModel` (page JS).

**Starting fit ("Apply Observed Fit")**: `fit_default_halo` does a grid search
over M₂₀₀ from 10¹¹ to 2×10¹² M☉ and r_s from 5 to 60 kpc. It minimises
χ² = Σ[(V_obs − V_model)/max(σ, 1 km/s)]², then snaps the result to the slider
steps (10¹⁰ M☉ and 0.5 kpc).

**Halo stability band**: in the sandbox, a warning shows when c is outside 3–25
(`STABLE_CONCENTRATION`). This range is a rule of thumb, not a published cut.

## Typical halo for a stellar mass (sandbox)

- **Stellar-to-halo mass relation** (Moster et al. 2013, z = 0):
  M*/M_h = 2N [(M_h/M₁)^(−β) + (M_h/M₁)^γ]⁻¹, with log M₁ = 11.59, N = 0.0351,
  β = 1.376, γ = 0.608. The code inverts it by bisection to get M₂₀₀ from M*.
  Page JS: `stellarMassOfHalo`, `haloMassFromStellar`.
- **Mass–concentration relation** (Dutton & Macciò 2014):
  log₁₀ c₂₀₀ = 0.905 − 0.101 log₁₀[M₂₀₀ / (10¹² h⁻¹ M☉)], with h = 0.7.
  Then r_s = r₂₀₀ / c. Page JS: `concentrationOf`, `baselineHalo`.

## MOND (Milgrom 1983)

μ(g/a₀) g = g_N, with the "simple-standard" interpolating function
μ(x) = x/√(1 + x²) (Famaey & McGaugh 2012). Solved exactly for g:

g = √{[g_N² + √(g_N⁴ + 4 g_N² a₀²)] / 2},  where g_N = V_bar²/r (SI), and V_MOND = √(g r)

- No halo is used in MOND mode (M₂₀₀ = 0).
- The limits are checked by the physics checks: g → g_N when g_N ≫ a₀, and
  V⁴ → G M a₀ when g_N ≪ a₀.
- **MOND best fit**: `findMondBestFit` grid-searches Υ_disk from 0.05 to 1.5
  and a₀ from 0.3 to 4 × 10⁻¹⁰ m/s², with Υ_bulge fixed, minimising χ².
- Implemented in `mondVelocity` and `extraVelocityModel` (page JS).

## Goodness of fit

χ² = Σᵢ [(V_obs,i − V_model,i) / σᵢ]²,  χ²_ν = χ² / (N − 2)

This is calculated only at measured SPARC radii and only for preset galaxies
(`updateChi`, page JS). `analysis/baryonic_baseline.py` uses dof = N instead
(there are no free parameters in a baryonic-only model). The two numbers are
therefore not directly comparable.

## Test-star orbits (website)

- **Massless test particles.** They feel the galaxy's gravity but not each other.
  The count is 800 for NGC 3198, and other galaxies scale it by
  (visible mass / template)^0.5, kept within 200–2400. Random numbers are
  seeded (`mulberry32(3198)`), so runs are repeatable.
- **Force:** spherical, a = −V_total²(r)/r² · **r** (r = 3D radius). Although the
  disk and bulge are not spherical, their combined pull is treated as a spherical
  approximation. The force is tabulated on 2048 radii out to 160 kpc
  (`fillForceTable`, `accelFactor`). Beyond 160 kpc, a ∝ r⁻².
- **Softening:** a cubic-spline kernel (`softeningFactor`), with
  ε = 0.5–1.0 kpc (`SOFTENING_RANGE_KPC`). It is exactly Newtonian beyond
  2.8 ε. The kernel is in the standard GADGET form; a source for it has not yet
  been added to the page (see Open questions).
- **Integrator:** velocity-Verlet (symplectic), 2 substeps per 1 Myr frame,
  over 1000 Myr (`integrateTrajectories`). A star counts as ejected beyond
  1.5 × R_max (`EJECTION_RADIUS_FACTOR`).
- **Initial conditions:** disk stars start on circular orbits in the *reference*
  gravity (the observed V_obs for a preset). Bulge stars are drawn from a Hernquist
  profile, with an isotropic velocity dispersion from the Jeans equation plus
  a mean rotation of 0.3 σ (`BULGE_ROTATION`). Changing the halo or MOND changes the gravity but not these
  starting speeds, so stars can drift or escape.
- **Star positions:** exponential disk with Gaussian heights (σ_z = 0.15 kpc); Hernquist
  bulge. The anatomy comes from the SPARC speeds (`anatomyFromRotation`): the
  disk scale length R_d = R_peak/2.16 and V²_peak = 0.3872 G M/R_d
  (Freeman 1970). These positions and the star colours are visual choices and
  do not affect the rotation curve.

## Offline orbit model (`analysis/simulate_baryonic_orbits.py`)

200 test particles in 2D (z = 0) move under **baryons only**, starting at
V_obs, for 300 Myr (301 frames). It uses `scipy.integrate.solve_ivp`. The script
shows that visible matter alone cannot hold stars moving at the observed speeds.
The output is `results/star_trajectories_baryonic.csv`, which is used by the
plotting scripts and the Blender model, **not** by the website.

## Open questions and known differences

These are documented here, not changed. Any change needs Will's decision.

1. χ²_ν uses dof = N − 2 on the page (also in MOND and visible-only modes) but
   dof = N in `baryonic_baseline.py`. The NFW starting fit uses a 1 km/s error
   floor and the page χ² does not.
2. The halo mass is truncated at r₂₀₀ (see NFW above).
3. In MOND mode, the site's CSV export puts the MOND boost in the `V_nfw_halo`
   column.
4. `anatomyFromRotation` estimates the gas mass from V_gas(last)²·r/G, which
   ignores the sign. This is visual only and does not affect the rotation curve.
5. The softening kernel has no source cited on the page yet.
6. The same baryonic and NFW formulas exist in five places: the generator
   Python, page JS, two analysis scripts and Blender. The website is the
   reference version. If you change one, change all of them, or record the
   difference here.
