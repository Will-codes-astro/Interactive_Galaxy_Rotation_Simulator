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
- Implemented in: `_chi_squared` (generator); `sparcBaryons` / `baryonsFor` / `applyBaryonState`
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

## Best realistic fits ("Apply Observed Fit", "MOND best fit")

Will's choice (15 Nov 2026): each best fit is the lowest χ² a model can reach
using only physically realistic values. It replaces the old starting fit (fixed
Υ, a 1 km/s error floor) and the page's MOND grid search (Υ_disk up to 1.5,
a₀ from 0.3 to 4 × 10⁻¹⁰). Every limit is `FIT_LIMIT_SIGMA` = 2 times the
published uncertainty either side of the typical value:

| Value | Limits | Source |
|---|---|---|
| Υ_disk | 0.5 × 10^(±0.2): 0.32 to 0.79 | 0.5 at 3.6 μm (McGaugh & Schombert 2014); 0.1 dex scatter, the prior width in Li et al. (2018, §2) |
| Υ_bulge (only galaxies with a bulge) | 0.7 × 10^(±0.2): 0.44 to 1.11 | as above |
| Distance D | D_SPARC ± 2 e_D | SPARC Table 1 (`table1.mrt`) |
| Inclination i | i_SPARC ± 2 e_i, at most 90° | SPARC Table 1 |
| NFW concentration c = r₂₀₀/r_s | c(M₂₀₀) × 10^(±0.22) | Dutton & Macciò (2014), eq. 8, intrinsic scatter 0.11 dex (§4.1) |
| NFW M₂₀₀ | 10¹⁰ to 10¹³ M☉, with r_s inside the slider's 5 to 60 kpc | search range only |
| MOND a₀ | 0.96 to 1.44 × 10⁻¹⁰ m/s² | 1.20 ± 0.24 (systematic) (McGaugh, Lelli & Schombert 2016). This is ±1 systematic uncertainty, not ±2σ: the value Will approved |

This follows Li et al. (2018), who fit each SPARC galaxy within "a realistic
region of the parameter space", except that they use Gaussian priors where this
site uses hard limits.

- **Distance and tilt.** A distance D scales every radius by f = D/D_SPARC and
  every squared visible-matter speed by f (mass ∝ D², radius ∝ D). An
  inclination i scales the measured speeds and their errors by
  s = sin i_SPARC / sin i (Li et al. 2018, eq. 5). Dividing the model by s
  instead gives the same χ², so the page leaves SPARC's points where they are
  and draws the model as V_model(r·f)/s (`frameOf`, `fitFrame`, `mondModelVelocity`).
- **Search** (`_minimise_in_box`): a grid of 5 values per parameter, then a
  pattern search from the 4 best grid points. On 15 Nov 2026 this matched a
  25⁴ brute-force grid for the MOND fits of UGC 128, NGC 3198 and NGC 2403 to
  within 0.1 in χ². D and i are partly degenerate, so their fitted values are
  less certain than the χ².
- **Snapping.** The NFW halo is snapped to the slider steps (10¹⁰ M☉, 0.5 kpc)
  and Υ, D and i are then refitted with the halo fixed, so the page shows the
  fit exactly. The physics checks confirm that the page's χ² matches the
  generator's (`chi2` in `galaxies.js`).
- **On the page.** Apply Observed Fit (and loading a preset) applies the NFW
  fit's halo, Υ, D and i; they stay while the halo sliders move, until Visible
  Matter Only returns to SPARC's standard values. The MOND best fit switch
  applies the MOND fit. Plain MOND and visible matter only always use SPARC's
  values.
- Implemented in: `fit_nfw`, `fit_mond`, `_realistic_settings`, `_chi_squared`
  (generator); `activeFit`, `applyBaryonState`, `scoreGalaxy` (page JS).

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
- **MOND best fit**: `fit_mond` (generator), within the realistic limits above.
- Implemented in `mondVelocity`, `mondModelVelocity` and `extraVelocityModel`
  (page JS), and `_mond_velocity` (generator, used only by the fit).

## Goodness of fit

χ² = Σᵢ [(V_obs,i − V_model,i) / σᵢ]²,  χ²_ν = χ² / (N − k)

k is the number of values tuned to the galaxy (`MODEL_TUNED_VALUES` in the
generator, read by `tunedCount` and `reducedChi` in page JS, and by
`tuned_count_text` for the k values quoted on the page): 2 for the NFW halo (M₂₀₀, r_s) or 1 for the MOND best fit (a₀), plus
the galaxy settings a best fit tunes (`tuned_settings`: Υ_disk, D, i, and
Υ_bulge where there is a bulge). Visible matter only and MOND at the default a₀
use SPARC's values, so k = 0. So the NFW best fit has k = 5 (6 with a bulge) and
the MOND best fit k = 4 (5). Will chose this rule on 14 Nov 2026 so that the two
models are compared fairly; before that the page used N − 2 in every mode. A
halo set by hand with the sliders still counts its 2, plus the fit's settings
while they are applied.

This is calculated only at measured SPARC radii and only for preset galaxies
(`updateChi` for the fit badge, `scoreGalaxy` for the scorecard).
`analysis/baryonic_baseline.py` also uses dof = N for its baryonic-only model.

**Fit verdict**: the words "Close fit" (χ²_ν < 2), "Rough fit" (2–10) and "Poor
fit" (> 10) come from `fitVerdict`. The cut-offs are a rule of thumb for this
site, not a published standard, and the page says so.

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

1. Resolved 15 Nov 2026: the NFW fit no longer uses a 1 km/s error floor or
   fixed Υ (see Best realistic fits). It still stays inside the scale-radius
   slider's range.
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
7. UGC 128 (added 15 Nov 2026, Will's choice from a shortlist). Resolved the
   same day: the old MOND best fit stopped at the top of its search range
   (Υ_disk = 1.5), with distance fixed. An unlimited search reaches χ²_ν ≈ 5.6,
   but only with Υ_disk = 2.7 and a₀ = 0.38 × 10⁻¹⁰, which are not realistic.
   The realistic MOND fit (above) varies distance and tilt too, and still
   scores χ²_ν ≈ 70, with Υ_disk, a₀ and the tilt all at their limits: within
   realistic values MOND cannot match this galaxy here. The NFW fit scores
   χ²_ν ≈ 3.5, with its distance at the lower limit (45 Mpc). Plain MOND
   scores χ²_ν ≈ 115. "Low surface brightness" is cited on the page to de Blok
   & McGaugh (1996), who use 60 Mpc for it.
8. Several best fits sit at a limit (for example the MOND a₀ for all four
   presets, and NGC 3198's NFW distance at 11.0 Mpc). The limits are hard
   edges; Li et al. (2018) use Gaussian priors, which would let a value move
   further at a cost in χ². A different FIT_LIMIT_SIGMA would change the
   scores, so it is Will's decision.
