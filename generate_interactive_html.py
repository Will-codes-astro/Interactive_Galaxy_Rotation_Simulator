"""Create a standalone Plotly dashboard from the local SPARC galaxy data."""

from __future__ import annotations

import csv
import itertools
import json
import math
from html import escape
from pathlib import Path
from typing import Any, NamedTuple

from analysis.sparc_data import SPARC_PATH, GalaxyInfo, galaxy_rows, load_galaxy_info


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_PATH = BASE_DIR / "index.html"
GALAXIES_JS_PATH = BASE_DIR / "galaxies.js"

G_KPC_KMS2_PER_MSUN = 4.30091e-6
H0_KMS_PER_KPC = 0.07
UPSILON_DISK = 0.50
UPSILON_BULGE = 0.70
M_PER_KPC = 3.08567758e19

# Best fits only use physically realistic values (Will's choice, 15 Nov 2026). Each limit is FIT_LIMIT_SIGMA
# times the published uncertainty either side of the typical value:
# - disk and bulge M/L: stellar-population models predict 0.5 and 0.7 at 3.6 um with about 0.1 dex scatter
#   (McGaugh & Schombert 2014; the prior width used by Li et al. 2018);
# - distance and inclination: SPARC Table 1's own errors (Lelli, McGaugh & Schombert 2016);
# - NFW concentration: the Dutton & Maccio (2014) relation, whose scatter is 0.11 dex;
# - MOND a0: the measured 1.20 +/- 0.02 (random, 1 sigma) x 10^-10 m/s^2 (McGaugh, Lelli & Schombert 2016, sec. III).
#   Their further +/-0.24 systematic is the 20% normalisation uncertainty in M/L, which the fits already vary,
#   so it is left out to avoid counting it twice; Li et al. (2018, sec. 5) use the same 1.20 +/- 0.02 (Will's choice, 15 Nov 2026).
FIT_LIMIT_SIGMA = 2.0
UPSILON_SCATTER_DEX = 0.1
HALO_CONCENTRATION_SCATTER_DEX = 0.11
MOND_A0_MEASURED_MS2 = 1.20e-10
MOND_A0_RANDOM_ERROR_MS2 = 0.02e-10
MOND_A0_RANGE_MS2 = (MOND_A0_MEASURED_MS2 - FIT_LIMIT_SIGMA * MOND_A0_RANDOM_ERROR_MS2,
                     MOND_A0_MEASURED_MS2 + FIT_LIMIT_SIGMA * MOND_A0_RANDOM_ERROR_MS2)
HALO_MASS_RANGE_MSUN = (1e10, 1e13)
HALO_SCALE_RANGE_KPC = (5.0, 60.0)  # the scale-radius slider's range
# Fit search: a grid of FIT_GRID_POINTS per value, then a pattern search from the best FIT_STARTS grid points,
# stopping when every step is below FIT_TOLERANCE of its range.
FIT_GRID_POINTS = 5
FIT_STARTS = 4
FIT_TOLERANCE = 1e-4
# A fitted value within LIMIT_TOLERANCE of its range from a limit is recorded as "at the limit" (the Scorecard marks it).
LIMIT_TOLERANCE = 1e-3
SETTING_NAMES = ("upsilon_disk", "upsilon_bulge", "distance_mpc", "inclination_deg")
# Values each model tunes to a galaxy, before the galaxy settings a best fit also tunes (its fit record's
# tuned_settings). k in chi2_r = chi2 / (N - k). The page's chi2_r and the Method text both come from this.
MODEL_TUNED_VALUES = {"visible": 0, "nfw": 2, "mond": 0, "mondBest": 1}


class GalaxyPreset(NamedTuple):
    """One entry of the galaxy dropdown."""

    sparc_id: str  # must match the Name column of data/sparc/table2.dat
    name: str
    summary: str


# The first preset is selected when the page opens. To offer another SPARC galaxy,
# add a preset here and run this script again.
GALAXY_PRESETS = (
    GalaxyPreset("NGC3198", "NGC 3198", "The classic test case: a calm spiral measured far beyond its stars"),
    GalaxyPreset("NGC2403", "NGC 2403", "A nearby, gas-rich spiral, measured in fine detail"),
    GalaxyPreset("NGC7331", "NGC 7331", "A big, bright spiral with a large central bulge"),
    GalaxyPreset("UGC00128", "UGC 128", "A faint, spread-out low-surface-brightness galaxy (de Blok & McGaugh 1996)"),
)
DEFAULT_PRESET = GALAXY_PRESETS[0]
ROTATION_COLUMNS = (
    "radius_kpc",
    "v_obs_kms",
    "v_obs_err_kms",
    "v_gas_kms",
    "v_disk_kms",
    "v_bulge_kms",
)
HALO_MASS_STEP_MSUN = 1e10  # resolution of the halo mass slider (0.01 x 10^12 Msun)
HALO_SCALE_STEP_KPC = 0.5  # resolution of the scale radius slider
# Concentrations c = r200 / rs of stable NFW halos; outside this range the sandbox shows its warning.
STABLE_CONCENTRATION = (3.0, 25.0)
# Test-star model. The default preset (NGC 3198) and the default sandbox galaxy are the BASE_STAR_COUNT
# template; other galaxies scale it by (visible mass / template mass) ** STAR_COUNT_EXPONENT, within the limits.
BASE_STAR_COUNT = 800
MIN_STAR_COUNT = 200
MAX_STAR_COUNT = 2400
STAR_COUNT_EXPONENT = 0.5
# A star is "ejected" once it is farther than EJECTION_RADIUS_FACTOR x R_max from the galactic centre.
EJECTION_RADIUS_FACTOR = 1.5
# Gravitational softening of the orbit integrator (kpc): the low end for pure disks; with a bulge it follows the
# bulge's scale radius within this range.
SOFTENING_RANGE_KPC = (0.5, 1.0)


def load_sparc_galaxy(sparc_id: str) -> dict[str, Any]:
    """One galaxy's rotation curve from the validated SPARC table (see analysis/sparc_data.py)."""
    rows = galaxy_rows(sparc_id, SPARC_PATH)
    rotation = {
        "radius_kpc": [row.rad_kpc for row in rows],
        "v_obs_kms": [abs(row.v_obs_kms) for row in rows],
        "v_obs_err_kms": [abs(row.e_v_obs_kms) for row in rows],
        # Vgas keeps its sign: the gas term is squared with its sign preserved.
        "v_gas_kms": [row.v_gas_kms for row in rows],
        "v_disk_kms": [abs(row.v_disk_kms) for row in rows],
        "v_bulge_kms": [abs(row.v_bulge_kms) for row in rows],
    }
    return {"distance_mpc": rows[-1].dist_mpc, "rotation": rotation}


def _nfw_shape(x: float) -> float:
    if x <= 0.0:
        return 0.0
    return math.log1p(x) - x / (1.0 + x)


def _nfw_velocity(radius_kpc: float, halo_mass_msun: float, scale_radius_kpc: float) -> float:
    critical_density = (
        3.0 * H0_KMS_PER_KPC**2 / (8.0 * math.pi * G_KPC_KMS2_PER_MSUN)
    )
    r200 = (
        3.0 * halo_mass_msun / (4.0 * math.pi * 200.0 * critical_density)
    ) ** (1.0 / 3.0)
    concentration = max(r200 / scale_radius_kpc, 1e-8)
    enclosed_mass = halo_mass_msun * _nfw_shape(radius_kpc / scale_radius_kpc)
    enclosed_mass /= max(_nfw_shape(concentration), 1e-12)
    if radius_kpc >= r200:
        enclosed_mass = halo_mass_msun
    return math.sqrt(G_KPC_KMS2_PER_MSUN * enclosed_mass / radius_kpc)


def _halo_concentration(m200_msun: float) -> float:
    """Typical NFW concentration for this halo mass at z = 0 (Dutton & Maccio 2014). Same as concentrationOf."""
    return 10 ** (0.905 - 0.101 * math.log10(m200_msun / (1e12 / 0.7)))


def _r200_kpc(m200_msun: float) -> float:
    critical_density = 3.0 * H0_KMS_PER_KPC**2 / (8.0 * math.pi * G_KPC_KMS2_PER_MSUN)
    return (3.0 * m200_msun / (4.0 * math.pi * 200.0 * critical_density)) ** (1.0 / 3.0)


def _mond_boost_rar(y: float) -> float:
    """g / g_N for y = g_N / a0: the radial acceleration relation, g = g_N / (1 - e^-sqrt(y)) (McGaugh, Lelli & Schombert 2016, eq. 4)."""
    return 1.0 / -math.expm1(-math.sqrt(y))


def _mond_boost_standard(y: float) -> float:
    """g / g_N for y = g_N / a0 with the "standard" mu(x) = x / sqrt(1 + x^2) (Famaey & McGaugh 2012). Only for the Method comparison."""
    return math.sqrt((1.0 + math.sqrt(1.0 + 4.0 / (y * y))) / 2.0)


def _mond_velocity(v_bar_kms: float, radius_kpc: float, a0_ms2: float, boost=_mond_boost_rar) -> float:
    """MOND speed, by default with the RAR interpolating function (Will's choice, 17 Nov 2026). Same as mondVelocity (page JS)."""
    r_m, v_ms = radius_kpc * M_PER_KPC, v_bar_kms * 1000.0
    if r_m <= 0.0 or v_ms <= 0.0:
        return 0.0
    g_n = v_ms * v_ms / r_m
    return math.sqrt(g_n * boost(g_n / a0_ms2) * r_m) / 1000.0


def mond_formula_table(galaxies: list[dict[str, Any]], galaxy_info: dict[str, GalaxyInfo]) -> str:
    """Method's table rows: plain MOND's reduced chi-squared (standard values, k = 0) with each interpolating function."""
    rows = []
    for galaxy in galaxies:
        info, n = galaxy_info[galaxy["id"].upper()], len(galaxy["rotation"]["radius_kpc"])
        standard = [lo for lo, _ in _standard_settings(info)]
        chi = [_chi_squared(galaxy["rotation"], info, standard, lambda v, r, b=boost: _mond_velocity(v, r, MOND_A0_MEASURED_MS2, b)) / n
               for boost in (_mond_boost_standard, _mond_boost_rar)]
        rows.append(f"<tr><td>{escape(galaxy['name'])}</td><td>{chi[0]:.3g}</td><td>{chi[1]:.3g}</td></tr>")
    return "".join(rows)


def _pattern_search(cost, x: list[float], fx: float, bounds, step: list[float]) -> tuple[float, list[float]]:
    """Move one value at a time while the cost falls, halving the steps when nothing helps."""
    while any(s > (hi - lo) * FIT_TOLERANCE for s, (lo, hi) in zip(step, bounds)):
        improved = False
        for d, (lo, hi) in enumerate(bounds):
            for sign in (1.0, -1.0):
                trial = x.copy()
                trial[d] = min(hi, max(lo, x[d] + sign * step[d]))
                f_trial = cost(trial)
                if f_trial < fx:
                    x, fx, improved = trial, f_trial, True
                    break
        if not improved:
            step = [s / 2.0 for s in step]
    return fx, x


def _minimise_in_box(cost, bounds) -> tuple[float, list[float]]:
    """Lowest cost with every value inside its limits: a coarse grid, then a pattern search from its best points."""
    axes = [
        [lo + (hi - lo) * k / (FIT_GRID_POINTS - 1) for k in range(FIT_GRID_POINTS)] if hi > lo else [lo]
        for lo, hi in bounds
    ]
    starts = sorted((cost(list(p)), list(p)) for p in itertools.product(*axes))[:FIT_STARTS]
    step = [(hi - lo) / (FIT_GRID_POINTS - 1) / 2.0 for lo, hi in bounds]
    return min(_pattern_search(cost, x, fx, bounds, step) for fx, x in starts)


def _realistic_settings(rotation: dict[str, list[float]], info: GalaxyInfo) -> list[tuple[float, float]]:
    """Limits for the galaxy's own settings: disk M/L, bulge M/L (only with a bulge), distance and inclination."""
    def ml_limits(typical: float) -> tuple[float, float]:
        spread = 10 ** (FIT_LIMIT_SIGMA * UPSILON_SCATTER_DEX)
        return typical / spread, typical * spread

    has_bulge = any(v > 0 for v in rotation["v_bulge_kms"])
    return [
        ml_limits(UPSILON_DISK),
        ml_limits(UPSILON_BULGE) if has_bulge else (UPSILON_BULGE, UPSILON_BULGE),
        (max(0.1, info.dist_mpc - FIT_LIMIT_SIGMA * info.e_dist_mpc), info.dist_mpc + FIT_LIMIT_SIGMA * info.e_dist_mpc),
        (max(5.0, info.inc_deg - FIT_LIMIT_SIGMA * info.e_inc_deg), min(90.0, info.inc_deg + FIT_LIMIT_SIGMA * info.e_inc_deg)),
    ]


def _chi_squared(rotation: dict[str, list[float]], info: GalaxyInfo, settings: list[float], model) -> float:
    """Chi-squared of a model at the given settings [disk M/L, bulge M/L, distance, inclination].

    A different distance D scales every radius by f = D / D_SPARC and every squared visible-matter speed by f
    (mass ~ D^2, radius ~ D). A different inclination i scales the measured speeds and their errors by
    s = sin(i_SPARC) / sin(i). Dividing the model by s instead leaves the same chi-squared, which is how the
    page draws the fitted model against SPARC's own points. model(v_bar, r) returns the model speed in km/s.
    """
    disk_ml, bulge_ml, distance_mpc, inclination_deg = settings
    f = distance_mpc / info.dist_mpc
    s = math.sin(math.radians(info.inc_deg)) / math.sin(math.radians(inclination_deg))
    total = 0.0
    for radius, observed, error, gas, disk, bulge in zip(*(rotation[name] for name in ROTATION_COLUMNS)):
        v_bar_sq = f * (gas * abs(gas) + disk_ml * disk**2 + bulge_ml * bulge**2)
        v_model = model(math.sqrt(max(0.0, v_bar_sq)), radius * f) / s
        total += ((observed - v_model) / error) ** 2
    return total


def _limits_hit(names, values, bounds) -> dict[str, str]:
    """Which fitted values stopped at a limit: {name: "low" or "high"}. Fixed values (no range) are skipped."""
    hits = {}
    for name, x, (lo, hi) in zip(names, values, bounds):
        if hi <= lo:
            continue
        tolerance = LIMIT_TOLERANCE * (hi - lo)
        if x <= lo + tolerance:
            hits[name] = "low"
        elif x >= hi - tolerance:
            hits[name] = "high"
    return hits


def _settings_record(settings: list[float], bounds) -> dict[str, Any]:
    disk_ml, bulge_ml, distance_mpc, inclination_deg = settings
    return {
        "upsilon_disk": round(disk_ml, 4),
        "upsilon_bulge": round(bulge_ml, 4),
        "distance_mpc": round(distance_mpc, 3),
        "inclination_deg": round(inclination_deg, 2),
        # Values tuned to the galaxy (a fixed bulge M/L is not tuned), counted in N - k.
        "tuned_settings": sum(hi > lo for lo, hi in bounds),
        "at_limit": _limits_hit(SETTING_NAMES, settings, bounds),
    }


def _standard_settings(info: GalaxyInfo) -> list[tuple[float, float]]:
    """SPARC's standard values held fixed (no range): disk and bulge M/L, distance and inclination."""
    return [(v, v) for v in (UPSILON_DISK, UPSILON_BULGE, info.dist_mpc, info.inc_deg)]


def fit_nfw(rotation: dict[str, list[float]], info: GalaxyInfo, settings_bounds=None) -> dict[str, Any]:
    """Best NFW halo: Fit halo to data and the scorecard's halo columns.

    Fitted: log M200, the concentration's offset from the typical value, and the galaxy's settings within
    settings_bounds (realistic limits by default; _standard_settings fixes them at SPARC's values, so only the halo
    is tuned). The halo is then snapped to the sliders' resolution and the settings refitted, so the page shows the
    result exactly.
    """
    settings_bounds = settings_bounds or _realistic_settings(rotation, info)
    log_c_limit = FIT_LIMIT_SIGMA * HALO_CONCENTRATION_SCATTER_DEX

    def halo_model(m200_msun: float, scale_radius_kpc: float):
        return lambda v_bar, radius: math.hypot(v_bar, _nfw_velocity(radius, m200_msun, scale_radius_kpc))

    def scale_radius(log_m200: float, log_c_offset: float) -> float:
        m200 = 10**log_m200
        return _r200_kpc(m200) / (_halo_concentration(m200) * 10**log_c_offset)

    def cost(p: list[float]) -> float:
        rs = scale_radius(p[0], p[1])
        if not HALO_SCALE_RANGE_KPC[0] <= rs <= HALO_SCALE_RANGE_KPC[1]:
            return math.inf
        return _chi_squared(rotation, info, p[2:], halo_model(10 ** p[0], rs))

    bounds = [tuple(math.log10(m) for m in HALO_MASS_RANGE_MSUN), (-log_c_limit, log_c_limit), *settings_bounds]
    _, best = _minimise_in_box(cost, bounds)
    m200 = round(10 ** best[0] / HALO_MASS_STEP_MSUN) * HALO_MASS_STEP_MSUN
    rs = round(scale_radius(best[0], best[1]) / HALO_SCALE_STEP_KPC) * HALO_SCALE_STEP_KPC
    chi, settings = _pattern_search(
        lambda p: _chi_squared(rotation, info, p, halo_model(m200, rs)),
        best[2:], _chi_squared(rotation, info, best[2:], halo_model(m200, rs)), settings_bounds, [(hi - lo) / 8.0 for lo, hi in settings_bounds],
    )
    record = _settings_record(settings, settings_bounds)
    # The halo is checked as shown, after snapping (which can move it slightly past a limit).
    halo_values = [math.log10(m200), math.log10(_r200_kpc(m200) / rs / _halo_concentration(m200)), rs]
    halo_bounds = [bounds[0], bounds[1], HALO_SCALE_RANGE_KPC]
    record["at_limit"] = {**_limits_hit(("m200_msun", "concentration", "scale_radius_kpc"), halo_values, halo_bounds), **record["at_limit"]}
    return {
        "m200_msun": m200,
        "scale_radius_kpc": rs,
        "concentration": round(_r200_kpc(m200) / rs, 3),
        **record,
        "chi2": round(chi, 3),
    }


def fit_mond(rotation: dict[str, list[float]], info: GalaxyInfo) -> dict[str, Any]:
    """Best MOND fit with realistic settings and a0 within its measured uncertainty: the Best realistic fits switch."""
    settings_bounds = _realistic_settings(rotation, info)

    def cost(p: list[float]) -> float:
        a0 = p[0] * 1e-10
        return _chi_squared(rotation, info, p[1:], lambda v_bar, radius: _mond_velocity(v_bar, radius, a0))

    bounds = [tuple(a / 1e-10 for a in MOND_A0_RANGE_MS2), *settings_bounds]
    chi, best = _minimise_in_box(cost, bounds)
    record = _settings_record(best[1:], settings_bounds)
    record["at_limit"] = {**_limits_hit(("a0_ms2",), best[:1], bounds[:1]), **record["at_limit"]}
    return {"a0_ms2": round(best[0], 4) * 1e-10, **record, "chi2": round(chi, 3)}


def build_galaxy_record(preset: GalaxyPreset, galaxy_info: dict[str, GalaxyInfo]) -> dict[str, Any]:
    """Collect everything the page needs to model one preset galaxy."""
    galaxy = load_sparc_galaxy(preset.sparc_id)
    info = galaxy_info[preset.sparc_id.upper()]
    if abs(info.dist_mpc - galaxy["distance_mpc"]) > 1e-6:
        raise ValueError(f"{preset.sparc_id}: Table 1 and table2 distances differ")
    return {
        "id": preset.sparc_id,
        "name": preset.name,
        "summary": preset.summary,
        "catalogue": {
            "distance_mpc": info.dist_mpc,
            "distance_err_mpc": info.e_dist_mpc,
            "inclination_deg": info.inc_deg,
            "inclination_err_deg": info.e_inc_deg,
            "hubble_type": int(info.hubble_type),
            "luminosity_36_1e9_lsun": info.lum_36_1e9_lsun,
            "hi_mass_1e9_msun": info.hi_mass_1e9_msun,
            "v_flat_kms": info.v_flat_kms,
        },
        "best_fit": {
            "nfw_standard": fit_nfw(galaxy["rotation"], info, _standard_settings(info)),
            "nfw": fit_nfw(galaxy["rotation"], info),
            "mond": fit_mond(galaxy["rotation"], info),
        },
        "rotation": galaxy["rotation"],
    }


GALAXIES_JS_HEADER = """/*
 * galaxies.js - SPARC rotation-curve presets for the Interactive Galaxy Rotation Simulator.
 *
 * Generated by generate_interactive_html.py from data/sparc/table2.dat and table1.mrt (SPARC mass
 * models and galaxy sample: Lelli, McGaugh & Schombert 2016, AJ 152, 157). Do not edit by hand: add or change presets in
 * GALAXY_PRESETS in the generator and run it again. The generator also inlines this code
 * into index.html, so the published page does not load this file separately.
 *
 * Radii are in kpc and velocities in km/s. v_gas_kms keeps SPARC's sign (the gas term can be
 * negative). catalogue holds SPARC Table 1's distance and inclination with their errors.
 * best_fit holds the generator's fits. nfw_standard tunes only the halo, with SPARC's standard
 * mass-to-light ratios, distance and inclination (the page's default). nfw and mond also tune those
 * settings, searched only over physically realistic values (and the halo concentration or a0): the
 * Best realistic fits switch. See fit_nfw and fit_mond. These are model results, not measurements.
 */
"""


def render_galaxies_js(galaxies: list[dict[str, Any]]) -> str:
    """Render the galaxies.js registry: a JavaScript lookup keyed by SPARC id."""

    def dump(value: Any) -> str:
        return json.dumps(value, allow_nan=False)

    lines = [GALAXIES_JS_HEADER.rstrip("\n"), "window.GALAXIES = {"]
    for galaxy in galaxies:
        lines += [
            f"  {dump(galaxy['id'])}: {{",
            f"    id: {dump(galaxy['id'])},",
            f"    name: {dump(galaxy['name'])},",
            f"    summary: {dump(galaxy['summary'])},",
            f"    catalogue: {dump(galaxy['catalogue'])},",
            "    best_fit: {",
            *(f"      {model}: {dump(fit)}," for model, fit in galaxy["best_fit"].items()),
            "    },",
            "    rotation: {",
            *(f"      {name}: {dump(galaxy['rotation'][name])}," for name in ROTATION_COLUMNS),
            "    },",
            "  },",
        ]
    lines.append("};")
    return "\n".join(lines) + "\n"


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Interactive Galaxy Rotation Simulator: explore SPARC rotation curves and star orbits for __GALAXY_NAMES__">
  <title>Interactive Galaxy Rotation Simulator</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&family=Source+Serif+4:opsz,wght@8..60,500;8..60,600&display=swap" rel="stylesheet">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/shaders/CopyShader.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/shaders/LuminosityHighPassShader.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/EffectComposer.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/RenderPass.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/ShaderPass.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/postprocessing/UnrealBloomPass.js"></script>
  <style>
    :root { color-scheme: dark; --bg:#050910; --panel:#0a111b; --line:#1b2a3a; --ink:#dbe7f4; --muted:#91a4b8; --cyan:#8ab4f8; --gold:#f2c879; }
    * { box-sizing:border-box; }
    sub, sup { line-height:0; } /* a subscript never makes its line taller than its neighbours */
    body { margin:0; background:var(--bg); color:var(--ink); font:14px/1.45 "Segoe UI",sans-serif; }
    header { display:flex; align-items:flex-end; justify-content:space-between; gap:20px; padding:22px 26px 16px; border-bottom:1px solid var(--line); background:#070d15; }
    .eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.12em; font-size:11px; }
    h1 { font-size:22px; line-height:1.15; margin:5px 0 0; font-weight:600; }
    .header-stats { display:flex; flex-wrap:nowrap; justify-content:flex-end; gap:8px; flex:0 0 auto; align-items:flex-start; }
    .stat { min-width:0; padding:7px 12px 8px; border:1px solid var(--line); border-radius:8px; background:linear-gradient(160deg,rgba(138,180,248,0.09),rgba(138,180,248,0.015)); white-space:nowrap; }
    .stat.ml { background:linear-gradient(160deg,rgba(242,200,121,0.12),rgba(242,200,121,0.02)); border-color:rgba(242,200,121,0.28); }
    .stat-label { color:var(--muted); font-size:9px; letter-spacing:.1em; text-transform:uppercase; }
    .stat-value { margin-top:1px; font-size:18px; font-weight:600; line-height:1.15; color:var(--cyan); font-variant-numeric:tabular-nums; }
    .stat-value small { font-size:12px; font-weight:400; color:var(--muted); margin-left:3px; }
    .stat.ml .stat-value { color:var(--gold); }
    .ml-col { display:flex; flex-direction:column; align-items:stretch; gap:6px; }
    .fit-toggle { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:5px 10px; border:1px solid rgba(138,180,248,0.28); border-radius:8px; background:rgba(138,180,248,0.05); font-size:11px; color:var(--muted); cursor:pointer; user-select:none; white-space:nowrap; }
    .fit-toggle[hidden] { display:none; }
    .fit-toggle input { position:absolute; opacity:0; width:0; height:0; }
    .fit-toggle .track { position:relative; width:30px; height:16px; border-radius:999px; background:#16222f; border:1px solid rgba(120,150,180,0.4); transition:background .2s,border-color .2s; flex:none; }
    .fit-toggle .track::after { content:""; position:absolute; top:2px; left:2px; width:10px; height:10px; border-radius:50%; background:#dbe7f4; transition:transform .2s; }
    .fit-toggle input:checked ~ .track { background:rgba(138,180,248,0.25); border-color:rgba(138,180,248,0.7); }
    .fit-toggle input:checked ~ .track::after { transform:translateX(14px); background:#8ab4f8; }
    .fit-toggle input:focus-visible ~ .track { outline:2px solid var(--cyan); outline-offset:2px; }
    .fit-toggle input:checked ~ .fit-label { color:var(--ink); }
    .ml-pair { display:flex; gap:16px; }
    .ml-source { display:block; margin-top:2px; color:var(--muted); font-size:10px; letter-spacing:.02em; }
    .ml-pair span { display:block; font-size:10px; font-weight:400; color:var(--muted); letter-spacing:.06em; }
    .ml-col { position:relative; }
    .ml-row { display:flex; gap:8px; }
    .ml-row .stat.ml { flex:1 1 auto; }
    .ml-col > .stat, .ml-col > .ml-row { flex:1 1 auto; } /* tiles fill the row's height, so they all match */
    @media (min-width:901px) { #tiltStat { min-width:166px; } } /* room for its longest source line, so the switch beside it stays put when it is flipped */
    .stat.ml .stat-label { display:flex; align-items:center; justify-content:space-between; gap:8px; }
    .ml-info { padding:0 6px; border:1px solid rgba(242,200,121,0.45); border-radius:999px; background:transparent; color:var(--gold); font:600 10px/1.5 "Segoe UI",sans-serif; letter-spacing:0; text-transform:none; cursor:pointer; }
    .ml-info:hover, .ml-info[aria-expanded="true"] { background:rgba(242,200,121,0.16); }
    .ml-info:focus-visible { outline:2px solid var(--cyan); outline-offset:2px; }
    .stat.ml.off { display:none; }
    [data-preset-only][hidden], .stability[hidden] { display:none !important; }
    .ml-explainer { position:absolute; top:calc(100% + 8px); right:0; z-index:40; width:340px; max-width:calc(100vw - 32px); padding:14px 16px; border:1px solid rgba(242,200,121,0.35); border-radius:8px; background:#0c1520; box-shadow:0 12px 32px rgba(0,0,0,0.55); font-size:12.5px; line-height:1.5; color:var(--ink); white-space:normal; }
    .ml-explainer { opacity:0; visibility:hidden; transform:translateY(-4px); pointer-events:none; transition:opacity .45s ease, transform .45s ease, visibility 0s linear .45s; }
    /* While the explainer is open, lift its section above the slider bar below it so the bar can't cover it. */
    .stats-row:has(.ml-explainer.open) { position:relative; z-index:160; }
    .stats-row { margin-bottom:14px; }
    .stats-row .header-stats { justify-content:flex-start; flex-wrap:wrap; }
    .stats-row .fit-toggle { align-self:center; padding:9px 14px; font-size:12px; }
    .sim-head .tutorial-launch { justify-self:end; }
    @media (max-width:900px) { .sim-head .tutorial-launch { justify-self:start; } }
    .ml-explainer.open { opacity:1; visibility:visible; transform:none; pointer-events:auto; transition:opacity .45s ease, transform .45s ease, visibility 0s; }
    @media (prefers-reduced-motion:reduce) { .ml-explainer, .ml-explainer.open { transform:none; transition:opacity .01s, visibility 0s; } }
    .ml-explainer h3 { margin:0 0 6px; font-size:13px; color:var(--gold); }
    .ml-explainer p { margin:0 0 8px; }
    .ml-explainer ul { margin:0 0 8px; padding-left:18px; }
    .ml-explainer .ml-source { font-size:11px; color:var(--muted); margin:0; }
    .ml-explainer .ml-close { position:absolute; top:8px; right:10px; border:0; background:none; color:var(--muted); font-size:16px; cursor:pointer; }
    main { padding:16px; max-width:1800px; margin:auto; }
    .controls { display:grid; grid-template-columns:minmax(230px,1fr) minmax(230px,1fr) minmax(280px,1.4fr); gap:14px; padding:13px 16px; border:1px solid var(--line); background:var(--panel); border-radius:5px; margin-bottom:14px; }
    .control label { display:flex; justify-content:space-between; gap:12px; color:var(--muted); font-size:12px; margin-bottom:7px; }
    .control output { color:var(--ink); font-variant-numeric:tabular-nums; }
    .control-help { color:var(--muted); font-size:11px; margin:6px 0 0; min-height:32px; }
    .fit-actions { display:flex; flex-wrap:wrap; gap:7px; margin-top:7px; }
    .fit-button { font-size:11px; }
    #visibleMatterOnly.tutorial-reveal { animation:tutorial-reveal .3s ease-out; }
    @keyframes tutorial-reveal { from { opacity:0; } to { opacity:1; } }
    .model-bar { display:flex; align-items:center; gap:16px; flex-wrap:wrap; padding:10px 16px; border:1px solid var(--line); background:var(--panel); border-radius:5px; margin-bottom:14px; }
    .model-bar .model-label { color:var(--muted); font-size:11px; letter-spacing:.08em; text-transform:uppercase; }
    .switch { display:inline-flex; align-items:center; gap:10px; cursor:pointer; user-select:none; }
    .switch input { position:absolute; opacity:0; width:0; height:0; }
    .switch .track { position:relative; width:44px; height:24px; border-radius:999px; background:#16222f; border:1px solid rgba(120,150,180,0.4); transition:background .2s,border-color .2s; flex:none; }
    .switch .track::after { content:""; position:absolute; top:2px; left:2px; width:18px; height:18px; border-radius:50%; background:#dbe7f4; transition:transform .2s; }
    .switch input:checked ~ .track { background:rgba(255,174,92,0.25); border-color:rgba(255,174,92,0.7); }
    .switch input:checked ~ .track::after { transform:translateX(20px); background:#ffae5c; }
    .switch input:focus-visible ~ .track { outline:2px solid var(--cyan); outline-offset:2px; }
    .switch .opt { color:var(--muted); font-size:13px; transition:color .2s; }
    .switch .opt.active { color:var(--ink); font-weight:600; }
    .control.disabled { opacity:.35; pointer-events:none; }
    .control.disabled input, .control.disabled button { cursor:not-allowed; }
    .mond-note { color:#ffc48a; font-size:12px; }
    .preset-picker { display:flex; align-items:center; gap:16px; flex-wrap:wrap; }
    .preset-picker[hidden], .sandbox[hidden], .stability-chip[hidden] { display:none; }
    /* Galaxy source: two big choice cards (real data or build your own), the first decision on the Simulator. Visual only. */
    .source-toggle { flex:1 1 100%; display:grid; grid-template-columns:1fr 1fr; gap:12px; }
    .source-toggle label { position:relative; display:flex; cursor:pointer; border-radius:12px; }
    .source-toggle input { position:absolute; inset:0; width:100%; height:100%; margin:0; opacity:0; cursor:pointer; }
    .source-card { flex:1; display:flex; align-items:flex-start; gap:14px; padding:14px 16px; border:1px solid var(--line-strong); border-radius:12px; background:rgba(255,255,255,.025); color:var(--muted); transition:background .2s,border-color .2s,box-shadow .2s,color .2s; }
    .source-icon { flex:none; width:40px; height:40px; padding:6px; border-radius:10px; background:rgba(138,180,248,.08); color:var(--muted); transition:color .2s,background .2s; }
    .source-text { display:flex; flex-direction:column; gap:2px; min-width:0; }
    .source-text small { font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; color:var(--muted); }
    .source-text b { font-size:17px; font-weight:600; color:var(--ink); }
    .source-text > span { font-size:12.5px; line-height:1.45; }
    .source-toggle label:hover input:not(:checked) + .source-card { border-color:rgba(138,180,248,.5); background:rgba(138,180,248,.05); color:var(--ink-2); }
    .source-toggle input:checked + .source-card { border-color:rgba(138,180,248,.75); background:linear-gradient(180deg,rgba(138,180,248,.16),rgba(138,180,248,.06)); color:var(--ink); box-shadow:0 0 0 1px rgba(138,180,248,.35), 0 0 22px rgba(120,170,255,.22); }
    .source-toggle input:checked + .source-card .source-icon { color:#cfe0ff; background:rgba(138,180,248,.22); }
    .source-toggle input:checked + .source-card .source-text small { color:var(--cyan); }
    .source-toggle input:focus-visible + .source-card { outline:2px solid var(--cyan); outline-offset:2px; }
    @media (max-width:640px) { .source-toggle { grid-template-columns:1fr; gap:8px; } .source-card { padding:11px 13px; gap:11px; } .source-icon { width:34px; height:34px; } .source-text b { font-size:15.5px; } }
    .sandbox { margin-bottom:14px; padding:13px 16px 14px; border:1px solid var(--line); background:var(--panel); border-radius:5px; }
    .galaxy-type { margin:0 0 12px; padding:10px 14px 11px; border:1px solid rgba(242,200,121,0.3); border-left:3px solid var(--gold); border-radius:6px; background:linear-gradient(160deg,rgba(242,200,121,0.08),rgba(242,200,121,0.01)); }
    .galaxy-type-kicker { margin:0; color:var(--muted); font-size:10px; letter-spacing:.1em; text-transform:uppercase; }
    .galaxy-type h3 { margin:2px 0 4px; color:var(--gold); font-size:16px; }
    .galaxy-type p { margin:0 0 5px; font-size:12.5px; line-height:1.45; }
    .galaxy-type .galaxy-type-example { color:var(--muted); font-size:12px; }
    .galaxy-type-tags { display:flex; flex-wrap:wrap; gap:6px; margin:6px 0; padding:0; list-style:none; }
    .galaxy-type-tags li { padding:2px 9px; border:1px solid var(--line); border-radius:999px; color:var(--cyan); font-size:11px; }
    .galaxy-type .galaxy-type-note { margin:0; color:var(--muted); font-size:10.5px; }
    .sandbox-head { display:flex; align-items:baseline; justify-content:space-between; gap:6px 16px; flex-wrap:wrap; margin-bottom:12px; }
    .sandbox-head p { margin:0; color:var(--muted); font-size:12px; }
    .sandbox-groups { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px 22px; }
    .sandbox-group { display:flex; flex-direction:column; gap:10px; min-width:0; }
    .sandbox-group h3 { margin:0; color:var(--cyan); font-size:11px; font-weight:600; letter-spacing:.1em; text-transform:uppercase; }
    .sandbox .control-help { min-height:0; }
    .sandbox-summary { margin:0; display:grid; grid-template-columns:1fr auto; gap:3px 12px; padding:9px 11px; border:1px solid rgba(120,150,180,0.2); border-radius:8px; background:rgba(120,150,180,0.05); font-size:12px; }
    .sandbox-stats-title { margin:10px 0 5px; }
    .sandbox-summary dt { color:var(--muted); white-space:nowrap; }
    .stat-q { position:relative; display:inline-grid; place-items:center; margin-left:5px; vertical-align:1px; width:14px; height:14px; flex:none; border:1px solid rgba(138,180,248,0.45); border-radius:3px; color:var(--cyan); font-size:9px; font-weight:600; line-height:1; cursor:help; }
    .stat-q::after { content:attr(data-tip); position:absolute; left:-4px; bottom:calc(100% + 6px); z-index:20; width:210px; padding:7px 9px; border:1px solid var(--line); border-radius:6px; background:#0d1825; color:var(--ink); font-size:11px; font-weight:400; line-height:1.4; text-align:left; white-space:normal; box-shadow:0 6px 18px rgba(0,0,0,0.5); opacity:0; visibility:hidden; pointer-events:none; transition:opacity .12s; }
    .stat-q:hover::after,.stat-q:focus-visible::after { opacity:1; visibility:visible; transition:opacity .3s ease .35s; }
    .sandbox-summary dd { margin:0; color:var(--ink); text-align:right; font-variant-numeric:tabular-nums; }
    .sandbox-footer { display:flex; align-items:center; justify-content:space-between; gap:10px 16px; flex-wrap:wrap; margin-top:12px; color:var(--muted); font-size:11px; }
    .sandbox-footer p { margin:0; }
    .stability { margin-top:14px; padding:12px 14px; border:1px solid rgba(95,208,224,0.3); border-radius:8px; background:rgba(95,208,224,0.05); }
    .stability[data-state=diffuse], .stability[data-state=dense] { border-color:rgba(255,138,128,0.55); background:rgba(255,138,128,0.07); }
    .stability[data-state=none] { border-color:rgba(120,150,180,0.25); background:rgba(120,150,180,0.05); }
    .stability-head { display:flex; align-items:center; justify-content:space-between; gap:8px 14px; flex-wrap:wrap; }
    .stability-title { color:#c3cfdb; font-size:11px; letter-spacing:.08em; text-transform:uppercase; }
    .stability-badge { padding:3px 11px; border-radius:999px; border:1px solid rgba(95,208,224,0.55); background:rgba(95,208,224,0.1); color:#8fe3ee; font-size:12px; font-weight:600; }
    .stability[data-state=diffuse] .stability-badge, .stability[data-state=dense] .stability-badge { border-color:rgba(255,138,128,0.7); background:rgba(255,138,128,0.14); color:#ffb4ac; }
    .stability[data-state=none] .stability-badge { border-color:rgba(120,150,180,0.35); background:rgba(120,150,180,0.08); color:var(--muted); font-weight:400; }
    .stability-reading { display:flex; align-items:baseline; flex-wrap:wrap; gap:4px 14px; margin:10px 0 2px; }
    .stability-c { font-size:20px; font-weight:600; color:var(--ink); font-variant-numeric:tabular-nums; }
    .stability-detail { color:var(--muted); font-size:12px; font-variant-numeric:tabular-nums; }
    .c-gauge { position:relative; height:8px; margin:12px 6px 0; border-radius:999px; background:rgba(255,138,128,0.28); }
    .c-band { position:absolute; top:0; bottom:0; border-radius:999px; background:rgba(95,208,224,0.55); }
    .c-marker { position:absolute; top:-4px; width:4px; height:16px; margin-left:-2px; border-radius:2px; background:#f4f7fb; box-shadow:0 0 0 2px rgba(5,9,16,0.8); transition:left .12s; }
    .stability[data-state=none] .c-marker { opacity:0; }
    .c-ticks { position:relative; height:16px; margin:5px 6px 0; color:var(--muted); font-size:10px; font-variant-numeric:tabular-nums; }
    .c-ticks span { position:absolute; transform:translateX(-50%); }
    .stability-message { margin:8px 0 0; color:#9fb0c2; font-size:12px; }
    .telemetry { position:absolute; top:8px; left:8px; z-index:100; padding:2px 6px; border:1px solid rgba(255,255,255,.07); border-radius:5px; background:rgba(8,14,24,.28); color:rgba(255,255,255,.5); font:9px/1.45 ui-monospace,Consolas,monospace; letter-spacing:.2px; pointer-events:none; }
    .tm-row[hidden], .halo-info[hidden] { display:none; }
    .tm-row { display:flex; justify-content:space-between; gap:9px; }
    .tm-row b { color:rgba(255,255,255,.78); font-weight:600; font-variant-numeric:tabular-nums; }
    /* χ² readout: fullscreen only (the curve panel shows it otherwise), and only for SPARC presets. */
    .tm-chi { display:none; }
    #orbitPlot:fullscreen .tm-chi:not([hidden]), #orbitPlot:-webkit-full-screen .tm-chi:not([hidden]), #orbitPlot.fullscreen-active .tm-chi:not([hidden]) { display:flex; }
    .tm-chi b.good { color:#93e3a0; } .tm-chi b.ok { color:#f2d675; } .tm-chi b.poor { color:#ff9d94; }
    .stability-chip { position:absolute; top:44px; left:12px; z-index:100; max-width:min(300px,calc(100% - 190px)); padding:5px 10px; border-radius:8px; border:1px solid rgba(255,138,128,0.6); background:rgba(40,10,10,0.6); -webkit-backdrop-filter:blur(10px); backdrop-filter:blur(10px); color:#ffb4ac; font-size:11px; font-weight:600; line-height:1.35; pointer-events:none; }
    @media (prefers-reduced-motion:reduce) { .c-marker { transition:none; } }
    @media (max-width:960px) { .sandbox-groups { grid-template-columns:1fr 1fr; } }
    @media (max-width:560px) { .sandbox-groups { grid-template-columns:1fr; } .preset-picker { gap:10px; } }
    .learn { margin-top:14px; display:grid; grid-template-columns:1fr 1fr; gap:14px; }
    .learn article { background:var(--panel); border:1px solid var(--line); border-radius:5px; padding:0 16px 6px; min-width:0; }
    .learn h3 { margin:0; padding:12px 0; font-size:14px; font-weight:600; }
    .learn article.nfw h3 { color:#b39ae8; } .learn article.mond h3 { color:#ffae5c; }
    .learn p { color:#9fb0c2; font-size:13px; margin:0 0 10px; }
    .learn p.tldr { color:#dbe7f4; }
    @media (max-width:900px) { .learn { grid-template-columns:1fr; } }
    input[type=range] { width:100%; accent-color:var(--cyan); }
    .panel-grid { display:grid; grid-template-columns:1.12fr 1fr; gap:14px; }
    .panel { min-width:0; background:var(--panel); border:1px solid var(--line); border-radius:5px; overflow:hidden; }
    .panel-heading { display:flex; align-items:baseline; justify-content:space-between; gap:12px; padding:13px 16px 0; }
    h2 { margin:0; font-size:14px; font-weight:600; }
    .panel-heading span { color:var(--muted); font-size:11px; }
    .key { margin:0 16px 14px; padding:0 14px; border:1px solid rgba(120,150,180,0.2); border-radius:8px; background:rgba(120,150,180,0.05); font-size:12px; color:#9fb0c2; }
    .key summary { cursor:pointer; padding:11px 0; color:#c3cfdb; font-size:11px; letter-spacing:.06em; text-transform:uppercase; }
    .key h3 { margin:14px 0 6px; color:#c3cfdb; font-size:12px; }
    .key p { margin:0 0 10px; }
    .key-content { padding-bottom:4px; }
    .key ul { list-style:none; margin:0; padding:0; display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr)); gap:5px 18px; }
    .key b { color:#e6edf5; font-weight:600; }
    .sw { display:inline-block; width:16px; height:3px; border-radius:2px; margin-right:7px; vertical-align:middle; }
    .sw.dot { width:8px; height:8px; border-radius:50%; margin:0 11px 0 4px; }
    .sw.dash { background:none; height:0; border-top:3px dotted; }
    .sw.bar { width:6px; height:10px; margin-right:2px; }
    .key .g { color:#93e3a0; } .key .o { color:#f2d675; } .key .r { color:#ff9d94; }
    #chiBadge { font-variant-numeric:tabular-nums; padding:3px 10px; border-radius:999px; border:1px solid rgba(120,150,180,0.3); background:rgba(120,150,180,0.08); color:#c3cfdb; font-size:12px; }
    #chiBadge.good { border-color:rgba(120,214,140,0.55); color:#93e3a0; } #chiBadge.ok { border-color:rgba(236,204,96,0.55); color:#f2d675; } #chiBadge.poor { border-color:rgba(255,138,128,0.6); color:#ff9d94; }
    .plot { width:100%; height:min(68vh,700px); min-height:420px; }
    /* Rotation-curve readout: a custom box, guide line and dots that glide between SPARC radii (visual only). */
    .curve-wrap { position:relative; }
    #curvePlot { touch-action:pan-y; -webkit-tap-highlight-color:transparent; }
    .curve-hover { position:absolute; inset:0; z-index:5; pointer-events:none; opacity:0; transition:opacity .18s ease; }
    .curve-hover.on { opacity:1; }
    .curve-hover .ch-line { position:absolute; left:0; top:0; width:0; border-left:1px dashed rgba(200,215,230,.45); will-change:transform; }
    .curve-hover .ch-dot { position:absolute; left:-5px; top:-5px; width:10px; height:10px; box-sizing:border-box; border-radius:50%; border:2px solid #070b14; box-shadow:0 0 0 1px rgba(255,255,255,.3),0 0 8px currentColor; background:currentColor; will-change:transform; }
    .curve-hover .ch-box { position:absolute; left:0; top:0; padding:4px 7px 5px; border:1px solid rgba(148,170,200,.25); border-radius:6px; background:rgba(10,15,26,.8); -webkit-backdrop-filter:blur(4px); backdrop-filter:blur(4px); color:#eef2f8; font:10.5px/1.4 Inter,"Segoe UI",system-ui,sans-serif; white-space:nowrap; box-shadow:0 4px 12px rgba(0,0,0,.35); will-change:transform; }
    .ch-box .ch-r { margin-bottom:1px; color:#8e9db3; font:500 10px "JetBrains Mono",ui-monospace,Consolas,monospace; }
    .ch-row { display:flex; align-items:center; gap:5px; }
    .ch-row i { flex:none; width:7px; height:2px; border-radius:1px; }
    .ch-row span { flex:1; margin-right:6px; color:#b8c4d6; }
    .ch-row b { font-weight:600; font-variant-numeric:tabular-nums; }
    .ch-row.ch-res { margin-top:2px; padding-top:2px; border-top:1px solid rgba(148,170,200,.16); }
    .curve-hint { display:none; margin:2px 16px 10px; color:var(--muted); font-size:12px; }
    .on-touch { display:none; }
    @media (hover:none) and (pointer:coarse) {
      .curve-hint { display:block; }
      .on-touch { display:inline; } .on-mouse { display:none; }
      .curve-hover .ch-box { padding:4px 7px 5px; font-size:11px; }
      .ch-box .ch-r { font-size:10.5px; }
      .stat-q { width:22px; height:22px; font-size:11px; }
      .stat-q:hover::after, .stat-q:focus-visible::after { opacity:0; visibility:hidden; }
      .stat-q.open::after { opacity:1; visibility:visible; width:min(210px,62vw); }
      .telemetry { padding:3px 7px; font-size:11px; }
      input[type=range] { height:34px; }
      input[type=range]::-webkit-slider-thumb { width:24px; height:24px; margin-top:-10px; }
      input[type=range]::-moz-range-thumb { width:20px; height:20px; }
      .fs-sliders input[type=range] { height:30px; }
      .fs-play { width:36px; height:36px; }
    }
    /* isolation: the readout, warning chip and star rings stack inside the 3D view, so they never show through the top menu. */
    #orbitPlot { position:relative; isolation:isolate; touch-action:none; overscroll-behavior:contain; user-select:none; -webkit-user-select:none; cursor:grab; }
    #orbitPlot.dragging { cursor:grabbing; }
    #orbitPlot canvas { display:block; width:100% !important; height:100% !important; touch-action:none; }
    #orbitPlot:fullscreen, #orbit-plot-container:fullscreen { width:100vw; height:100vh; max-height:none; border:0; border-radius:0; background:#020206; }
    #orbitPlot:-webkit-full-screen, #orbit-plot-container:-webkit-full-screen { width:100vw; height:100vh; max-height:none; border:0; border-radius:0; background:#020206; }
    .timeline-row { display:flex; align-items:center; gap:12px; }
    .timeline-row input[type=range] { flex:1; min-width:0; }
    .play-btn { display:inline-flex; align-items:center; justify-content:center; gap:6px; min-width:44px; min-height:44px; padding:0 14px; flex:none; font-size:13px; border-radius:22px; background:rgba(138,180,248,.12); border-color:var(--cyan); touch-action:manipulation; }
    .play-btn[aria-pressed=true] { background:rgba(138,180,248,.28); }
    @media(max-width:560px) { .play-btn .play-text { display:none; } .play-btn { width:44px; padding:0; } }
    .fullscreen-active { position:fixed !important; inset:0; z-index:9999; width:100vw !important; height:100dvh !important; max-height:none !important; min-height:0 !important; border:0 !important; border-radius:0 !important; background:#020206; }
    body.fullscreen-lock { overflow:hidden; }
    .fs-sliders { display:none; position:absolute; left:12px; bottom:12px; z-index:3; width:min(280px,calc(100% - 24px)); padding:8px 10px; gap:4px; flex-direction:column; border-radius:10px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.55); -webkit-backdrop-filter:blur(14px); backdrop-filter:blur(14px); color:var(--ink); font-size:11px; cursor:default; }
    #orbitPlot:fullscreen .fs-sliders, #orbitPlot:-webkit-full-screen .fs-sliders, #orbitPlot.fullscreen-active .fs-sliders { display:flex; }
    /* Outside fullscreen (laptops and tablets), the same box docks in the 3D view once the slider bar is off screen. */
    @media (min-width:641px) { body:not(.tutorial-active) #orbitPlot.slider-dock .fs-sliders { display:flex; } }
    .fs-sliders .fs-row { display:flex; flex-direction:column; gap:4px; }
    .fs-sliders .fs-row[hidden] { display:none; }
    .fs-sliders label { display:flex; justify-content:space-between; gap:8px; color:var(--muted); }
    .fs-sliders output { color:var(--ink); font-variant-numeric:tabular-nums; }
    .fs-sliders input[type=range] { width:100%; height:16px; margin:0; }
    .control[hidden] { display:none; }
    @media (min-width:961px) { .sim-stage .controls:has(#massControl[hidden]) { grid-template-columns:auto repeat(2,minmax(0,1fr)); } }
    .gravity-label { display:block; margin-bottom:9px; color:var(--ink-2); font-size:12.5px; font-weight:500; }
    /* Gravity model: two colour-coded cards. Purple = dark matter (the halo line and 3D halo), gold = MOND (its zone). */
    .model-cards { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr); gap:8px; margin:2px 0 0; }
    .model-card { position:relative; display:flex; cursor:pointer; border-radius:10px; }
    .model-card input { position:absolute; inset:0; width:100%; height:100%; margin:0; opacity:0; cursor:pointer; }
    .sim-stage .model-cards .model-card > span { flex:1 1 0; min-width:0; display:flex; flex-direction:column; gap:1px; padding:8px 11px 8px 26px; border:1px solid var(--line-strong); border-radius:10px; background:rgba(255,255,255,.025); color:var(--muted); transition:background .2s,border-color .2s,box-shadow .2s; }
    .model-card > span::before { content:''; position:absolute; left:10px; top:13px; width:9px; height:9px; border-radius:50%; border:1.5px solid currentColor; transition:background .2s; }
    .model-card b { font-size:14px; font-weight:600; color:var(--ink-2); }
    .model-card small { font-size:11.5px; line-height:1.3; }
    .model-card.is-dm { --mc:179,154,232; } .model-card.is-mond { --mc:242,200,121; }
    .model-card:hover input:not(:checked) + span { border-color:rgba(var(--mc),.55); }
    .model-card input:checked + span { border-color:rgba(var(--mc),.85); background:linear-gradient(180deg,rgba(var(--mc),.18),rgba(var(--mc),.06)); color:rgb(var(--mc)); box-shadow:0 0 0 1px rgba(var(--mc),.35), 0 0 18px rgba(var(--mc),.22); }
    .model-card input:checked + span::before { background:rgb(var(--mc)); }
    .model-card input:checked + span b { color:var(--ink); }
    .model-card input:focus-visible + span { outline:2px solid var(--cyan); outline-offset:2px; }
    .gravity-control .control-help { max-width:300px; }
    @media (max-width:960px) { .gravity-control { grid-column:1/-1; } }
    .fs-time { display:flex; align-items:center; justify-content:space-between; gap:8px; margin-bottom:4px; }
    .fs-play { flex:none; width:28px; height:28px; min-height:0; padding:0; border-radius:50%; display:inline-flex; align-items:center; justify-content:center; font-size:12px; line-height:1; }
    .fs-time output { margin-left:auto; }
    .track-input { width:100%; box-sizing:border-box; padding:6px 9px; border-radius:8px; border:1px solid rgba(255,255,255,.18); background:rgba(8,14,22,.6); color:var(--ink); font:inherit; font-size:13px; }
    .fs-time .track-input { width:58px; padding:3px 6px; font-size:11px; }
    .star-mark.track { color:#b8f58c; }
    .fs-time output { color:var(--muted); font-variant-numeric:tabular-nums; }
    /* Tutorial star marks: rings that follow chosen test stars, with a live label. Tutorial aids only. */
    .star-mark { position:absolute; left:0; top:0; z-index:4; width:0; height:0; color:#8ab4f8; pointer-events:none; }
    .star-mark.warm { color:#ffb27a; }
    .star-mark.edge::before { border-style:dashed; animation:none; }
    .star-mark.flip span { left:auto; right:17px; }
    .star-mark[hidden] { display:none; }
    .star-mark::before { content:''; position:absolute; left:-13px; top:-13px; width:22px; height:22px; border:2px solid currentColor; border-radius:50%; box-shadow:0 0 12px currentColor; animation:star-mark-pulse 1.6s ease-in-out infinite; }
    .star-mark span { position:absolute; left:17px; top:-10px; white-space:nowrap; padding:2px 7px; border-radius:6px; background:rgba(8,14,22,.8); color:#eef8ff; border:1px solid currentColor; font-size:11px; font-weight:600; font-variant-numeric:tabular-nums; }
    @keyframes star-mark-pulse { 0%,100% { transform:scale(1); opacity:1; } 50% { transform:scale(1.3); opacity:.6; } }
    @media (prefers-reduced-motion:reduce) { .star-mark::before { animation:none; } }
    /* Inside fullscreen the tutorial card stays bottom-right, clear of the corner controls. */
    #orbitPlot .tutorial-card, #orbitPlot .tutorial-card.on-left { left:auto; right:16px; bottom:16px; width:min(340px,calc(100% - 32px)); max-height:min(52vh,430px); }
    @media (max-width:700px) {
      #orbitPlot .tutorial-card, #orbitPlot .tutorial-card.on-left { right:10px; bottom:10px; width:calc(100% - 20px); max-height:38vh; padding:12px; }
      body.tutorial-active #orbitPlot .fs-sliders { bottom:calc(38vh + 20px); }
    }
    .view-panel { position:absolute; top:12px; right:12px; z-index:2; display:flex; flex-direction:column; gap:0; padding:10px; min-width:136px; border-radius:12px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.45); -webkit-backdrop-filter:blur(14px) saturate(140%); backdrop-filter:blur(14px) saturate(140%); box-shadow:0 8px 28px rgba(0,0,0,.35), inset 0 1px 0 rgba(255,255,255,.08); cursor:default; }
    .view-panel-title { color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; padding:0 2px 2px; }
    .view-panel button { width:100%; padding:6px 10px; font-size:12px; text-align:left; background:rgba(255,255,255,.06); border:1px solid rgba(255,255,255,.12); border-radius:8px; transition:background .15s,border-color .15s; }
    .view-panel button:hover, .view-panel button:focus-visible { background:rgba(138,180,248,.16); border-color:var(--cyan); outline:none; }
    .view-panel button:active { background:rgba(138,180,248,.28); border-color:var(--cyan); }
    .view-panel button[aria-pressed=true] { background:rgba(138,180,248,.22); border-color:var(--cyan); color:var(--ink); }
    .view-panel .view-row { display:grid; grid-template-columns:1fr 1fr; gap:6px; }
    .view-panel .view-row button { text-align:center; }
    .view-panel { min-width:188px; }
    .view-panel .view-primary button { font-weight:600; }
    .view-panel .view-label { color:var(--muted); font-size:9.5px; letter-spacing:.12em; text-transform:uppercase; padding:4px 2px 0; opacity:.8; }
    .view-panel [hidden] { display:none !important; }
    .view-panel .view-seg { display:grid; grid-auto-flow:column; grid-auto-columns:1fr; padding:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .view-seg button { padding:4px 2px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); }
    .view-panel .view-seg button:hover, .view-panel .view-seg button:focus-visible { color:var(--ink); background:rgba(138,180,248,.16); }
    .view-panel .view-switches { display:grid; grid-template-columns:1fr 1fr; gap:2px 8px; }
    .view-panel .view-switches button { display:flex; align-items:center; gap:7px; padding:4px 2px; font-size:11px; white-space:nowrap; border:0; background:transparent; color:var(--muted); }
    .view-panel .view-switches button:hover, .view-panel .view-switches button:focus-visible { background:transparent; color:var(--ink); }
    .view-panel .view-switches button[aria-checked=true] { background:transparent; color:var(--ink); }
    .view-panel .sw { position:relative; flex:none; width:24px; height:13px; border-radius:999px; background:rgba(255,255,255,.14); border:1px solid rgba(255,255,255,.18); transition:background .15s,border-color .15s; }
    .view-panel .sw::after { content:""; position:absolute; top:1px; left:1px; width:9px; height:9px; border-radius:50%; background:#cfd8e3; transition:transform .15s; }
    .view-panel [aria-checked=true] .sw { background:rgba(138,180,248,.45); border-color:var(--cyan); }
    .view-panel [aria-checked=true] .sw::after { transform:translateX(11px); background:#fff; }
    .view-panel .sw-eject[aria-checked=true] .sw { background:rgba(255,138,128,.45); border-color:#ff8a80; }
    .view-panel .view-switches .sw-mond { grid-column:1 / -1; } /* its longer label gets its own row, so the columns keep their width */
    .view-panel .sw-mond[aria-checked=true] .sw { background:rgba(242,200,121,.45); border-color:#f2c879; }
    /* With the 3D MOND zone on, the panels over the view go solid, so the gold behind never tints them or their text. */
    .mond-zone-on .view-panel, .mond-zone-on .fs-sliders, .mond-zone-on .telemetry, .mond-zone-on .stability-chip, .mond-zone-on .star-mark span { background:rgba(9,14,24,.94); -webkit-backdrop-filter:none; backdrop-filter:none; }
    /* Outside fullscreen the panel sits over a small plot, so it drops the section labels and tightens up. */
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel:not(.collapsed) { min-width:0; width:172px; padding:7px; top:8px; right:8px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-label { display:none; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel:not(.collapsed) > :not(.view-head) { margin-top:5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-row { gap:5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-row button { padding:4px 4px; font-size:11px; white-space:nowrap; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel :is(.view-seg,.speed-toggle,.mode-toggle) button { padding:3px 2px; font-size:10.5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-switches { grid-template-columns:auto auto; justify-content:space-between; gap:0 2px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-switches button { width:auto; gap:5px; padding:2px 0; font-size:10.5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .sw { width:20px; height:11px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .sw::after { width:7px; height:7px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel [aria-checked=true] .sw::after { transform:translateX(9px); }
    .view-panel .mode-toggle { display:grid; grid-template-columns:1fr 1fr; padding:2px; margin-top:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .mode-toggle button { padding:4px 8px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); }
    .view-panel .mode-toggle button[aria-checked=true] { background:rgba(138,180,248,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(138,180,248,.5); }
    .view-panel .speed-toggle { display:grid; grid-template-columns:repeat(4,1fr); padding:2px; margin-top:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .speed-toggle button { padding:4px 2px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); transition:background .15s,color .15s; }
    .view-panel .speed-toggle button:hover, .view-panel .speed-toggle button:focus-visible { color:var(--ink); outline:none; }
    .view-panel .speed-toggle button[aria-checked=true] { background:rgba(138,180,248,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(138,180,248,.5); }
    #orbitPlot.pan-mode { cursor:move; }
    #orbitPlot.pan-mode.dragging { cursor:grabbing; }
    .view-panel .view-head { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:0 2px 2px; background:transparent; border:0; border-radius:6px; color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; }
    .view-panel .view-head:hover, .view-panel .view-head:focus-visible { background:transparent; color:var(--ink); }
    .view-panel .chev { display:inline-block; font-size:12px; transition:transform .2s; }
    .view-panel.collapsed { min-width:0; padding:0; border-radius:999px; }
    .view-panel > :not(.view-head) { max-height:120px; margin-top:6px; opacity:1; visibility:visible; overflow:hidden; transition:max-height .2s ease, margin-top .2s ease, opacity .2s ease, visibility 0s; }
    .view-panel.collapsed > :not(.view-head) { width:0; min-width:0; max-height:0; margin-top:0; padding:0; border-width:0; opacity:0; visibility:hidden; pointer-events:none; transition:max-height .2s ease, margin-top .2s ease, opacity .2s ease, visibility 0s linear .2s; }
    .view-panel.collapsed .view-head { padding:7px 10px; border-radius:999px; }
    .view-panel.collapsed .chev { transform:rotate(90deg); }
    @media (hover:none) and (pointer:coarse) { .view-panel .mode-toggle { display:none; } .view-panel button { padding:9px 12px; font-size:13px; } .view-panel .speed-toggle button { padding:7px 4px; font-size:12px; } }
    button { background:#102130; color:var(--ink); border:1px solid #294255; border-radius:4px; padding:7px 12px; cursor:pointer; }
    button:hover { border-color:var(--cyan); }
    .halo-info { margin:0 16px 16px; padding:14px; border:1px solid rgba(179,154,232,.22); border-radius:8px; background:rgba(179,154,232,.05); }
    .halo-info h3 { margin:0; color:#d5c6f4; font-size:13px; }
    .halo-info-intro, .halo-empty { margin:5px 0 0; color:var(--muted); font-size:12px; }
    .halo-metrics { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:10px; margin:12px 0 0; }
    .halo-metrics[hidden] { display:none; }
    .halo-info:has(.halo-metrics[hidden]) .halo-info-intro { display:none; } /* zero halo mass: no values to hover over */
    .halo-metric { padding:10px; border:1px solid rgba(179,154,232,.14); border-radius:6px; background:rgba(5,9,16,.4); }
    .halo-metric button { display:block; width:100%; padding:0; border:0; background:none; text-align:left; color:#c3cfdb; font-size:11px; }
    .halo-metric button:focus-visible { outline:2px solid var(--cyan); outline-offset:3px; border-radius:2px; }
    .halo-metric .info-icon { color:#b39ae8; margin-left:3px; }
    .halo-metric .metric-value { display:block; margin-top:5px; color:#e6ddfa; font-size:18px; font-weight:600; font-variant-numeric:tabular-nums; }
    .halo-metric dd { display:grid; grid-template-rows:0fr; margin:0; color:var(--muted); font-size:11px; line-height:1.45; opacity:0; visibility:hidden; transition:grid-template-rows .45s ease, opacity .45s ease, visibility 0s linear .45s; }
    .halo-metric dd > span { min-height:0; overflow:hidden; }
    .halo-metric dd > span > span { display:block; padding-top:5px; }
    .halo-metric.expanded dd, .halo-metric.hovered dd { grid-template-rows:1fr; opacity:1; visibility:visible; transition-delay:0s; }
    @media (prefers-reduced-motion:reduce) { .halo-metric dd, .view-panel .sw, .view-panel .sw::after, .view-panel > :not(.view-head), .view-panel.collapsed > :not(.view-head) { transition:none; } }
    footer { margin:14px 2px 0; color:var(--muted); font-size:11px; display:flex; justify-content:space-between; gap:16px; }
    footer .export-btn { padding:0; border:0; background:none; color:var(--muted); font-size:11px; text-decoration:underline; text-underline-offset:3px; white-space:nowrap; }
    footer .export-btn:hover, footer .export-btn:focus-visible { color:var(--ink); }
    .tutorial-launch { margin-left:auto; background:transparent; color:var(--muted); font-size:12px; }
    .tutorial-launch:hover, .tutorial-launch:focus-visible { color:var(--ink); }
    .tutorial-welcome { width:min(460px,calc(100% - 32px)); max-height:calc(100% - 32px); padding:28px; border:1px solid rgba(138,180,248,.22); border-radius:18px; background:linear-gradient(145deg,#132739,#0b1522 65%); color:var(--ink); box-shadow:0 24px 80px rgba(0,0,0,.6); }
    .tutorial-welcome::backdrop { background:rgba(2,6,12,.8); }
    /* End of a tour: a short celebration (visual only). Three stars orbit a green tick, then the box settles. */
    .tutorial-done { position:fixed; inset:0; z-index:10002; display:grid; place-items:center; padding:16px; background:rgba(2,6,12,.6); animation:tutorial-done-in .45s ease both; }
    .tutorial-done[hidden] { display:none; }
    .tutorial-done.is-leaving { animation:tutorial-done-out .4s ease forwards; }
    .tutorial-done-box { width:min(400px,100%); padding:26px 24px 22px; text-align:center; border:1px solid rgba(165,227,194,.35); border-radius:18px; background:linear-gradient(145deg,#132739,#0b1522 65%); color:var(--ink); box-shadow:0 24px 80px rgba(0,0,0,.6),0 0 40px rgba(165,227,194,.12); animation:tutorial-done-pop .6s cubic-bezier(.2,1.4,.4,1) both; }
    .tutorial-done-box h2 { margin:14px 0 6px; font-size:22px; }
    .tutorial-done-box p { margin:0 0 18px; color:#b5c5d5; }
    .tutorial-done-orbit { position:relative; width:84px; height:84px; margin:0 auto; border:1px dashed rgba(165,227,194,.35); border-radius:50%; }
    .tutorial-done-orbit b { position:absolute; inset:22px; display:grid; place-items:center; border-radius:50%; background:#a5e3c2; color:#06131f; font-size:22px; box-shadow:0 0 22px rgba(165,227,194,.7); animation:tutorial-done-pop .6s .15s cubic-bezier(.2,1.4,.4,1) both; }
    .tutorial-done-orbit i { position:absolute; left:50%; top:50%; width:7px; height:7px; margin:-3.5px; border-radius:50%; background:#f4f7fb; box-shadow:0 0 8px #8ab4f8; animation:tutorial-done-spin 3s linear infinite; }
    .tutorial-done-orbit i:nth-child(2) { animation-delay:-1s; background:#f2c879; }
    .tutorial-done-orbit i:nth-child(3) { animation-delay:-2s; background:#a8ddff; }
    @keyframes tutorial-done-spin { from { transform:rotate(0) translateX(42px); } to { transform:rotate(360deg) translateX(42px); } }
    @keyframes tutorial-done-pop { from { opacity:0; transform:scale(.6); } }
    @keyframes tutorial-done-in { from { opacity:0; } }
    @keyframes tutorial-done-out { to { opacity:0; } }
    @media (prefers-reduced-motion:reduce) { .tutorial-done, .tutorial-done-box, .tutorial-done-orbit b, .tutorial-done-orbit i { animation:none; } .tutorial-done-orbit i { transform:rotate(calc(var(--k,0)*120deg)) translateX(42px); } .tutorial-done-orbit i:nth-child(2) { --k:1; } .tutorial-done-orbit i:nth-child(3) { --k:2; } }
    .tutorial-eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.1em; font-size:10px; font-weight:600; }
    .tutorial-welcome h2 { margin:8px 0; font-size:22px; }
    .tutorial-welcome p { color:#b5c5d5; margin:0 0 18px; }
    .tutorial-welcome-actions { display:flex; flex-wrap:wrap; gap:9px; }
    .tutorial-primary { background:rgba(138,180,248,.16); border-color:var(--cyan); color:#e4fbff; font-weight:600; }
    .tutorial-card { position:fixed; z-index:10000; right:16px; bottom:16px; display:flex; flex-direction:column; width:min(380px,calc(100% - 32px)); max-height:calc(100vh - 88px); overflow:hidden; padding:18px; border:1px solid rgba(138,180,248,.24); border-radius:16px; background:linear-gradient(155deg,#142838,#0b1623 70%); box-shadow:0 16px 48px rgba(0,0,0,.58); cursor:default; user-select:text; }
    .tutorial-card.on-left { left:16px; right:auto; }
    .tutorial-card[hidden] { display:none; }
    .tutorial-top { display:flex; flex:none; align-items:center; justify-content:space-between; gap:12px; }
    .tutorial-content { min-height:0; overflow-x:hidden; overflow-y:auto; overscroll-behavior:contain; }
    /* Hide the scrollbar while a step slides in or out, so it doesn't flash on Next. */
    .tutorial-card.is-switching .tutorial-content { overflow-y:hidden; }
    .tutorial-skip { border:0; padding:4px; background:transparent; color:var(--muted); font-size:11px; text-decoration:underline; text-underline-offset:3px; }
    .tutorial-skip:hover, .tutorial-skip:focus-visible { color:var(--ink); }
    .tutorial-progress-track { flex:none; height:3px; margin:10px 0 4px; background:rgba(138,180,248,.13); border-radius:4px; }
    .tutorial-progress-fill { height:100%; width:0; border-radius:inherit; background:var(--cyan); transition:width .25s ease; }
    .tutorial-card h2 { margin:10px 0 7px; font-size:18px; line-height:1.25; }
    .tutorial-card p { margin:0 0 10px; color:#c3d0de; font-size:13px; line-height:1.5; }
    .tutorial-card .tutorial-try { position:sticky; bottom:0; margin:2px 0 0; padding:10px 12px; border-left:2px solid var(--cyan); border-radius:0 7px 7px 0; background:#15283a; box-shadow:0 -10px 14px -6px #0f1e2d; color:#c7f5fa; }
    /* The task stays pinned in view at the bottom of the card while the explanation above it scrolls. */
    .tutorial-feedback { flex:none; display:flex; align-items:center; gap:10px; margin:9px 0 0; color:#a5e3c2; font-size:12px; animation:tutorial-pop .7s cubic-bezier(.2,.8,.2,1) both; }
    #tutorialFeedbackText { flex:1; min-width:0; }
    .tutorial-auto { flex:none; display:inline-flex; align-items:center; gap:6px; padding:3px 4px 3px 10px; border:1px solid rgba(165,227,194,.35); border-radius:999px; background:rgba(165,227,194,.08); color:#a5e3c2; font-size:11px; font-weight:600; cursor:pointer; transition:background .15s ease, opacity .2s ease; }
    .tutorial-auto[hidden] { display:none; }
    .tutorial-auto:hover, .tutorial-auto:focus-visible { background:rgba(165,227,194,.17); }
    .tutorial-auto.is-paused { opacity:.55; }
    .tutorial-clock { position:relative; width:20px; height:20px; }
    .tutorial-clock svg { display:block; width:100%; height:100%; transform:rotate(-90deg); }
    .tutorial-clock circle { fill:none; stroke-width:2; }
    .tutorial-clock .track { stroke:rgba(165,227,194,.2); }
    .tutorial-clock .ring { stroke:#a5e3c2; stroke-linecap:round; stroke-dasharray:43.98; }
    .tutorial-clock b { position:absolute; inset:0; display:grid; place-items:center; font-size:10px; font-variant-numeric:tabular-nums; }
    /* Moving between cards: the old text slides out the way you are going, the new text follows in, line by line. */
    .tutorial-card { transition:opacity .38s ease, transform .38s ease; }
    .tutorial-content, .tutorial-feedback { transition:opacity .42s ease, transform .42s ease; }
    .tutorial-try.part-in { animation:tutorial-pop .5s ease both; }
    .tutorial-card.leave :is(.tutorial-content,.tutorial-feedback) { opacity:0; transform:translateX(calc(var(--dir,1) * -12px)); }
    .tutorial-card.enter :is(#tutorialTitle,#tutorialText,#tutorialTry) { animation:tutorial-in .7s cubic-bezier(.2,.8,.2,1) both; }
    .tutorial-card.enter #tutorialText { animation-delay:.08s; }
    .tutorial-card.enter #tutorialTry { animation-delay:.16s; }
    /* When the card swaps sides of the page it fades out and rises back in, rather than jumping across. */
    .tutorial-card.leave-side { opacity:0; transform:translateY(10px) scale(.98); }
    .tutorial-card.enter-side { animation:tutorial-card-in .55s cubic-bezier(.2,.8,.2,1); }
    @keyframes tutorial-in { from { opacity:0; transform:translateX(calc(var(--dir,1) * 16px)); } }
    @keyframes tutorial-card-in { from { opacity:0; transform:translateY(10px) scale(.98); } }
    @keyframes tutorial-pop { from { opacity:0; transform:translateY(4px); } }
    /* The highlight settles onto each new control instead of appearing there. */
    .tutorial-spotlight.is-arriving { animation:tutorial-pulse 2.2s ease-in-out infinite, tutorial-arrive .7s cubic-bezier(.2,.8,.2,1); }
    @keyframes tutorial-arrive { from { opacity:0; transform:scale(1.06); } }
    .tutorial-feedback[hidden] { display:none; }
    .tutorial-actions { display:flex; flex:none; flex-wrap:wrap; justify-content:flex-end; gap:8px; margin-top:12px; }
    .tutorial-keys { align-self:center; margin-right:auto; color:var(--muted); font-size:11px; }
    @media (hover:none) and (pointer:coarse) { .tutorial-keys { display:none; } }
    .tutorial-actions button { min-height:36px; }
    .tutorial-actions button:disabled { opacity:.45; cursor:not-allowed; }
    .tutorial-spotlight { position:absolute; z-index:9998; box-sizing:border-box; pointer-events:none; border:3px solid #a8ddff; border-radius:9px; box-shadow:0 0 14px 4px rgba(120,210,255,.75), inset 0 0 12px rgba(120,210,255,.35); transition:opacity .45s ease; animation:tutorial-pulse 2.2s ease-in-out infinite; }
    /* Gentle pulse on the highlight ring (visual only). */
    @keyframes tutorial-pulse { 0%,100% { border-color:#a8ddff; box-shadow:0 0 14px 4px rgba(120,210,255,.75), inset 0 0 12px rgba(120,210,255,.35); } 50% { border-color:#ffffff; box-shadow:0 0 26px 9px rgba(120,210,255,.95), inset 0 0 16px rgba(120,210,255,.5); } }
    .tutorial-spotlight[hidden] { display:none; }
    /* One dim layer over the page with a hole for every highlighted control (main ring and the blue/yellow rings),
       so the things to use or watch stay at full brightness. No dim in full screen. */
    /* Rings and the dim layer are placed in page coordinates, so the browser scrolls them with the page and they stay
       locked to their controls (re-placing screen-fixed rings after each scroll made them jitter). A ring on something
       that does not scroll with the page (the top menu, the tour card, fullscreen) is fixed to the screen instead. */
    .tutorial-spotlight.is-fixed, .tutorial-also.is-fixed { position:fixed; }
    .tutorial-dim { position:absolute; top:0; left:0; z-index:9997; width:100%; height:100vh; pointer-events:none; fill:rgba(2,7,14,.42); opacity:0; transition:opacity .45s ease; }
    .tutorial-dim.is-on { opacity:1; }
    /* Every finished highlight (main ring, part rings, dim) holds green for 1.4 s then fades over 1.2 s. */
    .tutorial-dim.is-on.is-done { opacity:0; transition:opacity 1.2s ease-in-out 1.4s; }
    .tutorial-also { position:absolute; z-index:9998; box-sizing:border-box; pointer-events:none; border:2px solid #a8ddff; border-radius:9px; box-shadow:0 0 12px 3px rgba(120,210,255,.6); animation:tutorial-also-pulse 2.2s ease-in-out infinite; transition:opacity .3s ease; }
    .tutorial-also[hidden] { display:none; }
    /* Controls the current tour card does not use are locked (inert) and dimmed. */
    .tutorial-locked { opacity:.38; filter:saturate(.4); transition:opacity .25s ease; }
    /* Blue boxes = controls to use; yellow dashed boxes = things to watch. */
    /* Phones: live copies of watched values that are off screen, boxed like the yellow watch rings. */
    .tutorial-readout { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
    .tutorial-readout span { padding:4px 9px; border:1.5px dashed #ffcf7a; border-radius:8px; color:var(--ink); font-size:12px; font-variant-numeric:tabular-nums; }
    .tutorial-readout b { margin-right:5px; color:#ffcf7a; font-weight:600; }
    .tutorial-also.is-watch { border:2px dashed #ffcf7a; box-shadow:0 0 12px 3px rgba(255,207,122,.45); }
    /* A small numbered dot on the corner of each ring when a task has several steps (matches the card's list). */
    .tutorial-tag { position:absolute; left:-11px; top:-11px; display:grid; place-items:center; width:20px; height:20px; border-radius:50%; background:#a8ddff; color:#06131f; font-size:11px; font-weight:700; line-height:1; box-shadow:0 0 0 2px #06131f, 0 0 10px rgba(120,210,255,.8); pointer-events:none; }
    .tutorial-tag:empty { display:none; }
    .tutorial-try ol { margin:5px 0 0; padding-left:20px; }
    .tutorial-slider { margin-top:10px; padding:9px 12px 6px; border:1px solid rgba(138,180,248,.24); border-radius:7px; background:rgba(120,150,180,.06); }
    .tutorial-slider label { display:flex; justify-content:space-between; gap:10px; font-size:12px; color:var(--ink); }
    .tutorial-slider output { font-variant-numeric:tabular-nums; color:var(--cyan); }
    .tutorial-slider input { width:100%; margin:6px 0 0; }
    .tutorial-sliders { flex:none; }
    .tutorial-sliders[hidden] { display:none; }
    .tutorial-spotlight.in-card, .tutorial-also.in-card { z-index:10001; }
    .tutorial-try li { margin:2px 0; }
    .tutorial-try .tutorial-more-parts { display:block; margin-top:3px; color:var(--muted); font-size:11px; font-weight:400; }
    .tutorial-try .tutorial-first { display:block; margin-bottom:4px; color:#ffcf7a; }
    /* The three answers share one row; they only wrap on a very narrow card. */
    .tutorial-choices { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
    .tutorial-choices button { flex:1 1 auto; min-height:34px; padding:6px 8px; border-radius:999px; font-size:12px; white-space:nowrap; }
    .tutorial-try .tutorial-prediction { display:block; margin-bottom:4px; color:var(--muted); font-size:11px; }
    /* The current task's ring is the brightest; the other rings sit a little dimmer. */
    @keyframes tutorial-also-pulse { 0%,100% { opacity:.5; } 50% { opacity:.68; } }
    /* A finished part's ring turns green and fades away instead of vanishing (visual only). */
    .tutorial-also.is-leaving { animation:tutorial-ring-done 2.6s ease-in-out forwards; }
    /* The next part's ring hands over to the main ring with a short cross-fade instead of vanishing. */
    .tutorial-also.is-handover { animation:tutorial-handover .7s ease forwards; }
    @keyframes tutorial-handover { to { opacity:0; } }
    .tutorial-also.is-leaving .tutorial-tag { background:#a5e3c2; }
    @keyframes tutorial-ring-done { 0% { opacity:1; border-color:#a5e3c2; box-shadow:0 0 14px 4px rgba(165,227,194,.7); } 54% { opacity:1; border-color:#a5e3c2; box-shadow:0 0 18px 6px rgba(165,227,194,.6); } 100% { opacity:0; border-color:#a5e3c2; box-shadow:0 0 6px 1px rgba(165,227,194,0); } }
    .tutorial-pointer { position:absolute; left:0; top:0; width:25px; height:30px; filter:drop-shadow(0 2px 3px #000); animation:tutorial-gesture 2.8s ease-in-out infinite; will-change:transform; transition:opacity .5s ease; }
    .tutorial-spotlight.is-interacting .tutorial-pointer { opacity:.15; }
    .tutorial-spotlight.is-reached .tutorial-pointer { opacity:0; animation-play-state:paused; }
    /* When the task is done the highlight turns green with a soft glow and one small ripple, holds, then fades (visual only). */
    .tutorial-spotlight.is-complete { opacity:0; border-color:#a5e3c2; animation:tutorial-done-glow .7s ease-out; box-shadow:0 0 12px 3px rgba(165,227,194,.5); transition:opacity 1.2s ease-in-out 1.4s, border-color .5s ease, box-shadow .5s ease; }
    @keyframes tutorial-done-glow { 0% { box-shadow:0 0 18px 5px rgba(165,227,194,.7); } }
    /* A wrong prediction: the same fade, in orange instead of green, with no ripple. */
    .tutorial-feedback.is-missed { color:#ffb36b; }
    .tutorial-spotlight.is-complete.is-missed { border-color:#ffb36b; box-shadow:0 0 12px 3px rgba(255,179,107,.45); animation:none; }
    .tutorial-spotlight.is-complete.is-missed::after { content:none; }
    /* The ripple grows by a fixed distance, so it looks the same on a small button and on the whole graph. */
    .tutorial-spotlight.is-complete::after { content:''; position:absolute; inset:-3px; border:1.5px solid #a5e3c2; border-radius:11px; animation:tutorial-ripple .7s ease-out both; }
    .tutorial-spotlight.is-revisit, .tutorial-spotlight.is-revisit::after { animation:none; transition:none; opacity:0; }
    @keyframes tutorial-ripple { from { opacity:.5; } to { opacity:0; inset:-11px; border-radius:16px; } }
    /* The card's border tints green for a moment when a task is done. */
    .tutorial-card.is-celebrating { animation:tutorial-card-glow 1.2s ease-out; }
    @keyframes tutorial-card-glow { 0% { border-color:rgba(165,227,194,.45); } }
    #tutorialFeedbackText b { display:inline-block; font-weight:600; animation:tutorial-tick-pop .35s ease-out both; }
    @keyframes tutorial-tick-pop { from { opacity:0; transform:scale(.85); } }
    .tutorial-spotlight.is-complete .tutorial-pointer { opacity:0; animation-play-state:paused; }
    .tutorial-spotlight.is-click .tutorial-pointer { animation-name:tutorial-click; }
    .tutorial-trim { display:none !important; }
    /* look: a card that only asks you to read something gets the box but no pointer, since nothing needs hovering. */
    .tutorial-spotlight.is-look .tutorial-pointer { display:none; }
    /* groups: the pointer rests on each group's heading in turn (no drag), then pops back to the first. */
    @keyframes tutorial-gesture { 0%,15% { transform:translate(var(--from-x),var(--from-y)); } 65%,85% { transform:translate(var(--to-x),var(--to-y)); } 100% { transform:translate(var(--from-x),var(--from-y)); } }
    @keyframes tutorial-click { 0%,35%,100% { transform:translate(var(--to-x),var(--to-y)) scale(1); } 50%,65% { transform:translate(var(--to-x),var(--to-y)) scale(.8); } }
    @media (prefers-reduced-motion:reduce) { .tutorial-pointer, .tutorial-spotlight.is-click .tutorial-pointer, #visibleMatterOnly.tutorial-reveal { animation:none; } .tutorial-pointer { transform:translate(var(--to-x),var(--to-y)); } .tutorial-progress-fill, .tutorial-spotlight, .tutorial-card, .tutorial-content, .tutorial-feedback { transition:none; } .tutorial-spotlight, .tutorial-spotlight.is-arriving, .tutorial-spotlight.is-complete, .tutorial-card.is-celebrating, #tutorialFeedbackText b, .tutorial-also, .tutorial-feedback,.tutorial-card.enter :is(#tutorialTitle,#tutorialText,#tutorialTry), .tutorial-card.enter-side { animation:none; } .tutorial-spotlight.is-complete::after { display:none; } }
    #orbitPlot:fullscreen { touch-action:auto; }
    #orbitPlot:-webkit-full-screen { touch-action:auto; }
    #orbitPlot.fullscreen-active { touch-action:auto; }
    body.tutorial-active main { padding-bottom:calc(65vh + 32px); }
    @media(max-width:600px) {
      .tutorial-launch { margin-left:0; }
      .tutorial-card { right:10px; bottom:10px; width:calc(100% - 20px); max-height:45vh; padding:12px; }
      .tutorial-card.on-left { left:auto; right:10px; }
      .tutorial-card h2 { margin:5px 0; }
      .tutorial-card p { font-size:12px; margin-bottom:5px; }
      .tutorial-actions { margin-top:7px; }
      body.tutorial-active main { padding-bottom:calc(45vh + 20px); }
    }
    @media(max-width:960px) { .panel-grid { grid-template-columns:1fr; } .plot { height:62vh; } .controls { grid-template-columns:1fr 1fr; } .control:last-child { grid-column:1/-1; } }
    @media(max-width:560px) { header { align-items:flex-start; flex-direction:column; padding:17px; } .header-stats { justify-content:flex-start; flex-wrap:wrap; } main { padding:10px; } .controls { grid-template-columns:1fr; } .control:last-child { grid-column:auto; } .plot { min-height:360px; height:56vh; } footer { flex-direction:column; } }
    /* ===== Design system: tokens and component layer (visual only; no physics). ===== */
    :root { --bg:#060a12; --panel:#0c1322; --panel-2:#101929; --line:rgba(148,170,200,.14); --line-strong:rgba(148,170,200,.28); --ink:#e7edf6; --ink-2:#b8c4d6; --muted:#8392a8; --accent:#8ab4f8; --cyan:var(--accent); --gold:#e9c88a;
      --r:10px; --r-lg:14px; --f-ui:Inter,"Segoe UI",system-ui,sans-serif; --f-serif:"Source Serif 4",Georgia,serif; --f-mono:"JetBrains Mono",ui-monospace,Consolas,monospace; --shadow:inset 0 1px 0 rgba(255,255,255,.03),0 12px 32px rgba(0,0,0,.28); }
    @media (prefers-reduced-motion:no-preference) { html { scroll-behavior:smooth; } }
    body { font:14px/1.55 var(--f-ui); -webkit-font-smoothing:antialiased; background:radial-gradient(1100px 520px at 80% -8%,rgba(70,100,170,.16),transparent 62%),radial-gradient(800px 460px at -8% 6%,rgba(120,90,180,.08),transparent 60%),var(--bg); background-repeat:no-repeat; min-height:100vh; }
    .section-head[id], #curvePanel, #theoryGuide, #references { scroll-margin-top:72px; }
    .topbar { position:sticky; top:0; z-index:60; display:flex; align-items:center; justify-content:space-between; gap:16px; height:56px; padding:0 28px; border-bottom:1px solid var(--line); background:rgba(6,10,18,.74); -webkit-backdrop-filter:blur(14px) saturate(140%); backdrop-filter:blur(14px) saturate(140%); }
    .brand { display:inline-flex; align-items:center; gap:10px; color:var(--ink); font-weight:600; font-size:14px; letter-spacing:.01em; text-decoration:none; white-space:nowrap; }
    .brand svg { width:24px; height:24px; color:var(--accent); flex:none; }
    .topnav { display:flex; gap:2px; overflow-x:auto; scrollbar-width:none; }
    .topnav a { padding:6px 11px; border-radius:7px; color:var(--muted); font-size:13px; font-weight:500; text-decoration:none; white-space:nowrap; transition:color .15s,background .15s; }
    .topnav .nav-short { display:none; }
    @media (max-width:900px) { .topnav .nav-long { display:none; } .topnav .nav-short { display:inline; } } /* tablets: room for the menu and the gear */
    /* Phones: short labels, and the menu wraps onto a second row rather than running off the screen. */
    @media (max-width:640px) {
      .topbar { height:auto; min-height:56px; }
      .topnav { flex-wrap:wrap; justify-content:flex-end; row-gap:2px; overflow:visible; }
      .topnav .nav-long { display:none; } .topnav .nav-short { display:inline; }
    }
    @media (max-width:420px) { .topbar .brand { display:none; } .topbar { justify-content:center; } .topnav { margin-right:auto; } }
    .topnav a:hover, .topnav a:focus-visible { color:var(--ink); background:rgba(255,255,255,.05); outline:none; }
    /* Settings: a small gear at the far right of the top bar opens a panel of display choices. */
    .topnav { margin-left:auto; }
    .settings { position:relative; flex:none; margin-right:-10px; }
    .settings-btn { display:grid; place-items:center; width:32px; height:32px; padding:0; border:0; border-radius:8px; background:none; color:var(--muted); cursor:pointer; transition:color .15s,background .15s; }
    .settings-btn svg { width:18px; height:18px; transition:transform .4s ease; }
    .settings-btn:hover, .settings-btn:focus-visible, .settings-btn[aria-expanded=true] { color:var(--ink); background:rgba(255,255,255,.06); outline:none; }
    .settings-btn:focus-visible { box-shadow:0 0 0 2px var(--accent); }
    .settings-btn[aria-expanded=true] svg { transform:rotate(60deg); }
    .settings-panel { position:absolute; top:calc(100% + 10px); right:0; z-index:70; width:min(300px,calc(100vw - 24px)); padding:12px; border:1px solid var(--line-strong); border-radius:var(--r); background:#0d1625; box-shadow:0 14px 36px rgba(0,0,0,.5); }
    .settings-panel h2 { margin:0 2px 8px; color:var(--muted); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    .settings-row { display:flex; align-items:flex-start; gap:11px; padding:9px 10px; border-radius:8px; cursor:pointer; user-select:none; transition:background .15s; }
    .settings-row:hover { background:rgba(255,255,255,.045); }
    .settings-row input { appearance:none; -webkit-appearance:none; flex:none; display:grid; place-items:center; width:18px; height:18px; margin:1px 0 0; border:1.5px solid var(--line-strong); border-radius:5px; background:rgba(255,255,255,.03); cursor:pointer; transition:background .15s,border-color .15s; }
    .settings-row input::after { content:""; width:9px; height:5px; margin-top:-2px; border:solid #06101e; border-width:0 0 2px 2px; transform:rotate(-45deg) scale(0); transition:transform .15s ease; }
    .settings-row input:checked { background:var(--accent); border-color:var(--accent); }
    .settings-row input:checked::after { transform:rotate(-45deg) scale(1); }
    .settings-row input:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
    .settings-row input:disabled, .settings-row:has(input:disabled) { cursor:default; opacity:.55; }
    .settings-row b { display:block; color:var(--ink); font-size:13px; font-weight:500; }
    .settings-row small { display:block; margin-top:2px; color:var(--muted); font-size:11.5px; line-height:1.45; }
    .settings-note { display:flex; align-items:center; justify-content:space-between; gap:8px; margin:6px 10px 0; color:var(--muted); font-size:11px; }
    .settings-note button { padding:0; border:0; background:none; color:var(--accent); font:inherit; cursor:pointer; white-space:nowrap; }
    .settings-note button:hover, .settings-note button:focus-visible { text-decoration:underline; outline:none; }
    @media (max-width:640px) { .settings { margin-right:-4px; } }
    @media (prefers-reduced-motion:reduce) { .settings-btn svg, .settings-row input, .settings-row input::after { transition:none; } }
    /* Lite visuals (html.lite, see the first script): no animation, transition, blur or glow anywhere. Animations jump
       straight to their end state (after any delay, so a ring still shows green before it goes), so anything that fades
       out still ends hidden. Panels that relied on blur over the
       3D view get a near-solid background so their text stays readable. */
    html.lite { scroll-behavior:auto; }
    html.lite *, html.lite *::before, html.lite *::after { animation-duration:.01ms !important; animation-iteration-count:1 !important; transition-duration:.01ms !important; -webkit-backdrop-filter:none !important; backdrop-filter:none !important; }
    html.lite .tutorial-spotlight, html.lite .tutorial-also, html.lite .tutorial-spotlight.is-complete { box-shadow:none; }
    html.lite .tutorial-pointer { animation:none; transform:translate(var(--to-x),var(--to-y)); } /* as with reduce-motion: the pointer rests where to act */
    html.lite .tutorial-spotlight.is-complete::after { display:none; }
    html.lite .topbar { background:rgba(6,10,18,.97); }
    /* Open Settings sits above everything, the tour's rings, dim layer and card included (they are z-index 9997-10001). */
    .topbar.settings-open { z-index:10010; }
    html.lite .view-panel, html.lite .fs-sliders, html.lite .telemetry, html.lite .stability-chip, html.lite .star-mark span, html.lite .curve-hover .ch-box { background:rgba(9,14,24,.94); }
    a:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
    header { display:grid; grid-template-columns:minmax(0,1fr) auto; align-items:end; gap:28px 40px; max-width:1800px; margin:auto; padding:48px 28px 8px; border-bottom:0; background:none; }
    .eyebrow { color:var(--accent); font:500 11px/1.4 var(--f-mono); letter-spacing:.14em; }
    h1 { margin:12px 0; font:600 clamp(28px,3.4vw,40px)/1.08 var(--f-serif); letter-spacing:-.012em; color:var(--ink); }
    .lead { max-width:62ch; margin:0; color:var(--ink-2); font-size:15px; line-height:1.6; }
    .header-stats { gap:10px; align-items:stretch; }
    .stat, .stat.ml { padding:10px 14px 11px; border:1px solid var(--line); border-radius:var(--r); background:var(--panel); }
    .stat-label { font-size:10px; letter-spacing:.13em; color:var(--muted); }
    .stat-value, .stat.ml .stat-value { margin-top:4px; font:500 19px/1.15 var(--f-mono); color:var(--ink); }
    .stat-value small { font:400 11px var(--f-ui); }
    .ml-pair span { font-family:var(--f-ui); }
    .ml-info { border-color:var(--line-strong); color:var(--ink-2); font-family:var(--f-ui); }
    .ml-info:hover, .ml-info[aria-expanded="true"] { background:rgba(138,180,248,.14); color:var(--ink); }
    .ml-explainer { border-color:var(--line-strong); border-radius:var(--r); background:var(--panel-2); }
    .ml-explainer h3 { color:var(--ink); font-family:var(--f-serif); font-size:15px; }
    .fit-toggle { border-radius:var(--r); background:var(--panel); }
    main { padding:8px 28px 48px; }
    .section-head { display:grid; grid-template-columns:auto minmax(0,1fr); column-gap:14px; align-items:baseline; margin:36px 0 14px; padding-top:20px; border-top:1px solid var(--line); }
    main > .section-head:first-child { margin-top:20px; }
    .section-head .kicker { grid-row:span 2; color:var(--accent); font:500 12px var(--f-mono); letter-spacing:.06em; }
    .section-head h2 { font:600 22px/1.2 var(--f-serif); letter-spacing:-.005em; }
    .section-head p { grid-column:2; margin:3px 0 0; color:var(--muted); font-size:13px; }
    .model-bar, .controls, .sandbox, .panel, .learn article { border:1px solid var(--line); border-radius:var(--r-lg); background:linear-gradient(180deg,var(--panel-2),var(--panel)); box-shadow:var(--shadow); }
    .model-bar { padding:12px 18px; margin-bottom:12px; gap:14px 20px; }
    .model-bar .model-label { font:500 10.5px var(--f-mono); letter-spacing:.14em; }
    .controls { gap:22px 30px; padding:20px 22px; margin-bottom:0; }
    .control label { color:var(--ink-2); font-size:12.5px; font-weight:500; }
    .control output, .fs-sliders output { font-family:var(--f-mono); font-size:12px; color:var(--ink); white-space:nowrap; flex:none; }
    .control-help { font-size:12px; line-height:1.5; }
    /* Layout: the halo sliders sit directly above the plots they change (layout choice only). */
    .sim-stage .controls { grid-template-columns:auto minmax(0,1.35fr) repeat(3,minmax(0,1fr)); gap:14px 26px; padding:14px 20px; margin-bottom:14px; }
    .sim-stage .control-help { font-size:11.5px; line-height:1.45; min-height:0; margin-top:4px; }
    .sim-stage .control label { flex-wrap:wrap; row-gap:2px; } .sim-stage .control label > span { flex:1 0 auto; }
    @media (min-width:961px) and (max-width:1280px) { .sim-stage .play-btn .play-text { display:none; } .sim-stage .play-btn { width:44px; padding:0; } } /* room for the time slider in narrow columns */
    /* Phones: 3D view, then the rotation curve pinned under the menu, then the sliders, so a slider's effect on the
       curve stays in view while dragging (layout choice only). Unpinned during the tutorial, whose card needs the room. */
    @media (max-width:640px) {
      .sim-stage { display:flex; flex-direction:column; gap:14px; }
      .sim-stage .panel-grid { display:contents; }
      .sim-stage .panel-grid > .panel:first-child { order:1; }
      .sim-stage > .sandbox, #controlsAnchor { order:3; }
      .sim-stage > .controls { order:4; margin-bottom:0; }
      #curvePanel { order:2; position:sticky; top:56px; z-index:40; max-height:calc(100dvh - 56px - 180px); overflow-y:auto; overscroll-behavior:contain; box-shadow:0 12px 24px rgba(0,0,0,.6); }
      #curvePlot.plot { height:34dvh; min-height:220px; }
      body.tutorial-active #curvePanel { position:relative; max-height:none; overflow:hidden; }
    }
    /* Laptops and projectors: the slider bar scrolls with the page (not pinned), and each figure fits the screen below the
       menu, so both figures can be seen whole at once. Once the bar is scrolled away, the slider box (.fs-sliders) appears
       in the corner of the 3D view, so a slider can always be moved while watching both figures (layout choice only). */
    @media (min-width:961px) {
      .sim-stage .plot { height:clamp(360px, calc(100vh - 56px - 96px), 700px); min-height:0; }
    }
    /* Narrower laptops: the gravity model takes its own row, so the four slider columns have room for their values. */
    @media (min-width:961px) and (max-width:1279px) {
      .sim-stage .controls { grid-template-columns:repeat(4,minmax(0,1fr)); }
      .sim-stage .controls:has(#massControl[hidden]) { grid-template-columns:repeat(2,minmax(0,1fr)); }
      .sim-stage .gravity-control { grid-column:1/-1; display:flex; align-items:center; flex-wrap:wrap; gap:8px 18px; }
      .sim-stage .gravity-control .model-cards { flex:0 1 460px; }
      .sim-stage .gravity-control .control-help { flex:1 1 240px; margin:0; }
    }
    @media (max-width:960px) { .sim-stage .controls { grid-template-columns:1fr 1fr; } .sim-stage .control:last-child { grid-column:auto; } }
    @media (max-width:560px) { .sim-stage .controls { grid-template-columns:1fr; } }
    input[type=range] { -webkit-appearance:none; appearance:none; height:22px; margin:0; background:transparent; cursor:pointer; }
    input[type=range]::-webkit-slider-runnable-track { height:4px; border-radius:999px; background:rgba(148,170,200,.22); }
    input[type=range]::-moz-range-track { height:4px; border-radius:999px; background:rgba(148,170,200,.22); }
    input[type=range]::-moz-range-progress { height:4px; border-radius:999px; background:var(--accent); }
    input[type=range]::-webkit-slider-thumb { -webkit-appearance:none; width:16px; height:16px; margin-top:-6px; border-radius:50%; background:#eef3fa; border:2px solid var(--accent); box-shadow:0 1px 4px rgba(0,0,0,.5); transition:box-shadow .15s; }
    input[type=range]::-moz-range-thumb { width:12px; height:12px; border-radius:50%; background:#eef3fa; border:2px solid var(--accent); box-shadow:0 1px 4px rgba(0,0,0,.5); }
    input[type=range]:hover::-webkit-slider-thumb, input[type=range]:focus-visible::-webkit-slider-thumb { box-shadow:0 0 0 5px rgba(138,180,248,.2); }
    input[type=range]:focus-visible { outline:none; }
    input[type=range]:focus-visible::-moz-range-thumb { box-shadow:0 0 0 5px rgba(138,180,248,.25); }
    input[type=range]:disabled { cursor:not-allowed; }
    button { padding:7px 13px; border:1px solid var(--line-strong); border-radius:8px; background:rgba(255,255,255,.035); color:var(--ink); font:500 12.5px/1.2 var(--f-ui); transition:background .15s,border-color .15s,color .15s; }
    button:hover { border-color:rgba(138,180,248,.6); background:rgba(138,180,248,.08); }
    button:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
    #applyObservedFit, .tutorial-primary { background:rgba(138,180,248,.16); border-color:rgba(138,180,248,.6); color:#eef4ff; }
    #applyObservedFit:hover, .tutorial-primary:hover { background:rgba(138,180,248,.26); }
    .preset-picker { flex:1 1 420px; min-width:0; gap:12px 14px; }
    .galaxy-menu { position:relative; }
    .galaxy-menu-btn { display:flex; align-items:center; gap:10px; min-width:200px; padding:6px 36px 6px 7px; background-color:var(--panel-2); background-image:linear-gradient(45deg,transparent 50%,var(--muted) 50%),linear-gradient(135deg,var(--muted) 50%,transparent 50%); background-position:calc(100% - 19px) 52%,calc(100% - 14px) 52%; background-size:5px 5px,5px 5px; background-repeat:no-repeat; font:600 14px var(--f-ui); text-align:left; }
    .galaxy-menu-btn[aria-expanded="true"] { border-color:var(--accent); }
    .gm-code { display:inline-grid; place-items:center; flex:none; min-width:36px; height:24px; padding:0 6px; border:1px solid rgba(233,200,138,.42); border-radius:6px; background:rgba(233,200,138,.08); color:var(--gold); font:600 11.5px var(--f-mono); }
    .galaxy-menu-list { position:absolute; top:calc(100% + 6px); left:0; z-index:60; width:min(400px,calc(100vw - 32px)); margin:0; padding:6px; list-style:none; border:1px solid var(--line-strong); border-radius:var(--r); background:var(--panel-2); box-shadow:0 18px 44px rgba(0,0,0,.55); }
    .galaxy-menu-list[hidden] { display:none; }
    .galaxy-menu-list:focus { outline:none; }
    .galaxy-menu-list li { display:flex; align-items:flex-start; gap:11px; padding:9px 10px; border-radius:8px; cursor:pointer; }
    .galaxy-menu-list li.active { background:rgba(138,180,248,.12); }
    .galaxy-menu-list li[aria-selected="true"] .gm-text b::after { content:" ✓"; color:var(--accent); }
    .gm-text { display:flex; flex-direction:column; gap:2px; min-width:0; }
    .gm-text b { color:var(--ink); font:600 13.5px var(--f-ui); }
    .gm-text small { color:var(--muted); font-size:12px; line-height:1.4; }
    .galaxy-profile { flex:1 1 100%; display:grid; grid-template-columns:auto minmax(0,1fr); gap:4px 18px; padding:14px 18px 13px; border:1px solid rgba(233,200,138,.26); border-left:3px solid var(--gold); border-radius:var(--r); background:linear-gradient(160deg,rgba(233,200,138,.075),rgba(233,200,138,.01)); }
    .gp-badge { display:grid; place-items:center; width:66px; height:66px; margin-top:2px; border:1px solid rgba(233,200,138,.38); border-radius:50%; background:radial-gradient(circle,rgba(233,200,138,.3) 0,rgba(233,200,138,.06) 55%,transparent 72%); color:var(--gold); font:600 20px var(--f-serif); }
    .galaxy-profile h3 { margin:2px 0 3px; color:var(--ink); font:600 18px/1.25 var(--f-serif); }
    .gp-summary { margin:0 0 11px; color:var(--ink-2); font-size:13px; line-height:1.45; }
    .gp-stats { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:10px 16px; margin:0 0 12px; }
    .gp-stats dt { color:var(--muted); font:500 10px var(--f-mono); letter-spacing:.12em; text-transform:uppercase; }
    .gp-stats dd { margin:3px 0 0; }
    .gp-stats b { display:block; color:var(--ink); font:500 17px/1.2 var(--f-mono); }
    .gp-stats small { color:var(--muted); font-size:11.5px; }
    .gp-mass-head { display:flex; justify-content:space-between; flex-wrap:wrap; gap:3px 14px; margin-bottom:6px; color:var(--ink-2); font-size:12px; }
    .gp-mass-head b { color:var(--ink); font-family:var(--f-mono); font-weight:500; }
    .gp-dot { display:inline-block; width:8px; height:8px; margin-right:5px; border-radius:50%; }
    .gp-dot.stars { background:var(--gold); } .gp-dot.gas { background:rgba(138,180,248,.75); }
    .gp-bar { display:flex; height:8px; overflow:hidden; border-radius:999px; background:rgba(138,180,248,.75); }
    .gp-bar span { background:var(--gold); transition:flex-basis .45s ease; }
    .gp-note { margin:8px 0 0; color:var(--muted); font-size:10.5px; line-height:1.4; }
    @media (prefers-reduced-motion:reduce) { .gp-bar span { transition:none; } }
    @media (min-width:1000px) { .gp-body { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1.45fr); gap:0 30px; align-items:start; } .gp-summary { margin-bottom:0; } }
    @media (max-width:560px) { .galaxy-profile { grid-template-columns:1fr; } .gp-badge { width:48px; height:48px; font-size:16px; } .gp-stats { grid-template-columns:1fr 1fr; } }
    .switch .opt { font-size:13px; }
    .track-input { border-color:var(--line-strong); background:rgba(0,0,0,.25); font-family:var(--f-mono); }
    .sandbox-group h3 { font:500 10.5px var(--f-mono); letter-spacing:.14em; color:var(--accent); }
    .sandbox-summary, .key, .halo-info, .stability, .galaxy-type { border-radius:var(--r); }
    .sandbox-summary dd { font-family:var(--f-mono); }
    .galaxy-type h3 { font:600 18px var(--f-serif); }
    .panel-grid { gap:16px; }
    .panel-heading { align-items:center; flex-wrap:wrap; row-gap:6px; padding:16px 18px 8px; }
    .panel-heading h2 { font:600 15px var(--f-ui); color:var(--ink); }
    .panel-heading .fig { margin-right:10px; color:var(--accent); font:500 11px var(--f-mono); letter-spacing:.08em; text-transform:uppercase; }
    .panel-heading > span:not(#chiBadge) { color:var(--muted); font-size:12px; }
    #chiBadge { font:500 11.5px var(--f-mono); white-space:nowrap; }
    .key { background:rgba(148,170,200,.035); }
    .key summary { font:500 10.5px var(--f-mono); letter-spacing:.12em; color:var(--ink-2); }
    .sw.ldash { background:none; height:0; border-top:2.5px dashed; }
    .halo-info h3 { font:600 14px var(--f-ui); }
    .halo-metric .metric-value { font:500 17px var(--f-mono); }
    .halo-metric button, .halo-metric button:hover { background:none; border:0; cursor:default; }
    .stability-c { font-family:var(--f-mono); }
    .learn { gap:16px; margin-top:0; }
    .learn article { padding:0 22px 8px; }
    .learn h3 { padding:16px 0; font:600 18px var(--f-serif); }
    .learn p { color:var(--ink-2); font-size:13.5px; line-height:1.65; }
    .learn p.tldr { color:var(--ink); font-size:14px; }
    .schematic { margin:4px 0 16px; padding:12px 14px 10px; border:1px solid var(--line); border-radius:var(--r); background:rgba(0,0,0,.22); color:var(--muted); }
    .schematic svg { display:block; width:100%; max-width:440px; height:auto; margin:auto; font:10.5px var(--f-ui); }
    .schematic text:not([fill]) { fill:currentColor; }
    .schematic figcaption { margin-top:6px; font-size:11.5px; line-height:1.5; }
    .eq { display:block; margin:8px 0 12px; padding:10px 12px; border-left:2px solid var(--line-strong); border-radius:var(--r); background:rgba(0,0,0,.22); overflow-x:auto; color:var(--ink); }
    .eq math { font-family:"STIX Two Math","Cambria Math","Latin Modern Math",math,serif; font-size:1.25em; margin:0; }
    .method-grid .eq { margin:6px 0; padding:6px 10px; }
    .mini-table { border-collapse:collapse; margin:6px 0 12px; font-size:.92em; font-variant-numeric:tabular-nums; }
    .mini-table th, .mini-table td { padding:3px 12px 3px 0; text-align:left; border-bottom:1px solid var(--line); }
    .mini-table th { color:var(--muted); font-weight:500; }
    .method-grid .eq math { font-size:1.1em; }
    footer { flex-wrap:wrap; align-items:baseline; gap:8px 24px; margin:28px 0 0; padding:18px 22px; border:1px solid var(--line); border-radius:var(--r-lg); background:var(--panel); font-size:12px; }
    .footer-title { flex-basis:100%; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    footer { display:block; margin-top:0; }
    /* Two flowing columns (read down the left, then the right), so a short card leaves no hole beside a long one. */
    .method-grid { columns:2; column-gap:32px; margin:0 0 16px; }
    @media (max-width:900px) { .method-grid { columns:1; } }
    .method-grid > div { min-width:0; margin:0 0 22px; break-inside:avoid; } .method-grid dd p { margin:6px 0 0; } .method-grid dd ul { margin:6px 0 0; padding-left:18px; } .method-grid dd li { margin:3px 0; }
    .method-grid dt { margin-bottom:4px; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    .method-grid dd { margin:0; color:var(--ink-2); font-size:13.5px; line-height:1.6; }
    .tutorial-welcome, .tutorial-card { border-color:rgba(138,180,248,.24); background:linear-gradient(155deg,#121b2e,#0a101c 70%); }
    .tutorial-welcome h2, .tutorial-card h2 { font-family:var(--f-serif); }
    .tutorial-eyebrow { font-family:var(--f-mono); }
    .stat-q::after { background:var(--panel-2); border-color:var(--line-strong); }
    .telemetry { font-family:var(--f-mono); }
    @media (max-width:900px) { header { grid-template-columns:1fr; padding-top:36px; } .header-stats { justify-content:flex-start; flex-wrap:wrap; } }
    @media (max-width:640px) { .topbar { padding:6px 12px; gap:8px; } .brand span { display:none; } .topnav a { padding:6px 6px; font-size:12px; } header { padding:28px 16px 4px; } main { padding:0 12px 32px; } .section-head { margin-top:28px; } .section-head h2 { font-size:19px; } .controls { padding:16px; } .learn article { padding:0 16px 6px; } footer { padding:16px; } .stat { flex:1 1 140px; } .panel-heading > span:not(#chiBadge) { display:none; } }
    .key summary::before { content:""; display:inline-block; width:4px; }
    /* Small text: sizes and contrast raised for readability (visual only). */
    :root { --muted:#96a5ba; }
    .stat-label { font-size:10.5px; }
    .sim-stage .control-help { font-size:12px; }
    .schematic figcaption { font-size:12px; }
    footer .export-btn { font-size:12px; }
    /* The puzzle: the story before the controls. */
    .puzzle, .scorecard { border:1px solid var(--line); border-radius:var(--r-lg); background:linear-gradient(180deg,var(--panel-2),var(--panel)); box-shadow:var(--shadow); }
    .puzzle { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1.05fr); gap:8px 36px; align-items:center; padding:22px 24px; }
    .puzzle-steps { margin:0; padding:0; list-style:none; counter-reset:puzzle; }
    .puzzle-steps li { position:relative; padding:0 0 14px 40px; counter-increment:puzzle; }
    .puzzle-steps li:last-child { padding-bottom:0; }
    .puzzle-steps li::before { content:counter(puzzle); position:absolute; left:0; top:0; width:26px; height:26px; border:1px solid var(--line-strong); border-radius:50%; color:var(--accent); font:500 12px/24px var(--f-mono); text-align:center; }
    .puzzle-steps li:not(:last-child)::after { content:""; position:absolute; left:13px; top:32px; bottom:4px; width:1px; background:var(--line-strong); }
    .puzzle-steps h3 { margin:0 0 3px; padding-top:2px; color:var(--ink); font:600 16px/1.3 var(--f-serif); }
    .puzzle-steps p { margin:0; color:var(--ink-2); font-size:13.5px; line-height:1.6; }
    .puzzle-steps b { color:var(--ink); font-weight:600; }
    .puzzle-fig { margin:0; }
    .puzzle-fig svg { display:block; width:100%; max-width:520px; height:auto; margin:auto; font:13px var(--f-ui); }
    .puzzle-fig figcaption { margin-top:8px; color:var(--muted); font-size:12px; line-height:1.5; }
    .puzzle-cta { grid-column:1/-1; display:flex; flex-wrap:wrap; align-items:center; gap:10px 18px; margin-top:10px; padding-top:16px; border-top:1px solid var(--line); }
    .puzzle-cta button { padding:9px 16px; font-size:13px; }
    .puzzle-cta a, .refs a, .evidence-source a { color:var(--accent); text-underline-offset:3px; }
    .puzzle-cta a { font-size:13px; }
    @media (max-width:900px) { .puzzle { grid-template-columns:1fr; gap:18px; padding:18px 16px; } }
    /* Fit badge: a plain-language verdict in front of the numbers. */
    #chiBadge sub { font-size:.82em; }
    .chi-wrap { display:flex; align-items:center; flex-wrap:wrap; justify-content:flex-end; gap:4px 8px; }
    .chi-fit { color:var(--muted); font-size:11.5px; cursor:help; }
    #chiBadge b { margin-right:4px; font:600 11.5px var(--f-ui); letter-spacing:.01em; }
    /* Wide screens: a long fit label wraps inside its own space and both figure headings keep one fixed height,
       so the graph never jumps down when the label grows (e.g. "Visible matter only · 0 tuned · about 28 error bars off"). */
    @media (min-width:1200px) {
      .sim-stage .panel-heading { flex-wrap:nowrap; min-height:58px; }
      .sim-stage .panel-heading h2, #chiBadge { flex:none; }
      .chi-wrap { flex:1 1 0; min-width:0; flex-wrap:nowrap; }
      .chi-fit { min-width:0; text-align:right; line-height:1.3; }
    }
    @media (min-width:1200px) { .sim-stage .plot { height:clamp(360px, calc(100vh - 56px - 106px), 700px); } }
    /* Scorecard: every preset galaxy scored under each model. */
    .scorecard { padding:18px 22px 16px; }
    .scorecard > p { max-width:88ch; margin:0 0 14px; color:var(--ink-2); font-size:13.5px; line-height:1.65; }
    .score-scroll { overflow-x:auto; margin:0 -6px; padding:0 6px; }
    .score-table { width:100%; min-width:640px; border-collapse:collapse; font-size:13px; }
    .score-table th, .score-table td { padding:10px 12px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }
    .score-table thead th { color:var(--ink); font-weight:600; vertical-align:bottom; }
    .score-table thead th.score-group { color:var(--accent); border-bottom:1px solid var(--line-strong); }
    .score-table small { display:block; margin-top:2px; color:var(--muted); font-size:11.5px; font-weight:400; line-height:1.4; }
    .score-table tbody th { color:var(--ink); font-weight:600; white-space:nowrap; }
    .score-table tbody tr:last-child th, .score-table tbody tr:last-child td { border-bottom:0; }
    /* Phones: one stacked card per galaxy, each score labelled with its column, so nothing is off-screen. */
    @media (max-width:640px) {
      .score-table { min-width:0; }
      .score-table thead { display:none; }
      .score-table tbody, .score-table tr, .score-table tbody th, .score-table td { display:block; }
      .score-table tr { margin-bottom:14px; border:1px solid var(--line); border-radius:10px; }
      .score-table tbody tr:last-child th, .score-table tbody tr:last-child td { border-bottom:1px solid var(--line); }
      .score-table tr > :last-child, .score-table tbody tr:last-child td:last-child { border-bottom:0; }
      .score-table tbody th { white-space:normal; background:rgba(138,180,248,.06); }
      .score-table td[data-label]::before { content:attr(data-label); display:block; margin-bottom:4px; color:var(--ink); font-size:12px; font-weight:600; }
    }
    .score { display:inline-block; min-width:58px; padding:2px 9px; border:1px solid rgba(120,150,180,.3); border-radius:999px; font:500 13px var(--f-mono); font-variant-numeric:tabular-nums; text-align:center; }
    .score.good { border-color:rgba(120,214,140,.55); color:#93e3a0; } .score.ok { border-color:rgba(236,204,96,.55); color:#f2d675; } .score.poor { border-color:rgba(255,138,128,.6); color:#ff9d94; }
    .at-limit { padding:0 4px; border:1px solid rgba(255,174,92,.5); border-radius:4px; color:#ffc48a; font-size:10.5px; text-decoration:none; cursor:help; white-space:nowrap; }
    .score-show { margin-top:6px; padding:3px 9px; font-size:11.5px; }
    .glance-table td, .glance-table th { vertical-align:middle; }
    .glance-name { margin:0; font-weight:600; color:var(--ink); }
    .glance-pair { display:inline-flex; align-items:center; gap:6px; white-space:nowrap; }
    .glance-pair sup { color:#ffc48a; cursor:help; }
    .glance-arrow { color:var(--muted); }
    .glance-table .score { min-width:54px; }
    .glance-note { margin:6px 0 0; color:var(--muted); font-size:12px; } .glance-note sup { color:#ffc48a; }
    .score-details { margin-top:16px; }
    .score-details > summary { cursor:pointer; color:var(--accent); font-size:13px; font-weight:500; padding:6px 0; }
    .score-details[open] > summary { margin-bottom:10px; }
    .score-notes { margin:14px 0 0; padding:14px 0 0 18px; border-top:1px solid var(--line); color:var(--ink-2); font-size:13px; line-height:1.6; }
    .score-notes li { margin-bottom:6px; max-width:96ch; }
    .score-notes b { color:var(--ink); font-weight:600; }
    .score-notes-title { margin:24px 0 0; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    .score-notes-title + .score-notes { margin-top:8px; padding-top:0; border-top:0; }
    /* Wider evidence: a full-width box under the two theory boxes. */
    .learn article.evidence { grid-column:1/-1; }
    .evidence-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:14px; margin:6px 0 14px; }
    .evidence-card { padding:12px 14px; border:1px solid var(--line); border-radius:var(--r); background:rgba(0,0,0,.22); }
    .evidence-card h4 { margin:6px 0 6px; color:var(--ink); font:600 14.5px/1.3 var(--f-serif); }
    .learn .evidence-card p { margin:0 0 8px; font-size:13px; line-height:1.6; }
    .learn .evidence-card p.evidence-source { margin:0; color:var(--muted); font-size:11.5px; }
    .evidence-tag { display:inline-block; padding:2px 8px; border-radius:999px; font:500 10.5px var(--f-mono); letter-spacing:.06em; }
    .evidence-tag.dm { border:1px solid rgba(179,154,232,.5); color:#c9b8f0; } .evidence-tag.mond { border:1px solid rgba(255,174,92,.5); color:#ffc48a; }
    @media (max-width:1100px) { .evidence-grid { grid-template-columns:1fr 1fr; } } @media (max-width:560px) { .evidence-grid { grid-template-columns:1fr; } }
    /* Reference list. */
    .refs { margin:56px 0 18px; padding-top:28px; border-top:1px solid var(--line-strong); }
    .refs h3 { margin:0 0 14px; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    .refs ul { margin:0; padding:0; list-style:none; display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:12px; }
    .refs li { margin:0; padding:12px 14px; border:1px solid var(--line); border-left:2px solid var(--line-strong); border-radius:var(--r); background:rgba(0,0,0,.16); color:var(--ink-2); font-size:12.5px; line-height:1.6; overflow-wrap:anywhere; }
    .refs li a { display:block; width:fit-content; margin-top:6px; font:500 11px var(--f-mono); }
    .refs p { margin:8px 0 0; color:var(--muted); font-size:12px; }
    @media (max-width:900px) { .refs ul { grid-template-columns:1fr; } }
    @media (max-width:640px) { .scorecard { padding:16px; } }
    /* Pages: one .page shows at a time. A hidden page keeps its width (height 0, not display:none) so its plots stay sized. */
    .page:not(.active) { height:0; overflow:hidden; visibility:hidden; }
    .load-failed { margin:16px 0; padding:14px 16px; border:1px solid rgba(255,138,128,.6); border-radius:10px; background:rgba(255,138,128,.08); color:var(--ink); line-height:1.5; }
    body.sim-unavailable .page:not(.active) { height:auto; overflow:visible; visibility:visible; }
    body.sim-unavailable :is(.page[data-page="simulator"], .page[data-page="scorecard"], .page-next, .topnav a[href="#simulator"], .topnav a[href="#scorecard"]) { display:none; }
    .topnav a[aria-current="page"] { color:var(--ink); background:rgba(138,180,248,.16); }
    .page header { max-width:none; margin:0; padding:36px 0 4px; display:block; }
    .page > .section-head:first-child, .sim-head > .section-head { margin-top:20px; padding-top:8px; border-top:0; }
    .sim-head { display:grid; grid-template-columns:minmax(0,1fr) auto; align-items:end; gap:8px 40px; margin-bottom:14px; }
    .sim-head > .section-head { margin-bottom:0; }
    @media (max-width:900px) { .sim-head { grid-template-columns:1fr; } }
    .page-next { display:flex; justify-content:flex-end; margin-top:22px; }
    .next-link { display:inline-block; padding:9px 16px; border:1px solid rgba(138,180,248,.5); border-radius:8px; background:rgba(138,180,248,.08); color:#dbe7ff; font:500 13px/1.2 var(--f-ui); text-decoration:none; transition:background .15s,border-color .15s; }
    .next-link:hover, .next-link:focus-visible { background:rgba(138,180,248,.18); border-color:rgba(138,180,248,.8); }
    .puzzle-cta span { color:var(--ink-2); font-size:13.5px; }
    .answers { display:grid; grid-template-columns:1fr 1fr; gap:10px; margin-top:8px; }
    .answer { padding:10px 12px; border:1px solid var(--line); border-left-width:3px; border-radius:var(--r); background:rgba(0,0,0,.22); }
    .answer b { display:block; margin-bottom:2px; font:600 14px var(--f-ui); }
    .answer span { color:var(--ink-2); font-size:13px; line-height:1.5; }
    .answer.dm { border-left-color:#b39ae8; } .answer.dm b { color:#c9b8f0; }
    .answer.mond { border-left-color:#ffae5c; } .answer.mond b { color:#ffc48a; }
    @media (max-width:560px) { .answers { grid-template-columns:1fr; } }
    .score-legend { display:flex; flex-wrap:wrap; gap:8px 22px; margin:0 0 14px; padding:0; list-style:none; color:var(--ink-2); font-size:12.5px; }
    .score-legend .score { min-width:0; margin-right:8px; font-size:11.5px; }
    .score-notes-title + .score-notes { display:grid; grid-template-columns:repeat(3,1fr); gap:12px; padding:0; list-style:none; }
    /* Six cards: 3 + 3 on wide screens, 2 + 2 + 2 on tablets, one column on phones, so no row is left half empty. */
    @media (max-width:1000px) { .score-notes-title + .score-notes { grid-template-columns:repeat(2,1fr); } }
    @media (max-width:640px) { .score-notes-title + .score-notes { grid-template-columns:1fr; } }
    .score-notes-title + .score-notes li { max-width:none; margin:0; padding:12px 14px; border:1px solid var(--line); border-radius:var(--r); background:rgba(0,0,0,.22); font-size:12.5px; line-height:1.55; }
    .score-notes-title + .score-notes b { display:block; margin-bottom:4px; font:600 13.5px var(--f-serif); }
    .evidence-grid.two { grid-template-columns:1fr 1fr; margin-bottom:6px; }
    @media (max-width:700px) { .evidence-grid.two { grid-template-columns:1fr; } }
    .learn .evidence-card p { margin:0; } .learn .evidence-card p + p { margin-top:8px; }
    .evidence-open-title { margin:4px 0 0; color:var(--accent); font:500 10.5px var(--f-mono) !important; letter-spacing:.14em; text-transform:uppercase; }
    .maths { margin:0 0 12px; border:1px solid var(--line); border-radius:var(--r); background:rgba(0,0,0,.16); }
    .maths summary { padding:9px 12px; color:var(--accent); font:500 12.5px var(--f-ui); cursor:pointer; }
    .maths[open] { padding:0 12px 2px; } .maths[open] summary { padding:9px 0; }
  </style>
</head>
<body>
  <nav class="topbar" aria-label="Page sections">
    <a class="brand" href="#puzzle"><svg viewBox="0 0 32 32" aria-hidden="true"><ellipse cx="16" cy="16" rx="13" ry="5.5" transform="rotate(-24 16 16)" fill="none" stroke="currentColor" stroke-width="1.4" opacity=".55"/><ellipse cx="16" cy="16" rx="8" ry="3.2" transform="rotate(-24 16 16)" fill="none" stroke="currentColor" stroke-width="1.4" opacity=".85"/><circle cx="16" cy="16" r="2.4" fill="currentColor"/></svg><span>Galaxy Rotation Simulator</span></a>
    <div class="topnav"><a href="#puzzle" aria-current="page"><span class="nav-long">The puzzle</span><span class="nav-short">Puzzle</span></a><a href="#simulator">Simulator</a><a href="#scorecard">Scorecard</a><a href="#background">Background</a><a href="#method">Method<span class="nav-long"> &amp; sources</span></a></div>
    <div class="settings">
      <button class="settings-btn" id="settingsBtn" type="button" aria-label="Settings" title="Settings" aria-expanded="false" aria-controls="settingsPanel"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><circle cx="12" cy="12" r="3"/><path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 1 1-2.83 2.83l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 1 1-4 0v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 1 1-2.83-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 1 1 0-4h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 1 1 2.83-2.83l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 1 1 4 0v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 1 1 2.83 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 1 1 0 4h-.09a1.65 1.65 0 0 0-1.51 1z"/></svg></button>
      <div class="settings-panel" id="settingsPanel" role="dialog" aria-label="Settings" tabindex="-1" hidden>
        <h2>Settings</h2>
        <label class="settings-row"><input type="checkbox" id="qualityToggle" disabled><span><b>High-quality visuals</b><small>Glow, pulsing highlights, blur and motion, and a sharper 3D galaxy view. Untick it on a slow computer: everything gets simpler and runs smoother.</small></span></label>
        <p class="settings-note"><span id="qualityNote" aria-live="polite">The 3D view has not loaded.</span><button type="button" id="qualityAuto" hidden>Use automatic</button></p>
      </div>
    </div>
  </nav>
  <main>
    <div class="page active" data-page="puzzle">
    <header id="top">
    <div class="intro"><div class="eyebrow">SPARC galaxy dynamics · interactive model</div><h1>Interactive Galaxy Rotation Simulator</h1><p class="lead">Far from a galaxy’s centre, stars and gas orbit much faster than the gravity of everything we can see could explain. Test two explanations against real measurements: an unseen halo of dark matter, or a change to how gravity works.</p></div>
    </header>
    <div class="section-head" id="puzzle"><span class="kicker">01</span><h2>The puzzle</h2><p>Why galaxies seem to hold far more matter than we can see.</p></div>
    <section class="puzzle" aria-label="The missing mass problem">
      <ol class="puzzle-steps">
        <li><h3>A rotation curve is a speed chart</h3><p>Stars and gas circle the centre of their galaxy. Plot how fast they move against how far out they are, and you have a rotation curve.</p></li>
        <li><h3>What astronomers expected</h3><p>Gravity weakens with distance: Mercury orbits the Sun at 47 km/s, distant Neptune at only 5 km/s (Williams 2024). Most of a galaxy’s visible matter sits near its centre, so its outer stars and gas should be slow too: the dashed line.</p></li>
        <li><h3>What they measured</h3><p>The speeds stay high (Rubin, Ford &amp; Thonnard 1980). At the edge of __DEFAULT_GALAXY_NAME__, gas orbits at <b id="puzzleObs">–</b>. The gravity of the stars and gas we can see could only hold it in orbit at about <b id="puzzleBar">–</b>: any faster and it should fly away. To keep it circling that fast takes about <b id="puzzleRatio">–</b> times as much pull as the visible matter supplies.</p></li>
        <li><h3>Two possible answers</h3>
          <div class="answers">
            <div class="answer dm"><b>Dark matter</b><span>The galaxy sits inside a huge cloud of unseen mass.</span></div>
            <div class="answer mond"><b>MOND</b><span>Gravity itself works differently when it is very weak.</span></div>
          </div>
        </li>
      </ol>
      <figure class="puzzle-fig">
        <div id="puzzleChart"></div>
        <figcaption>Real measurements of __DEFAULT_GALAXY_NAME__ (SPARC; Lelli, McGaugh &amp; Schombert 2016). The dashed line is Newton’s gravity acting on the galaxy’s stars and gas alone.</figcaption>
      </figure>
      <div class="puzzle-cta"><span>Test both answers against real measurements:</span><button type="button" class="tutorial-primary" id="puzzleTour" title="The simulator tutorial, then the sandbox tutorial">Take the full guided tour</button><a class="next-link" href="#simulator">Open the simulator →</a></div>
    </section>
    </div>
    <div class="page" data-page="simulator">
    <div class="sim-head">
    <div class="section-head" id="simulator"><span class="kicker">02</span><h2>Model setup</h2><p>Choose a galaxy to model.</p></div>
    <button type="button" class="tutorial-launch" id="tutorialLaunch">Simulator tutorial</button>
    </div>
    <section class="model-bar galaxy-bar" id="galaxyBar" aria-label="Galaxy source">
      <div class="source-toggle" id="sourceToggle" role="radiogroup" aria-label="Galaxy source">
        <label><input type="radio" name="galaxyMode" value="preset" checked autocomplete="off"><span class="source-card"><svg class="source-icon" viewBox="0 0 40 40" aria-hidden="true"><path d="M20 20c3-1 4-5 1-7s-9 0-10 5 3 11 10 11 13-5 12-13S24 3 15 5" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round"/><circle cx="20" cy="20" r="2.4" fill="currentColor"/></svg><span class="source-text"><small>Real data</small><b>Preset Galaxies</b><span>Test both ideas against real galaxies whose speeds astronomers have measured (SPARC).</span></span></span></label>
        <label><input type="radio" name="galaxyMode" value="custom" autocomplete="off"><span class="source-card"><svg class="source-icon" viewBox="0 0 40 40" aria-hidden="true"><g stroke="currentColor" stroke-width="2.2" stroke-linecap="round"><path d="M8 11h24M8 20h24M8 29h24" opacity=".45"/></g><g fill="currentColor"><circle cx="15" cy="11" r="3.4"/><circle cx="26" cy="20" r="3.4"/><circle cx="12" cy="29" r="3.4"/></g></svg><span class="source-text"><small>Build your own</small><b>Custom Sandbox</b><span>Design a galaxy’s stars, gas and bulge, and see the dark matter halo it would need.</span></span></span></label>
      </div>
      <div class="preset-picker" id="presetPicker">
        <span class="model-label" id="galaxyLabel">Galaxy</span>
        <div class="galaxy-menu" id="galaxyMenu">
          <button type="button" class="galaxy-menu-btn" id="galaxyMenuBtn" aria-haspopup="listbox" aria-expanded="false" aria-labelledby="galaxyLabel galaxyMenuBtn" title="Choose which SPARC galaxy to model"><span class="gm-code" id="galaxyMenuCode"></span><span id="galaxyMenuName"></span></button>
          <ul class="galaxy-menu-list" id="galaxyMenuList" role="listbox" aria-labelledby="galaxyLabel" tabindex="-1" hidden></ul>
        </div>
        <div class="galaxy-profile" id="galaxyProfile" aria-live="polite">
          <div class="gp-badge" aria-hidden="true"><span id="gpCode"></span></div>
          <div class="gp-body">
            <div>
            <p class="galaxy-type-kicker" id="gpFamily"></p>
            <h3 id="gpShape"></h3>
            <p class="gp-summary" id="gpSummary"></p>
            </div>
            <div>
            <dl class="gp-stats">
              <div><dt>Distance</dt><dd><b id="gpDistance"></b><small id="gpDistanceLy"></small></dd></div>
              <div><dt>Flat rotation speed</dt><dd><b id="gpSpeed"></b><small id="gpSpeedKmh"></small></dd></div>
              <div><dt>Measured out to</dt><dd><b id="gpExtent"></b><small id="gpPoints"></small></dd></div>
            </dl>
            <div class="gp-mass-head"><span>Visible mass, in billions of Suns</span><span><i class="gp-dot stars"></i>Stars <b id="gpStars"></b> &nbsp; <i class="gp-dot gas"></i>Gas <b id="gpGas"></b></span></div>
            <div class="gp-bar" role="img" id="gpBar"><span id="gpStarsBar"></span></div>
            <p class="gp-note">Source: SPARC Table 1 (Lelli, McGaugh &amp; Schombert 2016).</p>
            </div>
          </div>
        </div>
      </div>
    </section>
    <div class="section-head" id="results"><span class="kicker">03</span><h2>Simulation</h2><p>Move the sliders: the orbits and the rotation curve respond straight away.</p></div>
    <div class="stats-row">
    <div class="header-stats" aria-label="Simulation summary">
      <div class="ml-col">
        <div class="stat ml" id="starsStat"><div class="stat-label"><span>Test stars</span><button type="button" class="ml-info" id="starsInfoBtn" aria-expanded="false" aria-controls="starsExplainer" title="What are test stars?">?</button></div><div class="stat-value"><span id="starCountStat">800</span></div></div>
        <div class="ml-explainer" id="starsExplainer" role="dialog" aria-labelledby="starsExplainerTitle" aria-hidden="true">
          <button type="button" class="ml-close" id="starsExplainerClose" aria-label="Close">×</button>
          <h3 id="starsExplainerTitle">What are test stars?</h3>
          <p>The dots in the 3D view are <b>test stars</b>: they feel the model’s gravity but do not pull on each other, so each one shows how that gravity alone moves a star.</p>
          <ul>
            <li><b>They start on balanced orbits</b>, set by the starting galaxy’s gravity (for a real galaxy, its measured speeds).</li>
            <li><b>Change the gravity and they respond.</b> Too little pull and they drift outwards; beyond 1.5 × the galaxy’s edge they count as ejected.</li>
          </ul>
          <p>The number is a display choice, not a real count: 800 for NGC 3198, scaled with each galaxy’s visible mass (200 to 2400). A real galaxy has billions of stars.</p>
        </div>
      </div>
      <div class="stat"><div class="stat-label">Timeline</div><div class="stat-value">1 Billion Years</div></div>
      <div class="ml-col">
        <div class="ml-row">
        <div class="stat ml" id="mlStat" title="Mass-to-light ratios used to turn the stars’ light into mass"><div class="stat-label"><span>Mass-to-light (M/L)</span><button type="button" class="ml-info" id="mlInfoBtn" aria-expanded="false" aria-controls="mlExplainer" title="What is a mass-to-light ratio?">?</button></div><div class="stat-value ml-pair"><div><b id="mlDisk">0.50</b><span>DISK</span></div><div><b id="mlBulge">0.70</b><span>BULGE</span></div></div><small class="ml-source" id="mlSource">standard SPARC values</small></div>
        <div class="stat ml" id="tiltStat" title="The galaxy’s tilt (inclination) in use"><div class="stat-label"><span>Tilt</span><button type="button" class="ml-info" id="tiltInfoBtn" aria-expanded="false" aria-controls="tiltExplainer" title="What is the tilt?">?</button></div><div class="stat-value"><b id="tiltValue">–</b></div><small class="ml-source" id="tiltSource">standard SPARC value</small></div>
        </div>
        <div class="ml-explainer" id="mlExplainer" role="dialog" aria-labelledby="mlExplainerTitle" aria-hidden="true">
          <button type="button" class="ml-close" id="mlExplainerClose" aria-label="Close">×</button>
          <h3 id="mlExplainerTitle">What is a mass-to-light ratio?</h3>
          <p>Telescopes measure a galaxy's <b>light</b>, not its mass. The mass-to-light ratio (M/L, or Υ) converts one into the other: Suns' worth of mass per Sun's worth of light.</p>
          <ul>
            <li><b>Low M/L means young stars.</b> Young, massive blue stars are very bright for their mass.</li>
            <li><b>High M/L means old stars.</b> Old populations are mostly small, dim red stars: lots of mass, not much light.</li>
          </ul>
          <p>SPARC’s standard values (in 3.6 μm infrared light) are <b>0.5 for the disk</b> and <b>0.7 for the older bulge</b>. They are averages, so every real galaxy differs a little. A higher M/L lifts the visible-matter curve and leaves less for dark matter to explain.</p>
          <p>The <b>Best realistic fits</b> switch lets each model scale them by 0.6 to 1.6 times, the range stellar models allow; the line under the tile says which values are in use. The Custom Sandbox sets masses directly, so needs no M/L.</p>
          <p class="ml-source">Source: Lelli, McGaugh &amp; Schombert (2016), SPARC; McGaugh &amp; Schombert (2014).</p>
        </div>
        <div class="ml-explainer" id="tiltExplainer" role="dialog" aria-labelledby="tiltExplainerTitle" aria-hidden="true">
          <button type="button" class="ml-close" id="tiltExplainerClose" aria-label="Close">×</button>
          <h3 id="tiltExplainerTitle">What is the tilt?</h3>
          <p>We see every galaxy at an angle. The <b>tilt</b> (or <b>inclination</b>, i) is that angle: 0° is face-on, 90° is edge-on.</p>
          <ul>
            <li><b>Telescopes only see motion towards or away from us</b> (the Doppler shift), not sideways.</li>
            <li><b>So a tilted disk looks slower than it really spins.</b> True speed = measured speed ÷ sin i.</li>
          </ul>
          <p>SPARC estimates each galaxy’s tilt and has already corrected its speeds for it, but the estimate has an error of a few degrees. Nearly face-on galaxies are the most sensitive: there a small change in tilt changes every speed a lot.</p>
          <p>The <b>Best realistic fits</b> switch lets each model change the tilt within twice SPARC’s error; the line under the tile says which value is in use.</p>
          <p class="ml-source">Source: Lelli, McGaugh &amp; Schombert (2016), SPARC; Li et al. (2018).</p>
        </div>
      </div>
        <label class="fit-toggle" id="bestFitWrap" title="Off: every model uses SPARC’s standard M/L, distance and tilt. On: the halo and MOND each tune them within the same realistic limits (and MOND its a₀) for their lowest χ²."><span class="fit-label">Best realistic fits</span><input type="checkbox" id="bestFitToggle" autocomplete="off"><span class="track" aria-hidden="true"></span></label>
    </div>
    </div>
    <div class="sim-stage" id="simStage">
    <section class="sandbox" id="sandboxPanel" aria-labelledby="sandboxTitle" hidden>
      <div class="sandbox-head"><h2 id="sandboxTitle">Custom galaxy builder</h2></div>
      <div class="galaxy-type" id="galaxyType">
        <p class="galaxy-type-kicker">What you have built</p>
        <h3 id="gtName" aria-live="polite"></h3>
        <p id="gtText"></p>
        <p class="galaxy-type-example" id="gtExample"></p>
        <ul class="galaxy-type-tags" id="gtTags" aria-label="Galaxy traits"></ul>
        <p class="galaxy-type-note">A rough guide, not a strict rule: the type comes from the star mass, bulge-to-total ratio, gas fraction and disk surface density (Hubble 1926; Freeman 1970; Graham &amp; Worley 2008).</p>
      </div>
      <div class="sandbox-groups" id="sandboxGroups">
        <div class="sandbox-group">
          <h3>Stellar disk</h3>
          <div class="control"><label for="sbDiskMass"><span>Disk mass, M<sub>disk</sub></span><output id="sbDiskMassValue" for="sbDiskMass"></output></label><input id="sbDiskMass" type="range" min="1e9" max="1e11" step="1e8" value="2e10" autocomplete="off"><p class="control-help">Mass of the stars in the flat disk.</p></div>
          <div class="control"><label for="sbDiskScale"><span>Scale length, R<sub>d</sub></span><output id="sbDiskScaleValue" for="sbDiskScale"></output></label><input id="sbDiskScale" type="range" min="1.5" max="6" step="0.1" value="3.1" autocomplete="off"><p class="control-help">How quickly the disk fades outwards. Bigger means more spread out.</p></div>
          <div class="control"><label for="sbTruncation"><span>Truncation radius, R<sub>max</sub></span><output id="sbTruncationValue" for="sbTruncation"></output></label><input id="sbTruncation" type="range" min="15" max="40" step="0.5" value="30" autocomplete="off"><p class="control-help">Outer edge of the visible galaxy. Stars beyond 1.5 × R<sub>max</sub> count as ejected.</p></div>
        </div>
        <div class="sandbox-group">
          <h3>Bulge</h3>
          <div class="control"><label for="sbBulgeMass"><span>Bulge mass, M<sub>bulge</sub></span><output id="sbBulgeMassValue" for="sbBulgeMass"></output></label><input id="sbBulgeMass" type="range" min="0" max="5e10" step="1e8" value="0" autocomplete="off"><p class="control-help">Mass of the round ball of stars at the centre. Zero means no bulge.</p></div>
          <div class="control" id="sbBulgeRadiusControl"><label for="sbBulgeRadius"><span>Hernquist scale radius, r<sub>b</sub></span><output id="sbBulgeRadiusValue" for="sbBulgeRadius"></output></label><input id="sbBulgeRadius" type="range" min="0.2" max="2" step="0.1" value="1" autocomplete="off"><p class="control-help" id="sbBulgeRadiusHelp">Size of the bulge (Hernquist profile; not the halo's scale radius). Half its mass lies within about 2.4 × r<sub>b</sub>.</p></div>
        </div>
        <div class="sandbox-group">
          <h3>Gas</h3>
          <div class="control"><label for="sbGasMass"><span>Gas mass, M<sub>gas</sub></span><output id="sbGasMassValue" for="sbGasMass"></output></label><input id="sbGasMass" type="range" min="5e8" max="2e10" step="1e8" value="1.4e10" autocomplete="off"><p class="control-help">Mass of cold hydrogen gas, spread three times wider than the stars.</p></div>
          <h3 class="sandbox-stats-title">Stats</h3>
          <dl class="sandbox-summary" id="sandboxSummary" aria-label="Sandbox galaxy stats">
            <dt>Visible mass<span class="stat-q" tabindex="0" role="img" aria-label="Disk stars + bulge stars + gas: everything you could see with a telescope." data-tip="Disk stars + bulge stars + gas: everything you could see with a telescope.">?</span></dt><dd id="sbTotalMass"></dd>
            <dt>Bulge-to-total<span class="stat-q" tabindex="0" role="img" aria-label="Share of the stars that sit in the central bulge. 0% is a pure disk." data-tip="Share of the stars that sit in the central bulge. 0% is a pure disk.">?</span></dt><dd id="sbBulgeFraction"></dd>
            <dt>Gas fraction<span class="stat-q" tabindex="0" role="img" aria-label="Share of the visible mass that is gas rather than stars." data-tip="Share of the visible mass that is gas rather than stars.">?</span></dt><dd id="sbGasFraction"></dd>
            <dt>Typical halo M<sub>200</sub><span class="stat-q" tabindex="0" role="img" aria-label="The dark matter halo mass usually found around this many stars (Moster et al. 2013). The halo mass slider starts here." data-tip="The dark matter halo mass usually found around this many stars (Moster et al. 2013). The halo mass slider starts here.">?</span></dt><dd id="sbHaloBaseline"></dd>
            <dt>Test stars<span class="stat-q" tabindex="0" role="img" aria-label="Stars in the 3D simulation, split between the disk and the bulge." data-tip="Stars in the 3D simulation, split between the disk and the bulge.">?</span></dt><dd id="sbStarCount"></dd>
            <dt>Ejection radius<span class="stat-q" tabindex="0" role="img" aria-label="A star that moves beyond this distance counts as ejected: gravity failed to hold it." data-tip="A star that moves beyond this distance counts as ejected: gravity failed to hold it.">?</span></dt><dd id="sbEjectRadius"></dd>
          </dl>
        </div>
      </div>
      <div class="sandbox-footer"><p id="sandboxNote">Changing the galaxy resets the halo to its typical size. This galaxy is imaginary, so there are no measured speeds to compare with.</p><button type="button" id="sandboxReset" title="Return the anatomy sliders to an NGC 3198-like galaxy">Reset anatomy</button></div>
    </section>
    <div id="controlsAnchor" aria-hidden="true"></div>
    <section class="controls" aria-label="Model controls">
    <div class="control gravity-control" id="modelBar" role="group" aria-labelledby="gravityLabel"><span class="gravity-label" id="gravityLabel">Gravity model</span><input type="checkbox" id="mondToggle" hidden autocomplete="off"><div class="model-cards" role="radiogroup" aria-labelledby="gravityLabel"><label class="model-card is-dm" id="optNfw"><input type="radio" name="gravityModel" value="nfw" checked autocomplete="off"><span><b>Dark matter halo</b><small>Extra unseen mass</small></span></label><label class="model-card is-mond" id="optMond"><input type="radio" name="gravityModel" value="mond" autocomplete="off"><span><b>MOND</b><small>Gravity changes, no extra mass</small></span></label></div><p class="control-help" id="nfwNote">Choose MOND to try changing gravity instead of adding dark matter.</p><p class="control-help mond-note" id="mondNote" hidden>MOND on: no dark matter. MOND boosts gravity everywhere: a little where it is strong, a lot where it is weaker than a<sub>0</sub> = <span id="mondA0Text">1.2×10⁻¹⁰</span> m/s².<span id="mondZoneNote" hidden> On the graph, the gold is in proportion to the boost: faint gold is a small boost, deep gold a big one.</span></p></div>
    <div class="control" id="massControl"><label for="massSlider"><span>Halo mass, M<sub>200</sub></span><output id="massValue"></output></label><input id="massSlider" type="range" min="0.00" max="2.00" step="0.01"><p class="control-help" id="massHelp">How much dark matter surrounds the galaxy. More mass, stronger gravity, faster orbits: the orange curve rises. 0 means visible matter only.</p><div class="fit-actions"><button id="applyObservedFit" class="fit-button" type="button">Fit halo to data</button><button id="visibleMatterOnly" class="fit-button" type="button" title="Set the dark matter halo mass to zero">Visible matter only</button></div></div>
    <div class="control" id="scaleControl"><label for="scaleSlider"><span>Scale radius, r<sub>s</sub></span><output id="scaleValue"></output></label><input id="scaleSlider" type="range" min="5" max="60" step="0.5"><p class="control-help" id="scaleHelp">How spread out the halo is, for the same mass. Small packs it near the centre and speeds up the inner stars; large spreads it out.</p></div>
      <div class="control"><label for="timeSlider"><span>Orbit time</span><output id="timeValue"></output></label><div class="timeline-row"><button type="button" id="playPauseBtn" class="play-btn" aria-pressed="false" aria-label="Play simulation" title="Play / pause (Space)"><span class="play-icon" aria-hidden="true">▶</span><span class="play-text">Play</span></button><input id="timeSlider" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"></div><p class="control-help">Moves the stars through one billion years of orbits. Drag it to any moment, or press Play.</p></div>
      <div class="control"><label for="trackInput"><span>Track a star</span><output id="trackMax"></output></label><input id="trackInput" class="track-input" type="number" min="1" step="1" inputmode="numeric" placeholder="Star number"><p class="control-help">Rings one star and shows its laps, speed and distance. Star 1 starts nearest the centre, star <span id="trackLast">800</span> furthest out. Clear the box to stop.</p></div>
    </section>
    <section class="panel-grid">
    <article class="panel"><div class="panel-heading"><h2><span class="fig">Fig. 1</span>3D galaxy simulation</h2><span>Test-star orbits over 1 billion years</span></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional __DEFAULT_GALAXY_NAME__ star trajectories"><div id="escape-counter" class="telemetry" aria-label="Live simulation readout"><div class="tm-row"><span>Ejected</span><b><span id="escape-count-val">0</span> / <span id="star-total-val">800</span></b></div><div class="tm-row" id="tmHaloRow"><span>Halo r<sub>200</sub></span><b id="tm-r200">—</b></div><div class="tm-row tm-chi" id="tmChiRow"><span>χ²ᵣ</span><b id="tm-chi">—</b></div></div><div class="stability-chip" id="stabilityChip" role="status" hidden></div><div class="view-panel" id="viewPanel" role="group" aria-label="3D view controls"><button type="button" class="view-head" id="viewToggle" aria-expanded="true" title="Collapse or expand the view controls"><span>View</span><span class="chev" aria-hidden="true">▾</span></button><div class="view-row view-primary"><button type="button" id="toggleFullscreen" aria-pressed="false" title="Toggle fullscreen view of 3D simulation">⛶ Fullscreen</button><button type="button" id="viewReset" title="Return to the default __DEFAULT_GALAXY_NAME__ view">↺ Reset</button></div><div class="view-label">Camera</div><div class="view-seg"><button type="button" id="viewFaceOn" title="Look straight down on the galactic plane">Face-on</button><button type="button" id="viewEdgeOn" title="Look along the galactic plane">Edge-on</button><button type="button" id="viewFullHalo" title="Pull back to see the full dark matter halo around the tiny disk">Full halo</button></div><div class="view-label">Show</div><div class="view-switches"><button type="button" role="switch" id="toggleOrbits" aria-checked="false" aria-pressed="false" title="Show or hide the orbit trajectory lines"><span class="sw" aria-hidden="true"></span>Orbits</button><button type="button" role="switch" id="toggleGas" aria-checked="true" aria-pressed="true" title="Show or hide the faint cold-gas cloud around the stars"><span class="sw" aria-hidden="true"></span>Gas</button><button type="button" role="switch" id="toggleHalo" aria-checked="true" aria-pressed="true" title="Show or hide the dark matter halo"><span class="sw" aria-hidden="true"></span>Halo</button><button type="button" role="switch" id="toggleEjection" class="sw-eject" aria-checked="false" aria-pressed="false" title="Show the red ejection line: stars beyond 1.5 × R_max count as ejected"><span class="sw" aria-hidden="true"></span>Ejection line</button><button type="button" role="switch" id="toggleMondZone" class="sw-mond" aria-checked="false" aria-pressed="false" hidden title="Show where MOND boosts gravity: the gold deepens where the boost is bigger, as on the graph (MOND only)"><span class="sw" aria-hidden="true"></span>Where MOND acts</button></div><div class="view-label">Speed</div><div class="speed-toggle" role="radiogroup" aria-label="Playback speed"><button type="button" role="radio" aria-checked="false" data-speed="0.25" title="0.25x speed">0.25x</button><button type="button" role="radio" aria-checked="false" data-speed="0.5" title="0.5x speed">0.5x</button><button type="button" role="radio" aria-checked="true" data-speed="1" title="1x speed">1x</button><button type="button" role="radio" aria-checked="false" data-speed="2" title="2x speed">2x</button></div><div class="mode-toggle" role="radiogroup" aria-label="Mouse drag mode" title="What a mouse drag does"><button type="button" role="radio" aria-checked="true" data-mode="rotate">Rotate</button><button type="button" role="radio" aria-checked="false" data-mode="pan">Pan</button></div></div><div class="fs-sliders" id="fsSliders" role="group" aria-label="Quick controls"><div class="fs-time"><button type="button" id="fsPlay" class="fs-play" aria-pressed="false" aria-label="Play simulation" title="Play / pause (Space)">▶</button><button type="button" id="fsRestart" class="fs-play" aria-label="Restart from 0 Myr" title="Back to 0 Myr and play">↺</button><input id="fsTrack" class="track-input" type="number" min="1" step="1" inputmode="numeric" placeholder="Star #" aria-label="Track a star by number" title="Type a star number to ring it and count its laps"><output id="fsTimeValue">0 / 1000 Myr</output></div><div class="fs-row" id="fsSandboxRows" hidden></div><div class="fs-row" id="fsHaloRows"><label for="fsMassSlider"><span>Halo mass</span><output id="fsMassValue"></output></label><input id="fsMassSlider" type="range" min="0.00" max="2.00" step="0.01"><label for="fsScaleSlider"><span>Scale radius</span><output id="fsScaleValue"></output></label><input id="fsScaleSlider" type="range" min="5" max="60" step="0.5"></div></div></div>
      <section class="halo-info" id="haloInfo" aria-labelledby="haloInfoTitle">
        <h3 id="haloInfoTitle">Your NFW halo</h3>
        <p class="halo-info-intro">The halo's size and shape update with the sliders. 1 kpc is about 3,260 light-years. <span class="on-mouse">Hover over a value to learn more.</span><span class="on-touch">Tap a value to learn more.</span></p>
        <dl class="halo-metrics" id="haloMetrics">
          <div class="halo-metric" id="haloRadiusCard"><dt><button type="button" aria-expanded="false" aria-controls="haloRadiusHelp">Halo radius, r<sub>200</sub><span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloRadius">—</span></button></dt><dd id="haloRadiusHelp"><span><span>Where the halo's average density falls to 200 times the Universe's critical density. A convention, not a sharp edge.</span></span></dd></div>
          <div class="halo-metric" id="haloScaleCard"><dt><button type="button" aria-expanded="false" aria-controls="haloScaleHelp">Scale radius, r<sub>s</sub><span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloScale">—</span></button></dt><dd id="haloScaleHelp"><span><span>Sets how the halo's density changes with distance from the centre.</span></span></dd></div>
          <div class="halo-metric" id="haloConcentrationCard"><dt><button type="button" aria-expanded="false" aria-controls="haloConcentrationHelp">Concentration, c<span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloConcentration">—</span></button></dt><dd id="haloConcentrationHelp"><span><span>Halo radius ÷ scale radius. Higher c means the dark matter is packed more towards the centre.</span></span></dd></div>
        </dl>
        <p class="halo-empty" id="haloEmpty" hidden></p>
      </section>
    </article>
      <article class="panel" id="curvePanel">
        <div class="panel-heading"><h2><span class="fig">Fig. 2</span>Rotation curve</h2><div class="chi-wrap"><small id="chiFit" class="chi-fit" title="Every model uses the same settings: SPARC’s standard values, or, with Best realistic fits on, values each model tunes within the same limits."></small><span id="chiBadge" title="Goodness of fit of visible matter + NFW halo to the SPARC points"><b id="chiVerdict"></b> <span id="chiNumbers">χ² –</span></span></div></div>
        <div class="curve-wrap"><div id="curvePlot" class="plot" role="img" aria-label="__DEFAULT_GALAXY_NAME__ observed, baryonic, and fitted rotation curves"></div></div>
        <p class="curve-hint">Tap the graph, or slide your finger along it, to read the speeds at each radius. Tap outside the graph to close.</p>
        <div class="stability" id="stabilityMonitor" data-state="stable" aria-label="Halo shape monitor" hidden>
          <div class="stability-head"><span class="stability-title">Halo shape monitor</span><span class="stability-badge" id="stabilityBadge" role="status"></span></div>
          <div class="stability-reading"><span class="stability-c" id="stabilityC">c = —</span><span class="stability-detail" id="stabilityDetail"></span></div>
          <div class="c-gauge" aria-hidden="true"><div class="c-band" id="cBand"></div><div class="c-marker" id="cMarker"></div></div>
          <div class="c-ticks" id="cTicks" aria-hidden="true"></div>
          <p class="stability-message" id="stabilityMessage"></p>
        </div>
        <details class="key">
          <summary>Understanding the colours<span data-preset-only> and χ²</span></summary>
          <div class="key-content">
            <ul><li data-preset-only><i class="sw dot" style="background:#f4f7fb"></i><b>White points</b> are measured speeds of <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> (SPARC); bars show uncertainty.</li><li><i class="sw ldash" style="border-color:#5fd0e0"></i><b>Cyan dashed</b> is the speed expected from visible matter only (stars + gas).</li><li><i class="sw dash" style="border-color:#b39ae8"></i><b>Purple</b> is the extra speed supplied by the dark matter halo (NFW mode).</li><li><i class="sw" style="background:#ffae5c"></i><b>Orange</b> is the combined model (or the MOND prediction in MOND mode).<span data-preset-only> The closer it hugs the white points, the better the fit.</span></li><li data-preset-only><i class="sw bar" style="background:#5fd0e0"></i><i class="sw bar" style="background:#ff8a80"></i><b>Residual bars:</b> grey bars show the gap at each point: up = model too slow, down = model too fast, no bar = perfect match.</li><li data-preset-only><b>Fit badge:</b> <span class="g">green</span> close fit (χ<sup>2</sup><sub>r</sub> below 2), <span class="o">yellow</span> rough fit (2 to 10), <span class="r">red</span> poor fit (above 10). These cut-offs are this site’s own rule of thumb, not a published standard.</li><li id="componentKey" hidden><b>Thin dotted lines</b> show your custom galaxy's stellar disk, bulge and gas contributions; their squared speeds add to give the cyan curve.</li></ul>
            <div data-preset-only>
            <h3>What do χ² and χ<sup>2</sup><sub>r</sub> mean?</h3>
            <p>For each point, take the gap between model and measurement, divide by that point’s error bar and square it. χ² is the total. χ<sup>2</sup><sub>r</sub> (reduced chi-squared) divides it by the number of points minus the values tuned to this galaxy, so galaxies and models compare fairly. Missing every point by one error bar gives about 1; by three, about 9.</p>
            <h3>Why not 0?</h3>
            <p>Measurements are noisy, so even a correct model misses by about one error bar. Both models also simplify real galaxies (a smooth round halo, an estimated star mass).</p>
            <h3>Why count tuned values?</h3>
            <p>More adjustable numbers can always bend a model closer to the data. With standard values, the halo tunes 2 (its mass and size) and MOND none: a<sub>0</sub> should be the same everywhere. Best realistic fits also tune the galaxy’s M/L, distance and tilt (and MOND’s a<sub>0</sub>): 4 to 6 in total.</p>
            </div>
          </div>
        </details>
      </article>
    </section>
    </div>
    <nav class="page-next" aria-label="Next page"><a class="next-link" href="#scorecard">Next: how both ideas score on every galaxy →</a></nav>
    </div>
    <div class="page" data-page="scorecard">
    <div class="section-head" id="scorecard"><span class="kicker">04</span><h2>Scorecard</h2><p>How closely each idea matches every galaxy on this site.</p></div>
    <section class="scorecard" aria-label="Fit scores for every preset galaxy">
      <p>Each fit gets a score, χ<sup>2</sup><sub>r</sub> (reduced chi-squared), which compares how far the model’s curve misses the measured points with the size of their error bars. Lower is better. About 1 means the misses are no bigger than the measurements’ own uncertainty: as close as the data allow.</p>
      <ul class="score-legend" aria-label="What the colours mean"><li><span class="score good">below 2</span>close fit</li><li><span class="score ok">2 to 10</span>rough fit</li><li><span class="score poor">above 10</span>poor fit</li></ul>
      <p class="score-notes-title">At a glance</p>
      <div class="score-scroll" tabindex="0" role="region" aria-label="Scorecard summary">
        <table class="score-table glance-table" id="glanceTable">
          <thead><tr><th scope="col">Galaxy</th><th scope="col">Visible matter only<small>Newton’s gravity, nothing added</small></th><th scope="col">Dark matter halo (NFW)<small>standard values → best realistic fit</small></th><th scope="col">MOND<small>standard values → best realistic fit</small></th></tr></thead>
          <tbody id="glanceBody"><tr><td colspan="4">Calculating the scores…</td></tr></tbody>
        </table>
      </div>
      <p class="glance-note"><sup>†</sup> The best fit needed at least one value at its realistic limit (see the full scores).</p>
      <details class="score-details" id="scoreDetails">
      <summary>Full scores, error bars and tuned values</summary>
      <div class="score-scroll" tabindex="0" role="region" aria-label="Scorecard table">
        <table class="score-table">
          <thead><tr><th scope="col" rowspan="2">Galaxy</th><th scope="colgroup" colspan="3" class="score-group">Standard values<small>SPARC’s M/L, distance and tilt for every model</small></th><th scope="colgroup" colspan="2" class="score-group">Best realistic fits<small>each model also tunes M/L, distance and tilt</small></th></tr><tr><th scope="col">Visible matter only<small>Newton’s gravity, nothing added</small></th><th scope="col">Dark matter halo (NFW)<small>halo mass and size tuned</small></th><th scope="col">MOND<small>nothing tuned</small></th><th scope="col">Dark matter halo (NFW)</th><th scope="col">MOND<small>a₀ also tuned</small></th></tr></thead>
          <tbody id="scoreBody"><tr><td colspan="6">Calculating the scores…</td></tr></tbody>
        </table>
      </div>
      </details>
      <p class="score-notes-title">Before drawing conclusions</p>
      <ul class="score-notes">
        <li><b>A fair test</b>With standard values, all three models start from the same visible matter: SPARC’s mass-to-light ratio (M/L, how much mass the stars have for each unit of their light), distance and tilt. Only the explanation for the extra pull differs. Best realistic fits let both models adjust those same values. The score allows for how many values each model tuned (listed in the full scores), so extra freedom does not come for free.</li>
        <li><b>Only realistic values</b>Best fits stay inside believable ranges: M/L values that models of stars allow (McGaugh &amp; Schombert 2014), distance and tilt within twice SPARC’s quoted errors, and a halo about as dense as is typical for its mass (Dutton &amp; Macciò 2014). A value marked <span class="at-limit">at limit</span> hit the edge of its range: the model could only have fitted that galaxy better with an unrealistic value.</li>
        <li><b>MOND’s constant</b>In best fits, MOND may also adjust a<sub>0</sub>, from __A0_MIN__ to __A0_MAX__ × 10⁻¹⁰ m/s² (within __FIT_LIMIT_SIGMA__× its measured error, like every other limit; McGaugh, Lelli &amp; Schombert 2016). Their quoted ±0.24 is wider because it includes the stars’ M/L, which the fits already vary. MOND says a<sub>0</sub> is the same in every galaxy (Milgrom 1983), so if galaxies need clearly different values, that is a warning sign even when the scores look good.</li>
        <li><b>Why MOND scores badly here</b>Plain MOND tunes nothing to each galaxy, and some of SPARC’s error bars are tiny (under 1 km/s across much of UGC 128), so a small but steady miss adds up to a large χ². The best fits here also use hard edges for M/L, distance, tilt and a<sub>0</sub>. Li et al. (2018) fitted the same formula to all 175 SPARC galaxies with soft limits, which let a value stray further at a cost, and found acceptable fits for most of them. A poor score here is a result for these galaxies under these limits, not a verdict on MOND.</li>
        <li><b>What the error bars cover</b>Scores use SPARC’s error bar on each measured speed. A wrong distance or tilt would shift the whole curve at once, not single points, so the best fits let distance and tilt vary instead.</li>
        <li><b>A small sample</b>These are <span id="scoreCount">a few</span> of SPARC’s 175 galaxies (Lelli, McGaugh &amp; Schombert 2016): enough to show how the method works, not to decide between the models.</li>
      </ul>
    </section>
    <nav class="page-next" aria-label="Next page"><a class="next-link" href="#background">Next: the two ideas and the wider evidence →</a></nav>
    </div>
    <div class="page" data-page="background">
    <div class="section-head" id="background"><span class="kicker">05</span><h2>Background</h2><p>The two explanations this simulator tests, and the wider evidence.</p></div>
    <section class="learn" id="theoryGuide" aria-label="Theory guide">
      <article class="nfw"><h3>What is a dark matter halo?</h3>
        <div class="learn-content">
        <p class="tldr">Dark matter is matter that gives out no light, so we only notice it through its gravity. A halo is a huge, roughly round cloud of it surrounding a galaxy and stretching far beyond its stars.</p>
        <p>The halo’s gravity supplies the missing pull, keeping the fast outer stars of <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> in orbit. Newton’s gravity stays exactly as it is; there is just more mass than we can see.</p>
        <p>This simulator uses an <b>NFW halo</b>, named after Navarro, Frenk &amp; White (1996). They did not invent its shape to fit rotation curves: it came out of computer simulations of how dark matter clumps together as the Universe grows. It is dense in the middle and thins out further away. Two sliders set it: its mass (M<sub>200</sub>) and its scale radius (r<sub>s</sub>, how spread out it is).</p>
        <figure class="schematic">
          <svg viewBox="0 0 320 150" role="img" aria-label="Schematic: the speed from visible matter falls with radius, the halo’s share grows, and together they stay level">
            <path d="M36 14V124H308" fill="none" stroke="currentColor" stroke-opacity=".35"/>
            <path d="M36.0 124.0L42.8 96.8L49.6 76.9L56.4 65.8L63.2 60.9L70.0 59.7L76.8 60.2L83.6 61.6L90.4 63.4L97.2 65.3L104.0 67.2L110.8 69.0L117.6 70.7L124.4 72.4L131.2 73.9L138.0 75.3L144.8 76.6L151.6 77.8L158.4 79.0L165.2 80.0L172.0 81.0L178.8 82.0L185.6 82.9L192.4 83.7L199.2 84.5L206.0 85.2L212.8 86.0L219.6 86.6L226.4 87.3L233.2 87.9L240.0 88.5L246.8 89.0L253.6 89.5L260.4 90.0L267.2 90.5L274.0 91.0L280.8 91.4L287.6 91.9L294.4 92.3L301.2 92.7L308.0 93.1" fill="none" stroke="#5fd0e0" stroke-width="2" stroke-dasharray="5 4"/>
            <path d="M36.0 124.0L42.8 101.0L49.6 92.3L56.4 86.1L63.2 81.3L70.0 77.4L76.8 74.1L83.6 71.3L90.4 68.9L97.2 66.7L104.0 64.9L110.8 63.2L117.6 61.7L124.4 60.4L131.2 59.2L138.0 58.1L144.8 57.1L151.6 56.2L158.4 55.4L165.2 54.7L172.0 54.1L178.8 53.4L185.6 52.9L192.4 52.4L199.2 51.9L206.0 51.5L212.8 51.1L219.6 50.8L226.4 50.4L233.2 50.1L240.0 49.9L246.8 49.6L253.6 49.4L260.4 49.2L267.2 49.0L274.0 48.8L280.8 48.7L287.6 48.6L294.4 48.4L301.2 48.3L308.0 48.2" fill="none" stroke="#b39ae8" stroke-width="2" stroke-dasharray="2 3"/>
            <path d="M36.0 124.0L42.8 88.3L49.6 67.2L56.4 54.6L63.2 47.9L70.0 44.5L76.8 43.0L83.6 42.3L90.4 42.1L97.2 42.0L308.0 42.0" fill="none" stroke="#ffae5c" stroke-width="2.4"/>
            <g fill="#f4f7fb"><circle cx="55" cy="57" r="3"/><circle cx="74" cy="43" r="3"/><circle cx="101" cy="42" r="3"/><circle cx="131" cy="42" r="3"/><circle cx="161" cy="42" r="3"/><circle cx="191" cy="42" r="3"/><circle cx="221" cy="42" r="3"/><circle cx="251" cy="42" r="3"/><circle cx="281" cy="42" r="3"/></g>
            <text x="300" y="30" text-anchor="end" fill="#ffae5c">Together: stays level</text>
            <text x="300" y="66" text-anchor="end" fill="#b39ae8">Halo: grows</text>
            <text x="300" y="110" text-anchor="end" fill="#5fd0e0">Visible matter: falls</text>
            <text x="172" y="142" text-anchor="middle">radius →</text>
            <text x="14" y="70" text-anchor="middle" transform="rotate(-90 14 70)">speed →</text>
          </svg>
          <figcaption>Schematic only, not real data. Colours match the simulator: the halo’s share of the speed (purple) grows as the visible matter’s (cyan) falls, so together (orange) they stay level.</figcaption>
        </figure>
        <details class="maths"><summary>Show the maths (A-level circular motion)</summary>
        <p><b>1. Orbits.</b> Gravity supplies the centripetal acceleration. For a round distribution of mass, only the mass M(r) inside the orbit pulls overall (Newton’s shell theorem):</p>
        <div class="eq"><math display="block"><mfrac><msup><mi>v</mi><mn>2</mn></msup><mi>r</mi></mfrac><mo>=</mo><mfrac><mrow><mi>G</mi><mi>M</mi><mo>(</mo><mi>r</mi><mo>)</mo></mrow><msup><mi>r</mi><mn>2</mn></msup></mfrac><mspace width="1em"/><mo>⇒</mo><mspace width="1em"/><mi>v</mi><mo>=</mo><msqrt><mfrac><mrow><mi>G</mi><mi>M</mi><mo>(</mo><mi>r</mi><mo>)</mo></mrow><mi>r</mi></mfrac></msqrt></math></div>
        <p><b>2. Visible matter alone.</b> Beyond most of the stars, M(r) hardly grows, so v ∝ 1/√r: the speed should fall, like the planets’ speeds in the Solar System.</p>
        <p><b>3. What a flat curve needs.</b> For v to stay constant, M(r) must grow in step with r. That happens if the density falls as 1/r²: each thin shell of thickness Δr then holds the same mass, 4πr²ρΔr, so every step outwards adds the same amount.</p>
        <p><b>4. The NFW density.</b> With ρ<sub>s</sub> setting how dense the halo is:</p>
        <div class="eq"><math display="block"><mi>ρ</mi><mo>(</mo><mi>r</mi><mo>)</mo><mo>=</mo><mfrac><msub><mi>ρ</mi><mi>s</mi></msub><mrow><mfrac><mi>r</mi><msub><mi>r</mi><mi>s</mi></msub></mfrac><msup><mrow><mo>(</mo><mn>1</mn><mo>+</mo><mfrac><mi>r</mi><msub><mi>r</mi><mi>s</mi></msub></mfrac><mo>)</mo></mrow><mn>2</mn></msup></mrow></mfrac></math></div>
        <p>Well inside r<sub>s</sub> the density falls as 1/r; far outside, as 1/r³. Around r<sub>s</sub> it falls roughly as 1/r², just what step 3 needs, so over that range the halo’s speed stays nearly level. It peaks at about 2.2 r<sub>s</sub>, then falls only slowly (Navarro, Frenk &amp; White 1996). Adding up the shells gives the mass formula used in the simulator (see Method).</p>
        </details>
        <p><b>A known weakness:</b> in many small, faint galaxies the speeds near the centre rise too gently for NFW’s dense middle. The data suggest a flatter centre instead, the “cusp–core problem” (de Blok 2010). Whether gas blown around by exploding stars could stir the dark matter enough to flatten it is still debated.</p>
        </div>
      </article>
      <article class="mond"><h3>What is MOND?</h3>
        <div class="learn-content">
        <p class="tldr">MOND (Modified Newtonian Dynamics) says there is no dark matter. Instead, gravity itself behaves differently when it is extremely weak.</p>
        <p>Mordehai Milgrom proposed it in 1983. MOND boosts gravity everywhere, but by how much depends on how strong gravity already is. Near a galaxy’s centre, where gravity is strong, the boost is small. The weaker gravity gets, the bigger the boost, so the outer stars keep moving fast.</p>
        <p><b>What a<sub>0</sub> is.</b> MOND’s constant, a<sub>0</sub> = 1.2×10⁻¹⁰ m/s² (about 10¹¹ times weaker than gravity at Earth’s surface), is a turning point, not a switch. Well above it the boost dies away fast; below it the boost keeps climbing. With the simulator’s formula, a star’s speed is boosted by about 2% where gravity is 10 × a<sub>0</sub>, 10% at 3 × a<sub>0</sub>, 26% at a<sub>0</sub>, 92% at a<sub>0</sub> ÷ 10 and 224% at a<sub>0</sub> ÷ 100. Even where gravity is strongest in the four galaxies here, MOND adds at least 4%.</p>
        <p><b>Where a<sub>0</sub> comes from.</b> It is not worked out from theory: it is measured from galaxies. Milgrom set it to match the galaxies known in 1983, and today it is fitted to SPARC’s rotation curves (McGaugh, Lelli &amp; Schombert 2016). Their study found that how much extra gravity a galaxy needs depends on its acceleration, whatever its size, which is why one number can work for every galaxy.</p>
        <p>There is no halo to adjust: once the visible matter is known, MOND predicts the whole rotation curve with one constant, a<sub>0</sub>, shared by every galaxy.</p>
        <figure class="schematic">
          <svg viewBox="0 0 320 150" role="img" aria-label="Schematic: Newtonian gravity from visible matter gives falling speeds, but MOND levels off where the acceleration drops below a0">
            <path d="M36 14V124H308" fill="none" stroke="currentColor" stroke-opacity=".35"/>
            <path d="M118 18V124" fill="none" stroke="currentColor" stroke-opacity=".45" stroke-dasharray="2 3"/>
            <path d="M36 124C58 52 78 38 108 42S196 78 302 98" fill="none" stroke="#5fd0e0" stroke-width="2" stroke-dasharray="5 4"/>
            <path d="M36 124C58 52 78 38 108 42S170 45 302 43" fill="none" stroke="#ffae5c" stroke-width="2.4"/>
            <g fill="#f4f7fb"><circle cx="58" cy="70" r="3"/><circle cx="84" cy="49" r="3"/><circle cx="110" cy="42" r="3"/><circle cx="140" cy="40" r="3"/><circle cx="170" cy="39" r="3"/><circle cx="200" cy="40" r="3"/><circle cx="230" cy="38" r="3"/><circle cx="260" cy="40" r="3"/><circle cx="290" cy="39" r="3"/></g>
            <text x="122" y="24">g ≈ a₀</text>
            <text x="300" y="30" text-anchor="end" fill="#ffae5c">MOND: levels off</text>
            <text x="300" y="114" text-anchor="end" fill="#5fd0e0">Newton, visible matter: falls</text>
            <text x="172" y="142" text-anchor="middle">radius →</text>
            <text x="14" y="70" text-anchor="middle" transform="rotate(-90 14 70)">speed →</text>
          </svg>
          <figcaption>Schematic only, not real data. The dotted line marks where gravity falls to a₀. MOND already boosts a little inside it; beyond it the boost grows large and the speed levels off instead of falling (Milgrom 1983).</figcaption>
        </figure>
        <details class="maths"><summary>Show the maths (A-level circular motion)</summary>
        <p><b>1. Orbits.</b> Gravity supplies the centripetal acceleration, and Newton gives the pull g<sub>N</sub> of the visible mass M:</p>
        <div class="eq"><math display="block"><mi>g</mi><mo>=</mo><mfrac><msup><mi>v</mi><mn>2</mn></msup><mi>r</mi></mfrac><mo>,</mo><mspace width="1.5em"/><msub><mi>g</mi><mi>N</mi></msub><mo>=</mo><mfrac><mrow><mi>G</mi><mi>M</mi></mrow><msup><mi>r</mi><mn>2</mn></msup></mfrac></math></div>
        <p><b>2. MOND’s rule.</b> Strong gravity is almost unchanged; weak gravity is boosted:</p>
        <div class="eq"><math display="block"><mi>g</mi><mo>=</mo><msub><mi>g</mi><mi>N</mi></msub><mtext>&#8195;if&#8194;</mtext><msub><mi>g</mi><mi>N</mi></msub><mo>≫</mo><msub><mi>a</mi><mn>0</mn></msub></math></div>
        <div class="eq"><math display="block"><mi>g</mi><mo>=</mo><msqrt><msub><mi>g</mi><mi>N</mi></msub><msub><mi>a</mi><mn>0</mn></msub></msqrt><mtext>&#8195;if&#8194;</mtext><msub><mi>g</mi><mi>N</mi></msub><mo>≪</mo><msub><mi>a</mi><mn>0</mn></msub></math></div>
        <p><b>3. Far out,</b> put the weak rule into step 1:</p>
        <div class="eq"><math display="block"><mfrac><msup><mi>v</mi><mn>2</mn></msup><mi>r</mi></mfrac><mo>=</mo><mfrac><msqrt><mi>G</mi><mi>M</mi><msub><mi>a</mi><mn>0</mn></msub></msqrt><mi>r</mi></mfrac></math></div>
        <div class="eq"><math display="block"><mo>⇒</mo><mspace width="0.5em"/><msup><mi>v</mi><mn>4</mn></msup><mo>=</mo><mi>G</mi><mi>M</mi><msub><mi>a</mi><mn>0</mn></msub></math></div>
        <p>r cancels, so the speed stops falling: a flat rotation curve set by the visible mass alone (Milgrom 1983). Between the two limits the simulator blends smoothly, using a formula fitted to SPARC galaxies (McGaugh, Lelli &amp; Schombert 2016; see Method).</p>
        </details>
        <p>Compare both on the Simulator page; one galaxy alone cannot settle which is right.</p>
        </div>
      </article>
      <article class="evidence"><h3>Beyond rotation curves: the wider evidence</h3>
        <div class="learn-content">
        <p class="tldr">Rotation curves are one clue among several. A winning explanation has to account for all of them.</p>
        <div class="evidence-grid">
          <div class="evidence-card"><span class="evidence-tag dm">Points to dark matter</span><h4>Galaxy clusters</h4><p><b>Seen:</b> galaxies in the Coma cluster move too fast for the cluster’s visible mass to hold them together (Zwicky 1933).</p><p><b>Why it matters:</b> the first hint of missing mass. MOND shrinks the gap but still needs some unseen matter (Famaey &amp; McGaugh 2012).</p></div>
          <div class="evidence-card"><span class="evidence-tag dm">Points to dark matter</span><h4>The Bullet Cluster</h4><p><b>Seen:</b> two clusters collided. Their hot gas, most of the visible matter, slowed down in the middle, but the bending of background light shows most of the mass carried on with the galaxies (Clowe et al. 2006).</p><p><b>Why it matters:</b> the mass and the visible matter ended up in different places.</p></div>
          <div class="evidence-card"><span class="evidence-tag dm">Points to dark matter</span><h4>The oldest light</h4><p><b>Seen:</b> the cosmic microwave background, light left over from the early Universe, has a pattern of slightly hotter and cooler patches.</p><p><b>Why it matters:</b> the pattern fits a universe with about five times as much dark matter as ordinary matter (Planck Collaboration 2020).</p></div>
          <div class="evidence-card"><span class="evidence-tag mond">Points to MOND</span><h4>Visible matter predicts the pull</h4><p><b>Seen:</b> across 153 galaxies, the pull worked out from the measured speeds closely follows one curve set only by the visible matter (McGaugh, Lelli &amp; Schombert 2016).</p><p><b>Why it matters:</b> MOND predicts exactly this. Dark matter can match it only if every halo is closely tied to its galaxy’s visible matter, which still needs explaining.</p></div>
        </div>
        <p class="evidence-open-title">Still open</p>
        <div class="evidence-grid two">
          <div class="evidence-card"><h4>No confirmed dark matter particle</h4><p>Decades of searches have found none (Bertone &amp; Hooper 2018). In September 2026 the LZ experiment reported one unexplained event that might be one, but at 2.6σ it is far short of the 5σ usually needed to claim a discovery, and it is not yet peer reviewed (LZ Collaboration 2026).</p></div>
          <div class="evidence-card"><h4>MOND beyond galaxies</h4><p>MOND has not matched galaxy clusters or the cosmic microwave background without adding unseen matter of its own (Famaey &amp; McGaugh 2012).</p></div>
        </div>
        </div>
      </article>
    </section>
    <nav class="page-next" aria-label="Next page"><a class="next-link" href="#method">Next: method and sources →</a></nav>
    </div>
    <div class="page" data-page="method">
    <div class="section-head" id="method"><span class="kicker">06</span><h2>Method &amp; sources</h2><p>What the simulator calculates, and where the numbers and ideas come from.</p></div>
    <footer>
      <dl class="method-grid">
        <div data-preset-only><dt>Data</dt><dd><p>SPARC (Lelli, McGaugh &amp; Schombert 2016) gives, at each radius r of each galaxy:</p><ul><li><b>V<sub>obs</sub> ± σ</b>: the measured orbital speed and its uncertainty (error bar).</li><li><b>V<sub>gas</sub>, V<sub>disk</sub>, V<sub>bulge</sub></b>: the speed the gas, the stellar disk or the bulge would give on its own, worked out from infrared images of the stars and radio maps of the hydrogen gas.</li></ul><p>Read from the project’s copy of SPARC’s table2.dat.</p></dd></div>
        <div data-preset-only><dt>Tilt (inclination)</dt><dd><p>We see each disk at a tilt i: 0° face-on, 90° edge-on. A telescope only measures motion along our line of sight (the Doppler shift), so it sees just part of each star’s orbital speed:</p><span class="eq"><math display="block"><msub><mi>V</mi><mi>obs</mi></msub><mo>=</mo><mfrac><msub><mi>V</mi><mi>line of sight</mi></msub><mrow><mi>sin</mi><mi>i</mi></mrow></mfrac></math></span><p><b>In words:</b> divide the speed we see by sin i to get the true speed. SPARC has already done this using its own tilt, i<sub>SPARC</sub>, which has an error of a few degrees.</p><p><b>If the true tilt is different</b>, every point moves by the same factor (Li et al. 2018, eq. 5):</p><span class="eq"><math display="block"><mi>s</mi><mo>=</mo><mfrac><mrow><mi>sin</mi><msub><mi>i</mi><mi>SPARC</mi></msub></mrow><mrow><mi>sin</mi><mi>i</mi></mrow></mfrac></math></span><p><b>Example:</b> SPARC says 60° but the truth is 55°: s = 0.866 / 0.819 ≈ 1.06, so the real speeds are 6% higher. At 30° instead of 25°, s ≈ 1.18: nearly face-on galaxies are far more sensitive.</p><p><b>On the page,</b> Best realistic fits may change i within __FIT_LIMIT_SIGMA__× SPARC’s error (at most 90°). The points stay where SPARC put them and the model is divided by s instead, which gives the same χ². Tilt and distance can partly make up for each other, so the fitted tilt is less certain than the score.</p></dd></div>
        <div><dt>Visible matter</dt><dd><div id="footerModel"><span class="eq"><math display="block"><msub><mi>V</mi><mi>bar</mi></msub><mo>=</mo><msqrt><msub><mi>V</mi><mi>gas</mi></msub><mo stretchy="false" lspace="0" rspace="0">|</mo><msub><mi>V</mi><mi>gas</mi></msub><mo stretchy="false" lspace="0" rspace="0">|</mo><mo>+</mo><mn id="ftDisk">0.50</mn><msubsup><mi>V</mi><mi>disk</mi><mn>2</mn></msubsup><mo>+</mo><mn id="ftBulge">0.70</mn><msubsup><mi>V</mi><mi>bulge</mi><mn>2</mn></msubsup></msqrt></math></span><p><b>In words:</b> the speed a star would orbit at if only the visible matter pulled on it.</p><p><b>Why add squares?</b> For a circular orbit, the pull g = V²/r, so V² = gr. Pulls from different masses add, so V² values add; then take the square root.</p><p><b>The numbers</b> are mass-to-light ratios (in use: <b id="ftMlSource">standard SPARC values</b>; the standard 0.50 and 0.70 come from McGaugh &amp; Schombert 2014). SPARC gives V<sub>disk</sub> and V<sub>bulge</sub> for a ratio of 1, and V² is proportional to mass, so scaling the mass scales V².</p><p><b>Why V<sub>gas</sub>|V<sub>gas</sub>|?</b> Much of a galaxy’s gas lies in a ring further out. The shell theorem (only the mass inside you pulls) holds for a round cloud, not a flat disk: inside a flat ring, the nearer side of the ring pulls harder than the far side, so its overall pull points outwards and SPARC lists V<sub>gas</sub> as negative. V<sub>gas</sub>|V<sub>gas</sub>| keeps that minus sign (squaring would lose it), so the gas then lowers the total.</p></div><div id="footerCustom" hidden><span class="eq"><math display="block"><msub><mi>V</mi><mi>bar</mi></msub><mo>=</mo><msqrt><msubsup><mi>V</mi><mi>disk</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>bulge</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>gas</mi><mn>2</mn></msubsup></msqrt></math></span><p><b>In words:</b> the speed from your sandbox galaxy’s visible matter alone, added as squares for the same reason as in a real galaxy (V² = gr, and pulls add).</p><p>Each speed comes from the masses you set: an exponential disk (density fading steadily outwards) for the stars and the gas, and a Hernquist bulge (a standard ball-shaped profile).</p></div></dd></div>
        <div><dt>Dark matter (NFW halo)</dt><dd><span class="eq"><math display="block"><msub><mi>V</mi><mi>model</mi></msub><mo>=</mo><msqrt><msubsup><mi>V</mi><mi>bar</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>halo</mi><mn>2</mn></msubsup></msqrt></math></span><p><b>In words:</b> the model’s speed is visible matter plus the halo, added as squares like the parts above.</p><span class="eq"><math display="block"><msub><mi>V</mi><mi>halo</mi></msub><mo>=</mo><msqrt><mfrac><mrow><mi>G</mi><msub><mi>M</mi><mi>halo</mi></msub><mo>(</mo><mi>r</mi><mo>)</mo></mrow><mi>r</mi></mfrac></msqrt><mo>,</mo><mspace width="1em"/><msub><mi>M</mi><mi>halo</mi></msub><mo>(</mo><mi>r</mi><mo>)</mo><mo>=</mo><msub><mi>M</mi><mn>200</mn></msub><mfrac><mrow><mi>f</mi><mo>(</mo><mi>r</mi><mo>/</mo><msub><mi>r</mi><mi>s</mi></msub><mo>)</mo></mrow><mrow><mi>f</mi><mo>(</mo><mi>c</mi><mo>)</mo></mrow></mfrac></math></span><span class="eq"><math display="block"><mi>f</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mi>ln</mi><mo>(</mo><mn>1</mn><mo>+</mo><mi>x</mi><mo>)</mo><mo>−</mo><mfrac><mi>x</mi><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfrac><mo>,</mo><mspace width="1em"/><mi>c</mi><mo>=</mo><mfrac><msub><mi>r</mi><mn>200</mn></msub><msub><mi>r</mi><mi>s</mi></msub></mfrac></math></span><p><b>The halo’s speed</b> is the A-level orbit speed V = √(GM/r). For a round cloud, only the mass inside the orbit pulls overall (Newton’s shell theorem), so M<sub>halo</sub>(r) is the dark matter inside radius r.</p><p><b>f</b> describes how that mass builds up with radius in an NFW halo (Navarro, Frenk &amp; White 1996); dividing by f(c) makes the mass inside r<sub>200</sub> exactly M<sub>200</sub>. <b>r<sub>200</sub></b> is the halo’s conventional edge, where its average density falls to 200 times the Universe’s critical density; the model adds no mass beyond it. <b>r<sub>s</sub></b> (the scale radius slider) sets how concentrated the halo is.</p></dd></div>
        <div><dt>MOND</dt><dd><p>No halo (Milgrom 1983). Three steps, with a<sub>0</sub> = 1.2×10⁻¹⁰ m/s²:</p><p><b>1. Newton’s pull</b> from the visible matter:</p><span class="eq"><math display="block"><msub><mi>g</mi><mi>N</mi></msub><mo>=</mo><mfrac><msubsup><mi>V</mi><mi>bar</mi><mn>2</mn></msubsup><mi>r</mi></mfrac></math></span><p><b>2. MOND’s boost.</b> Divide by a switch-on factor:</p><span class="eq"><math display="block"><mi>g</mi><mo>=</mo><mfrac><msub><mi>g</mi><mi>N</mi></msub><mrow><mn>1</mn><mo>−</mo><msup><mi>e</mi><mrow><mo>−</mo><msqrt><msub><mi>g</mi><mi>N</mi></msub><mo>/</mo><msub><mi>a</mi><mn>0</mn></msub></msqrt></mrow></msup></mrow></mfrac></math></span><p>When g<sub>N</sub> is much bigger than a<sub>0</sub>, e<sup>−…</sup> is almost 0, so g = g<sub>N</sub>: plain Newton. When g<sub>N</sub> is much smaller, 1 − e<sup>−x</sup> ≈ x, which gives g ≈ √(g<sub>N</sub>a<sub>0</sub>), the stronger pull that flattens the curve. This is the radial acceleration relation’s formula (McGaugh, Lelli &amp; Schombert 2016, eq. 4).</p><p><b>Why this formula?</b> MOND fixes the two limits but not how gravity switches between them. Older work often used the “standard” formula, μ(x) = x/√(1+x²) (Famaey &amp; McGaugh 2012). This site uses the RAR formula because a<sub>0</sub> was measured with it, from the same SPARC data, and the realistic best fits follow Li et al. (2018), who use it too. So MOND’s formula has already been fitted to SPARC data (one curve shared by 153 galaxies, not tuned to each), while the NFW halo is fitted to each galaxy on its own: neither model meets this data completely fresh. The choice of formula matters for single galaxies. χ<sup>2</sup><sub>r</sub> for MOND at standard values:</p><table class="mini-table"><thead><tr><th>Galaxy</th><th>Standard</th><th>RAR (used)</th></tr></thead><tbody>__MOND_FORMULA_ROWS__</tbody></table><p><b>3. Back to a speed</b>, since g = V²/r:</p><span class="eq"><math display="block"><msub><mi>V</mi><mi>model</mi></msub><mo>=</mo><msqrt><mi>g</mi><mi>r</mi></msqrt></math></span></dd></div>
        <div><dt>Simulation</dt><dd><p>The 3D view follows massless <b>test stars</b> (800 for NGC 3198, scaled with each galaxy’s visible mass from 200 to 2400): they feel the model’s gravity but do not pull on each other, so each moves on its own. Disk stars start on circular orbits, bulge stars with random motions. Each star’s acceleration is a = V<sub>model</sub>(r)²/r, towards the centre.</p><p><b>Simplifications.</b> The pull is worked out as if the galaxy’s mass were a round ball, although the disk is flat, so a star’s up-and-down bobbing through the disk is only approximate. Very near the centre the pull is softened over 0.5–1 kpc (the spline softening of the GADGET-2 code; Springel 2005), so a star passing close to the centre is not flung out by an unrealistically huge pull.</p><p>An orbit cannot be solved in one go, so the computer moves each star in small time steps, Δt = 0.5 Myr, 2000 times over 1 Gyr. Each step uses the <b>Velocity-Verlet method</b>:</p><span class="eq"><math display="block"><msub><mi>x</mi><mi>new</mi></msub><mo>=</mo><mi>x</mi><mo>+</mo><mi>v</mi><mi>Δ</mi><mi>t</mi><mo>+</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mi>a</mi><msup><mrow><mi>Δ</mi><mi>t</mi></mrow><mn>2</mn></msup></math></span><p>This is SUVAT’s s = ut + ½at². Next, it finds the new acceleration a<sub>new</sub> at the new position, then updates the velocity with the <i>average</i> of the old and new accelerations:</p><span class="eq"><math display="block"><msub><mi>v</mi><mi>new</mi></msub><mo>=</mo><mi>v</mi><mo>+</mo><mfrac><mrow><mi>a</mi><mo>+</mo><msub><mi>a</mi><mi>new</mi></msub></mrow><mn>2</mn></mfrac><mi>Δ</mi><mi>t</mi></math></span><p><b>Why this method?</b> Using only the old acceleration would make each step slightly too big or small, and the error would build up, slowly spiralling stars in or out. Averaging keeps the energy error small and stops it building up over time, so a star that leaves does so because of the physics, not the size of the time step.</p></dd></div>
        <div data-preset-only><dt>Goodness of fit</dt><dd><span class="eq"><math display="block"><msup><mi>χ</mi><mn>2</mn></msup><mo>=</mo><mo>∑</mo><msup><mrow><mo>(</mo><mfrac><mrow><msub><mi>V</mi><mi>obs</mi></msub><mo>−</mo><msub><mi>V</mi><mi>model</mi></msub></mrow><mi>σ</mi></mfrac><mo>)</mo></mrow><mn>2</mn></msup></math></span><p><b>In words:</b> for each measured point, find how far the model misses, divide by that point’s error bar σ, square it, then add them all up. Squaring makes misses above and below both count, and big misses count much more.</p><span class="eq"><math display="block"><msubsup><mi>χ</mi><mi>r</mi><mn>2</mn></msubsup><mo>=</mo><mfrac><msup><mi>χ</mi><mn>2</mn></msup><mrow><mi>N</mi><mo>−</mo><mi>k</mi></mrow></mfrac></math></span><p><b>N</b> is the number of points and <b>k</b> the number of values tuned to fit this galaxy, so χ<sup>2</sup><sub>r</sub> is roughly the squared miss per point. About 1 means misses the size of the error bars.</p><p>With standard values: k = 0 for visible matter only and MOND, and k = 2 for the halo (M<sub>200</sub> and r<sub>s</sub>, also when set by hand). Best realistic fits also tune disk M/L (and bulge M/L), plus distance and tilt within __FIT_LIMIT_SIGMA__× SPARC’s errors: k = __K_NFW_FIT__ for NFW and __K_MOND_FIT__ for MOND (which adds a<sub>0</sub>).</p></dd></div>
      </dl>
      <section class="refs" id="references" aria-labelledby="refsTitle">
        <h3 id="refsTitle">References</h3>
        <ul>
          <li>Bertone, G. and Hooper, D. (2018) ‘History of dark matter’, <i>Reviews of Modern Physics</i>, 90(4), 045002. <a href="https://doi.org/10.1103/RevModPhys.90.045002" target="_blank" rel="noopener">doi:10.1103/RevModPhys.90.045002</a></li>
          <li>Clowe, D., Bradač, M., Gonzalez, A.H., Markevitch, M., Randall, S.W., Jones, C. and Zaritsky, D. (2006) ‘A direct empirical proof of the existence of dark matter’, <i>The Astrophysical Journal</i>, 648(2), pp. L109–L113. <a href="https://doi.org/10.1086/508162" target="_blank" rel="noopener">doi:10.1086/508162</a></li>
          <li>de Blok, W.J.G. (2010) ‘The core-cusp problem’, <i>Advances in Astronomy</i>, 2010, 789293. <a href="https://doi.org/10.1155/2010/789293" target="_blank" rel="noopener">doi:10.1155/2010/789293</a></li>
          <li>de Blok, W.J.G. and McGaugh, S.S. (1996) ‘Does low surface brightness mean low density?’, <i>The Astrophysical Journal</i>, 469(2), pp. L89–L92. <a href="https://doi.org/10.1086/310266" target="_blank" rel="noopener">doi:10.1086/310266</a></li>
          <li>Dutton, A.A. and Macciò, A.V. (2014) ‘Cold dark matter haloes in the Planck era: evolution of structural parameters for Einasto and NFW profiles’, <i>Monthly Notices of the Royal Astronomical Society</i>, 441(4), pp. 3359–3374. <a href="https://doi.org/10.1093/mnras/stu742" target="_blank" rel="noopener">doi:10.1093/mnras/stu742</a></li>
          <li>Famaey, B. and McGaugh, S.S. (2012) ‘Modified Newtonian Dynamics (MOND): observational phenomenology and relativistic extensions’, <i>Living Reviews in Relativity</i>, 15, 10. <a href="https://doi.org/10.12942/lrr-2012-10" target="_blank" rel="noopener">doi:10.12942/lrr-2012-10</a></li>
          <li>Freeman, K.C. (1970) ‘On the disks of spiral and S0 galaxies’, <i>The Astrophysical Journal</i>, 160, p. 811. <a href="https://doi.org/10.1086/150474" target="_blank" rel="noopener">doi:10.1086/150474</a></li>
          <li>Graham, A.W. and Worley, C.C. (2008) ‘Inclination- and dust-corrected galaxy parameters: bulge-to-disc ratios and size–luminosity relations’, <i>Monthly Notices of the Royal Astronomical Society</i>, 388(4), pp. 1708–1728. <a href="https://doi.org/10.1111/j.1365-2966.2008.13506.x" target="_blank" rel="noopener">doi:10.1111/j.1365-2966.2008.13506.x</a></li>
          <li>Hubble, E.P. (1926) ‘Extragalactic nebulae’, <i>The Astrophysical Journal</i>, 64, p. 321. <a href="https://doi.org/10.1086/143018" target="_blank" rel="noopener">doi:10.1086/143018</a></li>
          <li>Lelli, F., McGaugh, S.S. and Schombert, J.M. (2016) ‘SPARC: mass models for 175 disk galaxies with Spitzer photometry and accurate rotation curves’, <i>The Astronomical Journal</i>, 152(6), 157. <a href="https://doi.org/10.3847/0004-6256/152/6/157" target="_blank" rel="noopener">doi:10.3847/0004-6256/152/6/157</a></li>
          <li>Li, P., Lelli, F., McGaugh, S. and Schombert, J. (2018) ‘Fitting the radial acceleration relation to individual SPARC galaxies’, <i>Astronomy &amp; Astrophysics</i>, 615, A3. <a href="https://doi.org/10.1051/0004-6361/201732547" target="_blank" rel="noopener">doi:10.1051/0004-6361/201732547</a></li>
          <li>LZ Collaboration (Akerib, D.S. et al.) (2026) ‘Search for dark matter particle interactions in an extended nuclear recoil energy window with the LUX-ZEPLIN (LZ) experiment’. Preprint, not yet peer reviewed. <a href="https://arxiv.org/abs/2609.02823" target="_blank" rel="noopener">arXiv:2609.02823</a></li>
          <li>McGaugh, S.S., Lelli, F. and Schombert, J.M. (2016) ‘Radial acceleration relation in rotationally supported galaxies’, <i>Physical Review Letters</i>, 117(20), 201101. <a href="https://doi.org/10.1103/PhysRevLett.117.201101" target="_blank" rel="noopener">doi:10.1103/PhysRevLett.117.201101</a></li>
          <li>McGaugh, S.S. and Schombert, J.M. (2014) ‘Color–mass-to-light-ratio relations for disk galaxies’, <i>The Astronomical Journal</i>, 148(5), 77. <a href="https://doi.org/10.1088/0004-6256/148/5/77" target="_blank" rel="noopener">doi:10.1088/0004-6256/148/5/77</a></li>
          <li>Milgrom, M. (1983) ‘A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis’, <i>The Astrophysical Journal</i>, 270, p. 365. <a href="https://doi.org/10.1086/161130" target="_blank" rel="noopener">doi:10.1086/161130</a></li>
          <li>Moster, B.P., Naab, T. and White, S.D.M. (2013) ‘Galactic star formation and accretion histories from matching galaxies to dark matter haloes’, <i>Monthly Notices of the Royal Astronomical Society</i>, 428(4), pp. 3121–3138. <a href="https://doi.org/10.1093/mnras/sts261" target="_blank" rel="noopener">doi:10.1093/mnras/sts261</a></li>
          <li>Navarro, J.F., Frenk, C.S. and White, S.D.M. (1996) ‘The structure of cold dark matter halos’, <i>The Astrophysical Journal</i>, 462, p. 563. <a href="https://doi.org/10.1086/177173" target="_blank" rel="noopener">doi:10.1086/177173</a></li>
          <li>Planck Collaboration (2020) ‘Planck 2018 results. VI. Cosmological parameters’, <i>Astronomy &amp; Astrophysics</i>, 641, A6. <a href="https://doi.org/10.1051/0004-6361/201833910" target="_blank" rel="noopener">doi:10.1051/0004-6361/201833910</a></li>
          <li>Rubin, V.C., Ford, W.K. and Thonnard, N. (1980) ‘Rotational properties of 21 Sc galaxies with a large range of luminosities and radii, from NGC 4605 (R = 4 kpc) to UGC 2885 (R = 122 kpc)’, <i>The Astrophysical Journal</i>, 238, p. 471. <a href="https://doi.org/10.1086/158003" target="_blank" rel="noopener">doi:10.1086/158003</a></li>
          <li>Springel, V. (2005) ‘The cosmological simulation code GADGET-2’, <i>Monthly Notices of the Royal Astronomical Society</i>, 364(4), pp. 1105–1134. <a href="https://doi.org/10.1111/j.1365-2966.2005.09655.x" target="_blank" rel="noopener">doi:10.1111/j.1365-2966.2005.09655.x</a></li>
          <li>Williams, D.R. (2024) <i>Planetary fact sheet</i>. NASA Goddard Space Flight Center. <a href="https://nssdc.gsfc.nasa.gov/planetary/factsheet/" target="_blank" rel="noopener">nssdc.gsfc.nasa.gov/planetary/factsheet</a></li>
          <li>Zwicky, F. (1933) ‘Die Rotverschiebung von extragalaktischen Nebeln’, <i>Helvetica Physica Acta</i>, 6, pp. 110–127. English translation (2009): ‘The redshift of extragalactic nebulae’, <i>General Relativity and Gravitation</i>, 41(1), pp. 207–224. <a href="https://doi.org/10.1007/s10714-008-0707-4" target="_blank" rel="noopener">doi:10.1007/s10714-008-0707-4</a></li>
        </ul>
        <p>Each link opens the published source in a new tab.</p>
      </section>
      <button type="button" id="exportCsv" class="export-btn" title="Download the time-series log and the rotation-curve model as two CSV files">Export raw data (CSV)</button>
    </footer>
    </div>
  </main>
  <dialog class="tutorial-welcome" id="tutorialWelcome" aria-labelledby="welcomeTitle" aria-describedby="welcomeDescription">
    <div class="tutorial-eyebrow">Explore <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span></div>
    <h2 id="welcomeTitle">Welcome to the simulator</h2>
    <p id="welcomeDescription">New here? The simulator tutorial is a short, hands-on tour of this page. You will find the puzzle in a real galaxy’s rotation curve, then test two answers to it: dark matter and modified gravity. Skip it at any time. The Custom Sandbox, where you build your own galaxy, has a tutorial of its own.</p>
    <div class="tutorial-welcome-actions"><button type="button" class="tutorial-primary" id="tutorialStart">Start simulator tutorial</button><button type="button" id="tutorialWelcomeSkip">I've been here before</button></div>
  </dialog>
  <div class="tutorial-done" id="tutorialDone" role="status" aria-labelledby="tutorialDoneTitle" hidden>
    <div class="tutorial-done-box">
      <div class="tutorial-done-orbit" aria-hidden="true"><i></i><i></i><i></i><b>✓</b></div>
      <h2 id="tutorialDoneTitle">Tutorial complete!</h2>
      <p>You’re all set. Now it’s yours to play with: try other galaxies, switch between dark matter and MOND, or build a galaxy of your own.</p>
      <button type="button" class="tutorial-primary" id="tutorialDoneClose">Start exploring</button>
    </div>
  </div>
  <div class="tutorial-spotlight" id="tutorialSpotlight" aria-hidden="true" hidden>
    <span class="tutorial-tag" id="tutorialTag"></span>
    <svg class="tutorial-pointer" viewBox="0 0 28 34" aria-hidden="true"><path d="M3 2v23l6-5 4 10 5-2-4-10h8z" fill="#eaffff" stroke="#09212c" stroke-width="2" stroke-linejoin="round"/></svg>
  </div>
  <svg class="tutorial-dim" id="tutorialDim" aria-hidden="true"><mask id="tutorialDimMask"><rect width="100%" height="100%" fill="#fff"/><g id="tutorialDimHoles" fill="#000"></g></mask><rect width="100%" height="100%" mask="url(#tutorialDimMask)"/></svg>
  <aside class="tutorial-card" id="tutorialCard" aria-labelledby="tutorialTitle" hidden>
    <div class="tutorial-top"><span class="tutorial-eyebrow" id="tutorialProgress"></span><button type="button" class="tutorial-skip" id="tutorialSkip">Skip tutorial</button></div>
    <div class="tutorial-progress-track" aria-hidden="true"><div class="tutorial-progress-fill" id="tutorialProgressFill"></div></div>
    <div class="tutorial-content" id="tutorialContent" tabindex="0"><div aria-live="polite" aria-atomic="true"><h2 id="tutorialTitle" tabindex="-1"></h2><p id="tutorialText"></p><p class="tutorial-try" id="tutorialTry"></p></div></div>
    <div class="tutorial-sliders" id="tutorialSliders" hidden></div>
    <div class="tutorial-readout" id="tutorialReadout" aria-live="polite" hidden></div>
    <div class="tutorial-feedback" id="tutorialFeedback" hidden><span id="tutorialFeedbackText" role="status" aria-live="polite"></span><button type="button" class="tutorial-auto" id="tutorialAuto" aria-label="Go to the next slide now" title="Moving on in a few seconds. Press to go now." hidden>Next slide<span class="tutorial-clock" aria-hidden="true"><svg viewBox="0 0 20 20"><circle class="track" cx="10" cy="10" r="7"/><circle class="ring" id="tutorialAutoRing" cx="10" cy="10" r="7"/></svg><b id="tutorialAutoNum">5</b></span></button></div>
    <div class="tutorial-actions"><span class="tutorial-keys">← → to navigate</span><button type="button" id="tutorialBack">Back</button><button type="button" id="tutorialMore" hidden>Sandbox tutorial</button><button type="button" class="tutorial-primary" id="tutorialNext">Next</button></div>
  </aside>
  <script>
__GALAXIES_JS__
  </script>
  <script>
    /* Lite visuals (Settings' "High-quality visuals" unticked, or chosen automatically on a slow device): the html
       element gets the class "lite", which turns off the page's glows, pulsing highlights, blur and motion, and the 3D
       view's bloom. Set here, before anything draws, from the choice remembered in this browser. calmMotion() is the
       one test the page's scripts use for "skip the animation": lite visuals or the system's reduce-motion setting. */
    try{if(localStorage.getItem('galaxySimQuality')==='low')document.documentElement.classList.add('lite');}catch(e){}
    window.calmMotion=()=>document.documentElement.classList.contains('lite')||matchMedia('(prefers-reduced-motion: reduce)').matches;
    /* Runs before the simulator. If a library did not load (a network that blocks CDNs, or no internet) or the
       browser has WebGL turned off, say so plainly, show every page as one long scroll so the text still reads,
       and stop: the simulator below needs all of them. */
    (()=>{
      const canvas=document.createElement('canvas');
      let webgl=false;
      try{webgl=!!(window.WebGLRenderingContext&&(canvas.getContext('webgl2')||canvas.getContext('webgl')));}catch(err){}
      const missing=[];
      if(typeof Plotly==='undefined')missing.push('Plotly (graphs)');
      if(typeof THREE==='undefined'||!THREE.EffectComposer||!THREE.UnrealBloomPass)missing.push('Three.js (3D view)');
      if(!missing.length&&webgl)return;
      window.SIM_UNAVAILABLE=missing.length?`could not load ${missing.join(' and ')}`:'WebGL is turned off in this browser';
      const note=document.createElement('div');
      note.className='load-failed';note.setAttribute('role','alert');
      note.innerHTML=`<b>The interactive simulator could not start:</b> ${window.SIM_UNAVAILABLE}. `+(missing.length
        ?'These load from the internet, so a school or work network may be blocking them. Try another network, or reload.'
        :'Try another browser, or turn on hardware acceleration in this one.')+' The explanations below still work.';
      document.body.classList.add('sim-unavailable');
      document.querySelector('main').prepend(note);
    })();
  </script>
  <script>
    if(window.SIM_UNAVAILABLE)throw new Error('Simulator not started: '+window.SIM_UNAVAILABLE); // see the check above
    document.querySelectorAll('details').forEach(details=>{
      const summary=details.querySelector('summary');
      const content=summary.nextElementSibling;
      let animation=null, targetOpen=details.open;
      function updateDisclosure(){
        const opening=!targetOpen;
        targetOpen=opening;
        const height=details.open?content.getBoundingClientRect().height:0;
        const opacity=details.open?Number(getComputedStyle(content).opacity):0;
        if(animation){
          animation.cancel();
          animation=null;
          content.style.overflow='';
        }
        details.open=true;
        content.inert=!opening;
        if(calmMotion()){
          details.open=opening;
          content.inert=false;
          return;
        }
        content.style.overflow='hidden';
        animation=content.animate(
          [{height:height+'px',opacity},{height:(opening?content.scrollHeight:0)+'px',opacity:opening?1:0}],
          {duration:opening?350:200,easing:'ease',fill:'both'}
        );
        animation.onfinish=()=>{
          details.open=opening;
          content.inert=false;
          content.style.overflow='';
          animation.cancel();
          animation=null;
        };
      }
      summary.addEventListener('click',event=>{
        event.preventDefault();
        updateDisclosure();
      });
    });
    const DATA = __DATA_JSON__;
    // GALAXIES (the inlined galaxies.js above) holds one SPARC rotation curve per preset galaxy.
    // Best realistic fits switch: off, every model uses SPARC's standard values (the halo fit tunes only the halo);
    // on, the NFW and MOND best fits also tune M/L, distance and tilt within the same realistic limits.
    let bestFits = false;
    const nfwFitOf = galaxy => galaxy.best_fit[bestFits ? 'nfw' : 'nfw_standard'];
    const haloDefaultsOf = galaxy => ({massT: nfwFitOf(galaxy).m200_msun / 1e12, scaleRadius: nfwFitOf(galaxy).scale_radius_kpc});
    let activeGalaxy = GALAXIES[DATA.default_galaxy_id];
    let DEFAULTS = haloDefaultsOf(activeGalaxy);
    let customActive = false; // true while the Custom Sandbox, not a SPARC preset, drives the physics
    let customModel = null; // PhysicsEngine circular-speed functions of the sandbox galaxy
    let customGridRmax = null; // R_max of the sandbox's radius grid (null while a preset's SPARC data are loaded)
    const G = 4.30091e-6;
    const H0 = 0.07;
    const BG = '#050910';
    const GRID = '#1c2936';
    const textColor = '#b8c4d6';
    const SIM_TOTAL_TIME_MYR = 1000.0;
    const SIM_DT_MYR = 1.0;
    const frameCount = 1001; // 0 to 1000 Myr (1 Gyr)
    const SUBSTEPS = 2; // 2 Velocity-Verlet substeps per 1 Myr frame
    const DT_SUB = SIM_DT_MYR / SUBSTEPS;
        const KM_S_TO_KPC_MYR = 0.0010227121650537077;

        /* ---------- Physics engine: gravitational potentials of a custom galaxy ----------
           Pure functions of the anatomy sliders (no DOM). Each component contributes its circular speed
           V_c^2 = R dPhi/dR in (km/s)^2, and the components add in quadrature:
             stellar disk, gas  exponential disks, Sigma(R) = Sigma_0 exp(-R/R_s) cut at R_max. Inside R_max: Freeman's (1970)
                                Bessel-function solution for the infinite disk, minus the annulus beyond R_max (a sum of
                                softened elliptic-integral rings). Outside R_max: the rings out to R_max.
             bulge              Hernquist sphere, Phi = -G M / (r + r_b), so V_c^2 = G M r / (r + r_b)^2.
           The halo baseline couples to the stellar mass through the Moster et al. (2013) stellar-to-halo mass relation
           (inverted for M200) and the Dutton & Maccio (2014) mass-concentration relation (giving the NFW scale radius). */
        const PhysicsEngine=(()=>{
            const TWO_PI=2*Math.PI;
            const RHO_CRIT=3*H0*H0/(8*Math.PI*G);
            const GAS_SCALE_RATIO=3; // the gas disk is this many times more extended than the stellar disk
            let ellipK=0,ellipE=0;
            // Complete elliptic integrals K(m) and E(m), m = k^2, by the arithmetic-geometric mean.
            function ellipticKE(m){
                let a=1,b=Math.sqrt(1-m),sum=m/2,weight=1;
                for(let i=0;i<40;i++){
                    const c=(a-b)/2,next=(a+b)/2;
                    b=Math.sqrt(a*b);a=next;
                    sum+=weight*c*c;weight*=2;
                    if(c<1e-16*a)break;
                }
                ellipK=Math.PI/(2*a);ellipE=ellipK*(1-sum);
            }
            // Modified Bessel functions (Abramowitz & Stegun 9.8.1-9.8.8, relative error ~1e-7).
            const poly=(t,c)=>{let s=0;for(let i=c.length-1;i>=0;i--)s=s*t+c[i];return s;};
            const I0_SMALL=[1,3.5156229,3.0899424,1.2067492,0.2659732,0.0360768,0.0045813];
            const I0_LARGE=[0.39894228,0.01328592,0.00225319,-0.00157565,0.00916281,-0.02057706,0.02635537,-0.01647633,0.00392377];
            const I1_SMALL=[0.5,0.87890594,0.51498869,0.15084934,0.02658733,0.00301532,0.00032411];
            const I1_LARGE=[0.39894228,-0.03988024,-0.00362018,0.00163801,-0.01031555,0.02282967,-0.02895312,0.01787654,-0.00420059];
            const K0_SMALL=[-0.57721566,0.42278420,0.23069756,0.03488590,0.00262698,0.00010750,0.00000740];
            const K0_LARGE=[1.25331414,-0.07832358,0.02189568,-0.01062446,0.00587872,-0.00251540,0.00053208];
            const K1_SMALL=[1,0.15443144,-0.67278579,-0.18156897,-0.01919402,-0.00110404,-0.00004686];
            const K1_LARGE=[1.25331414,0.23498619,-0.03655620,0.01504268,-0.00780353,0.00325614,-0.00068245];
            const besselI0=x=>x<3.75?poly((x/3.75)**2,I0_SMALL):Math.exp(x)/Math.sqrt(x)*poly(3.75/x,I0_LARGE);
            const besselI1=x=>x<3.75?x*poly((x/3.75)**2,I1_SMALL):Math.exp(x)/Math.sqrt(x)*poly(3.75/x,I1_LARGE);
            const besselK0=x=>x<=2?-Math.log(x/2)*besselI0(x)+poly((x/2)**2,K0_SMALL):Math.exp(-x)/Math.sqrt(x)*poly(2/x,K0_LARGE);
            const besselK1=x=>x<=2?Math.log(x/2)*besselI1(x)+poly((x/2)**2,K1_SMALL)/x:Math.exp(-x)/Math.sqrt(x)*poly(2/x,K1_LARGE);
            // Gauss-Legendre nodes and weights on [-1, 1].
            function gaussLegendre(n){
                const x=new Float64Array(n),w=new Float64Array(n);
                for(let i=0;i<n;i++){
                    let z=Math.cos(Math.PI*(i+0.75)/(n+0.5)),dp=1;
                    for(let it=0;it<100;it++){
                        let p1=1,p2=0;
                        for(let j=0;j<n;j++){const p3=p2;p2=p1;p1=((2*j+1)*z*p2-j*p3)/(j+1);}
                        dp=n*(z*p1-p2)/(z*z-1);
                        const dz=p1/dp;z-=dz;
                        if(Math.abs(dz)<1e-15)break;
                    }
                    x[i]=z;w[i]=2/((1-z*z)*dp*dp);
                }
                return{x,w};
            }
            const GL=gaussLegendre(10),PANELS=7;
            // (dPhi/dR)/G for a thin ring of unit mass and radius a, seen at in-plane radius R with Plummer softening eps.
            // Positive pulls inwards. Phi_ring = -(2/pi) K(m) / sqrt(s), s = (a+R)^2 + eps^2, m = 4aR/s.
            function ringGradient(a,R,eps2){
                const s=(a+R)*(a+R)+eps2,m=4*a*R/s;
                ellipticKE(m);
                const mc=1-m,dm=4*a*(a*a-R*R+eps2)/(s*s);
                const dK=m<1e-6?Math.PI/8*(1+9*m/16):(ellipE-mc*ellipK)/(2*m*mc);
                const root=Math.sqrt(s);
                return -(2/Math.PI)*(dK*dm/root-ellipK*(a+R)/(s*root));
            }
            // Exponential disk of the given mass inside r_max (so the slider is the mass actually present).
            function exponentialDisk(mass,scale,rMax,eps){
                const xm=rMax/scale,sigma0=mass/(TWO_PI*scale*scale*(1-(1+xm)*Math.exp(-xm)));
                const eps2=eps*eps;
                // integral of 2 pi a Sigma(a) (dPhi_ring/dR)/G over [lo, hi]: composite Gauss-Legendre whose panels
                // shrink geometrically towards the end nearest R, where the ring kernel is sharpest.
                function ringSum(R,lo,hi,sharpEnd){
                    let sum=0,far=hi-lo;
                    for(let p=0;p<PANELS;p++){
                        const near=p===PANELS-1?0:far/2;
                        const from=sharpEnd==='high'?hi-far:lo+near,to=sharpEnd==='high'?hi-near:lo+far;
                        const mid=(from+to)/2,half=(to-from)/2;
                        for(let i=0;i<GL.x.length;i++){
                            const a=mid+half*GL.x[i];
                            sum+=GL.w[i]*half*a*Math.exp(-a/scale)*ringGradient(a,R,eps2);
                        }
                        far=near;
                    }
                    return TWO_PI*sigma0*sum;
                }
                function vSq(R){
                    if(!(R>0))return 0;
                    if(R<rMax){
                        const y=R/(2*scale);
                        const freeman=4*Math.PI*G*sigma0*scale*y*y*(besselI0(y)*besselK0(y)-besselI1(y)*besselK1(y));
                        return Math.max(0,freeman-G*R*ringSum(R,rMax,rMax+16*scale,'low'));
                    }
                    return Math.max(0,G*R*ringSum(R,0,rMax,'high'));
                }
                return{sigma0,vSq};
            }
            // Circular-speed functions (km/s squared) of a galaxy with the given anatomy.
            function buildGalaxy(anatomy,eps){
                const stars=exponentialDisk(anatomy.diskMass,anatomy.diskScale,anatomy.rMax,eps);
                const gas=exponentialDisk(anatomy.gasMass,GAS_SCALE_RATIO*anatomy.diskScale,anatomy.rMax,eps);
                const bulgeMass=anatomy.bulgeMass>0?anatomy.bulgeMass:0,rb=anatomy.bulgeRadius;
                const bulgeSq=r=>bulgeMass>0?G*bulgeMass*r/((r+rb)*(r+rb)):0;
                return{diskSq:stars.vSq,gasSq:gas.vSq,bulgeSq,vSq:r=>stars.vSq(r)+gas.vSq(r)+bulgeSq(r)};
            }
            // Moster et al. (2013), z = 0: M*/M_h = 2 N [(M_h/M_1)^-beta + (M_h/M_1)^gamma]^-1.
            const SHMR={logM1:11.59,norm:0.0351,beta:1.376,gamma:0.608};
            function stellarMassOfHalo(m200){
                const x=m200/10**SHMR.logM1;
                return m200*2*SHMR.norm/(x**-SHMR.beta+x**SHMR.gamma);
            }
            function haloMassFromStellar(mStar){
                let lo=7,hi=16;
                for(let i=0;i<60;i++){const mid=(lo+hi)/2;if(stellarMassOfHalo(10**mid)<mStar)lo=mid;else hi=mid;}
                return 10**((lo+hi)/2);
            }
            // Dutton & Maccio (2014): log10 c200 = 0.905 - 0.101 log10(M200 / (10^12 h^-1 Msun)), with h = 0.7.
            const concentrationOf=m200=>10**(0.905-0.101*Math.log10(m200/(1e12/0.7)));
            const virialRadius=m200=>Math.cbrt(3*m200/(4*Math.PI*200*RHO_CRIT));
            function baselineHalo(stellarMass){
                const m200=haloMassFromStellar(stellarMass);
                return{m200,scaleRadius:virialRadius(m200)/concentrationOf(m200)};
            }
            // Rough anatomy of a SPARC galaxy from its mass-model speeds; it only shapes the 3D star population.
            function anatomyFromRotation(rotation,upsilonDisk,upsilonBulge){
                const r=rotation.radius_kpc,n=r.length,rLast=r[n-1];
                const vDisk=rotation.v_disk_kms.map(Math.abs),vBulge=rotation.v_bulge_kms.map(Math.abs);
                let peak=0;
                for(let i=1;i<n;i++)if(vDisk[i]>vDisk[peak])peak=i;
                // An exponential disk's speed peaks at 2.16 R_d, where V^2 = 0.3872 G M / R_d (Freeman 1970).
                const diskScale=r[peak]/2.16,diskMass=upsilonDisk*vDisk[peak]**2*diskScale/(0.3872*G);
                // A spherical bulge's enclosed mass is V^2 r / G; its Hernquist half-mass radius is (1 + sqrt 2) r_b.
                const enclosed=r.map((x,i)=>upsilonBulge*vBulge[i]**2*x/G),bulgeMass=enclosed[n-1];
                let halfRadius=r[0];
                for(let i=0;i<n;i++){
                    if(enclosed[i]>=bulgeMass/2){
                        const previous=i>0?i-1:0,span=enclosed[i]-enclosed[previous]||1;
                        halfRadius=i>0?r[previous]+(r[i]-r[previous])*(bulgeMass/2-enclosed[previous])/span:r[0];
                        break;
                    }
                }
                const gasMass=rotation.v_gas_kms[n-1]**2*rLast/G;
                return{diskMass,diskScale,rMax:rLast,bulgeMass:bulgeMass>1e6?bulgeMass:0,bulgeRadius:Math.max(0.05,halfRadius/(1+Math.SQRT2)),gasMass};
            }
            return{buildGalaxy,baselineHalo,anatomyFromRotation,GAS_SCALE_RATIO};
        })();

        /* ---------- Test-star population ----------
           The 800-star NGC 3198 disk is the template: other galaxies scale it with their visible mass (starCountFor).
           Each star owns fixed random numbers (STAR_SEEDS), so moving a slider slides the stars smoothly instead of
           reshuffling them. Stars are ordered bulge first, then disk. */
        const PM=DATA.particle_model;
        const MAX_STARS=PM.max_count,MIN_BULGE_STARS=12,MIN_DISK_STARS=24,BULGE_ROTATION=0.3;
        const GALAXY={scaleHeight:0.15,diskMin:0.45,bulgeMin:0.1,blackHoleRadius:0.16};
        function mulberry32(seed){return()=>{seed|=0;seed=seed+0x6D2B79F5|0;let t=Math.imul(seed^seed>>>15,1|seed);
            t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
        // Box-Muller transform: converts pairs of uniform random numbers to standard normal distribution
        function boxMuller(mean,stdDev,rng){
            let u1=rng();
            while(u1<=1e-15)u1=rng();
            const u2=rng();
            return mean+stdDev*Math.sqrt(-2.0*Math.log(u1))*Math.cos(2.0*Math.PI*u2);
        }
        // Fixed uniform and Gaussian random numbers for every star the page can ever show.
        const STAR_SEEDS=(()=>{
            const rng=mulberry32(3198),uniform=new Float32Array(MAX_STARS*4),gauss=new Float32Array(MAX_STARS*4);
            for(let i=0;i<MAX_STARS*4;i++){uniform[i]=rng();gauss[i]=boxMuller(0,1,rng);}
            return{uniform,gauss};
        })();
        const starCountFor=(visibleMass,templateMass)=>Math.max(PM.min_count,Math.min(PM.max_count,
            Math.round(PM.base_count*Math.pow(visibleMass/templateMass,PM.count_exponent))));
        const visibleMassOf=anatomy=>anatomy.diskMass+anatomy.bulgeMass+anatomy.gasMass;
        const presetAnatomy=galaxy=>PhysicsEngine.anatomyFromRotation(galaxy.rotation,DATA.upsilon_disk,DATA.upsilon_bulge);
        const PRESET_TEMPLATE_MASS=visibleMassOf(presetAnatomy(GALAXIES[DATA.default_galaxy_id]));
        // The sandbox template is its own slider defaults, an NGC 3198-like galaxy.
        const sandboxTemplateMass=()=>['sbDiskMass','sbBulgeMass','sbGasMass'].reduce((sum,id)=>sum+Number(document.getElementById(id).defaultValue),0);
        // Radius enclosing the fraction u of an exponential disk's stars (cumulative 1-(1+x)e^-x, x = R/R_d) between rMin and rMax.
        function diskRadius(u,scale,rMin,rMax){
            const cdf=x=>1-(1+x)*Math.exp(-x),c0=cdf(rMin/scale),target=c0+u*(cdf(rMax/scale)-c0);
            let lo=rMin/scale,hi=rMax/scale;
            for(let i=0;i<40;i++){const mid=(lo+hi)/2;if(cdf(mid)<target)lo=mid;else hi=mid;}
            return scale*(lo+hi)/2;
        }
        // Radius enclosing the fraction u of a Hernquist bulge's stars (cumulative (r/(r+a))^2) between rMin and rMax.
        function hernquistRadius(u,a,rMin,rMax){
            const fraction=r=>(r/(r+a))**2,f0=fraction(rMin),s=Math.sqrt(f0+u*(fraction(rMax)-f0));
            return a*s/(1-s);
        }
        // Positions of `count` stars for an anatomy: an exponential disk to R_max and, when there is a bulge, a
        // spheroidal Hernquist cluster (its share of the stars follows its share of the stellar mass).
        function generatePopulation(anatomy,count){
            let bulgeCount=0;
            if(anatomy.bulgeMass>0){
                const share=anatomy.bulgeMass/(anatomy.diskMass+anatomy.bulgeMass);
                bulgeCount=Math.min(count-MIN_DISK_STARS,Math.max(MIN_BULGE_STARS,Math.round(count*share)));
            }
            const bulgeCut=Math.min(anatomy.rMax,25*anatomy.bulgeRadius);
            const u=STAR_SEEDS.uniform,g=STAR_SEEDS.gauss,list=new Array(count);
            for(let k=0;k<count;k++){
                const phi=2*Math.PI*u[4*k+1],kick=[g[4*k+1],g[4*k+2],g[4*k+3]];
                if(k<bulgeCount){
                    const r=hernquistRadius(u[4*k],anatomy.bulgeRadius,GALAXY.bulgeMin,bulgeCut);
                    const cosT=2*u[4*k+2]-1,sinT=Math.sqrt(1-cosT*cosT);
                    list[k]={kind:'bulge',position:[r*sinT*Math.cos(phi),r*sinT*Math.sin(phi),r*cosT],radius:r,kick,colorSeed:u[4*k+3]};
                }else{
                    const R=diskRadius(u[4*k],anatomy.diskScale,GALAXY.diskMin,anatomy.rMax);
                    list[k]={kind:'disk',position:[R*Math.cos(phi),R*Math.sin(phi),g[4*k]*GALAXY.scaleHeight],radius:R,kick,colorSeed:u[4*k+3]};
                }
            }
            // Number the stars from the centre outwards (by starting distance), so star 1 is innermost and the last is outermost.
            list.sort((a,b)=>Math.hypot(...a.position)-Math.hypot(...b.position));
            return{particles:list,bulgeCount,bulgeCut};
        }
        /* ---------- Active galaxy: SPARC rotation curve ----------
           galaxies.js stores the raw SPARC columns. loadGalaxyData() rebuilds every array below for the chosen
           galaxy, so the rest of the page reads radii / vObs / vBary without knowing which galaxy is active. */
        const mlState={disk:DATA.upsilon_disk,bulge:DATA.upsilon_bulge};
        /* A best fit may also use a distance and tilt (inclination) within SPARC's errors. As factors on SPARC's values:
           f = D / D_SPARC scales the radii and the squared visible-matter speeds; s = sin(i_SPARC) / sin(i) scales the
           measured speeds. The page keeps SPARC's points where they are and divides the model by s instead, which gives
           the same chi-squared (see _chi_squared in the generator). SPARC_FRAME (f = s = 1) is SPARC's own values. */
        const SPARC_FRAME={f:1,s:1};
        let fitFrame=SPARC_FRAME;
        const frameOf=(galaxy,fit)=>fit?{f:fit.distance_mpc/galaxy.catalogue.distance_mpc,
            s:Math.sin(galaxy.catalogue.inclination_deg*Math.PI/180)/Math.sin(fit.inclination_deg*Math.PI/180)}:SPARC_FRAME;
        let rotationRows=[],radii=[],vObs=[],obsV=[],obsErr=[],vBary=[];
        function sparcRowsOf(galaxy){
            const t=galaxy.rotation;
            return t.radius_kpc.map((radius,i)=>{
                const gas=t.v_gas_kms[i];
                return{radius_kpc:radius,observed_kms:Math.abs(t.v_obs_kms[i]),error_kms:Math.abs(t.v_obs_err_kms[i]),
                    gas_kms:Math.abs(gas),gas_sq:gas*Math.abs(gas),disk_kms:Math.abs(t.v_disk_kms[i]),bulge_kms:Math.abs(t.v_bulge_kms[i])};
            });
        }
        // Visible-matter speed at each SPARC row for the given disk and bulge mass-to-light ratios (the gas term keeps its sign).
        const sparcBaryons=(rows,ud,ub,frame=SPARC_FRAME)=>rows.map(p=>Math.sqrt(frame.f*Math.max(0,p.gas_sq+ud*p.disk_kms**2+ub*p.bulge_kms**2))/frame.s);
        const baryonsFor=(ud,ub)=>customActive&&customModel
            ?rotationRows.map(p=>Math.sqrt(customModel.vSq(p.radius_kpc)))
            :sparcBaryons(rotationRows,ud,ub,fitFrame);
        function loadGalaxyData(galaxy){
            rotationRows=sparcRowsOf(galaxy);
            radii=rotationRows.map(p=>p.radius_kpc);
            vObs=rotationRows.map(p=>p.observed_kms);
            obsV=vObs.slice();
            obsErr=rotationRows.map(p=>p.error_kms);
            vBary=baryonsFor(mlState.disk,mlState.bulge);
        }
        loadGalaxyData(activeGalaxy);
        // The Custom Sandbox galaxy is imaginary, so it has no measured speeds: its curves are drawn on an even
        // grid of 60 radii out to 1.25 x R_max instead of at a SPARC galaxy's data radii. A few extra radii inside
        // the first grid step (CUSTOM_GRID_INNER, as fractions of that step) let the curve start right at the centre,
        // where every speed goes to zero, instead of jumping in at the first step. r = 0 itself is left out, since
        // some speed and gravity formulas divide by r.
        const CUSTOM_GRID_INNER=[0.05,0.25,0.5,0.75];
        function loadCustomGrid(rMax){
            const n=60,rEnd=1.25*rMax,step=rEnd/n;
            const gridRadii=[...CUSTOM_GRID_INNER.map(f=>f*step),...Array.from({length:n},(_,i)=>step*(i+1))];
            rotationRows=gridRadii.map(r=>({radius_kpc:r,observed_kms:0,error_kms:0,gas_kms:0,gas_sq:0,disk_kms:0,bulge_kms:0}));
            radii=rotationRows.map(p=>p.radius_kpc);
            vObs=radii.map(()=>0);obsV=vObs.slice();obsErr=vObs.slice();
            vBary=radii.map(()=>0);
        }
        function observedVelocity(radius){
            const last=radii.length-1;
            if(radius<=radii[0])return vObs[0]*radius/radii[0];
            if(radius>=radii[last])return vObs[last];
            let low=0,high=last;
            while(high-low>1){const m=(low+high)>>1;if(radii[m]<=radius)low=m;else high=m;}
            return vObs[low]+(radius-radii[low])/(radii[high]-radii[low])*(vObs[high]-vObs[low]);
        }
        function nfwVelocityModel(massSolar,scaleRadius,frame=SPARC_FRAME) {
            if(massSolar<=0)return ()=>0;
            const rhoCrit=3*H0*H0/(8*Math.PI*G);
            const r200=Math.cbrt(3*massSolar/(4*Math.PI*200*rhoCrit));
            const shape=x=>x<=0?0:Math.log1p(x)-x/(1+x);
            const fc=Math.max(shape(r200/scaleRadius),1e-12);
            return radius=>{
                const r=radius*frame.f;
                if(r<=0)return 0;
                const enclosed=r>=r200?massSolar:massSolar*shape(r/scaleRadius)/fc;
                return Math.sqrt(G*enclosed/r)/frame.s;
            };
        }
        function haloCurve(massT,scaleRadius) {return radii.map(nfwVelocityModel(massT*1e12,scaleRadius,fitFrame));}
        const M_PER_KPC = 3.08567758e19;
        const MOND_A0_DEFAULT = 1.2e-10; // m/s^2
        let MOND_A0 = MOND_A0_DEFAULT;
        let mondMode = false;
        // SI calculation with the RAR interpolating function, g = g_N / (1 − e^−√(g_N/a₀)) (McGaugh, Lelli & Schombert 2016,
        // eq. 4), the form a₀ was measured with. −expm1(−x) is 1 − e^−x without rounding loss when x is tiny.
        function mondVelocity(vBarKms, rKpc, a0 = MOND_A0) {
            const rM = rKpc * M_PER_KPC, vMs = vBarKms * 1000;
            if (rM <= 0 || vMs <= 0) return 0;
            const gN = vMs * vMs / rM;
            const g = gN / -Math.expm1(-Math.sqrt(gN / a0));
            return Math.sqrt(g * rM) / 1000;
        }
        // MOND prediction drawn against SPARC's points, for a drawn visible-matter speed and a best fit's distance and tilt.
        const mondModelVelocity = (vBarKms, rKpc, a0 = MOND_A0, frame = fitFrame) => mondVelocity(vBarKms * frame.s, rKpc * frame.f, a0) / frame.s;
        // Velocity added on top of the baryons (in quadrature): NFW halo, or the MOND boost.
        function extraVelocityModel(massSolar, scaleRadius, vbOf) {
            if (!mondMode) return nfwVelocityModel(massSolar, scaleRadius, fitFrame);
            return radius => {
                if (radius <= 0) return 0;
                const vb = vbOf(radius), vm = mondModelVelocity(vb, radius);
                return Math.sqrt(Math.max(vm * vm - vb * vb, 0));
            };
        }
        // Baryonic circular speed: the custom galaxy's engine table, or the SPARC mass model interpolated at the data radii.
        let customBaryonTable=null;
        function baryonicVelocity(radius) {
            if(customActive&&customBaryonTable){
                const t=radius*INV_DR_TABLE;
                if(t>=ACCEL_TABLE_SIZE-1)return Math.sqrt(customBaryonTable[ACCEL_TABLE_SIZE-1]*R_MAX_TABLE/radius);
                const i=t|0,f=t-i;
                return Math.sqrt(customBaryonTable[i]+f*(customBaryonTable[i+1]-customBaryonTable[i]));
            }
            if(radius<=radii[0])return vBary[0]*radius/radii[0];
            const last=radii.length-1;
            if(radius>=radii[last])return vBary[last]*Math.sqrt(radii[last]/radius);
            let low=0,high=last;
            while(high-low>1){const middle=Math.floor((low+high)/2);if(radii[middle]<=radius)low=middle;else high=middle;}
            const fraction=(radius-radii[low])/(radii[high]-radii[low]);
            return vBary[low]+fraction*(vBary[high]-vBary[low]);
        }

        // Radial acceleration tables: a = -k(r) r_vec (k in 1/Myr^2) for the total potential of baryons + halo.
        const ACCEL_TABLE_SIZE = 2048;
        const R_MAX_TABLE = 160.0;
        const INV_DR_TABLE = (ACCEL_TABLE_SIZE - 1) / R_MAX_TABLE;
        const DR_TABLE = R_MAX_TABLE / (ACCEL_TABLE_SIZE - 1);
        const SOFTENING_SPLINE = 2.8; // the cubic-spline force is exactly Newtonian beyond 2.8 epsilon (Plummer-equivalent epsilon)
        // Gravitational softening: the fraction of the Newtonian force left at u = r / (2.8 epsilon), which smooths the
        // dense centre so no star meets a force singularity or a time step too coarse for its orbit.
        function softeningFactor(u) {
            if (u >= 1) return 1;
            const u3 = u * u * u;
            return (u < 0.5
                ? 10.666666666667 + u * u * (32 * u - 38.4)
                : 21.333333333333 - 48 * u + 38.4 * u * u - 10.666666666667 * u3 - 0.066666666667 / u3) * u3;
        }
        // A pure disk gets the smallest softening; with a bulge it follows the bulge's scale radius within the allowed range.
        const softeningFor = anatomy => anatomy.bulgeMass > 0
            ? Math.min(PM.softening_kpc[1], Math.max(PM.softening_kpc[0], anatomy.bulgeRadius))
            : PM.softening_kpc[0];
        let softening = PM.softening_kpc[0];
        const forceTable = new Float64Array(ACCEL_TABLE_SIZE), referenceTable = new Float64Array(ACCEL_TABLE_SIZE);
        function fillForceTable(table, speedSquaredAt) {
            const scale = SOFTENING_SPLINE * softening;
            for (let i = 0; i < ACCEL_TABLE_SIZE; i++) {
                const r = Math.max(i * DR_TABLE, 1e-4);
                table[i] = speedSquaredAt(r) * (KM_S_TO_KPC_MYR ** 2) * softeningFactor(r / scale) / (r * r);
            }
        }
        function accelFactor(table, r) {
            const idxFloat = r * INV_DR_TABLE;
            if (idxFloat >= ACCEL_TABLE_SIZE - 1) {
                const ratio = R_MAX_TABLE / r;
                return table[ACCEL_TABLE_SIZE - 1] * ratio * ratio * ratio;
            }
            const idx = idxFloat | 0;
            return table[idx] + (idxFloat - idx) * (table[idx + 1] - table[idx]);
        }

        /* ---------- Active population and its initial orbits ----------
           Stars start on the orbits that balance a reference galaxy: the observed SPARC curve of a preset, or the
           baryons plus the baseline halo of a custom galaxy. Moving the halo sliders or switching to MOND then changes
           the gravity under them. Disk stars begin on circular orbits; bulge stars get an isotropic velocity dispersion
           from the Jeans equation, so the spheroid is held up by pressure instead of rotation. */
        const presetReferenceSpeedSquared = r => observedVelocity(r) ** 2;
        let referenceSpeedSquared = presetReferenceSpeedSquared;
        let particles = [], starCount = 0, bulgeStarCount = 0, bulgeCutRadius = 1, currentAnatomy = null, ejectionRadius = 1;
        const JEANS_POINTS = 240, JEANS_R_MIN = 0.01, JEANS_R_MAX = 1500;
        // sigma^2(r) = (1/rho) integral_r^inf rho g dr' for a Hernquist tracer (rho ~ 1/(r (r+a)^3)) in the reference gravity g.
        function bulgeDispersionSquared(scale) {
            const lnStep = Math.log(JEANS_R_MAX / JEANS_R_MIN) / (JEANS_POINTS - 1), sigma2 = new Float64Array(JEANS_POINTS);
            let tail = 0, previous = 0;
            for (let i = JEANS_POINTS - 1; i >= 0; i--) {
                const r = JEANS_R_MIN * Math.exp(i * lnStep), rho = 1 / (r * (r + scale) ** 3);
                const flux = rho * accelFactor(referenceTable, r) * r * r;
                if (i < JEANS_POINTS - 1) tail += 0.5 * (flux + previous) * lnStep;
                previous = flux;
                sigma2[i] = tail / rho;
            }
            return r => {
                const t = Math.log(Math.max(r, JEANS_R_MIN) / JEANS_R_MIN) / lnStep, i = Math.min(Math.floor(t), JEANS_POINTS - 2);
                return Math.max(0, sigma2[i] + Math.min(1, t - i) * (sigma2[i + 1] - sigma2[i]));
            };
        }
        function assignInitialVelocities() {
            fillForceTable(referenceTable, referenceSpeedSquared);
            const dispersion = bulgeStarCount > 0 ? bulgeDispersionSquared(currentAnatomy.bulgeRadius) : null;
            const streaming = Math.sqrt(1 - BULGE_ROTATION ** 2);
            for (const p of particles) {
                const [x, y, z] = p.position;
                if (p.kind === 'disk') {
                    const rho = Math.max(Math.hypot(x, y), 1e-6);
                    const speed = Math.sqrt(accelFactor(referenceTable, Math.hypot(x, y, z))) * rho;
                    p.velocity = [-y / rho * speed, x / rho * speed, 0];
                } else {
                    const sigma = Math.sqrt(dispersion(p.radius)), rho = Math.hypot(x, y);
                    const c = rho > 1e-9 ? x / rho : 1, s = rho > 1e-9 ? y / rho : 0;
                    const vR = sigma * p.kick[0], vPhi = sigma * (streaming * p.kick[1] + BULGE_ROTATION), vZ = sigma * p.kick[2];
                    p.velocity = [vR * c - vPhi * s, vR * s + vPhi * c, vZ];
                }
            }
        }
        // Builds the test stars for an anatomy. referenceSpeedSquared must already describe the galaxy.
        function setPopulation(anatomy, count) {
            currentAnatomy = anatomy;
            softening = softeningFor(anatomy);
            ejectionRadius = PM.ejection_factor * anatomy.rMax;
            const population = generatePopulation(anatomy, count);
            particles = population.particles;
            bulgeStarCount = population.bulgeCount;
            bulgeCutRadius = population.bulgeCut;
            starCount = particles.length;
            assignInitialVelocities();
        }

        let sortedEscapeTimes = [], escapeTimeByStar = new Float64Array(0);
        let activeHaloMassT = DEFAULTS.massT;

        // One shared buffer holds every star's 1 Gyr path (frameCount x 3 floats each); it only grows.
        let trajectoryBuffer = null, trajectoryViews = [];
        function trajectoryViewsFor(count) {
            if (trajectoryViews.length < count) {
                const stride = frameCount * 3;
                trajectoryBuffer = new Float32Array(count * stride);
                trajectoryViews = Array.from({length: count}, (_, k) => trajectoryBuffer.subarray(k * stride, (k + 1) * stride));
            }
            return trajectoryViews;
        }

        // Velocity-Verlet symplectic integration loop across the 1 Gyr (1000 Myr) timeline
        function integrateTrajectories(massT, scaleRadius) {
            const haloVelocity = extraVelocityModel(massT * 1e12, scaleRadius, baryonicVelocity);
            fillForceTable(forceTable, r => {
                const vb = baryonicVelocity(r), vh = haloVelocity(r);
                return vb * vb + vh * vh;
            });
            const count = particles.length, trajectories = trajectoryViewsFor(count);
            const dtSub = DT_SUB;
            const dtSub2Half = 0.5 * dtSub * dtSub;
            const dtSubHalf = 0.5 * dtSub;
            const limit = ejectionRadius, limit2 = limit * limit;
            const escapeTimes = [];
            if (escapeTimeByStar.length < count) escapeTimeByStar = new Float64Array(count);

            for (let p = 0; p < count; p++) {
                const particle = particles[p], positions = trajectories[p];
                let x = particle.position[0], y = particle.position[1], z = particle.position[2];
                let vx = particle.velocity[0], vy = particle.velocity[1], vz = particle.velocity[2];
                let k = accelFactor(forceTable, Math.max(Math.sqrt(x * x + y * y + z * z), 1e-4));
                let ax = -k * x, ay = -k * y, az = -k * z;
                positions[0] = x; positions[1] = y; positions[2] = z;

                // A star is ejected once it is beyond the ejection radius (1.5 x R_max); the crossing time is kept.
                let escTime = x * x + y * y + z * z > limit2 ? 0.0 : Infinity;

                for (let frame = 1; frame < frameCount; frame++) {
                    const prevR2 = x * x + y * y + z * z;
                    for (let s = 0; s < SUBSTEPS; s++) {
                        // Velocity-Verlet step 1: advance position
                        x += vx * dtSub + ax * dtSub2Half;
                        y += vy * dtSub + ay * dtSub2Half;
                        z += vz * dtSub + az * dtSub2Half;

                        // Velocity-Verlet step 2: evaluate acceleration at new position
                        k = accelFactor(forceTable, Math.max(Math.sqrt(x * x + y * y + z * z), 1e-4));
                        const nax = -k * x, nay = -k * y, naz = -k * z;

                        // Velocity-Verlet step 3: advance velocity
                        vx += (ax + nax) * dtSubHalf;
                        vy += (ay + nay) * dtSubHalf;
                        vz += (az + naz) * dtSubHalf;

                        ax = nax; ay = nay; az = naz;
                    }

                    const currR2 = x * x + y * y + z * z;
                    if (escTime === Infinity && currR2 > limit2) {
                        const rPrev = Math.sqrt(prevR2), rCurr = Math.sqrt(currR2);
                        const frac = Math.max(0, Math.min(1, (limit - rPrev) / Math.max(rCurr - rPrev, 1e-6)));
                        escTime = (frame - 1 + frac) * SIM_DT_MYR;
                    }

                    const idx = frame * 3;
                    positions[idx] = x; positions[idx + 1] = y; positions[idx + 2] = z;
                }

                escapeTimeByStar[p] = escTime; // kept per star so the tutorial can follow the first one out
                if (escTime < Infinity) {
                    escapeTimes.push(escTime);
                }
            }

            sortedEscapeTimes = escapeTimes.sort((a, b) => a - b);
            activeHaloMassT = mondMode ? 0 : massT;
            return trajectories;
        }
        setPopulation(presetAnatomy(activeGalaxy), starCountFor(visibleMassOf(presetAnatomy(activeGalaxy)), PRESET_TEMPLATE_MASS));
        let activeTrajectories = integrateTrajectories(DEFAULTS.massT, DEFAULTS.scaleRadius);
        let currentFrame=0,orbitModelTimer=null;
        const orbitPlot=document.getElementById('orbitPlot');

        /* ---------- Camera / input module: damped orbit controls (mouse, pen, touch) ---------- */
        class SmoothOrbitControls {
            // Z-up spherical orbit. Input only edits `goal`; update(dt) eases `view` toward it
            // with a frame-rate independent exponential filter, so nothing ever snaps.
            constructor(camera,element,options={}) {
                this.camera=camera;this.element=element;this.dragMode='rotate';
                this.o=Object.assign({damping:11,touchDamping:5,rotateSpeed:1,panSpeed:1,zoomSpeed:0.0016,inertia:0.14,
                    minDistance:5,maxDistance:1000,minPhi:0.04,maxPhi:Math.PI-0.04,panLimit:80,
                    initial:{theta:-0.82,phi:1.01,radius:72,target:[0,0,0]},onChange:()=>{}},options);
                const i=this.o.initial;
                this.goal={theta:i.theta,phi:i.phi,radius:i.radius,target:new THREE.Vector3(...i.target)};
                this.view={theta:i.theta,phi:i.phi,radius:i.radius,target:this.goal.target.clone()};
                this.pointers=new Map();this.pinch=null;this.lastMove=0;this.velocity={theta:0,phi:0};
                this._right=new THREE.Vector3();this._up=new THREE.Vector3();
                this.enabled=true; // false while the guided tour is open: the camera only moves when a card asks
                this._bind();this._apply();
            }
            _bind() {
                const el=this.element;
                this._h={
                    down:e=>this._down(e),move:e=>this._move(e),up:e=>this._up_(e),
                    wheel:e=>{if(this.enabled)this._wheel(e);},menu:e=>e.preventDefault(),dbl:()=>{if(this.enabled)this.reset();}
                };
                el.addEventListener('pointerdown',this._h.down);
                el.addEventListener('pointermove',this._h.move);
                el.addEventListener('pointerup',this._h.up);
                el.addEventListener('pointercancel',this._h.up);
                el.addEventListener('lostpointercapture',this._h.up);
                el.addEventListener('wheel',this._h.wheel,{passive:false});
                el.addEventListener('contextmenu',this._h.menu);
                el.addEventListener('dblclick',this._h.dbl);
            }
            dispose() {
                const el=this.element,h=this._h;
                el.removeEventListener('pointerdown',h.down);el.removeEventListener('pointermove',h.move);
                el.removeEventListener('pointerup',h.up);el.removeEventListener('pointercancel',h.up);
                el.removeEventListener('lostpointercapture',h.up);el.removeEventListener('wheel',h.wheel);
                el.removeEventListener('contextmenu',h.menu);el.removeEventListener('dblclick',h.dbl);
            }
            // Eases to a preset; theta takes the shortest way round so the camera never spins the long way.
            flyTo({theta,phi,radius,target,rate}={}) {
                const g=this.goal;
                this.flyRate=rate||0;
                if(theta!==undefined){
                    const turn=Math.PI*2;
                    g.theta=theta+Math.round((this.view.theta-theta)/turn)*turn;
                }
                if(phi!==undefined)g.phi=this._clampPhi(phi);
                if(radius!==undefined)g.radius=this._clampRadius(radius);
                if(target)g.target.set(...target);
                this.velocity.theta=this.velocity.phi=0;
                this.o.onChange();
            }
            reset() {this.flyTo(Object.assign({rate:this.view.radius>200?3:0},this.o.initial));}
            _down(e) {
                if(!this.enabled||(e.pointerType==='mouse'&&e.button>2))return;
                try{this.element.setPointerCapture(e.pointerId);}catch(_){}
                this.pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,button:e.button,shift:e.shiftKey||e.ctrlKey||e.metaKey,type:e.pointerType});
                this.touchActive=e.pointerType==='touch';
                this.flyRate=0;
                this.velocity.theta=this.velocity.phi=0;
                if(this.pointers.size===2)this.pinch=this._pinchState();
                this.element.classList.add('dragging');
                this.o.onChange();
            }
            _pinchState() {
                const [a,b]=[...this.pointers.values()];
                return{dist:Math.hypot(a.x-b.x,a.y-b.y)||1,x:(a.x+b.x)/2,y:(a.y+b.y)/2};
            }
            _move(e) {
                const p=this.pointers.get(e.pointerId);
                if(!p)return;
                // Coalesced events keep fast swipes smooth on high-polling-rate devices.
                const events=e.getCoalescedEvents?e.getCoalescedEvents():[];
                for(const ev of (events.length?events:[e])) {
                    const dx=ev.clientX-p.x,dy=ev.clientY-p.y;
                    p.x=ev.clientX;p.y=ev.clientY;
                    if(this.pointers.size===1) {
                        const mousePan=p.type!=='touch'&&p.button===0&&this.dragMode==='pan';
                        if(p.button===2||p.button===1||p.shift||mousePan)this._pan(dx,dy);
                        else this._rotate(dx,dy,ev.timeStamp);
                    } else if(this.pointers.size===2&&this.pinch) {
                        const s=this._pinchState();
                        this.goal.radius=this._clampRadius(this.goal.radius*this.pinch.dist/s.dist);
                        this._pan(s.x-this.pinch.x,s.y-this.pinch.y);
                        this.pinch=s;
                    }
                }
            }
            _up_(e) {
                if(!this.pointers.delete(e.pointerId))return;
                try{this.element.releasePointerCapture(e.pointerId);}catch(_){}
                this.pinch=null;
                if(this.pointers.size===1)this.pinch=null;
                if(this.pointers.size===0) {
                    this.element.classList.remove('dragging');
                    // Release momentum: only if the pointer was still moving when lifted.
                    if(performance.now()-this.lastMove<80) {
                        const k=this.o.inertia*1000;
                        this.goal.theta+=this.velocity.theta*k;
                        this.goal.phi=this._clampPhi(this.goal.phi+this.velocity.phi*k);
                    }
                    this.velocity.theta=this.velocity.phi=0;
                } else if(this.pointers.size===1&&this.pinch===null) {
                    // Remaining finger continues as a rotate without a jump.
                    const rest=[...this.pointers.values()][0];rest.button=0;rest.shift=false;
                }
                this.o.onChange();
            }
            _height() {return this.element.clientHeight||1;}
            _rotate(dx,dy,t) {
                const f=2*Math.PI*this.o.rotateSpeed/this._height();
                const dTheta=-dx*f,dPhi=-dy*f;
                this.goal.theta+=dTheta;
                this.goal.phi=this._clampPhi(this.goal.phi+dPhi);
                const now=performance.now(),dt=Math.max(now-this.lastMove,1);
                this.velocity.theta=0.6*this.velocity.theta+0.4*dTheta/dt;
                this.velocity.phi=0.6*this.velocity.phi+0.4*dPhi/dt;
                this.lastMove=now;
                this.o.onChange();
            }
            _pan(dx,dy) {
                const cam=this.camera,k=2*this.goal.radius*Math.tan(THREE.MathUtils.degToRad(cam.fov)/2)/this._height()*this.o.panSpeed;
                cam.updateMatrixWorld();
                this._right.setFromMatrixColumn(cam.matrixWorld,0);
                this._up.setFromMatrixColumn(cam.matrixWorld,1);
                this.goal.target.addScaledVector(this._right,-dx*k).addScaledVector(this._up,dy*k);
                const L=this.o.panLimit;
                this.goal.target.set(...this.goal.target.toArray().map(v=>Math.max(-L,Math.min(L,v))));
                this.o.onChange();
            }
            _wheel(e) {
                this.flyRate=0;
                e.preventDefault();
                const unit=e.deltaMode===1?16:e.deltaMode===2?100:1;
                const delta=Math.max(-400,Math.min(400,e.deltaY*unit))*(e.ctrlKey?2.5:1);
                this.goal.radius=this._clampRadius(this.goal.radius*Math.exp(delta*this.o.zoomSpeed));
                this.o.onChange();
            }
            _clampPhi(v) {return Math.max(this.o.minPhi,Math.min(this.o.maxPhi,v));}
            _clampRadius(v) {return Math.max(this.o.minDistance,Math.min(this.o.maxDistance,v));}
            // Returns true while the view is still easing toward the goal (caller keeps rendering).
            update(dt) {
                const g=this.goal,v=this.view,rate=this.flyRate||(this.touchActive?this.o.touchDamping:this.o.damping),a=1-Math.exp(-rate*Math.min(dt,0.1));
                const dTheta=g.theta-v.theta,dPhi=g.phi-v.phi,dR=g.radius-v.radius;
                const dT=g.target.distanceTo(v.target);
                if(Math.abs(dTheta)<1e-5&&Math.abs(dPhi)<1e-5&&Math.abs(dR)<1e-3&&dT<1e-3) {
                    this.flyRate=0;
                    if(dTheta||dPhi||dR||dT){v.theta=g.theta;v.phi=g.phi;v.radius=g.radius;v.target.copy(g.target);this._apply();return true;}
                    return false;
                }
                v.theta+=dTheta*a;v.phi+=dPhi*a;v.radius+=dR*a;v.target.lerp(g.target,a);
                this._apply();
                return true;
            }
            _apply() {
                const v=this.view,s=Math.sin(v.phi);
                this.camera.up.set(0,0,1);
                this.camera.position.set(v.target.x+v.radius*s*Math.cos(v.theta),v.target.y+v.radius*s*Math.sin(v.theta),v.target.z+v.radius*Math.cos(v.phi));
                this.camera.lookAt(v.target);
            }
        }

        /* ---------- Three.js scene ---------- */
        // With bloom, the scene is drawn into the composer's own target and the canvas only receives a full-screen copy,
        // so the canvas needs no antialiasing of its own (the composer target does that where it is worth it).
        const BLOOM_AVAILABLE=!!(THREE.EffectComposer&&THREE.UnrealBloomPass);
        const renderer=new THREE.WebGLRenderer({antialias:!BLOOM_AVAILABLE,powerPreference:'high-performance'});
        renderer.setClearColor(0x020206, 1.0); // Deep space dark navy/black.
        orbitPlot.appendChild(renderer.domElement);
        const scene=new THREE.Scene();
        const camera=new THREE.PerspectiveCamera(40,1,1,3000);
        // Bloom (visual choice, like a camera exposure): bright light bleeds softly into its surroundings, then a soft
        // highlight shoulder stops piled-up light clipping to flat white. Falls back to a plain render if the add-ons fail to load.
        let bloomPass=null; // switched off with lite visuals (Settings)
        const composer=(()=>{
            if(!BLOOM_AVAILABLE)return null;
            try{
                const gl2=renderer.capabilities&&renderer.capabilities.isWebGL2;
                const opts={type:THREE.HalfFloatType,format:THREE.RGBAFormat};
                const target=gl2&&THREE.WebGLMultisampleRenderTarget?new THREE.WebGLMultisampleRenderTarget(1,1,opts):new THREE.WebGLRenderTarget(1,1,opts);
                const fx=new THREE.EffectComposer(renderer,target); // MSAA samples are set per resolution level in resizeOrbit
                fx.addPass(new THREE.RenderPass(scene,camera));
                bloomPass=new THREE.UnrealBloomPass(new THREE.Vector2(256,256),0.55,0.45,0.62);fx.addPass(bloomPass);
                fx.addPass(new THREE.ShaderPass({uniforms:{tDiffuse:{value:null}},
                    vertexShader:`varying vec2 vUv;void main(){vUv=uv;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
                    fragmentShader:`uniform sampler2D tDiffuse;varying vec2 vUv;
                        void main(){vec4 c=texture2D(tDiffuse,vUv);vec3 x=c.rgb,k=0.7+0.3*(1.0-exp(-(x-0.7)/0.3));
                            gl_FragColor=vec4(mix(x,k,step(0.7,x)),1.0);}`}));
                return fx;
            }catch(e){console.warn('Bloom unavailable',e);return null;}
        })();
        const distantStarCount=900,distantStarPositions=new Float32Array(distantStarCount*3);
        let starSeed=3198;
        const randomStar=()=>{starSeed=(starSeed*1664525+1013904223)>>>0;return starSeed/4294967296;};
        for(let i=0;i<distantStarCount;i++){
            const y=randomStar()*2-1,r=Math.sqrt(1-y*y),angle=randomStar()*Math.PI*2;
            const distance=820+randomStar()*120;
            distantStarPositions[i*3]=Math.cos(angle)*r*distance;
            distantStarPositions[i*3+1]=y*distance;
            distantStarPositions[i*3+2]=Math.sin(angle)*r*distance;
        }
        const distantStarGeometry=new THREE.BufferGeometry();
        distantStarGeometry.setAttribute('position',new THREE.BufferAttribute(distantStarPositions,3));
        const distantStars=new THREE.Points(distantStarGeometry,new THREE.PointsMaterial({
            color:0x68708a,size:1.1,sizeAttenuation:false,transparent:true,opacity:0.16,depthWrite:false,toneMapped:false
        }));
        distantStars.frustumCulled=false;scene.add(distantStars);
        // Three-quarter view: ~39 degrees above the disk plane, framed to fit the galaxy.
        const GALAXY_VIEW=(()=>{const e=[0.55,-1.0,0.8],r=Math.hypot(...e);return{theta:Math.atan2(e[1],e[0]),phi:Math.acos(e[2]/r),radius:r*40,target:[0,0,0]};})();
        // Every simulated star (bulge + disk) is one point of light: a bright core with soft wings, like a star seen through
        // a telescope, and diffraction spikes on the brightest few (visual choice). The buffers hold the largest possible
        // population and draw only the first starCount; syncStarMeshes() fills colour and brightness, updateStarTransforms() positions.
        const starGeometry=new THREE.BufferGeometry();
        starGeometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(MAX_STARS*3),3).setUsage(THREE.DynamicDrawUsage));
        starGeometry.setAttribute('aColor',new THREE.BufferAttribute(new Float32Array(MAX_STARS*3),3));
        starGeometry.setAttribute('aLum',new THREE.BufferAttribute(new Float32Array(MAX_STARS),1));
        const starUniforms={uScale:{value:600},uDpr:{value:1},uWorld:{value:1.1}};
        const starPoints=new THREE.Points(starGeometry,new THREE.ShaderMaterial({
            uniforms:starUniforms,transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,
            vertexShader:`attribute vec3 aColor;attribute float aLum;uniform float uScale;uniform float uDpr;uniform float uWorld;
                varying vec3 vCol;varying float vInt;varying float vSpike;
                void main(){vCol=aColor;vSpike=step(0.88,aLum);vInt=0.3+0.7*aLum;
                    vec4 mv=modelViewMatrix*vec4(position,1.0);
                    float px=clamp(uWorld*(0.55+0.9*aLum)*uScale/max(-mv.z,1.0),2.2*uDpr,22.0*uDpr);
                    gl_PointSize=px*mix(1.0,2.6,vSpike);gl_Position=projectionMatrix*mv;}`,
            fragmentShader:`varying vec3 vCol;varying float vInt;varying float vSpike;
                void main(){vec2 d=gl_PointCoord*2.0-1.0;float r=length(d);if(r>1.0)discard;
                    vec2 q=d*mix(1.0,2.6,vSpike);float qr=length(q),q2=qr*qr;
                    float core=exp(-q2*22.0);
                    float wing=0.28/(1.0+q2*40.0)*(1.0-smoothstep(0.45,1.0,qr));
                    float spike=vSpike*(exp(-abs(d.x)*70.0)+exp(-abs(d.y)*70.0))*pow(1.0-r,2.0)*0.55;
                    vec3 col=mix(vCol,vec3(1.0),core*0.35);
                    gl_FragColor=vec4(col*(core+wing+spike)*vInt,1.0);}`
        }));
        starPoints.frustumCulled=false;starPoints.renderOrder=1;
        scene.add(starPoints);

        // Star colour and brightness from a stellar mass drawn per star (its colorSeed):
        //   masses follow the Kroupa (2001) IMF slope dN/dM ~ M^-2.3 above 0.5 solar masses,
        //   main-sequence luminosity L ~ M^3.5 and radius R ~ M^0.8, so temperature T ~ (L/R^2)^(1/4) ~ M^0.475 (T_sun = 5772 K),
        //   and the colour is that temperature's blackbody colour. Brightness is shown on a log scale, like magnitudes.
        // Old bulge stars stop at ~1 solar mass (heavier ones have died) plus a few red giants; the disk gets more young,
        // massive blue stars further out, so discs look bluer outwards (de Jong 1996). Simplified, visual only.
        const IMF_SLOPE=2.3,M_LO=0.5,M_HI=8,LOG_L_RANGE=3.5*Math.log10(M_HI/M_LO);
        function blackbodyRGB(T,c){ // Approximate blackbody colour (fit to the CIE colour-matching data, 1000-40000 K)
            const t=T/100,k=v=>Math.min(255,Math.max(0,v))/255;
            if(t<=66)c.setRGB(1,k(99.4708025861*Math.log(t)-161.1195681661),t<=19?0:k(138.5177312231*Math.log(t-10)-305.0447927307));
            else c.setRGB(k(329.698727446*Math.pow(t-60,-0.1332047592)),k(288.1221695283*Math.pow(t-60,-0.0755148492)),1);
            return c;
        }
        function starAppearance(pt,c){ // returns brightness 0..1 and fills c with the star's colour
            if(pt.kind==='bulge'&&pt.colorSeed>0.92){blackbodyRGB(4000+4000*(pt.colorSeed-0.92),c);return 0.62;} // red giant
            const mHi=pt.kind==='bulge'?1.0:1.2+(M_HI-1.2)*Math.min(1,pt.radius/(3*currentAnatomy.diskScale));
            const a=1-IMF_SLOPE,lo=M_LO**a,u=pt.kind==='bulge'?pt.colorSeed/0.92:pt.colorSeed;
            const M=(lo+u*(mHi**a-lo))**(1/a);
            blackbodyRGB(5772*M**0.475,c);
            const grey=(c.r+c.g+c.b)/3,sat=1.7; // colours boosted like an astrophoto (visual choice)
            c.setRGB(Math.max(0,grey+(c.r-grey)*sat),Math.max(0,grey+(c.g-grey)*sat),Math.max(0,grey+(c.b-grey)*sat));
            return Math.min(1,3.5*Math.log10(M/M_LO)/LOG_L_RANGE);
        }

        // Diffuse starlight (visual layer drawn from the real model): the billions of stars too faint to see one by one.
        // The disk glow follows the exponential disk the physics uses, I(R) ~ exp(-R/R_d) (Freeman 1970), with the galaxy's own
        // scale length; the bulge glow is a de Vaucouleurs profile (Hernquist 1990 links it to the bulge model), always facing the camera.
        // Their relative brightness comes from the disk and bulge masses (same mass-to-light ratio assumed).
        const diskGlowUniforms={uRd:{value:3},uRmax:{value:20},uGain:{value:0.3},
            uInner:{value:new THREE.Color(0xffd9a8)},uOuter:{value:new THREE.Color(0xc4d4ff)}};
        const diskGlow=new THREE.Mesh(new THREE.PlaneGeometry(2,2),new THREE.ShaderMaterial({
            uniforms:diskGlowUniforms,transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,side:THREE.DoubleSide,
            vertexShader:`varying vec2 vP;uniform float uRmax;void main(){vP=position.xy*uRmax;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
            fragmentShader:`varying vec2 vP;uniform float uRd;uniform float uRmax;uniform float uGain;uniform vec3 uInner;uniform vec3 uOuter;
                void main(){float R=length(vP);float I=exp(-R/uRd)*(1.0-smoothstep(0.8*uRmax,uRmax,R));
                    vec3 col=mix(uInner,uOuter,smoothstep(0.0,3.0*uRd,R));I=uGain*I/(1.0+0.5*I);
                    gl_FragColor=vec4(col*I,1.0);}`
        }));
        diskGlow.frustumCulled=false;diskGlow.renderOrder=-1;scene.add(diskGlow);
        const bulgeGlowUniforms={uRe:{value:0.5},uAmp:{value:0},uGain:{value:0.3},uColor:{value:new THREE.Color(0xffcc88)}};
        const bulgeGlow=new THREE.Mesh(new THREE.PlaneGeometry(2,2),new THREE.ShaderMaterial({
            uniforms:bulgeGlowUniforms,transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,
            vertexShader:`varying vec2 vP;uniform float uRe;void main(){vP=position.xy*6.0*uRe;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
            fragmentShader:`varying vec2 vP;uniform float uRe;uniform float uAmp;uniform float uGain;uniform vec3 uColor;
                void main(){float R=length(vP);float I=uAmp*exp(-7.669*(pow(max(R,0.03*uRe)/uRe,0.25)-1.0))*(1.0-smoothstep(0.7,1.0,R/(6.0*uRe)));
                    I=uGain*I/(1.0+0.5*I);gl_FragColor=vec4(uColor*I,1.0);}`
        }));
        bulgeGlow.frustumCulled=false;bulgeGlow.renderOrder=-1;scene.add(bulgeGlow);
        function updateDiffuseGlow(anatomy){
            const Rd=anatomy.diskScale;
            diskGlowUniforms.uRd.value=Rd;diskGlowUniforms.uRmax.value=anatomy.rMax;
            diskGlow.scale.set(anatomy.rMax,anatomy.rMax,1);
            const Re=1.8153*anatomy.bulgeRadius; // projected half-light radius of a Hernquist bulge
            const amp=anatomy.bulgeMass>0&&anatomy.diskMass>0?anatomy.bulgeMass*2*Rd*Rd/(7.215*anatomy.diskMass*Re*Re):0;
            bulgeGlowUniforms.uRe.value=Re;bulgeGlowUniforms.uAmp.value=amp;
            bulgeGlow.visible=amp>0;bulgeGlow.scale.set(6*Re,6*Re,1);
        }

        // Black hole (visual choice only, not a physics result): shadow, thin photon ring, an accretion disk
        // that tilts with the camera, the far side of the disk lensed over and under the shadow, Doppler
        // beaming (the side coming towards you is brighter and whiter) and inner gas swirling faster than outer gas.
        const BH_R=GALAXY.blackHoleRadius,BH_EXTENT=3.2;
        const blackHole=new THREE.Mesh(
            new THREE.PlaneGeometry(2*BH_R*BH_EXTENT,2*BH_R*BH_EXTENT),
            new THREE.ShaderMaterial({
                uniforms:{uView:{value:new THREE.Vector3(0,0,1)},uTime:{value:0}},
                transparent:true,depthWrite:false,blending:THREE.CustomBlending,
                blendEquation:THREE.AddEquation,blendSrc:THREE.OneFactor,blendDst:THREE.OneMinusSrcAlphaFactor,
                vertexShader:`varying vec2 vP;void main(){vP=position.xy;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
                fragmentShader:`
                    varying vec2 vP;
                    uniform vec3 uView;   // galaxy spin axis in camera space
                    uniform float uTime;
                    float hash(vec2 q){return fract(sin(dot(q,vec2(127.1,311.7)))*43758.5453);}
                    float noise(vec2 q){vec2 i=floor(q),f=fract(q);f=f*f*(3.0-2.0*f);
                        return mix(mix(hash(i),hash(i+vec2(1,0)),f.x),mix(hash(i+vec2(0,1)),hash(i+vec2(1,1)),f.x),f.y);}
                    float streaks(float r,float phi){
                        // Keplerian-style shear: angular speed ~ r^-1.5, so inner gas laps the outer gas.
                        float a=phi+uTime*pow(r,-1.5);
                        vec2 q=vec2(a*2.5,log(r)*9.0);
                        return 0.55*noise(q)+0.3*noise(q*2.1+3.7)+0.15*noise(q*4.3+9.1);
                    }
                    vec3 heat(float t){ // 0 = cool deep orange, 1 = hot white
                        return mix(mix(vec3(0.75,0.16,0.02),vec3(1.0,0.55,0.12),smoothstep(0.0,0.5,t)),vec3(1.0,0.93,0.80),smoothstep(0.5,1.0,t));
                    }
                    void main(){
                        vec2 p=vP/${BH_R.toFixed(3)};
                        float d=length(p);
                        float cosI=clamp(abs(uView.z),0.0,1.0),sinI=sqrt(1.0-cosI*cosI);
                        float side=uView.z>=0.0?1.0:-1.0;
                        vec2 ax=length(uView.xy)>1e-4?normalize(uView.xy):vec2(0.0,1.0);
                        // q.x along the disk's long axis, q.y along its short axis (positive = far side)
                        vec2 q=vec2(dot(p,vec2(ax.y,-ax.x)),dot(p,ax)*side);
                        const float RIN=1.25,ROUT=3.0;

                        // Direct image of the disk, squashed by the viewing angle, with a little thickness edge-on.
                        float c=max(cosI,0.04);
                        float rD=length(vec2(q.x,q.y/c));
                        float phiD=atan(q.y/c,q.x);
                        float thick=exp(-pow(q.y/(0.05+0.95*cosI*ROUT),2.0));
                        float band=smoothstep(RIN-0.05,RIN+0.12,rD)*(1.0-smoothstep(ROUT-0.7,ROUT,rD));
                        float prof=pow(RIN/max(rD,RIN),1.8);
                        float vlos=sinI*cos(phiD)*side*inversesqrt(max(rD,RIN));   // + = coming towards you
                        float beam=pow(max(1.0+0.75*vlos,0.15),3.0);
                        float diskI=band*prof*(0.45+0.9*streaks(rD,phiD))*beam*mix(1.0,thick,sinI);
                        float diskT=clamp(prof*0.9+0.35*vlos,0.0,1.0);
                        // The far half sits behind the hole, so the shadow hides it.
                        float behind=step(0.0,q.y)*(1.0-smoothstep(0.96,1.02,d));
                        diskI*=1.0-behind;

                        // Lensed image of the far disk, bent over the top (and faintly under the bottom).
                        float ang=atan(q.y,q.x);
                        float lr=1.0+0.42*max(sinI,0.25);
                        float ringW=0.10+0.10*cosI;
                        float lensRing=exp(-pow((d-lr)/ringW,2.0));
                        float over=mix(0.55,1.0,smoothstep(-0.2,0.8,q.y/max(d,1e-3)));
                        float under=0.35*smoothstep(0.2,-0.8,q.y/max(d,1e-3));
                        float lBeam=pow(max(1.0+0.6*sinI*cos(ang)*side,0.2),2.0);
                        float lensI=lensRing*(over+under)*(0.5+0.8*streaks(lr*1.3,ang))*lBeam*mix(0.6,1.0,sinI);

                        // Photon ring: a razor-thin bright circle just outside the shadow.
                        float photon=exp(-pow((d-1.03)/0.022,2.0))*1.6+exp(-pow((d-1.07)/0.07,2.0))*0.35;
                        float shadow=1.0-smoothstep(0.95,1.01,d);
                        float halo=exp(-pow((d-1.2)/0.9,2.0))*0.10;
                        float fade=1.0-smoothstep(2.6,3.15,d);

                        vec3 col=heat(diskT)*diskI*2.4+heat(0.75)*lensI*2.0+vec3(1.0,0.9,0.75)*photon+vec3(0.9,0.35,0.08)*halo;
                        col=(vec3(1.0)-exp(-col*1.3))*fade*(1.0-shadow);   // soft tone map so nothing clips
                        // A soft dark backdrop round the hole so the disk stands out against the bright galactic centre.
                        float backdrop=0.6*(1.0-smoothstep(1.0,2.9,d));
                        float a=max(max(shadow*0.97,backdrop),clamp(max(col.r,max(col.g,col.b)),0.0,1.0)*0.85)*fade;
                        gl_FragColor=vec4(col,a);
                    }
                `
            })
        );
        blackHole.position.set(0,0,0);
        blackHole.renderOrder=3;scene.add(blackHole);
        const bhAxis=new THREE.Vector3(),bhInv=new THREE.Quaternion();
        function updateBlackHole(t){
            blackHole.quaternion.copy(camera.quaternion);
            bhAxis.set(0,0,1).applyQuaternion(bhInv.copy(camera.quaternion).invert());
            blackHole.material.uniforms.uView.value.copy(bhAxis);
            blackHole.material.uniforms.uTime.value=t*0.04;
        }

        // Optional orbit trajectory lines: full paths of an evenly spaced subset of the stars (mostly disk, a few bulge).
        const ORBIT_LINE_STARS=60,ORBIT_LINE_STEP=5;
        let orbitLineIdx=[];
        const orbitSegs=Math.ceil((frameCount-1)/ORBIT_LINE_STEP);
        const orbitLineGeom=new THREE.BufferGeometry();
        orbitLineGeom.setAttribute('position',new THREE.BufferAttribute(new Float32Array(ORBIT_LINE_STARS*orbitSegs*6),3).setUsage(THREE.DynamicDrawUsage));
        const orbitLines=new THREE.LineSegments(orbitLineGeom,new THREE.LineBasicMaterial({color:0x9ff0ff,transparent:true,opacity:0.6,depthWrite:false,blending:THREE.AdditiveBlending}));
        orbitLines.frustumCulled=false;orbitLines.visible=false;orbitLines.renderOrder=0;scene.add(orbitLines);
        // A faint circle in the disk plane at the ejection radius (1.5 x R_max), shown by the Ejection line switch or the tutorial.
        const ejectionRing=new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints(Array.from({length:128},(_,i)=>new THREE.Vector3(Math.cos(i*Math.PI/64),Math.sin(i*Math.PI/64),0))),
            new THREE.LineBasicMaterial({color:0xff8a80,transparent:true,opacity:0.55,depthWrite:false}));
        ejectionRing.frustumCulled=false;ejectionRing.visible=false;scene.add(ejectionRing);
        let ejectionLineShown=false;
        // MOND boost in 3D (MOND mode, "Where MOND acts" switch): the graph's gold shading laid into the disk plane. Each
        // radius glows with the same shade as the graph at that radius (mondZoneLayout hands over mondZoneLevel), so
        // the gold is in proportion to MOND's boost: faint near the centre and deeper outwards; past the last data
        // radius it holds the last shade and fades softly to nothing by MOND_PLANE_EXTENT x that radius. No line or
        // ring: MOND has no switch-on point. It is drawn in the graph's gold
        // (MOND_GOLD, #f2c879) so the two read as one thing. A haze of faint points spread a little
        // above and below the plane, brighter where the shade is deeper, gives it depth like the stars. All a visual
        // choice: the plane and haze carry no physics beyond the graph's shading.
        const MOND_PLANE_SAMPLES=256,MOND_PLANE_EXTENT=1.6,MOND_PLANE_STRENGTH=0.16,MOND_GOLD=[242/255,200/255,121/255];
        const MOND_HAZE_POINTS=3000,MOND_HAZE_THICKNESS=0.012,MOND_HAZE_STRENGTH=0.65;
        let mondZoneShown=false,mondZoneReady=false;
        const mondZone=new THREE.Group();mondZone.visible=false;scene.add(mondZone);
        const mondLevelData=new Uint8Array(MOND_PLANE_SAMPLES*4);
        const mondLevelTexture=new THREE.DataTexture(mondLevelData,MOND_PLANE_SAMPLES,1,THREE.RGBAFormat);
        mondLevelTexture.magFilter=mondLevelTexture.minFilter=THREE.LinearFilter;
        const mondPlaneUniforms={uLevels:{value:mondLevelTexture},uDataEdge:{value:1/MOND_PLANE_EXTENT},
            uStrength:{value:MOND_PLANE_STRENGTH},uGold:{value:new THREE.Vector3(...MOND_GOLD)}};
        const mondPlane=new THREE.Mesh(new THREE.CircleGeometry(1,256),new THREE.ShaderMaterial({
            uniforms:mondPlaneUniforms,transparent:true,depthWrite:false,side:THREE.DoubleSide,blending:THREE.AdditiveBlending,
            vertexShader:`varying float vU;void main(){vU=length(position.xy);gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
            fragmentShader:`varying float vU;uniform sampler2D uLevels;uniform float uDataEdge;uniform float uStrength;uniform vec3 uGold;
                void main(){float level=texture2D(uLevels,vec2(min(vU,uDataEdge)/uDataEdge*(1.0-1.0/${MOND_PLANE_SAMPLES}.0)+0.5/${MOND_PLANE_SAMPLES}.0,0.5)).r;
                float fade=1.0-smoothstep(uDataEdge,1.0,vU);
                gl_FragColor=vec4(mix(uGold,vec3(1.0),0.3),uStrength*level*fade);}`
        }));
        mondPlane.frustumCulled=false;mondPlane.renderOrder=1;mondZone.add(mondPlane);
        // Haze: fixed points spread evenly over the plane's area; their brightness is set from the shade where they sit.
        const mondHazeU=new Float32Array(MOND_HAZE_POINTS),mondHazeGeom=new THREE.BufferGeometry();
        (()=>{const rng=mulberry32(5310),pos=new Float32Array(MOND_HAZE_POINTS*3),tw=new Float32Array(MOND_HAZE_POINTS);
            for(let i=0;i<MOND_HAZE_POINTS;i++){
                const u=Math.sqrt(rng()),phi=2*Math.PI*rng();mondHazeU[i]=u;
                pos[3*i]=u*Math.cos(phi);pos[3*i+1]=u*Math.sin(phi);pos[3*i+2]=MOND_HAZE_THICKNESS*boxMuller(0,1,rng);
                tw[i]=0.4+0.6*rng();
            }
            mondHazeGeom.setAttribute('position',new THREE.BufferAttribute(pos,3));
            mondHazeGeom.setAttribute('aTwinkle',new THREE.BufferAttribute(tw,1));
            mondHazeGeom.setAttribute('aAlpha',new THREE.BufferAttribute(new Float32Array(MOND_HAZE_POINTS),1));})();
        const mondHazeUniforms={uSize:{value:1},uScale:{value:600}};
        const mondHaze=new THREE.Points(mondHazeGeom,new THREE.ShaderMaterial({
            uniforms:mondHazeUniforms,transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,
            vertexShader:`attribute float aAlpha;attribute float aTwinkle;uniform float uSize;uniform float uScale;varying float vA;
                void main(){vA=aAlpha*aTwinkle;vec4 mv=modelViewMatrix*vec4(position,1.0);
                gl_PointSize=clamp(uSize*uScale/max(-mv.z,1.0),1.5,9.0);gl_Position=projectionMatrix*mv;}`,
            fragmentShader:`varying float vA;void main(){vec2 d=gl_PointCoord-0.5;float r2=dot(d,d)*4.0;if(r2>1.0)discard;
                gl_FragColor=vec4(${MOND_GOLD.join(',')},vA*exp(-3.0*r2));}`
        }));
        mondHaze.frustumCulled=false;mondHaze.renderOrder=1;mondZone.add(mondHaze);
        const mondZoneToggle=document.getElementById('toggleMondZone');
        // While the gold shows, the panels over the 3D view turn solid (CSS .mond-zone-on) so it never tints them.
        const showMondZone=()=>{
            mondZone.visible=mondMode&&mondZoneShown&&mondZoneReady;
            orbitPlot.classList.toggle('mond-zone-on',mondZone.visible);
            requestRender();
        };
        // levelAt(rKpc): the graph's shade at that radius, 0 to 1; lastKpc: the last data radius. levelAt null clears
        // it (NFW mode).
        function setMondZone3D(levelAt,lastKpc){
            mondZoneToggle.hidden=!mondMode; // each model shows only its own controls (no greyed-out ones)
            mondZoneReady=!!levelAt&&lastKpc>0;
            if(mondZoneReady){
                const extent=MOND_PLANE_EXTENT*lastKpc,levels=new Float32Array(MOND_PLANE_SAMPLES);
                for(let i=0;i<MOND_PLANE_SAMPLES;i++){levels[i]=levelAt(lastKpc*i/(MOND_PLANE_SAMPLES-1));mondLevelData[4*i]=Math.round(255*levels[i]);}
                mondLevelTexture.needsUpdate=true;
                mondZone.scale.set(extent,extent,extent);
                const alpha=mondHazeGeom.attributes.aAlpha.array,dataEdge=1/MOND_PLANE_EXTENT;
                for(let i=0;i<MOND_HAZE_POINTS;i++){
                    const u=mondHazeU[i],k=Math.min(1,u/dataEdge)*(MOND_PLANE_SAMPLES-1),j=Math.min(Math.floor(k),MOND_PLANE_SAMPLES-2);
                    alpha[i]=MOND_HAZE_STRENGTH*(levels[j]+(levels[j+1]-levels[j])*(k-j))*(1-smoothStep(dataEdge,1,u));
                }
                mondHazeGeom.attributes.aAlpha.needsUpdate=true;
                mondHazeUniforms.uSize.value=0.012*extent;
            }
            showMondZone();
        }
        mondZoneToggle.addEventListener('click',()=>{
            if(!mondMode)return;
            mondZoneShown=!mondZoneShown;
            setSwitch(mondZoneToggle,mondZoneShown);
            showMondZone();
        });

        // Cold gas: a faint, thicker, more extended exponential cloud (PhysicsEngine.GAS_SCALE_RATIO x the stellar scale
        // length, cut at R_max) drawn as soft points; its brightness follows the gas mass.
        const GAS_POINTS=2400,GAS_SCALE_HEIGHT=0.35,GAS_REFERENCE_MASS=2e10;
        const gasSeeds=(()=>{
            const rng=mulberry32(2403),uniform=new Float32Array(GAS_POINTS*2),gauss=new Float32Array(GAS_POINTS*3);
            for(let i=0;i<GAS_POINTS*2;i++)uniform[i]=rng();
            for(let i=0;i<GAS_POINTS*3;i++)gauss[i]=boxMuller(0,1,rng);
            return{uniform,gauss};
        })();
        function createSoftTexture(){
            const canvas=document.createElement('canvas');
            canvas.width=canvas.height=64;
            const ctx=canvas.getContext('2d'),grad=ctx.createRadialGradient(32,32,0,32,32,32);
            grad.addColorStop(0,'rgba(255,255,255,0.9)');grad.addColorStop(0.45,'rgba(255,255,255,0.3)');grad.addColorStop(1,'rgba(255,255,255,0)');
            ctx.fillStyle=grad;ctx.fillRect(0,0,64,64);
            return new THREE.CanvasTexture(canvas);
        }
        const gasGeometry=new THREE.BufferGeometry();
        gasGeometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(GAS_POINTS*3),3).setUsage(THREE.DynamicDrawUsage));
        const gasMaterial=new THREE.PointsMaterial({map:createSoftTexture(),color:0x5fc8e6,size:1.6,sizeAttenuation:true,transparent:true,opacity:0.2,blending:THREE.AdditiveBlending,depthWrite:false});
        const gasCloud=new THREE.Points(gasGeometry,gasMaterial);
        gasCloud.frustumCulled=false;gasCloud.renderOrder=0;scene.add(gasCloud);
        let gasBaseOpacity=0.2; // set from the gas mass; the render loop dims it when zoomed in (visual choice)
        // Visual choice: gas is dim close up (star view), brightest at the reset view, easing back to normal far out.
        const smoothStep=(a,b,x)=>{const t=Math.min(1,Math.max(0,(x-a)/(b-a)));return t*t*(3-2*t);};
        const GAS_RESET_RADIUS=72,GAS_RESET_BOOST=1.22; // the reset view is the reference look; zooming either way works off it
        const gasZoomOpacity=r=>gasBaseOpacity*(r<=GAS_RESET_RADIUS
            ?0.15+(GAS_RESET_BOOST-0.15)*smoothStep(28,GAS_RESET_RADIUS,r)
            :1+(GAS_RESET_BOOST-1)*(1-smoothStep(GAS_RESET_RADIUS,200,r)));
        function fillGasCloud(anatomy){
            const a=gasGeometry.attributes.position.array,scale=PhysicsEngine.GAS_SCALE_RATIO*anatomy.diskScale;
            for(let i=0;i<GAS_POINTS;i++){
                const R=diskRadius(gasSeeds.uniform[2*i],scale,0,anatomy.rMax),phi=2*Math.PI*gasSeeds.uniform[2*i+1];
                a[3*i]=R*Math.cos(phi);a[3*i+1]=R*Math.sin(phi);a[3*i+2]=gasSeeds.gauss[3*i]*GAS_SCALE_HEIGHT;
            }
            gasGeometry.attributes.position.needsUpdate=true;
            gasBaseOpacity=0.1+0.2*Math.min(1,Math.sqrt(anatomy.gasMass/GAS_REFERENCE_MASS));
        }

        // Refills the star meshes (position, size, colour) and the orbit-line selection from the current population.
        function syncStarMeshes(){
            fillGasCloud(currentAnatomy);
            // Camera framing for Reset / Face-on / Edge-on: NGC 3198 (extent 22 kpc) keeps its original distance.
            const extent=Math.min(55,Math.max(11,Math.max(6*currentAnatomy.diskScale,0.5*currentAnatomy.rMax)));
            GALAXY_VIEW.radius=Math.hypot(0.55,-1.0,0.8)*40*extent/22;
            updateDiffuseGlow(currentAnatomy);
            starGeometry.setDrawRange(0,starCount);
            const c=new THREE.Color(),pos=starGeometry.attributes.position.array,
                cols=starGeometry.attributes.aColor.array,lums=starGeometry.attributes.aLum.array;
            for(let k=0;k<starCount;k++){
                const pt=particles[k],p=pt.position;
                lums[k]=starAppearance(pt,c);
                cols[k*3]=c.r;cols[k*3+1]=c.g;cols[k*3+2]=c.b;
                pos[k*3]=p[0];pos[k*3+1]=p[1];pos[k*3+2]=p[2];
            }
            starGeometry.attributes.position.needsUpdate=true;
            starGeometry.attributes.aColor.needsUpdate=true;starGeometry.attributes.aLum.needsUpdate=true;
            const diskStars=starCount-bulgeStarCount;
            const bulgeLines=bulgeStarCount>0?Math.min(bulgeStarCount,Math.max(6,Math.round(ORBIT_LINE_STARS*bulgeStarCount/starCount))):0;
            const diskLines=Math.min(diskStars,ORBIT_LINE_STARS-bulgeLines);
            orbitLineIdx=[];
            const bulgeIdx=[],diskIdx=[];
            for(let k=0;k<starCount;k++)(particles[k].kind==='bulge'?bulgeIdx:diskIdx).push(k);
            for(let j=0;j<bulgeLines;j++)orbitLineIdx.push(bulgeIdx[Math.floor(j*bulgeStarCount/bulgeLines)]);
            for(let j=0;j<diskLines;j++)orbitLineIdx.push(diskIdx[Math.floor(j*diskStars/diskLines)]);
            orbitLineGeom.setDrawRange(0,orbitLineIdx.length*orbitSegs*2);
        }

        /* ---------- NFW dark matter halo visual ----------
           Points are drawn from the true NFW enclosed-mass profile M(<r) ~ ln(1+x)-x/(1+x), x=r/rs,
           truncated at r200. Fixed uniform samples are re-mapped whenever M200 or rs changes, so the
           cloud's size (r200 ~ M200^(1/3)) and central concentration (c=r200/rs) track the sliders. */
        const HALO_POINTS=16000,HALO_TABLE=512;
        const haloGroup=new THREE.Group();haloGroup.visible=false;scene.add(haloGroup);
        const haloU=new Float32Array(HALO_POINTS),haloDirs=new Float32Array(HALO_POINTS*3);
        (()=>{const rng=mulberry32(8200);for(let i=0;i<HALO_POINTS;i++){
            haloU[i]=rng();const cz=2*rng()-1,sz=Math.sqrt(1-cz*cz),ph=2*Math.PI*rng();
            haloDirs[i*3]=sz*Math.cos(ph);haloDirs[i*3+1]=sz*Math.sin(ph);haloDirs[i*3+2]=cz;}})();
        const haloGeom=new THREE.BufferGeometry();
        haloGeom.setAttribute('position',new THREE.BufferAttribute(new Float32Array(HALO_POINTS*3),3));
        haloGeom.setAttribute('aFrac',new THREE.BufferAttribute(new Float32Array(HALO_POINTS),1));
        const haloUniforms={uScale:{value:600},uAlpha:{value:0.2},uSize:{value:2.5}};
        const haloPoints=new THREE.Points(haloGeom,new THREE.ShaderMaterial({
            uniforms:haloUniforms,transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,
            vertexShader:`            attribute float aFrac;uniform float uScale;uniform float uSize;varying float vFrac;varying float vFade;
                            void main(){vFrac=aFrac;vec4 mv=modelViewMatrix*vec4(position,1.0);vFade=0.08+0.92*smoothstep(40.0,200.0,-mv.z);
                            gl_PointSize=clamp(uSize*uScale/max(-mv.z,1.0),2.0,10.0);gl_Position=projectionMatrix*mv;}`,
                        fragmentShader:`varying float vFrac;varying float vFade;uniform float uAlpha;
                void main(){vec2 d=gl_PointCoord-0.5;float r2=dot(d,d)*4.0;if(r2>1.0)discard;
                vec3 inner=vec3(0.55,0.35,1.0),outer=vec3(0.15,0.35,0.9);
                float a=exp(-3.0*r2)*uAlpha*vFade*mix(1.0,0.55,vFrac);
                gl_FragColor=vec4(mix(inner,outer,sqrt(vFrac)),a);}`
        }));
        haloPoints.frustumCulled=false;haloPoints.renderOrder=1;haloGroup.add(haloPoints);
        // Faint rim marking the virial radius r200.
        const haloShell=new THREE.Mesh(new THREE.IcosahedronGeometry(1,4),new THREE.ShaderMaterial({
            uniforms:{uFade:{value:1}},
            transparent:true,depthWrite:false,side:THREE.DoubleSide,blending:THREE.AdditiveBlending,
            vertexShader:`varying vec3 vN;varying vec3 vV;void main(){vec4 mv=modelViewMatrix*vec4(position,1.0);
                vN=normalize(normalMatrix*normal);vV=normalize(-mv.xyz);gl_Position=projectionMatrix*mv;}`,
            fragmentShader:`varying vec3 vN;varying vec3 vV;uniform float uFade;void main(){float f=1.0-abs(dot(normalize(vN),normalize(vV)));
                gl_FragColor=vec4(0.3,0.45,1.0,uFade*(0.05+0.4*pow(f,2.5)));}`
        }));
        haloShell.frustumCulled=false;haloShell.renderOrder=1;haloGroup.add(haloShell);
        // Geodesic (icosahedral) grid: uniform cells with no pole pinching.
                const haloGrid=new THREE.LineSegments(new THREE.WireframeGeometry(new THREE.IcosahedronGeometry(1,2)),
            new THREE.LineBasicMaterial({color:0x5a7dff,transparent:true,opacity:0.13,depthWrite:false,blending:THREE.AdditiveBlending}));
        haloGrid.frustumCulled=false;haloGrid.renderOrder=1;haloGroup.add(haloGrid);
        let haloR200=158;
        const haloToggle=document.getElementById('toggleHalo');
        const haloMetrics=document.getElementById('haloMetrics'),haloEmpty=document.getElementById('haloEmpty');
        let haloHidden=false;
        const fullHaloBtn=document.getElementById('viewFullHalo');
        fullHaloBtn.addEventListener('click',()=>{
            if(haloHidden)haloToggle.click();
            controls.flyTo({theta:GALAXY_VIEW.theta,phi:1.2,radius:haloR200*4.2,target:[0,0,0],rate:2.4});
        });
        const nfwShape=x=>Math.log1p(x)-x/(1+x);
        function updateHaloVisual(massT,rs){
            const show=!mondMode&&!haloHidden&&massT>0;
            document.getElementById('tm-r200').textContent=mondMode||massT<=0?'no halo':`${Math.cbrt(3*massT*1e12/(4*Math.PI*200*(3*H0*H0/(8*Math.PI*G)))).toFixed(0)} kpc`;
            haloGroup.visible=show;
            // MOND has no halo, so every halo control and readout goes (and Full halo goes at zero halo mass too).
            haloToggle.hidden=mondMode;fullHaloBtn.hidden=mondMode||massT<=0;
            document.getElementById('haloInfo').hidden=document.getElementById('tmHaloRow').hidden=mondMode;
            haloMetrics.hidden=mondMode||massT<=0;
            haloEmpty.hidden=!haloMetrics.hidden;
            haloEmpty.textContent='Set the halo mass above zero to see its size and shape.';
            if(mondMode||massT<=0){requestRender();return;}
            const rhoCrit=3*H0*H0/(8*Math.PI*G);
            const r200=Math.cbrt(3*massT*1e12/(4*Math.PI*200*rhoCrit)),c=r200/rs,mc=nfwShape(c);
            // Invert the enclosed-mass fraction -> x table by bisection.
            const table=new Float32Array(HALO_TABLE+1);
            for(let k=1;k<=HALO_TABLE;k++){
                const target=k/HALO_TABLE*mc;let lo=0,hi=c;
                for(let it=0;it<40;it++){const mid=(lo+hi)/2;if(nfwShape(mid)<target)lo=mid;else hi=mid;}
                table[k]=(lo+hi)/2*rs;
            }
            const pos=haloGeom.attributes.position.array,frac=haloGeom.attributes.aFrac.array;
            for(let i=0;i<HALO_POINTS;i++){
                const t=haloU[i]*HALO_TABLE,k=Math.min(Math.floor(t),HALO_TABLE-1),r=table[k]+(t-k)*(table[k+1]-table[k]);
                pos[i*3]=haloDirs[i*3]*r;pos[i*3+1]=haloDirs[i*3+1]*r;pos[i*3+2]=haloDirs[i*3+2]*r;frac[i]=r/r200;
            }
            haloGeom.attributes.position.needsUpdate=true;haloGeom.attributes.aFrac.needsUpdate=true;
            haloShell.scale.setScalar(r200);haloGrid.scale.setScalar(r200);haloR200=r200;
            haloUniforms.uAlpha.value=0.05+0.28*Math.min(1,massT/2);
            document.getElementById('haloRadius').textContent=`${r200.toFixed(0)} kpc`;
            document.getElementById('haloScale').textContent=`${rs.toFixed(1)} kpc`;
            document.getElementById('haloConcentration').textContent=c.toFixed(1);
            requestRender();
        }
        haloToggle.addEventListener('click',()=>{
            haloHidden=!haloHidden;
            setSwitch(haloToggle,!haloHidden);
            updateHaloVisual(Number(massSlider.value),Number(scaleSlider.value));
        });
        function fillOrbitLines(){
            const a=orbitLineGeom.attributes.position.array;let o=0;
            for(const idx of orbitLineIdx){
                const t=activeTrajectories[idx];
                for(let s=0;s<orbitSegs;s++){
                    const f0=s*ORBIT_LINE_STEP;
                    const f1=Math.min((s+1)*ORBIT_LINE_STEP,frameCount-1);
                    const i0=f0*3,i1=f1*3;
                    a[o++]=t[i0];a[o++]=t[i0+1];a[o++]=t[i0+2];
                    a[o++]=t[i1];a[o++]=t[i1+1];a[o++]=t[i1+2];
                }
            }
            orbitLineGeom.attributes.position.needsUpdate=true;
        }
        syncStarMeshes();
        fillOrbitLines();

        // A star is ejected once it passes the ejection radius, PM.ejection_factor x R_max of the active galaxy.
        const escapeCountVal=document.getElementById('escape-count-val');
        let lastEscaped=-1;
        // The corner readout follows the active population: its star total and where the ejection radius lies.
        function updateStarTelemetry(){
            document.getElementById('star-total-val').textContent=starCount;
            document.getElementById('starCountStat').textContent=starCount;
            // The star locator's range follows the population, so a bigger galaxy offers more stars to find.
            ['trackInput','fsTrack'].forEach(id=>{const el=document.getElementById(id);if(el)el.max=starCount;});
            const trackMaxEl=document.getElementById('trackMax');if(trackMaxEl)trackMaxEl.textContent=`1–${starCount}`;const trackLastEl=document.getElementById('trackLast');if(trackLastEl)trackLastEl.textContent=starCount;
            if(typeof ejectionRing!=='undefined')ejectionRing.scale.setScalar(ejectionRadius);
            document.getElementById('escape-counter').title=`Stars farther than ${PM.ejection_factor} × R_max = ${ejectionRadius.toFixed(1)} kpc from the centre count as ejected`;
            lastEscaped=-1;
        }
        updateStarTelemetry();

        // Synchronized escape counter using bisection on exact particle ejection times
        function updateEscapeCounter(timeMyr){
            const t=Math.max(0,Math.min(SIM_TOTAL_TIME_MYR,timeMyr));
            let low=0,high=sortedEscapeTimes.length;
            while(low<high){
                const mid=(low+high)>>1;
                if(sortedEscapeTimes[mid]<=t)low=mid+1;
                else high=mid;
            }
            const escaped=low;
            if(escaped!==lastEscaped){lastEscaped=escaped;escapeCountVal.textContent=escaped;}
        }

        // Sub-frame continuous orbit updates for locked 60 FPS animation
        function updateStarTransforms(timeMyr){
            const fFloat=Math.max(0,Math.min(frameCount-1,timeMyr/SIM_DT_MYR));
            const f0=Math.floor(fFloat);
            const f1=Math.min(f0+1,frameCount-1);
            const alpha=fFloat-f0;
            const posArray=starGeometry.attributes.position.array;
            const i0=f0*3,i1=f1*3;
            for(let k=0;k<starCount;k++){
                const t=activeTrajectories[k];
                posArray[k*3]=t[i0]+(t[i1]-t[i0])*alpha;
                posArray[k*3+1]=t[i0+1]+(t[i1+1]-t[i0+1])*alpha;
                posArray[k*3+2]=t[i0+2]+(t[i1+2]-t[i0+2])*alpha;
            }
            starGeometry.attributes.position.needsUpdate=true;
        }

        // Render loop: continuous delta-time playback spanning 0 to 1000 Myr (1 Gyr)
        let needsRender=true,lastTime=performance.now(),simTime=0.0,playing=false,lastReportedFrame=-1,playbackMultiplier=1.0;
        const PLAYBACK_SPEED_MYR_PER_SEC=25.0; // half of the nominal 50 Myr/s; UI labels stay 0.25x-2x
        const controls=new SmoothOrbitControls(camera,orbitPlot,{initial:GALAXY_VIEW,onChange:()=>{needsRender=true;}});
        function requestRender(){needsRender=true;}
        // Panel clicks must not reach the controls' drag / wheel / double-click handlers.
        const viewPanel=document.getElementById('viewPanel');
        function setSwitch(btn,on){btn.setAttribute('aria-checked',String(on));btn.setAttribute('aria-pressed',String(on));}
        ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu'].forEach(type=>viewPanel.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
        const viewToggle=document.getElementById('viewToggle');
        viewPanel.classList.add('collapsed');
        viewToggle.setAttribute('aria-expanded','false');
        viewToggle.addEventListener('click',()=>{
            const collapsed=viewPanel.classList.toggle('collapsed');
            viewToggle.setAttribute('aria-expanded',String(!collapsed));
        });
        // The "?" tips open on hover with a mouse; on a touch screen a tap opens one (and closes the others).
        document.addEventListener('click',e=>{
            const tip=e.target.closest&&e.target.closest('.stat-q');
            document.querySelectorAll('.stat-q.open').forEach(q=>{if(q!==tip)q.classList.remove('open');});
            if(tip)tip.classList.toggle('open');
        });
        document.querySelectorAll('.halo-metric').forEach(metric=>{
            const button=metric.querySelector('button');
            // A mouse has to rest on the card for a moment before it opens, so sweeping past does not open it.
            let openTimer=null;
            metric.addEventListener('pointerenter',event=>{
                if(event.pointerType==='touch')return;
                clearTimeout(openTimer);openTimer=setTimeout(()=>{metric.classList.add('hovered');sync();},HOVER_OPEN_DELAY_MS);
            });
            metric.addEventListener('pointerleave',event=>{
                if(event.pointerType==='touch')return;
                clearTimeout(openTimer);metric.classList.remove('hovered');sync();
            });
            // A mouse only hovers (a click does not pin the card open); touch and keyboard toggle it instead, one card at a time.
            const sync=()=>{
                const open=metric.classList.contains('hovered')||metric.classList.contains('expanded');
                button.setAttribute('aria-expanded',String(open));
                trackExplainerRead(metric,open);
            };
            let lastPointer='';
            button.addEventListener('pointerdown',event=>{lastPointer=event.pointerType;});
            button.addEventListener('click',()=>{
                const mouse=lastPointer==='mouse';lastPointer='';
                if(mouse)return;
                // Opening one card closes the others, so only one explainer is open at a time.
                document.querySelectorAll('.halo-metric.expanded').forEach(other=>{
                    if(other===metric)return;
                    other.classList.remove('expanded');
                    const open=other.classList.contains('hovered');
                    other.querySelector('button').setAttribute('aria-expanded',String(open));
                    trackExplainerRead(other,open);
                });
                metric.classList.toggle('expanded');sync();
            });
        });
        function resetView(){controls.reset();requestRender();}
        document.getElementById('viewReset').addEventListener('click',resetView);
        const nativeFsElement=()=>document.fullscreenElement||document.webkitFullscreenElement||null;
        function setCssFullscreen(on){
            orbitPlot.classList.toggle('fullscreen-active',on);
            document.body.classList.toggle('fullscreen-lock',on);
            if(fsBtn)fsBtn.setAttribute('aria-pressed',String(on));
            syncTutorialFullscreen();
            resizeOrbit();requestAnimationFrame(resizeOrbit);
        }
        function toggleFullscreen() {
            const container = orbitPlot;
            if (nativeFsElement()) {
                (document.exitFullscreen||document.webkitExitFullscreen).call(document);
                return;
            }
            if (container.classList.contains('fullscreen-active')) { setCssFullscreen(false); return; }
            const request=container.requestFullscreen||container.webkitRequestFullscreen;
            if (!request) { setCssFullscreen(true); return; }
            try {
                const result=request.call(container);
                if (result&&result.catch) result.catch(()=>setCssFullscreen(true));
            } catch (err) { setCssFullscreen(true); }
            // iOS Safari silently ignores the request on non-video elements.
            setTimeout(()=>{ if(!nativeFsElement()&&!container.classList.contains('fullscreen-active'))setCssFullscreen(true); },250);
        }
        const fsBtn = document.getElementById('toggleFullscreen');
        fsBtn.addEventListener('click', toggleFullscreen);
        // Builds the Time vs. Escaped Count log from the exact escape times of the active run (independent of scrubbing).
        const LOG_INTERVAL_MYR=5;
        function buildSimulationLog(){
            const log=[];
            let n=0;
            for(let t=0;t<=SIM_TOTAL_TIME_MYR+1e-9;t+=LOG_INTERVAL_MYR){
                while(n<sortedEscapeTimes.length&&sortedEscapeTimes[n]<=t)n++;
                log.push({time:t,escaped:n,total:particles.length,haloMass:activeHaloMassT});
            }
            return log;
        }
        function interpolateAt(values,radius){
            const last=radii.length-1;
            if(radius<=radii[0])return values[0];
            if(radius>=radii[last])return values[last];
            let low=0,high=last;
            while(high-low>1){const m=(low+high)>>1;if(radii[m]<=radius)low=m;else high=m;}
            return values[low]+(radius-radii[low])/(radii[high]-radii[low])*(values[high]-values[low]);
        }
        // Radial profile on a 1 kpc grid, limited to the range covered by the active galaxy's SPARC data.
        function buildRotationCurveCsv(){
            const gas=customActive?componentCurves[2].slice():rotationRows.map(p=>p.gas_kms),gasSq=customActive?componentCurves[2].map(v=>v*v):rotationRows.map(p=>p.gas_sq);
            const disk=customActive?componentCurves[0].map((v,i)=>Math.hypot(v,componentCurves[1][i])):rotationRows.map(p=>Math.sqrt(mlState.disk*p.disk_kms**2+mlState.bulge*p.bulge_kms**2)),err=rotationRows.map(p=>p.error_kms);
            const rs=Number(scaleSlider.value),halo=extraVelocityModel(activeHaloMassT*1e12,rs,r=>interpolateAt(vBary,r));
            const f=v=>v.toFixed(3),frameScale=Math.sqrt(fitFrame.f)/fitFrame.s; // a best fit's distance and tilt (1 otherwise)
            // The sandbox galaxy has no measurements, so its file leaves out the observed and residual columns.
            let csv=customActive?"Radius_kpc,V_gas,V_stars,V_bar,V_nfw_halo,V_total_model\n":"Radius_kpc,V_obs,V_obs_err,V_gas,V_stars,V_bar,V_nfw_halo,V_total_model,Residual_V_obs_minus_V_bar\n";
            for(let r=Math.max(1,Math.ceil(radii[0]-1e-9));r<=radii[radii.length-1]+1e-9;r++){
                const vo=interpolateAt(vObs,r),eo=interpolateAt(err,r),vg=interpolateAt(gas,r)*frameScale,vd=interpolateAt(disk,r)*frameScale;
                const vb=Math.sqrt(Math.max(0,interpolateAt(gasSq,r)*frameScale**2+vd*vd)),vh=mondMode?Math.sqrt(Math.max(mondModelVelocity(vb,r)**2-vb*vb,0)):halo(r);
                const row=customActive?[r,f(vg),f(vd),f(vb),f(vh),f(Math.hypot(vb,vh))]:[r,f(vo),f(eo),f(vg),f(vd),f(vb),f(vh),f(Math.hypot(vb,vh)),f(vo-vb)];
                csv+=row.join(',')+"\n";
            }
            return csv;
        }
        function downloadCsv(filename,content){
            const link=document.createElement('a');
            link.href='data:text/csv;charset=utf-8,'+encodeURIComponent(content);
            link.download=filename;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
        function exportSimulationData(){
            let csvContent="Time_Myr,Escaped_Count,Total_Particles,Halo_Mass_M200\n";
            buildSimulationLog().forEach(row=>{
                csvContent+=`${row.time},${row.escaped},${row.total},${(row.haloMass*1e12).toExponential(4)}\n`;
            });
            const filePrefix=customActive?'custom_galaxy':activeGalaxy.id.toLowerCase();
            downloadCsv(`${filePrefix}_timeseries_log.csv`,csvContent);
            // Browsers may block back-to-back downloads, so stagger the second file.
            setTimeout(()=>downloadCsv(`${filePrefix}_rotation_curve_model.csv`,buildRotationCurveCsv()),300);
        }
        document.getElementById('exportCsv').addEventListener('click',exportSimulationData);
        const orbitToggle=document.getElementById('toggleOrbits');
        orbitToggle.addEventListener('click',()=>{
            orbitLines.visible=!orbitLines.visible;
            setSwitch(orbitToggle,orbitLines.visible);
            requestRender();
        });
        const gasToggle=document.getElementById('toggleGas');
        gasToggle.addEventListener('click',()=>{
            gasCloud.visible=!gasCloud.visible;
            setSwitch(gasToggle,gasCloud.visible);
            requestRender();
        });
        // Ejection line: the red ring at PM.ejection_factor x R_max. Off by default; the ejection tour cards ask for it.
        const ejectionToggle=document.getElementById('toggleEjection');
        ejectionToggle.addEventListener('click',()=>{
            ejectionLineShown=!ejectionLineShown;
            setSwitch(ejectionToggle,ejectionLineShown);
            ejectionRing.scale.setScalar(ejectionRadius);ejectionRing.visible=ejectionLineShown;
            requestRender();
        });
        document.getElementById('viewFaceOn').addEventListener('click',()=>controls.flyTo({phi:0.04,target:[0,0,0],radius:GALAXY_VIEW.radius,rate:3}));
        document.getElementById('viewEdgeOn').addEventListener('click',()=>controls.flyTo({phi:Math.PI/2,target:[0,0,0],radius:GALAXY_VIEW.radius,rate:3}));
        // Rotate/Pan only changes what a primary mouse drag does; touch gestures are unaffected.
        const modeButtons=[...viewPanel.querySelectorAll('[data-mode]')];
        function setDragMode(mode){
            controls.dragMode=mode;
            orbitPlot.classList.toggle('pan-mode',mode==='pan');
            modeButtons.forEach(b=>b.setAttribute('aria-checked',String(b.dataset.mode===mode)));
        }
        modeButtons.forEach(b=>b.addEventListener('click',()=>setDragMode(b.dataset.mode)));
        // Playback speed multiplier (0.25x, 0.5x, 1x, 2x)
        const speedButtons=[...viewPanel.querySelectorAll('[data-speed]')];
        function setPlaybackMultiplier(speed){
            playbackMultiplier=parseFloat(speed)||1.0;
            speedButtons.forEach(b=>b.setAttribute('aria-checked',String(parseFloat(b.dataset.speed)===playbackMultiplier)));
        }
        speedButtons.forEach((b,idx)=>{
            b.addEventListener('click',()=>setPlaybackMultiplier(b.dataset.speed));
            b.addEventListener('keydown',e=>{
                let nextIdx=-1;
                if(e.key==='ArrowRight'||e.key==='ArrowDown')nextIdx=(idx+1)%speedButtons.length;
                else if(e.key==='ArrowLeft'||e.key==='ArrowUp')nextIdx=(idx-1+speedButtons.length)%speedButtons.length;
                if(nextIdx>=0){
                    e.preventDefault();
                    speedButtons[nextIdx].focus();
                    setPlaybackMultiplier(speedButtons[nextIdx].dataset.speed);
                }
            });
        });
        let containerOriginalWidth = orbitPlot.clientWidth || 800;
        let containerOriginalHeight = orbitPlot.clientHeight || 600;
        /* Adaptive resolution (performance only; the physics is untouched). The GPU cost of the 3D view grows with its pixel
           count, so a large or high-DPI screen (fullscreen especially) can drop well below 60 fps. While the view redraws every
           frame, resolution levels step down when most frames run slow and back up after a spell of smooth frames; a level that
           proves too slow waits longer each time before it is tried again, so the resolution does not flicker between two.
           Level 0 is the screen's own pixel ratio (capped at 2), with 4x MSAA only on low-DPI screens, where edges show most. */
        const RES_SCALES=[1,0.85,0.7,0.6,0.5],RES_MIN_PIXEL_RATIO=0.6;
        const RES_WINDOW_FRAMES=20,RES_SLOW_MS=1000/45,RES_FAST_MS=1000/55,RES_UPGRADE_WINDOWS=8,RES_IGNORE_MS=250;
        let resLevels=[],resLevel=0,resFrames=0,resSlow=0,resFastWindows=0,resUpgradeWindows=RES_UPGRADE_WINDOWS,resLastDowngradeFrom=-1;
        /* Settings' "High-quality visuals" tick box. Until the visitor uses it, quality is automatic (the adaptive steps above)
           and the box shows whether the view is at full resolution. Ticking fixes level 0; unticking fixes the lowest level
           and lite visuals for the whole page (no glow, pulsing, blur or motion; no bloom in the 3D view). In automatic
           mode, lite visuals also come on, for the rest of the visit, once the resolution has had to drop to its lowest
           level. The choice is remembered in this browser (localStorage) so it can be set once before a presentation. */
        const QUALITY_STORAGE_KEY='galaxySimQuality';
        let qualityChoice=null,autoLite=false; // qualityChoice: 'high', 'low' or null (automatic)
        try{const saved=localStorage.getItem(QUALITY_STORAGE_KEY);if(saved==='high'||saved==='low')qualityChoice=saved;}catch(e){}
        function buildResolutionLevels(){
            const native=Math.min(window.devicePixelRatio||1,2);
            resLevels=native<1.5?[{pixelRatio:native,msaa:true}]:[];
            RES_SCALES.forEach(s=>{const pixelRatio=Math.round(native*s*100)/100;if(pixelRatio>=RES_MIN_PIXEL_RATIO||s===1)resLevels.push({pixelRatio,msaa:false});});
            resLevel=qualityChoice==='high'?0:qualityChoice==='low'?resLevels.length-1:Math.min(resLevel,resLevels.length-1);
        }
        buildResolutionLevels();
        const qualityToggle=document.getElementById('qualityToggle'),qualityNote=document.getElementById('qualityNote'),qualityAuto=document.getElementById('qualityAuto');
        function syncQualitySetting(){
            if(!qualityChoice&&resLevel===resLevels.length-1&&resLevels.length>1)autoLite=true;
            const lite=qualityChoice==='low'||!qualityChoice&&autoLite;
            if(document.documentElement.classList.contains('lite')!==lite){document.documentElement.classList.toggle('lite',lite);needsRender=true;}
            if(bloomPass)bloomPass.enabled=!lite;
            const full=resLevel===0&&!lite;
            qualityToggle.checked=qualityChoice?qualityChoice==='high':full;
            qualityNote.textContent=qualityChoice?'Your choice, remembered on this device.'
                :full?'Automatic: high quality on this device.':'Automatic: lowered to keep the 3D view smooth.';
            qualityAuto.hidden=!qualityChoice;
        }
        function setQualityChoice(choice){
            qualityChoice=choice;autoLite=false;
            try{if(choice)localStorage.setItem(QUALITY_STORAGE_KEY,choice);else localStorage.removeItem(QUALITY_STORAGE_KEY);}catch(e){}
            resUpgradeWindows=RES_UPGRADE_WINDOWS;resLastDowngradeFrom=-1; // automatic starts afresh
            setResolutionLevel(choice==='high'?0:choice==='low'?resLevels.length-1:resLevel);
        }
        qualityToggle.disabled=false;
        qualityToggle.addEventListener('change',()=>setQualityChoice(qualityToggle.checked?'high':'low'));
        qualityAuto.addEventListener('click',()=>{setQualityChoice(null);qualityToggle.focus();});
        function setResolutionLevel(level){
            resLevel=level;resFrames=resSlow=resFastWindows=0;
            resizeOrbit();syncQualitySetting();
        }
        syncQualitySetting();
        // Called once per drawn frame with the time since the previous one (only when the previous frame was drawn too).
        function trackFrameTime(ms){
            if(qualityChoice)return; // the visitor fixed the quality in Settings
            if(ms>RES_IGNORE_MS)return; // a background tab or a one-off stall, not the steady frame rate
            resFrames++;if(ms>RES_SLOW_MS)resSlow++;
            if(resFrames<RES_WINDOW_FRAMES)return;
            const slow=resSlow>RES_WINDOW_FRAMES/2,fast=resSlow===0&&ms<RES_FAST_MS;
            resFrames=resSlow=0;
            if(slow&&resLevel<resLevels.length-1){
                // Falling straight back from a level just tried means it is too slow: wait twice as long before the next try.
                resUpgradeWindows=resLastDowngradeFrom===resLevel?Math.min(resUpgradeWindows*2,64*RES_UPGRADE_WINDOWS):RES_UPGRADE_WINDOWS;
                resLastDowngradeFrom=resLevel;setResolutionLevel(resLevel+1);
            }else if(fast&&resLevel>0&&++resFastWindows>=resUpgradeWindows){
                resLastDowngradeFrom=resLevel-1;setResolutionLevel(resLevel-1);
            }else if(!fast)resFastWindows=0;
        }
        function resizeOrbit(){
            if (!nativeFsElement()) {
                containerOriginalWidth = orbitPlot.clientWidth;
                containerOriginalHeight = orbitPlot.clientHeight;
            }
            if(Math.min(window.devicePixelRatio||1,2)!==(resLevels[0]&&resLevels[0].pixelRatio))buildResolutionLevels(); // moved to another screen
            const w=Math.max(orbitPlot.clientWidth,1),h=Math.max(orbitPlot.clientHeight,1),{pixelRatio,msaa}=resLevels[resLevel];
            renderer.setPixelRatio(pixelRatio);renderer.setSize(w,h,false);
            if(composer){
                const samples=msaa?4:0;
                [composer.renderTarget1,composer.renderTarget2].forEach(t=>{if(t.samples!==undefined&&t.samples!==samples){t.samples=samples;t.dispose();}});
                composer.setPixelRatio(pixelRatio);composer.setSize(w,h);
            }
            starUniforms.uDpr.value=pixelRatio;
            camera.aspect=w/h;camera.updateProjectionMatrix();
            needsRender=true;
        }
        if(window.ResizeObserver)new ResizeObserver(resizeOrbit).observe(orbitPlot);
        window.addEventListener('resize',resizeOrbit);
        function onFullscreenChange() {
            const isFullscreen = Boolean(nativeFsElement() || orbitPlot.classList.contains('fullscreen-active'));
            if (fsBtn) fsBtn.setAttribute('aria-pressed', String(isFullscreen));
            // Layout settles after the change event, so resize now and on the next frame.
            resizeOrbit();
            requestAnimationFrame(resizeOrbit);
        }
        document.addEventListener('fullscreenchange', onFullscreenChange);
        document.addEventListener('webkitfullscreenchange', onFullscreenChange);
        let drewLastFrame=false;
        function renderLoop(now){
            const frameMs=now-lastTime,dt=Math.min(frameMs/1000,0.1);
            lastTime=now;
            if(controls.update(dt))needsRender=true;
            if(playing){
                simTime+=dt*PLAYBACK_SPEED_MYR_PER_SEC*playbackMultiplier;
                if(simTime>=SIM_TOTAL_TIME_MYR){
                    simTime=SIM_TOTAL_TIME_MYR;
                    setPlaying(false);
                }
                updateStarTransforms(simTime);
                updateEscapeCounter(simTime);
                const intMyr=Math.floor(simTime);
                if(intMyr!==lastReportedFrame){
                    lastReportedFrame=intMyr;
                    syncTimelineUI(simTime);
                }
                needsRender=true;
            }
            if(needsRender&&drewLastFrame)trackFrameTime(frameMs); // only back-to-back frames measure the frame rate
            drewLastFrame=needsRender;
            if(needsRender){needsRender=false;haloUniforms.uScale.value=starUniforms.uScale.value=mondHazeUniforms.uScale.value=renderer.domElement.height/(2*Math.tan(camera.fov*Math.PI/360));
                // Diffuse glow: dimmer close up so single stars stand out; the disk plane brightens edge-on (longer line of sight).
                const glowZoom=0.35+0.65*smoothStep(20,GAS_RESET_RADIUS,controls.view.radius);
                const tilt=Math.abs(Math.cos(controls.view.phi));
                diskGlowUniforms.uGain.value=0.6*glowZoom*Math.min(3,1/Math.max(tilt,0.33));
                bulgeGlowUniforms.uGain.value=0.6*glowZoom;bulgeGlow.quaternion.copy(camera.quaternion);
                // Fade the boundary grid as the camera closes in on the visible galaxy.
                const haloFade=0.3+0.7*Math.min(1,Math.max(0,(controls.view.radius-30)/120));
                haloShell.material.uniforms.uFade.value=haloFade;gasMaterial.opacity=gasZoomOpacity(controls.view.radius);haloGrid.material.opacity=0.13*haloFade;updateBlackHole(simTime);if(composer)composer.render();else renderer.render(scene,camera);updateStarMarks();}
            requestAnimationFrame(renderLoop);
        }
        resizeOrbit();requestAnimationFrame(renderLoop);

    const initialHalo=haloCurve(DEFAULTS.massT,DEFAULTS.scaleRadius);
    const totalCurve=vBary.map((v,i)=>Math.hypot(v,initialHalo[i]));
    // χ² needs measured speeds, so the badge is hidden for the imaginary sandbox galaxy.
    // Reduced chi-squared divides by N - k, where k counts the values tuned to the galaxy: the halo's mass and scale
    // radius (NFW) or a0 (MOND best fit), plus the galaxy settings a best fit tunes (disk M/L, bulge M/L where there is
    // a bulge, distance and tilt). Visible matter only and plain MOND use SPARC's standard values and tune nothing.
    const MODEL_VALUES=DATA.model_tuned_values; // MODEL_TUNED_VALUES in the generator
    const tunedCount=(model,fit)=>MODEL_VALUES[model]+(fit?fit.tuned_settings:0);
    const reducedChi=(chi,pointCount,k)=>chi/Math.max(1,pointCount-k);
    // The fit in use: the halo fit (Fit halo to data, until Visible matter only) in NFW mode, or the MOND best fit
    // when Best realistic fits is on. Which halo fit depends on the switch (nfwFitOf). With the halo mass at zero
    // the model is visible matter only, so SPARC's standard values apply however the slider got there.
    let nfwFitApplied=true,appliedFit=null; // appliedFit: the fit applyBaryonState last applied
    function activeFit(){
      if(customActive)return null;
      if(mondMode)return bestFits?activeGalaxy.best_fit.mond:null;
      return nfwFitApplied&&Number(document.getElementById('massSlider').value)>0?nfwFitOf(activeGalaxy):null;
    }
    const activeFitModel=()=>mondMode?(bestFits?'mondBest':'mond'):(Number(document.getElementById('massSlider').value)>0?'nfw':'visible');
    // Plain-language verdict for a reduced chi-squared. The cut-offs are this site's rule of thumb, not a published standard.
    const CLOSE_FIT_RED_CHI=2,ROUGH_FIT_RED_CHI=10;
    const fitVerdict=red=>red<CLOSE_FIT_RED_CHI?{cls:'good',word:'Close fit'}:red<ROUGH_FIT_RED_CHI?{cls:'ok',word:'Rough fit'}:{cls:'poor',word:'Poor fit'};
    let shownRedChi=Infinity;  // the reduced χ² on the fit badge now (the tutorial's fitting goal reads it)
    let shownModelKms=null;    // the model speeds behind it, one per measured point
    // Names the model and fit the badge is scoring, so like-for-like pairs are easy to spot.
    function fitLabel(model){
      if(model==='visible')return 'Visible matter only';
      if(model==='mond')return 'Plain MOND';
      if(model==='mondBest')return 'MOND best fit';
      const fit=activeFit(),massT=Number(document.getElementById('massSlider').value),rs=Number(document.getElementById('scaleSlider').value);
      const isFit=fit&&Math.abs(massT-fit.m200_msun/1e12)<1e-6&&Math.abs(rs-fit.scale_radius_kpc)<1e-6;
      return isFit?(bestFits?'NFW best fit':'NFW halo fit'):fit&&bestFits?'Your halo + fit’s M/L, distance, tilt':'Your halo';
    }
    function updateChi(total){
      const b=document.getElementById('chiBadge');
      b.hidden=customActive;document.getElementById('chiFit').hidden=customActive;document.getElementById('tmChiRow').hidden=customActive;
      if(customActive)return;
      const model=activeFitModel(),k=tunedCount(model,activeFit());
      const chi=obsV.reduce((s,v,i)=>s+((v-total[i])/obsErr[i])**2,0),red=reducedChi(chi,obsV.length,k),verdict=fitVerdict(red);
      shownRedChi=red;shownModelKms=total;
      const miss=Math.sqrt(chi/obsV.length); // root-mean-square gap in error bars, as on the Scorecard
      document.getElementById('chiFit').textContent=`${fitLabel(model)} · ${k} tuned · about ${miss<10?miss.toFixed(1):miss.toFixed(0)} error bars off`;
      document.getElementById('chiVerdict').textContent=verdict.word;
      const nums=document.getElementById('chiNumbers');
      nums.innerHTML='χ²<sub>r</sub> '+red.toFixed(2);
      nums.title=`χ² = ${chi.toFixed(0)} over ${obsV.length} measured speeds, ${k} values tuned. Lower is better; about 1 is as close as the measurements allow.`;
      b.className=verdict.cls;
      const fc=document.getElementById('tm-chi');fc.textContent=red.toFixed(2);fc.className=verdict.cls;
    }
    const residualsOf=total=>obsV.map((v,i)=>v-total[i]);
    // Residual bars are one neutral grey (direction shows the sign), so they never share a colour with a curve or a fit verdict.
    const RESIDUAL_COLOR='rgba(170,184,204,0.6)',RESIDUAL_SWATCH='#aab8cc';
    const AXIS={gridcolor:'rgba(148,170,200,0.08)',gridwidth:1,zerolinecolor:'rgba(148,170,200,0.3)',linecolor:'rgba(148,170,200,0.32)',showline:true,mirror:false,tickcolor:'rgba(148,170,200,0.32)',ticks:'outside',ticklen:4,automargin:true,tickfont:{size:11,color:'#8e9db3',family:'JetBrains Mono,ui-monospace,Consolas,monospace'},title:{font:{size:12,color:'#b8c4d6'}}};
    // Residual bars are about half the mean spacing of the SPARC radii, so dense and sparse curves both read clearly.
    const residualBarWidth=()=>Math.round(53*(radii[radii.length-1]-radii[0])/(radii.length-1))/100;
    // The custom galaxy's stellar disk, bulge and gas speeds are drawn as thin dotted lines (hidden for presets).
    const COMPONENT_STYLE=[['Stellar disk','#ffd166'],['Bulge','#ff8fab'],['Gas','#7bd88f']];
    let componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
    const componentVisible=i=>customActive&&(i!==1||sandbox.bulgeMass>0); // the sliders' anatomy: the stars may still be catching up
    // "Before" line (a visual aid only): while a slider that reshapes the curve is moved, and for CURVE_GHOST_HOLD_MS
    // after it is let go, a faint wide band in the same colour stays where the changed line was, labelled "before".
    // Anatomy sliders leave the cyan visible-matter line behind; the halo sliders leave the orange total.
    const CURVE_GHOST_HOLD_MS=4000,CURVE_GHOST_FADE_IN_MS=300,CURVE_GHOST_FADE_OUT_MS=900; // it fades in and out, never pops
    const curveGhost={x:[],y:[],rgb:'95,208,224',visible:false,fade:0,leaving:false};
    const curveGhostTrace=()=>({type:'scatter',mode:'lines+text',name:'Before',showlegend:false,hoverinfo:'skip',visible:curveGhost.visible,opacity:curveGhost.fade,x:curveGhost.x,y:curveGhost.y,
      text:curveGhost.y.map((_,i,a)=>i===a.length-1?'before':''),textposition:'top left',textfont:{size:11,color:`rgba(${curveGhost.rgb},0.8)`},
      line:{color:`rgba(${curveGhost.rgb},0.2)`,width:7,shape:'spline',smoothing:0.6}});
    const buildCurveData=(halo,total)=>[
      {type:'scatter',mode:'markers',name:LEGEND_NAMES.observed,visible:!customActive,x:radii,y:obsV,error_y:{type:'data',array:obsErr,visible:true,color:'rgba(232,238,246,0.5)',thickness:1.1,width:3},marker:{size:6.5,color:'#f4f7fb',line:{color:'#070b14',width:1.4}},customdata:obsErr,hovertemplate:'Observed <b>%{y:.0f}</b> ± %{customdata:.0f} km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:LEGEND_NAMES.visible,x:radii,y:vBary.slice(),line:{color:'#5fd0e0',width:2.2,dash:'6px,4px',shape:'spline',smoothing:0.6},hovertemplate:'Visible <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Dark matter halo',x:radii,y:halo,fill:'tozeroy',fillcolor:'rgba(179,154,232,0.08)',line:{color:'#b39ae8',width:1.8,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:'Halo <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Total model',x:radii,y:total,line:{color:'#ffae5c',width:3,shape:'spline',smoothing:0.6},hovertemplate:'Model <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'bar',name:'Residual',visible:!customActive,x:radii,y:residualsOf(total),yaxis:'y2',showlegend:false,width:residualBarWidth(),marker:{color:RESIDUAL_COLOR},hovertemplate:'Residual <b>%{y:+.0f}</b><extra></extra>'},
      ...COMPONENT_STYLE.map(([name,color],i)=>({type:'scatter',mode:'lines',name,x:radii,y:componentCurves[i].slice(),visible:componentVisible(i),line:{color,width:1.6,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:name+' <b>%{y:.0f}</b><extra></extra>'})),
      curveGhostTrace() // last trace: updateHalo's Plotly.update leaves it alone
    ];
    // Narrow graphs (phones, and laptops under 1200px where it shares the row with the 3D view) use shorter legend names
    // so no entry is cut off.
    const NARROW_LEGEND=innerWidth<600||(innerWidth>960&&innerWidth<1200);
    const LEGEND_NAMES=NARROW_LEGEND?{observed:'Observed',visible:'Visible matter',noHalo:'No halo'}:{observed:'Observed (SPARC)',visible:'Visible matter (baryons)',noHalo:'Dark matter halo (none)'};
    const curveLayout={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',margin:{l:66,r:20,t:44,b:54},font:{color:textColor,size:12,family:'Inter,Segoe UI,system-ui,sans-serif'},
      xaxis:{...AXIS,fixedrange:true,title:{text:'Radius, R (kpc)',standoff:10,font:{size:12,color:'#b8c4d6'}},anchor:'y2',showgrid:true,zeroline:false,hoverformat:'.1f'},
      yaxis:{...AXIS,title:{text:'Orbital speed, V (km/s)',standoff:10,font:{size:12,color:'#b8c4d6'}},range:[0,190],fixedrange:true,domain:[0.32,1],zeroline:false},
      yaxis2:{...AXIS,title:{text:'Obs − model',standoff:6,font:{size:11,color:'#b8c4d6'}},domain:[0,0.24],anchor:'x',zeroline:true,zerolinewidth:1.5,nticks:4,fixedrange:true},
      bargap:0.2,legend:{orientation:'h',yanchor:'bottom',y:1.02,xanchor:'left',x:0,font:{size:11.5,color:'#cdd6e3'},itemsizing:'constant',itemwidth:36,tracegroupgap:4},
      hovermode:false,dragmode:false};
    // On a phone-width screen the legend sits two entries per row (with shorter names, see LEGEND_NAMES) and the
    // speed axis has a shorter title, leaving more height for the graph.
    if(NARROW_LEGEND){
      curveLayout.legend={...curveLayout.legend,entrywidth:0.5,entrywidthmode:'fraction',itemwidth:30,font:{...curveLayout.legend.font,size:11}};
      curveLayout.yaxis.title={...curveLayout.yaxis.title,text:'Speed, V (km/s)'};
    }
    const curveConfig={responsive:true,displaylogo:false,displayModeBar:false,scrollZoom:false,doubleClick:false};
    // The sandbox has no residuals, so its velocity axis fills the whole plot.
    // A preset galaxy's axis starts at its first measured radius (nothing is drawn inside it, so nothing is extrapolated)
    // and hugs the y axis: CURVE_X_PAD_FRACTION of the span either side keeps the end points and bars from being cut.
    const CURVE_X_PAD_FRACTION=0.012;
    const presetXRange=()=>{const a=radii[0],b=radii[radii.length-1],pad=(b-a)*CURVE_X_PAD_FRACTION;return[Math.max(0,a-pad),b+pad];};
    const curveLayoutNow=()=>({...(customActive?{...curveLayout,xaxis:{...curveLayout.xaxis,anchor:'y',range:[0,radii[radii.length-1]]},yaxis:{...curveLayout.yaxis,domain:[0,1]},yaxis2:{...curveLayout.yaxis2,visible:false}}:{...curveLayout,xaxis:{...curveLayout.xaxis,range:presetXRange()}}),...mondZoneLayout()});
    // MOND boost shading (MOND mode, every galaxy). MOND has no switch-on point: it boosts gravity at every radius, by
    // an amount that grows as the visible matter's own pull g_N = V_bar²/R weakens (Milgrom 1983; the size of the
    // boost follows the interpolating function in mondVelocity, McGaugh, Lelli & Schombert 2016). So the gold's
    // opacity is simply proportional to MOND's speed boost, V_MOND/V_bar − 1: MOND_ZONE_ALPHA_PER_BOOST × boost, so a
    // 1 % boost is all but invisible and 100 % is clear, capped at MOND_ZONE_ALPHA_MAX (a 200 % boost). There is no
    // line: the visitor judges where the boost becomes significant. g_N is taken in the same distance and tilt frame
    // as the MOND curve (mondModelVelocity), so with Best realistic fits on it follows the fit's distance, tilt and
    // a₀. Only the speed panel is shaded (y domain), not the residuals. The shading is a visual aid. The graph is cut
    // into MOND_ZONE_STRIPS strips, each shaded at its midpoint's boost, rounded to MOND_ZONE_LEVELS steps (fine
    // enough that a 1 % boost still gets its own faint step) so neighbouring strips of the same shade merge.
    const MOND_ZONE_ALPHA_PER_BOOST=0.2,MOND_ZONE_ALPHA_MAX=0.4,MOND_ZONE_STRIPS=120,MOND_ZONE_LEVELS=400;
    const mondBoostOfGN=gN=>{const vb=Math.sqrt(gN*M_PER_KPC)/1000;return mondVelocity(vb,1)/vb-1;};
    // ln g_N at any radius, linear between the radii where it is finite; held flat beyond them.
    function mondLogGNAt(logG,r){
      const pts=radii.map((x,i)=>[x,logG[i]]).filter(p=>isFinite(p[1]));
      if(!pts.length)return -Infinity;
      if(r<=pts[0][0])return pts[0][1];
      for(let i=1;i<pts.length;i++)if(r<=pts[i][0]){const[a,ga]=pts[i-1],[b,gb]=pts[i];return ga+(gb-ga)*(r-a)/(b-a);}
      return pts[pts.length-1][1];
    }
    // The shade at radius r, 0 (no boost) to 1 (MOND_ZONE_ALPHA_MAX), proportional to the boost; the graph and the
    // 3D plane both use it.
    function mondZoneLevel(logG,r){
      const lg=mondLogGNAt(logG,r),boost=isFinite(lg)?mondBoostOfGN(Math.exp(lg)):Infinity;
      return Math.min(1,Math.max(0,boost*MOND_ZONE_ALPHA_PER_BOOST/MOND_ZONE_ALPHA_MAX));
    }
    function mondZoneStrips([x0,x1],logG){
      const width=(radii[radii.length-1]-radii[0])/MOND_ZONE_STRIPS,n=Math.max(1,Math.ceil((x1-x0)/width)),strips=[];
      for(let k=0;k<n;k++){
        const a=x0+(x1-x0)*k/n,b=x0+(x1-x0)*(k+1)/n,level=Math.round(mondZoneLevel(logG,(a+b)/2)*MOND_ZONE_LEVELS);
        const last=strips[strips.length-1];
        if(last&&last.level===level)last.x1=b;else strips.push({x0:a,x1:b,level});
      }
      return strips.filter(s=>s.level>0).map(({x0,x1,level})=>({type:'rect',layer:'below',xref:'x',yref:'y domain',x0,x1,y0:0,y1:1,line:{width:0},
        fillcolor:`rgba(242,200,121,${(MOND_ZONE_ALPHA_MAX*level/MOND_ZONE_LEVELS).toFixed(4)})`}));
    }
    const mondLogGN=()=>radii.map((r,i)=>r>0&&vBary[i]>0?Math.log((vBary[i]*fitFrame.s*1000)**2/(r*fitFrame.f*M_PER_KPC)):-Infinity);
    // True when the visible pull is above a₀ somewhere, where MOND's boost is under ≈26 % and the gold is faint.
    const mondFaintSomewhere=()=>mondLogGN().some(g=>g>=Math.log(MOND_A0));
    function mondZoneLayout(){
      if(!mondMode){setMondZone3D(null);return{shapes:[],annotations:[]};}
      const logG=mondLogGN(),last=radii[radii.length-1];
      setMondZone3D(r=>mondZoneLevel(logG,r),last); // the 3D view's gold plane uses the same shade at every radius
      return{
        shapes:mondZoneStrips([radii[0],last],logG),
        // The label sits in the bottom-right corner, below the curves, where the graph is emptiest.
        annotations:[{xref:'x',yref:'y domain',x:last,xanchor:'right',xshift:-6,y:0.02,yanchor:'bottom',showarrow:false,
          text:'Deeper gold = bigger MOND boost',font:{size:11,color:'#f2c879'}}]
      };
    }
    Plotly.newPlot('curvePlot',buildCurveData(initialHalo,totalCurve),curveLayoutNow(),curveConfig);
    updateChi(totalCurve);
    setupCurveReadout();
    /* ---------- Rotation-curve readout ----------
       Replaces Plotly's hover box, which jumped from point to point. The numbers are always the values at a real
       SPARC radius (the nearest one to the pointer); only the box, guide line and dots glide between radii, so the
       motion is smooth without inventing values between measurements. Mouse: hover. Touch: tap or slide along the
       graph; the readout stays put (and updates live as sliders move) until you tap outside the graph. */
    function setupCurveReadout(){
      const gd=document.getElementById('curvePlot'),wrap=gd.parentElement;
      const layer=document.createElement('div');
      layer.className='curve-hover';layer.setAttribute('aria-hidden','true');
      layer.innerHTML='<div class="ch-line"></div><div class="ch-dots"></div><div class="ch-box"><div class="ch-r"></div><div class="ch-rows"></div></div>';
      wrap.appendChild(layer);
      const line=layer.querySelector('.ch-line'),dotsEl=layer.querySelector('.ch-dots'),box=layer.querySelector('.ch-box');
      const rEl=layer.querySelector('.ch-r'),rowsEl=layer.querySelector('.ch-rows');
      // [trace index, short label (null = 'Model' or 'MOND'), colour]; trace 4 (residual) is added as its own row.
      const ROWS=[[0,'Observed','#f4f7fb'],[3,null,'#ffae5c'],[1,'Visible','#5fd0e0'],[2,'Halo','#b39ae8'],[5,'Disk','#ffd166'],[6,'Bulge','#ff8fab'],[7,'Gas','#7bd88f']];
      const isOn=t=>!!t&&t.visible!==false&&t.visible!=='legendonly';
      let active=false,pinned=false,touchDown=false,idx=-1,target=0,shown=0,boxX=0,boxY=0,boxGoal=0,lastT=0,raf=0,dots=[];
      const axes=()=>{const fl=gd._fullLayout;return fl&&fl.xaxis&&fl.yaxis&&fl.xaxis.p2l?fl:null;};
      const radiiNow=()=>(gd.data&&gd.data[3]&&gd.data[3].x)||[];
      // Nearest SPARC radius to a screen point, or null if the point is outside the plotting area (e.g. on the legend).
      function indexAt(clientX,clientY){
        const fl=axes(),xs=radiiNow();if(!fl||!xs.length)return null;
        const rect=gd.getBoundingClientRect(),px=clientX-rect.left-fl.xaxis._offset,py=clientY-rect.top;
        if(px<-24||px>fl.xaxis._length+24||py<fl.yaxis._offset-8)return null;
        const x=fl.xaxis.p2l(px);let best=0;
        for(let i=1;i<xs.length;i++)if(Math.abs(xs[i]-x)<Math.abs(xs[best]-x))best=i;
        return best;
      }
      function fill(i){
        const d=gd.data,xs=radiiNow();if(!xs.length)return;
        rEl.textContent=`R = ${Number(xs[i]).toFixed(1)} kpc`;
        let html='';dots=[];
        ROWS.forEach(([k,label,colour])=>{
          const t=d[k];if(!isOn(t)||!t.y||t.y[i]==null)return;
          const v=Number(t.y[i]);
          const value=k===0?`${v.toFixed(0)} ± ${Number(t.customdata[i]).toFixed(0)}`:v.toFixed(0);
          html+=`<div class="ch-row"><i style="background:${colour}"></i><span>${label||(mondMode?'MOND':'Model')}</span><b>${value} km/s</b></div>`;
          dots.push([k,colour]);
        });
        if(isOn(d[4])&&d[4].y){const r=Number(d[4].y[i]);html+=`<div class="ch-row ch-res"><i style="background:${RESIDUAL_SWATCH}"></i><span>Residual</span><b>${r>=0?'+':'−'}${Math.abs(r).toFixed(0)} km/s</b></div>`;}
        rowsEl.innerHTML=html;
        dotsEl.innerHTML=dots.map(([,c])=>`<div class="ch-dot" style="color:${c}"></div>`).join('');
      }
      // Draw at a fractional radius index f (between two SPARC points while gliding).
      function place(f){
        const fl=axes(),xs=radiiNow(),d=gd.data;if(!fl||!xs.length)return;
        const i0=Math.max(0,Math.min(xs.length-1,Math.floor(f))),i1=Math.min(xs.length-1,i0+1),t=f-i0;
        const at=arr=>Number(arr[i0])+(Number(arr[i1])-Number(arr[i0]))*t;
        const xa=fl.xaxis,ya=fl.yaxis,ya2=fl.yaxis2;
        const x=xa._offset+xa.l2p(at(xs)),top=ya._offset;
        const bottom=isOn(d[4])&&ya2&&ya2._length?ya2._offset+ya2._length:ya._offset+ya._length;
        line.style.transform=`translate(${x}px,${top}px)`;line.style.height=`${bottom-top}px`;
        [...dotsEl.children].forEach((el,j)=>{
          const k=dots[j]&&dots[j][0],tr=d[k];if(!tr||!tr.y)return;
          const y=Math.max(top,Math.min(ya._offset+ya._length,ya._offset+ya.l2p(at(tr.y))));
          el.style.transform=`translate(${x}px,${y}px)`;
        });
        // The box sits in the lower part of the speed panel (usually empty: the curves rise to the top), right of
        // the line and flipping left near the edge; boxX eases so the flip is smooth too.
        const w=box.offsetWidth,W=gd.clientWidth;
        boxGoal=x+10+w>W-4?x-10-w:x+10;boxGoal=Math.max(4,Math.min(W-w-4,boxGoal));
        boxY=Math.max(top+6,ya._offset+ya._length-box.offsetHeight-8);
        box.style.transform=`translate(${boxX}px,${boxY}px)`;
        return x;
      }
      function tick(now){
        raf=0;
        const dt=lastT?Math.min(0.1,(now-lastT)/1000):1/60;lastT=now;
        const k=calmMotion()?1:1-Math.exp(-16*dt);
        shown+=(target-shown)*k;if(Math.abs(target-shown)<0.002)shown=target;
        place(shown);
        boxX+=(boxGoal-boxX)*k;if(Math.abs(boxGoal-boxX)<0.3)boxX=boxGoal;
        box.style.transform=`translate(${boxX}px,${boxY}px)`;
        if(active&&(shown!==target||boxX!==boxGoal))raf=requestAnimationFrame(tick);else lastT=0;
      }
      const kick=()=>{if(!raf)raf=requestAnimationFrame(tick);};
      function show(i){
        if(i!==idx){idx=i;fill(i);}
        target=i;
        if(!active){active=true;shown=i;lastT=0;place(i);boxX=boxGoal;layer.classList.add('on');}
        kick();
      }
      function hide(){active=pinned=false;idx=-1;layer.classList.remove('on');}
      wrap.addEventListener('pointermove',e=>{
        if(e.pointerType==='touch'&&!touchDown)return;
        const i=indexAt(e.clientX,e.clientY);
        if(i===null){if(!pinned)hide();return;}
        show(i);
      },{passive:true});
      wrap.addEventListener('pointerleave',e=>{if(e.pointerType!=='touch'&&!pinned)hide();},{passive:true});
      wrap.addEventListener('pointerdown',e=>{
        if(e.pointerType!=='touch')return;
        const i=indexAt(e.clientX,e.clientY);if(i===null)return;
        touchDown=true;pinned=true;show(i);
      },{passive:true});
      ['pointerup','pointercancel'].forEach(type=>window.addEventListener(type,()=>{touchDown=false;},{passive:true}));
      document.addEventListener('pointerdown',e=>{if(pinned&&!wrap.contains(e.target))hide();},{passive:true});
      // Slider moves, MOND, legend clicks and galaxy changes redraw the plot: refresh the numbers in place.
      gd.on('plotly_afterplot',()=>{
        if(!active)return;
        const n=radiiNow().length;if(!n){hide();return;}
        if(idx>n-1){target=shown=n-1;}
        const keep=Math.min(idx,n-1);idx=-1;fill(keep);idx=keep;place(shown);kick();
      });
      window.addEventListener('resize',()=>{if(active){place(shown);boxX=boxGoal;kick();}});
    }

    const massSlider=document.getElementById('massSlider'),scaleSlider=document.getElementById('scaleSlider');
    let sliderFrame=0,curveSettleTimer=null,sliderDragging=false;
    const fsSliders=document.getElementById('fsSliders'),fsMass=document.getElementById('fsMassSlider'),fsScale=document.getElementById('fsScaleSlider');
    // The mass slider grows for massive halos so a fitted value is never pinned to the end of its track. In the Custom
    // Sandbox it spans 0 to HALO_FINE_RANGE x the stellar-to-halo baseline, so it fine-tunes that baseline: haloFine is
    // M200 as a multiple of it, reset to 1 whenever the anatomy changes. rsManual records that the user took over r_s.
    const HALO_FINE_RANGE=5;
    let haloFine=1,rsManual=false;
    function configureHaloSliders(){
      if(customActive){
        const maxT=DEFAULTS.massT*HALO_FINE_RANGE;
        massSlider.max=fsMass.max=String(maxT);
        massSlider.step=fsMass.step=String(maxT/500);
        massSlider.value=haloFine*DEFAULTS.massT;
        if(!rsManual)scaleSlider.value=DEFAULTS.scaleRadius;
        return;
      }
      massSlider.step=fsMass.step='0.01';
      const maxT=String(Math.max(2,Math.ceil(2*DEFAULTS.massT/0.5)*0.5));
      massSlider.max=maxT;fsMass.max=maxT;
      massSlider.value=DEFAULTS.massT;scaleSlider.value=DEFAULTS.scaleRadius;
    }
    configureHaloSliders();
    const haloMassText=massT=>customActive?formatMass(massT*1e12):`${massT.toFixed(2)} × 10¹² M☉`;
    function refreshOrbits(){
      clearTimeout(orbitModelTimer);
      if(pendingCustomStars){if(customActive)rebuildCustomStars();else pendingCustomStars=null;}
      activeTrajectories=integrateTrajectories(Number(massSlider.value),Number(scaleSlider.value));
      fillOrbitLines();
      updateStarTransforms(simTime);
      updateEscapeCounter(simTime);
      requestRender();
    }
    function scheduleOrbitUpdate(){
      clearTimeout(orbitModelTimer);
      if(sliderDragging)return;
      orbitModelTimer=setTimeout(refreshOrbits,120);
    }
    function syncFsSliders(){
      fsMass.value=massSlider.value;fsScale.value=scaleSlider.value;
      document.getElementById('fsMassValue').textContent=haloMassText(Number(massSlider.value));
      document.getElementById('fsScaleValue').textContent=`${Number(scaleSlider.value).toFixed(1)} kpc`;
      document.getElementById('fsHaloRows').hidden=mondMode; // no halo in MOND, so no halo sliders
      document.getElementById('fsSandboxRows').hidden=!customActive;
    }
    [[fsMass,massSlider],[fsScale,scaleSlider]].forEach(([mini,main])=>{
      mini.addEventListener('input',()=>{main.value=mini.value;main.dispatchEvent(new Event('input'));});
      mini.addEventListener('change',()=>{main.value=mini.value;main.dispatchEvent(new Event('change'));});
      mini.addEventListener('pointerdown',()=>{sliderDragging=true;clearTimeout(orbitModelTimer);});
    });
    ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu','touchstart'].forEach(type=>fsSliders.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
    // Halo and total-model curves for the active galaxy under the current gravity model.
    function modelCurves(massT,rs){
      const halo=mondMode?vBary.map(()=>0):haloCurve(massT,rs);
      const total=mondMode?vBary.map((v,i)=>mondModelVelocity(v,radii[i])):vBary.map((v,i)=>Math.hypot(v,halo[i]));
      return{halo,total};
    }
    /* ---------- Custom sandbox ----------
       "Custom Sandbox" swaps the galaxy dropdown for six anatomy sliders (disk, bulge, gas) and a halo stability
       monitor. sandbox holds the current anatomy; a 'sandboxchange' event announces every change. The monitor
       never blocks anything: out-of-range halos only raise a warning so unphysical "monster galaxies" can be run. */
    const [STABLE_C_MIN,STABLE_C_MAX]=DATA.stable_concentration;
    const superscript=n=>String(n).replace('-','⁻').replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);
    const formatMass=m=>{
      if(!(m>0))return '0 M☉';
      const e=Math.floor(Math.log10(m));
      return `${(m/10**e).toFixed(2)} × 10${superscript(e)} M☉`;
    };
    const formatKpc=v=>`${v.toFixed(1)} kpc`;
    const formatPercent=f=>`${(100*f).toFixed(1)}%`;
    const SANDBOX_SLIDERS=[
      {key:'diskMass',id:'sbDiskMass',format:formatMass},
      {key:'diskScale',id:'sbDiskScale',format:formatKpc},
      {key:'truncation',id:'sbTruncation',format:formatKpc},
      {key:'bulgeMass',id:'sbBulgeMass',format:formatMass},
      {key:'bulgeRadius',id:'sbBulgeRadius',format:formatKpc},
      {key:'gasMass',id:'sbGasMass',format:formatMass},
    ].map(slider=>({...slider,input:document.getElementById(slider.id),output:document.getElementById(slider.id+'Value')}));
    const sandbox={};
    const sandboxPanel=document.getElementById('sandboxPanel');
    const bulgeRadiusControl=document.getElementById('sbBulgeRadiusControl');
    const bulgeRadiusHelp=document.getElementById('sbBulgeRadiusHelp');
    const BULGE_RADIUS_HELP=bulgeRadiusHelp.innerHTML;
    // Galaxy type card: names the kind of galaxy the anatomy sliders describe. Thresholds are rules of thumb:
    // Hubble-type bulge-to-total ratios (Graham & Worley 2008) and low surface brightness meaning a central stellar
    // surface density well below Freeman's (1970) typical disk value.
    function classifyGalaxy(g){
      const stellar=g.diskMass+g.bulgeMass,bt=g.bulgeMass/stellar,gasFrac=g.gasMass/(stellar+g.gasMass);
      const sigma0=g.diskMass/(2*Math.PI*(g.diskScale*1e3)**2); // central disk surface density, M☉/pc²
      const tags=[];
      if(stellar>=8e10)tags.push('Giant');
      if(stellar<3e9)tags.push('Dark matter dominated');
      if(gasFrac>=0.4)tags.push('Gas-rich: lots of fuel for new stars');else if(gasFrac<0.1)tags.push('Gas-poor: few new stars');
      if(g.bulgeMass<=0)tags.push('No bulge');else if(bt>=0.4)tags.push('Big bulge');
      if(sigma0<100)tags.push('Faint, spread-out disk');else if(g.diskScale<=2)tags.push('Compact disk');
      let type;
      if(stellar>=1.2e11&&g.gasMass>=1.5e10)type={name:'Monster galaxy',
        text:'Every dial turned up. This is heavier than almost any spiral known: a bright, gas-stuffed giant that needs a huge halo to keep its outer stars from flying off.',
        example:'Closest real cousin: UGC 2885 (“Rubin’s galaxy”), one of the largest spirals known.'};
      else if(stellar<3e9&&gasFrac>=0.5)type={name:'Gas-rich dwarf irregular',
        text:'A small, scruffy galaxy with more gas than stars. Its gravity is so weak that dark matter (or MOND) does almost all the work, which makes dwarfs some of the best tests in this project.',
        example:'Real example: DDO 154, where the gas outweighs the stars.'};
      else if(stellar<3e9)type={name:'Dwarf galaxy',
        text:'A lightweight with only a few billion suns of stars. Dwarfs are the most common galaxies in the Universe, but they are faint, and their slow-spinning stars are mostly held in by dark matter.',
        example:'Real example: the Large Magellanic Cloud, a dwarf orbiting the Milky Way.'};
      else if(sigma0<100&&bt<0.25)type={name:'Low surface brightness galaxy',
        text:'A ghostly disk: its stars are spread so thinly that it barely shows up against the night sky. These galaxies are packed with dark matter compared with their stars.',
        example:'Real examples: UGC 128, and the giant Malin 1.'};
      else if(bt>=0.6&&gasFrac<0.1)type={name:'Lenticular (S0) galaxy',
        text:'A big, bright bulge with a smooth disk and almost no gas left. With no fuel, it has stopped making new stars: halfway between a spiral and an elliptical.',
        example:'Real example: NGC 5866, the Spindle Galaxy.'};
      else if(bt>=0.6)type={name:'Bulge-heavy spiral (Sa)',
        text:'The central bulge outshines the disk, so the spiral arms are tightly wound and faint. Push the bulge any further and it starts to look like an elliptical.',
        example:'Real example: M104, the Sombrero Galaxy.'};
      else if(bt>=0.25)type={name:'Spiral with a bright bulge (Sb)',
        text:'A classic grand spiral: a glowing central bulge wrapped in a disk of stars and gas. The bulge makes the inner stars orbit fast, so the curve rises steeply.',
        example:'Real example: Andromeda (M31), our nearest big neighbour.'};
      else if(bt>0)type={name:'Late-type spiral (Sbc–Sc)',
        text:'Mostly disk with a small bulge at its heart and open, loose spiral arms. Plenty of gas keeps new stars forming.',
        example:'Real example: the Milky Way, our home galaxy.'};
      else type={name:'Bulgeless disk galaxy (Sc–Sd)',
        text:'A flat, pure disk with no central bulge. Its rotation curve stays flat far beyond the stars, one of the clearest signs that something unseen is there.',
        example:'Real examples: NGC 3198 (this project’s main galaxy) and M33.'};
      if(stellar>=8e10&&type.name!=='Monster galaxy'&&!type.name.startsWith('Giant'))type.name='Giant '+type.name.charAt(0).toLowerCase()+type.name.slice(1);
      return{...type,tags};
    }
    function updateGalaxyType(g){
      const t=classifyGalaxy(g);
      const nameEl=document.getElementById('gtName');
      if(nameEl.textContent!==t.name)nameEl.textContent=t.name;
      document.getElementById('gtText').textContent=t.text;
      document.getElementById('gtExample').textContent=t.example;
      document.getElementById('gtTags').replaceChildren(...t.tags.map(tag=>{const li=document.createElement('li');li.textContent=tag;return li;}));
    }
    function updateSandbox(){
      SANDBOX_SLIDERS.forEach(({key,input,output,format})=>{
        sandbox[key]=Number(input.value);
        output.textContent=format(sandbox[key]);
      });
      const noBulge=sandbox.bulgeMass<=0;
      bulgeRadiusControl.classList.toggle('disabled',noBulge);
      bulgeRadiusControl.setAttribute('aria-disabled',String(noBulge));
      document.getElementById('sbBulgeRadius').disabled=noBulge;
      bulgeRadiusHelp.innerHTML=noBulge?'No bulge: this is a pure disk galaxy.':BULGE_RADIUS_HELP;
      const stellar=sandbox.diskMass+sandbox.bulgeMass,visible=stellar+sandbox.gasMass;
      document.getElementById('sbTotalMass').textContent=formatMass(visible);
      document.getElementById('sbBulgeFraction').textContent=formatPercent(sandbox.bulgeMass/stellar);
      document.getElementById('sbGasFraction').textContent=formatPercent(sandbox.gasMass/visible);
      updateGalaxyType(sandbox);
      fsSandboxSliders.forEach(({mini,input,output,miniOutput})=>{mini.value=input.value;miniOutput.textContent=output.textContent;});
      document.dispatchEvent(new CustomEvent('sandboxchange',{detail:{...sandbox}}));
    }
    // The quick-controls box (fullscreen, or docked in the 3D view) also carries the sandbox's three mass sliders.
    const FS_SANDBOX_KEYS=['diskMass','bulgeMass','gasMass'];
    const fsSandboxSliders=SANDBOX_SLIDERS.filter(s=>FS_SANDBOX_KEYS.includes(s.key)).map(s=>{
      const mini=document.createElement('input'),label=document.createElement('label'),name=document.createElement('span'),miniOutput=document.createElement('output');
      mini.type='range';mini.id='fs-'+s.id;mini.min=s.input.min;mini.max=s.input.max;mini.step=s.input.step;mini.value=s.input.value;
      label.htmlFor=mini.id;name.textContent=s.input.labels[0].querySelector('span').textContent.split(',')[0];
      label.append(name,miniOutput);document.getElementById('fsSandboxRows').append(label,mini);
      mini.addEventListener('input',()=>{s.input.value=mini.value;s.input.dispatchEvent(new Event('input'));});
      return {...s,mini,miniOutput};
    });
    SANDBOX_SLIDERS.forEach(({input})=>input.addEventListener('input',updateSandbox));
    document.getElementById('sandboxReset').addEventListener('click',()=>{
      SANDBOX_SLIDERS.forEach(({input})=>{input.value=input.defaultValue;});
      updateSandbox();
    });
    updateSandbox();

    // Halo shape monitor: concentration c = r200 / rs of the current NFW halo.
    const stabilityMonitor=document.getElementById('stabilityMonitor'),stabilityBadge=document.getElementById('stabilityBadge');
    const stabilityChip=document.getElementById('stabilityChip');
    const C_GAUGE_MIN=1,C_GAUGE_MAX=100;
    const gaugePosition=c=>100*Math.min(1,Math.max(0,Math.log(c/C_GAUGE_MIN)/Math.log(C_GAUGE_MAX/C_GAUGE_MIN)))+'%';
    const cBand=document.getElementById('cBand');
    cBand.style.left=gaugePosition(STABLE_C_MIN);
    cBand.style.width=`${parseFloat(gaugePosition(STABLE_C_MAX))-parseFloat(gaugePosition(STABLE_C_MIN))}%`;
    [C_GAUGE_MIN,STABLE_C_MIN,10,STABLE_C_MAX,C_GAUGE_MAX].forEach(c=>{
      const tick=document.createElement('span');
      tick.textContent=c;tick.style.left=gaugePosition(c);
      document.getElementById('cTicks').append(tick);
    });
    const UNSTABLE_BADGE='⚠️ Unusual halo shape';
    function stabilityReport(){
      if(mondMode)return{state:'none',badge:'No halo to assess',message:'MOND has no dark matter halo, so there is no concentration to check. Choose Dark matter halo to see it.'};
      const massT=Number(massSlider.value),rs=Number(scaleSlider.value);
      if(massT<=0)return{state:'none',badge:'No halo to assess',message:'With the halo mass at zero there is no halo to check. Raise it to see its concentration.'};
      const r200=Math.cbrt(3*massT*1e12/(4*Math.PI*200*(3*H0*H0/(8*Math.PI*G))));
      const c=r200/rs,range=`${STABLE_C_MIN} and ${STABLE_C_MAX}`;
      const report={c,r200,rs};
      if(c<STABLE_C_MIN)return{...report,state:'diffuse',badge:UNSTABLE_BADGE,
        message:`c is below ${STABLE_C_MIN}: this halo is far more spread out than typical NFW halos, whose concentrations mostly lie between ${range} (a rough guide, not a strict limit). You can still run the simulation to see how the galaxy behaves.`};
      if(c>STABLE_C_MAX)return{...report,state:'dense',badge:UNSTABLE_BADGE,
        message:`c is above ${STABLE_C_MAX}: this halo packs its mass into an extremely dense centre, beyond typical NFW halos, whose concentrations mostly lie between ${range} (a rough guide, not a strict limit). You can still run the simulation to see how the galaxy behaves.`};
      return{...report,state:'stable',badge:'✓ Typical halo shape',message:`c is inside the ${STABLE_C_MIN}–${STABLE_C_MAX} range typical of NFW halos. The monitor follows the halo mass and scale radius sliders.`};
    }
    function updateStability(){
      const report=stabilityReport(),unstable=report.state==='diffuse'||report.state==='dense';
      stabilityMonitor.dataset.state=report.state;
      if(stabilityBadge.textContent!==report.badge)stabilityBadge.textContent=report.badge;
      document.getElementById('stabilityC').textContent=report.c===undefined?'c = —':`c = ${report.c.toFixed(1)}`;
      document.getElementById('stabilityDetail').innerHTML=report.c===undefined?'':`r<sub>200</sub> ${report.r200.toFixed(0)} kpc ÷ r<sub>s</sub> ${report.rs.toFixed(1)} kpc`;
      document.getElementById('cMarker').style.left=report.c===undefined?'0%':gaugePosition(report.c);
      document.getElementById('stabilityMessage').textContent=report.message;
      const showChip=galaxyMode==='custom'&&unstable;
      stabilityChip.hidden=!showChip;
      if(showChip&&stabilityChip.textContent!==UNSTABLE_BADGE)stabilityChip.textContent=UNSTABLE_BADGE;
    }

    let galaxyMode='preset';
    const galaxyModeRadios=[...document.querySelectorAll('input[name=galaxyMode]')];
    function setGalaxyMode(mode){
      const changed=mode!==galaxyMode;
      galaxyMode=mode;
      const custom=mode==='custom';
      galaxyModeRadios.forEach(radio=>{radio.checked=radio.value===mode;});
      sandboxPanel.hidden=!custom;
      stabilityMonitor.hidden=!custom;
      document.getElementById('presetPicker').hidden=custom; // nothing real to compare an imaginary galaxy with
      if(changed){if(custom)activateCustom();else loadPreset(activeGalaxy);}
      updateStability();
      document.dispatchEvent(new CustomEvent('galaxymodechange',{detail:{mode}}));
    }
    galaxyModeRadios.forEach(radio=>radio.addEventListener('change',()=>{if(radio.checked)setGalaxyMode(radio.value);}));
    setGalaxyMode(galaxyModeRadios.find(radio=>radio.checked).value);

    // A full Plotly.update re-lays out the whole chart (~40 ms), too slow to follow a slider. While a slider moves, the
    // curves are redrawn through Plotly's transition path instead (about 1 ms): only the line, fill and bar shapes move,
    // and the speed axis only grows if a curve would leave it. Once the slider settles, a full update tidies everything
    // (axis fit, names), CURVE_SETTLE_MS after the last move or at once when it is let go.
    const CURVE_SETTLE_MS=150;
    const LIVE_CURVE_ANIMATION={transition:{duration:0},frame:{duration:0,redraw:false},mode:'immediate'};
    function updateHalo(live){
      live=live===true; // event listeners pass an Event here
      if(activeFit()!==appliedFit)applyBaryonState(); // the halo slider reached or left zero
      syncFsSliders();
      updateStability();
      cancelAnimationFrame(sliderFrame);sliderFrame=0;clearTimeout(curveSettleTimer);curveSettleTimer=null;
      const massT=Number(massSlider.value),rs=Number(scaleSlider.value);
      const{halo,total}=modelCurves(massT,rs);
      document.getElementById('massValue').textContent=haloMassText(massT);
      document.getElementById('scaleValue').textContent=`${rs.toFixed(1)} kpc`;
      const res=residualsOf(total);updateChi(total);
      document.getElementById('mondZoneNote').hidden=!mondMode;
      const fit=activeFit();
      document.getElementById('chiBadge').title=(mondMode?'Goodness of fit of the MOND prediction to the SPARC points':'Goodness of fit of visible matter + NFW halo to the SPARC points')
        +(fit&&bestFits?`, using the best fit's distance (${fit.distance_mpc.toFixed(1)} Mpc) and tilt (${fit.inclination_deg.toFixed(0)}°)`:'');
      const yMax=Math.max(...obsV.map((v,i)=>v+obsErr[i]),...vBary,...halo,...total,...(curveGhost.visible?curveGhost.y:[]))*1.1;
      // One Plotly call redraws the model lines, the three custom-galaxy component lines and the residual bars (trace 4).
      const curvePlot=document.getElementById('curvePlot'),visible=[true,!mondMode,true,!customActive,...COMPONENT_STYLE.map((_,i)=>componentVisible(i))];
      // A line appearing or vanishing (e.g. the bulge set to zero) needs the full update even mid-drag.
      if(live&&visible.every((v,k)=>curvePlot.data[k+1].visible===v)){
        const range=curvePlot.layout.yaxis.range;
        Plotly.animate(curvePlot,{data:[vBary.slice(),halo,total,res,...componentCurves.map(c=>c.slice())].map(y=>({y})),
          traces:[1,2,3,4,5,6,7],layout:yMax>range[1]?{yaxis:{range:[0,yMax]}}:{}},LIVE_CURVE_ANIMATION).catch(()=>{}); // a newer move interrupts it
        curveSettleTimer=setTimeout(updateHalo,CURVE_SETTLE_MS);
        updateHaloVisual(massT,rs);
        clearTimeout(orbitModelTimer); // the settling full update schedules the orbits once
        return;
      }
      Plotly.update('curvePlot',{y:[vBary.slice(),halo,total,res,...componentCurves.map(c=>c.slice())],name:[LEGEND_NAMES.visible,mondMode?LEGEND_NAMES.noHalo:'Dark matter halo',mondMode?'MOND model':'Total model','Residual',...COMPONENT_STYLE.map(([name])=>name)],visible},{'yaxis.range':[0,yMax],...mondZoneLayout()},[1,2,3,4,5,6,7]);
      updateHaloVisual(massT,rs);
      scheduleOrbitUpdate();
    }
    const mondToggle=document.getElementById('mondToggle'),bestFitToggle=document.getElementById('bestFitToggle');
    const formatA0=a0=>{const e=Math.floor(Math.log10(a0));return `${(a0/10**e).toFixed(2)}×10${String(e).replace('-','⁻').replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`;};
    function applyBaryonState(){
      // The best fit in use sets the M/L ratios, a0, distance and tilt; otherwise SPARC's standard values apply.
      // A custom galaxy is built from masses, so none of these apply to it (activeFit is null).
      const fit=activeFit();
      appliedFit=fit;
      mlState.disk=fit?fit.upsilon_disk:DATA.upsilon_disk;
      mlState.bulge=fit?fit.upsilon_bulge:DATA.upsilon_bulge;
      MOND_A0=fit&&mondMode?fit.a0_ms2:MOND_A0_DEFAULT;
      fitFrame=frameOf(activeGalaxy,fit);
      baryonsFor(mlState.disk,mlState.bulge).forEach((v,i)=>{vBary[i]=v;});
      document.getElementById('bestFitWrap').hidden=customActive;
      const mlStat=document.getElementById('mlStat');
      mlStat.classList.toggle('off',customActive); // hidden in the sandbox, which sets stellar masses directly
      document.getElementById('tiltStat').classList.toggle('off',customActive); // a sandbox galaxy has no measured tilt
      document.getElementById('mlDisk').textContent=mlState.disk.toFixed(2);
      document.getElementById('mlBulge').textContent=mlState.bulge.toFixed(2);
      document.getElementById('ftDisk').textContent=mlState.disk.toFixed(2);
      document.getElementById('ftBulge').textContent=mlState.bulge.toFixed(2);
      const mlSource=fit&&bestFits?(mondMode?'MOND best fit':'NFW best fit'):'standard SPARC values';
      document.getElementById('mlSource').textContent=mlSource;
      document.getElementById('ftMlSource').textContent=mlSource;
      mlStat.title=`Stellar mass-to-light ratios in use (${mlSource}). Best realistic fits may change them within realistic limits; otherwise SPARC's standard 0.50 and 0.70 apply.`;
      document.getElementById('mondA0Text').textContent=formatA0(MOND_A0);
      if(!customActive){
        const sparcTilt=activeGalaxy.catalogue.inclination_deg,tilt=fit?fit.inclination_deg:sparcTilt;
        document.getElementById('tiltValue').textContent=`${tilt.toFixed(0)}°`;
        document.getElementById('tiltSource').textContent=fit&&bestFits?`${mlSource} (SPARC ${sparcTilt.toFixed(0)}°)`:'standard SPARC value';
      }
    }
    bestFitToggle.addEventListener('change',()=>{
      bestFits=bestFitToggle.checked;
      DEFAULTS=haloDefaultsOf(activeGalaxy);
      if(!customActive){ // the halo sliders follow the new fit unless the halo was set to zero (Visible matter only)
        const m=massSlider.value,rs=scaleSlider.value;
        configureHaloSliders();
        if(!nfwFitApplied){massSlider.value=m;scaleSlider.value=rs;}
      }
      applyBaryonState();updateHalo();
    });
    applyBaryonState();updateHalo(); // the page opens with the preset's halo fit at SPARC's standard values
    // "?" on the M/L and tilt tiles: a mouse hover fades the explainer in and moving away fades it out
    // (a short grace period lets the pointer travel into the box). Touch and keyboard users tap/press
    // the button instead, and close it with ×, Escape or a tap elsewhere.
    // Hover explainers (halo cards, the ? boxes): a mouse rests this long before one opens, and it counts as read
    // (for the tutorial) only once it has stayed open this long, so a brushed-past box is neither opened nor ticked.
    const HOVER_OPEN_DELAY_MS=350;
    const EXPLAINER_READ_MS=1500;
    // Fires an 'explainerread' event on el once it has stayed open for EXPLAINER_READ_MS; closing early cancels it.
    function trackExplainerRead(el,open){
      if(!open){clearTimeout(el._readTimer);el._readTimer=null;return;}
      if(!el._readTimer)el._readTimer=setTimeout(()=>el.dispatchEvent(new Event('explainerread')),EXPLAINER_READ_MS);
    }
    function setupExplainer(btnId,boxId,closeId){
      const btn=document.getElementById(btnId),box=document.getElementById(boxId);
      let hover=false,pinned=false,closeTimer=null,lastPointer='';
      const isOpen=()=>box.classList.contains('open');
      const render=()=>{
        const open=hover||pinned;
        if(open&&!isOpen()){
          // Right-aligned under the tile by default; flip to left-aligned if that would run off a narrow screen.
          box.style.left='';box.style.right='';
          if(box.getBoundingClientRect().left<16){box.style.left='0';box.style.right='auto';}
        }
        box.classList.toggle('open',open);
        box.setAttribute('aria-hidden',String(!open));
        btn.setAttribute('aria-expanded',String(open));
        trackExplainerRead(btn,open);
      };
      const close=()=>{hover=false;pinned=false;clearTimeout(closeTimer);render();};
      [btn,box].forEach(el=>{
        // Already open (moving from the ? into its box): stay open. Otherwise rest on it for a moment first.
        el.addEventListener('pointerenter',e=>{if(e.pointerType!=='mouse')return;clearTimeout(closeTimer);closeTimer=setTimeout(()=>{hover=true;render();},isOpen()?0:HOVER_OPEN_DELAY_MS);});
        el.addEventListener('pointerleave',e=>{if(e.pointerType!=='mouse')return;clearTimeout(closeTimer);closeTimer=setTimeout(()=>{hover=false;render();},150);});
      });
      btn.addEventListener('pointerdown',e=>{lastPointer=e.pointerType;});
      btn.addEventListener('click',e=>{
        e.stopPropagation();
        if(lastPointer==='mouse'){lastPointer='';return;}  // mouse users already see it on hover
        lastPointer='';pinned=!pinned;render();
      });
      document.getElementById(closeId).addEventListener('click',()=>{close();btn.focus();});
      document.addEventListener('click',e=>{if(isOpen()&&!box.contains(e.target)&&e.target!==btn)close();});
      document.addEventListener('keydown',e=>{if(e.key==='Escape'&&isOpen()){close();btn.focus();}});
    }
    setupExplainer('mlInfoBtn','mlExplainer','mlExplainerClose');
    setupExplainer('tiltInfoBtn','tiltExplainer','tiltExplainerClose');
    setupExplainer('starsInfoBtn','starsExplainer','starsExplainerClose');
    // The two model cards set the hidden mondToggle, which every other listener (tutorial, "before" line) watches.
    const gravityRadios=[...document.querySelectorAll('input[name=gravityModel]')];
    gravityRadios.forEach(radio=>radio.addEventListener('change',()=>{
      if(!radio.checked||(radio.value==='mond')===mondToggle.checked)return;
      mondToggle.checked=radio.value==='mond';mondToggle.dispatchEvent(new Event('change',{bubbles:true}));
    }));
    mondToggle.addEventListener('change',()=>{
      mondMode=mondToggle.checked;
      gravityRadios.forEach(radio=>{radio.checked=(radio.value==='mond')===mondMode;});
      document.getElementById('mondNote').hidden=!mondMode;
      document.getElementById('nfwNote').hidden=mondMode;
      applyBaryonState();
      // MOND has no halo, so its sliders are hidden rather than greyed out.
      ['massControl','scaleControl'].forEach(id=>{document.getElementById(id).hidden=mondMode;});
      updateHalo();
    });
    function queueSliderUpdate(){
      updateStability();
      document.getElementById('massValue').textContent=haloMassText(Number(massSlider.value));
      document.getElementById('scaleValue').textContent=`${Number(scaleSlider.value).toFixed(1)} kpc`;
      if(!sliderFrame)sliderFrame=requestAnimationFrame(()=>{sliderFrame=0;updateHalo(true);}); // at most one redraw per frame
    }
    [massSlider,scaleSlider].forEach(slider=>{
      slider.addEventListener('input',queueSliderUpdate);
      slider.addEventListener('change',updateHalo);
      slider.addEventListener('pointerdown',()=>{sliderDragging=true;clearTimeout(orbitModelTimer);});
    });
    // Only the user's own drags change haloFine / rsManual; the sandbox moves the sliders in code without input events.
    ['input','change'].forEach(type=>{
      massSlider.addEventListener(type,()=>{if(customActive)haloFine=Number(massSlider.value)/DEFAULTS.massT;});
      scaleSlider.addEventListener(type,()=>{if(customActive)rsManual=true;});
    });
    function finishSliderDrag(){
      if(!sliderDragging)return;
      sliderDragging=false;
      scheduleOrbitUpdate();
    }
    // Capture phase: the fullscreen panel stops pointer events bubbling (so drags don't spin the camera), and a
    // bubbling listener here never heard its sliders being released, so the orbits were never rebuilt.
    window.addEventListener('pointerup',finishSliderDrag,true);
    window.addEventListener('pointercancel',finishSliderDrag,true);
    updateHalo();
    document.getElementById('applyObservedFit').addEventListener('click',()=>{
        massSlider.value=DEFAULTS.massT;
        scaleSlider.value=DEFAULTS.scaleRadius;
        haloFine=1;rsManual=false;
        nfwFitApplied=true;applyBaryonState();
        updateHalo();
    });
    document.getElementById('visibleMatterOnly').addEventListener('click',()=>{
        massSlider.value=0;
        haloFine=0;
        nfwFitApplied=false;applyBaryonState(); // visible matter with SPARC's standard values
        updateHalo();
    });

    /* ---------- Galaxy presets and the Custom Sandbox physics ----------
       The dropdown swaps the active SPARC galaxy: its baryonic baseline, observed halo fit, slider ranges,
       rotation-curve plot, chi-squared and orbits are all rebuilt from that galaxy's galaxies.js entry.
       The Custom Sandbox hides the dropdown: its imaginary galaxy has no measurements to compare with, and the
       anatomy sliders drive every output through PhysicsEngine. */
    // SPARC Table 1's Hubble type codes (its note 1), each with the family and a plain description of the shape.
    // SPARC Table 1's Hubble type codes (its note 1), each with the family and how its arms look. Whether there is a
    // bulge comes from the galaxy's own SPARC mass model (any non-zero V_bulge), not from the type, so the card never
    // claims a bulge the model leaves out (NGC 3198 is type Sc but SPARC fits it with no bulge).
    const HUBBLE_TYPES=[['S0','Lenticular','No spiral arms'],['Sa','Spiral','Tightly wound arms'],
      ['Sab','Spiral','Fairly tight arms'],['Sb','Spiral','Medium arms'],['Sbc','Spiral','Fairly open arms'],
      ['Sc','Spiral','Loose arms'],['Scd','Spiral','Very loose arms'],['Sd','Spiral','Loose, patchy arms'],
      ['Sdm','Late-type spiral','Faint, patchy arms'],['Sm','Magellanic spiral','Barely-there arms'],
      ['Im','Irregular','No clear shape'],['BCD','Blue compact dwarf','Small, bright and forming stars']];
    const LIGHT_YEARS_PER_MPC=3.2616e6;
    const GAS_PER_HI=1.33; // SPARC counts helium by scaling the hydrogen mass (Lelli, McGaugh & Schombert 2016)
    const hubbleOf=galaxy=>HUBBLE_TYPES[galaxy.catalogue.hubble_type];
    /* Galaxy menu: a button and a listbox, so each option can carry its type badge and a one-line description.
       Arrow keys, Home/End, Enter/Space and Escape work as in a native menu. */
    const galaxyMenu=document.getElementById('galaxyMenu'),galaxyMenuBtn=document.getElementById('galaxyMenuBtn'),
      galaxyMenuList=document.getElementById('galaxyMenuList');
    const galaxyOptions=Object.values(GALAXIES).map(galaxy=>{
      const li=document.createElement('li'),code=document.createElement('span'),text=document.createElement('span'),
        name=document.createElement('b'),blurb=document.createElement('small');
      li.id=`galaxyOption-${galaxy.id}`;li.setAttribute('role','option');li.dataset.galaxy=galaxy.id;
      code.className='gm-code';code.textContent=hubbleOf(galaxy)[0];
      text.className='gm-text';name.textContent=galaxy.name;blurb.textContent=galaxy.summary;
      text.append(name,blurb);li.append(code,text);galaxyMenuList.append(li);
      return li;
    });
    let galaxyMenuActive=0;
    function setGalaxyMenuActive(index){
      galaxyMenuActive=(index+galaxyOptions.length)%galaxyOptions.length;
      galaxyOptions.forEach((li,i)=>li.classList.toggle('active',i===galaxyMenuActive));
      galaxyMenuList.setAttribute('aria-activedescendant',galaxyOptions[galaxyMenuActive].id);
    }
    function openGalaxyMenu(){
      galaxyMenuList.hidden=false;galaxyMenuBtn.setAttribute('aria-expanded','true');
      // Below the button by default; nudged left if that would run off a narrow screen.
      galaxyMenuList.style.left='';
      const overflow=galaxyMenuList.getBoundingClientRect().right-(document.documentElement.clientWidth-16);
      if(overflow>0)galaxyMenuList.style.left=`${-overflow}px`;
      setGalaxyMenuActive(galaxyOptions.findIndex(li=>li.dataset.galaxy===activeGalaxy.id));
      galaxyMenuList.focus();
    }
    function closeGalaxyMenu(returnFocus){
      if(galaxyMenuList.hidden)return;
      galaxyMenuList.hidden=true;galaxyMenuBtn.setAttribute('aria-expanded','false');
      if(returnFocus)galaxyMenuBtn.focus();
    }
    function pickGalaxy(id){
      closeGalaxyMenu(true);
      if(id===activeGalaxy.id)return;
      selectGalaxy(id);
      completeTutorialAction('galaxies','New galaxy loaded: its own data, visible matter and starting halo. Its profile is under the menu.');
    }
    galaxyMenuBtn.addEventListener('click',()=>{if(galaxyMenuList.hidden)openGalaxyMenu();else closeGalaxyMenu(true);});
    galaxyMenuBtn.addEventListener('keydown',e=>{if(e.key==='ArrowDown'||e.key==='ArrowUp'){e.preventDefault();openGalaxyMenu();}});
    galaxyMenuList.addEventListener('keydown',e=>{
      if(e.key==='ArrowDown'||e.key==='ArrowUp'){e.preventDefault();setGalaxyMenuActive(galaxyMenuActive+(e.key==='ArrowDown'?1:-1));}
      else if(e.key==='Home'||e.key==='End'){e.preventDefault();setGalaxyMenuActive(e.key==='Home'?0:-1);}
      else if(e.key==='Enter'||e.key===' '){e.preventDefault();pickGalaxy(galaxyOptions[galaxyMenuActive].dataset.galaxy);}
      else if(e.key==='Escape'){e.preventDefault();e.stopPropagation();closeGalaxyMenu(true);}
      else if(e.key==='Tab')closeGalaxyMenu(false);
    });
    galaxyMenuList.addEventListener('click',e=>{const li=e.target.closest('[role=option]');if(li)pickGalaxy(li.dataset.galaxy);});
    galaxyMenuList.addEventListener('pointermove',e=>{const li=e.target.closest('[role=option]');if(li)setGalaxyMenuActive(galaxyOptions.indexOf(li));});
    document.addEventListener('click',e=>{if(!galaxyMenu.contains(e.target))closeGalaxyMenu(false);});
    // Profile card: SPARC Table 1's catalogue values, with everyday units beside them.
    function renderGalaxyProfile(galaxy){
      const c=galaxy.catalogue,[code,family,arms]=hubbleOf(galaxy),radii=galaxy.rotation.radius_kpc;
      const shape=`${arms}, ${galaxy.rotation.v_bulge_kms.some(v=>v!==0)?'with a central bulge':'no central bulge'}`;
      const stars=DATA.upsilon_disk*c.luminosity_36_1e9_lsun,gas=GAS_PER_HI*c.hi_mass_1e9_msun;
      const billions=v=>v>=10?v.toFixed(0):v.toFixed(1),set=(id,text)=>{document.getElementById(id).textContent=text;};
      set('gpCode',code);set('gpFamily',`${family} · Hubble type ${code}`);set('gpShape',shape);set('gpSummary',galaxy.summary);
      set('gpDistance',`${c.distance_mpc} Mpc`);set('gpDistanceLy',`${Math.round(c.distance_mpc*LIGHT_YEARS_PER_MPC/1e6)} million light-years`);
      set('gpSpeed',`${c.v_flat_kms.toFixed(0)} km/s`);set('gpSpeedKmh',`${(Math.round(c.v_flat_kms*3.6)*1000).toLocaleString('en-GB')} km/h`);
      set('gpExtent',`${radii[radii.length-1].toFixed(0)} kpc`);set('gpPoints',`${radii.length} measured speeds`);
      set('gpStars',billions(stars));set('gpGas',billions(gas));
      document.getElementById('gpStarsBar').style.flexBasis=`${(100*stars/(stars+gas)).toFixed(1)}%`;
      document.getElementById('gpBar').setAttribute('aria-label',`Stars ${Math.round(100*stars/(stars+gas))}% and gas ${Math.round(100*gas/(stars+gas))}% of the visible mass`);
    }
    function applyGalaxyLabels(galaxy){
      document.querySelectorAll('[data-galaxy-name]').forEach(el=>{el.textContent=galaxy.name;});
      document.getElementById('galaxyMenuCode').textContent=hubbleOf(galaxy)[0];
      document.getElementById('galaxyMenuName').textContent=galaxy.name;
      galaxyOptions.forEach(li=>li.setAttribute('aria-selected',String(li.dataset.galaxy===galaxy.id)));
      renderGalaxyProfile(galaxy);
      orbitPlot.setAttribute('aria-label',`Interactive three-dimensional ${galaxyMode==='custom'?'custom galaxy':galaxy.name} star trajectories`);
      document.getElementById('curvePlot').setAttribute('aria-label',galaxyMode==='custom'?'Custom galaxy baryonic, halo and model rotation curves':`${galaxy.name} observed, baryonic, and fitted rotation curves`);
      document.getElementById('viewReset').title=`Return to the default ${galaxyMode==='custom'?'custom galaxy':galaxy.name} view`;
    }
    applyGalaxyLabels(activeGalaxy);
    // Wording that differs between a preset galaxy and the sandbox.
    function applyModeLabels(){
      const custom=customActive,fitButton=document.getElementById('applyObservedFit');
      document.querySelectorAll('[data-preset-only]').forEach(el=>{el.hidden=custom;});
      document.getElementById('massHelp').textContent=custom
        ?'How much dark matter surrounds your galaxy: 0 to 5 times the typical halo for these stars (Moster et al. 2013).'
        :'How much dark matter surrounds the galaxy. More mass, stronger gravity, faster orbits: the orange curve rises. 0 means visible matter only.';
      document.getElementById('scaleHelp').textContent=custom
        ?'How spread out the halo is. Small packs the dark matter near the centre. Starts at the typical value for this halo mass (Dutton & Macciò 2014).'
        :'How spread out the halo is, for the same mass. Small packs it near the centre and speeds up the inner stars; large spreads it out.';
      fitButton.textContent=custom?'Typical halo':'Fit halo to data';
      if(custom)fitButton.title='Put both halo sliders back to the typical halo for this mass of stars (Moster et al. 2013; Dutton & Macciò 2014)';else fitButton.removeAttribute('title');
      document.getElementById('footerModel').hidden=custom;
      document.getElementById('footerCustom').hidden=!custom;
      document.getElementById('componentKey').hidden=!custom;
      applyGalaxyLabels(activeGalaxy);
    }
    // Makes a SPARC preset the active galaxy: the page then follows its data, with stars on its observed orbits.
    function loadPreset(galaxy){
      customActive=false;customGridRmax=null;
      activeGalaxy=galaxy;
      DEFAULTS=haloDefaultsOf(galaxy);
      loadGalaxyData(galaxy);
      nfwFitApplied=true;applyBaryonState();
      configureHaloSliders();
      applyModeLabels();
      componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
      referenceSpeedSquared=presetReferenceSpeedSquared;
      const anatomy=presetAnatomy(galaxy);
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),PRESET_TEMPLATE_MASS));
      syncStarMeshes();updateStarTelemetry();controls.reset();
      // Restart the orbit clock and integrate the new population in the new galaxy's potential.
      simTime=0;syncTimelineUI(0);
      refreshOrbits();
      const{halo,total}=modelCurves(DEFAULTS.massT,DEFAULTS.scaleRadius);
      Plotly.react('curvePlot',buildCurveData(halo,total),curveLayoutNow(),curveConfig);
      updateHalo();clearTimeout(orbitModelTimer);
      document.dispatchEvent(new CustomEvent('galaxychange',{detail:{galaxy}}));
    }
    // Custom Sandbox: everything below is recomputed from the anatomy sliders, as fast as the sliders move.
    const COMPONENT_SQ=['diskSq','bulgeSq','gasSq']; // same order as COMPONENT_STYLE
    const snapScaleRadius=rs=>Math.min(60,Math.max(5,Math.round(rs/0.5)*0.5)); // the scale-radius slider's range and step
    function updateCustomSummary(baseline){
      document.getElementById('sbHaloBaseline').textContent=formatMass(baseline.m200);
      document.getElementById('sbStarCount').textContent=bulgeStarCount>0?`${starCount} (${starCount-bulgeStarCount} disk + ${bulgeStarCount} bulge)`:`${starCount} (all disk)`;
      document.getElementById('sbEjectRadius').textContent=`${ejectionRadius.toFixed(1)} kpc (${PM.ejection_factor} × R_max)`;
    }
    // V^2 of an exponential disk is proportional to its mass, so the costly per-radius tables are kept per unit mass and
    // only rebuilt when a disk's shape (R_d or R_max) changes. Nodes are exact near the centre, where V^2 changes fastest,
    // and every fourth node further out, interpolated linearly in between.
    const diskShape={key:'',stars:new Float64Array(ACCEL_TABLE_SIZE),gas:new Float64Array(ACCEL_TABLE_SIZE)};
    function fillCoarseTable(table,valueAt){
      const EXACT_NODES=64,STRIDE=4;
      let previous=0,previousValue=0;
      table[0]=0;
      for(let i=1;i<ACCEL_TABLE_SIZE;i++){
        if(i<EXACT_NODES||i%STRIDE===0||i===ACCEL_TABLE_SIZE-1){
          const value=valueAt(i*DR_TABLE);
          for(let j=previous+1;j<i;j++)table[j]=previousValue+(value-previousValue)*(j-previous)/(i-previous);
          table[i]=value;previous=i;previousValue=value;
        }
      }
    }
    // New test stars (and their 1 Gyr orbits, ~100 ms) are too slow to rebuild on every slider move, so a live rebuild
    // only records the anatomy here; refreshOrbits builds the stars from it once the slider pauses (scheduleOrbitUpdate).
    let pendingCustomStars=null;
    function rebuildCustomStars(){
      const{anatomy,baseline}=pendingCustomStars;pendingCustomStars=null;
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),sandboxTemplateMass()));
      syncStarMeshes();updateStarTelemetry();updateCustomSummary(baseline);
    }
    function rebuildCustomGalaxy({replot=false,live=false}={}){
      const anatomy={diskMass:sandbox.diskMass,diskScale:sandbox.diskScale,rMax:sandbox.truncation,bulgeMass:sandbox.bulgeMass,bulgeRadius:sandbox.bulgeRadius,gasMass:sandbox.gasMass};
      if(anatomy.rMax!==customGridRmax){loadCustomGrid(anatomy.rMax);customGridRmax=anatomy.rMax;replot=true;}
      // Rotation curve: the baryonic speed of the new anatomy on the force-table grid and at the data radii.
      customModel=PhysicsEngine.buildGalaxy(anatomy,PM.softening_kpc[0]);
      const shapeKey=`${anatomy.diskScale}|${anatomy.rMax}`;
      if(shapeKey!==diskShape.key){
        fillCoarseTable(diskShape.stars,r=>customModel.diskSq(r)/anatomy.diskMass);
        fillCoarseTable(diskShape.gas,r=>customModel.gasSq(r)/anatomy.gasMass);
        diskShape.key=shapeKey;
      }
      if(!customBaryonTable)customBaryonTable=new Float64Array(ACCEL_TABLE_SIZE);
      for(let i=0;i<ACCEL_TABLE_SIZE;i++)customBaryonTable[i]=diskShape.stars[i]*anatomy.diskMass+diskShape.gas[i]*anatomy.gasMass+customModel.bulgeSq(i*DR_TABLE);
      applyBaryonState();
      componentCurves=COMPONENT_SQ.map(key=>radii.map(r=>Math.sqrt(customModel[key](r))));
      // Halo: the stellar-to-halo mass relation sets the baseline; haloFine and rsManual hold the user's tuning until the anatomy changes.
      const baseline=PhysicsEngine.baselineHalo(anatomy.diskMass+anatomy.bulgeMass);
      DEFAULTS={massT:baseline.m200/1e12,scaleRadius:snapScaleRadius(baseline.scaleRadius)};
      configureHaloSliders();
      // Stars start on the orbits of the baryons plus baseline halo, then feel whatever halo the sliders leave.
      const baselineHalo=nfwVelocityModel(baseline.m200,DEFAULTS.scaleRadius);
      referenceSpeedSquared=r=>{const vb=baryonicVelocity(r),vh=baselineHalo(r);return vb*vb+vh*vh;};
      pendingCustomStars={anatomy,baseline};
      if(replot){
        const{halo,total}=modelCurves(Number(massSlider.value),Number(scaleSlider.value));
        Plotly.react('curvePlot',buildCurveData(halo,total),curveLayoutNow(),curveConfig);
      }
      if(live){updateCustomSummary(baseline);updateHalo(true);return;} // updateHalo schedules the stars and orbits
      refreshOrbits();
      updateHalo();clearTimeout(orbitModelTimer);
    }
    function activateCustom(){
      customActive=true;haloFine=1;rsManual=false;
      applyModeLabels();
      simTime=0;syncTimelineUI(0);
      rebuildCustomGalaxy({replot:true});
      controls.reset();
      document.dispatchEvent(new CustomEvent('galaxychange',{detail:{galaxy:activeGalaxy}}));
    }
    let customUpdateFrame=0;
    document.addEventListener('sandboxchange',()=>{
      if(!customActive)return;
      haloFine=1;rsManual=false; // a new galaxy gets its own typical halo: earlier halo tuning is dropped
      if(!customUpdateFrame)customUpdateFrame=requestAnimationFrame(()=>{customUpdateFrame=0;rebuildCustomGalaxy({live:true});});
    });
    document.getElementById('sandboxReset').addEventListener('click',()=>{haloFine=1;rsManual=false;});
    function selectGalaxy(id){
      const galaxy=GALAXIES[id];
      if(!galaxy||galaxy===activeGalaxy||customActive)return;
      loadPreset(galaxy);
    }

    /* ---------- Scorecard ----------
       Every preset galaxy scored under each model with the page's own functions. Standard values: visible matter only,
       the halo fit (halo only tuned) and plain MOND. Best realistic fits: the NFW and MOND best fits (realistic limits,
       from the generator's fit_nfw and fit_mond). The numbers match the fit badge. */
    function scoreGalaxy(galaxy){
      const rows=sparcRowsOf(galaxy),vBar=sparcBaryons(rows,DATA.upsilon_disk,DATA.upsilon_bulge);
      const chiOf=model=>rows.reduce((s,p,i)=>s+((p.observed_kms-model[i])/p.error_kms)**2,0);
      const{nfw_standard:std,nfw,mond}=galaxy.best_fit,nfwFrame=frameOf(galaxy,nfw),mondFrame=frameOf(galaxy,mond);
      const haloStd=nfwVelocityModel(std.m200_msun,std.scale_radius_kpc);
      const vBarNfw=sparcBaryons(rows,nfw.upsilon_disk,nfw.upsilon_bulge,nfwFrame),vBarMond=sparcBaryons(rows,mond.upsilon_disk,mond.upsilon_bulge,mondFrame);
      const halo=nfwVelocityModel(nfw.m200_msun,nfw.scale_radius_kpc,nfwFrame);
      return{pointCount:rows.length,
        k:{visible:0,nfwStd:tunedCount('nfw',std),mond:0,nfw:tunedCount('nfw',nfw),mondBest:tunedCount('mondBest',mond)},
        chi:{
        visible:chiOf(vBar),
        nfwStd:chiOf(rows.map((p,i)=>Math.hypot(vBar[i],haloStd(p.radius_kpc)))),
        nfw:chiOf(rows.map((p,i)=>Math.hypot(vBarNfw[i],halo(p.radius_kpc)))),
        mond:chiOf(rows.map((p,i)=>mondVelocity(vBar[i],p.radius_kpc,MOND_A0_DEFAULT))),
        mondBest:chiOf(rows.map((p,i)=>mondModelVelocity(vBarMond[i],p.radius_kpc,mond.a0_ms2,mondFrame)))}};
    }
    // The galaxy settings a best fit chose, beside SPARC's own values.
    // A fitted value that stopped at one of its realistic limits (fit.at_limit, recorded by the generator).
    function limitMark(fit,key){
      const side=fit.at_limit&&fit.at_limit[key];
      return side?` <abbr class="at-limit" title="At its ${side==='low'?'lower':'upper'} realistic limit: a closer fit would need a value outside the realistic range">at limit</abbr>`:'';
    }
    function fitSettingsHtml(galaxy,fit){
      const c=galaxy.catalogue,bulge=fit.tuned_settings>3?` · bulge M/L ${fit.upsilon_bulge.toFixed(2)}${limitMark(fit,'upsilon_bulge')}`:'';
      return `M/L ${fit.upsilon_disk.toFixed(2)}${limitMark(fit,'upsilon_disk')}${bulge} · ${fit.distance_mpc.toFixed(1)} Mpc${limitMark(fit,'distance_mpc')} (SPARC ${c.distance_mpc}) · tilt ${fit.inclination_deg.toFixed(0)}°${limitMark(fit,'inclination_deg')} (SPARC ${c.inclination_deg.toFixed(0)}°)`;
    }
    const SCORE_COLUMNS={visible:'Visible matter only',nfwStd:'Dark matter halo (NFW), standard values',mond:'MOND, standard values',nfw:'NFW, best realistic fit',mondBest:'MOND, best realistic fit'};
    function scoreRowHtml(galaxy){
      const score=scoreGalaxy(galaxy),{nfw_standard:std,nfw,mond}=galaxy.best_fit;
      const haloText=fit=>`M<sub>200</sub> ${(fit.m200_msun/1e12).toFixed(2)}×10¹² M☉${limitMark(fit,'m200_msun')} · concentration ${fit.concentration.toFixed(1)}${limitMark(fit,'concentration')}${limitMark(fit,'scale_radius_kpc')}`;
      const stdTuned=`<small>${score.k.nfwStd} tuned: ${haloText(std)}</small>`;
      const cell=(model,note='')=>{ // data-label names the column on phones, where the table stacks
        const red=reducedChi(score.chi[model],score.pointCount,score.k[model]),verdict=fitVerdict(red);
        const miss=Math.sqrt(score.chi[model]/score.pointCount); // root-mean-square gap, in error bars
        return `<td data-label="${SCORE_COLUMNS[model]}"><span class="score ${verdict.cls}" title="χ² = ${score.chi[model].toFixed(0)} over ${score.pointCount} measured speeds, ${score.k[model]} values tuned">${red<100?red.toFixed(2):red.toFixed(0)}</span><small>${verdict.word} · about ${miss<10?miss.toFixed(1):miss.toFixed(0)} error bars off</small>${note}</td>`;
      };
      const nfwTuned=`<small>${score.k.nfw} tuned: ${haloText(nfw)} · ${fitSettingsHtml(galaxy,nfw)}</small>`;
      const mondTuned=`<small>${score.k.mondBest} tuned: a₀ ${(mond.a0_ms2/1e-10).toFixed(2)}×10⁻¹⁰${limitMark(mond,'a0_ms2')} (standard ${(MOND_A0_DEFAULT/1e-10).toFixed(2)}) · ${fitSettingsHtml(galaxy,mond)}</small>`;
      return `<tr><th scope="row">${galaxy.name}<small>${score.pointCount} measured speeds</small><button type="button" class="score-show" data-galaxy="${galaxy.id}">Show in simulator</button></th>${cell('visible')}${cell('nfwStd',stdTuned)}${cell('mond')}${cell('nfw',nfwTuned)}${cell('mondBest',mondTuned)}</tr>`;
    }
    // At a glance: one verdict per model (standard values → best realistic fit), from the same scores as the full table.
    function glanceRowHtml(galaxy){
      const score=scoreGalaxy(galaxy),{nfw,mond}=galaxy.best_fit;
      const pill=model=>{const v=fitVerdict(reducedChi(score.chi[model],score.pointCount,score.k[model]));return `<span class="score ${v.cls}">${v.word.replace(' fit','')}</span>`;};
      const limited=fit=>fit.at_limit&&Object.values(fit.at_limit).some(Boolean)?'<sup title="Needed a value at its realistic limit">†</sup>':'';
      const pair=(std,best,fit)=>`<span class="glance-pair">${pill(std)}<span class="glance-arrow" aria-label="then, with the best realistic fit,">→</span>${pill(best)}${limited(fit)}</span>`;
      return `<tr><th scope="row"><button type="button" class="score-show glance-name" data-galaxy="${galaxy.id}" title="Show in simulator">${galaxy.name}</button></th><td data-label="${SCORE_COLUMNS.visible}">${pill('visible')}</td><td data-label="Dark matter halo (NFW): standard → best fit">${pair('nfwStd','nfw',nfw)}</td><td data-label="MOND: standard → best fit">${pair('mond','mondBest',mond)}</td></tr>`;
    }
    const scoreBody=document.getElementById('scoreBody');
    function renderScorecard(){
      const galaxies=Object.values(GALAXIES);
      document.getElementById('scoreCount').textContent=String(galaxies.length);
      scoreBody.innerHTML=galaxies.map(scoreRowHtml).join('');
      document.getElementById('glanceBody').innerHTML=galaxies.map(glanceRowHtml).join('');
    }
    setTimeout(renderScorecard,400);
    document.querySelector('.scorecard').addEventListener('click',e=>{
      const button=e.target.closest('.score-show');
      if(!button)return;
      if(customActive)setGalaxyMode('preset');
      selectGalaxy(button.dataset.galaxy);
      showPage('simulator',{scroll:false});
      document.getElementById('results').scrollIntoView({behavior:'instant'});
    });

    /* ---------- The puzzle figure ----------
       The default galaxy's real SPARC speeds against the speed its visible matter predicts, drawn once as an SVG.
       The colours match the rotation-curve plot; the shaded gap is a visual aid only. */
    function drawPuzzleFigure(){
      const galaxy=GALAXIES[DATA.default_galaxy_id],rows=sparcRowsOf(galaxy),vBar=sparcBaryons(rows,DATA.upsilon_disk,DATA.upsilon_bulge);
      const n=rows.length,outer=rows[n-1],W=400,H=250,LEFT=44,RIGHT=12,TOP=14,BOTTOM=38;
      const rMax=outer.radius_kpc*1.04,vMax=Math.max(...rows.map(p=>p.observed_kms+p.error_kms))*1.14;
      const x=r=>LEFT+(W-LEFT-RIGHT)*r/rMax,y=v=>H-BOTTOM-(H-TOP-BOTTOM)*v/vMax,f=value=>value.toFixed(1);
      const tickStep=(max,target)=>{const raw=max/target,pow=10**Math.floor(Math.log10(raw));return [1,2,5,10].map(m=>m*pow).find(step=>step>=raw);};
      const ticks=(max,target)=>{const step=tickStep(max,target),out=[];for(let t=0;t<=max;t+=step)out.push(t);return out;};
      const axis='#96a5ba',white='#f4f7fb',cyan='#5fd0e0',gold='#e9c88a';
      let svg=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Rotation curve of ${galaxy.name}: measured speeds stay near ${Math.round(outer.observed_kms)} km/s out to ${Math.round(outer.radius_kpc)} kpc, while the speed expected from visible matter falls to ${Math.round(vBar[n-1])} km/s">`;
      ticks(vMax,4).forEach(v=>{svg+=`<path d="M${LEFT} ${f(y(v))}H${W-RIGHT}" stroke="${axis}" stroke-opacity="${v?0.14:0.5}"/><text x="${LEFT-6}" y="${f(y(v)+4)}" text-anchor="end" fill="${axis}">${v}</text>`;});
      ticks(rMax,5).forEach(r=>{svg+=`<path d="M${f(x(r))} ${H-BOTTOM}v4" stroke="${axis}" stroke-opacity=".5"/><text x="${f(x(r))}" y="${H-BOTTOM+17}" text-anchor="middle" fill="${axis}">${r}</text>`;});
      svg+=`<text x="${(LEFT+W-RIGHT)/2}" y="${H-4}" text-anchor="middle" fill="${axis}">Distance from the centre (kpc)</text><text x="12" y="${(TOP+H-BOTTOM)/2}" text-anchor="middle" fill="${axis}" transform="rotate(-90 12 ${(TOP+H-BOTTOM)/2})">Orbital speed (km/s)</text>`;
      // Shade the gap from the last radius where the visible matter still keeps up with the measured speed.
      let from=0;
      rows.forEach((p,i)=>{if(p.observed_kms<=vBar[i])from=Math.min(i+1,n-1);});
      const top=rows.slice(from).map(p=>`${f(x(p.radius_kpc))},${f(y(p.observed_kms))}`),bottom=rows.slice(from).map((p,i)=>`${f(x(p.radius_kpc))},${f(y(vBar[from+i]))}`).reverse();
      svg+=`<polygon points="${top.concat(bottom).join(' ')}" fill="${gold}" fill-opacity=".13"/>`;
      svg+=`<polyline points="${rows.map((p,i)=>`${f(x(p.radius_kpc))},${f(y(vBar[i]))}`).join(' ')}" fill="none" stroke="${cyan}" stroke-width="2" stroke-dasharray="5 4"/>`;
      rows.forEach(p=>{svg+=`<path d="M${f(x(p.radius_kpc))} ${f(y(p.observed_kms-p.error_kms))}V${f(y(p.observed_kms+p.error_kms))}" stroke="${white}" stroke-opacity=".55"/><circle cx="${f(x(p.radius_kpc))}" cy="${f(y(p.observed_kms))}" r="2.4" fill="${white}"/>`;});
      const mid=rows[Math.max(from,Math.round(0.7*(n-1)))],midBar=vBar[rows.indexOf(mid)];
      svg+=`<text x="${f(x(mid.radius_kpc))}" y="${f(y((mid.observed_kms+midBar)/2)+4)}" text-anchor="middle" fill="${gold}" font-weight="600">the gap</text>`;
      svg+=`<text x="${W-RIGHT}" y="${f(y(outer.observed_kms+outer.error_kms)-9)}" text-anchor="end" fill="${white}">measured</text>`;
      svg+=`<text x="${W-RIGHT}" y="${f(y(vBar[n-1])+17)}" text-anchor="end" fill="${cyan}">expected from visible matter</text></svg>`;
      document.getElementById('puzzleChart').innerHTML=svg;
      const ratio=(outer.observed_kms/vBar[n-1])**2;
      document.getElementById('puzzleObs').textContent=`${Math.round(outer.observed_kms)} km/s`;
      document.getElementById('puzzleBar').textContent=`${Math.round(vBar[n-1])} km/s`;
      document.getElementById('puzzleRatio').textContent=ratio<3?ratio.toFixed(1):String(Math.round(ratio));
    }
    drawPuzzleFigure();

    /* ---------- Settings panel ----------
       The gear opens it; a click outside, Escape or the gear again closes it. The settings themselves are wired where they act
       (High-quality visuals: beside the 3D view's resolution levels). */
    const settingsBtn=document.getElementById('settingsBtn'),settingsPanel=document.getElementById('settingsPanel');
    function setSettingsOpen(open){
      settingsPanel.hidden=!open;settingsBtn.setAttribute('aria-expanded',String(open));
      document.querySelector('.topbar').classList.toggle('settings-open',open); // lifts the menu over the tour's rings and dim
      if(open)(settingsPanel.querySelector('input:not(:disabled)')||settingsPanel).focus();
    }
    settingsBtn.addEventListener('click',()=>setSettingsOpen(settingsPanel.hidden));
    document.addEventListener('pointerdown',e=>{if(!settingsPanel.hidden&&!e.target.closest('.settings'))setSettingsOpen(false);});
    document.addEventListener('keydown',e=>{if(e.key==='Escape'&&!settingsPanel.hidden){setSettingsOpen(false);settingsBtn.focus();}});

    /* ---------- Pages ----------
       The site is one file split into pages. The top menu, or any #link, shows one .page and hides the rest; a #link
       to something inside a page (such as #curvePanel) opens that page and scrolls to it. */
    const pageEls=[...document.querySelectorAll('.page')];
    let activePage='puzzle';
    const pageOf=el=>{const page=el&&el.closest?el.closest('.page'):null;return page?page.dataset.page:null;};
    const pageLink=page=>document.querySelector(`.topnav a[href="#${page}"]`);
    function showPage(page,{scroll=true}={}){
      if(!pageEls.some(el=>el.dataset.page===page))return;
      const changed=page!==activePage;
      activePage=page;
      pageEls.forEach(el=>el.classList.toggle('active',el.dataset.page===page));
      document.querySelectorAll('.topnav a').forEach(a=>{if(a.getAttribute('href')==='#'+page)a.setAttribute('aria-current','page');else a.removeAttribute('aria-current');});
      const link=pageLink(page);
      if(link)link.scrollIntoView({block:'nearest',inline:'nearest',behavior:'instant'}); // keeps the menu entry in view on phones
      if(scroll)window.scrollTo({top:0,behavior:'instant'});
      if(changed)document.dispatchEvent(new CustomEvent('pagechange',{detail:{page}}));
    }
    function followHash(hash){
      const id=String(hash).replace(/^#/,'');
      if(!id)return false;
      if(pageEls.some(el=>el.dataset.page===id)){showPage(id);return true;}
      const target=document.getElementById(id),page=pageOf(target);
      if(!page)return false;
      showPage(page);target.scrollIntoView({behavior:'instant'});
      return true;
    }
    document.addEventListener('click',e=>{
      const link=e.target.closest?e.target.closest('a[href^="#"]'):null;
      if(!link||!followHash(link.getAttribute('href')))return;
      e.preventDefault();
      history.pushState(null,'',link.getAttribute('href'));
    });
    window.addEventListener('popstate',()=>{if(!followHash(location.hash))showPage('puzzle');});

    const timeSlider=document.getElementById('timeSlider');
    function setPlaying(isPlaying){
        playing=isPlaying;
        const b=document.getElementById('playPauseBtn');
        if(b){
            b.setAttribute('aria-pressed',String(playing));
            b.setAttribute('aria-label',playing?'Pause simulation':'Play simulation');
            b.querySelector('.play-icon').textContent=playing?'❚❚':'▶';
            b.querySelector('.play-text').textContent=playing?'Pause':'Play';
        }
        const fb=document.getElementById('fsPlay');
        fb.textContent=playing?'❚❚':'▶';fb.setAttribute('aria-pressed',String(playing));fb.setAttribute('aria-label',playing?'Pause simulation':'Play simulation');
    }
    document.getElementById('playPauseBtn').addEventListener('click',()=>togglePlayback());
    // Fullscreen hides the timeline, so the corner panel has its own Play button and clock.
    document.getElementById('fsPlay').addEventListener('click',()=>togglePlayback());

    function togglePlayback(){
        if(playing){
            setPlaying(false);
        }else{
            if(simTime>=SIM_TOTAL_TIME_MYR){
                simTime=0.0;
                updateStarTransforms(simTime);
                updateEscapeCounter(simTime);
                syncTimelineUI(simTime);
            }
            setPlaying(true);
        }
        requestRender();
    }
    // After a mouse/touch click, drop focus from the button so Space / arrow keys act globally instead of re-activating it.
    document.addEventListener('click',e=>{
        if(e.detail===0)return;
        const btn=e.target instanceof Element?e.target.closest('button'):null;
        if(btn&&!btn.disabled)setTimeout(()=>{if(document.activeElement===btn)btn.blur();},0);
    });
    // Same for sliders and toggle switches once the mouse/finger lets go, so the arrow keys go back to the tour
    // (or the fullscreen time step) instead of nudging the last slider touched. Keyboard focus via Tab is untouched.
    document.addEventListener('pointerup',()=>setTimeout(()=>{
        const a=document.activeElement;
        if(a&&a.matches&&a.matches('input[type=range],input[type=checkbox],input[type=radio]'))a.blur();
    },0));
    window.addEventListener('keydown',e=>{
        if(e.code==='Escape'&&orbitPlot.classList.contains('fullscreen-active')){setCssFullscreen(false);return;}
        // Space plays/pauses even when a 3D-view button or the tour card has focus (e.g. just after pressing Next or
        // Orbits), so the tour's 'press Space' works; elsewhere it keeps its normal job of pressing the focused button.
        const t=e.target&&e.target.closest?e.target:null;
        const busy=t&&(t.closest('input:not([type=range]),select,textarea,dialog,[contenteditable="true"]')||(t.closest('button,a,summary')&&!t.closest('#orbitPlot,#tutorialCard')));
        if(e.code==='Space'&&!busy&&!document.getElementById('fsPlay').inert){
            e.preventDefault();
            if(!e.repeat)togglePlayback();
        }
    });

    function syncTimelineUI(timeMyr){
        const val=Math.round(Math.max(0,Math.min(SIM_TOTAL_TIME_MYR,timeMyr)));
        timeSlider.value=val;
        document.getElementById('timeValue').textContent=`${val} Myr / 1000 Myr`;
        document.getElementById('fsTimeValue').textContent=`${val} / 1000 Myr`;
    }

    function onSliderInput(val){
        simTime=Math.max(0,Math.min(SIM_TOTAL_TIME_MYR,Number(val)));
        updateStarTransforms(simTime);
        updateEscapeCounter(simTime);
        syncTimelineUI(simTime);
        requestRender();
    }
    timeSlider.addEventListener('input',()=>onSliderInput(timeSlider.value));
    updateStarTransforms(0);updateEscapeCounter(0);syncTimelineUI(0);

    /* ---------- Star marks (tutorial aid) ----------
       Rings that follow chosen test stars across the 3D view with a live label, so the tour can point out things
       that are easy to miss: an inner star lapping an outer one, a star bobbing through the disk, the first star to
       leave. They only read the trajectories already computed; they never change the simulation. */
    const KPC_PER_MYR_TO_KMS=977.8;
    const starMarkEls=[0,1,2].map(()=>{
      const el=document.createElement('div');el.className='star-mark';el.hidden=true;
      el.append(document.createElement('span'));orbitPlot.append(el);return el;
    });
    let starMarks=[];
    const markVec=new THREE.Vector3();
    function starPositionAt(idx,timeMyr){
      const t=activeTrajectories[idx],f=Math.max(0,Math.min(frameCount-1,timeMyr/SIM_DT_MYR)),f0=Math.floor(f),f1=Math.min(f0+1,frameCount-1),a=f-f0;
      return [0,1,2].map(c=>t[f0*3+c]+(t[f1*3+c]-t[f0*3+c])*a);
    }
    function starSpeedAt(idx,timeMyr){
      const t=activeTrajectories[idx],f=Math.min(frameCount-2,Math.floor(Math.max(0,timeMyr)/SIM_DT_MYR));
      return Math.hypot(t[f*3+3]-t[f*3],t[f*3+4]-t[f*3+1],t[f*3+5]-t[f*3+2])/SIM_DT_MYR*KPC_PER_MYR_TO_KMS;
    }
    // Laps round the centre between fromMyr and timeMyr.
    function lapsAt(idx,timeMyr,fromMyr=0){
      const t=activeTrajectories[idx],last=Math.min(frameCount-1,Math.floor(Math.max(0,timeMyr)/SIM_DT_MYR));
      let angle=0;
      for(let f=Math.floor(Math.max(0,fromMyr)/SIM_DT_MYR)+1;f<=last;f++){
        let d=Math.atan2(t[f*3+1],t[f*3])-Math.atan2(t[f*3-2],t[f*3-3]);
        if(d>Math.PI)d-=2*Math.PI;else if(d<-Math.PI)d+=2*Math.PI;
        angle+=d;
      }
      return Math.abs(angle)/(2*Math.PI);
    }
    // Stars are numbered by radius, so bulge and disk stars are mixed in the list: filter by kind.
    function diskStarNear(radius,tilted){
      let best=-1,score=Infinity;
      for(let k=0;k<starCount;k++){
        if(particles[k].kind!=='disk')continue;
        const p=particles[k].position,R=Math.hypot(p[0],p[1]),off=Math.abs(R-radius)/radius;
        const s=tilted?(off<0.25?-Math.abs(p[2]):Infinity):off+Math.abs(p[2]);
        if(s<score){score=s;best=k;}
      }
      return best;
    }
    function firstEscaper(){
      let best=-1,time=Infinity;
      for(let k=0;k<starCount;k++)if(escapeTimeByStar[k]<time){time=escapeTimeByStar[k];best=k;}
      return best;
    }
    function bulgeStar(){
      if(!bulgeStarCount)return -1;
      let best=0,score=Infinity;
      for(let k=0;k<starCount;k++){if(particles[k].kind!=='bulge')continue;const s=Math.abs(particles[k].radius-1.5*currentAnatomy.bulgeRadius);if(s<score){score=s;best=k;}}
      return best;
    }
    // The race card counts laps from the moment the race starts (raceStartMyr), not from 0 Myr, so the counts never
    // give the answer away before the viewer has predicted it. Display only: it reads the computed orbits.
    let raceStartMyr=0;
    const lapText=(name,idx)=>`${name} · ${lapsAt(idx,simTime,raceStartMyr).toFixed(1)} laps · ${starSpeedAt(idx,simTime).toFixed(0)} km/s`;
    const STAR_MARKS={
      race:[
        {pick:()=>diskStarNear(0.1*currentAnatomy.rMax,false),label:idx=>lapText('Inner star',idx)},
        {pick:()=>diskStarNear(0.4*currentAnatomy.rMax,false),warm:true,label:idx=>lapText('Outer star',idx)}
      ],
      tilt:[{pick:()=>diskStarNear(0.4*currentAnatomy.rMax,true),label:idx=>{const z=starPositionAt(idx,simTime)[2];return `Height ${z>=0?'+':'−'}${Math.abs(z).toFixed(2)} kpc`;}}],
      escaper:[{pick:firstEscaper,warm:true,label:idx=>{
        const out=escapeTimeByStar[idx];
        return simTime<out?`First to leave · ${Math.round(out)} Myr`:`Ejected · ${Math.hypot(...starPositionAt(idx,simTime)).toFixed(0)} kpc out`;
      }}],
      bulge:[{pick:bulgeStar,warm:true,label:()=>'Bulge star'}]
    };
    // Star tracker (main view): the viewer types a star number; it gets its own ring alongside any tutorial marks.
    let trackedStar=-1;
    const trackInputs=['trackInput','fsTrack'].map(id=>document.getElementById(id));
    const trackMark={track:true,pick:()=>trackedStar<starCount?trackedStar:-1,label:idx=>{
      const r=Math.hypot(...starPositionAt(idx,simTime));
      return `Star ${idx+1} · ${lapsAt(idx,simTime).toFixed(1)} laps · ${starSpeedAt(idx,simTime).toFixed(0)} km/s · ${r.toFixed(1)} kpc`;
    }};
    function syncTrackLimits(){
      trackInputs.forEach(el=>el.max=starCount);
      document.getElementById('trackMax').textContent=`1–${starCount}`;document.getElementById('trackLast').textContent=starCount;
    }
    trackInputs.forEach(el=>el.addEventListener('input',()=>{
      const n=parseInt(el.value,10);
      trackedStar=n>=1&&n<=starCount?n-1:-1;
      trackInputs.forEach(o=>{if(o!==el)o.value=el.value;});
      requestRender();
    }));
    trackInputs.forEach(el=>el.addEventListener('focus',syncTrackLimits));
    syncTrackLimits();
    function setStarMarks(names){
      starMarks=(names||[]).flatMap(name=>STAR_MARKS[name]);
      requestRender();
    }
    const STAR_LABEL_CHAR_PX=6.6,STAR_LABEL_PAD_PX=34; // rough label width per character, plus the ring and padding
    function updateStarMarks(){
      // Read the view size once: reading it between style writes would force a fresh page layout for every mark.
      const viewW=orbitPlot.clientWidth,viewH=orbitPlot.clientHeight;
      starMarkEls.forEach((el,i)=>{
        const marks=trackedStar>=0?[...starMarks.slice(0,2),trackMark]:starMarks;
        const mark=marks[i],idx=mark?mark.pick():-1;
        if(idx<0){el.hidden=true;return;}
        markVec.set(...starPositionAt(idx,simTime)).project(camera);
        if(markVec.z>1){el.hidden=true;return;}
        // A star that has left the screen stays marked at the edge, on the side it went.
        const edge=Math.max(Math.abs(markVec.x),Math.abs(markVec.y))>0.92;
        if(edge){const f=0.92/Math.max(Math.abs(markVec.x),Math.abs(markVec.y));markVec.x*=f;markVec.y*=f;}
        el.classList.toggle('edge',edge);
        // Keep the label on screen near the right edge: flip it to the left when its estimated width (from its length,
        // so nothing is measured mid-frame) would run past the edge, which on a phone happens well short of it.
        const text=mark.label(idx),px=(markVec.x+1)/2*viewW;
        el.classList.toggle('flip',markVec.x>0.4||px+text.length*STAR_LABEL_CHAR_PX+STAR_LABEL_PAD_PX>viewW);
        el.hidden=false;
        el.style.transform=`translate(${px.toFixed(1)}px,${((1-markVec.y)/2*viewH).toFixed(1)}px)`;
        el.classList.toggle('warm',!!mark.warm);el.classList.toggle('track',!!mark.track);
        if(el.firstChild.textContent!==text)el.firstChild.textContent=text;
      });
      tutorialWatchTick();
    }

    /* ---------- Guided tours ----------
       Two tours share one card: 'main' (the rotation-curve puzzle, the halo, MOND and the 3D view) and 'sandbox'
       (building a custom galaxy). Steps are looked up by id, never by position, so they can be added or reordered.
       The main ring sits on target; also lists other controls the card names, which get dashed rings.
       A step can need NFW mode (nfw), a halo mass above zero (halo) or one galaxy source (mode). When the page is
       not in that state, the card says how to get there and the pointer moves to the control that does it.
       The tours never change the simulator's settings themselves. While a card is open, every control it does not name
       is locked (see tutorialAllowed). A guided tour from the puzzle page runs the main tour, then the sandbox tour. */
    // Leave (~0.45 s) plus the staggered enter animation (~0.85 s): the card's scrollbar stays hidden for this long.
    const TUTORIAL_SWITCH_MS=1350;let tutorialScrollTimer=0;
    const TUTORIAL_RING_MIN_RADIUS_PX=8;  // softest corner any highlight ring gets
    const TUTORIAL_BESIDE_PX=16;          // an opened box this close to its control is drawn inside the control's ring
    const TUTORIAL_RING_FADE_MS=2650;     // a done part's green ring holds, then fades (matches tutorial-ring-done)
    const TUTORIAL_HANDOVER_MS=700;       // the next part's ring cross-fades into the main ring (matches tutorial-handover)
    const TUTORIAL_NEAR_PX=160;           // a ringed thing this close above the control is scrolled into view with it
    const TOURS={
      main:{label:'Simulator tutorial',steps:[
        {id:'puzzle',mode:'preset',target:'curvePlot',gesture:'scan',tag:'Move along the right half',title:'Spot the puzzle',
         checks:{scanNear:()=>tutorialFlags.scan0,scanMid:()=>tutorialFlags.scan1,scanFar:()=>tutorialFlags.scan2},
         doneText:'Spotted: all the way out, the white points stay well above the cyan line. Press Next to find out why.',
         text:'Each white point is a measured orbit speed. The cyan line is how fast they would move if only the stars and gas we can see were pulling on them.',
         action:'Slowly move your mouse (or slide your finger) across the whole right half of the graph, from the middle to the far edge, reading the speeds as you go. Do the white points stay above the cyan line?'},
        {id:'removeHalo',nfw:true,target:'massSlider',cardSlider:{massSlider:'curvePlot'},gesture:'zero',allow:['visibleMatterOnly'],watch:['curvePlot:Watch the orange line'],title:'Take the dark matter away',
         text:'That gap is the puzzle. The orange line is the model: visible matter plus a dark matter halo.',
         action:'Drag Halo mass (blue box) all the way left, to 0. Left lowers it.'},
        {id:'addHalo',nfw:true,mode:'preset',target:'massSlider',cardSlider:{massSlider:'curvePlot'},gesture:'increase',watch:['curvePlot:Watch the orange line'],title:'Add the halo back',
         checks:{massSlider:()=>Number(massSlider.value)>0&&tutorialOuterReached()},doneText:'Gap closed: enough dark matter lifts the orange line up to the measured speeds.',
         text:'A halo is a huge cloud of unseen matter. Its gravity makes the stars orbit faster.',
         action:'Drag Halo mass to the right until the orange line climbs up to the white points on the right of the graph.'},
        {id:'haloShape',nfw:true,target:'scaleSlider',cardSlider:{scaleSlider:'curvePlot'},gesture:'adjust',tag:'Drag both ways',watch:['curvePlot:Watch the left of the orange line'],title:'Pack the halo in or spread it out',
         text:'The scale radius sets how tightly the halo is packed. Its total mass stays the same.',
         action:'1. Drag Scale radius left: the left of the orange line rises. 2. Drag it right: that part drops.'},
        {id:'haloNumbers',nfw:true,halo:true,cardRight:true,scrollPad:130,target:'haloRadiusCard',gesture:'click',tag:'1 · Open',also:['haloScaleCard:2 · Open','haloConcentrationCard:3 · Open'],title:'Read the halo’s numbers',
         checks:{haloRadiusCard:()=>tutorialFlags.haloRadiusCard,haloScaleCard:()=>tutorialFlags.haloScaleCard,haloConcentrationCard:()=>tutorialFlags.haloConcentrationCard},doneText:'All three read: these numbers follow the halo sliders.',
         text:'These cards describe your halo. r₂₀₀ marks its outer edge. Concentration c (r₂₀₀ ÷ scale radius) says how tightly packed it is.',
         action:'1. Hover over (or tap) Halo radius. 2. Now Scale radius. 3. Now Concentration.'},
        {id:'fit',nfw:true,mode:'preset',target:'massSlider',cardSlider:{massSlider:'curvePlot',scaleSlider:'curvePlot'},gesture:'adjust',tag:'1 · Drag',also:['scaleSlider:2 · Drag'],watch:['chiBadge:χ² (lower is better)'],title:'Score your fit',
         checks:{massSlider:()=>Number(massSlider.value)!==tutorialStart.mass&&shownRedChi>=tutorialChiGoal(),scaleSlider:()=>Number(scaleSlider.value)!==tutorialStart.scale,
           chiGoal:()=>tutorialLatched.massSlider&&tutorialLatched.scaleSlider&&shownRedChi<tutorialChiGoal()},doneText:'Fitted by hand: the halo now matches the measured speeds.',
         text:'χ² (yellow box) scores the fit: it adds up the gaps between the orange line and the white points. Lower is better.',
         action:'1. Drag Halo mass until the yellow box no longer reads Close fit. 2. Drag Scale radius and watch χ². 3. Fine-tune both until it reads Close fit again.'},
        {id:'mond',target:'modelBar',gesture:'switch',tag:'Choose MOND',watch:['chiBadge:Compare χ²'],keepWatch:true,title:'Test a rival idea',
         text:'MOND has no dark matter. Instead, it changes the law of gravity, so it pulls harder than Newton’s law says: a little where gravity is strong, much more where it is weak. The deeper the gold on the graph, the bigger the boost.',
         action:'Click the MOND card in the blue box. Compare χ² in the yellow box.'},
        {id:'mondZone3d',mond:true,viewPanel:true,given:['toggleMondZone'],target:'toggleMondZone',gesture:'click',tag:'Switch on',title:'See where MOND acts',
         checks:{toggleMondZone:()=>mondMode&&mondZoneShown},doneText:'Faint gold means a small boost; deep gold means MOND is doing most of the work.',
         text:'The same gold as the graph, laid on the galaxy itself. Faint near the centre, where MOND’s boost is small, and deeper further out, always in proportion to the boost.',
         action:'Switch on Where MOND acts in the View panel.'},
        {id:'bestFit',mond:true,mode:'preset',target:'bestFitWrap',gesture:'click',tag:'Turn on',watch:['chiBadge:Watch χ²','mlStat:M/L','tiltStat:Tilt'],keepWatch:true,title:'Give both their best shot',
         text:'Nobody knows exactly how heavy a galaxy’s stars are for their brightness (the mass-to-light ratio), how far away it is, or how tilted it is. This switch lets each model adjust them within realistic limits.',
         action:'Turn on Best realistic fits (blue box). Watch χ², M/L and tilt change.'},
        // Watching the simulation: these steps run inside fullscreen (fs), where the tour card follows.
        {id:'fsEnter',given:['toggleFullscreen'],nfw:true,viewPanel:true,target:'toggleFullscreen',gesture:'click',tag:'1 · Fullscreen',also:['fsPlay:2 · Play'],freeCamera:true,title:'Now watch it happen',
         checks:{toggleFullscreen:tutorialIsFullscreen,fsPlay:()=>playing},doneText:'Fullscreen and playing: watch for as long as you like, then press Next.',
         text:'The 3D view moves hundreds of test stars under the model’s gravity for a billion years. Jerky on this computer? Untick High-quality visuals in Settings (the gear, top right).',
         action:'1. Press Fullscreen: the tour comes too. 2. Press ▶ Play (or Space).'},
        {id:'fsClean',given:['toggleGas','toggleHalo'],fs:true,viewPanel:true,target:'toggleGas',tag:'1 · Gas off',also:['toggleHalo:2 · Halo off'],gesture:'click',freeCamera:true,title:'Clear the stage',
         text:'The blue haze is gas. The purple cloud is the halo, drawn so you can see it: real dark matter is invisible.',
         action:'1. Switch off Gas. 2. Switch off Halo.',
         checks:{toggleGas:()=>!gasCloud.visible,toggleHalo:()=>haloHidden||mondMode||Number(massSlider.value)===0},doneText:'Clear view: everything left is a test star.'},
        {id:'fsOrbits',given:['toggleOrbits','fsPlay'],fs:true,viewPanel:true,target:'viewFaceOn',tag:'1 · Face-on',also:['toggleOrbits:2 · Orbits on','fsPlay:3 · Play'],gesture:'click',freeCamera:()=>!tutorialLatched.viewFaceOn,title:'Look down on the disk',
         text:'Face-on, you look straight down on the disk. Each line is one star’s path: near-perfect circles.',
         action:'1. Press Face-on. 2. Switch on Orbits. 3. Press ▶ Play.',
         checks:{viewFaceOn:()=>controls.goal.phi<0.3,toggleOrbits:()=>orbitLines.visible,fsPlay:()=>playing||tutorialPlayedMyr>0},doneText:'Circles everywhere: a galaxy in balance.'},
        {id:'fsRace',fs:true,nfw:true,marks:['race'],target:'fsPlay',gesture:'click',tag:'2 · Play',checks:{fsPlay:()=>tutorialFlags.prediction!=null&&playing},title:'Which star laps more?',
         predict:{question:'Predict first: over the next 400 million years, which ringed star will finish more laps?',options:['Inner star','Outer star','About the same'],answer:0},
         text:'The two ringed stars move at similar speeds, but the outer one is about four times further out.',
         action:'Press ▶ Play and watch the lap counts on the two ringed stars.',
         done:()=>tutorialFlags.prediction!=null&&tutorialPlayedMyr>=400,
         doneText:()=>tutorialFlags.prediction===0?'You called it: at similar speeds, the inner star’s much shorter lap means it laps far more often.':'The inner star wins. At similar speeds, the outer star’s lap (2πr) is about four times longer, so the inner star finishes several laps for each of the outer star’s.'},
        {id:'fsTilt',fs:true,viewPanel:true,marks:['tilt'],target:'viewEdgeOn',tag:'1 · Edge-on',also:['fsPlay:2 · Keep playing'],gesture:'click',title:'Turn it on its side',
         text:'Edge-on, the disk is a thin line. The ringed star bobs above and below it once a lap.',
         action:'1. Press Edge-on. 2. Press ▶ Play and watch the ringed star’s height.',
         checks:{viewEdgeOn:()=>Math.abs(controls.goal.phi-Math.PI/2)<0.2,fsPlay:()=>playing},doneText:'Up, down, up: one bob per lap.'},
        {id:'fsTrack',fs:true,target:'fsTrack',gesture:'click',tag:'Type a number',title:'Follow any star',
         checks:{fsTrack:()=>trackedStar>=0},doneText:'Tracked: its ring shows its laps, speed and distance.',
         text:'Every star has a number: 1 starts nearest the centre.',
         action:'Type a star number in the Star # box.'},
        {id:'fsEscape',given:['toggleEjection'],fs:true,nfw:true,viewPanel:true,marks:['escaper'],ring:true,target:'toggleEjection',tag:'1 · Ejection line on',also:['fsMassSlider:2 · Drag to 0','fsRestart:3 · Restart'],
         checks:{toggleEjection:()=>ejectionLineShown,fsMassSlider:()=>Number(massSlider.value)===0,fsRestart:()=>tutorialFlags.restart0},gesture:'click',watch:['escape-counter:Ejected count'],title:'Pull the halo away',
         text:'A star that crosses the red circle (1.5 × the galaxy’s visible edge) counts as ejected. Without the halo, the outer stars move too fast for the visible matter’s gravity to hold them.',
         action:'1. Switch on Ejection line. 2. Drag Halo mass to 0. 3. Press Restart and watch the ringed star.',
         done:prev=>{const k=firstEscaper();return Number(massSlider.value)===0&&k>=0&&prev<escapeTimeByStar[k]&&simTime>=escapeTimeByStar[k];},doneText:'Gone: visible matter alone cannot hold the outer stars.'},
        {id:'fsFullHalo',fs:true,nfw:true,halo:true,viewPanel:true,target:'fsRestart',gesture:'click',tag:'1 · Restart',also:['viewFullHalo:2 · Full halo'],title:'Put the halo back, and step back',
         checks:{fsRestart:()=>tutorialFlags.restartHalo,viewFullHalo:()=>tutorialFlags.fullHalo},doneText:'Zoomed out: the whole visible galaxy is that bright speck.',
         text:'With the halo back, the stars stay. The sphere marks the halo’s edge, r₂₀₀; the bright speck is the whole visible galaxy.',
         action:'1. Press Restart: the stars stay put. 2. Press Full halo to zoom out.'},
        {id:'galaxies',mode:'preset',target:'galaxyMenu',gesture:'click',tag:'Choose a galaxy',watch:['galaxyProfile:Its profile'],keepWatch:true,title:'Try another galaxy',
         text:'Each galaxy has its own measured speeds. The card in the yellow box says what kind of galaxy it is.',
         action:'Open the Galaxy menu (blue box) and pick another galaxy.'},
        {id:'ml',mode:'preset',cardRight:true,target:'mlInfoBtn',gesture:'click',tag:'Hover or tap',title:'How heavy are the stars?',
         text:'Telescopes measure light, not mass. The mass-to-light ratio converts one into the other, and it is an estimate.',
         action:'Hover over (or tap) the ? in the blue box.'},
        {id:'scorecard',closing:true,target:'glanceTable',gesture:'look',tag:'Read along the rows',title:'How both ideas score',
         text:'Each row is a galaxy, scored under each idea: close (blue), rough (orange) or poor (red). An arrow goes from the standard values to the best realistic fit; a † means that fit needed a value at the edge of what is realistic.',
         action:'Read along the rows in the blue box, then press Next.'},
        {id:'finish',closing:true,noRing:true,target:'theoryGuide',gesture:'scan',tag:'Read on',title:'What this does and doesn’t show',
         text:'A close fit does not prove either idea: both models are simplified and the measurements have errors.',
         action:'Read on below to learn more.'}
      ]},
      sandbox:{label:'Sandbox tutorial',steps:[
        {id:'sbIntro',mode:'custom',target:'sandboxGroups',gesture:'look',tag:'Your galaxy’s sliders',title:'Build your own galaxy',
         text:'Here you design the galaxy: how much its stars and gas weigh, and how they are spread out.',
         action:'Look at the three groups of sliders in the blue box, then press Next.'},
        {id:'sbDisk',mode:'custom',target:'sbDiskMass',cardSlider:{sbDiskMass:'curvePlot'},gesture:'adjust',tag:'Drag both ways',watch:['massValue:Halo mass','curvePlot:Watch the cyan line'],title:'Make the disk heavier or lighter',
         text:'The disk is the flat, spinning part where most stars are. Galaxies with more stars sit in bigger halos.',
         action:'Drag Disk mass up and down. Watch the cyan line and Halo mass.'},
        {id:'sbSize',mode:'custom',partPause:true,target:'sbDiskScale',cardSlider:{sbDiskScale:'curvePlot',sbTruncation:'orbitPlot'},tag:'1 · Drag',also:['sbTruncation:2 · Drag'],gesture:'adjust',
         checks:{sbDiskScale:()=>sbVal('sbDiskScale')!==tutorialStart.sbDiskScale,sbTruncation:()=>sbVal('sbTruncation')!==tutorialStart.sbTruncation},doneText:'Disk reshaped and edge moved.',title:'Spread the stars out',
         text:'Scale length sets how quickly the disk fades outwards. Truncation radius is where it ends.',
         action:'1. Drag Scale length and watch the peak of the cyan line. 2. Drag Truncation radius and watch the 3D galaxy’s size.'},
        {id:'sbBulge',mode:'custom',target:'sbBulgeMass',cardSlider:{sbBulgeMass:'curvePlot'},gesture:'increase',watch:['curvePlot:Watch the left of the cyan line'],allow:['sbBulgeRadius'],title:'Add a bulge',
         text:'A bulge is a ball of stars at the centre. It speeds up the innermost orbits.',
         action:'Drag Bulge mass above 0. Watch the left of the cyan line climb.'},
        {id:'sbGas',mode:'custom',target:'sbGasMass',cardSlider:{sbGasMass:'curvePlot'},gesture:'adjust',tag:'Drag both ways',watch:['curvePlot:Watch the right of the cyan line'],title:'Add the gas',
         text:'Gas spreads about three times wider than the stars, so it matters most far out.',
         action:'Drag Gas mass up and down. Watch the right of the cyan line.'},
        {id:'sbSummary',mode:'custom',target:'sandboxSummary',cardSlider:{sbDiskMass:'sandboxSummary'},cardSliderPhone:true,gesture:'look',watch:['sbHaloBaseline:Typical halo mass'],tag:'Watch the totals',also:['sbDiskMass:Drag this'],allow:['sandboxGroups'],title:'Check the totals',
         text:'You never set the dark matter directly: your galaxy gets the halo typical for its stars (Moster et al. 2013).',
         action:'Drag Disk mass and watch Typical halo mass change.'},
        {id:'sbMond',mode:'custom',partPause:true,given:['optMond'],target:'optMond',gesture:'click',tag:'1 · Choose MOND',also:['sbDiskMass:2 · Drag right'],cardSlider:{sbDiskMass:'curvePlot'},watch:['curvePlot:Watch the gold zone'],title:'Where MOND steps in',
         // The default disk is already heavy enough for faint gold, so the part also needs a real drag to the right
         // (or the slider already at its top end).
         checks:{optMond:()=>mondMode,sbDiskMass:()=>mondMode&&mondFaintSomewhere()&&(sbVal('sbDiskMass')>tutorialStart.sbDiskMass||sbVal('sbDiskMass')>=Number(document.getElementById('sbDiskMass').max))},doneText:'Near a heavy disk’s centre gravity is strong, so MOND changes it only a little. Its boost grows further out.',
         text:'MOND boosts gravity everywhere, but only a little where gravity is strong. The weaker it gets, the bigger the boost, and the deeper the gold on the graph.',
         action:'1. Choose MOND. 2. Drag Disk mass to the right until the gold near the centre of the graph turns very faint.'},
        {id:'sbMondZone3d',mode:'custom',mond:true,viewPanel:true,given:['toggleMondZone'],target:'toggleMondZone',gesture:'click',tag:'Switch on',title:'MOND in your galaxy',
         checks:{toggleMondZone:()=>mondMode&&mondZoneShown},doneText:'Dim in the middle, where your disk pulls hard, and deeper gold outwards, where MOND takes over.',
         text:'The graph’s gold, laid on your galaxy: the deeper the gold, the more MOND boosts gravity there.',
         action:'Switch on Where MOND acts in the View panel.'},
        {id:'sbMonitor',mode:'custom',nfw:true,halo:true,target:'stabilityMonitor',gesture:'look',tag:'Look at the gauge',title:'Meet the halo shape monitor',
         text:`This monitor checks your halo has a believable shape, using its concentration c. The shaded band is the realistic range, c from ${STABLE_C_MIN} to ${STABLE_C_MAX}.`,
         action:'Look at the monitor (blue box) and find its marker. Is it inside the shaded band? Then press Next.'},
        {id:'sbStability',mode:'custom',nfw:true,halo:true,target:'scaleSlider',cardSlider:{scaleSlider:'stabilityMonitor'},tag:'Drag to each end',watch:['stabilityMonitor:Watch the marker','stabilityChip:Warning shows here'],gesture:'adjust',title:'Break the halo',
         text:'A tiny or huge scale radius makes a halo unlike any real one. The monitor turns red, but it still runs.',
         action:'Drag Scale radius to one end until the monitor turns red.'},
        {id:'sbHalo',mode:'custom',nfw:true,target:'massSlider',cardSlider:{massSlider:'curvePlot'},tag:'1 · Drag',also:['applyObservedFit:2 · Reset'],watch:['curvePlot:Watch the orange line'],gesture:'decrease',
         checks:{massSlider:()=>Number(massSlider.value)!==tutorialStart.mass,applyObservedFit:()=>tutorialFlags.typical},doneText:'Halo back to typical. Try a heavy or light halo in fullscreen after the tour.',title:'Tune the halo yourself',
         text:'Halo mass runs from 0 to 5 × the typical halo. The stars start balanced for the typical one.',
         action:'1. Drag Halo mass low, then high. 2. Press Typical halo to reset it.'},
        {id:'sbFs',given:['toggleFullscreen'],mode:'custom',nfw:true,viewPanel:true,target:'toggleFullscreen',gesture:'click',tag:'1 · Fullscreen',also:['fsPlay:2 · Play'],freeCamera:true,title:'Watch your galaxy',
         checks:{toggleFullscreen:tutorialIsFullscreen,fsPlay:()=>playing},doneText:'Fullscreen and playing: watch for as long as you like, then press Next.',
         text:'Now watch your own galaxy run.',
         action:'1. Press Fullscreen: the tour comes too. 2. Press ▶ Play (or Space).'},
        {id:'sbOrbits',given:['toggleOrbits','toggleGas','toggleHalo','fsPlay'],fs:true,viewPanel:true,marks:['bulge'],target:'toggleOrbits',tag:'1 · Orbits on',also:['toggleGas:2 · Gas off','toggleHalo:3 · Halo off','fsPlay:4 · Play'],gesture:'click',freeCamera:true,title:'Two kinds of orbit',
         checks:{toggleOrbits:()=>orbitLines.visible,toggleGas:()=>!gasCloud.visible,toggleHalo:()=>haloHidden||Number(massSlider.value)===0,fsPlay:()=>playing||tutorialPlayedMyr>0},
         text:'Disk stars trace circles. Bulge stars (the ringed one) trace tangled loops: random motions hold a bulge up.',
         action:'1. Switch on Orbits. 2. Switch off Gas. 3. Switch off Halo. 4. Press ▶ Play.',
         doneText:'Circles in the disk, a tangle in the bulge.'},
        {id:'sbFsEscape',given:['toggleEjection'],fs:true,nfw:true,viewPanel:true,marks:['escaper'],ring:true,target:'toggleEjection',tag:'1 · Ejection line on',also:['fsMassSlider:2 · Drag to 0','fsRestart:3 · Restart'],
         checks:{toggleEjection:()=>ejectionLineShown,fsMassSlider:()=>Number(massSlider.value)===0,fsRestart:()=>tutorialFlags.restart0},gesture:'click',watch:['escape-counter:Ejected count'],title:'Starve it of dark matter',
         text:'Take the halo away and the ringed star is the first your galaxy loses.',
         action:'1. Switch on Ejection line. 2. Drag Halo mass to 0. 3. Press Restart and watch the ringed star.',
         done:prev=>{const k=firstEscaper();return Number(massSlider.value)===0&&k>=0&&prev<escapeTimeByStar[k]&&simTime>=escapeTimeByStar[k];},doneText:'Lost: without dark matter, the outer stars fly off.'}
      ]}
    };
    const tutorialWelcome=document.getElementById('tutorialWelcome'),tutorialCard=document.getElementById('tutorialCard');
    const tutorialTitle=document.getElementById('tutorialTitle'),tutorialProgress=document.getElementById('tutorialProgress');
    const tutorialText=document.getElementById('tutorialText'),tutorialTry=document.getElementById('tutorialTry');
    const tutorialBack=document.getElementById('tutorialBack'),tutorialNext=document.getElementById('tutorialNext');
    const tutorialMore=document.getElementById('tutorialMore');
    const tutorialLaunch=document.getElementById('tutorialLaunch');
    const tutorialSpotlight=document.getElementById('tutorialSpotlight');
    const tutorialProgressFill=document.getElementById('tutorialProgressFill');
    const tutorialFeedback=document.getElementById('tutorialFeedback');
    const visibleMatterOnly=document.getElementById('visibleMatterOnly');
    const sourceToggle=document.getElementById('sourceToggle');
    const tutorialTopbar=document.querySelector('.topbar'),tutorialBar=document.querySelector('.sim-stage .controls'),tutorialPlots=document.querySelector('.sim-stage .panel-grid');
    let tutorialTour='main',tutorialSteps=TOURS.main.steps,tutorialIndex=0,tutorialTarget=null,tutorialCompleted={};
    // The full guided tour (from the puzzle page) runs every tour in this order, one after the other.
    // Its closing cards (Scorecard, then Background) are held back until after the sandbox tour, so it ends on them.
    TOURS.closing={label:'Wrapping up',steps:TOURS.main.steps.filter(s=>s.closing)};
    const FULL_TOUR=['main','sandbox','closing'];let tutorialFull=false;
    const tourSteps=(tour,full)=>full&&tour==='main'?TOURS.main.steps.filter(s=>!s.closing):TOURS[tour].steps;
    let tutorialLastMass=0,tutorialLastScale=0,tutorialLastTime=0,tutorialDrag=null,tutorialShortcutUnlocked=false,tutorialBlockTarget=null,tutorialAimTarget=null,tutorialCarryTarget=null;
    const tutorialStep=()=>tutorialSteps[tutorialIndex];
    // Multi-part steps: checks maps each ringed control to a test of whether its part is done. A part, once done,
    // stays done for that card; its ring disappears and the pointer moves on to the next part.
    // tutorialStale: when a card opens with every part already true (an earlier card was skipped, so the page is
    // already in this card's end state), a part only counts once it has been undone and redone on this card, so a
    // skipped card never earns praise for nothing. Holds the ids seen false since the card opened (false = not stale).
    let tutorialFlags={},tutorialLatched={},tutorialStart={},tutorialCheckSig='',tutorialStale=null;
    // step.given: parts that count as done when already true as the card opens (MOND already chosen, orbit lines
    // already on, already playing), so a card never asks for a switch to be undone just to be redone. step.partPause: when a part is done, its ring turns green and
    // fades before the next part appears (tutorialPartHold holds the done part's id until then).
    let tutorialPartHold=null,tutorialPartHoldTimer=0;
    function holdTutorialPart(id){
      clearTimeout(tutorialPartHoldTimer);tutorialPartHold=id;
      if(id)tutorialPartHoldTimer=setTimeout(()=>{tutorialPartHold=null;refreshTutorialGuidance();},TUTORIAL_RING_FADE_MS);
    }
    const sbVal=id=>{const el=document.getElementById(id);return el?Number(el.value):0;};
    function tutorialResetChecks(){
      tutorialFlags={};tutorialLatched={};tutorialCheckSig='';tutorialStale=null;holdTutorialPart(null);
      tutorialStart={mass:Number(massSlider.value),scale:Number(scaleSlider.value),sbDiskMass:sbVal('sbDiskMass'),sbDiskScale:sbVal('sbDiskScale'),sbTruncation:sbVal('sbTruncation')};
    }
    function tutorialUpdateChecks(){
      const step=tutorialStep();
      if(!step.checks||tutorialCard.hidden||tutorialBlocker())return;
      if(tutorialStale===null){
        (step.given||[]).forEach(id=>{if(step.checks[id]())tutorialLatched[id]=true;});
        const rest=Object.keys(step.checks).filter(id=>!tutorialLatched[id]);
        tutorialStale=rest.length&&rest.every(id=>step.checks[id]())?{}:false;
      }
      const before=step.partPause&&!tutorialPartHold?tutorialCurrentAim(step):null;
      for(const id in step.checks){
        const ok=step.checks[id]();
        if(tutorialStale&&!ok)tutorialStale[id]=true;
        // A part only counts once the slider (or finger) is let go, so its ring never moves on mid-drag;
        // this check runs every frame, so it latches on the first frame after release.
        if(!tutorialHolding&&!tutorialLatched[id]&&ok&&(!tutorialStale||tutorialStale[id]))tutorialLatched[id]=true;
      }
      const sig=Object.keys(tutorialLatched).join();
      if(sig===tutorialCheckSig)return;
      tutorialCheckSig=sig;
      const allDone=Object.keys(step.checks).every(id=>tutorialLatched[id]);
      if(allDone)holdTutorialPart(null);else if(before&&tutorialLatched[before.id])holdTutorialPart(before.id);
      if(!step.done&&allDone)completeTutorialAction(step.id,step.doneText,!!step.fs||tutorialIsFullscreen());
      renderTutorialTry();
      updateTutorialSpotlight();
    }
    // The control to point at now: the first part not yet done (null once every part is done).
    function tutorialCurrentAim(step){
      if(tutorialCompleted[step.id]&&step.checks||tutorialPartHold)return null;
      // A predict-first card points at its answer buttons (on the card) before anything else.
      if(step.predict&&tutorialFlags.prediction==null)return{id:'tutorialChoices',label:'1 · Predict',gesture:'click'};
      const first={id:step.target,label:step.tag,gesture:step.gesture};
      if(!step.checks)return first;
      const parts=[first,...(step.also||[]).map(e=>{const k=e.indexOf(':'),id=k<0?e:e.slice(0,k),el=document.getElementById(id);
        // A slider part labelled 'Drag to 0' points left to 0, 'Drag right' points right; any other slider part either way.
        const label=k<0?'':e.slice(k+1),range=el&&el.matches('input[type=range]');
        return{id,label,gesture:range?(/\bto 0\b/.test(label)?'zero':/\bright\b/.test(label)?'increase':'adjust'):'click'};})];
      return parts.find(p=>!tutorialLatched[p.id])||null;
    }
    // Slider steps whose result shows far from the slider (the sandbox sliders sit above the plots): the card carries a
    // copy of every slider the card uses, and the page scrolls to what the current part's slider changes.
    // step.cardSlider maps a slider id to the id of the view to watch. tutorialSliderId: the current part's slider.
    const tutorialSliders=document.getElementById('tutorialSliders');
    let tutorialSliderId=null,tutorialSliderSig='';
    // The card's copy of a real slider (null when the card does not carry one).
    const tutorialCopy=id=>id&&!tutorialSliders.hidden?tutorialSliders.querySelector(`input[data-real="${id}"]`):null;
    // cardSliderPhone: the copy is only needed on a phone, where the slider and what it changes are far apart.
    const tutorialCarriesSliders=step=>!tutorialCard.hidden&&!!step.cardSlider&&!tutorialBlocker()&&!(step.cardSliderPhone&&!tutorialPhone());
    function tutorialCardSliderId(){
      const step=tutorialStep();
      if(!tutorialCarriesSliders(step))return null;
      // A card with no parts whose ring is on something else (the totals) still carries its slider.
      const aim=tutorialCurrentAim(step),keys=Object.keys(step.cardSlider);
      const id=tutorialPartHold||(aim?(step.cardSlider[aim.id]||step.checks?aim.id:keys[0]):keys.pop()),real=step.cardSlider[id]&&document.getElementById(id);
      if(!real)return null;
      return id;
    }
    // scroll: glide to the new view when the part (and so the slider) changes.
    function syncTutorialSlider(scroll){
      const step=tutorialStep(),ids=tutorialCarriesSliders(step)?Object.keys(step.cardSlider).filter(id=>document.getElementById(id)):[];
      // The copies are rebuilt only when the card's set of sliders changes, and never mid-drag.
      if(ids.join()!==tutorialSliderSig&&!(tutorialHolding&&tutorialSliderSig)){
        tutorialSliderSig=ids.join();tutorialSliders.hidden=!ids.length;
        tutorialSliders.replaceChildren(...ids.map(id=>{
          const real=document.getElementById(id),box=document.createElement('div');
          box.className='tutorial-slider';
          box.innerHTML=`<label for="tutorialCopy-${id}"><span>${real.labels[0].querySelector('span').innerHTML}</span><output for="tutorialCopy-${id}"></output></label><input id="tutorialCopy-${id}" type="range" autocomplete="off" data-real="${id}">`;
          const input=box.querySelector('input');['max','min','step'].forEach(k=>{input[k]=real[k];});
          return box;
        }));
      }
      const id=tutorialCardSliderId();
      if(id!==tutorialSliderId&&!(tutorialHolding&&tutorialSliderId)){
        tutorialSliderId=id; // never swap mid-drag
        if(scroll)scrollTutorialTarget(true);
      }
      tutorialSliders.querySelectorAll('input').forEach(copy=>{
        const real=document.getElementById(copy.dataset.real),shown=real.labels[0].querySelector('output').textContent,out=copy.labels[0].querySelector('output');
        if(copy.value!==real.value)copy.value=real.value;
        if(out.textContent!==shown)out.textContent=shown;
      });
    }
    // A copy drives its real slider, so every listener (the model, the tour's checks) runs as if it moved itself.
    const tutorialRealOf=copy=>copy&&copy.dataset&&copy.dataset.real?document.getElementById(copy.dataset.real):null;
    tutorialSliders.addEventListener('input',e=>{
      const real=tutorialRealOf(e.target);if(!real)return;
      real.value=e.target.value;real.dispatchEvent(new Event('input',{bubbles:true}));
      syncTutorialSlider();
    });
    tutorialSliders.addEventListener('change',e=>{const real=tutorialRealOf(e.target);if(real)real.dispatchEvent(new Event('change',{bubbles:true}));});
    tutorialSliders.addEventListener('pointerdown',()=>tutorialSpotlight.classList.add('is-interacting'));
    // Watch steps finish by themselves: each render checks the step's done() test (camera angle, layers, play time).
    let tutorialPlayedMyr=0,tutorialWatchLast=0,tutorialScaleSeen=[0,0,0];
    function tutorialWatchTick(){
      if(tutorialCard.hidden)return;
      const prev=tutorialWatchLast;
      if(simTime>prev)tutorialPlayedMyr+=simTime-prev;
      else if(simTime<prev)raceStartMyr=simTime; // the clock went back (Restart, or a drag): the race's lap counts start again from here
      tutorialWatchLast=simTime;
      const step=tutorialStep();
      // Watch steps never count down: the viewer keeps watching and presses Next when ready.
      if(step.done&&!tutorialCompleted[step.id]&&!tutorialBlocker()&&step.done(prev))completeTutorialAction(step.id,typeof step.doneText==='function'?step.doneText():step.doneText,true,!!step.predict&&tutorialFlags.prediction!==step.predict.answer);
    }
    // What has to change on the page before this step's own action can be done (null when nothing is in the way).
    function tutorialBlocker(){
      const blocker=tutorialPageBlocker();
      if(!blocker)return null;
      // Fullscreen hides the page, so anything to fix there means leaving fullscreen first (Esc, or the button in View).
      if(tutorialIsFullscreen()&&!orbitPlot.contains(blocker.target)){
        const viewClosed=viewPanel.classList.contains('collapsed');
        return{target:viewClosed?viewToggle:fsBtn,gesture:'click',say:`Leave fullscreen first: press Esc${viewClosed?' (or open View and press Fullscreen)':' or Fullscreen'}. Then: ${blocker.say}`};
      }
      return blocker;
    }
    function tutorialPageBlocker(){
      const step=tutorialStep(),fullscreen=tutorialIsFullscreen();
      if(!step.fs&&fullscreen&&!orbitPlot.contains(document.getElementById(step.target))&&!/^(fsEnter|sbFs)$/.test(step.id)){
        const viewClosed=viewPanel.classList.contains('collapsed');
        return{target:viewClosed?viewToggle:fsBtn,gesture:'click',say:'Leave fullscreen first: press Esc. '};
      }
      const stepPage=pageOf(document.getElementById(step.target));
      if(stepPage&&stepPage!==activePage)return{target:pageLink(stepPage),gesture:'click',say:`Go back to the ${pageLink(stepPage).textContent} page first (top menu). `};
      if(step.mode&&step.mode!==galaxyMode){
        return{target:sourceToggle.querySelector(`input[value=${step.mode}]`).parentElement,gesture:'click',
          say:step.mode==='custom'?'Switch to Custom Sandbox first. ':'Switch back to Preset Galaxies first. '};
      }
      if(step.nfw&&mondMode){
        return{target:document.getElementById('optNfw'),gesture:'click',
          say:step.halo?'Choose Dark matter halo first; MOND does not use one. ':'Choose Dark matter halo first. '};
      }
      if(step.mond&&!mondMode)return{target:document.getElementById('optMond'),gesture:'click',say:'Choose MOND first. '};
      // Page settings come first: fixing them needs the page, so fullscreen is asked for last.
      if(step.fs&&!fullscreen){
        const viewClosed=viewPanel.classList.contains('collapsed');
        return{target:viewClosed?viewToggle:fsBtn,gesture:'click',say:viewClosed?'Open View in the 3D panel and press Fullscreen first. ':'Press Fullscreen first. '};
      }
      if(step.halo&&Number(massSlider.value)===0)return{target:fullscreen?fsMass:massSlider,gesture:'increase',say:'Increase halo mass above 0 first so there is a halo to see. '};
      if(step.viewPanel&&viewPanel.classList.contains('collapsed'))return{target:viewToggle,gesture:'click',say:'Open View first. '};
      return null;
    }
    // The fitting goal: a Close fit, or for a galaxy whose published halo fit can't reach one, within
    // TUTORIAL_FIT_SLACK of that fit's reduced χ², so the goal is always reachable with the two sliders.
    const TUTORIAL_FIT_SLACK=1.25;
    let tutorialGoalCache=null;
    function tutorialChiGoal(){
      if(tutorialGoalCache&&tutorialGoalCache.galaxy===activeGalaxy)return tutorialGoalCache.goal;
      const score=scoreGalaxy(activeGalaxy),fitRed=reducedChi(score.chi.nfwStd,score.pointCount,score.k.nfwStd);
      tutorialGoalCache={galaxy:activeGalaxy,goal:Math.max(CLOSE_FIT_RED_CHI,fitRed*TUTORIAL_FIT_SLACK)};
      return tutorialGoalCache.goal;
    }
    // Add the halo back is done once the model reaches TUTORIAL_REACH_FRACTION of the measured speeds,
    // averaged over the outer third of the points (where the gap is).
    const TUTORIAL_REACH_FRACTION=0.9;
    function tutorialOuterReached(){
      if(!shownModelKms)return false;
      let model=0,measured=0;
      for(let i=Math.floor(obsV.length*2/3);i<obsV.length;i++){model+=shownModelKms[i];measured+=obsV[i];}
      return model>=TUTORIAL_REACH_FRACTION*measured;
    }
    function tutorialInstruction(){
      const step=tutorialStep();
      if(step.id==='puzzle'&&!tutorialCompleted.puzzle){const n=[0,1,2].filter(i=>tutorialFlags['scan'+i]).length;if(n)return step.action+` (${n} of 3 stretches covered)`;}
      if(step.id==='fit'&&tutorialChiGoal()>CLOSE_FIT_RED_CHI){
        const goal=tutorialChiGoal().toFixed(1);
        return step.action.replace('the yellow box no longer reads Close fit',`χ²ᵣ (the second number in the yellow box) rises above ${goal}`)
          .replace('the yellow box reads Close fit',`χ²ᵣ drops below ${goal}, close to the best the published fit manages for this galaxy`);
      }
      if(step.id==='removeHalo'){
        if(tutorialShortcutUnlocked)return 'Compare the orange line with the white points on the right. (The new Visible matter only button is a shortcut back to 0.)';
        if(Number(massSlider.value)===0)return 'Halo mass is already 0: drag it right, then back to 0.';
      }
      return step.action;
    }
    // Which numbered part of the action the viewer is on: the first part not yet done (the last part once all are).
    function tutorialPartIndex(step,count){
      if(tutorialCompleted[step.id])return count-1;
      if(step.checks){const i=tutorialPartHold?Object.keys(step.checks).indexOf(tutorialPartHold):Object.keys(step.checks).findIndex(id=>!tutorialLatched[id]);return i<0?count-1:Math.min(i,count-1);}
      if(step.id==='haloShape'){const[start,min]=tutorialScaleSeen;return min<start||start<=Number(scaleSlider.min)?1:0;}
      return 0;
    }
    // "First:" line for anything in the way, then the action; "1. … 2. …" actions show one part at a time,
    // keeping its number, and the next part appears once that one is done.
    let tutorialPartKey='';
    // Splits a "1. … 2. … 3. …" action at its own numbers, taken in order, so a value in the text
    // ("Drag Halo mass to 0. 3. Press Restart") is never mistaken for the start of a part.
    function tutorialParts(action){
      if(!/^1\.\s/.test(action))return[action];
      const parts=[];let rest=action;
      for(let n=2;;n++){const k=rest.indexOf(` ${n}. `);if(k<0)break;parts.push(rest.slice(0,k));rest=rest.slice(k+1);}
      parts.push(rest);
      return parts;
    }
    function renderTutorialTry(){
      const blocker=tutorialBlocker(),action=tutorialInstruction();
      tutorialTry.replaceChildren();
      if(blocker){const first=document.createElement('strong');first.className='tutorial-first';first.textContent='First: '+blocker.say.replace(/\s+first\b/gi,'').trim();tutorialTry.appendChild(first);}
      const step=tutorialStep();
      if(step.predict&&!blocker){
        const pick=tutorialFlags.prediction;
        if(pick==null){tutorialTry.append(tutorialPredictionPrompt(step.predict));return;}
        const note=document.createElement('span');note.className='tutorial-prediction';note.textContent='Your prediction: '+step.predict.options[pick];tutorialTry.appendChild(note);
      }
      const head=document.createElement('strong');head.textContent=blocker?'Then: ':'Instructions: ';tutorialTry.appendChild(head);
      const parts=tutorialParts(action);
      if(parts.length>1){
        const i=tutorialPartIndex(tutorialStep(),parts.length),list=document.createElement('ol'),li=document.createElement('li');
        // A new part fades in, so the instruction never just jumps.
        const partKey=tutorialStep().id+':'+i;
        if(tutorialPartKey&&partKey!==tutorialPartKey&&tutorialPartKey.startsWith(tutorialStep().id+':')){tutorialTry.classList.remove('part-in');void tutorialTry.offsetWidth;tutorialTry.classList.add('part-in');}
        tutorialPartKey=partKey;
        list.start=i+1;li.textContent=parts[i].replace(/^\d\.\s/,'');list.appendChild(li);
        if(i<parts.length-1){const more=document.createElement('span');more.className='tutorial-more-parts';more.textContent=`Part ${i+1} of ${parts.length}: the next part appears when this one is done.`;tutorialTry.appendChild(list);tutorialTry.appendChild(more);return;}
        tutorialTry.appendChild(list);
      }else tutorialTry.appendChild(document.createTextNode(action));
    }
    // A predict-first question: the step's controls stay locked (and the race paused) until it is answered.
    function tutorialPredictionPrompt(predict){
      const box=document.createElement('div'),q=document.createElement('strong'),row=document.createElement('div');
      q.textContent=predict.question;row.className='tutorial-choices';row.id='tutorialChoices';
      predict.options.forEach((label,i)=>{const b=document.createElement('button');b.type='button';b.textContent=label;
        b.addEventListener('click',()=>{tutorialFlags.prediction=i;tutorialPlayedMyr=0;tutorialWatchLast=raceStartMyr=simTime;refreshTutorialGuidance();});row.appendChild(b);});
      box.append(q,row);return box;
    }
    function refreshTutorialGuidance(){
      if(tutorialCard.hidden)return;
      renderTutorialTry();
      syncTutorialSlider(true);
      updateTutorialSpotlight();
    }
    // A gentle eased glide (our own animation, so it is smooth in every browser), about 0.6 to 1.1 s depending on
    // distance. A wheel, touch or key press from the viewer stops it at once.
    let tutorialGlide=0;
    function tutorialGlideTo(y){
      cancelAnimationFrame(tutorialGlide);
      const from=window.scrollY,dist=y-from;if(Math.abs(dist)<2)return;
      const dur=Math.min(1100,600+Math.abs(dist)*.25),t0=performance.now();
      const ease=t=>t<.5?4*t*t*t:1-Math.pow(-2*t+2,3)/2;
      const stop=()=>{cancelAnimationFrame(tutorialGlide);['wheel','touchstart','keydown'].forEach(t=>window.removeEventListener(t,stop));};
      ['wheel','touchstart','keydown'].forEach(t=>window.addEventListener(t,stop,{passive:true}));
      (function step(now){
        const t=Math.min(1,(now-t0)/dur);
        window.scrollTo({top:from+dist*ease(t),behavior:'instant'});
        if(t<1)tutorialGlide=requestAnimationFrame(step);else stop();
      })(t0);
    }
    function scrollTutorialTarget(smooth){
      const step=tutorialStep(),blocker=tutorialBlocker();
      const view=!blocker&&tutorialSliderId&&document.getElementById(step.cardSlider[tutorialSliderId]);
      const target=blocker?blocker.target:view||(step.scroll&&innerWidth>960?document.getElementById(step.scroll):tutorialTarget);
      if(tutorialTopbar.contains(target))return; // the top menu is always in view
      const y=window.scrollY,clear=tutorialTopbar.offsetHeight+16; // the gap kept between the top menu and the control
      // A slider comes into view with its label and value (its whole .control), not just its track.
      const boxTop=el=>(el.closest('.control')||el).getBoundingClientRect().top;
      // The other things the card rings (to use or to watch), in the order the step lists them.
      const extras=blocker?[]:tutorialRingList(step).map(r=>document.getElementById(r.id)).filter(el=>el&&el.getClientRects().length);
      // Anything ringed just above the control (the M/L and tilt tiles beside the best-fit switch) comes with it.
      const targetTop=boxTop(target);
      let top=Math.min(targetTop,...extras.map(boxTop).filter(t=>t<targetTop&&t>targetTop-TUTORIAL_NEAR_PX))+y-clear;
      // Sliders and plots: line the slider bar up under the top bar so as much as possible of the plots it changes is in view.
      const anchor=document.getElementById('controlsAnchor'),inPlots=tutorialPlots.contains(target);
      if(tutorialBar&&innerWidth>960&&!document.fullscreenElement&&(tutorialBar.contains(target)||inPlots)){
        top=anchor.getBoundingClientRect().top+y-56;
        // A target in the plots (e.g. the halo numbers): scroll on until it is in view, but no further than keeps the bar
        // on screen. scrollPad leaves room below for boxes that open on hover (the halo cards' explainers).
        const low=inPlots?target:null;
        if(low){const r=low.getBoundingClientRect(),need=r.bottom+y-innerHeight+16+(step.scrollPad||0);if(need>top)top=Math.min(need,r.top+y-64-tutorialBar.offsetHeight);}
      }
      // Phones: when the control (or what its card slider changes) and everything else the card rings fit between the
      // top menu and the tour card together, scroll so they are all in view at once.
      if(tutorialPhone()&&!blocker){
        const els=[view||target,...extras.filter(el=>!(step.cardSlider&&step.cardSlider[el.id]))],box=el=>(el.closest('.control')||el).getBoundingClientRect();
        const spanTop=Math.min(...els.map(el=>box(el).top))+y,spanBottom=Math.max(...els.map(el=>box(el).bottom))+y;
        const room=innerHeight-tutorialCard.offsetHeight-20-clear;
        if(spanBottom-spanTop<=room)top=spanTop-clear;
      }
      if(smooth)tutorialGlideTo(Math.max(0,top));else window.scrollTo({top:Math.max(0,top),behavior:'auto'});
    }
    // Once whatever was in the way has been sorted out (or something new is), bring the right control into view.
    function tutorialRecheck(){
      if(tutorialCard.hidden)return;
      const blocker=tutorialBlocker(),target=blocker?blocker.target:null;
      if(target!==tutorialBlockTarget){
        tutorialBlockTarget=target;
        if(!tutorialIsFullscreen())scrollTutorialTarget();
      }
      refreshTutorialGuidance();
    }
    // A new galaxy resets the halo sliders and the timeline, so resync what the tutorial last saw.
    document.addEventListener('galaxychange',()=>{
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialLastTime=Number(timeSlider.value);
      refreshTutorialGuidance();
    });
    // Opening or closing the sandbox panel moves the page content, so re-aim the tutorial pointer.
    document.addEventListener('galaxymodechange',()=>{
      tutorialLaunch.textContent=TOURS[galaxyMode==='custom'?'sandbox':'main'].label;
      tutorialRecheck();
    });
    // noAuto: no countdown, so the viewer can watch the simulation and press Next when ready.
    // A word of praise before each 'done' message, taken in turn so it varies from card to card.
    const TUTORIAL_PRAISE=['Nice work!','Great!','Spot on!','Well done!','Brilliant!','Perfect!'];
    // A wrong prediction is the one place the tour gives negative feedback: orange 'Not quite.', no tick, no green.
    const tutorialMissed={};
    function showTutorialFeedback(text,missed){
      tutorialFeedbackText.textContent='';
      tutorialFeedback.classList.toggle('is-missed',!!missed);
      if(missed){tutorialFeedbackText.append(Object.assign(document.createElement('b'),{textContent:'Not quite.'}),' '+text);return;}
      const cut=text.indexOf('!')+1;
      tutorialFeedbackText.append(Object.assign(document.createElement('b'),{textContent:'✓ '+text.slice(0,cut)}),text.slice(cut));
    }
    // Done while a slider is still held: wait until it is let go, so the viewer can keep dragging and watching.
    let tutorialPending=null;
    // missed: a wrong prediction (no praise, no tick, orange instead of green).
    function completeTutorialAction(id,message,noAuto,missed){
      if(tutorialCard.hidden||tutorialStep().id!==id||tutorialCompleted[id])return;
      if(tutorialHolding){tutorialPending=[id,message,noAuto,missed];return;}
      tutorialPending=null;
      tutorialMissed[id]=!!missed;
      tutorialCompleted[id]=missed?message:TUTORIAL_PRAISE[tutorialIndex%TUTORIAL_PRAISE.length]+' '+message;
      tutorialFeedback.hidden=false;
      showTutorialFeedback(tutorialCompleted[id],missed);
      tutorialCard.classList.remove('is-celebrating');void tutorialCard.offsetWidth;if(!missed)tutorialCard.classList.add('is-celebrating');
      if(!noAuto)startTutorialAuto();
      tutorialSpotlight.classList.toggle('is-missed',!!missed);
      tutorialSpotlight.classList.remove('is-revisit');tutorialSpotlight.classList.add('is-complete');
      updateTutorialNextGate();
    }
    // The 'press Fullscreen' steps cannot be skipped with Next: the steps after them need fullscreen.
    function updateTutorialNextGate(){
      const step=tutorialStep(),gated=/^(fsEnter|sbFs)$/.test(step.id)&&!tutorialIsFullscreen();
      tutorialNext.disabled=gated;
      tutorialNext.title=gated?'Press Fullscreen to carry on':'';
    }
    const tutorialPhone=()=>matchMedia('(max-width:600px)').matches;
    function tutorialIsFullscreen(){return !!nativeFsElement()||orbitPlot.classList.contains('fullscreen-active');}
    // Keeps the highlight glued to its control when the control or the page changes size without a scroll or window
    // resize: a hover box opening inside it, a <details> section expanding above it, the sandbox panel appearing.
    // ResizeObserver runs after layout and before paint, so a ring is re-fitted in the same frame its control grows
    // (the χ² badge widening as its numbers get longer), never a frame or more later.
    const spotlightResize=window.ResizeObserver?new ResizeObserver(()=>{if(!tutorialCard.hidden){updateTutorialSpotlight();positionTutorialAlso();}}):null;
    if(spotlightResize)spotlightResize.observe(document.body);
    // When a part of a card is done and the main ring moves on, a green copy stays where it was and fades slowly.
    // Each ghost is kept with its control, so it follows the control if the page shifts while it fades.
    let tutorialGhosts=[],tutorialAimId=null; // tutorialAimId: the part the main ring is on
    function leaveTutorialGhost(){
      const prev=tutorialAimTarget;
      if(!prev||!tutorialLatched[tutorialAimId]||tutorialSpotlight.hidden||calmMotion())return;
      const ghost=document.createElement('div');
      ghost.className='tutorial-also is-leaving'+(tutorialSpotlight.classList.contains('is-fixed')?' is-fixed':'');ghost.setAttribute('aria-hidden','true');
      ['left','top','width','height','borderRadius','clipPath'].forEach(k=>{ghost.style[k]=tutorialSpotlight.style[k];});
      tutorialSpotlight.parentElement.appendChild(ghost);
      tutorialGhosts.push([ghost,prev]);
      setTimeout(()=>ghost.remove(),TUTORIAL_RING_FADE_MS);
    }
    function updateTutorialSpotlight(){
      if(tutorialCard.hidden){tutorialSpotlight.hidden=true;positionTutorialAlso(true);return;}
      positionTutorialAlso();
      const step=tutorialStep(),blocker=tutorialBlocker(),aim=blocker?null:tutorialCurrentAim(step);
      // noRing: a card about a whole section (the last one) gets no box or dimming, which would cover most of the page.
      if(step.noRing&&!blocker){tutorialSpotlight.hidden=true;return;}
      // A many-part card just finished: the ring stays on the last part so it can turn green and fade out slowly.
      const finished=!blocker&&!aim&&tutorialCompleted[step.id]&&tutorialAimTarget&&tutorialAimTarget.isConnected&&!tutorialSpotlight.classList.contains('is-revisit');
      // Every ringed part is done but the card still waits on something unringed (Score your fit's χ² goal): the last
      // part still turns green and fades like the others, then the main ring goes.
      if(!blocker&&!aim&&!finished){
        if(tutorialAimTarget){leaveTutorialGhost();tutorialAimTarget=null;tutorialAimId=null;}
        tutorialSpotlight.hidden=true;return;
      }
      let target=finished?tutorialAimTarget:blocker?blocker.target:tutorialCopy(aim.id)||document.getElementById(aim.id),gesture=finished?'click':blocker?blocker.gesture:aim.gesture;
      if(gesture==='switch'&&!finished){
        target=document.getElementById('optMond');gesture='click';
      }
      // A new control to aim at brings the pointer back; it fades again once the user's mouse reaches it.
      if(target!==tutorialAimTarget||!finished&&!blocker&&aim.id!==tutorialAimId){
        leaveTutorialGhost();
        if(spotlightResize){if(tutorialAimTarget)spotlightResize.unobserve(tutorialAimTarget);spotlightResize.observe(target);}
        const carried=target===tutorialCarryTarget;
        tutorialAimTarget=target;tutorialAimId=!blocker&&aim?aim.id:null;
        if(!carried){tutorialSpotlight.classList.remove('is-reached','is-arriving');void tutorialSpotlight.offsetWidth;tutorialSpotlight.classList.add('is-arriving');}
      }
      tutorialCarryTarget=null;
      tutorialSpotlight.classList.toggle('in-card',tutorialCard.contains(target));
      dodgeTutorialBox(target);
      const rect=tutorialSpotRect(target),pad=4;
      tutorialSpotLast=rectKey(rect);
      if(!fitTutorialRing(tutorialSpotlight,target)){tutorialSpotlight.hidden=true;return;}
      tutorialSpotlight.hidden=false;
      tutorialSpotlight.classList.toggle('is-click',gesture==='click');
      tutorialSpotlight.classList.toggle('is-look',gesture==='look');
      const tag=document.getElementById('tutorialTag'),label=blocker||!aim?'':tutorialTagNumber(aim.label);
      if(tag.textContent!==label)tag.textContent=label;
      let fromX=rect.width*.5,toX=fromX,fromY=rect.height*.5,toY=fromY;
      if(gesture==='scan'){fromX=rect.width*.55;toX=rect.width*.82;fromY=toY=rect.height*.43;}
      if(gesture==='drag'){fromX=rect.width*.35;toX=rect.width*.62;fromY=rect.height*.5;toY=rect.height*.42;}
      // zero: slide to the left end; decrease / increase: slide left / right; adjust: slide whichever way has more room.
      // The ring also holds the slider's label row, so the pointer is placed on the slider track itself.
      if(gesture==='zero'||gesture==='decrease'||gesture==='increase'||gesture==='adjust'){
        const fraction=(Number(target.value)-Number(target.min))/(Number(target.max)-Number(target.min));
        const track=target.getBoundingClientRect(),x0=track.left-rect.left;
        fromX=x0+8+(track.width-16)*fraction;
        const destination=gesture==='zero'?0:gesture==='decrease'?Math.max(0,Math.min(.45,fraction-.2)):gesture==='increase'?Math.min(1,Math.max(.55,fraction+.2)):(fraction>.5?.2:.8);
        toX=x0+8+(track.width-16)*destination;
        fromY=toY=track.top-rect.top+track.height/2;
      }
      tutorialSpotlight.style.setProperty('--from-x',(fromX+pad-2)+'px');
      tutorialSpotlight.style.setProperty('--to-x',(toX+pad-2)+'px');
      tutorialSpotlight.style.setProperty('--from-y',(fromY+pad-2)+'px');
      tutorialSpotlight.style.setProperty('--to-y',(toY+pad-2)+'px');
    }
    // An info box the highlighted control opens (the M/L ? box) must never sit under the card: if it would, the card
    // moves to the other side of the screen and stays there for the rest of that card.
    let tutorialDodged=false; // moved aside once on this card already (never flips back and forth)
    function dodgeTutorialBox(target){
      if(tutorialDodged||tutorialIsFullscreen()||innerWidth<=600||tutorialCard.contains(target))return;
      const card=tutorialCard.getBoundingClientRect();
      const hit=[...target.querySelectorAll('[aria-expanded="true"][aria-controls]'),...(target.matches('[aria-expanded="true"][aria-controls]')?[target]:[])]
        .map(btn=>document.getElementById(btn.getAttribute('aria-controls'))).filter(pop=>pop&&!pop.hidden)
        .some(pop=>{const p=pop.getBoundingClientRect();return p.width&&p.left<card.right&&p.right>card.left&&p.top<card.bottom&&p.bottom>card.top;});
      if(hit){tutorialDodged=true;tutorialCard.classList.toggle('on-left');}
    }
    // The control plus anything it has opened (an info box from a ? button, a halo card's explainer), so the
    // highlight grows to cover a box that expands out of it.
    // A slider's <output> stretches to fill its label row, so measure the number text itself.
    // A slider is ringed together with its label row (its name and value) and its description, so the box shows what is being changed;
    // the tour card's own slider is ringed round its whole panel.
    const tutorialSliderLabel=el=>el.matches('input[type=range]')&&el.id?document.querySelector(`label[for="${el.id}"]`):null;
    const tutorialRingBox=el=>el&&el.closest?el.closest('.tutorial-slider')||el:el;
    function tutorialBaseRect(target){
      if(tutorialRingBox(target)!==target)return tutorialRingBox(target).getBoundingClientRect();
      let r=target.getBoundingClientRect();
      if(target.tagName==='OUTPUT'&&target.textContent.trim()){const range=document.createRange();range.selectNodeContents(target);r=range.getBoundingClientRect();}
      const label=tutorialSliderLabel(target);
      if(!label)return r;
      // The slider's own description lines (.control-help) in its .control join the box too. They are measured every
      // frame like the rest, so the box grows or shrinks with any change to their wording.
      const control=target.closest('.control');
      const parts=[label,...(control?[...control.children].filter(el=>el.matches('.control-help')):[])].filter(el=>el.getClientRects().length).map(el=>el.getBoundingClientRect());
      if(!parts.length)return r;
      const left=Math.min(r.left,...parts.map(p=>p.left)),top=Math.min(r.top,...parts.map(p=>p.top)),right=Math.max(r.right,...parts.map(p=>p.right)),bottom=Math.max(r.bottom,...parts.map(p=>p.bottom));
      return{left,top,right,bottom,width:right-left,height:bottom-top};
    }
    function tutorialSpotRect(target){
      const r=tutorialBaseRect(target);
      let left=r.left,top=r.top,right=r.right,bottom=r.bottom;
      const openers=[...target.querySelectorAll('[aria-expanded="true"][aria-controls]')];
      if(target.matches('[aria-expanded="true"][aria-controls]'))openers.push(target);
      const grow=pop=>{
        const p=pop.getBoundingClientRect();
        if(!p.width||!p.height)return;
        left=Math.min(left,p.left);top=Math.min(top,p.top);right=Math.max(right,p.right);bottom=Math.max(bottom,p.bottom);
      };
      // Only a box that opens right beside the control joins its ring. A box that opens further away (the M/L ? box,
      // which drops below both stat tiles) would stretch the ring over controls it has nothing to do with.
      const beside=pop=>{const p=pop.getBoundingClientRect();return p.left<=r.right+TUTORIAL_BESIDE_PX&&p.right>=r.left-TUTORIAL_BESIDE_PX&&p.top<=r.bottom+TUTORIAL_BESIDE_PX&&p.bottom>=r.top-TUTORIAL_BESIDE_PX;};
      openers.forEach(btn=>{const pop=document.getElementById(btn.getAttribute('aria-controls'));if(pop&&!pop.hidden&&beside(pop))grow(pop);});
      // An open drop-down list inside the control (the Galaxy menu) is part of it too.
      target.querySelectorAll('[role="listbox"]:not([hidden])').forEach(grow);
      return{left,top,right,bottom,width:right-left,height:bottom-top};
    }
    // A ring whose control has scrolled under the top menu is cut off at that edge, so no highlight is ever drawn over
    // the menu. Returns how far down the page is covered.
    function tutorialCoverBottom(el){
      if(tutorialIsFullscreen()||tutorialTopbar.contains(el))return 0;
      return tutorialTopbar.getBoundingClientRect().bottom;
    }
    // Parts of the page that stay put on screen while it scrolls; a ring on them is fixed too (see .is-fixed).
    const TUTORIAL_SCREEN_FIXED='.topbar,.tutorial-card,.tutorial-done,.fullscreen-active';
    // Fits a ring round its control (with any box the control has opened); false when the control is off screen.
    function fitTutorialRing(ring,el){
      const r=el&&el.getClientRects().length?tutorialSpotRect(el):null;
      if(!r||!r.width||r.bottom<0||r.top>innerHeight||r.right<0||r.left>innerWidth)return false;
      // The ring is border-box, so it sits exactly pad px outside the control on every side, with matching corners.
      const pad=4,radius=parseFloat(getComputedStyle(tutorialRingBox(el)).borderTopLeftRadius)||0;
      const fixed=tutorialIsFullscreen()||!!el.closest(TUTORIAL_SCREEN_FIXED),dx=fixed?0:scrollX,dy=fixed?0:scrollY;
      ring.classList.toggle('is-fixed',fixed);
      ring.style.left=(r.left-pad+dx)+'px';ring.style.top=(r.top-pad+dy)+'px';
      ring.style.width=(r.width+2*pad)+'px';ring.style.height=(r.height+2*pad)+'px';
      // Corners follow the control itself, not a box it has opened (a round ? button with its square info box).
      ring.style.borderRadius=Math.max(TUTORIAL_RING_MIN_RADIUS_PX,Math.min(radius,el.getBoundingClientRect().height/2)+pad)+'px';
      clipTutorialRing(ring,el,{left:r.left-pad,top:r.top-pad,right:r.right+pad,bottom:r.bottom+pad});
      return true;
    }
    /* Pop-ups the visitor opens (info boxes, the galaxy list, the View and slider panels) always sit above the tour's
       rings: a ring is cut away wherever one covers it, and wherever the top menu does. A pop-up that is part of the
       ring's own control (inside it, or taken into its box) is left alone. Settings is lifted above it all instead. */
    const TUTORIAL_POPUPS='#mlExplainer,#tiltExplainer,#starsExplainer,#haloRadiusHelp,#haloScaleHelp,#haloConcentrationHelp,#galaxyMenuList,#viewPanel,#fsSliders';
    const TUTORIAL_RING_GLOW_PX=40; // room kept round a ring for its glow
    const tutorialPopupShown=el=>{if(!el.getClientRects().length)return false;const cs=getComputedStyle(el);return cs.visibility!=='hidden'&&parseFloat(cs.opacity)>0.05;};
    function tutorialPopupRects(el,box){
      return[...document.querySelectorAll(TUTORIAL_POPUPS)].filter(pop=>!pop.contains(el)&&!el.contains(pop)&&tutorialPopupShown(pop))
        .map(pop=>pop.getBoundingClientRect()).filter(p=>!(p.left>=box.left-1&&p.right<=box.right+1&&p.top>=box.top-1&&p.bottom<=box.bottom+1));
    }
    function clipTutorialRing(ring,el,box){
      const g=TUTORIAL_RING_GLOW_PX,L=box.left-g,T=box.top-g,R=box.right+g,B=box.bottom+g,cover=tutorialCoverBottom(el);
      const cuts=tutorialPopupRects(el,box).concat(cover>T?[{left:L,top:T,right:R,bottom:cover}]:[])
        .map(p=>({l:Math.max(L,p.left),t:Math.max(T,p.top),r:Math.min(R,p.right),b:Math.min(B,p.bottom)})).filter(p=>p.l<p.r&&p.t<p.b);
      let clip='';
      if(cuts.length){
        // The covered area as non-overlapping grid cells, so the even-odd rule cuts each place once even where pop-ups overlap.
        const xs=[...new Set(cuts.flatMap(c=>[c.l,c.r]))].sort((a,b)=>a-b),ys=[...new Set(cuts.flatMap(c=>[c.t,c.b]))].sort((a,b)=>a-b);
        const rect=(l,t,r,b)=>`M${(l-box.left).toFixed(1)} ${(t-box.top).toFixed(1)}H${(r-box.left).toFixed(1)}V${(b-box.top).toFixed(1)}H${(l-box.left).toFixed(1)}Z`;
        let d=rect(L,T,R,B);
        for(let j=0;j<ys.length-1;j++)for(let i=0;i<xs.length-1;i++){
          const cx=(xs[i]+xs[i+1])/2,cy=(ys[j]+ys[j+1])/2;
          if(cuts.some(c=>cx>c.l&&cx<c.r&&cy>c.t&&cy<c.b))d+=rect(xs[i],ys[j],xs[i+1],ys[j+1]);
        }
        clip=`path(evenodd,'${d}')`;
      }
      if(ring.style.clipPath!==clip)ring.style.clipPath=clip;
    }
    // Dashed rings on the other controls a card names (hidden while something else has to be fixed first).
    // One ring per control id, so a ring whose part is done can fade out where it is instead of vanishing.
    const tutorialAlsoRings=new Map();
    function positionTutorialAlso(hideAll){
      const step=tutorialSteps[tutorialIndex];
      const aim=step&&!hideAll&&!tutorialCard.hidden?tutorialCurrentAim(step):null;
      // Rings go once their part is done, or once the whole card is done; the part being pointed at has the main ring.
      // A watched value inside the main ring (the number in the slider's own label row) needs no ring of its own.
      const aimLabel=aim&&!tutorialCopy(aim.id)&&document.getElementById(aim.id)?tutorialSliderLabel(document.getElementById(aim.id)):null;
      // step.keepWatch: the watch rings stay after the card is done, as the values they ring change only then.
      const done=step&&tutorialCompleted[step.id],off=hideAll||tutorialCard.hidden||!step||tutorialBlocker();
      const list=off||done&&!step.keepWatch?[]:tutorialRingList(step).filter(r=>(done?r.watch:!tutorialLatched[r.id])&&!(aim&&aim.id===r.id)&&!(r.watch&&aimLabel&&aimLabel.contains(document.getElementById(r.id))));
      // A slider the card carries is ringed twice, on the card's copy and on the real slider, so it is clear they are
      // the same control. Each ring has a key (the real one's is 'real:' + its id) and remembers its part (id) and element.
      const ids=[],withReal=(r,key)=>{const copy=tutorialCopy(r.id),real=document.getElementById(r.id);
        ids.push({...r,key:key+r.id,el:copy||real});if(copy)ids.push({...r,key:'real:'+key+r.id,el:real});};
      list.forEach(r=>withReal(r,''));
      if(aim&&tutorialCopy(aim.id))ids.push({id:aim.id,key:'real:'+aim.id,el:document.getElementById(aim.id),label:tutorialTagNumber(aim.label),watch:false});
      // Every slider part done but the card still waits (Score your fit's 'fine-tune both'): all its sliders stay ringed.
      if(!off&&!done&&step.cardSlider&&step.checks&&!tutorialPartHold&&Object.keys(step.cardSlider).every(id=>tutorialLatched[id]))
        Object.keys(step.cardSlider).forEach(id=>withReal({id,label:'',watch:false},'tune:'));
      const keep=new Set(ids.map(item=>item.key)),cardKey=tutorialTour+tutorialIndex;
      const still=hideAll||tutorialCard.hidden||calmMotion();
      tutorialAlsoRings.forEach((ring,key)=>{
        if(keep.has(key))return;
        const id=ring.dataset.part;
        // Only a part actually done on this card turns green and fades where it is (earned), and the next part's ring
        // cross-fades into the main ring (handover). Anything else (a new card, a skipped card, something in the way,
        // the tour ending) clears at once, so no ring is left floating over a spot the page has scrolled away from.
        const earned=ring.dataset.card===cardKey&&step&&(tutorialLatched[id]||tutorialCompleted[step.id]);
        const handover=!earned&&aim&&aim.id===id&&ring.dataset.card===cardKey;
        // A ring on its way out still follows its control, so it never hangs where the control used to be when the
        // page shifts under it (the slider bar shrinking in MOND mode, a heading wrapping).
        if(!(earned||handover)||still||ring.hidden||!fitTutorialRing(ring,ring.tutorialEl)){ring.remove();tutorialAlsoRings.delete(key);return;}
        const cls=handover?'is-handover':'is-leaving';
        if(!ring.classList.contains(cls)){ring.classList.add(cls);setTimeout(()=>{ring.remove();if(tutorialAlsoRings.get(key)===ring)tutorialAlsoRings.delete(key);},handover?TUTORIAL_HANDOVER_MS:TUTORIAL_RING_FADE_MS);}
      });
      tutorialGhosts=tutorialGhosts.filter(([ghost,el])=>ghost.isConnected&&(!still&&fitTutorialRing(ghost,el)||(ghost.remove(),false)));
      // While the control being pointed at has a box open over the page (the Galaxy menu's list), rings on anything
      // the box covers are hidden, so no ring is drawn across it.
      const aimEl=aim?document.getElementById(aim.id):null,open=aimEl?tutorialSpotRect(aimEl):null,shut=aimEl?tutorialBaseRect(aimEl):null;
      const box=open&&(open.width>shut.width+1||open.height>shut.height+1)?open:null;
      const underBox=el=>{if(!box||aimEl.contains(el))return false;const r=el.getBoundingClientRect();return r.left<box.right&&r.right>box.left&&r.top<box.bottom&&r.bottom>box.top;};
      ids.forEach(item=>{
        let ring=tutorialAlsoRings.get(item.key);
        const el=item.el;
        // A box the last card also had on this control carries over as it is (no restart); one already fading
        // green from that card is replaced by a fresh one.
        if(ring&&ring.dataset.card!==cardKey){if(/is-leaving|is-handover/.test(ring.className)){ring.remove();ring=null;}else ring.dataset.card=cardKey;}
        if(!ring){ring=document.createElement('div');ring.dataset.card=cardKey;ring.setAttribute('aria-hidden','true');ring.appendChild(document.createElement('span')).className='tutorial-tag';tutorialAlsoRings.set(item.key,ring);if(spotlightResize&&el)spotlightResize.observe(el);}
        ring.dataset.part=item.id;ring.tutorialEl=el;
        const cls='tutorial-also'+(item.watch?' is-watch':'')+(el&&tutorialCard.contains(el)?' in-card':'');if(ring.className!==cls)ring.className=cls;if(ring.firstChild.textContent!==item.label)ring.firstChild.textContent=item.label;
        if(ring.parentElement!==tutorialSpotlight.parentElement)tutorialSpotlight.parentElement.appendChild(ring);
        ring.hidden=!fitTutorialRing(ring,el)||underBox(el);
      });
    }
    // A step's extra rings: also = other controls to use (blue), watch = things to watch (yellow).
    // Each entry is 'elementId' or 'elementId:Label'.
    function tutorialRingList(step){
      const parse=(entry,watch)=>{const k=entry.indexOf(':');return{id:k<0?entry:entry.slice(0,k),label:k<0||watch?'':tutorialTagNumber(entry.slice(k+1)),watch};};
      return[...(step.also||[]).map(e=>parse(e,false)),...(step.watch||[]).map(e=>parse(e,true))];
    }
    // Only the step number of a tag such as '2 · Halo off' is shown on screen; the words are notes for this file.
    const tutorialTagNumber=t=>((t||'').match(/^(\d)\s·/)||[])[1]||'';
    const rectKey=r=>[r.left,r.top,r.width,r.height].map(v=>Math.round(v*2)).join();
    // Every frame, re-measure the control and re-fit the ring if it has moved or changed size (an expanding box,
    // a layout shift above it, a CSS transition), so the highlight always hugs it.
    let tutorialSpotLast='';
    // Locks: while a card is open, only the controls it names (or the one that clears what is in the way) can be
    // used, so a stray click cannot undo the step. Always free: playback speed, info (?) buttons and their boxes,
    // and the tour card. Locked controls are made inert (no clicks, keys or focus) and dimmed, then set free again
    // when the tour ends. inert is used instead of disabled, so the page's own disabled states are left alone.
    const TUTORIAL_LOCKABLE='input,select,button,textarea,a[href]';
    const TUTORIAL_ALWAYS_FREE='[data-speed],[aria-controls][aria-expanded],.ml-close,summary';
    const tutorialLockScopes=[document.querySelector('.page[data-page="simulator"]'),document.querySelector('.topbar')].filter(Boolean);
    // On a freeCamera card the camera buttons are free too, so the view can be reset, turned or pulled back.
    const TUTORIAL_CAMERA_BUTTONS=['viewReset','viewFaceOn','viewEdgeOn','viewFullHalo'];
    let tutorialLockKey=null;
    // step.freeCamera: the 3D view can be dragged, zoomed and turned as normal on this card (a function: only while true).
    const tutorialFreeCamera=step=>!!step&&(typeof step.freeCamera==='function'?step.freeCamera():!!step.freeCamera);
    function tutorialAllowed(){
      if(tutorialCard.hidden)return null;
      const step=tutorialStep(),blocker=tutorialBlocker();
      if(!blocker&&step.predict&&tutorialFlags.prediction==null)return[];
      const ids=blocker?[]:[step.target,...tutorialRingList(step).filter(r=>!r.watch).map(r=>r.id),...(step.allow||[]),...(tutorialFreeCamera(step)?TUTORIAL_CAMERA_BUTTONS:[])];
      return[...(blocker?[blocker.target]:[]),...ids.map(id=>document.getElementById(id))].filter(Boolean).map(el=>el.closest('label')||el);
    }
    function applyTutorialLocks(){
      const freeCamera=tutorialFreeCamera(tutorialCard.hidden?null:tutorialStep());
      const roots=tutorialAllowed(),key=roots&&roots.map(el=>el.id||el.className).join('|')+'|'+tutorialIndex+tutorialTour+freeCamera;
      if(key===tutorialLockKey)return;
      tutorialLockKey=key;
      controls.enabled=!roots||freeCamera;
      tutorialLockScopes.forEach(scope=>scope.querySelectorAll(TUTORIAL_LOCKABLE).forEach(el=>{
        const unit=el.closest('label')||el;
        if(tutorialCard.contains(unit))return;
        // The tour's own launch button (beside the Model setup heading) is always locked while a tour runs.
        const locked=!!roots&&(unit===tutorialLaunch||!unit.matches(TUTORIAL_ALWAYS_FREE)&&!unit.closest('#mlExplainer,#tiltExplainer,#starsExplainer,#settingsPanel')&&!roots.some(r=>r.contains(unit)||unit.contains(r)));
        if(unit.inert!==locked)unit.inert=locked;
        unit.classList.toggle('tutorial-locked',locked);
      }));
    }
    // Cuts a rounded hole in the dim layer for the main ring, each visible extra ring, and every control the card
    // names (with any info box it has opened), so parts already done and text to read stay bright too. Holes are
    // drawn into a mask, so where two overlap (an opened info box over the next card) they still stay clear.
    const tutorialDim=document.getElementById('tutorialDim'),tutorialDimHoles=document.getElementById('tutorialDimHoles');
    let tutorialDimLast='';
    function updateTutorialDim(){
      const on=!tutorialCard.hidden&&(!tutorialSpotlight.hidden||tutorialPartHold)&&!tutorialSpotlight.classList.contains('is-revisit')&&!tutorialIsFullscreen();
      tutorialDim.classList.toggle('is-on',on);
      tutorialDim.classList.toggle('is-done',tutorialSpotlight.classList.contains('is-complete'));
      if(!on)return;
      // Like the rings, a hole stops at the top menu when its control has scrolled under it.
      // The dim layer covers the whole page and scrolls with it, so holes are in page coordinates too.
      const pageHeight=Math.max(innerHeight,Math.ceil(document.body.getBoundingClientRect().bottom+scrollY))+'px';
      if(tutorialDim.style.height!==pageHeight)tutorialDim.style.height=pageHeight;
      const holes=[],step=tutorialStep(),pad=4;
      const cutHole=(el,r,k)=>{const cut=Math.max(0,tutorialCoverBottom(el)-r.top);if(cut<r.height)holes.push([{left:r.left+scrollX,top:r.top+cut+scrollY,width:r.width,height:r.height-cut},k]);};
      cutHole(tutorialAimTarget,tutorialSpotlight.getBoundingClientRect(),parseFloat(tutorialSpotlight.style.borderRadius)||0);
      tutorialAlsoRings.forEach(ring=>{if(!ring.hidden&&ring.isConnected&&ring.tutorialEl)cutHole(ring.tutorialEl,ring.getBoundingClientRect(),parseFloat(ring.style.borderRadius)||0);});
      if(!tutorialBlocker())[step.target,...tutorialRingList(step).map(r=>r.id)].forEach(id=>{
        const el=document.getElementById(id);
        if(!el||!el.getClientRects().length)return;
        const r=tutorialSpotRect(el);
        cutHole(el,{left:r.left-pad,top:r.top-pad,width:r.width+2*pad,height:r.height+2*pad},TUTORIAL_RING_MIN_RADIUS_PX);
      });
      // Open pop-ups stay bright too (the always-open View and slider panels follow the dim like the rest of the page).
      document.querySelectorAll('#mlExplainer,#tiltExplainer,#starsExplainer,#haloRadiusHelp,#haloScaleHelp,#haloConcentrationHelp,#galaxyMenuList').forEach(pop=>{
        if(tutorialPopupShown(pop))cutHole(pop,pop.getBoundingClientRect(),parseFloat(getComputedStyle(pop).borderTopLeftRadius)||0);
      });
      const html=holes.map(([r,k])=>`<rect x="${r.left}" y="${r.top}" width="${r.width}" height="${r.height}" rx="${Math.min(k,r.width/2,r.height/2)}"/>`).join('');
      if(html!==tutorialDimLast){tutorialDimLast=html;tutorialDimHoles.innerHTML=html;}
    }
    // Phones in fullscreen: the View panel and the halo slider panel would cover most of the 3D view, so while the tour
    // is open they show only the rows holding a control the card names (a row's heading goes with it). Everything
    // comes back when the card moves on, fullscreen ends or the tour closes.
    let tutorialTrimKey='';
    function trimTutorialPanels(){
      const on=!tutorialCard.hidden&&tutorialPhone()&&tutorialIsFullscreen();
      const step=on?tutorialStep():null,blocker=on?tutorialBlocker():null;
      const named=on?[blocker&&blocker.target,...(blocker?[]:[step.target,...tutorialRingList(step).map(r=>r.id)].map(id=>document.getElementById(id)))].filter(Boolean):[];
      const key=on?named.map(el=>el.id).join():'';
      if(key===tutorialTrimKey)return;
      tutorialTrimKey=key;
      const fs=document.getElementById('fsSliders');
      const rows=[...viewPanel.children].filter(el=>el!==viewToggle&&!el.matches('.view-label')).concat([...fs.querySelectorAll(':scope > :not(.fs-row), .fs-row > *')].filter(el=>el.tagName!=='LABEL'));
      rows.forEach(row=>{
        const hide=on&&!named.some(el=>row===el||row.contains(el));
        row.classList.toggle('tutorial-trim',hide);
        // A row's heading (the View panel's labels, the slider's label) sits just before it.
        const head=row.previousElementSibling;
        if(head&&(head.matches('.view-label')||head.tagName==='LABEL'))head.classList.toggle('tutorial-trim',hide);
      });
    }
    // Phones: a value the card says to watch that is off screen (χ² above the graph while the model cards sit below it;
    // the M/L and tilt tiles far up the page) is copied live onto the card, so it can be read without scrolling.
    const TUTORIAL_READOUTS={
      chiBadge:{label:'Fit',value:el=>el.innerText.replace(/\s+/g,' ').trim()},
      mlStat:{label:'M/L',value:el=>[...el.querySelectorAll('.stat-value b')].map(b=>b.textContent).join(' / ')},
      tiltStat:{label:'Tilt',value:()=>document.getElementById('tiltValue').textContent}
    };
    const tutorialReadout=document.getElementById('tutorialReadout');let tutorialReadoutKey='';
    function syncTutorialReadout(){
      const step=tutorialStep(),card=tutorialCard.getBoundingClientRect(),top=tutorialCoverBottom(tutorialCard);
      const offScreen=el=>{if(!el.getClientRects().length)return false;const r=el.getBoundingClientRect();return r.top<top||r.bottom>card.top;};
      const items=tutorialPhone()&&!tutorialBlocker()?tutorialRingList(step).filter(r=>r.watch&&TUTORIAL_READOUTS[r.id]).map(r=>[r.id,document.getElementById(r.id)]).filter(([,el])=>el&&offScreen(el)):[];
      const key=items.map(([id,el])=>id+'='+TUTORIAL_READOUTS[id].value(el)).join('|');
      if(key===tutorialReadoutKey)return;
      tutorialReadoutKey=key;tutorialReadout.hidden=!items.length;
      tutorialReadout.replaceChildren(...items.map(([id,el])=>{const chip=document.createElement('span');chip.append(Object.assign(document.createElement('b'),{textContent:TUTORIAL_READOUTS[id].label}),TUTORIAL_READOUTS[id].value(el));return chip;}));
    }
    (function trackTutorialSpot(){
      if(!tutorialCard.hidden){tutorialUpdateChecks();syncTutorialSlider(true);syncTutorialReadout();}
      else if(!tutorialReadout.hidden){tutorialReadout.hidden=true;tutorialReadoutKey='';}
      applyTutorialLocks();
      trimTutorialPanels();
      positionTutorialAlso();
      if(!tutorialCard.hidden&&tutorialAimTarget&&rectKey(tutorialSpotRect(tutorialAimTarget))!==tutorialSpotLast)updateTutorialSpotlight();
      else if(!tutorialSpotlight.hidden&&tutorialAimTarget){
        const st=tutorialSpotlight.style,fixed=tutorialSpotlight.classList.contains('is-fixed'),left=parseFloat(st.left)-(fixed?0:scrollX),top=parseFloat(st.top)-(fixed?0:scrollY);
        clipTutorialRing(tutorialSpotlight,tutorialAimTarget,{left,top,right:left+parseFloat(st.width),bottom:top+parseFloat(st.height)});
      }
      updateTutorialDim();
      requestAnimationFrame(trackTutorialSpot);
    })();
    function syncTutorialFullscreen(){
      const fullscreen=tutorialIsFullscreen();
      const parent=fullscreen?orbitPlot:document.body;
      if(tutorialCard.parentElement!==parent)parent.appendChild(tutorialCard);
      if(tutorialSpotlight.parentElement!==parent)parent.appendChild(tutorialSpotlight);
      if(tutorialCard.hidden)return;
      updateTutorialNextGate();
      tutorialRecheck();
    }
    function showTutorialStep(){
      const step=tutorialStep(),lastStep=tutorialIndex===tutorialSteps.length-1;
      const backHadFocus=document.activeElement===tutorialBack;
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialScaleSeen=[tutorialLastScale,tutorialLastScale,tutorialLastScale];
      tutorialLastTime=Number(timeSlider.value);
      tutorialDrag=null;
      // A box still showing on the control the new card starts on stays put instead of replaying its arrival.
      tutorialCarryTarget=!tutorialSpotlight.hidden&&!tutorialSpotlight.classList.contains('is-complete')?tutorialAimTarget:null;
      tutorialAimTarget=null;tutorialAimId=null;
      tutorialResetChecks();
      tutorialTarget=document.getElementById(step.target);
      // Phones: a card that uses the View panel opens it, so there is no extra toggle to find on a small screen.
      if(tutorialPhone()&&(step.viewPanel||step.fs)&&viewPanel.classList.contains('collapsed'))viewToggle.click();
      if(pageOf(tutorialTarget))showPage(pageOf(tutorialTarget),{scroll:false});
      const blocker=tutorialBlocker();
      tutorialBlockTarget=blocker?blocker.target:null;
      const part=tutorialFull?`Part ${FULL_TOUR.indexOf(tutorialTour)+1} of ${FULL_TOUR.length} · `:'';
      tutorialProgress.textContent=`${part}${TOURS[tutorialTour].label} · ${tutorialIndex+1} / ${tutorialSteps.length}`;
      tutorialProgressFill.style.width=`${100*(tutorialIndex+1)/tutorialSteps.length}%`;
      tutorialTitle.textContent=step.title;
      tutorialText.textContent=step.text;
      renderTutorialTry();
      syncTutorialSlider();
      tutorialFeedback.hidden=!tutorialCompleted[step.id];
      showTutorialFeedback(tutorialCompleted[step.id]||'',tutorialMissed[step.id]);
      tutorialCard.classList.remove('is-celebrating');
      stopTutorialAuto();
      // Back on a card already done: no green glow replayed; the ring is simply not shown.
      tutorialSpotlight.classList.toggle('is-complete',!!tutorialCompleted[step.id]);
      tutorialSpotlight.classList.toggle('is-revisit',!!tutorialCompleted[step.id]);
      document.getElementById('tutorialContent').scrollTop=0;
      tutorialBack.disabled=tutorialIndex===0;
      const nextTour=tutorialFull&&FULL_TOUR[FULL_TOUR.indexOf(tutorialTour)+1];
      tutorialNext.textContent=lastStep?(nextTour?`Next: ${TOURS[nextTour].label}`:'Finish'):'Next';
      updateTutorialNextGate();
      tutorialMore.hidden=!(lastStep&&tutorialTour==='main'&&!nextTour);
      tutorialCard.classList.toggle('on-left',innerWidth>960&&!step.cardRight);tutorialDodged=false;
      setStarMarks(step.marks);
      // The ejection steps pull the camera back far enough to frame the ring, so the escaping star stays in view.
      if(step.ring&&tutorialIsFullscreen())controls.flyTo({phi:0.35,radius:ejectionRadius*4,target:[0,0,0],rate:2.4});
      tutorialPlayedMyr=0;tutorialWatchLast=raceStartMyr=simTime;
      if(step.predict&&playing)setPlaying(false);
      syncTutorialFullscreen();
      if(tutorialBack.disabled&&backHadFocus)tutorialNext.focus({preventScroll:true});
      if(!tutorialIsFullscreen()){
        scrollTutorialTarget();
      }
      updateTutorialSpotlight();
    }
    // full: run every tour in FULL_TOUR from this one on (the guided tour from the puzzle page).
    function startTutorial(tour,full){
      if(tutorialWelcome.open)tutorialWelcome.close();
      tourOffered=true;
      tutorialTour=TOURS[tour]?tour:'main';
      tutorialFull=!!full;
      if(tutorialTour==='sandbox'&&galaxyMode!=='custom')setGalaxyMode('custom');
      tutorialSteps=tourSteps(tutorialTour,tutorialFull);
      stopTutorialAuto();clearTimeout(tutorialSwitchTimer);tutorialCard.classList.remove('leave','leave-side');
      tutorialIndex=0;
      tutorialCompleted={};
      // In the main tour the Visible matter only shortcut is revealed by the step that sets the halo to zero.
      visibleMatterOnly.hidden=tutorialTour==='main'&&!tutorialShortcutUnlocked;
      document.body.classList.add('tutorial-active');
      tutorialCard.hidden=false;
      showTutorialStep();
      tutorialTitle.focus({preventScroll:true});
    }
    // After the step's action is done (and any slider let go), a 10 s countdown moves to the next card. It pauses while a mouse button or
    // finger is held down (mid-drag), and pressing it goes straight on.
    const tutorialFeedbackText=document.getElementById('tutorialFeedbackText'),tutorialAuto=document.getElementById('tutorialAuto');
    const tutorialAutoNum=document.getElementById('tutorialAutoNum'),tutorialAutoRing=document.getElementById('tutorialAutoRing');
    const TUTORIAL_AUTO_MS=10000,TUTORIAL_RING=2*Math.PI*7;
    let tutorialAutoLeft=0,tutorialAutoLast=0,tutorialAutoFrame=0,tutorialHolding=false,tutorialSwitchTimer=0;
    function startTutorialAuto(){
      stopTutorialAuto();
      if(tutorialIndex>=tutorialSteps.length-1)return;
      tutorialAutoLeft=TUTORIAL_AUTO_MS;tutorialAutoLast=performance.now();
      tutorialAutoNum.textContent=String(Math.ceil(TUTORIAL_AUTO_MS/1000));tutorialAutoRing.style.strokeDashoffset='0';
      tutorialAuto.hidden=false;
      tutorialAutoFrame=requestAnimationFrame(tickTutorialAuto);
    }
    function stopTutorialAuto(){cancelAnimationFrame(tutorialAutoFrame);tutorialAutoFrame=0;tutorialAuto.hidden=true;}
    function tickTutorialAuto(now){
      const dt=Math.min(100,Math.max(0,now-tutorialAutoLast));tutorialAutoLast=now;
      if(!tutorialHolding&&!document.hidden)tutorialAutoLeft-=dt;
      tutorialAuto.classList.toggle('is-paused',tutorialHolding);
      if(tutorialAutoLeft<=0){stopTutorialAuto();goTutorialStep(tutorialIndex+1,1);return;}
      tutorialAutoNum.textContent=String(Math.ceil(tutorialAutoLeft/1000));
      tutorialAutoRing.style.strokeDashoffset=String(TUTORIAL_RING*(1-tutorialAutoLeft/TUTORIAL_AUTO_MS));
      tutorialAutoFrame=requestAnimationFrame(tickTutorialAuto);
    }
    tutorialAuto.addEventListener('click',()=>goTutorialStep(tutorialIndex+1,1));
    window.addEventListener('pointerdown',e=>{tutorialHolding=!tutorialCard.contains(e.target)||!!tutorialRealOf(e.target);},{capture:true,passive:true});
    // Still using the controls (another drag, a key press on a slider): the countdown starts again from 10 s,
    // so a card never moves on while the viewer is mid-experiment.
    ['pointerdown','input'].forEach(type=>window.addEventListener(type,e=>{
      if(tutorialAutoFrame&&!tutorialCard.contains(e.target)&&e.target.closest&&e.target.closest('input,button,select'))tutorialAutoLeft=TUTORIAL_AUTO_MS;
    },{capture:true,passive:true}));
    ['pointerup','pointercancel'].forEach(type=>window.addEventListener(type,()=>{
      tutorialHolding=false;
      if(tutorialPending){const pending=tutorialPending;tutorialPending=null;setTimeout(()=>completeTutorialAction(...pending),250);}
    },{capture:true,passive:true}));
    // Card-to-card animation: dir is 1 going forward, -1 going back.
    function goTutorialStep(index,dir){
      if(index<0||index>=tutorialSteps.length||tutorialCard.hidden)return;
      stopTutorialAuto();clearTimeout(tutorialSwitchTimer);
      if(calmMotion()){tutorialIndex=index;showTutorialStep();return;}
      const sideChange=!tutorialIsFullscreen()&&(innerWidth>960&&!tutorialSteps[index].cardRight)!==tutorialCard.classList.contains('on-left');
      tutorialCard.classList.remove('enter','enter-side','leave','leave-side');
      tutorialCard.style.setProperty('--dir',dir);
      tutorialCard.classList.add(sideChange?'leave-side':'leave','is-switching');
      clearTimeout(tutorialScrollTimer);tutorialScrollTimer=setTimeout(()=>tutorialCard.classList.remove('is-switching'),TUTORIAL_SWITCH_MS);
      tutorialSwitchTimer=setTimeout(()=>{
        tutorialCard.classList.remove('leave','leave-side');
        tutorialIndex=index;showTutorialStep();
        void tutorialCard.offsetWidth;
        tutorialCard.classList.add(sideChange?'enter-side':'enter');
      },sideChange?460:440);
    }
    function endTutorial(){
      stopTutorialAuto();clearTimeout(tutorialSwitchTimer);tutorialCard.classList.remove('enter','enter-side','leave','leave-side');
      tutorialTarget=null;
      tutorialCard.hidden=true;
      syncTutorialSlider();
      tutorialSpotlight.hidden=true;
      visibleMatterOnly.hidden=false;
      tutorialSpotlight.classList.remove('is-interacting');
      document.body.classList.remove('tutorial-active');
      setStarMarks([]);applyTutorialLocks();
      (tutorialIsFullscreen()?fsBtn:tutorialLaunch).focus({preventScroll:true});
    }
    ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu'].forEach(type=>tutorialCard.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
    document.addEventListener('fullscreenchange',syncTutorialFullscreen);
    document.addEventListener('webkitfullscreenchange',syncTutorialFullscreen);
    // The animated pointer only shows the way. Once the mouse (or a finger) reaches the highlighted control it fades
    // out and stays hidden for that step, so it never keeps sliding about while the user is working.
    // Listening on window in the capture phase sees the pointer even over the graphs, which stop their own events.
    function tutorialPointerReach(e){
      if(tutorialCard.hidden||tutorialSpotlight.hidden||!tutorialAimTarget)return;
      const r=tutorialAimTarget.getBoundingClientRect(),x=e.clientX,y=e.clientY;
      if(x<r.left-6||x>r.right+6||y<r.top-6||y>r.bottom+6)return;
      tutorialSpotlight.classList.add('is-reached');
      // Spot the puzzle needs all three stretches of the graph's right half (where the gap is) to be visited.
      if(tutorialStep().id==='puzzle'&&x>r.left+r.width*.5)tutorialFlags['scan'+Math.min(2,Math.floor((x-r.left-r.width*.5)/(r.width/6)))]=true;
      if(e.type==='pointerdown')completeTutorialAction('sbMonitor','Found it. Next, push the marker out of the band.');
    }
    window.addEventListener('pointermove',tutorialPointerReach,{capture:true,passive:true});
    window.addEventListener('pointerdown',tutorialPointerReach,{capture:true,passive:true});
    // Capture also catches scrolling inside panels, not just the page.
    window.addEventListener('scroll',updateTutorialSpotlight,{capture:true,passive:true});
    window.addEventListener('resize',updateTutorialSpotlight);
    // Every card should show all its text without scrolling. When it does not fit (mostly in fullscreen, where the
    // view buttons above it limit its height), the card widens step by step so the text wraps into fewer lines.
    // Phones keep their full-width card and scroll as a last resort.
    const TUTORIAL_CARD_MAX_WIDTH_PX=540,TUTORIAL_CARD_WIDEN_STEP_PX=40,TUTORIAL_CARD_SIDE_CLEARANCE_PX=290,TUTORIAL_CARD_PANEL_GAP_PX=12;
    const tutorialContentBox=document.getElementById('tutorialContent');
    function fitTutorialCard(){
      if(tutorialCard.hidden)return;
      tutorialCard.style.width='';tutorialCard.style.maxHeight='';
      if(innerWidth<=700)return;
      const overflows=()=>tutorialContentBox.scrollHeight>tutorialContentBox.clientHeight+1;
      const maxWidth=Math.min(TUTORIAL_CARD_MAX_WIDTH_PX,innerWidth-TUTORIAL_CARD_SIDE_CLEARANCE_PX);
      let width=tutorialCard.offsetWidth;
      while(overflows()&&width<maxWidth){
        width=Math.min(maxWidth,width+TUTORIAL_CARD_WIDEN_STEP_PX);
        tutorialCard.style.width=width+'px';
      }
      // Still too long in fullscreen (short screens): grow upwards into the gap below the view buttons.
      if(overflows()&&tutorialIsFullscreen()){
        const room=tutorialCard.getBoundingClientRect().bottom-tutorialViewPanel.getBoundingClientRect().bottom-TUTORIAL_CARD_PANEL_GAP_PX;
        if(room>tutorialCard.offsetHeight)tutorialCard.style.maxHeight=room+'px';
      }
    }
    const tutorialViewPanel=document.getElementById('viewPanel');
    if(window.ResizeObserver){const fitWatch=new ResizeObserver(fitTutorialCard);[tutorialContentBox,tutorialContentBox.firstElementChild,tutorialFeedback,tutorialViewPanel].forEach(el=>fitWatch.observe(el));}
    window.addEventListener('resize',fitTutorialCard);
    // Scrolled right away from the highlighted control (and anything else the card names)? Once the page stops
    // moving, glide back to it so the viewer never loses their place. Never while a finger or button is held down.
    let tutorialSnapTimer=0;
    const tutorialOnScreen=el=>{if(!el||!el.getClientRects().length)return false;const r=el.getBoundingClientRect();return r.bottom>40&&r.top<innerHeight-40;};
    function tutorialSnapBack(){
      if(tutorialCard.hidden||tutorialIsFullscreen())return;
      if(tutorialHolding){tutorialSnapTimer=setTimeout(tutorialSnapBack,600);return;}
      const step=tutorialStep(),blocker=tutorialBlocker();
      const els=blocker?[blocker.target]:[tutorialTarget,...tutorialRingList(step).map(r=>document.getElementById(r.id))];
      if(!els[0]||!els[0].getClientRects().length||els.some(tutorialOnScreen))return;
      scrollTutorialTarget(true);
      tutorialSpotlight.classList.remove('is-arriving','is-reached');void tutorialSpotlight.offsetWidth;tutorialSpotlight.classList.add('is-arriving');
    }
    window.addEventListener('scroll',()=>{
      if(tutorialCard.hidden)return;
      clearTimeout(tutorialSnapTimer);tutorialSnapTimer=setTimeout(tutorialSnapBack,1200);
    },{passive:true});
    // Slider dock: while less than SLIDER_DOCK_MIN_VISIBLE of the slider bar is on screen, the quick-controls box shows in
    // the 3D view (CSS: laptops and tablets only, never during a tour), so the sliders stay in reach beside both figures.
    const SLIDER_DOCK_MIN_VISIBLE=0.6;
    (()=>{const bar=document.querySelector('.sim-stage .controls'),plot=document.getElementById('orbitPlot');
      if(!bar||!plot||!window.IntersectionObserver)return;
      new IntersectionObserver(([entry])=>plot.classList.toggle('slider-dock',entry.intersectionRatio<SLIDER_DOCK_MIN_VISIBLE),
        {threshold:[0,SLIDER_DOCK_MIN_VISIBLE,1],rootMargin:'-56px 0px 0px 0px'}).observe(bar);})();
    const sandboxInputs=SANDBOX_SLIDERS.map(slider=>slider.input);
    // The "before" line: the curve is copied as a slider is grabbed (or first nudged from the keyboard), shown once it
    // moves, and cleared CURVE_GHOST_HOLD_MS after it is let go, or at once when the galaxy or model changes.
    const curveGhostSliders=[massSlider,scaleSlider,...sandboxInputs];
    let curveGhostArmed=false,curveGhostTimer=0,curveGhostFadeFrame=0;
    // Eases the band's opacity to `to` over ms, then calls done.
    function fadeCurveGhost(to,ms,done){
      cancelAnimationFrame(curveGhostFadeFrame);
      const plot=document.getElementById('curvePlot'),from=curveGhost.fade,t0=performance.now();
      if(calmMotion())ms=0;
      (function step(now){
        const t=ms?Math.min(1,(now-t0)/ms):1;
        curveGhost.fade=from+(to-from)*t*(2-t);
        Plotly.restyle(plot,{opacity:curveGhost.fade},[plot.data.length-1]);
        if(t<1)curveGhostFadeFrame=requestAnimationFrame(step);else if(done)done();
      })(t0);
    }
    function setCurveGhost(on){
      const plot=document.getElementById('curvePlot');
      curveGhost.leaving=false;
      if(!on){cancelAnimationFrame(curveGhostFadeFrame);curveGhost.visible=false;curveGhost.fade=0;Plotly.restyle(plot,{visible:false,opacity:0},[plot.data.length-1]);return;}
      curveGhost.visible=true;
      const t=curveGhostTrace(); // restyle reads an array as one value per trace, so each data array is wrapped
      Plotly.restyle(plot,{visible:true,opacity:curveGhost.fade,x:[t.x],y:[t.y],text:[t.text],'line.color':t.line.color,'textfont.color':t.textfont.color},[plot.data.length-1]);
      fadeCurveGhost(1,CURVE_GHOST_FADE_IN_MS);
    }
    // gently: fade out (after the hold); otherwise (new galaxy or model) clear at once.
    function clearCurveGhost(gently){
      clearTimeout(curveGhostTimer);curveGhostArmed=false;
      if(!curveGhost.visible)return;
      if(gently===true){curveGhost.leaving=true;fadeCurveGhost(0,CURVE_GHOST_FADE_OUT_MS,()=>setCurveGhost(false));}
      else setCurveGhost(false);
    }
    function armCurveGhost(e){
      const slider=tutorialRealOf(e.target)||e.target;
      if(!curveGhostSliders.includes(slider)||e.type==='keydown'&&curveGhost.visible&&!curveGhost.leaving)return;
      const halo=slider===massSlider||slider===scaleSlider,line=document.getElementById('curvePlot').data[halo?3:1];
      curveGhost.x=line.x.slice();curveGhost.y=line.y.slice();curveGhost.rgb=halo?'255,174,92':'95,208,224';
      clearTimeout(curveGhostTimer);curveGhostArmed=true;
    }
    window.addEventListener('pointerdown',armCurveGhost,{capture:true,passive:true});
    window.addEventListener('keydown',armCurveGhost,{capture:true,passive:true});
    curveGhostSliders.forEach(slider=>{
      slider.addEventListener('input',()=>{if(curveGhostArmed&&(!curveGhost.visible||curveGhost.leaving))setCurveGhost(true);clearTimeout(curveGhostTimer);});
      slider.addEventListener('change',()=>{clearTimeout(curveGhostTimer);curveGhostTimer=setTimeout(()=>clearCurveGhost(true),CURVE_GHOST_HOLD_MS);});
    });
    ['galaxychange','galaxymodechange'].forEach(type=>document.addEventListener(type,clearCurveGhost));
    document.getElementById('mondToggle').addEventListener('change',clearCurveGhost);
    [massSlider,scaleSlider,...sandboxInputs].forEach(slider=>{
      slider.addEventListener('pointerdown',()=>tutorialSpotlight.classList.add('is-interacting'));
      slider.addEventListener('change',updateTutorialSpotlight);
    });
    // The stability step is done once the halo sliders have pushed the concentration outside the realistic range.
    function checkStabilityStep(){
      const state=stabilityReport().state;
      if(state==='diffuse'||state==='dense')completeTutorialAction('sbStability','Warning shown: this halo is outside the range seen for real halos.');
    }
    massSlider.addEventListener('input',()=>{
      const value=Number(massSlider.value);
      const crossedZero=(value===0)!==(tutorialLastMass===0);
      if(!mondMode){
        if(!tutorialCard.hidden&&tutorialTour==='main'&&value===0&&!tutorialShortcutUnlocked){
          tutorialShortcutUnlocked=true;
          visibleMatterOnly.hidden=false;
          visibleMatterOnly.classList.add('tutorial-reveal');
        }
        if(value===0)completeTutorialAction('removeHalo','Halo removed: the orange line now shows visible matter only.');
        checkStabilityStep();
      }
      tutorialLastMass=value;
      if(crossedZero)refreshTutorialGuidance();
    });
    // Scrolling waits until the slider is released, so the page never jumps mid-drag.
    massSlider.addEventListener('change',tutorialRecheck);
    scaleSlider.addEventListener('input',()=>{
      const value=Number(scaleSlider.value);
      if(!mondMode&&value!==tutorialLastScale){
        // 'Move the scale radius both ways': done once it has gone below and above where the card started
        // (or reached an end of the track on the other side).
        const[start]=tutorialScaleSeen,lo=Number(scaleSlider.min),hi=Number(scaleSlider.max);
        tutorialScaleSeen[1]=Math.min(tutorialScaleSeen[1],value);tutorialScaleSeen[2]=Math.max(tutorialScaleSeen[2],value);
        const down=tutorialScaleSeen[1]<start||start<=lo,up=tutorialScaleSeen[2]>start||start>=hi;
        if(down&&up)completeTutorialAction('haloShape','Both ways: the same mass, arranged differently.');
        else if(down&&tutorialStep().id==='haloShape')renderTutorialTry();
        checkStabilityStep();
      }
      tutorialLastScale=value;
    });
    document.querySelectorAll('#haloMetrics .halo-metric').forEach(metric=>{
      metric.addEventListener('explainerread',()=>{tutorialFlags[metric.id]=true;});
    });
    visibleMatterOnly.addEventListener('click',()=>{
      if(!mondMode)completeTutorialAction('removeHalo','Halo removed: the orange line now shows visible matter only.');
      tutorialLastMass=Number(massSlider.value);
      tutorialRecheck();
    });
    document.getElementById('applyObservedFit').addEventListener('click',()=>{
      tutorialFlags.typical=true;
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialRecheck();
    });
    timeSlider.addEventListener('input',()=>{tutorialLastTime=Number(timeSlider.value);});
    // Opening View, or a layer or camera button, can clear what a step was waiting for.
    viewToggle.addEventListener('click',tutorialRecheck);
    ['toggleGas','toggleHalo','toggleOrbits','toggleEjection','toggleMondZone','viewFaceOn','viewEdgeOn','viewReset'].forEach(id=>document.getElementById(id).addEventListener('click',()=>{requestRender();refreshTutorialGuidance();}));
    fullHaloBtn.addEventListener('click',()=>{tutorialFlags.fullHalo=true;});
    // Restart (fullscreen corner panel): back to 0 Myr and play, so a new halo can be watched from the start.
    document.getElementById('fsRestart').addEventListener('click',()=>{tutorialFlags[Number(massSlider.value)===0?'restart0':'restartHalo']=true;onSliderInput(0);setPlaying(true);});
    document.getElementById('mlInfoBtn').addEventListener('explainerread',()=>completeTutorialAction('ml','Explainer read.'));
    const SANDBOX_STEP_FEEDBACK={
      sbDiskMass:['sbDisk','Disk mass changed: the typical halo followed it.'],
      sbGasMass:['sbGas','Gas mass changed.']
    };
    sandboxInputs.forEach(input=>input.addEventListener('input',()=>{
      const feedback=SANDBOX_STEP_FEEDBACK[input.id];
      if(feedback)completeTutorialAction(feedback[0],feedback[1]);
      if(input.id==='sbBulgeMass'&&Number(input.value)>0)completeTutorialAction('sbBulge','Bulge added: the inner orbits speed up.');
      completeTutorialAction('sbSummary','Totals updated.');
    }));
    window.addEventListener('pointerup',()=>tutorialSpotlight.classList.remove('is-interacting'));
    window.addEventListener('pointercancel',()=>tutorialSpotlight.classList.remove('is-interacting'));
    mondToggle.addEventListener('change',()=>{
      if(mondMode)completeTutorialAction('mond','MOND selected: no halo, modified gravity instead.');
      tutorialRecheck();
    });
    bestFitToggle.addEventListener('change',()=>{
      if(bestFitToggle.checked)completeTutorialAction('bestFit','Best fits on: each model has re-tuned the mass-to-light ratio, distance and tilt (and MOND its a₀) within realistic limits.');
    });
    document.getElementById('tutorialStart').addEventListener('click',()=>startTutorial('main'));
    document.getElementById('tutorialWelcomeSkip').addEventListener('click',()=>tutorialWelcome.close());
    tutorialLaunch.addEventListener('click',()=>startTutorial(galaxyMode==='custom'?'sandbox':'main'));
    document.getElementById('puzzleTour').addEventListener('click',()=>{if(galaxyMode==='custom')setGalaxyMode('preset');startTutorial('main',true);});
    document.getElementById('tutorialSkip').addEventListener('click',endTutorial);
    tutorialBack.addEventListener('click',()=>{if(tutorialIndex>0)goTutorialStep(tutorialIndex-1,-1);});
    tutorialNext.addEventListener('click',()=>{
      if(tutorialIndex===tutorialSteps.length-1){
        const nextTour=tutorialFull&&FULL_TOUR[FULL_TOUR.indexOf(tutorialTour)+1];
        if(!nextTour){endTutorial();showTutorialDone();return;}
        // The sandbox tour ends in fullscreen but Wrapping up is on the page: leave fullscreen first, and start the
        // next tour only once the exit has finished, so it never opens half-way through the switch.
        if(!tutorialIsFullscreen()||TOURS[nextTour].steps[0].fs){startTutorial(nextTour,true);return;}
        let started=false;const go=()=>{if(started)return;started=true;startTutorial(nextTour,true);};
        if(nativeFsElement()){document.addEventListener('fullscreenchange',()=>requestAnimationFrame(go),{once:true});setTimeout(go,600);}
        toggleFullscreen();
        if(!nativeFsElement())go();
        return;
      }
      goTutorialStep(tutorialIndex+1,1);
    });
    // Finishing a tour (Finish on its last card, not Skip) shows a short "Tutorial complete" celebration. It sits inside
    // the 3D view when that is fullscreen, so it is never hidden behind it, and closes on its button, a click outside
    // it, Esc, or by itself after TUTORIAL_DONE_MS.
    const tutorialDone=document.getElementById('tutorialDone'),TUTORIAL_DONE_MS=9000;let tutorialDoneTimer=0;
    function showTutorialDone(){
      (tutorialIsFullscreen()?orbitPlot:document.body).appendChild(tutorialDone);
      tutorialDone.classList.remove('is-leaving');tutorialDone.hidden=false;
      document.getElementById('tutorialDoneClose').focus({preventScroll:true});
      clearTimeout(tutorialDoneTimer);tutorialDoneTimer=setTimeout(closeTutorialDone,TUTORIAL_DONE_MS);
    }
    function closeTutorialDone(){
      if(tutorialDone.hidden||tutorialDone.classList.contains('is-leaving'))return;
      clearTimeout(tutorialDoneTimer);tutorialDone.classList.add('is-leaving');
      setTimeout(()=>{tutorialDone.hidden=true;tutorialDone.classList.remove('is-leaving');},calmMotion()?0:400);
    }
    document.getElementById('tutorialDoneClose').addEventListener('click',closeTutorialDone);
    tutorialDone.addEventListener('click',e=>{if(e.target===tutorialDone)closeTutorialDone();});
    document.addEventListener('keydown',e=>{if(e.key==='Escape'&&!tutorialDone.hidden){e.stopPropagation();closeTutorialDone();}},{capture:true});
    // The main tour's last step offers the sandbox tour; pressing it opens the Custom Sandbox.
    tutorialMore.addEventListener('click',()=>startTutorial('sandbox'));
    // Left/Right arrows: move between tour cards while the tour is open; otherwise, in fullscreen only, step the
    // simulation time by ARROW_TIME_STEP_MYR (hold to scrub).
    // While the tour is open, Left/Right always turn the card, even when the last control used (a slider, the Galaxy
    // menu, a button) still has focus; only a box being typed in keeps them. Up/Down still move a focused slider.
    // Runs in the capture phase so the focused control never sees the key as well.
    const ARROW_TIME_STEP_MYR=50;
    document.addEventListener('keydown',e=>{
      if(e.altKey||e.ctrlKey||e.metaKey||e.shiftKey)return;
      if(e.key!=='ArrowLeft'&&e.key!=='ArrowRight')return;
      const el=e.target instanceof Element?e.target:null;
      const dir=e.key==='ArrowRight'?1:-1;
      if(!tutorialCard.hidden){
        if(el&&el.closest('textarea,[contenteditable]:not([contenteditable="false"]),input:not([type=range],[type=checkbox],[type=radio],[type=button],[type=submit])'))return;
        e.preventDefault();e.stopPropagation();
        if(!e.repeat)(dir>0?tutorialNext:tutorialBack).click();
        return;
      }
      if(el&&el.closest('input,select,textarea,button,a,summary,[contenteditable]:not([contenteditable="false"]),[role="slider"],[role="radio"]'))return;
      if(tutorialIsFullscreen()){
        e.preventDefault();
        onSliderInput(simTime+dir*ARROW_TIME_STEP_MYR);
      }
    },true);
    // The tour is offered the first time the Simulator page opens in each page load (a refresh or new tab offers it
    // again; Will's choice), not over the opening story. Nothing is stored between visits.
    let tourOffered=false;
    document.addEventListener('pagechange',e=>{
      tutorialRecheck();
      if(e.detail.page==='simulator'&&!tourOffered&&tutorialCard.hidden){tourOffered=true;tutorialWelcome.showModal();}
    });
    followHash(location.hash);
  </script>
</body>
</html>
"""


def join_names(names: list[str]) -> str:
    """Join names for prose: 'A', 'A and B', 'A, B and C'."""
    if len(names) < 2:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def tuned_count_text(galaxies: list[dict[str, Any]], model: str, fit_key: str) -> str:
    """k for one best fit across the presets, e.g. '5 (6 for NGC 7331)'. Matches tunedCount (page JS)."""
    counts = {g["name"]: MODEL_TUNED_VALUES[model] + g["best_fit"][fit_key]["tuned_settings"] for g in galaxies}
    usual = max(set(counts.values()), key=list(counts.values()).count)
    others = sorted({k for k in counts.values() if k != usual})
    extra = "; ".join(f"{k} for {join_names([n for n, c in counts.items() if c == k])}" for k in others)
    return f"{usual} ({extra})" if extra else str(usual)


def main() -> None:
    if not SPARC_PATH.is_file():
        raise FileNotFoundError(f"Missing SPARC table: {SPARC_PATH}")

    galaxy_info = load_galaxy_info()
    galaxies = [build_galaxy_record(preset, galaxy_info) for preset in GALAXY_PRESETS]
    galaxies_js = render_galaxies_js(galaxies)
    GALAXIES_JS_PATH.write_text(galaxies_js, encoding="utf-8")
    payload = {
        "upsilon_disk": UPSILON_DISK,
        "upsilon_bulge": UPSILON_BULGE,
        "default_galaxy_id": DEFAULT_PRESET.sparc_id,
        "stable_concentration": list(STABLE_CONCENTRATION),
        "model_tuned_values": MODEL_TUNED_VALUES,
        "particle_model": {
            "base_count": BASE_STAR_COUNT,
            "min_count": MIN_STAR_COUNT,
            "max_count": MAX_STAR_COUNT,
            "count_exponent": STAR_COUNT_EXPONENT,
            "ejection_factor": EJECTION_RADIUS_FACTOR,
            "softening_kpc": list(SOFTENING_RANGE_KPC),
        },
    }
    html = HTML_TEMPLATE.replace(
        "__DATA_JSON__",
        json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/"),
    )
    # The deployed site is index.html alone, so the galaxies.js code is inlined into it.
    html = html.replace("__GALAXIES_JS__", galaxies_js.replace("</", "<\\/"))
    html = html.replace("__DEFAULT_GALAXY_NAME__", escape(DEFAULT_PRESET.name))
    html = html.replace(
        "__GALAXY_NAMES__", escape(join_names([preset.name for preset in GALAXY_PRESETS]))
    )
    html = html.replace("__FRAME_MAX__", "1000")
    html = html.replace("__MOND_FORMULA_ROWS__", mond_formula_table(galaxies, galaxy_info))
    html = html.replace("__K_NFW_FIT__", escape(tuned_count_text(galaxies, "nfw", "nfw")))
    html = html.replace("__K_MOND_FIT__", escape(tuned_count_text(galaxies, "mondBest", "mond")))
    html = html.replace("__FIT_LIMIT_SIGMA__", f"{FIT_LIMIT_SIGMA:g}")
    html = html.replace("__A0_MIN__", f"{MOND_A0_RANGE_MS2[0] / 1e-10:.2f}")
    html = html.replace("__A0_MAX__", f"{MOND_A0_RANGE_MS2[1] / 1e-10:.2f}")
    OUTPUT_PATH.write_text(html, encoding="utf-8")
    size_mb = OUTPUT_PATH.stat().st_size / (1024 * 1024)
    print(f"Generated {OUTPUT_PATH.name} ({size_mb:.2f} MiB) and {GALAXIES_JS_PATH.name}")
    print(f"Embedded {len(galaxies)} SPARC galaxies")
    for galaxy in galaxies:
        std, nfw, mond = (galaxy["best_fit"][key] for key in ("nfw_standard", "nfw", "mond"))
        print(
            f"  {galaxy['name']}: {len(galaxy['rotation']['radius_kpc'])} points; "
            f"NFW standard M200={std['m200_msun']:.3g} Msun, rs={std['scale_radius_kpc']:.1f} kpc, chi2={std['chi2']:.1f}; "
            f"NFW best M200={nfw['m200_msun']:.3g} Msun, rs={nfw['scale_radius_kpc']:.1f} kpc, chi2={nfw['chi2']:.1f}; "
            f"MOND a0={mond['a0_ms2']:.3g}, chi2={mond['chi2']:.1f}"
        )


if __name__ == "__main__":
    main()