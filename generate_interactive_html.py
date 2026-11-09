"""Create a standalone Plotly dashboard from the local SPARC galaxy data."""

from __future__ import annotations

import csv
import json
import math
from html import escape
from pathlib import Path
from typing import Any, NamedTuple


BASE_DIR = Path(__file__).resolve().parent
TRAJECTORY_PATH = BASE_DIR / "results" / "star_trajectories_baryonic.csv"
SPARC_PATH = BASE_DIR / "table2.dat"
OUTPUT_PATH = BASE_DIR / "index.html"
GALAXIES_JS_PATH = BASE_DIR / "galaxies.js"

G_KPC_KMS2_PER_MSUN = 4.30091e-6
H0_KMS_PER_KPC = 0.07
UPSILON_DISK = 0.50
UPSILON_BULGE = 0.70
STAR_COUNT = 200


class GalaxyPreset(NamedTuple):
    """One entry of the galaxy dropdown."""

    sparc_id: str  # must match the Name column of table2.dat
    name: str
    summary: str


# The first preset is selected when the page opens. To offer another SPARC galaxy,
# add a preset here and run this script again.
GALAXY_PRESETS = (
    GalaxyPreset("NGC3198", "NGC 3198", "Classic late-type spiral benchmark"),
    GalaxyPreset("NGC2403", "NGC 2403", "Gas-rich late-type spiral"),
    GalaxyPreset("NGC7331", "NGC 7331", "Massive early-type spiral with a prominent bulge"),
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


def load_trajectories(path: Path) -> list[dict[str, Any]]:
    by_star: dict[int, dict[str, Any]] = {}
    with path.open("r", newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            star_id = int(row["star_id"])
            star = by_star.setdefault(
                star_id,
                {
                    "id": star_id,
                    "v_rot_kms": float(row["v_obs_kms"]),
                    "initial_r_kpc": float(row["initial_r_kpc"]),
                    "frames": [],
                },
            )
            star["frames"].append(
                (
                    int(row["frame"]),
                    float(row["x_kpc"]),
                    float(row["y_kpc"]),
                    float(row["z_kpc"]),
                )
            )

    if len(by_star) < STAR_COUNT:
        raise ValueError(f"Expected at least {STAR_COUNT} stars in {path}")

    stars = []
    expected_frames: int | None = None
    for star_id in sorted(by_star):
        star = by_star[star_id]
        star["frames"].sort(key=lambda frame: frame[0])
        if expected_frames is None:
            expected_frames = len(star["frames"])
        if len(star["frames"]) != expected_frames:
            raise ValueError(f"Star {star_id} has an incomplete trajectory in {path}")
        star["positions"] = [list(frame[1:]) for frame in star.pop("frames")]
        stars.append(star)

    return stars[:STAR_COUNT]


def load_sparc_galaxy(path: Path, sparc_id: str) -> dict[str, Any]:
    """Read one galaxy's rotation curve from a SPARC table (the table2.dat column layout)."""
    distance_mpc = None
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as source:
        for line in source:
            columns = line.split()
            if len(columns) < 10 or columns[0].upper() != sparc_id.upper():
                continue
            try:
                distance, radius, observed, error, gas, disk, bulge = map(float, columns[1:8])
            except ValueError:
                continue
            distance_mpc = distance
            # Vgas keeps its sign: the gas term is squared with its sign preserved.
            rows.append((radius, abs(observed), abs(error), gas, abs(disk), abs(bulge)))

    if distance_mpc is None or len(rows) < 3:
        raise ValueError(f"No usable {sparc_id} SPARC points found in {path}")
    rows.sort(key=lambda row: row[0])
    rotation = {name: list(values) for name, values in zip(ROTATION_COLUMNS, zip(*rows))}
    return {"distance_mpc": distance_mpc, "rotation": rotation}


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


def fit_default_halo(curve: list[dict[str, float]]) -> tuple[float, float]:
    """Grid-fit a usable M200 and scale-radius starting point to SPARC."""
    critical_density = (
        3.0 * H0_KMS_PER_KPC**2 / (8.0 * math.pi * G_KPC_KMS2_PER_MSUN)
    )

    def score(mass_billion_solar: float, scale_radius: float) -> float:
        mass = mass_billion_solar * 1e9
        r200 = (3.0 * mass / (4.0 * math.pi * 200.0 * critical_density)) ** (1.0 / 3.0)
        concentration = max(r200 / scale_radius, 1e-8)
        shape_concentration = max(_nfw_shape(concentration), 1e-12)
        chi_squared = 0.0
        for point in curve:
            radius = point["radius_kpc"]
            enclosed_mass = mass * _nfw_shape(radius / scale_radius) / shape_concentration
            if radius >= r200:
                enclosed_mass = mass
            dark_velocity = math.sqrt(G_KPC_KMS2_PER_MSUN * enclosed_mass / radius)
            total_velocity = math.hypot(point["baryonic_kms"], dark_velocity)
            uncertainty = max(point["error_kms"], 1.0)
            chi_squared += ((point["observed_kms"] - total_velocity) / uncertainty) ** 2
        return chi_squared

    coarse = min(
        (
            (score(mass, scale_radius), mass, scale_radius)
            for mass in range(100, 2001, 20)
            for scale_radius in range(5, 61, 1)
        ),
        key=lambda item: item[0],
    )
    _, coarse_mass, coarse_radius = coarse
    fine = min(
        (
            (score(mass, scale_radius), mass, scale_radius)
            for mass in range(max(100, coarse_mass - 40), min(2000, coarse_mass + 40) + 1, 2)
            for scale_radius in range(max(5, coarse_radius - 3), min(60, coarse_radius + 3) + 1)
        ),
        key=lambda item: item[0],
    )
    return fine[1] * 1e9, float(fine[2])


def fit_galaxy_halo(rotation: dict[str, list[float]]) -> tuple[float, float]:
    """Return one galaxy's starting NFW halo: M200 in Msun and scale radius in kpc.

    The fit is snapped to the sliders' resolution so the page can display it exactly.
    """
    curve = []
    for radius, observed, error, gas, disk, bulge in zip(
        *(rotation[name] for name in ROTATION_COLUMNS)
    ):
        baryonic_squared = gas * abs(gas) + UPSILON_DISK * disk**2 + UPSILON_BULGE * bulge**2
        curve.append(
            {
                "radius_kpc": radius,
                "observed_kms": observed,
                "error_kms": error,
                "baryonic_kms": math.sqrt(max(0.0, baryonic_squared)),
            }
        )
    mass_msun, scale_radius_kpc = fit_default_halo(curve)
    return (
        round(mass_msun / HALO_MASS_STEP_MSUN) * HALO_MASS_STEP_MSUN,
        round(scale_radius_kpc / HALO_SCALE_STEP_KPC) * HALO_SCALE_STEP_KPC,
    )


def build_galaxy_record(preset: GalaxyPreset, sparc_path: Path) -> dict[str, Any]:
    """Collect everything the page needs to model one preset galaxy."""
    galaxy = load_sparc_galaxy(sparc_path, preset.sparc_id)
    mass_msun, scale_radius_kpc = fit_galaxy_halo(galaxy["rotation"])
    return {
        "id": preset.sparc_id,
        "name": preset.name,
        "summary": preset.summary,
        "distance_mpc": galaxy["distance_mpc"],
        "default_halo": {"m200_msun": mass_msun, "scale_radius_kpc": scale_radius_kpc},
        "rotation": galaxy["rotation"],
    }


GALAXIES_JS_HEADER = """/*
 * galaxies.js - SPARC rotation-curve presets for the Interactive Galaxy Rotation Simulator.
 *
 * Generated by generate_interactive_html.py from table2.dat (SPARC mass models: Lelli,
 * McGaugh & Schombert 2016, AJ 152, 157). Do not edit by hand: add or change presets in
 * GALAXY_PRESETS in the generator and run it again. The generator also inlines this code
 * into index.html, so the published page does not load this file separately.
 *
 * Radii are in kpc and velocities in km/s. v_gas_kms keeps SPARC's sign (the gas term can be
 * negative). default_halo is the generator's NFW fit to the observed curve with the page's
 * fixed mass-to-light ratios: the starting point of the Apply Observed Fit button.
 */
"""


def render_galaxies_js(galaxies: list[dict[str, Any]]) -> str:
    """Render the galaxies.js registry: a JavaScript lookup keyed by SPARC id."""

    def dump(value: Any) -> str:
        return json.dumps(value, allow_nan=False)

    lines = [GALAXIES_JS_HEADER.rstrip("\n"), "window.GALAXIES = {"]
    for galaxy in galaxies:
        halo = galaxy["default_halo"]
        lines += [
            f"  {dump(galaxy['id'])}: {{",
            f"    id: {dump(galaxy['id'])},",
            f"    name: {dump(galaxy['name'])},",
            f"    summary: {dump(galaxy['summary'])},",
            f"    distance_mpc: {dump(galaxy['distance_mpc'])},",
            "    default_halo: {",
            f"      m200_msun: {dump(halo['m200_msun'])},",
            f"      scale_radius_kpc: {dump(halo['scale_radius_kpc'])},",
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
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
  <style>
    :root { color-scheme: dark; --bg:#050910; --panel:#0a111b; --line:#1b2a3a; --ink:#dbe7f4; --muted:#91a4b8; --cyan:#7be4f1; --gold:#f2c879; }
    * { box-sizing:border-box; }
    body { margin:0; background:var(--bg); color:var(--ink); font:14px/1.45 "Segoe UI",sans-serif; }
    header { display:flex; align-items:flex-end; justify-content:space-between; gap:20px; padding:22px 26px 16px; border-bottom:1px solid var(--line); background:#070d15; }
    .eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.12em; font-size:11px; }
    h1 { font-size:22px; line-height:1.15; margin:5px 0 0; font-weight:600; }
    .header-stats { display:flex; flex-wrap:nowrap; justify-content:flex-end; gap:8px; flex:0 0 auto; align-items:flex-start; }
    .stat { min-width:0; padding:7px 12px 8px; border:1px solid var(--line); border-radius:8px; background:linear-gradient(160deg,rgba(123,228,241,0.09),rgba(123,228,241,0.015)); white-space:nowrap; }
    .stat.ml { background:linear-gradient(160deg,rgba(242,200,121,0.12),rgba(242,200,121,0.02)); border-color:rgba(242,200,121,0.28); }
    .stat-label { color:var(--muted); font-size:9px; letter-spacing:.1em; text-transform:uppercase; }
    .stat-value { margin-top:1px; font-size:18px; font-weight:600; line-height:1.15; color:var(--cyan); font-variant-numeric:tabular-nums; }
    .stat-value small { font-size:12px; font-weight:400; color:var(--muted); margin-left:3px; }
    .stat.ml .stat-value { color:var(--gold); }
    .ml-col { display:flex; flex-direction:column; align-items:stretch; gap:6px; }
    .fit-toggle { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:5px 10px; border:1px solid rgba(255,174,92,0.28); border-radius:8px; background:rgba(255,174,92,0.05); font-size:11px; color:var(--muted); cursor:pointer; user-select:none; white-space:nowrap; }
    .fit-toggle[hidden] { display:none; }
    .fit-toggle input { position:absolute; opacity:0; width:0; height:0; }
    .fit-toggle .track { position:relative; width:30px; height:16px; border-radius:999px; background:#16222f; border:1px solid rgba(120,150,180,0.4); transition:background .2s,border-color .2s; flex:none; }
    .fit-toggle .track::after { content:""; position:absolute; top:2px; left:2px; width:10px; height:10px; border-radius:50%; background:#dbe7f4; transition:transform .2s; }
    .fit-toggle input:checked ~ .track { background:rgba(255,174,92,0.25); border-color:rgba(255,174,92,0.7); }
    .fit-toggle input:checked ~ .track::after { transform:translateX(14px); background:#ffae5c; }
    .fit-toggle input:focus-visible ~ .track { outline:2px solid var(--cyan); outline-offset:2px; }
    .fit-toggle input:checked ~ .fit-label { color:var(--ink); }
    .ml-pair { display:flex; gap:16px; }
    .ml-pair span { display:block; font-size:10px; font-weight:400; color:var(--muted); letter-spacing:.06em; }
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
    .galaxy-bar select { min-width:150px; padding:6px 32px 6px 12px; font:inherit; font-size:13px; color:var(--ink); background-color:#102130; background-image:linear-gradient(45deg,transparent 50%,var(--muted) 50%),linear-gradient(135deg,var(--muted) 50%,transparent 50%); background-position:calc(100% - 17px) 55%,calc(100% - 12px) 55%; background-size:5px 5px,5px 5px; background-repeat:no-repeat; border:1px solid #294255; border-radius:4px; -webkit-appearance:none; appearance:none; cursor:pointer; }
    .galaxy-bar select:hover { border-color:var(--cyan); }
    .galaxy-bar select:focus-visible { outline:2px solid var(--cyan); outline-offset:2px; }
    .galaxy-bar option { background:#0a111b; color:var(--ink); }
    .galaxy-summary { color:var(--muted); font-size:12px; }
    .preset-picker { display:flex; align-items:center; gap:16px; flex-wrap:wrap; }
    .preset-picker[hidden], .sandbox[hidden], .stability-chip[hidden] { display:none; }
    .source-toggle { display:inline-grid; grid-template-columns:1fr 1fr; padding:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .source-toggle label { position:relative; cursor:pointer; }
    .source-toggle input { position:absolute; inset:0; width:100%; height:100%; margin:0; opacity:0; cursor:pointer; }
    .source-toggle span { display:block; padding:5px 14px; border-radius:999px; color:var(--muted); font-size:12px; text-align:center; white-space:nowrap; transition:background .15s,color .15s; }
    .source-toggle label:hover input:not(:checked) + span { color:var(--ink); }
    .source-toggle input:checked + span { background:rgba(123,228,241,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(123,228,241,.5); }
    .source-toggle input:focus-visible + span { outline:2px solid var(--cyan); outline-offset:2px; }
    .sandbox { margin-bottom:14px; padding:13px 16px 14px; border:1px solid var(--line); background:var(--panel); border-radius:5px; }
    .sandbox-head { display:flex; align-items:baseline; justify-content:space-between; gap:6px 16px; flex-wrap:wrap; margin-bottom:12px; }
    .sandbox-head p { margin:0; color:var(--muted); font-size:12px; }
    .sandbox-groups { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px 22px; }
    .sandbox-group { display:flex; flex-direction:column; gap:10px; min-width:0; }
    .sandbox-group h3 { margin:0; color:var(--cyan); font-size:11px; font-weight:600; letter-spacing:.1em; text-transform:uppercase; }
    .sandbox .control-help { min-height:0; }
    .sandbox-summary { margin:0; display:grid; grid-template-columns:1fr auto; gap:3px 12px; padding:9px 11px; border:1px solid rgba(120,150,180,0.2); border-radius:8px; background:rgba(120,150,180,0.05); font-size:12px; }
    .sandbox-summary dt { color:var(--muted); }
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
    .stability-chip { position:absolute; top:64px; left:12px; z-index:100; max-width:min(300px,calc(100% - 190px)); padding:5px 10px; border-radius:8px; border:1px solid rgba(255,138,128,0.6); background:rgba(40,10,10,0.6); -webkit-backdrop-filter:blur(10px); backdrop-filter:blur(10px); color:#ffb4ac; font-size:11px; font-weight:600; line-height:1.35; pointer-events:none; }
    @media (prefers-reduced-motion:reduce) { .c-marker { transition:none; } }
    @media (max-width:960px) { .sandbox-groups { grid-template-columns:1fr 1fr; } }
    @media (max-width:560px) { .sandbox-groups { grid-template-columns:1fr; } .preset-picker { gap:10px; } }
    .learn { margin-top:14px; display:grid; grid-template-columns:1fr 1fr; gap:14px; }
    .learn details { background:var(--panel); border:1px solid var(--line); border-radius:5px; padding:0 16px; }
    .learn summary { cursor:pointer; padding:12px 0; font-size:14px; font-weight:600; list-style:none; display:flex; justify-content:space-between; align-items:center; }
    .learn summary::-webkit-details-marker { display:none; }
    .learn summary::after { content:"▾"; color:var(--muted); transition:transform .2s; }
    .learn details[open] summary::after { transform:rotate(180deg); }
    .learn details.nfw summary { color:#b39ae8; } .learn details.mond summary { color:#ffae5c; }
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
    .key .g { color:#8fe3ee; } .key .o { color:#ffc48a; } .key .r { color:#ff9d94; }
    #chiBadge { font-variant-numeric:tabular-nums; padding:3px 10px; border-radius:999px; border:1px solid rgba(120,150,180,0.3); background:rgba(120,150,180,0.08); color:#c3cfdb; font-size:12px; }
    #chiBadge.good { border-color:rgba(95,208,224,0.55); color:#8fe3ee; } #chiBadge.ok { border-color:rgba(255,174,92,0.55); color:#ffc48a; } #chiBadge.poor { border-color:rgba(255,138,128,0.6); color:#ff9d94; }
    .plot { width:100%; height:min(68vh,700px); min-height:420px; }
    #orbitPlot { position:relative; touch-action:none; overscroll-behavior:contain; user-select:none; -webkit-user-select:none; cursor:grab; }
    #orbitPlot.dragging { cursor:grabbing; }
    #orbitPlot canvas { display:block; width:100% !important; height:100% !important; touch-action:none; }
    #orbitPlot:fullscreen, #orbit-plot-container:fullscreen { width:100vw; height:100vh; max-height:none; border:0; border-radius:0; background:#020206; }
    #orbitPlot:-webkit-full-screen, #orbit-plot-container:-webkit-full-screen { width:100vw; height:100vh; max-height:none; border:0; border-radius:0; background:#020206; }
    .timeline-row { display:flex; align-items:center; gap:12px; }
    .timeline-row input[type=range] { flex:1; min-width:0; }
    .play-btn { display:inline-flex; align-items:center; justify-content:center; gap:6px; min-width:44px; min-height:44px; padding:0 14px; flex:none; font-size:13px; border-radius:22px; background:rgba(123,228,241,.12); border-color:var(--cyan); touch-action:manipulation; }
    .play-btn[aria-pressed=true] { background:rgba(123,228,241,.28); }
    @media(max-width:560px) { .play-btn .play-text { display:none; } .play-btn { width:44px; padding:0; } }
    .fullscreen-active { position:fixed !important; inset:0; z-index:9999; width:100vw !important; height:100dvh !important; max-height:none !important; min-height:0 !important; border:0 !important; border-radius:0 !important; background:#020206; }
    body.fullscreen-lock { overflow:hidden; }
    .fs-sliders { display:none; position:absolute; left:12px; bottom:12px; z-index:3; width:min(240px,calc(100% - 24px)); padding:8px 10px; gap:4px; flex-direction:column; border-radius:10px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.55); -webkit-backdrop-filter:blur(14px); backdrop-filter:blur(14px); color:var(--ink); font-size:11px; cursor:default; }
    #orbitPlot:fullscreen .fs-sliders, #orbitPlot:-webkit-full-screen .fs-sliders, #orbitPlot.fullscreen-active .fs-sliders { display:flex; }
    .fs-sliders label { display:flex; justify-content:space-between; gap:8px; color:var(--muted); }
    .fs-sliders output { color:var(--ink); font-variant-numeric:tabular-nums; }
    .fs-sliders input[type=range] { width:100%; height:16px; margin:0; }
    .fs-sliders.disabled { opacity:.35; pointer-events:none; }
    .view-panel { position:absolute; top:12px; right:12px; z-index:2; display:flex; flex-direction:column; gap:0; padding:10px; min-width:136px; border-radius:12px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.45); -webkit-backdrop-filter:blur(14px) saturate(140%); backdrop-filter:blur(14px) saturate(140%); box-shadow:0 8px 28px rgba(0,0,0,.35), inset 0 1px 0 rgba(255,255,255,.08); cursor:default; }
    .view-panel-title { color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; padding:0 2px 2px; }
    .view-panel button { width:100%; padding:6px 10px; font-size:12px; text-align:left; background:rgba(255,255,255,.06); border:1px solid rgba(255,255,255,.12); border-radius:8px; transition:background .15s,border-color .15s; }
    .view-panel button:hover, .view-panel button:focus-visible { background:rgba(123,228,241,.16); border-color:var(--cyan); outline:none; }
    .view-panel button:active { background:rgba(123,228,241,.28); border-color:var(--cyan); }
    .view-panel button[aria-pressed=true] { background:rgba(123,228,241,.22); border-color:var(--cyan); color:var(--ink); }
    body.tutorial-view-intro .view-panel :is(#toggleHalo,#viewFullHalo):not(:focus-visible) { opacity:.4; }
    .view-panel .view-row { display:grid; grid-template-columns:1fr 1fr; gap:6px; }
    .view-panel .view-row button { text-align:center; }
    .view-panel .mode-toggle { display:grid; grid-template-columns:1fr 1fr; padding:2px; margin-top:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .mode-toggle button { padding:4px 8px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); }
    .view-panel .mode-toggle button[aria-checked=true] { background:rgba(123,228,241,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(123,228,241,.5); }
    .view-panel .speed-toggle { display:grid; grid-template-columns:repeat(4,1fr); padding:2px; margin-top:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .speed-toggle button { padding:4px 2px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); transition:background .15s,color .15s; }
    .view-panel .speed-toggle button:hover, .view-panel .speed-toggle button:focus-visible { color:var(--ink); outline:none; }
    .view-panel .speed-toggle button[aria-checked=true] { background:rgba(123,228,241,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(123,228,241,.5); }
    #orbitPlot.pan-mode { cursor:move; }
    #orbitPlot.pan-mode.dragging { cursor:grabbing; }
    .view-panel .view-head { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:0 2px 2px; background:transparent; border:0; border-radius:6px; color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; }
    .view-panel .view-head:hover, .view-panel .view-head:focus-visible { background:transparent; color:var(--ink); }
    .view-panel .chev { display:inline-block; font-size:12px; transition:transform .2s; }
    .view-panel.collapsed { min-width:0; padding:0; border-radius:999px; }
    .view-panel > :not(.view-head) { max-height:80px; margin-top:6px; opacity:1; visibility:visible; overflow:hidden; transition:max-height .2s ease, margin-top .2s ease, opacity .2s ease, visibility 0s; }
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
    .halo-metric { padding:10px; border:1px solid rgba(179,154,232,.14); border-radius:6px; background:rgba(5,9,16,.4); }
    .halo-metric button { display:block; width:100%; padding:0; border:0; background:none; text-align:left; color:#c3cfdb; font-size:11px; }
    .halo-metric button:focus-visible { outline:2px solid var(--cyan); outline-offset:3px; border-radius:2px; }
    .halo-metric .info-icon { color:#b39ae8; margin-left:3px; }
    .halo-metric .metric-value { display:block; margin-top:5px; color:#e6ddfa; font-size:18px; font-weight:600; font-variant-numeric:tabular-nums; }
    .halo-metric dd { display:grid; grid-template-rows:0fr; margin:0; color:var(--muted); font-size:11px; line-height:1.45; opacity:0; visibility:hidden; transition:grid-template-rows .2s ease, opacity .2s ease, visibility 0s linear .2s; }
    .halo-metric dd > span { min-height:0; overflow:hidden; }
    .halo-metric dd > span > span { display:block; padding-top:5px; }
    .halo-metric.expanded dd, .halo-metric.hovered dd { grid-template-rows:1fr; opacity:1; visibility:visible; transition-delay:0s; }
    @media (prefers-reduced-motion:reduce) { .halo-metric dd, .view-panel > :not(.view-head), .view-panel.collapsed > :not(.view-head) { transition:none; } }
    footer { margin:14px 2px 0; color:var(--muted); font-size:11px; display:flex; justify-content:space-between; gap:16px; }
    footer .export-btn { padding:0; border:0; background:none; color:var(--muted); font-size:11px; text-decoration:underline; text-underline-offset:3px; white-space:nowrap; }
    footer .export-btn:hover, footer .export-btn:focus-visible { color:var(--ink); }
    .tutorial-launch { margin-left:auto; background:transparent; color:var(--muted); font-size:12px; }
    .tutorial-launch:hover, .tutorial-launch:focus-visible { color:var(--ink); }
    .tutorial-welcome { width:min(460px,calc(100% - 32px)); max-height:calc(100% - 32px); padding:28px; border:1px solid rgba(123,228,241,.22); border-radius:18px; background:linear-gradient(145deg,#132739,#0b1522 65%); color:var(--ink); box-shadow:0 24px 80px rgba(0,0,0,.6); }
    .tutorial-welcome::backdrop { background:rgba(2,6,12,.8); }
    .tutorial-eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.1em; font-size:10px; font-weight:600; }
    .tutorial-welcome h2 { margin:8px 0; font-size:22px; }
    .tutorial-welcome p { color:#b5c5d5; margin:0 0 18px; }
    .tutorial-welcome-actions { display:flex; flex-wrap:wrap; gap:9px; }
    .tutorial-primary { background:rgba(123,228,241,.16); border-color:var(--cyan); color:#e4fbff; font-weight:600; }
    .tutorial-card { position:fixed; z-index:10000; right:16px; bottom:16px; display:flex; flex-direction:column; width:min(380px,calc(100% - 32px)); max-height:min(65vh,480px); overflow:hidden; padding:18px; border:1px solid rgba(123,228,241,.24); border-radius:16px; background:linear-gradient(155deg,#142838,#0b1623 70%); box-shadow:0 16px 48px rgba(0,0,0,.58); cursor:default; user-select:text; }
    .tutorial-card.on-left { left:16px; right:auto; }
    .tutorial-card[hidden] { display:none; }
    .tutorial-top { display:flex; flex:none; align-items:center; justify-content:space-between; gap:12px; }
    .tutorial-content { min-height:0; overflow-y:auto; overscroll-behavior:contain; }
    .tutorial-skip { border:0; padding:4px; background:transparent; color:var(--muted); font-size:11px; text-decoration:underline; text-underline-offset:3px; }
    .tutorial-skip:hover, .tutorial-skip:focus-visible { color:var(--ink); }
    .tutorial-progress-track { flex:none; height:3px; margin:10px 0 4px; background:rgba(123,228,241,.13); border-radius:4px; }
    .tutorial-progress-fill { height:100%; width:0; border-radius:inherit; background:var(--cyan); transition:width .25s ease; }
    .tutorial-card h2 { margin:10px 0 7px; font-size:18px; line-height:1.25; }
    .tutorial-card p { margin:0 0 10px; color:#c3d0de; font-size:13px; line-height:1.5; }
    .tutorial-card .tutorial-try { margin:2px 0 0; padding:10px 12px; border-left:2px solid var(--cyan); border-radius:0 7px 7px 0; background:rgba(123,228,241,.07); color:#c7f5fa; }
    .tutorial-feedback { flex:none; margin:9px 0 0; color:#a5e3c2; font-size:12px; }
    .tutorial-feedback[hidden] { display:none; }
    .tutorial-card .tutorial-fullscreen { color:#ffc48a; }
    .tutorial-actions { display:flex; flex:none; justify-content:flex-end; gap:8px; margin-top:12px; }
    .tutorial-keys { align-self:center; margin-right:auto; color:var(--muted); font-size:11px; }
    @media (hover:none) and (pointer:coarse) { .tutorial-keys { display:none; } }
    .tutorial-actions button { min-height:36px; }
    .tutorial-actions button:disabled { opacity:.45; cursor:not-allowed; }
    .tutorial-spotlight { position:fixed; z-index:9998; pointer-events:none; border:2px solid #8ceaf5; border-radius:9px; box-shadow:0 0 0 9999px rgba(2,7,14,.30), 0 0 0 4px rgba(123,228,241,.16); transition:opacity .45s ease; }
    .tutorial-spotlight[hidden] { display:none; }
    .tutorial-pointer { position:absolute; left:0; top:0; width:25px; height:30px; filter:drop-shadow(0 2px 3px #000); animation:tutorial-gesture 2.8s ease-in-out infinite; will-change:transform; }
    .tutorial-spotlight.is-interacting .tutorial-pointer { opacity:.15; }
    .tutorial-spotlight.is-complete { opacity:0; }
    .is-complete .tutorial-pointer { animation-play-state:paused; }
    .tutorial-spotlight.is-click .tutorial-pointer { animation-name:tutorial-click; }
    @keyframes tutorial-gesture { 0%,15% { transform:translate(var(--from-x),var(--from-y)); } 65%,85% { transform:translate(var(--to-x),var(--to-y)); } 100% { transform:translate(var(--from-x),var(--from-y)); } }
    @keyframes tutorial-click { 0%,35%,100% { transform:translate(var(--to-x),var(--to-y)) scale(1); } 50%,65% { transform:translate(var(--to-x),var(--to-y)) scale(.8); } }
    @media (prefers-reduced-motion:reduce) { .tutorial-pointer, .tutorial-spotlight.is-click .tutorial-pointer, #visibleMatterOnly.tutorial-reveal { animation:none; } .tutorial-pointer { transform:translate(var(--to-x),var(--to-y)); } .tutorial-progress-fill, .tutorial-spotlight { transition:none; } }
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
  </style>
</head>
<body>
  <header>
    <div><div class="eyebrow">SPARC galaxy dynamics · interactive model</div><h1>Interactive Galaxy Rotation Simulator</h1></div>
    <div class="header-stats" aria-label="Simulation summary">
      <div class="stat"><div class="stat-label">Test stars</div><div class="stat-value"><span id="starCountStat">800</span><small>massless</small></div></div>
      <div class="stat"><div class="stat-label">Timeline</div><div class="stat-value">1 Gyr<small>1000 Myr</small></div></div>
      <div class="ml-col">
        <div class="stat ml" title="Stellar mass-to-light ratios used for the baryonic (visible-matter) trajectories"><div class="stat-label">Mass-to-light (M/L)</div><div class="stat-value ml-pair"><div><b id="mlDisk">0.50</b><span>DISK</span></div><div><b id="mlBulge">0.70</b><span>BULGE</span></div></div></div>
        <label class="fit-toggle" id="mondFitWrap" hidden title="Re-fit a₀ and the disk M/L to the SPARC data to minimise χ² (MOND mode only)"><span class="fit-label">MOND best fit</span><input type="checkbox" id="mondFitToggle" autocomplete="off"><span class="track" aria-hidden="true"></span></label>
      </div>
    </div>
  </header>
  <main>
    <section class="model-bar galaxy-bar" id="galaxyBar" aria-label="Galaxy source">
      <div class="source-toggle" role="radiogroup" aria-label="Galaxy source">
        <label><input type="radio" name="galaxyMode" value="preset" checked autocomplete="off"><span>Preset Galaxies</span></label>
        <label><input type="radio" name="galaxyMode" value="custom" autocomplete="off"><span>Custom Sandbox</span></label>
      </div>
      <div class="preset-picker" id="presetPicker">
        <label class="model-label" id="galaxyLabel" for="galaxySelect">Galaxy</label>
        <select id="galaxySelect" autocomplete="off" aria-describedby="galaxySummary" title="Choose which SPARC galaxy to model"></select>
        <span class="galaxy-summary" id="galaxySummary" role="status" aria-live="polite"></span>
      </div>
    </section>
    <section class="sandbox" id="sandboxPanel" aria-labelledby="sandboxTitle" hidden>
      <div class="sandbox-head"><h2 id="sandboxTitle">Custom galaxy builder</h2><p>Build a galaxy from its visible anatomy. Nothing is blocked: push the sliders to extremes to build a “monster galaxy”.</p></div>
      <div class="sandbox-groups">
        <div class="sandbox-group">
          <h3>Stellar disk</h3>
          <div class="control"><label for="sbDiskMass"><span>Disk mass, M<sub>disk</sub></span><output id="sbDiskMassValue" for="sbDiskMass"></output></label><input id="sbDiskMass" type="range" min="1e9" max="1e11" step="1e8" value="2e10" autocomplete="off"><p class="control-help">Total mass of the stars in the flat disk (exponential profile).</p></div>
          <div class="control"><label for="sbDiskScale"><span>Scale length, R<sub>d</sub></span><output id="sbDiskScaleValue" for="sbDiskScale"></output></label><input id="sbDiskScale" type="range" min="1.5" max="6" step="0.1" value="3.1" autocomplete="off"><p class="control-help">How quickly the disk fades outwards: the brightness falls by a factor of e every R<sub>d</sub>.</p></div>
          <div class="control"><label for="sbTruncation"><span>Truncation radius, R<sub>max</sub></span><output id="sbTruncationValue" for="sbTruncation"></output></label><input id="sbTruncation" type="range" min="15" max="40" step="0.5" value="30" autocomplete="off"><p class="control-help">Outer edge of the galaxy's visible matter. Stars that stray beyond 1.5 × R<sub>max</sub> count as ejected.</p></div>
        </div>
        <div class="sandbox-group">
          <h3>Bulge</h3>
          <div class="control"><label for="sbBulgeMass"><span>Bulge mass, M<sub>bulge</sub></span><output id="sbBulgeMassValue" for="sbBulgeMass"></output></label><input id="sbBulgeMass" type="range" min="0" max="5e10" step="1e8" value="0" autocomplete="off"><p class="control-help">Mass of the round central cluster of stars. Zero gives a pure disk galaxy.</p></div>
          <div class="control" id="sbBulgeRadiusControl"><label for="sbBulgeRadius"><span>Effective radius, r<sub>b</sub></span><output id="sbBulgeRadiusValue" for="sbBulgeRadius"></output></label><input id="sbBulgeRadius" type="range" min="0.2" max="2" step="0.1" value="1" autocomplete="off"><p class="control-help" id="sbBulgeRadiusHelp">Scale radius of the bulge's Hernquist profile. Half of the bulge's mass lies within about 2.4 times this radius.</p></div>
        </div>
        <div class="sandbox-group">
          <h3>Gas</h3>
          <div class="control"><label for="sbGasMass"><span>Gas mass, M<sub>gas</sub></span><output id="sbGasMassValue" for="sbGasMass"></output></label><input id="sbGasMass" type="range" min="5e8" max="2e10" step="1e8" value="1.4e10" autocomplete="off"><p class="control-help">Mass of the cold gas, mostly hydrogen, in a disk three times more extended than the stars.</p></div>
          <dl class="sandbox-summary" aria-label="Sandbox galaxy summary">
            <dt>Total visible mass</dt><dd id="sbTotalMass"></dd>
            <dt>Bulge-to-total (stars)</dt><dd id="sbBulgeFraction"></dd>
            <dt>Gas fraction</dt><dd id="sbGasFraction"></dd>
            <dt>Baseline halo M<sub>200</sub> (SHMR)</dt><dd id="sbHaloBaseline"></dd>
            <dt>Test stars</dt><dd id="sbStarCount"></dd>
            <dt>Ejection radius</dt><dd id="sbEjectRadius"></dd>
          </dl>
        </div>
      </div>
      <div class="stability" id="stabilityMonitor" data-state="stable" aria-label="Halo stability monitor">
        <div class="stability-head"><span class="stability-title">Halo stability monitor</span><span class="stability-badge" id="stabilityBadge" role="status"></span></div>
        <div class="stability-reading"><span class="stability-c" id="stabilityC">c = —</span><span class="stability-detail" id="stabilityDetail"></span></div>
        <div class="c-gauge" aria-hidden="true"><div class="c-band" id="cBand"></div><div class="c-marker" id="cMarker"></div></div>
        <div class="c-ticks" id="cTicks" aria-hidden="true"></div>
        <p class="stability-message" id="stabilityMessage"></p>
      </div>
      <div class="sandbox-footer"><p id="sandboxNote">The anatomy sliders drive the rotation curve, χ², halo and 3D stars live. The halo follows the stellar-to-halo mass relation: fine-tune its mass and scale radius with the halo sliders below. The white points and χ² come from the SPARC galaxy chosen under Compare with.</p><button type="button" id="sandboxReset" title="Return the anatomy sliders to an NGC 3198-like galaxy">Reset anatomy</button></div>
    </section>
    <section class="model-bar" id="modelBar" aria-label="Gravity model">
      <span class="model-label">Gravity model</span>
      <label class="switch" title="Switch between a dark matter halo and modified gravity (MOND)"><input type="checkbox" id="mondToggle" role="switch" autocomplete="off" aria-label="Use MOND theory instead of NFW dark matter halo"><span class="opt active" id="optNfw">NFW Dark Matter Halo</span><span class="track" aria-hidden="true"></span><span class="opt" id="optMond">MOND Theory</span></label>
      <span class="mond-note" id="mondNote" hidden>MOND active: no dark matter (M<sub>200</sub> = 0). Gravity is boosted below a<sub>0</sub> = <span id="mondA0Text">1.2×10⁻¹⁰</span> m/s².</span>
      <button type="button" class="tutorial-launch" id="tutorialLaunch">Tutorial</button>
    </section>
    <section class="controls" aria-label="Dark matter controls">
    <div class="control" id="massControl"><label for="massSlider"><span>Halo mass, M<sub>200</sub></span><output id="massValue"></output></label><input id="massSlider" type="range" min="0.00" max="2.00" step="0.01"><p class="control-help" id="massHelp">Updates the rotation curve and star orbits. Set it to 0 for visible matter only.</p><div class="fit-actions"><button id="applyObservedFit" class="fit-button" type="button">Apply Observed Fit</button><button id="visibleMatterOnly" class="fit-button" type="button" title="Set the dark matter halo mass to zero">Visible Matter Only</button></div></div>
    <div class="control" id="scaleControl"><label for="scaleSlider"><span>Scale radius, r<sub>s</sub></span><output id="scaleValue"></output></label><input id="scaleSlider" type="range" min="5" max="60" step="0.5"><p class="control-help" id="scaleHelp">Controls halo concentration: a smaller radius packs more dark matter near the center.</p></div>
      <div class="control"><label for="timeSlider"><span>Orbit time</span><output id="timeValue"></output></label><div class="timeline-row"><button type="button" id="playPauseBtn" class="play-btn" aria-pressed="false" aria-label="Play simulation" title="Play / pause (Space)"><span class="play-icon" aria-hidden="true">▶</span><span class="play-text">Play</span></button><input id="timeSlider" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"></div></div>
    </section>
    <section class="panel-grid">
    <article class="panel"><div class="panel-heading"><h2>Stellar orbit trajectories</h2></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional __DEFAULT_GALAXY_NAME__ star trajectories"><div id="escape-counter" style="position: absolute; top: 12px; left: 12px; color: rgba(255,255,255,0.75); font-family: monospace; font-size: 11px; background: rgba(0,0,0,0.45); padding: 4px 8px; border-radius: 4px; pointer-events: none; z-index: 100; letter-spacing: 0.5px;"><div>Ejected: <span id="escape-count-val">0</span> / <span id="star-total-val">800</span></div><div id="star-detail" style="font-size: 10px; opacity: 0.7;"></div></div><div class="stability-chip" id="stabilityChip" role="status" hidden></div><div class="view-panel" id="viewPanel" role="group" aria-label="3D view controls"><button type="button" class="view-head" id="viewToggle" aria-expanded="true" title="Collapse or expand the view controls"><span>View</span><span class="chev" aria-hidden="true">▾</span></button><button type="button" id="viewReset" title="Return to the default __DEFAULT_GALAXY_NAME__ view">Reset view</button><button type="button" id="toggleFullscreen" aria-pressed="false" title="Toggle fullscreen view of 3D simulation">Fullscreen</button><button type="button" id="toggleOrbits" aria-pressed="false" title="Show or hide the orbit trajectory lines">Orbit lines</button><button type="button" id="toggleHalo" aria-pressed="false" title="Hide the dark matter halo (the halo is not drawn in MOND mode)">Hide halo</button><button type="button" id="viewFullHalo" title="Pull back to see the full dark matter halo around the tiny disk">View Full Halo</button><div class="view-row"><button type="button" id="viewFaceOn" title="Look straight down on the galactic plane">Face-on</button><button type="button" id="viewEdgeOn" title="Look along the galactic plane">Edge-on</button></div><div class="mode-toggle" role="radiogroup" aria-label="Mouse drag mode"><button type="button" role="radio" aria-checked="true" data-mode="rotate">Rotate</button><button type="button" role="radio" aria-checked="false" data-mode="pan">Pan</button></div><div class="speed-toggle" role="radiogroup" aria-label="Playback speed"><button type="button" role="radio" aria-checked="false" data-speed="0.25" title="0.25x speed">0.25x</button><button type="button" role="radio" aria-checked="false" data-speed="0.5" title="0.5x speed">0.5x</button><button type="button" role="radio" aria-checked="true" data-speed="1" title="1x speed">1x</button><button type="button" role="radio" aria-checked="false" data-speed="2" title="2x speed">2x</button></div></div><div class="fs-sliders" id="fsSliders" role="group" aria-label="Dark matter parameters"><label for="fsMassSlider"><span>Halo mass</span><output id="fsMassValue"></output></label><input id="fsMassSlider" type="range" min="0.00" max="2.00" step="0.01"><label for="fsScaleSlider"><span>Scale radius</span><output id="fsScaleValue"></output></label><input id="fsScaleSlider" type="range" min="5" max="60" step="0.5"></div></div>
      <section class="halo-info" aria-labelledby="haloInfoTitle">
        <h3 id="haloInfoTitle">Your NFW halo</h3>
        <p class="halo-info-intro">The halo's size and shape update with the sliders. 1 kpc is about 3,260 light-years. Hover or select a value to learn more.</p>
        <dl class="halo-metrics" id="haloMetrics">
          <div class="halo-metric"><dt><button type="button" aria-expanded="false" aria-controls="haloRadiusHelp">Halo radius, r<sub>200</sub><span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloRadius">—</span></button></dt><dd id="haloRadiusHelp"><span><span>The model's outer halo radius; not a sharp physical edge.</span></span></dd></div>
          <div class="halo-metric"><dt><button type="button" aria-expanded="false" aria-controls="haloScaleHelp">Scale radius, r<sub>s</sub><span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloScale">—</span></button></dt><dd id="haloScaleHelp"><span><span>Sets how the halo's density changes with distance from the centre.</span></span></dd></div>
          <div class="halo-metric"><dt><button type="button" aria-expanded="false" aria-controls="haloConcentrationHelp">Concentration, c<span class="info-icon" aria-hidden="true">ⓘ</span><span class="metric-value" id="haloConcentration">—</span></button></dt><dd id="haloConcentrationHelp"><span><span>Halo radius divided by scale radius: a value of 10 means the halo radius is 10 times the scale radius. Higher values mean a more centrally concentrated halo.</span></span></dd></div>
        </dl>
        <p class="halo-empty" id="haloEmpty" hidden></p>
      </section>
    </article>
      <article class="panel">
        <div class="panel-heading"><h2>Rotation curve</h2><span id="chiBadge" title="Goodness of fit of baryons + NFW halo to the SPARC points">χ² –</span></div>
        <div id="curvePlot" class="plot" role="img" aria-label="__DEFAULT_GALAXY_NAME__ observed, baryonic, and fitted rotation curves"></div>
        <details class="key">
          <summary>Understanding the colours and χ²</summary>
          <div class="key-content">
            <ul><li><i class="sw dot" style="background:#f4f7fb"></i><b>White points</b> are measured speeds of <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> (SPARC); bars show uncertainty.</li><li><i class="sw" style="background:#5fd0e0"></i><b>Cyan</b> is the speed expected from visible matter only (stars + gas).</li><li><i class="sw dash" style="border-color:#b39ae8"></i><b>Purple</b> is the extra speed supplied by the dark matter halo (NFW mode).</li><li><i class="sw" style="background:#ffae5c"></i><b>Orange</b> is the combined model (or the MOND prediction in MOND mode). The closer it hugs the white points, the better the fit.</li><li><i class="sw bar" style="background:#5fd0e0"></i><i class="sw bar" style="background:#ff8a80"></i><b>Residual bars:</b> teal = model is too slow, coral = model is too fast, no bar = perfect match.</li><li><b>χ² badge:</b> <span class="g">teal</span> good fit, <span class="o">orange</span> fair, <span class="r">red</span> poor.</li><li id="componentKey" hidden><b>Thin dotted lines</b> show your custom galaxy's stellar disk, bulge and gas contributions; they add in quadrature to the cyan baryonic curve.</li></ul>
            <h3>Why isn't χ² exactly 0?</h3>
            <p>χ² adds up how far the model is from each measured point, in units of that point's uncertainty. A score of 0 would mean every point matched perfectly, which real data never allows: measurements carry noise, so even a correct model scatters around them by roughly one error bar per point (χ²ᵣ ≈ 1 is a good fit).</p>
            <p>Models also simplify reality. NFW assumes a smooth spherical halo, MOND uses a characteristic acceleration a<sub>0</sub>, and both rely on an estimated stellar mass for the visible matter. Their small mismatches add to the score, which is why neither theory reaches 0.</p>
          </div>
        </details>
      </article>
    </section>
    <section class="learn" id="theoryGuide" aria-label="Theory guide">
      <details class="nfw"><summary>What is a dark matter halo?</summary>
        <div class="learn-content">
        <p class="tldr">A dark matter halo is an extended region of invisible matter whose gravity adds to that of a galaxy's stars and gas.</p>
        <p>Stars and gas far from the centre of <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> orbit faster than we would expect from the visible matter alone. A halo provides extra gravity that can help explain these speeds. Dark matter does not emit light, but its gravity can reveal its presence.</p>
        <p>This simulator uses an NFW halo: a simple model that is denser near the centre and more diffuse farther out. Adjust the halo mass (M<sub>200</sub>, the mass within r<sub>200</sub>) and scale radius (r<sub>s</sub>) to see how its gravity changes the predicted rotation curve.</p>
        </div>
      </details>
      <details class="mond"><summary>What is MOND?</summary>
        <div class="learn-content">
        <p class="tldr">MOND (Modified Newtonian Dynamics) is an alternative idea: it changes the usual prediction for gravity at very low accelerations, rather than adding a dark matter halo.</p>
        <p>In a galaxy's outskirts, the gravitational acceleration from visible matter is very small. MOND predicts a stronger effect than Newtonian gravity there, helping explain why orbital speeds remain high. The change becomes important around an acceleration called a<sub>0</sub>; the default here is 1.2×10⁻¹⁰ m/s².</p>
        <p>Switch between the two models and compare their curves with the white observed points. A close match is informative, but one galaxy alone cannot settle which explanation is correct.</p>
        </div>
      </details>
    </section>
    <footer><span id="footerModel">Baryons: V<sub>bar</sub>² = V<sub>gas</sub>² + <span id="ftDisk">0.50</span>V<sub>disk</sub>² + <span id="ftBulge">0.70</span>V<sub>bulge</sub>². Halo: NFW, M<sub>200</sub> and r<sub>s</sub> update locally.</span><span id="footerCustom" hidden>Baryons: V<sub>bar</sub>² = V<sub>disk</sub>² + V<sub>bulge</sub>² + V<sub>gas</sub>², from exponential disks and a Hernquist bulge. Halo: NFW, starting from the stellar-to-halo mass relation.</span><span>SPARC values from the project’s local table2.dat.</span><button type="button" id="exportCsv" class="export-btn" title="Download the time-series log and the rotation-curve model as two CSV files">Export raw data (CSV)</button></footer>
  </main>
  <dialog class="tutorial-welcome" id="tutorialWelcome" aria-labelledby="welcomeTitle" aria-describedby="welcomeDescription">
    <div class="tutorial-eyebrow">Explore <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span></div>
    <h2 id="welcomeTitle">Welcome to the simulator</h2>
    <p id="welcomeDescription">New here? Take a short, hands-on tour of the galaxy, its rotation curve and the ideas behind it. You can skip the tour at any time.</p>
    <div class="tutorial-welcome-actions"><button type="button" class="tutorial-primary" id="tutorialStart">Start tutorial</button><button type="button" id="tutorialWelcomeSkip">I've been here before</button></div>
  </dialog>
  <div class="tutorial-spotlight" id="tutorialSpotlight" aria-hidden="true" hidden>
    <svg class="tutorial-pointer" viewBox="0 0 28 34" aria-hidden="true"><path d="M3 2v23l6-5 4 10 5-2-4-10h8z" fill="#eaffff" stroke="#09212c" stroke-width="2" stroke-linejoin="round"/></svg>
  </div>
  <aside class="tutorial-card" id="tutorialCard" aria-labelledby="tutorialTitle" hidden>
    <div class="tutorial-top"><span class="tutorial-eyebrow" id="tutorialProgress"></span><button type="button" class="tutorial-skip" id="tutorialSkip">Skip tutorial</button></div>
    <div class="tutorial-progress-track" aria-hidden="true"><div class="tutorial-progress-fill" id="tutorialProgressFill"></div></div>
    <div class="tutorial-content" id="tutorialContent" tabindex="0"><div aria-live="polite" aria-atomic="true"><h2 id="tutorialTitle" tabindex="-1"></h2><p id="tutorialText"></p><p class="tutorial-try" id="tutorialTry"></p></div></div>
    <div class="tutorial-feedback" id="tutorialFeedback" role="status" aria-live="polite" hidden></div>
    <p class="tutorial-fullscreen" id="tutorialFullscreen" hidden>Exit fullscreen to continue the tutorial.</p>
    <div class="tutorial-actions"><span class="tutorial-keys">← → to navigate</span><button type="button" id="tutorialBack">Back</button><button type="button" class="tutorial-primary" id="tutorialNext">Next</button></div>
  </aside>
  <script>
__GALAXIES_JS__
  </script>
  <script>
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
        if(window.matchMedia('(prefers-reduced-motion: reduce)').matches){
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
    const haloDefaultsOf = galaxy => ({massT: galaxy.default_halo.m200_msun / 1e12, scaleRadius: galaxy.default_halo.scale_radius_kpc});
    let activeGalaxy = GALAXIES[DATA.default_galaxy_id];
    let DEFAULTS = haloDefaultsOf(activeGalaxy);
    let customActive = false; // true while the Custom Sandbox, not a SPARC preset, drives the physics
    let customModel = null; // PhysicsEngine circular-speed functions of the sandbox galaxy
    const G = 4.30091e-6;
    const H0 = 0.07;
    const BG = '#050910';
    const GRID = '#1c2936';
    const textColor = '#bac9d8';
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
            return{buildGalaxy,baselineHalo,anatomyFromRotation};
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
            return{particles:list,bulgeCount,bulgeCut};
        }
        /* ---------- Active galaxy: SPARC rotation curve ----------
           galaxies.js stores the raw SPARC columns. loadGalaxyData() rebuilds every array below for the chosen
           galaxy, so the rest of the page reads radii / vObs / vBary without knowing which galaxy is active. */
        const mlState={disk:DATA.upsilon_disk,bulge:DATA.upsilon_bulge};
        let rotationRows=[],radii=[],vObs=[],obsV=[],obsErr=[],vBary=[];
        const baryonsFor=(ud,ub)=>customActive&&customModel
            ?rotationRows.map(p=>Math.sqrt(customModel.vSq(p.radius_kpc)))
            :rotationRows.map(p=>Math.sqrt(Math.max(0,p.gas_sq+ud*p.disk_kms**2+ub*p.bulge_kms**2)));
        function loadGalaxyData(galaxy){
            const t=galaxy.rotation;
            rotationRows=t.radius_kpc.map((radius,i)=>{
                const gas=t.v_gas_kms[i];
                return{radius_kpc:radius,observed_kms:Math.abs(t.v_obs_kms[i]),error_kms:Math.abs(t.v_obs_err_kms[i]),
                    gas_kms:Math.abs(gas),gas_sq:gas*Math.abs(gas),disk_kms:Math.abs(t.v_disk_kms[i]),bulge_kms:Math.abs(t.v_bulge_kms[i])};
            });
            radii=rotationRows.map(p=>p.radius_kpc);
            vObs=rotationRows.map(p=>p.observed_kms);
            obsV=vObs.slice();
            obsErr=rotationRows.map(p=>p.error_kms);
            vBary=baryonsFor(mlState.disk,mlState.bulge);
        }
        loadGalaxyData(activeGalaxy);
        function observedVelocity(radius){
            const last=radii.length-1;
            if(radius<=radii[0])return vObs[0]*radius/radii[0];
            if(radius>=radii[last])return vObs[last];
            let low=0,high=last;
            while(high-low>1){const m=(low+high)>>1;if(radii[m]<=radius)low=m;else high=m;}
            return vObs[low]+(radius-radii[low])/(radii[high]-radii[low])*(vObs[high]-vObs[low]);
        }
        function nfwVelocityModel(massSolar,scaleRadius) {
            if(massSolar<=0)return ()=>0;
            const rhoCrit=3*H0*H0/(8*Math.PI*G);
            const r200=Math.cbrt(3*massSolar/(4*Math.PI*200*rhoCrit));
            const shape=x=>x<=0?0:Math.log1p(x)-x/(1+x);
            const fc=Math.max(shape(r200/scaleRadius),1e-12);
            return radius=>{
                if(radius<=0)return 0;
                const enclosed=radius>=r200?massSolar:massSolar*shape(radius/scaleRadius)/fc;
                return Math.sqrt(G*enclosed/radius);
            };
        }
        function haloCurve(massT,scaleRadius) {return radii.map(nfwVelocityModel(massT*1e12,scaleRadius));}
        const M_PER_KPC = 3.08567758e19;
        const MOND_A0_DEFAULT = 1.2e-10; // m/s^2
        let MOND_A0 = MOND_A0_DEFAULT;
        let mondMode = false;
        // SI calculation with the standard interpolating function mu(x) = x / sqrt(1 + x^2)
        function mondVelocity(vBarKms, rKpc, a0 = MOND_A0) {
            const rM = rKpc * M_PER_KPC, vMs = vBarKms * 1000;
            if (rM <= 0 || vMs <= 0) return 0;
            const gN = vMs * vMs / rM;
            const g = Math.sqrt((gN * gN + Math.sqrt(gN ** 4 + 4 * gN * gN * a0 * a0)) / 2);
            return Math.sqrt(g * rM) / 1000;
        }
        // Velocity added on top of the baryons (in quadrature): NFW halo, or the MOND boost.
        function extraVelocityModel(massSolar, scaleRadius, vbOf) {
            if (!mondMode) return nfwVelocityModel(massSolar, scaleRadius);
            return radius => {
                if (radius <= 0) return 0;
                const vb = vbOf(radius), vm = mondVelocity(vb, radius);
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

        let sortedEscapeTimes = [];
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
                this._bind();this._apply();
            }
            _bind() {
                const el=this.element;
                this._h={
                    down:e=>this._down(e),move:e=>this._move(e),up:e=>this._up_(e),
                    wheel:e=>this._wheel(e),menu:e=>e.preventDefault(),dbl:()=>this.reset()
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
                if(e.pointerType==='mouse'&&e.button>2)return;
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
        const renderer=new THREE.WebGLRenderer({antialias:true,powerPreference:'high-performance'});
        renderer.setClearColor(0x020206, 1.0); // Deep space dark navy/black.
        orbitPlot.appendChild(renderer.domElement);
        const scene=new THREE.Scene();
        const camera=new THREE.PerspectiveCamera(40,1,1,3000);
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
        // Every star (bulge + disk) is an InstancedMesh sphere; no glowing orbs or halo clutter. The mesh holds the largest
        // possible population and draws only the first starCount instances; syncStarMeshes() refills it.
        const STAR_BASE_RADIUS=0.13;
        const starMaterial = new THREE.MeshBasicMaterial({
            color: 0xffffff, // Neutral base so per-instance colors show unchanged
            transparent: true,
            opacity: 0.6,
            blending: THREE.AdditiveBlending,
            depthWrite: false // Crucial for additive blending so overlapping stars glow properly
        });
        const starMesh=new THREE.InstancedMesh(new THREE.SphereGeometry(STAR_BASE_RADIUS,10,8),starMaterial,MAX_STARS);
        starMesh.setColorAt(0,new THREE.Color(0xffffff)); // allocates the per-instance colours while count is still the capacity
        starMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
        starMesh.frustumCulled=false;starMesh.renderOrder=1;
        scene.add(starMesh);

        function createGlowTexture(){
            const canvas=document.createElement('canvas');
            canvas.width=128;canvas.height=128;
            const ctx=canvas.getContext('2d');
            const grad=ctx.createRadialGradient(64,64,0,64,64,64);
            grad.addColorStop(0.0,'rgba(255, 255, 255, 1.0)');
            grad.addColorStop(0.15,'rgba(255, 235, 200, 0.9)');
            grad.addColorStop(0.4,'rgba(255, 170, 60, 0.4)');
            // Fade to transparent white (not black) to preserve color vibrancy
            grad.addColorStop(1.0,'rgba(255, 255, 255, 0.0)');
            ctx.fillStyle=grad;ctx.fillRect(0,0,128,128);
            return new THREE.CanvasTexture(canvas);
        }
        // Soft glow halo per star, tinted by the star's instance color and following the instanced spheres.
        const glowGeometry=new THREE.BufferGeometry();
        glowGeometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(MAX_STARS*3),3).setUsage(THREE.DynamicDrawUsage));
        glowGeometry.setAttribute('color',new THREE.BufferAttribute(new Float32Array(MAX_STARS*3),3).setUsage(THREE.DynamicDrawUsage));
        const glowPoints=new THREE.Points(glowGeometry,new THREE.PointsMaterial({
            map:createGlowTexture(),size:1.6,sizeAttenuation:true,vertexColors:true,transparent:true,
            opacity:0.55,blending:THREE.AdditiveBlending,depthWrite:false
        }));
        glowPoints.frustumCulled=false;glowPoints.renderOrder=0;
        scene.add(glowPoints);

        // Black hole: compact core shadow with gravitationally lensed photon ring and warm amber/orange glowing arcs.
        const BH_R=GALAXY.blackHoleRadius,BH_EXTENT=3.2;
        const blackHole=new THREE.Mesh(
            new THREE.PlaneGeometry(2*BH_R*BH_EXTENT,2*BH_R*BH_EXTENT),
            new THREE.ShaderMaterial({
                transparent:true,depthWrite:false,blending:THREE.CustomBlending,
                blendEquation:THREE.AddEquation,blendSrc:THREE.OneFactor,blendDst:THREE.OneMinusSrcAlphaFactor,
                vertexShader:`varying vec2 vP;void main(){vP=position.xy;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
                fragmentShader:`
                    varying vec2 vP;
                    void main(){
                        vec2 p=vP/${BH_R.toFixed(3)};
                        float d=length(p);
                        float shadow=1.0-smoothstep(0.93,0.99,d);

                        // Relativistic Doppler beaming across the disk (counter-clockwise spin)
                        float beam=1.0+0.38*(-p.x/max(d,0.25));

                        // 1. Primary circular photon ring hugging the event horizon
                        float ringDist=abs(d-1.035);
                        float photonRing=exp(-pow(ringDist/0.040,2.0))*3.6;
                        float photonShoulder=exp(-pow((d-1.06)/0.10,2.0))*1.2;

                        // 2. Gravitational lensing arcs: upper arc bent over top, lower arc bent under bottom
                        float rUp=length(vec2(p.x/1.80,(p.y+0.25)/1.50));
                        float arcUp=exp(-pow(abs(rUp-1.0)/0.075,2.0))*2.4*smoothstep(-0.20,0.15,p.y);

                        float rDown=length(vec2(p.x/1.65,(p.y-0.20)/1.30));
                        float arcDown=exp(-pow(abs(rDown-1.0)/0.070,2.0))*1.9*smoothstep(0.20,-0.15,p.y);

                        // 3. Equatorial accretion disk extending horizontally
                        float rDisk=length(vec2(p.x/2.15,p.y/0.28));
                        float disk=exp(-pow(abs(rDisk-1.0)/0.095,2.0))*2.0;

                        // Combined emission with Doppler boosting and diffuse amber halo
                        float totalRings=(photonRing+photonShoulder+arcUp+arcDown+disk*0.85)*beam;
                        float glow=exp(-pow((d-1.15)/0.60,2.0))*0.30;
                        float fade=1.0-smoothstep(2.3,3.1,d);

                        // Warm amber/orange palette
                        vec3 cDeepOrange=vec3(0.92,0.24,0.02);
                        vec3 cWarmAmber=vec3(1.0,0.58,0.10);
                        vec3 cGold=vec3(1.0,0.85,0.38);
                        vec3 cHotWhite=vec3(1.0,0.96,0.85);

                        vec3 col=cDeepOrange*glow;
                        col+=cWarmAmber*totalRings;
                        col+=cGold*(pow(photonRing,1.2)*0.45+(arcUp+arcDown)*0.35);
                        col+=cHotWhite*pow(photonRing,2.0)*0.16;

                        // Translucent rendering: allow background stars to subtly show through
                        const float OPACITY = 0.55;
                        col = col * fade * (1.0 - shadow * 0.70) * OPACITY;
                        gl_FragColor = vec4(col, shadow * 0.50);
                    }
                `
            })
        );
        blackHole.position.set(0,0,0);
        blackHole.renderOrder=3;scene.add(blackHole);

        // Optional orbit trajectory lines: full paths of an evenly spaced subset of the stars (mostly disk, a few bulge).
        const ORBIT_LINE_STARS=60,ORBIT_LINE_STEP=5;
        let orbitLineIdx=[];
        const orbitSegs=Math.ceil((frameCount-1)/ORBIT_LINE_STEP);
        const orbitLineGeom=new THREE.BufferGeometry();
        orbitLineGeom.setAttribute('position',new THREE.BufferAttribute(new Float32Array(ORBIT_LINE_STARS*orbitSegs*6),3).setUsage(THREE.DynamicDrawUsage));
        const orbitLines=new THREE.LineSegments(orbitLineGeom,new THREE.LineBasicMaterial({color:0x9ff0ff,transparent:true,opacity:0.6,depthWrite:false,blending:THREE.AdditiveBlending}));
        orbitLines.frustumCulled=false;orbitLines.visible=false;orbitLines.renderOrder=0;scene.add(orbitLines);

        // Refills the star meshes (position, size, colour) and the orbit-line selection from the current population.
        function syncStarMeshes(){
            starMesh.count=starCount;
            glowGeometry.setDrawRange(0,starCount);
            const m=new THREE.Matrix4(),c=new THREE.Color(),glowColors=glowGeometry.attributes.color.array;
            for(let k=0;k<starCount;k++){
                const pt=particles[k],p=pt.position;
                // ~85% warm gold/amber, ~15% crisp ice-blue (seeded for reproducibility)
                if(pt.colorSeed<0.15)c.set(0x88ccff);
                else c.set((pt.colorSeed-0.15)/0.85>0.5?0xffeedd:0xffaa33);
                const sc=pt.kind==='bulge'?0.30+0.08*Math.min(1,pt.radius/bulgeCutRadius):0.32+0.12*(pt.radius/currentAnatomy.rMax);
                m.makeScale(sc,sc,sc).setPosition(p[0],p[1],p[2]);
                starMesh.setMatrixAt(k,m);
                starMesh.setColorAt(k,c);
                glowColors[k*3]=c.r;glowColors[k*3+1]=c.g;glowColors[k*3+2]=c.b;
            }
            starMesh.instanceMatrix.needsUpdate=true;starMesh.instanceColor.needsUpdate=true;
            glowGeometry.attributes.color.needsUpdate=true;
            const diskStars=starCount-bulgeStarCount;
            const bulgeLines=bulgeStarCount>0?Math.min(bulgeStarCount,Math.max(6,Math.round(ORBIT_LINE_STARS*bulgeStarCount/starCount))):0;
            const diskLines=Math.min(diskStars,ORBIT_LINE_STARS-bulgeLines);
            orbitLineIdx=[];
            for(let j=0;j<bulgeLines;j++)orbitLineIdx.push(Math.floor(j*bulgeStarCount/bulgeLines));
            for(let j=0;j<diskLines;j++)orbitLineIdx.push(bulgeStarCount+Math.floor(j*diskStars/diskLines));
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
            haloGroup.visible=show;
            haloToggle.disabled=mondMode;fullHaloBtn.disabled=mondMode||massT<=0;
            haloToggle.setAttribute('aria-disabled',String(mondMode));
            haloMetrics.hidden=mondMode||massT<=0;
            haloEmpty.hidden=!haloMetrics.hidden;
            haloEmpty.textContent=mondMode?'MOND does not use a dark matter halo.':'Set the halo mass above zero to see its size and shape.';
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
            haloToggle.setAttribute('aria-pressed',String(haloHidden));
            haloToggle.textContent=haloHidden?'Show halo':'Hide halo';
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
        // Live telemetry follows the active population: total stars, the disk / bulge split and the ejection radius.
        function updateStarTelemetry(){
            const diskStars=starCount-bulgeStarCount;
            document.getElementById('star-total-val').textContent=starCount;
            document.getElementById('starCountStat').textContent=starCount;
            document.getElementById('star-detail').textContent=`${diskStars} disk${bulgeStarCount>0?` + ${bulgeStarCount} bulge`:''} · ejected beyond ${ejectionRadius.toFixed(0)} kpc`;
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
            const matArray=starMesh.instanceMatrix.array;
            const glowArray=glowGeometry.attributes.position.array;
            const i0=f0*3,i1=f1*3;
            for(let k=0;k<starCount;k++){
                const t=activeTrajectories[k];
                const x=t[i0]+(t[i1]-t[i0])*alpha;
                const y=t[i0+1]+(t[i1+1]-t[i0+1])*alpha;
                const z=t[i0+2]+(t[i1+2]-t[i0+2])*alpha;
                const off=k*16;
                matArray[off+12]=x;matArray[off+13]=y;matArray[off+14]=z;
                glowArray[k*3]=x;glowArray[k*3+1]=y;glowArray[k*3+2]=z;
            }
            starMesh.instanceMatrix.needsUpdate=true;
            glowGeometry.attributes.position.needsUpdate=true;
        }

        // Render loop: continuous delta-time playback spanning 0 to 1000 Myr (1 Gyr)
        let needsRender=true,lastTime=performance.now(),simTime=0.0,playing=false,lastReportedFrame=-1,playbackMultiplier=1.0;
        window.simTime=simTime;
        const PLAYBACK_SPEED_MYR_PER_SEC=25.0; // half of the nominal 50 Myr/s; UI labels stay 0.25x-2x
        const controls=new SmoothOrbitControls(camera,orbitPlot,{initial:GALAXY_VIEW,onChange:()=>{needsRender=true;}});
        function requestRender(){needsRender=true;}
        // Panel clicks must not reach the controls' drag / wheel / double-click handlers.
        const viewPanel=document.getElementById('viewPanel');
        ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu'].forEach(type=>viewPanel.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
        const viewToggle=document.getElementById('viewToggle');
        viewPanel.classList.add('collapsed');
        viewToggle.setAttribute('aria-expanded','false');
        viewToggle.addEventListener('click',()=>{
            const collapsed=viewPanel.classList.toggle('collapsed');
            viewToggle.setAttribute('aria-expanded',String(!collapsed));
        });
        document.querySelectorAll('.halo-metric').forEach(metric=>{
            const button=metric.querySelector('button');
            metric.addEventListener('pointerenter',event=>{
                if(event.pointerType!=='touch')metric.classList.add('hovered');
            });
            metric.addEventListener('pointerleave',event=>{
                if(event.pointerType!=='touch')metric.classList.remove('hovered');
            });
            button.addEventListener('click',()=>{
                const expanded=button.getAttribute('aria-expanded')!=='true';
                button.setAttribute('aria-expanded',String(expanded));
                metric.classList.toggle('expanded',expanded);
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
        window.toggleFullscreen = toggleFullscreen;
        const fsBtn = document.getElementById('toggleFullscreen') || document.getElementById('viewFullscreen');
        if (fsBtn) {
            fsBtn.addEventListener('click', toggleFullscreen);
        }
        const viewFsBtn = document.getElementById('viewFullscreen');
        if (viewFsBtn && viewFsBtn !== fsBtn) {
            viewFsBtn.addEventListener('click', toggleFullscreen);
        }
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
            const f=v=>v.toFixed(3);
            let csv="Radius_kpc,V_obs,V_obs_err,V_gas,V_disk,V_bar,V_nfw_halo,V_total_model,Residual_V_obs_minus_V_bar\n";
            for(let r=Math.max(1,Math.ceil(radii[0]-1e-9));r<=radii[radii.length-1]+1e-9;r++){
                const vo=interpolateAt(vObs,r),eo=interpolateAt(err,r),vg=interpolateAt(gas,r),vd=interpolateAt(disk,r);
                const vb=Math.sqrt(Math.max(0,interpolateAt(gasSq,r)+vd*vd)),vh=mondMode?Math.sqrt(Math.max(mondVelocity(vb,r)**2-vb*vb,0)):halo(r);
                csv+=[r,f(vo),f(eo),f(vg),f(vd),f(vb),f(vh),f(Math.hypot(vb,vh)),f(vo-vb)].join(',')+"\n";
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
            orbitToggle.setAttribute('aria-pressed',String(orbitLines.visible));
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
        function resizeOrbit(){
            if (!nativeFsElement()) {
                containerOriginalWidth = orbitPlot.clientWidth;
                containerOriginalHeight = orbitPlot.clientHeight;
            }
            const w=Math.max(orbitPlot.clientWidth,1),h=Math.max(orbitPlot.clientHeight,1),dpr=Math.min(window.devicePixelRatio||1,2);
            renderer.setPixelRatio(dpr);renderer.setSize(w,h,false);
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
        function renderLoop(now){
            const dt=Math.min((now-lastTime)/1000,0.1);
            lastTime=now;
            if(controls.update(dt))needsRender=true;
            if(playing){
                simTime+=dt*PLAYBACK_SPEED_MYR_PER_SEC*playbackMultiplier;
                window.simTime=simTime;
                if(simTime>=SIM_TOTAL_TIME_MYR){
                    simTime=SIM_TOTAL_TIME_MYR;
                    window.simTime=simTime;
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
            if(needsRender){needsRender=false;haloUniforms.uScale.value=renderer.domElement.height/(2*Math.tan(camera.fov*Math.PI/360));
                // Fade the boundary grid as the camera closes in on the visible galaxy.
                const haloFade=0.3+0.7*Math.min(1,Math.max(0,(controls.view.radius-30)/120));
                haloShell.material.uniforms.uFade.value=haloFade;haloGrid.material.opacity=0.13*haloFade;blackHole.quaternion.copy(camera.quaternion);renderer.render(scene,camera);}
            requestAnimationFrame(renderLoop);
        }
        resizeOrbit();requestAnimationFrame(renderLoop);

    const initialHalo=haloCurve(DEFAULTS.massT,DEFAULTS.scaleRadius);
    const totalCurve=vBary.map((v,i)=>Math.hypot(v,initialHalo[i]));
    function updateChi(total){const chi=obsV.reduce((s,v,i)=>s+((v-total[i])/obsErr[i])**2,0),dof=Math.max(1,obsV.length-2),red=chi/dof,b=document.getElementById('chiBadge');b.textContent='χ² '+chi.toFixed(0)+' · χ²ᵣ '+red.toFixed(2);b.className=red<2?'good':red<10?'ok':'poor';}
    const residualsOf=total=>obsV.map((v,i)=>v-total[i]);
    const AXIS={gridcolor:'rgba(120,150,180,0.10)',zerolinecolor:'rgba(120,150,180,0.22)',linecolor:'rgba(120,150,180,0.25)',tickcolor:'rgba(120,150,180,0.25)',ticks:'outside',ticklen:4,automargin:true,tickfont:{size:11,color:'#9fb0c2'}};
    // Residual bars are about half the mean spacing of the SPARC radii, so dense and sparse curves both read clearly.
    const residualBarWidth=()=>Math.round(53*(radii[radii.length-1]-radii[0])/(radii.length-1))/100;
    // The custom galaxy's stellar disk, bulge and gas speeds are drawn as thin dotted lines (hidden for presets).
    const COMPONENT_STYLE=[['Stellar disk','#ffd166'],['Bulge','#ff8fab'],['Gas','#7bd88f']];
    let componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
    const componentVisible=i=>customActive&&(i!==1||currentAnatomy.bulgeMass>0);
    const buildCurveData=(halo,total)=>[
      {type:'scatter',mode:'markers',name:'Observed (SPARC)',x:radii,y:obsV,error_y:{type:'data',array:obsErr,visible:true,color:'rgba(225,235,245,0.35)',thickness:1,width:2},marker:{size:7,color:'#f4f7fb',line:{color:'rgba(5,9,16,0.9)',width:1.2}},hovertemplate:'V_obs  <b>%{y:.1f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Baryonic',x:radii,y:vBary.slice(),line:{color:'#5fd0e0',width:2.5,shape:'spline',smoothing:0.6},hovertemplate:'V_baryonic  <b>%{y:.1f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Dark matter halo',x:radii,y:halo,fill:'tozeroy',fillcolor:'rgba(179,154,232,0.10)',line:{color:'#b39ae8',width:2,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:'V_halo  <b>%{y:.1f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Total model',x:radii,y:total,line:{color:'#ffae5c',width:3.5,shape:'spline',smoothing:0.6},hovertemplate:'V_model  <b>%{y:.1f}</b> km/s<extra></extra>'},
      {type:'bar',name:'Residual',x:radii,y:residualsOf(total),yaxis:'y2',showlegend:false,width:residualBarWidth(),marker:{color:residualsOf(total).map(r=>r>=0?'rgba(95,208,224,0.75)':'rgba(255,138,128,0.75)')},hovertemplate:'Residual  <b>%{y:+.1f}</b> km/s<extra></extra>'},
      ...COMPONENT_STYLE.map(([name,color],i)=>({type:'scatter',mode:'lines',name,x:radii,y:componentCurves[i].slice(),visible:componentVisible(i),line:{color,width:1.6,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:name+'  <b>%{y:.1f}</b> km/s<extra></extra>'}))
    ];
    const curveLayout={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',margin:{l:70,r:22,t:46,b:56},font:{color:textColor,size:12,family:'Inter,Segoe UI,system-ui,sans-serif'},
      xaxis:{...AXIS,fixedrange:true,title:{text:'Radius (kpc)',standoff:10},anchor:'y2',showgrid:false,zeroline:false,spikemode:'across',spikecolor:'rgba(200,215,230,0.35)',spikethickness:1,spikedash:'dot',spikesnap:'cursor'},
      yaxis:{...AXIS,title:{text:'Orbital velocity (km/s)',standoff:10},range:[0,190],fixedrange:true,domain:[0.32,1],zeroline:false},
      yaxis2:{...AXIS,title:{text:'Residual',standoff:6},domain:[0,0.24],anchor:'x',zeroline:true,zerolinewidth:1.5,nticks:4,fixedrange:true},
      bargap:0.2,legend:{orientation:'h',yanchor:'bottom',y:1.03,xanchor:'left',x:0,font:{size:11,color:'#c3cfdb'},itemsizing:'constant'},
      hovermode:'x unified',hoverlabel:{bgcolor:'rgba(12,17,24,0.96)',bordercolor:'rgba(120,150,180,0.35)',font:{color:'#f1f5f9',size:12}},dragmode:false};
    const curveConfig={responsive:true,displaylogo:false,displayModeBar:false,scrollZoom:false,doubleClick:false};
    Plotly.newPlot('curvePlot',buildCurveData(initialHalo,totalCurve),curveLayout,curveConfig);
    updateChi(totalCurve);

    const massSlider=document.getElementById('massSlider'),scaleSlider=document.getElementById('scaleSlider');
    let sliderUpdateTimer=null,sliderDragging=false;
    const fsSliders=document.getElementById('fsSliders'),fsMass=document.getElementById('fsMassSlider'),fsScale=document.getElementById('fsScaleSlider');
    // The mass slider grows for massive halos so a fitted value is never pinned to the end of its track. In the Custom
    // Sandbox it spans 0 to HALO_FINE_RANGE x the stellar-to-halo baseline, so it fine-tunes that baseline: haloFine is
    // M200 as a multiple of it, kept when the anatomy moves the baseline. rsManual records that the user took over r_s.
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
      fsSliders.classList.toggle('disabled',mondMode);
      fsMass.disabled=fsScale.disabled=mondMode;
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
      const total=mondMode?vBary.map((v,i)=>mondVelocity(v,radii[i])):vBary.map((v,i)=>Math.hypot(v,halo[i]));
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
    const BULGE_RADIUS_HELP=bulgeRadiusHelp.textContent;
    function updateSandbox(){
      SANDBOX_SLIDERS.forEach(({key,input,output,format})=>{
        sandbox[key]=Number(input.value);
        output.textContent=format(sandbox[key]);
      });
      const noBulge=sandbox.bulgeMass<=0;
      bulgeRadiusControl.classList.toggle('disabled',noBulge);
      bulgeRadiusControl.setAttribute('aria-disabled',String(noBulge));
      document.getElementById('sbBulgeRadius').disabled=noBulge;
      bulgeRadiusHelp.textContent=noBulge?'No bulge: this is a pure disk galaxy.':BULGE_RADIUS_HELP;
      const stellar=sandbox.diskMass+sandbox.bulgeMass,visible=stellar+sandbox.gasMass;
      document.getElementById('sbTotalMass').textContent=formatMass(visible);
      document.getElementById('sbBulgeFraction').textContent=formatPercent(sandbox.bulgeMass/stellar);
      document.getElementById('sbGasFraction').textContent=formatPercent(sandbox.gasMass/visible);
      document.dispatchEvent(new CustomEvent('sandboxchange',{detail:{...sandbox}}));
    }
    SANDBOX_SLIDERS.forEach(({input})=>input.addEventListener('input',updateSandbox));
    document.getElementById('sandboxReset').addEventListener('click',()=>{
      SANDBOX_SLIDERS.forEach(({input})=>{input.value=input.defaultValue;});
      updateSandbox();
    });
    updateSandbox();

    // Halo stability monitor: concentration c = r200 / rs of the current NFW halo.
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
    const UNSTABLE_BADGE='⚠️ Unphysical / Unstable Galaxy Profile';
    function stabilityReport(){
      if(mondMode)return{state:'none',badge:'No halo to assess',message:'MOND has no dark matter halo, so there is no concentration to check. Switch back to the NFW halo to see it.'};
      const massT=Number(massSlider.value),rs=Number(scaleSlider.value);
      if(massT<=0)return{state:'none',badge:'No halo to assess',message:'With the halo mass at zero there is no halo to check. Raise it to see its concentration.'};
      const r200=Math.cbrt(3*massT*1e12/(4*Math.PI*200*(3*H0*H0/(8*Math.PI*G))));
      const c=r200/rs,range=`${STABLE_C_MIN} and ${STABLE_C_MAX}`;
      const report={c,r200,rs};
      if(c<STABLE_C_MIN)return{...report,state:'diffuse',badge:UNSTABLE_BADGE,
        message:`c is below ${STABLE_C_MIN}: this halo is far more spread out than real NFW halos, whose concentrations stay between ${range}. You can still run the simulation to see how the galaxy behaves.`};
      if(c>STABLE_C_MAX)return{...report,state:'dense',badge:UNSTABLE_BADGE,
        message:`c is above ${STABLE_C_MAX}: this halo packs its mass into an extremely dense centre, beyond real NFW halos, whose concentrations stay between ${range}. You can still run the simulation to see how the galaxy behaves.`};
      return{...report,state:'stable',badge:'✓ Stable halo profile',message:`c is inside the ${STABLE_C_MIN}–${STABLE_C_MAX} range seen for real NFW halos. The monitor follows the halo mass and scale radius sliders.`};
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
      if(changed){if(custom)activateCustom();else loadPreset(activeGalaxy);}
      updateStability();
      document.dispatchEvent(new CustomEvent('galaxymodechange',{detail:{mode}}));
    }
    galaxyModeRadios.forEach(radio=>radio.addEventListener('change',()=>{if(radio.checked)setGalaxyMode(radio.value);}));
    setGalaxyMode(galaxyModeRadios.find(radio=>radio.checked).value);

    function updateHalo(){
      syncFsSliders();
      updateStability();
      clearTimeout(sliderUpdateTimer);sliderUpdateTimer=null;
      const massT=Number(massSlider.value),rs=Number(scaleSlider.value);
      const{halo,total}=modelCurves(massT,rs);
      document.getElementById('massValue').textContent=haloMassText(massT);
      document.getElementById('scaleValue').textContent=`${rs.toFixed(1)} kpc`;
      const res=residualsOf(total);updateChi(total);
      document.getElementById('chiBadge').title=mondMode?'Goodness of fit of the MOND prediction to the SPARC points':'Goodness of fit of baryons + NFW halo to the SPARC points';
      const yMax=Math.max(...obsV.map((v,i)=>v+obsErr[i]),...vBary,...halo,...total)*1.1;
      // One Plotly call redraws the model lines, the three custom-galaxy component lines and the residual bars (trace 4).
      const residualColors=res.map(r=>r>=0?'rgba(95,208,224,0.75)':'rgba(255,138,128,0.75)');
      Plotly.update('curvePlot',{y:[vBary.slice(),halo,total,res,...componentCurves.map(c=>c.slice())],name:['Baryonic',mondMode?'Dark matter halo (none)':'Dark matter halo',mondMode?'MOND model':'Total model','Residual',...COMPONENT_STYLE.map(([name])=>name)],visible:[true,!mondMode,true,true,...COMPONENT_STYLE.map((_,i)=>componentVisible(i))],'marker.color':[undefined,undefined,undefined,residualColors,undefined,undefined,undefined]},{'yaxis.range':[0,yMax]},[1,2,3,4,5,6,7]);
      updateHaloVisual(massT,rs);
      scheduleOrbitUpdate();
    }
    const mondToggle=document.getElementById('mondToggle'),mondFitToggle=document.getElementById('mondFitToggle');
    // Grid search over a0 and the disk M/L (bulge M/L stays fixed) for the lowest chi-squared.
    let mondBest=null;
    function findMondBestFit(){
      if(mondBest)return mondBest;
      let best={chi:Infinity,disk:mlState.disk,a0:MOND_A0_DEFAULT};
      for(let ud=0.05;ud<=1.5001;ud+=0.01){
        const vb=baryonsFor(ud,DATA.upsilon_bulge);
        for(let a0=0.3e-10;a0<=4e-10;a0+=0.01e-10){
          let chi=0;
          for(let i=0;i<vb.length;i++){const d=(obsV[i]-mondVelocity(vb[i],radii[i],a0))/obsErr[i];chi+=d*d;}
          if(chi<best.chi)best={chi,disk:ud,a0};
        }
      }
      return mondBest=best;
    }
    function applyBaryonState(){
      // A custom galaxy is built from masses, so mass-to-light ratios and the MOND best fit (which tunes them) do not apply.
      const fit=mondMode&&!customActive&&mondFitToggle.checked,best=fit?findMondBestFit():null;
      mlState.disk=fit?best.disk:DATA.upsilon_disk;
      MOND_A0=fit?best.a0:MOND_A0_DEFAULT;
      baryonsFor(mlState.disk,mlState.bulge).forEach((v,i)=>{vBary[i]=v;});
      document.getElementById('mondFitWrap').hidden=!mondMode||customActive;
      const ratioText=ratio=>customActive?'—':ratio.toFixed(2);
      document.getElementById('mlDisk').textContent=ratioText(mlState.disk);
      document.getElementById('mlBulge').textContent=ratioText(mlState.bulge);
      document.getElementById('ftDisk').textContent=mlState.disk.toFixed(2);
      document.getElementById('ftBulge').textContent=mlState.bulge.toFixed(2);
      const e=Math.floor(Math.log10(MOND_A0));
      document.getElementById('mondA0Text').textContent=`${(MOND_A0/10**e).toFixed(2)}×10${String(e).replace('-','⁻').replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`;
    }
    mondFitToggle.addEventListener('change',()=>{applyBaryonState();updateHalo();});
    mondToggle.addEventListener('change',()=>{
      mondMode=mondToggle.checked;
      if(!mondMode)mondFitToggle.checked=false;
      document.getElementById('optNfw').classList.toggle('active',!mondMode);
      document.getElementById('optMond').classList.toggle('active',mondMode);
      document.getElementById('mondNote').hidden=!mondMode;
      applyBaryonState();
      ['massControl','scaleControl'].forEach(id=>{
        const c=document.getElementById(id);
        c.classList.toggle('disabled',mondMode);
        c.setAttribute('aria-disabled',String(mondMode));
        c.querySelectorAll('input,button').forEach(el=>{el.disabled=mondMode;});
      });
      updateHalo();
    });
    function queueSliderUpdate(){
      updateStability();
      document.getElementById('massValue').textContent=haloMassText(Number(massSlider.value));
      document.getElementById('scaleValue').textContent=`${Number(scaleSlider.value).toFixed(1)} kpc`;
      if(sliderUpdateTimer===null)sliderUpdateTimer=setTimeout(updateHalo,100);
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
    window.addEventListener('pointerup',finishSliderDrag);
    window.addEventListener('pointercancel',finishSliderDrag);
    updateHalo();
    document.getElementById('applyObservedFit').addEventListener('click',()=>{
        massSlider.value=DEFAULTS.massT;
        scaleSlider.value=DEFAULTS.scaleRadius;
        haloFine=1;rsManual=false;
        updateHalo();
    });
    document.getElementById('visibleMatterOnly').addEventListener('click',()=>{
        massSlider.value=0;
        haloFine=0;
        updateHalo();
    });

    /* ---------- Galaxy presets and the Custom Sandbox physics ----------
       The dropdown swaps the active SPARC galaxy: its baryonic baseline, observed halo fit, slider ranges,
       rotation-curve plot, chi-squared and orbits are all rebuilt from that galaxy's galaxies.js entry.
       In the Custom Sandbox the same dropdown ("Compare with") only chooses which galaxy's measured speeds the
       custom galaxy is compared with; the anatomy sliders then drive every output through PhysicsEngine. */
    const galaxySelect=document.getElementById('galaxySelect'),galaxySummary=document.getElementById('galaxySummary');
    Object.values(GALAXIES).forEach(galaxy=>galaxySelect.add(new Option(galaxy.name,galaxy.id)));
    function applyGalaxyLabels(galaxy){
      document.querySelectorAll('[data-galaxy-name]').forEach(el=>{el.textContent=galaxy.name;});
      galaxySummary.textContent=`${galaxy.summary} · ${galaxy.distance_mpc} Mpc · ${galaxy.rotation.radius_kpc.length} SPARC points${galaxyMode==='custom'?' · data points and χ² for your galaxy':''}`;
      orbitPlot.setAttribute('aria-label',`Interactive three-dimensional ${galaxyMode==='custom'?'custom galaxy':galaxy.name} star trajectories`);
      document.getElementById('curvePlot').setAttribute('aria-label',`${galaxyMode==='custom'?'Custom galaxy compared with '+galaxy.name:galaxy.name} observed, baryonic, and fitted rotation curves`);
      document.getElementById('viewReset').title=`Return to the default ${galaxyMode==='custom'?'custom galaxy':galaxy.name} view`;
    }
    galaxySelect.value=activeGalaxy.id;
    applyGalaxyLabels(activeGalaxy);
    // Wording that differs between a preset galaxy and the sandbox.
    function applyModeLabels(){
      const custom=customActive,fitButton=document.getElementById('applyObservedFit');
      document.getElementById('galaxyLabel').textContent=custom?'Compare with':'Galaxy';
      galaxySelect.title=custom?'Choose the SPARC galaxy whose measured speeds your custom galaxy is compared with':'Choose which SPARC galaxy to model';
      document.getElementById('massHelp').textContent=custom
        ?'Fine-tunes the halo around the stellar-to-halo mass relation baseline for your galaxy: the slider spans 0 to 5 times that baseline. Set it to 0 for visible matter only.'
        :'Updates the rotation curve and star orbits. Set it to 0 for visible matter only.';
      document.getElementById('scaleHelp').textContent=custom
        ?'Controls halo concentration. It follows the mass–concentration relation for your galaxy until you move it.'
        :'Controls halo concentration: a smaller radius packs more dark matter near the center.';
      fitButton.textContent=custom?'Apply SHMR Baseline':'Apply Observed Fit';
      if(custom)fitButton.title='Return the halo to the stellar-to-halo mass relation baseline';else fitButton.removeAttribute('title');
      document.getElementById('footerModel').hidden=custom;
      document.getElementById('footerCustom').hidden=!custom;
      document.getElementById('componentKey').hidden=!custom;
      applyGalaxyLabels(activeGalaxy);
    }
    // Makes a SPARC preset the active galaxy: the page then follows its data, with stars on its observed orbits.
    function loadPreset(galaxy){
      customActive=false;
      activeGalaxy=galaxy;
      DEFAULTS=haloDefaultsOf(galaxy);
      mondBest=null; // a cached MOND fit belongs to the previous galaxy
      loadGalaxyData(galaxy);
      applyBaryonState();
      configureHaloSliders();
      applyModeLabels();
      componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
      referenceSpeedSquared=presetReferenceSpeedSquared;
      const anatomy=presetAnatomy(galaxy);
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),PRESET_TEMPLATE_MASS));
      syncStarMeshes();updateStarTelemetry();
      // Restart the orbit clock and integrate the new population in the new galaxy's potential.
      simTime=0;window.simTime=0;syncTimelineUI(0);
      refreshOrbits();
      const{halo,total}=modelCurves(DEFAULTS.massT,DEFAULTS.scaleRadius);
      Plotly.react('curvePlot',buildCurveData(halo,total),curveLayout,curveConfig);
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
    function rebuildCustomGalaxy({replot=false}={}){
      const anatomy={diskMass:sandbox.diskMass,diskScale:sandbox.diskScale,rMax:sandbox.truncation,bulgeMass:sandbox.bulgeMass,bulgeRadius:sandbox.bulgeRadius,gasMass:sandbox.gasMass};
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
      // Halo: the stellar-to-halo mass relation sets the baseline; haloFine and rsManual keep the user's own tuning.
      const baseline=PhysicsEngine.baselineHalo(anatomy.diskMass+anatomy.bulgeMass);
      DEFAULTS={massT:baseline.m200/1e12,scaleRadius:snapScaleRadius(baseline.scaleRadius)};
      configureHaloSliders();
      // Stars start on the orbits of the baryons plus baseline halo, then feel whatever halo the sliders leave.
      const baselineHalo=nfwVelocityModel(baseline.m200,DEFAULTS.scaleRadius);
      referenceSpeedSquared=r=>{const vb=baryonicVelocity(r),vh=baselineHalo(r);return vb*vb+vh*vh;};
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),sandboxTemplateMass()));
      syncStarMeshes();updateStarTelemetry();updateCustomSummary(baseline);
      if(replot){
        const{halo,total}=modelCurves(Number(massSlider.value),Number(scaleSlider.value));
        Plotly.react('curvePlot',buildCurveData(halo,total),curveLayout,curveConfig);
      }
      refreshOrbits();
      updateHalo();clearTimeout(orbitModelTimer);
    }
    function activateCustom(){
      customActive=true;haloFine=1;rsManual=false;mondBest=null;
      applyModeLabels();
      simTime=0;window.simTime=0;syncTimelineUI(0);
      rebuildCustomGalaxy({replot:true});
      document.dispatchEvent(new CustomEvent('galaxychange',{detail:{galaxy:activeGalaxy}}));
    }
    let customUpdateTimer=null;
    document.addEventListener('sandboxchange',()=>{
      if(!customActive||customUpdateTimer!==null)return;
      customUpdateTimer=setTimeout(()=>{customUpdateTimer=null;rebuildCustomGalaxy();},40);
    });
    document.getElementById('sandboxReset').addEventListener('click',()=>{haloFine=1;rsManual=false;});
    function selectGalaxy(id){
      const galaxy=GALAXIES[id];
      if(!galaxy||galaxy===activeGalaxy)return;
      if(!customActive){loadPreset(galaxy);return;}
      // Custom Sandbox: only the measured data change; the anatomy and halo stay as built.
      activeGalaxy=galaxy;
      mondBest=null;
      loadGalaxyData(galaxy);
      applyGalaxyLabels(galaxy);
      rebuildCustomGalaxy({replot:true});
      document.dispatchEvent(new CustomEvent('galaxychange',{detail:{galaxy}}));
    }
    galaxySelect.addEventListener('change',()=>selectGalaxy(galaxySelect.value));

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
    }
    document.getElementById('playPauseBtn').addEventListener('click',()=>togglePlayback());

    function togglePlayback(){
        if(playing){
            setPlaying(false);
        }else{
            if(simTime>=SIM_TOTAL_TIME_MYR){
                simTime=0.0;
                window.simTime=simTime;
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
    window.addEventListener('keydown',e=>{
        if(e.code==='Escape'&&orbitPlot.classList.contains('fullscreen-active')){setCssFullscreen(false);return;}
        if(e.code==='Space'&&!(e.target&&e.target.closest&&e.target.closest('button,input,select,textarea,a,summary,dialog,#tutorialCard,[contenteditable="true"]'))){
            e.preventDefault();
            if(!e.repeat)togglePlayback();
        }
    });

    function syncTimelineUI(timeMyr){
        const val=Math.round(Math.max(0,Math.min(SIM_TOTAL_TIME_MYR,timeMyr)));
        timeSlider.value=val;
        document.getElementById('timeValue').textContent=`${val} Myr / 1000 Myr`;
    }

    function onSliderInput(val){
        simTime=Math.max(0,Math.min(SIM_TOTAL_TIME_MYR,Number(val)));
        window.simTime=simTime;
        updateStarTransforms(simTime);
        updateEscapeCounter(simTime);
        syncTimelineUI(simTime);
        requestRender();
    }
    timeSlider.addEventListener('input',()=>onSliderInput(timeSlider.value));
    updateStarTransforms(0);updateEscapeCounter(0);syncTimelineUI(0);

    const tutorialSteps=[
      {target:'curvePlot',gesture:'scan',title:'Read the rotation curve',
       text:'The horizontal axis is distance from the galaxy’s centre; the vertical axis is orbital speed. White points are measured speeds, with bars showing uncertainty. Cyan predicts the speed from visible matter alone.',
       action:'Look towards the right of the graph: how does the cyan curve compare with the white points?'},
      {target:'massSlider',gesture:'zero',title:'Try visible matter alone',
       text:'The orange line is the current model prediction. In NFW mode it combines visible matter with a dark matter halo; the purple line shows the halo contribution.',
       action:'Drag halo mass all the way left to 0. Compare the orange line with the observed points.'},
      {target:'massSlider',gesture:'increase',title:'Add a dark matter halo',
       text:'A halo adds gravity without adding visible light. The mass slider changes the halo mass within its model radius, r₂₀₀.',
       action:'Increase halo mass and watch the orange curve. Can you bring it closer to the white points?'},
      {target:'scaleSlider',gesture:'adjust',title:'Change the halo’s shape',
       text:'The scale radius affects how the halo’s density changes with distance. At the same mass, a smaller scale radius makes the halo more concentrated.',
       action:'Move the scale radius slider and compare the curve. The halo cards below the 3D view show the resulting size and concentration.'},
      {target:'applyObservedFit',gesture:'click',title:'How close is the fit?',
       text:'Residual bars show observed speed minus model speed. χ² combines the differences while accounting for measurement uncertainties; a lower value means a closer fit to these points.',
       action:'Try different settings, or press Apply Observed Fit to return to the starting fit. A close fit alone does not prove a theory.'},
      {target:'modelBar',gesture:'switch',title:'Compare another explanation',
       text:'MOND predicts a stronger gravitational effect at very low accelerations, without a dark matter halo. The orange curve becomes its prediction; halo sliders are unavailable in this mode.',
       action:'Switch to MOND and compare the orange curve with the observations. Switch back to NFW whenever you like.'},
      {target:'playPauseBtn',gesture:'click',title:'Watch a galaxy evolve',
       text:'The 3D view shows simulated, massless test stars, not observations of individual stars. The models calculate their motion over one billion years.',
       action:'Press Play or move the time slider. Try changing the model and see how the orbits respond.'},
      {target:'orbitPlot',gesture:'drag',title:'Interpret the animation',
       text:'“Ejected” counts stars that have travelled beyond 1.5 times the galaxy’s outer radius; it does not necessarily mean they are gravitationally unbound. The halo is a visual representation of the model.',
       action:'Drag to rotate the view. Use the mouse wheel or a pinch gesture to zoom.'},
      {target:'viewToggle',gesture:'click',title:'Explore the View controls',
       text:'View opens the camera and display controls. Face-on and Edge-on change your angle; Reset view returns to the start. With a mouse, Rotate or Pan changes how dragging moves the view. Hide halo and Orbit lines change what is drawn; you can also adjust speed or enter fullscreen.',
       action:'Press View to open or close the controls. Try a camera angle or turn on orbit lines. We’ll explore the dimmed halo controls next.'},
      {target:'viewFullHalo',gesture:'click',title:'See the whole halo',
       text:'In the NFW model, a vast dark matter halo surrounds the much smaller visible galaxy. The blue sphere is a model illustration, not a photograph; its outer radius r₂₀₀ is a conventional boundary, not a sharp edge.',
       action:'Press View Full Halo to zoom out. The galaxy becomes tiny at this scale. Use Reset view to return.'},
      {target:'theoryGuide',gesture:'guide',title:'Keep exploring',
       text:'Both models simplify a real galaxy. The result also depends on measurement uncertainties and the estimated stellar mass. One galaxy alone cannot settle which explanation is correct.',
       action:'Open the theory guide below for more. You can restart this tutorial from the button beside the gravity model switch.'}
    ];
    const tutorialWelcome=document.getElementById('tutorialWelcome'),tutorialCard=document.getElementById('tutorialCard');
    const tutorialTitle=document.getElementById('tutorialTitle'),tutorialProgress=document.getElementById('tutorialProgress');
    const tutorialText=document.getElementById('tutorialText'),tutorialTry=document.getElementById('tutorialTry');
    const tutorialBack=document.getElementById('tutorialBack'),tutorialNext=document.getElementById('tutorialNext');
    const tutorialLaunch=document.getElementById('tutorialLaunch'),tutorialFullscreen=document.getElementById('tutorialFullscreen');
    const tutorialSpotlight=document.getElementById('tutorialSpotlight');
    const tutorialProgressFill=document.getElementById('tutorialProgressFill');
    const tutorialFeedback=document.getElementById('tutorialFeedback');
    const visibleMatterOnly=document.getElementById('visibleMatterOnly');
    let tutorialIndex=0,tutorialTarget=null,tutorialCompleted=[],tutorialLastMass=0,tutorialLastScale=0,tutorialLastTime=0,tutorialDrag=null,tutorialShortcutUnlocked=false;
    // A new galaxy resets the halo sliders and the timeline, so resync what the tutorial last saw.
    document.addEventListener('galaxychange',()=>{
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialLastTime=Number(timeSlider.value);
      refreshTutorialGuidance();
    });
    // Opening or closing the sandbox panel moves the page content, so re-aim the tutorial pointer.
    document.addEventListener('galaxymodechange',refreshTutorialGuidance);
    function tutorialInstruction(){
      const step=tutorialSteps[tutorialIndex];
      if(tutorialIndex>=1&&tutorialIndex<=4&&mondMode)return 'Switch back to NFW first. '+step.action;
      if(tutorialIndex===1){
        if(tutorialShortcutUnlocked)return 'The Visible Matter Only button is a shortcut to set halo mass to 0 again. Compare the orange line with the observed points.';
        if(Number(massSlider.value)===0)return 'The mass is already 0. Move it up, then drag it all the way left to 0.';
      }
      if(tutorialIndex===9){
        if(mondMode)return 'Switch back to NFW first; MOND does not use a dark matter halo. '+step.action;
        if(Number(massSlider.value)===0)return 'Increase halo mass above 0 first so there is a halo to see. '+step.action;
        if(viewPanel.classList.contains('collapsed'))return 'Open View first. '+step.action;
      }
      return step.action;
    }
    function tutorialAttentionTarget(){
      if(mondMode&&((tutorialIndex>=1&&tutorialIndex<=4)||tutorialIndex===9))return document.querySelector('#modelBar .track');
      if(tutorialIndex===9){
        if(Number(massSlider.value)===0)return massSlider;
        if(viewPanel.classList.contains('collapsed'))return viewToggle;
      }
      return tutorialTarget;
    }
    function refreshTutorialGuidance(){
      if(tutorialCard.hidden)return;
      tutorialTry.textContent='Try it: '+tutorialInstruction();
      updateTutorialSpotlight();
    }
    function scrollTutorialTarget(){
      const step=tutorialSteps[tutorialIndex],target=tutorialAttentionTarget();
      const top=target.getBoundingClientRect().top+window.scrollY-(step.gesture==='scan'||step.gesture==='drag'?16:75);
      window.scrollTo({top:Math.max(0,top),behavior:'auto'});
    }
    function completeTutorialAction(step,message){
      if(tutorialCard.hidden||tutorialIndex!==step||tutorialCompleted[step])return;
      tutorialCompleted[step]=message;
      tutorialFeedback.hidden=false;
      tutorialFeedback.textContent='✓ '+message;
      tutorialSpotlight.classList.add('is-complete');
    }
    function tutorialIsFullscreen(){return !!nativeFsElement()||orbitPlot.classList.contains('fullscreen-active');}
    function updateTutorialSpotlight(){
      if(tutorialCard.hidden||tutorialIsFullscreen()){tutorialSpotlight.hidden=true;return;}
      const step=tutorialSteps[tutorialIndex];
      let target=tutorialAttentionTarget(),gesture=step.gesture;
      if(target!==tutorialTarget){
        gesture=target===massSlider?'increase':'click';
      }else if(gesture==='switch'){
        target=document.querySelector('#modelBar .track');gesture='click';
      }else if(gesture==='guide'){
        target=document.querySelector('#theoryGuide summary');gesture='click';
      }
      const rect=target.getBoundingClientRect();
      if(rect.bottom<0||rect.top>innerHeight||rect.right<0||rect.left>innerWidth){tutorialSpotlight.hidden=true;return;}
      tutorialSpotlight.hidden=false;
      const pad=gesture==='scan'||gesture==='drag'?4:6;
      tutorialSpotlight.style.left=(rect.left-pad)+'px';
      tutorialSpotlight.style.top=(rect.top-pad)+'px';
      tutorialSpotlight.style.width=(rect.width+2*pad)+'px';
      tutorialSpotlight.style.height=(rect.height+2*pad)+'px';
      tutorialSpotlight.classList.toggle('is-click',gesture==='click');
      let fromX=rect.width*.5,toX=fromX,fromY=rect.height*.5,toY=fromY;
      if(gesture==='scan'){fromX=rect.width*.55;toX=rect.width*.82;fromY=toY=rect.height*.43;}
      if(gesture==='drag'){fromX=rect.width*.35;toX=rect.width*.62;fromY=rect.height*.5;toY=rect.height*.42;}
      if(gesture==='zero'||gesture==='increase'||gesture==='adjust'){
        const fraction=(Number(target.value)-Number(target.min))/(Number(target.max)-Number(target.min));
        fromX=8+(rect.width-16)*fraction;
        const destination=gesture==='zero'?0:gesture==='increase'?Math.min(1,Math.max(.55,fraction+.2)):(fraction>.5?.2:.8);
        toX=8+(rect.width-16)*destination;
      }
      tutorialSpotlight.style.setProperty('--from-x',(fromX+pad)+'px');
      tutorialSpotlight.style.setProperty('--to-x',(toX+pad)+'px');
      tutorialSpotlight.style.setProperty('--from-y',(fromY+pad)+'px');
      tutorialSpotlight.style.setProperty('--to-y',(toY+pad)+'px');
    }
    function syncTutorialFullscreen(){
      const fullscreen=tutorialIsFullscreen();
      const parent=fullscreen?orbitPlot:document.body;
      if(tutorialCard.parentElement!==parent)parent.appendChild(tutorialCard);
      if(tutorialSpotlight.parentElement!==parent)parent.appendChild(tutorialSpotlight);
      tutorialFullscreen.hidden=!fullscreen;
      tutorialBack.disabled=fullscreen||tutorialIndex===0;
      tutorialNext.disabled=fullscreen;
      updateTutorialSpotlight();
    }
    function showTutorialStep(){
      const step=tutorialSteps[tutorialIndex];
      const backHadFocus=document.activeElement===tutorialBack;
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialLastTime=Number(timeSlider.value);
      tutorialDrag=null;
      tutorialTarget=document.getElementById(step.target);
      tutorialProgress.textContent=`Tutorial · ${tutorialIndex+1} / ${tutorialSteps.length}`;
      tutorialProgressFill.style.width=`${100*(tutorialIndex+1)/tutorialSteps.length}%`;
      tutorialTitle.textContent=step.title;
      tutorialText.textContent=step.text;
      tutorialTry.textContent='Try it: '+tutorialInstruction();
      tutorialFeedback.hidden=!tutorialCompleted[tutorialIndex];
      tutorialFeedback.textContent=tutorialCompleted[tutorialIndex]?'✓ '+tutorialCompleted[tutorialIndex]:'';
      tutorialSpotlight.classList.toggle('is-complete',!!tutorialCompleted[tutorialIndex]);
      document.getElementById('tutorialContent').scrollTop=0;
      tutorialBack.disabled=tutorialIndex===0;
      tutorialNext.textContent=tutorialIndex===tutorialSteps.length-1?'Finish':'Next';
      tutorialCard.classList.toggle('on-left',innerWidth>960&&tutorialIndex!==7);
      document.body.classList.toggle('tutorial-view-intro',tutorialIndex===8);
      syncTutorialFullscreen();
      if(tutorialBack.disabled&&backHadFocus)tutorialNext.focus({preventScroll:true});
      if(!tutorialIsFullscreen()){
        scrollTutorialTarget();
      }
      updateTutorialSpotlight();
    }
    function startTutorial(){
      if(tutorialWelcome.open)tutorialWelcome.close();
      tutorialIndex=0;
      tutorialCompleted=[];
      visibleMatterOnly.hidden=!tutorialShortcutUnlocked;
      document.body.classList.add('tutorial-active');
      tutorialCard.hidden=false;
      showTutorialStep();
      tutorialTitle.focus({preventScroll:true});
    }
    function endTutorial(){
      tutorialTarget=null;
      tutorialCard.hidden=true;
      tutorialSpotlight.hidden=true;
      visibleMatterOnly.hidden=false;
      tutorialSpotlight.classList.remove('is-interacting');
      document.body.classList.remove('tutorial-active');
      document.body.classList.remove('tutorial-view-intro');
      (tutorialIsFullscreen()?fsBtn:tutorialLaunch).focus({preventScroll:true});
    }
    ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu'].forEach(type=>tutorialCard.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
    document.addEventListener('fullscreenchange',syncTutorialFullscreen);
    document.addEventListener('webkitfullscreenchange',syncTutorialFullscreen);
    window.addEventListener('scroll',updateTutorialSpotlight,{passive:true});
    window.addEventListener('resize',updateTutorialSpotlight);
    [massSlider,scaleSlider].forEach(slider=>{
      slider.addEventListener('pointerdown',()=>tutorialSpotlight.classList.add('is-interacting'));
      slider.addEventListener('change',updateTutorialSpotlight);
    });
    massSlider.addEventListener('input',()=>{
      const value=Number(massSlider.value);
      const crossedZero=(value===0)!==(tutorialLastMass===0);
      if(!mondMode){
        if(!tutorialCard.hidden&&value===0&&!tutorialShortcutUnlocked){
          tutorialShortcutUnlocked=true;
          visibleMatterOnly.hidden=false;
          visibleMatterOnly.classList.add('tutorial-reveal');
        }
        if(value===0)completeTutorialAction(1,'Visible matter only selected.');
        if(value>tutorialLastMass)completeTutorialAction(2,'Halo mass increased.');
        if(value!==tutorialLastMass)completeTutorialAction(4,'Model settings adjusted.');
      }
      tutorialLastMass=value;
      if(crossedZero&&(tutorialIndex===1||tutorialIndex===9))refreshTutorialGuidance();
    });
    massSlider.addEventListener('change',()=>{
      if(!tutorialCard.hidden&&tutorialIndex===9&&!mondMode&&Number(massSlider.value)>0)scrollTutorialTarget();
    });
    scaleSlider.addEventListener('input',()=>{
      const value=Number(scaleSlider.value);
      if(!mondMode&&value!==tutorialLastScale){
        completeTutorialAction(3,'Halo shape adjusted.');
        completeTutorialAction(4,'Model settings adjusted.');
      }
      tutorialLastScale=value;
    });
    document.getElementById('visibleMatterOnly').addEventListener('click',()=>{
      if(!mondMode)completeTutorialAction(1,'Visible matter only selected.');
      tutorialLastMass=Number(massSlider.value);
    });
    document.getElementById('applyObservedFit').addEventListener('click',()=>{
      if(!mondMode)completeTutorialAction(4,'Starting fit restored.');
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
    });
    document.getElementById('playPauseBtn').addEventListener('click',()=>{
      if(playing)completeTutorialAction(6,'Simulation playing.');
    });
    timeSlider.addEventListener('input',()=>{
      const value=Number(timeSlider.value);
      if(value!==tutorialLastTime)completeTutorialAction(6,'Timeline adjusted.');
      tutorialLastTime=value;
    });
    orbitPlot.addEventListener('pointerdown',e=>{
      if(!tutorialCard.hidden&&tutorialIndex===7&&!e.target.closest('#viewPanel'))tutorialDrag={id:e.pointerId,x:e.clientX,y:e.clientY};
    });
    orbitPlot.addEventListener('pointermove',e=>{
      if(tutorialDrag&&e.pointerId===tutorialDrag.id&&Math.hypot(e.clientX-tutorialDrag.x,e.clientY-tutorialDrag.y)>6){
        completeTutorialAction(7,'View explored.');
        tutorialDrag=null;
      }
    });
    window.addEventListener('pointerup',()=>{tutorialDrag=null;});
    window.addEventListener('pointercancel',()=>{tutorialDrag=null;});
    ['viewReset','viewFaceOn','viewEdgeOn','toggleOrbits','viewFullHalo'].forEach(id=>{
      document.getElementById(id).addEventListener('click',()=>completeTutorialAction(7,'View explored.'));
    });
    viewToggle.addEventListener('click',()=>{
      completeTutorialAction(8,'View controls toggled.');
      refreshTutorialGuidance();
      if(!tutorialCard.hidden&&tutorialIndex===9)scrollTutorialTarget();
    });
    fullHaloBtn.addEventListener('click',()=>completeTutorialAction(9,'Full halo view selected.'));
    document.querySelectorAll('#theoryGuide details').forEach(details=>{
      details.addEventListener('toggle',()=>{if(details.open)completeTutorialAction(10,'Theory guide opened.');});
    });
    window.addEventListener('pointerup',()=>tutorialSpotlight.classList.remove('is-interacting'));
    window.addEventListener('pointercancel',()=>tutorialSpotlight.classList.remove('is-interacting'));
    [visibleMatterOnly,document.getElementById('applyObservedFit')].forEach(button=>button.addEventListener('click',refreshTutorialGuidance));
    mondToggle.addEventListener('change',()=>{
      if(mondMode)completeTutorialAction(5,'MOND selected.');
      refreshTutorialGuidance();
      if(!tutorialCard.hidden&&tutorialIndex===9)scrollTutorialTarget();
    });
    document.getElementById('tutorialStart').addEventListener('click',startTutorial);
    document.getElementById('tutorialWelcomeSkip').addEventListener('click',()=>tutorialWelcome.close());
    tutorialLaunch.addEventListener('click',startTutorial);
    document.getElementById('tutorialSkip').addEventListener('click',endTutorial);
    tutorialBack.addEventListener('click',()=>{if(tutorialIndex>0){tutorialIndex--;showTutorialStep();}});
    tutorialNext.addEventListener('click',()=>{
      if(tutorialIndex===tutorialSteps.length-1){endTutorial();return;}
      tutorialIndex++;showTutorialStep();
    });
    document.addEventListener('keydown',e=>{
      if(tutorialCard.hidden||tutorialIsFullscreen()||e.repeat||e.altKey||e.ctrlKey||e.metaKey||e.shiftKey)return;
      if(e.key!=='ArrowLeft'&&e.key!=='ArrowRight')return;
      const control=e.target instanceof Element?e.target.closest('input,select,textarea,button,a,summary,[contenteditable]:not([contenteditable="false"]),[role="slider"],[role="radio"]'):null;
      if(control&&control!==tutorialBack&&control!==tutorialNext)return;
      e.preventDefault();
      (e.key==='ArrowRight'?tutorialNext:tutorialBack).click();
    });
    tutorialWelcome.showModal();
  </script>
</body>
</html>
"""


def join_names(names: list[str]) -> str:
    """Join names for prose: 'A', 'A and B', 'A, B and C'."""
    if len(names) < 2:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def main() -> None:
    if not TRAJECTORY_PATH.is_file():
        raise FileNotFoundError(f"Missing trajectory CSV: {TRAJECTORY_PATH}")
    if not SPARC_PATH.is_file():
        raise FileNotFoundError(f"Missing SPARC table: {SPARC_PATH}")

    stars = load_trajectories(TRAJECTORY_PATH)
    galaxies = [build_galaxy_record(preset, SPARC_PATH) for preset in GALAXY_PRESETS]
    galaxies_js = render_galaxies_js(galaxies)
    GALAXIES_JS_PATH.write_text(galaxies_js, encoding="utf-8")
    payload = {
        "stars": stars,
        "upsilon_disk": UPSILON_DISK,
        "upsilon_bulge": UPSILON_BULGE,
        "frame_step_myr": 1.0,
        "total_time_myr": 1000.0,
        "default_galaxy_id": DEFAULT_PRESET.sparc_id,
        "stable_concentration": list(STABLE_CONCENTRATION),
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
    OUTPUT_PATH.write_text(html, encoding="utf-8")
    size_mb = OUTPUT_PATH.stat().st_size / (1024 * 1024)
    print(f"Generated {OUTPUT_PATH.name} ({size_mb:.2f} MiB) and {GALAXIES_JS_PATH.name}")
    print(f"Embedded {len(stars)} stars and {len(galaxies)} SPARC galaxies")
    for galaxy in galaxies:
        halo = galaxy["default_halo"]
        print(
            f"  {galaxy['name']}: {len(galaxy['rotation']['radius_kpc'])} points, "
            f"NFW starting fit M200={halo['m200_msun']:.3g} Msun, "
            f"rs={halo['scale_radius_kpc']:.1f} kpc"
        )



if __name__ == "__main__":
    main()