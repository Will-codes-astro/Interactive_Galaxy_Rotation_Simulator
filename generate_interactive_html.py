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
    .ml-col { position:relative; }
    .stat.ml .stat-label { display:flex; align-items:center; justify-content:space-between; gap:8px; }
    .ml-info { padding:0 6px; border:1px solid rgba(242,200,121,0.45); border-radius:999px; background:transparent; color:var(--gold); font:600 10px/1.5 "Segoe UI",sans-serif; letter-spacing:0; text-transform:none; cursor:pointer; }
    .ml-info:hover, .ml-info[aria-expanded="true"] { background:rgba(242,200,121,0.16); }
    .ml-info:focus-visible { outline:2px solid var(--cyan); outline-offset:2px; }
    .stat.ml.off { display:none; }
    [data-preset-only][hidden], .stability[hidden] { display:none !important; }
    .ml-explainer { position:absolute; top:calc(100% + 8px); right:0; z-index:40; width:340px; max-width:calc(100vw - 32px); padding:14px 16px; border:1px solid rgba(242,200,121,0.35); border-radius:8px; background:#0c1520; box-shadow:0 12px 32px rgba(0,0,0,0.55); font-size:12.5px; line-height:1.5; color:var(--ink); white-space:normal; }
    .ml-explainer { opacity:0; visibility:hidden; transform:translateY(-4px); pointer-events:none; transition:opacity .2s ease, transform .2s ease, visibility 0s linear .2s; }
    .ml-explainer.open { opacity:1; visibility:visible; transform:none; pointer-events:auto; transition:opacity .2s ease, transform .2s ease, visibility 0s; }
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
    .source-toggle input:checked + span { background:rgba(138,180,248,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(138,180,248,.5); }
    .source-toggle input:focus-visible + span { outline:2px solid var(--cyan); outline-offset:2px; }
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
    .stat-q:hover::after,.stat-q:focus-visible::after { opacity:1; visibility:visible; }
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
    .tm-row { display:flex; justify-content:space-between; gap:9px; }
    .tm-row b { color:rgba(255,255,255,.78); font-weight:600; font-variant-numeric:tabular-nums; }
    /* χ² readout: fullscreen only (the curve panel shows it otherwise), and only for SPARC presets. */
    .tm-chi { display:none; }
    #orbitPlot:fullscreen .tm-chi:not([hidden]), #orbitPlot:-webkit-full-screen .tm-chi:not([hidden]), #orbitPlot.fullscreen-active .tm-chi:not([hidden]) { display:flex; }
    .tm-chi b.good { color:#8fe3ee; } .tm-chi b.ok { color:#ffc48a; } .tm-chi b.poor { color:#ff9d94; }
    .stability-chip { position:absolute; top:44px; left:12px; z-index:100; max-width:min(300px,calc(100% - 190px)); padding:5px 10px; border-radius:8px; border:1px solid rgba(255,138,128,0.6); background:rgba(40,10,10,0.6); -webkit-backdrop-filter:blur(10px); backdrop-filter:blur(10px); color:#ffb4ac; font-size:11px; font-weight:600; line-height:1.35; pointer-events:none; }
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
    .play-btn { display:inline-flex; align-items:center; justify-content:center; gap:6px; min-width:44px; min-height:44px; padding:0 14px; flex:none; font-size:13px; border-radius:22px; background:rgba(138,180,248,.12); border-color:var(--cyan); touch-action:manipulation; }
    .play-btn[aria-pressed=true] { background:rgba(138,180,248,.28); }
    @media(max-width:560px) { .play-btn .play-text { display:none; } .play-btn { width:44px; padding:0; } }
    .fullscreen-active { position:fixed !important; inset:0; z-index:9999; width:100vw !important; height:100dvh !important; max-height:none !important; min-height:0 !important; border:0 !important; border-radius:0 !important; background:#020206; }
    body.fullscreen-lock { overflow:hidden; }
    .fs-sliders { display:none; position:absolute; left:12px; bottom:12px; z-index:3; width:min(240px,calc(100% - 24px)); padding:8px 10px; gap:4px; flex-direction:column; border-radius:10px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.55); -webkit-backdrop-filter:blur(14px); backdrop-filter:blur(14px); color:var(--ink); font-size:11px; cursor:default; }
    #orbitPlot:fullscreen .fs-sliders, #orbitPlot:-webkit-full-screen .fs-sliders, #orbitPlot.fullscreen-active .fs-sliders { display:flex; }
    .fs-sliders label { display:flex; justify-content:space-between; gap:8px; color:var(--muted); }
    .fs-sliders output { color:var(--ink); font-variant-numeric:tabular-nums; }
    .fs-sliders input[type=range] { width:100%; height:16px; margin:0; }
    .fs-sliders.disabled :is(label,input) { opacity:.35; pointer-events:none; }
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
    .view-panel .view-seg { display:grid; grid-template-columns:repeat(3,1fr); padding:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
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
    /* Outside fullscreen the panel sits over a small plot, so it drops the section labels and tightens up. */
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel:not(.collapsed) { min-width:0; width:172px; padding:7px; top:8px; right:8px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-label { display:none; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel:not(.collapsed) > :not(.view-head) { margin-top:5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-row { gap:5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-row button { padding:4px 4px; font-size:11px; white-space:nowrap; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel :is(.view-seg,.speed-toggle,.mode-toggle) button { padding:3px 2px; font-size:10.5px; }
    #orbitPlot:not(.fullscreen-active):not(:fullscreen) .view-panel .view-switches { grid-template-columns:auto auto; justify-content:space-between; gap:0 6px; }
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
    .halo-metric { padding:10px; border:1px solid rgba(179,154,232,.14); border-radius:6px; background:rgba(5,9,16,.4); }
    .halo-metric button { display:block; width:100%; padding:0; border:0; background:none; text-align:left; color:#c3cfdb; font-size:11px; }
    .halo-metric button:focus-visible { outline:2px solid var(--cyan); outline-offset:3px; border-radius:2px; }
    .halo-metric .info-icon { color:#b39ae8; margin-left:3px; }
    .halo-metric .metric-value { display:block; margin-top:5px; color:#e6ddfa; font-size:18px; font-weight:600; font-variant-numeric:tabular-nums; }
    .halo-metric dd { display:grid; grid-template-rows:0fr; margin:0; color:var(--muted); font-size:11px; line-height:1.45; opacity:0; visibility:hidden; transition:grid-template-rows .2s ease, opacity .2s ease, visibility 0s linear .2s; }
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
    .tutorial-eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.1em; font-size:10px; font-weight:600; }
    .tutorial-welcome h2 { margin:8px 0; font-size:22px; }
    .tutorial-welcome p { color:#b5c5d5; margin:0 0 18px; }
    .tutorial-welcome-actions { display:flex; flex-wrap:wrap; gap:9px; }
    .tutorial-primary { background:rgba(138,180,248,.16); border-color:var(--cyan); color:#e4fbff; font-weight:600; }
    .tutorial-card { position:fixed; z-index:10000; right:16px; bottom:16px; display:flex; flex-direction:column; width:min(380px,calc(100% - 32px)); max-height:min(65vh,480px); overflow:hidden; padding:18px; border:1px solid rgba(138,180,248,.24); border-radius:16px; background:linear-gradient(155deg,#142838,#0b1623 70%); box-shadow:0 16px 48px rgba(0,0,0,.58); cursor:default; user-select:text; }
    .tutorial-card.on-left { left:16px; right:auto; }
    .tutorial-card[hidden] { display:none; }
    .tutorial-top { display:flex; flex:none; align-items:center; justify-content:space-between; gap:12px; }
    .tutorial-content { min-height:0; overflow-y:auto; overscroll-behavior:contain; }
    .tutorial-skip { border:0; padding:4px; background:transparent; color:var(--muted); font-size:11px; text-decoration:underline; text-underline-offset:3px; }
    .tutorial-skip:hover, .tutorial-skip:focus-visible { color:var(--ink); }
    .tutorial-progress-track { flex:none; height:3px; margin:10px 0 4px; background:rgba(138,180,248,.13); border-radius:4px; }
    .tutorial-progress-fill { height:100%; width:0; border-radius:inherit; background:var(--cyan); transition:width .25s ease; }
    .tutorial-card h2 { margin:10px 0 7px; font-size:18px; line-height:1.25; }
    .tutorial-card p { margin:0 0 10px; color:#c3d0de; font-size:13px; line-height:1.5; }
    .tutorial-card .tutorial-try { margin:2px 0 0; padding:10px 12px; border-left:2px solid var(--cyan); border-radius:0 7px 7px 0; background:rgba(138,180,248,.07); color:#c7f5fa; }
    .tutorial-feedback { flex:none; margin:9px 0 0; color:#a5e3c2; font-size:12px; }
    .tutorial-feedback[hidden] { display:none; }
    .tutorial-actions { display:flex; flex:none; flex-wrap:wrap; justify-content:flex-end; gap:8px; margin-top:12px; }
    .tutorial-keys { align-self:center; margin-right:auto; color:var(--muted); font-size:11px; }
    @media (hover:none) and (pointer:coarse) { .tutorial-keys { display:none; } }
    .tutorial-actions button { min-height:36px; }
    .tutorial-actions button:disabled { opacity:.45; cursor:not-allowed; }
    .tutorial-spotlight { position:fixed; z-index:9998; pointer-events:none; border:2px solid #8ab4f8; border-radius:9px; box-shadow:0 0 0 9999px rgba(2,7,14,.30), 0 0 0 4px rgba(138,180,248,.16); transition:opacity .45s ease; }
    .tutorial-spotlight[hidden] { display:none; }
    .tutorial-pointer { position:absolute; left:0; top:0; width:25px; height:30px; filter:drop-shadow(0 2px 3px #000); animation:tutorial-gesture 2.8s ease-in-out infinite; will-change:transform; transition:opacity .5s ease; }
    .tutorial-spotlight.is-interacting .tutorial-pointer { opacity:.15; }
    .tutorial-spotlight.is-reached .tutorial-pointer { opacity:0; animation-play-state:paused; }
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
    /* ===== Design system: tokens and component layer (visual only; no physics). ===== */
    :root { --bg:#060a12; --panel:#0c1322; --panel-2:#101929; --line:rgba(148,170,200,.14); --line-strong:rgba(148,170,200,.28); --ink:#e7edf6; --ink-2:#b8c4d6; --muted:#8392a8; --accent:#8ab4f8; --cyan:var(--accent); --gold:#e9c88a;
      --r:10px; --r-lg:14px; --f-ui:Inter,"Segoe UI",system-ui,sans-serif; --f-serif:"Source Serif 4",Georgia,serif; --f-mono:"JetBrains Mono",ui-monospace,Consolas,monospace; --shadow:inset 0 1px 0 rgba(255,255,255,.03),0 12px 32px rgba(0,0,0,.28); }
    @media (prefers-reduced-motion:no-preference) { html { scroll-behavior:smooth; } }
    body { font:14px/1.55 var(--f-ui); -webkit-font-smoothing:antialiased; background:radial-gradient(1100px 520px at 80% -8%,rgba(70,100,170,.16),transparent 62%),radial-gradient(800px 460px at -8% 6%,rgba(120,90,180,.08),transparent 60%),var(--bg); background-repeat:no-repeat; min-height:100vh; }
    #top, #simulator, #results, #curvePanel, #theoryGuide, #method { scroll-margin-top:72px; }
    .topbar { position:sticky; top:0; z-index:60; display:flex; align-items:center; justify-content:space-between; gap:16px; height:56px; padding:0 28px; border-bottom:1px solid var(--line); background:rgba(6,10,18,.74); -webkit-backdrop-filter:blur(14px) saturate(140%); backdrop-filter:blur(14px) saturate(140%); }
    .brand { display:inline-flex; align-items:center; gap:10px; color:var(--ink); font-weight:600; font-size:14px; letter-spacing:.01em; text-decoration:none; white-space:nowrap; }
    .brand svg { width:24px; height:24px; color:var(--accent); flex:none; }
    .topnav { display:flex; gap:2px; overflow-x:auto; scrollbar-width:none; }
    .topnav a { padding:6px 11px; border-radius:7px; color:var(--muted); font-size:13px; font-weight:500; text-decoration:none; white-space:nowrap; transition:color .15s,background .15s; }
    .topnav a:hover, .topnav a:focus-visible { color:var(--ink); background:rgba(255,255,255,.05); outline:none; }
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
    .model-bar, .controls, .sandbox, .panel, .learn details { border:1px solid var(--line); border-radius:var(--r-lg); background:linear-gradient(180deg,var(--panel-2),var(--panel)); box-shadow:var(--shadow); }
    .model-bar { padding:12px 18px; margin-bottom:12px; gap:14px 20px; }
    .model-bar .model-label { font:500 10.5px var(--f-mono); letter-spacing:.14em; }
    .controls { gap:22px 30px; padding:20px 22px; margin-bottom:0; }
    .control label { color:var(--ink-2); font-size:12.5px; font-weight:500; }
    .control output, .fs-sliders output { font-family:var(--f-mono); font-size:12px; color:var(--ink); }
    .control-help { font-size:12px; line-height:1.5; }
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
    .galaxy-bar select { padding:7px 34px 7px 12px; border:1px solid var(--line-strong); border-radius:8px; background-color:var(--panel-2); font:500 13px var(--f-ui); }
    .galaxy-summary { color:var(--ink-2); font-size:12.5px; }
    .source-toggle span { padding:6px 16px; font-weight:500; }
    .switch .opt { font-size:13px; }
    .track-input { border-color:var(--line-strong); background:rgba(0,0,0,.25); font-family:var(--f-mono); }
    .sandbox-group h3 { font:500 10.5px var(--f-mono); letter-spacing:.14em; color:var(--accent); }
    .sandbox-summary, .key, .halo-info, .stability, .galaxy-type { border-radius:var(--r); }
    .sandbox-summary dd { font-family:var(--f-mono); }
    .galaxy-type h3 { font:600 18px var(--f-serif); }
    .panel-grid { gap:16px; }
    .panel-heading { align-items:center; padding:16px 18px 8px; }
    .panel-heading h2 { font:600 15px var(--f-ui); color:var(--ink); }
    .panel-heading .fig { margin-right:10px; color:var(--accent); font:500 11px var(--f-mono); letter-spacing:.08em; text-transform:uppercase; }
    .panel-heading > span:not(#chiBadge) { color:var(--muted); font-size:12px; }
    #chiBadge { font:500 11.5px var(--f-mono); }
    .key { background:rgba(148,170,200,.035); }
    .key summary { font:500 10.5px var(--f-mono); letter-spacing:.12em; color:var(--ink-2); }
    .sw.ldash { background:none; height:0; border-top:2.5px dashed; }
    .halo-info h3 { font:600 14px var(--f-ui); }
    .halo-metric .metric-value { font:500 17px var(--f-mono); }
    .stability-c { font-family:var(--f-mono); }
    .learn { gap:16px; margin-top:0; }
    .learn details { padding:0 22px; }
    .learn summary { padding:16px 0; font:600 18px var(--f-serif); }
    .learn p { color:var(--ink-2); font-size:13.5px; line-height:1.65; }
    .learn p.tldr { color:var(--ink); font-size:14px; }
    .schematic { margin:4px 0 16px; padding:12px 14px 10px; border:1px solid var(--line); border-radius:var(--r); background:rgba(0,0,0,.22); color:var(--muted); }
    .schematic svg { display:block; width:100%; max-width:440px; height:auto; margin:auto; font:10.5px var(--f-ui); }
    .schematic text:not([fill]) { fill:currentColor; }
    .schematic figcaption { margin-top:6px; font-size:11.5px; line-height:1.5; }
    .eq { display:block; margin:8px 0 12px; padding:10px 12px; border-left:2px solid var(--line-strong); border-radius:var(--r); background:rgba(0,0,0,.22); overflow-x:auto; color:var(--ink); }
    .eq math { font-family:"STIX Two Math","Cambria Math","Latin Modern Math",math,serif; font-size:1.25em; margin:0; }
    .method-grid .eq { margin:6px 0; padding:6px 10px; }
    .method-grid .eq math { font-size:1.1em; }
    footer { flex-wrap:wrap; align-items:baseline; gap:8px 24px; margin:28px 0 0; padding:18px 22px; border:1px solid var(--line); border-radius:var(--r-lg); background:var(--panel); font-size:12px; }
    .footer-title { flex-basis:100%; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    footer { display:block; margin-top:0; }
    @media (max-width:900px) { .method-grid { grid-template-columns:1fr 1fr !important; } } @media (max-width:560px) { .method-grid { grid-template-columns:1fr !important; } }
    .method-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:16px 28px; margin:0 0 16px; }
    .method-grid dt { margin-bottom:4px; color:var(--accent); font:500 10.5px var(--f-mono); letter-spacing:.14em; text-transform:uppercase; }
    .method-grid dd { margin:0; color:var(--ink-2); font-size:12.5px; line-height:1.55; }
    .tutorial-welcome, .tutorial-card { border-color:rgba(138,180,248,.24); background:linear-gradient(155deg,#121b2e,#0a101c 70%); }
    .tutorial-welcome h2, .tutorial-card h2 { font-family:var(--f-serif); }
    .tutorial-eyebrow { font-family:var(--f-mono); }
    .stat-q::after { background:var(--panel-2); border-color:var(--line-strong); }
    .telemetry { font-family:var(--f-mono); }
    @media (max-width:900px) { header { grid-template-columns:1fr; padding-top:36px; } .header-stats { justify-content:flex-start; flex-wrap:wrap; } }
    @media (max-width:640px) { .topbar { padding:0 12px; gap:8px; } .brand span { display:none; } .topnav a { padding:6px 8px; font-size:12.5px; } header { padding:28px 16px 4px; } main { padding:0 12px 32px; } .section-head { margin-top:28px; } .section-head h2 { font-size:19px; } .controls { padding:16px; } .learn details { padding:0 16px; } footer { padding:16px; } .stat { flex:1 1 140px; } .panel-heading > span:not(#chiBadge) { display:none; } }
    .key summary::before { content:""; display:inline-block; width:4px; }
  </style>
</head>
<body>
  <nav class="topbar" aria-label="Page sections">
    <a class="brand" href="#top"><svg viewBox="0 0 32 32" aria-hidden="true"><ellipse cx="16" cy="16" rx="13" ry="5.5" transform="rotate(-24 16 16)" fill="none" stroke="currentColor" stroke-width="1.4" opacity=".55"/><ellipse cx="16" cy="16" rx="8" ry="3.2" transform="rotate(-24 16 16)" fill="none" stroke="currentColor" stroke-width="1.4" opacity=".85"/><circle cx="16" cy="16" r="2.4" fill="currentColor"/></svg><span>Galaxy Rotation Simulator</span></a>
    <div class="topnav"><a href="#simulator">Model Setup</a><a href="#results">Simulation</a><a href="#curvePanel">Rotation curve</a><a href="#theoryGuide">Theory</a><a href="#method">Method &amp; data</a></div>
  </nav>
  <header id="top">
    <div class="intro"><div class="eyebrow">SPARC galaxy dynamics · interactive model</div><h1>Interactive Galaxy Rotation Simulator</h1><p class="lead">Compare a dark matter halo (NFW) with modified gravity (MOND) against measured rotation curves from the SPARC survey (Lelli, McGaugh &amp; Schombert 2016).</p></div>
    <div class="header-stats" aria-label="Simulation summary">
      <div class="stat"><div class="stat-label">Test stars</div><div class="stat-value"><span id="starCountStat">800</span><small>massless</small></div></div>
      <div class="stat"><div class="stat-label">Timeline</div><div class="stat-value">1 Gyr<small>1000 Myr</small></div></div>
      <div class="ml-col">
        <div class="stat ml" id="mlStat" title="Stellar mass-to-light ratios used for the baryonic (visible-matter) trajectories"><div class="stat-label"><span>Mass-to-light (M/L)</span><button type="button" class="ml-info" id="mlInfoBtn" aria-expanded="false" aria-controls="mlExplainer" title="What is a mass-to-light ratio?">?</button></div><div class="stat-value ml-pair"><div><b id="mlDisk">0.50</b><span>DISK</span></div><div><b id="mlBulge">0.70</b><span>BULGE</span></div></div></div>
        <div class="ml-explainer" id="mlExplainer" role="dialog" aria-labelledby="mlExplainerTitle" aria-hidden="true">
          <button type="button" class="ml-close" id="mlExplainerClose" aria-label="Close">×</button>
          <h3 id="mlExplainerTitle">What is a mass-to-light ratio?</h3>
          <p>Telescopes measure a galaxy's <b>light</b>, not its mass. The mass-to-light ratio (M/L, or the Greek letter upsilon, Υ) converts one into the other: it is how many Suns' worth of mass the stars hold for every Sun's worth of light they give off.</p>
          <ul>
            <li><b>Low M/L means young stars.</b> Young, massive blue stars are very bright for their mass.</li>
            <li><b>High M/L means old stars.</b> Old populations are mostly small, dim red stars: lots of mass, not much light.</li>
          </ul>
          <p>The presets use <b>0.5 for the disk</b> and <b>0.7 for the bulge</b>, which is older (both measured in 3.6 μm infrared light). These are standard average values, not measured for each galaxy, and every real galaxy is a little different. Changing them raises or lowers the visible-matter curve, which changes how much dark matter is needed.</p>
          <p>The Custom Sandbox sets the masses directly, so it does not need an M/L ratio.</p>
          <p class="ml-source">Source: Lelli, McGaugh &amp; Schombert (2016), SPARC; McGaugh &amp; Schombert (2014).</p>
        </div>
        <label class="fit-toggle" id="mondFitWrap" hidden title="Re-fit a₀ and the disk M/L to the SPARC data to minimise χ² (MOND mode only)"><span class="fit-label">MOND best fit</span><input type="checkbox" id="mondFitToggle" autocomplete="off"><span class="track" aria-hidden="true"></span></label>
      </div>
    </div>
  </header>
  <main>
    <div class="section-head" id="simulator"><span class="kicker">01</span><h2>Model setup</h2><p>Choose a galaxy, a theory of gravity and the halo parameters.</p></div>
    <section class="model-bar galaxy-bar" id="galaxyBar" aria-label="Galaxy source">
      <div class="source-toggle" id="sourceToggle" role="radiogroup" aria-label="Galaxy source">
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
          <div class="control" id="sbBulgeRadiusControl"><label for="sbBulgeRadius"><span>Effective radius, r<sub>b</sub></span><output id="sbBulgeRadiusValue" for="sbBulgeRadius"></output></label><input id="sbBulgeRadius" type="range" min="0.2" max="2" step="0.1" value="1" autocomplete="off"><p class="control-help" id="sbBulgeRadiusHelp">Size of the bulge (Hernquist profile). Half its mass lies within about 2.4 × r<sub>b</sub>.</p></div>
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
    <section class="model-bar" id="modelBar" aria-label="Gravity model">
      <span class="model-label">Gravity model</span>
      <label class="switch" title="Switch between a dark matter halo and modified gravity (MOND)"><input type="checkbox" id="mondToggle" role="switch" autocomplete="off" aria-label="Use MOND theory instead of NFW dark matter halo"><span class="opt active" id="optNfw">NFW Dark Matter Halo</span><span class="track" aria-hidden="true"></span><span class="opt" id="optMond">MOND Theory</span></label>
      <span class="mond-note" id="mondNote" hidden>MOND active: no dark matter (M<sub>200</sub> = 0). Gravity is boosted below a<sub>0</sub> = <span id="mondA0Text">1.2×10⁻¹⁰</span> m/s².</span>
      <button type="button" class="tutorial-launch" id="tutorialLaunch">Tutorial</button>
    </section>
    <section class="controls" aria-label="Dark matter controls">
    <div class="control" id="massControl"><label for="massSlider"><span>Halo mass, M<sub>200</sub></span><output id="massValue"></output></label><input id="massSlider" type="range" min="0.00" max="2.00" step="0.01"><p class="control-help" id="massHelp">How much dark matter surrounds the galaxy. More mass means stronger gravity, so stars orbit faster and the orange curve rises. Set it to 0 for visible matter only.</p><div class="fit-actions"><button id="applyObservedFit" class="fit-button" type="button">Apply Observed Fit</button><button id="visibleMatterOnly" class="fit-button" type="button" title="Set the dark matter halo mass to zero">Visible Matter Only</button></div></div>
    <div class="control" id="scaleControl"><label for="scaleSlider"><span>Scale radius, r<sub>s</sub></span><output id="scaleValue"></output></label><input id="scaleSlider" type="range" min="5" max="60" step="0.5"><p class="control-help" id="scaleHelp">How spread out the halo is, for the same total mass. A small radius packs the dark matter near the centre and speeds up the inner stars; a large radius spreads it out.</p></div>
      <div class="control"><label for="timeSlider"><span>Orbit time</span><output id="timeValue"></output></label><div class="timeline-row"><button type="button" id="playPauseBtn" class="play-btn" aria-pressed="false" aria-label="Play simulation" title="Play / pause (Space)"><span class="play-icon" aria-hidden="true">▶</span><span class="play-text">Play</span></button><input id="timeSlider" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"></div><p class="control-help">Moves the stars through one billion years of orbits. Drag it to any moment, or press Play.</p></div>
      <div class="control"><label for="trackInput"><span>Track a star</span><output id="trackMax"></output></label><input id="trackInput" class="track-input" type="number" min="1" step="1" inputmode="numeric" placeholder="Star number, e.g. 120"><p class="control-help">Rings one star and counts its laps round the centre since 0 Myr, with its speed and distance. Clear the box to stop.</p><p class="control-help">Stars are numbered from the centre outwards: star 1 starts closest to the centre and star <span id="trackLast">800</span> starts furthest out.</p></div>
    </section>
    <div class="section-head" id="results"><span class="kicker">02</span><h2>Simulation</h2><p>Star orbits in 3D, and the model's rotation curve against the measured speeds.</p></div>
    <section class="panel-grid">
    <article class="panel"><div class="panel-heading"><h2><span class="fig">Fig. 1</span>3D galaxy simulation</h2><span>Test-star orbits over 1 Gyr</span></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional __DEFAULT_GALAXY_NAME__ star trajectories"><div id="escape-counter" class="telemetry" aria-label="Live simulation readout"><div class="tm-row"><span>Ejected</span><b><span id="escape-count-val">0</span> / <span id="star-total-val">800</span></b></div><div class="tm-row"><span>Halo r<sub>200</sub></span><b id="tm-r200">—</b></div><div class="tm-row tm-chi" id="tmChiRow"><span>χ²ᵣ</span><b id="tm-chi">—</b></div></div><div class="stability-chip" id="stabilityChip" role="status" hidden></div><div class="view-panel" id="viewPanel" role="group" aria-label="3D view controls"><button type="button" class="view-head" id="viewToggle" aria-expanded="true" title="Collapse or expand the view controls"><span>View</span><span class="chev" aria-hidden="true">▾</span></button><div class="view-row view-primary"><button type="button" id="toggleFullscreen" aria-pressed="false" title="Toggle fullscreen view of 3D simulation">⛶ Fullscreen</button><button type="button" id="viewReset" title="Return to the default __DEFAULT_GALAXY_NAME__ view">↺ Reset</button></div><div class="view-label">Camera</div><div class="view-seg"><button type="button" id="viewFaceOn" title="Look straight down on the galactic plane">Face-on</button><button type="button" id="viewEdgeOn" title="Look along the galactic plane">Edge-on</button><button type="button" id="viewFullHalo" title="Pull back to see the full dark matter halo around the tiny disk">Full halo</button></div><div class="view-label">Show</div><div class="view-switches"><button type="button" role="switch" id="toggleOrbits" aria-checked="false" aria-pressed="false" title="Show or hide the orbit trajectory lines"><span class="sw" aria-hidden="true"></span>Orbits</button><button type="button" role="switch" id="toggleGas" aria-checked="true" aria-pressed="true" title="Show or hide the faint cold-gas cloud around the stars"><span class="sw" aria-hidden="true"></span>Gas</button><button type="button" role="switch" id="toggleHalo" aria-checked="true" aria-pressed="true" title="Show or hide the dark matter halo (the halo is not drawn in MOND mode)"><span class="sw" aria-hidden="true"></span>Halo</button><button type="button" role="switch" id="toggleEjection" class="sw-eject" aria-checked="false" aria-pressed="false" title="Show the red ejection line: stars beyond 1.5 × R_max count as ejected"><span class="sw" aria-hidden="true"></span>Ejection line</button></div><div class="view-label">Speed</div><div class="speed-toggle" role="radiogroup" aria-label="Playback speed"><button type="button" role="radio" aria-checked="false" data-speed="0.25" title="0.25x speed">0.25x</button><button type="button" role="radio" aria-checked="false" data-speed="0.5" title="0.5x speed">0.5x</button><button type="button" role="radio" aria-checked="true" data-speed="1" title="1x speed">1x</button><button type="button" role="radio" aria-checked="false" data-speed="2" title="2x speed">2x</button></div><div class="mode-toggle" role="radiogroup" aria-label="Mouse drag mode" title="What a mouse drag does"><button type="button" role="radio" aria-checked="true" data-mode="rotate">Rotate</button><button type="button" role="radio" aria-checked="false" data-mode="pan">Pan</button></div></div><div class="fs-sliders" id="fsSliders" role="group" aria-label="Fullscreen controls"><div class="fs-time"><button type="button" id="fsPlay" class="fs-play" aria-pressed="false" aria-label="Play simulation" title="Play / pause (Space)">▶</button><button type="button" id="fsRestart" class="fs-play" aria-label="Restart from 0 Myr" title="Back to 0 Myr and play">↺</button><input id="fsTrack" class="track-input" type="number" min="1" step="1" inputmode="numeric" placeholder="Star #" aria-label="Track a star by number" title="Type a star number to ring it and count its laps"><output id="fsTimeValue">0 / 1000 Myr</output></div><label for="fsMassSlider"><span>Halo mass</span><output id="fsMassValue"></output></label><input id="fsMassSlider" type="range" min="0.00" max="2.00" step="0.01"><label for="fsScaleSlider"><span>Scale radius</span><output id="fsScaleValue"></output></label><input id="fsScaleSlider" type="range" min="5" max="60" step="0.5"></div></div>
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
      <article class="panel" id="curvePanel">
        <div class="panel-heading"><h2><span class="fig">Fig. 2</span>Rotation curve</h2><span id="chiBadge" title="Goodness of fit of baryons + NFW halo to the SPARC points">χ² –</span></div>
        <div id="curvePlot" class="plot" role="img" aria-label="__DEFAULT_GALAXY_NAME__ observed, baryonic, and fitted rotation curves"></div>
        <div class="stability" id="stabilityMonitor" data-state="stable" aria-label="Halo stability monitor" hidden>
          <div class="stability-head"><span class="stability-title">Halo stability monitor</span><span class="stability-badge" id="stabilityBadge" role="status"></span></div>
          <div class="stability-reading"><span class="stability-c" id="stabilityC">c = —</span><span class="stability-detail" id="stabilityDetail"></span></div>
          <div class="c-gauge" aria-hidden="true"><div class="c-band" id="cBand"></div><div class="c-marker" id="cMarker"></div></div>
          <div class="c-ticks" id="cTicks" aria-hidden="true"></div>
          <p class="stability-message" id="stabilityMessage"></p>
        </div>
        <details class="key">
          <summary>Understanding the colours<span data-preset-only> and χ²</span></summary>
          <div class="key-content">
            <ul><li data-preset-only><i class="sw dot" style="background:#f4f7fb"></i><b>White points</b> are measured speeds of <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> (SPARC); bars show uncertainty.</li><li><i class="sw ldash" style="border-color:#5fd0e0"></i><b>Cyan dashed</b> is the speed expected from visible matter only (stars + gas).</li><li><i class="sw dash" style="border-color:#b39ae8"></i><b>Purple</b> is the extra speed supplied by the dark matter halo (NFW mode).</li><li><i class="sw" style="background:#ffae5c"></i><b>Orange</b> is the combined model (or the MOND prediction in MOND mode).<span data-preset-only> The closer it hugs the white points, the better the fit.</span></li><li data-preset-only><i class="sw bar" style="background:#5fd0e0"></i><i class="sw bar" style="background:#ff8a80"></i><b>Residual bars:</b> teal = model is too slow, coral = model is too fast, no bar = perfect match.</li><li data-preset-only><b>χ² badge:</b> <span class="g">teal</span> good fit, <span class="o">orange</span> fair, <span class="r">red</span> poor.</li><li id="componentKey" hidden><b>Thin dotted lines</b> show your custom galaxy's stellar disk, bulge and gas contributions; they add in quadrature to the cyan baryonic curve.</li></ul>
            <div data-preset-only>
            <h3>Why isn't χ² exactly 0?</h3>
            <p>χ² adds up how far the model is from each measured point, in units of that point's uncertainty. A score of 0 would mean every point matched perfectly, which real data never allows: measurements carry noise, so even a correct model scatters around them by roughly one error bar per point (χ²ᵣ ≈ 1 is a good fit).</p>
            <p>Models also simplify reality. NFW assumes a smooth spherical halo, MOND uses a characteristic acceleration a<sub>0</sub>, and both rely on an estimated stellar mass for the visible matter. Their small mismatches add to the score, which is why neither theory reaches 0.</p>
            </div>
          </div>
        </details>
      </article>
    </section>
    <div class="section-head"><span class="kicker">03</span><h2>Background</h2><p>The two explanations this simulator tests.</p></div>
    <section class="learn" id="theoryGuide" aria-label="Theory guide">
      <details class="nfw"><summary>What is a dark matter halo?</summary>
        <div class="learn-content">
        <p class="tldr">Dark matter is matter we cannot see; we only notice it through its gravity. A halo is a huge cloud of it surrounding a galaxy.</p>
        <figure class="schematic">
          <svg viewBox="0 0 320 150" role="img" aria-label="Schematic: the speed expected from visible matter falls with radius, but measured speeds stay roughly flat">
            <path d="M36 14V124H308" fill="none" stroke="currentColor" stroke-opacity=".35"/>
            <path d="M36 124C58 52 78 38 108 42S196 78 302 98" fill="none" stroke="#5fd0e0" stroke-width="2" stroke-dasharray="5 4"/>
            <g fill="#f4f7fb"><circle cx="58" cy="70" r="3"/><circle cx="84" cy="49" r="3"/><circle cx="110" cy="42" r="3"/><circle cx="140" cy="40" r="3"/><circle cx="170" cy="39" r="3"/><circle cx="200" cy="40" r="3"/><circle cx="230" cy="38" r="3"/><circle cx="260" cy="40" r="3"/><circle cx="290" cy="39" r="3"/></g>
            <text x="300" y="30" text-anchor="end" fill="#f4f7fb">measured: stays flat</text>
            <text x="300" y="114" text-anchor="end" fill="#5fd0e0">visible matter only: falls</text>
            <text x="172" y="142" text-anchor="middle">radius →</text>
            <text x="14" y="70" text-anchor="middle" transform="rotate(-90 14 70)">speed →</text>
          </svg>
          <figcaption>Schematic only, not real data. Far from most of the visible mass, Newtonian gravity predicts speeds falling roughly as 1/√r, yet measured curves stay roughly flat (Rubin, Ford &amp; Thonnard 1980).</figcaption>
        </figure>
        <p>The outer stars in <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span> move too fast for the visible matter alone to hold them in orbit. A halo adds the extra gravity needed to keep them there.</p>
        <p>This simulator uses an NFW halo (Navarro, Frenk &amp; White 1996): packed tightly in the middle and thinning out further away. Two sliders control it: the halo's mass (M<sub>200</sub>) and its scale radius (r<sub>s</sub>), which sets how spread out it is.</p>
        </div>
      </details>
      <details class="mond"><summary>What is MOND?</summary>
        <div class="learn-content">
        <p class="tldr">MOND (Modified Newtonian Dynamics) says there is no dark matter. Instead, gravity behaves a little differently when it is extremely weak.</p>
        <p>Near a galaxy's centre gravity is strong, and MOND gives the same answer as Newton. Out in the edges gravity becomes tiny, weaker than a threshold called a<sub>0</sub> (1.2×10⁻¹⁰ m/s² here). Below that, MOND makes gravity stronger than Newton predicts, so the outer stars keep moving fast.</p>
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
          <figcaption>Schematic only, not real data. Inside the dotted line gravity is strong and MOND matches Newton. Outside it, gravity from visible matter drops below a₀ and MOND boosts it, so the speed levels off instead of falling (Milgrom 1983).</figcaption>
        </figure>
        <p>The rule in symbols:</p>
        <div class="eq"><math display="block"><mi>μ</mi><mo>(</mo><mfrac><mi>g</mi><msub><mi>a</mi><mn>0</mn></msub></mfrac><mo>)</mo><mi>g</mi><mo>=</mo><msub><mi>g</mi><mi>N</mi></msub></math></div>
        <p>Here g is the real pull of gravity and g<sub>N</sub> is what Newton predicts from the visible matter. μ acts like a dial: it is about 1 when gravity is strong (plain Newton) and drops below 1 when gravity is weak, which boosts g. This simulator uses (Famaey &amp; McGaugh 2012):</p>
        <div class="eq"><math display="block"><mi>μ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mfrac><mi>x</mi><msqrt><mn>1</mn><mo>+</mo><msup><mi>x</mi><mn>2</mn></msup></msqrt></mfrac></math></div>
        <p>Far out, the speed stops falling and settles at a fixed value, where M is the visible mass:</p>
        <div class="eq"><math display="block"><msup><mi>v</mi><mn>4</mn></msup><mo>=</mo><mi>G</mi><mi>M</mi><msub><mi>a</mi><mn>0</mn></msub></math></div>
        <p>So the flat outer speed comes from visible matter alone. The only things to tune are a<sub>0</sub> and the mass-to-light ratio, which is why the <b>MOND best fit</b> switch changes those two.</p>
        <p>Switch between the two models and compare each curve with the white measured points. One galaxy alone cannot settle which idea is right.</p>
        </div>
      </details>
    </section>
    <div class="section-head" id="method"><span class="kicker">04</span><h2>Method &amp; data</h2><p>What the simulator calculates, and where the numbers come from.</p></div>
    <footer>
      <dl class="method-grid">
        <div data-preset-only><dt>Data</dt><dd>Measured speeds V<sub>obs</sub> ± σ and the gas, disk and bulge speeds at each radius, from SPARC (Lelli, McGaugh &amp; Schombert 2016), read from the project’s table2.dat.</dd></div>
        <div><dt>Visible matter</dt><dd><span id="footerModel"><span class="eq"><math display="block"><msubsup><mi>V</mi><mi>bar</mi><mn>2</mn></msubsup><mo>=</mo><msub><mi>V</mi><mi>gas</mi></msub><mo>|</mo><msub><mi>V</mi><mi>gas</mi></msub><mo>|</mo><mo>+</mo><mn id="ftDisk">0.50</mn><msubsup><mi>V</mi><mi>disk</mi><mn>2</mn></msubsup><mo>+</mo><mn id="ftBulge">0.70</mn><msubsup><mi>V</mi><mi>bulge</mi><mn>2</mn></msubsup></math></span>The numbers are the disk and bulge M/L ratios.</span><span id="footerCustom" hidden><span class="eq"><math display="block"><msubsup><mi>V</mi><mi>bar</mi><mn>2</mn></msubsup><mo>=</mo><msubsup><mi>V</mi><mi>disk</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>bulge</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>gas</mi><mn>2</mn></msubsup></math></span>From exponential disks and a Hernquist bulge.</span></dd></div>
        <div><dt>Dark matter</dt><dd>NFW halo (Navarro, Frenk &amp; White 1996), set by M<sub>200</sub> and r<sub>s</sub>. Model speed:<span class="eq"><math display="block"><msubsup><mi>V</mi><mi>model</mi><mn>2</mn></msubsup><mo>=</mo><msubsup><mi>V</mi><mi>bar</mi><mn>2</mn></msubsup><mo>+</mo><msubsup><mi>V</mi><mi>halo</mi><mn>2</mn></msubsup></math></span></dd></div>
        <div><dt>MOND</dt><dd>No halo (Milgrom 1983). a<sub>0</sub> = 1.2×10⁻¹⁰ m/s².<span class="eq"><math display="block"><mi>μ</mi><mo>(</mo><mfrac><mi>g</mi><msub><mi>a</mi><mn>0</mn></msub></mfrac><mo>)</mo><mi>g</mi><mo>=</mo><msub><mi>g</mi><mi>N</mi></msub><mo>,</mo><mspace width="1em"/><mi>μ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mfrac><mi>x</mi><msqrt><mn>1</mn><mo>+</mo><msup><mi>x</mi><mn>2</mn></msup></msqrt></mfrac></math></span></dd></div>
        <div><dt>Simulation</dt><dd>Massless test stars (800 for the presets) moved by Velocity-Verlet steps, 2 per Myr, over 1 Gyr.</dd></div>
        <div data-preset-only><dt>Goodness of fit</dt><dd><span class="eq"><math display="block"><msup><mi>χ</mi><mn>2</mn></msup><mo>=</mo><mo>∑</mo><msup><mrow><mo>(</mo><mfrac><mrow><msub><mi>V</mi><mi>obs</mi></msub><mo>−</mo><msub><mi>V</mi><mi>model</mi></msub></mrow><mi>σ</mi></mfrac><mo>)</mo></mrow><mn>2</mn></msup><mo>,</mo><mspace width="1em"/><msubsup><mi>χ</mi><mi>r</mi><mn>2</mn></msubsup><mo>=</mo><mfrac><msup><mi>χ</mi><mn>2</mn></msup><mrow><mi>N</mi><mo>−</mo><mn>2</mn></mrow></mfrac></math></span>for N measured points.</dd></div>
      </dl>
      <button type="button" id="exportCsv" class="export-btn" title="Download the time-series log and the rotation-curve model as two CSV files">Export raw data (CSV)</button>
    </footer>
  </main>
  <dialog class="tutorial-welcome" id="tutorialWelcome" aria-labelledby="welcomeTitle" aria-describedby="welcomeDescription">
    <div class="tutorial-eyebrow">Explore <span data-galaxy-name>__DEFAULT_GALAXY_NAME__</span></div>
    <h2 id="welcomeTitle">Welcome to the simulator</h2>
    <p id="welcomeDescription">New here? Take a short, hands-on tour. You will find the puzzle hidden in a real galaxy’s rotation curve, then test two explanations for it: dark matter and modified gravity. You can skip the tour at any time. The Custom Sandbox, where you build your own galaxy, has a tour of its own.</p>
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
    <div class="tutorial-actions"><span class="tutorial-keys">← → to navigate</span><button type="button" id="tutorialBack">Back</button><button type="button" id="tutorialMore" hidden>Sandbox tour</button><button type="button" class="tutorial-primary" id="tutorialNext">Next</button></div>
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
        // The Custom Sandbox galaxy is imaginary, so it has no measured speeds: its curves are drawn on an even
        // grid of 60 radii out to 1.25 x R_max instead of at a SPARC galaxy's data radii.
        function loadCustomGrid(rMax){
            const n=60,rEnd=1.25*rMax;
            rotationRows=Array.from({length:n},(_,i)=>({radius_kpc:rEnd*(i+1)/n,observed_kms:0,error_kms:0,gas_kms:0,gas_sq:0,disk_kms:0,bulge_kms:0}));
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
        // Bloom (visual choice, like a camera exposure): bright light bleeds softly into its surroundings, then a soft
        // highlight shoulder stops piled-up light clipping to flat white. Falls back to a plain render if the add-ons fail to load.
        const composer=(()=>{
            if(!THREE.EffectComposer||!THREE.UnrealBloomPass)return null;
            try{
                const gl2=renderer.capabilities&&renderer.capabilities.isWebGL2;
                const opts={type:THREE.HalfFloatType,format:THREE.RGBAFormat};
                const target=gl2&&THREE.WebGLMultisampleRenderTarget?new THREE.WebGLMultisampleRenderTarget(1,1,opts):new THREE.WebGLRenderTarget(1,1,opts);
                if(gl2)target.samples=4;
                const fx=new THREE.EffectComposer(renderer,target);
                fx.addPass(new THREE.RenderPass(scene,camera));
                fx.addPass(new THREE.UnrealBloomPass(new THREE.Vector2(256,256),0.55,0.45,0.62));
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
        window.simTime=simTime;
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
            // The sandbox galaxy has no measurements, so its file leaves out the observed and residual columns.
            let csv=customActive?"Radius_kpc,V_gas,V_disk,V_bar,V_nfw_halo,V_total_model\n":"Radius_kpc,V_obs,V_obs_err,V_gas,V_disk,V_bar,V_nfw_halo,V_total_model,Residual_V_obs_minus_V_bar\n";
            for(let r=Math.max(1,Math.ceil(radii[0]-1e-9));r<=radii[radii.length-1]+1e-9;r++){
                const vo=interpolateAt(vObs,r),eo=interpolateAt(err,r),vg=interpolateAt(gas,r),vd=interpolateAt(disk,r);
                const vb=Math.sqrt(Math.max(0,interpolateAt(gasSq,r)+vd*vd)),vh=mondMode?Math.sqrt(Math.max(mondVelocity(vb,r)**2-vb*vb,0)):halo(r);
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
        // Ejection line: the red ring at PM.ejection_factor x R_max. Off by default; the tutorial can also show it.
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
        function resizeOrbit(){
            if (!nativeFsElement()) {
                containerOriginalWidth = orbitPlot.clientWidth;
                containerOriginalHeight = orbitPlot.clientHeight;
            }
            const w=Math.max(orbitPlot.clientWidth,1),h=Math.max(orbitPlot.clientHeight,1),dpr=Math.min(window.devicePixelRatio||1,2);
            renderer.setPixelRatio(dpr);renderer.setSize(w,h,false);
            if(composer){composer.setPixelRatio(dpr);composer.setSize(w,h);}
            starUniforms.uDpr.value=dpr;
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
            if(needsRender){needsRender=false;haloUniforms.uScale.value=starUniforms.uScale.value=renderer.domElement.height/(2*Math.tan(camera.fov*Math.PI/360));
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
    function updateChi(total){const b=document.getElementById('chiBadge');b.hidden=customActive;document.getElementById('tmChiRow').hidden=customActive;if(customActive)return;const chi=obsV.reduce((s,v,i)=>s+((v-total[i])/obsErr[i])**2,0),dof=Math.max(1,obsV.length-2),red=chi/dof;b.textContent='χ² '+chi.toFixed(0)+' · χ²ᵣ '+red.toFixed(2);b.className=red<2?'good':red<10?'ok':'poor';const fc=document.getElementById('tm-chi');fc.textContent=red.toFixed(2);fc.className=b.className;}
    const residualsOf=total=>obsV.map((v,i)=>v-total[i]);
    const AXIS={gridcolor:'rgba(148,170,200,0.08)',gridwidth:1,zerolinecolor:'rgba(148,170,200,0.3)',linecolor:'rgba(148,170,200,0.32)',showline:true,mirror:false,tickcolor:'rgba(148,170,200,0.32)',ticks:'outside',ticklen:4,automargin:true,tickfont:{size:11,color:'#8e9db3',family:'JetBrains Mono,ui-monospace,Consolas,monospace'},title:{font:{size:12,color:'#b8c4d6'}}};
    // Residual bars are about half the mean spacing of the SPARC radii, so dense and sparse curves both read clearly.
    const residualBarWidth=()=>Math.round(53*(radii[radii.length-1]-radii[0])/(radii.length-1))/100;
    // The custom galaxy's stellar disk, bulge and gas speeds are drawn as thin dotted lines (hidden for presets).
    const COMPONENT_STYLE=[['Stellar disk','#ffd166'],['Bulge','#ff8fab'],['Gas','#7bd88f']];
    let componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
    const componentVisible=i=>customActive&&(i!==1||currentAnatomy.bulgeMass>0);
    const buildCurveData=(halo,total)=>[
      {type:'scatter',mode:'markers',name:'Observed (SPARC)',visible:!customActive,x:radii,y:obsV,error_y:{type:'data',array:obsErr,visible:true,color:'rgba(232,238,246,0.5)',thickness:1.1,width:3},marker:{size:6.5,color:'#f4f7fb',line:{color:'#070b14',width:1.4}},customdata:obsErr,hovertemplate:'Observed <b>%{y:.0f}</b> ± %{customdata:.0f} km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Visible matter (baryons)',x:radii,y:vBary.slice(),line:{color:'#5fd0e0',width:2.2,dash:'6px,4px',shape:'spline',smoothing:0.6},hovertemplate:'Baryons <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Dark matter halo',x:radii,y:halo,fill:'tozeroy',fillcolor:'rgba(179,154,232,0.08)',line:{color:'#b39ae8',width:1.8,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:'Halo <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Total model',x:radii,y:total,line:{color:'#ffae5c',width:3,shape:'spline',smoothing:0.6},hovertemplate:'Model <b>%{y:.0f}</b> km/s<extra></extra>'},
      {type:'bar',name:'Residual',visible:!customActive,x:radii,y:residualsOf(total),yaxis:'y2',showlegend:false,width:residualBarWidth(),marker:{color:residualsOf(total).map(r=>r>=0?'rgba(95,208,224,0.6)':'rgba(255,138,128,0.6)')},hovertemplate:'Residual <b>%{y:+.0f}</b><extra></extra>'},
      ...COMPONENT_STYLE.map(([name,color],i)=>({type:'scatter',mode:'lines',name,x:radii,y:componentCurves[i].slice(),visible:componentVisible(i),line:{color,width:1.6,dash:'dot',shape:'spline',smoothing:0.6},hovertemplate:name+' <b>%{y:.0f}</b><extra></extra>'}))
    ];
    const curveLayout={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',margin:{l:66,r:20,t:44,b:54},font:{color:textColor,size:12,family:'Inter,Segoe UI,system-ui,sans-serif'},
      xaxis:{...AXIS,fixedrange:true,title:{text:'Radius, R (kpc)',standoff:10,font:{size:12,color:'#b8c4d6'}},anchor:'y2',showgrid:true,zeroline:false,spikemode:'across',spikecolor:'rgba(200,215,230,0.35)',spikethickness:1,spikedash:'dot',spikesnap:'cursor',hoverformat:'.1f'},
      yaxis:{...AXIS,title:{text:'Orbital speed, V (km/s)',standoff:10,font:{size:12,color:'#b8c4d6'}},range:[0,190],fixedrange:true,domain:[0.32,1],zeroline:false},
      yaxis2:{...AXIS,title:{text:'Obs − model',standoff:6,font:{size:11,color:'#b8c4d6'}},domain:[0,0.24],anchor:'x',zeroline:true,zerolinewidth:1.5,nticks:4,fixedrange:true},
      bargap:0.2,legend:{orientation:'h',yanchor:'bottom',y:1.02,xanchor:'left',x:0,font:{size:11.5,color:'#cdd6e3'},itemsizing:'constant',itemwidth:36,tracegroupgap:4},
      hovermode:'x unified',hoverlabel:{bgcolor:'rgba(10,15,26,0.94)',bordercolor:'rgba(148,170,200,0.3)',font:{color:'#eef2f8',size:11,family:'Inter,Segoe UI,sans-serif'}},dragmode:false};
    const curveConfig={responsive:true,displaylogo:false,displayModeBar:false,scrollZoom:false,doubleClick:false};
    // The sandbox has no residuals, so its velocity axis fills the whole plot.
    const curveLayoutNow=()=>customActive?{...curveLayout,xaxis:{...curveLayout.xaxis,anchor:'y'},yaxis:{...curveLayout.yaxis,domain:[0,1]},yaxis2:{...curveLayout.yaxis2,visible:false}}:curveLayout;
    Plotly.newPlot('curvePlot',buildCurveData(initialHalo,totalCurve),curveLayoutNow(),curveConfig);
    updateChi(totalCurve);

    const massSlider=document.getElementById('massSlider'),scaleSlider=document.getElementById('scaleSlider');
    let sliderUpdateTimer=null,sliderDragging=false;
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
      stabilityMonitor.hidden=!custom;
      document.getElementById('presetPicker').hidden=custom; // nothing real to compare an imaginary galaxy with
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
      const residualColors=res.map(r=>r>=0?'rgba(95,208,224,0.6)':'rgba(255,138,128,0.6)');
      Plotly.update('curvePlot',{y:[vBary.slice(),halo,total,res,...componentCurves.map(c=>c.slice())],name:['Visible matter (baryons)',mondMode?'Dark matter halo (none)':'Dark matter halo',mondMode?'MOND model':'Total model','Residual',...COMPONENT_STYLE.map(([name])=>name)],visible:[true,!mondMode,true,!customActive,...COMPONENT_STYLE.map((_,i)=>componentVisible(i))],'marker.color':[undefined,undefined,undefined,residualColors,undefined,undefined,undefined]},{'yaxis.range':[0,yMax]},[1,2,3,4,5,6,7]);
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
      const mlStat=document.getElementById('mlStat');
      mlStat.classList.toggle('off',customActive); // hidden in the sandbox, which sets stellar masses directly
      document.getElementById('mlDisk').textContent=mlState.disk.toFixed(2);
      document.getElementById('mlBulge').textContent=mlState.bulge.toFixed(2);
      document.getElementById('ftDisk').textContent=mlState.disk.toFixed(2);
      document.getElementById('ftBulge').textContent=mlState.bulge.toFixed(2);
      const e=Math.floor(Math.log10(MOND_A0));
      document.getElementById('mondA0Text').textContent=`${(MOND_A0/10**e).toFixed(2)}×10${String(e).replace('-','⁻').replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`;
    }
    mondFitToggle.addEventListener('change',()=>{applyBaryonState();updateHalo();});
    // "?" on the M/L tile: a mouse hover fades the explainer in and moving away fades it out
    // (a short grace period lets the pointer travel into the box). Touch and keyboard users tap/press
    // the button instead, and close it with ×, Escape or a tap elsewhere.
    const mlInfoBtn=document.getElementById('mlInfoBtn'),mlExplainer=document.getElementById('mlExplainer');
    let mlHover=false,mlPinned=false,mlCloseTimer=null,mlLastPointer='';
    const mlIsOpen=()=>mlExplainer.classList.contains('open');
    const renderMlExplainer=()=>{
      const open=mlHover||mlPinned;
      if(open&&!mlIsOpen()){
        // Right-aligned under the tile by default; flip to left-aligned if that would run off a narrow screen.
        mlExplainer.style.left='';mlExplainer.style.right='';
        if(mlExplainer.getBoundingClientRect().left<16){mlExplainer.style.left='0';mlExplainer.style.right='auto';}
      }
      mlExplainer.classList.toggle('open',open);
      mlExplainer.setAttribute('aria-hidden',String(!open));
      mlInfoBtn.setAttribute('aria-expanded',String(open));
    };
    const closeMlExplainer=()=>{mlHover=false;mlPinned=false;clearTimeout(mlCloseTimer);renderMlExplainer();};
    [mlInfoBtn,mlExplainer].forEach(el=>{
      el.addEventListener('pointerenter',e=>{if(e.pointerType!=='mouse')return;clearTimeout(mlCloseTimer);mlHover=true;renderMlExplainer();});
      el.addEventListener('pointerleave',e=>{if(e.pointerType!=='mouse')return;clearTimeout(mlCloseTimer);mlCloseTimer=setTimeout(()=>{mlHover=false;renderMlExplainer();},150);});
    });
    mlInfoBtn.addEventListener('pointerdown',e=>{mlLastPointer=e.pointerType;});
    mlInfoBtn.addEventListener('click',e=>{
      e.stopPropagation();
      if(mlLastPointer==='mouse'){mlLastPointer='';return;}  // mouse users already see it on hover
      mlLastPointer='';mlPinned=!mlPinned;renderMlExplainer();
    });
    document.getElementById('mlExplainerClose').addEventListener('click',()=>{closeMlExplainer();mlInfoBtn.focus();});
    document.addEventListener('click',e=>{if(mlIsOpen()&&!mlExplainer.contains(e.target)&&e.target!==mlInfoBtn)closeMlExplainer();});
    document.addEventListener('keydown',e=>{if(e.key==='Escape'&&mlIsOpen()){closeMlExplainer();mlInfoBtn.focus();}});
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
    // Capture phase: the fullscreen panel stops pointer events bubbling (so drags don't spin the camera), and a
    // bubbling listener here never heard its sliders being released, so the orbits were never rebuilt.
    window.addEventListener('pointerup',finishSliderDrag,true);
    window.addEventListener('pointercancel',finishSliderDrag,true);
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
       The Custom Sandbox hides the dropdown: its imaginary galaxy has no measurements to compare with, and the
       anatomy sliders drive every output through PhysicsEngine. */
    const galaxySelect=document.getElementById('galaxySelect'),galaxySummary=document.getElementById('galaxySummary');
    Object.values(GALAXIES).forEach(galaxy=>galaxySelect.add(new Option(galaxy.name,galaxy.id)));
    function applyGalaxyLabels(galaxy){
      document.querySelectorAll('[data-galaxy-name]').forEach(el=>{el.textContent=galaxy.name;});
      galaxySummary.textContent=`${galaxy.summary} · ${galaxy.rotation.radius_kpc.length} measured speeds (SPARC)`;
      orbitPlot.setAttribute('aria-label',`Interactive three-dimensional ${galaxyMode==='custom'?'custom galaxy':galaxy.name} star trajectories`);
      document.getElementById('curvePlot').setAttribute('aria-label',galaxyMode==='custom'?'Custom galaxy baryonic, halo and model rotation curves':`${galaxy.name} observed, baryonic, and fitted rotation curves`);
      document.getElementById('viewReset').title=`Return to the default ${galaxyMode==='custom'?'custom galaxy':galaxy.name} view`;
    }
    galaxySelect.value=activeGalaxy.id;
    applyGalaxyLabels(activeGalaxy);
    // Wording that differs between a preset galaxy and the sandbox.
    function applyModeLabels(){
      const custom=customActive,fitButton=document.getElementById('applyObservedFit');
      document.querySelectorAll('[data-preset-only]').forEach(el=>{el.hidden=custom;});
      document.getElementById('massHelp').textContent=custom
        ?'How much dark matter surrounds your galaxy. Starts at the typical halo for these stars (Moster et al. 2013); runs from 0 to 5 times that.'
        :'How much dark matter surrounds the galaxy. More mass means stronger gravity, so stars orbit faster and the orange curve rises. Set it to 0 for visible matter only.';
      document.getElementById('scaleHelp').textContent=custom
        ?'How spread out the halo is. Small packs the dark matter near the centre. Starts at the typical value for this halo mass (Dutton & Macciò 2014).'
        :'How spread out the halo is, for the same total mass. A small radius packs the dark matter near the centre and speeds up the inner stars; a large radius spreads it out.';
      fitButton.textContent=custom?'Typical Halo':'Apply Observed Fit';
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
      mondBest=null; // a cached MOND fit belongs to the previous galaxy
      loadGalaxyData(galaxy);
      applyBaryonState();
      configureHaloSliders();
      applyModeLabels();
      componentCurves=COMPONENT_STYLE.map(()=>radii.map(()=>0));
      referenceSpeedSquared=presetReferenceSpeedSquared;
      const anatomy=presetAnatomy(galaxy);
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),PRESET_TEMPLATE_MASS));
      syncStarMeshes();updateStarTelemetry();controls.reset();
      // Restart the orbit clock and integrate the new population in the new galaxy's potential.
      simTime=0;window.simTime=0;syncTimelineUI(0);
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
    function rebuildCustomGalaxy({replot=false}={}){
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
      setPopulation(anatomy,starCountFor(visibleMassOf(anatomy),sandboxTemplateMass()));
      syncStarMeshes();updateStarTelemetry();updateCustomSummary(baseline);
      if(replot){
        const{halo,total}=modelCurves(Number(massSlider.value),Number(scaleSlider.value));
        Plotly.react('curvePlot',buildCurveData(halo,total),curveLayoutNow(),curveConfig);
      }
      refreshOrbits();
      updateHalo();clearTimeout(orbitModelTimer);
    }
    function activateCustom(){
      customActive=true;haloFine=1;rsManual=false;mondBest=null;
      applyModeLabels();
      simTime=0;window.simTime=0;syncTimelineUI(0);
      rebuildCustomGalaxy({replot:true});
      controls.reset();
      document.dispatchEvent(new CustomEvent('galaxychange',{detail:{galaxy:activeGalaxy}}));
    }
    let customUpdateTimer=null;
    document.addEventListener('sandboxchange',()=>{
      if(!customActive)return;
      haloFine=1;rsManual=false; // a new galaxy gets its own typical halo: earlier halo tuning is dropped
      if(customUpdateTimer!==null)return;
      customUpdateTimer=setTimeout(()=>{customUpdateTimer=null;rebuildCustomGalaxy();},40);
    });
    document.getElementById('sandboxReset').addEventListener('click',()=>{haloFine=1;rsManual=false;});
    function selectGalaxy(id){
      const galaxy=GALAXIES[id];
      if(!galaxy||galaxy===activeGalaxy||customActive)return;
      loadPreset(galaxy);
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
        document.getElementById('fsTimeValue').textContent=`${val} / 1000 Myr`;
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
    function lapsAt(idx,timeMyr){
      const t=activeTrajectories[idx],last=Math.min(frameCount-1,Math.floor(Math.max(0,timeMyr)/SIM_DT_MYR));
      let angle=0;
      for(let f=1;f<=last;f++){
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
    const lapText=(name,idx)=>`${name} · ${lapsAt(idx,simTime).toFixed(1)} laps · ${starSpeedAt(idx,simTime).toFixed(0)} km/s`;
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
    function updateStarMarks(){
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
        el.classList.toggle('flip',markVec.x>0.4); // keep the label on screen near the right edge
        el.hidden=false;
        el.style.transform=`translate(${((markVec.x+1)/2*orbitPlot.clientWidth).toFixed(1)}px,${((1-markVec.y)/2*orbitPlot.clientHeight).toFixed(1)}px)`;
        el.classList.toggle('warm',!!mark.warm);el.classList.toggle('track',!!mark.track);
        const text=mark.label(idx);
        if(el.firstChild.textContent!==text)el.firstChild.textContent=text;
      });
      tutorialWatchTick();
    }

    /* ---------- Guided tours ----------
       Two tours share one card: 'main' (the rotation-curve puzzle, the halo, MOND and the 3D view) and 'sandbox'
       (building a custom galaxy). Steps are looked up by id, never by position, so they can be added or reordered.
       A step can need NFW mode (nfw), a halo mass above zero (halo) or one galaxy source (mode). When the page is
       not in that state, the card says how to get there and the pointer moves to the control that does it.
       The tours never change the simulator's settings themselves. */
    const TOURS={
      main:{label:'Tutorial',steps:[
        {id:'puzzle',target:'curvePlot',gesture:'scan',title:'Spot the puzzle',
         text:'This graph shows how fast things orbit at each distance from the galaxy’s centre. Each white point is a measured speed. The cyan line is the speed you would expect if the stars and gas we can see were all the matter there is.',
         action:'Move your mouse (or tap) along the right of the graph. Do the white points fall with the cyan line, or stay above it?'},
        {id:'removeHalo',nfw:true,target:'massSlider',gesture:'zero',title:'Take the dark matter away',
         text:'That gap is the puzzle: the outer galaxy spins faster than its visible matter can explain. The orange line is the model’s prediction, visible matter plus a dark matter halo. The dashed purple line is the halo’s share.',
         action:'Drag halo mass all the way left to 0. The orange line drops onto the cyan one, because only visible matter is left.'},
        {id:'addHalo',nfw:true,target:'massSlider',gesture:'increase',title:'Add the halo back',
         text:'A dark matter halo is a huge cloud of invisible matter around the galaxy. Its gravity makes the stars orbit faster. This slider sets its mass.',
         action:'Drag halo mass to the right and watch the orange line rise. Can you make it run through the white points?'},
        {id:'haloShape',nfw:true,target:'scaleSlider',gesture:'adjust',title:'Pack the halo in or spread it out',
         text:'The scale radius changes the halo’s shape without changing its total mass. A small value packs the dark matter close to the centre. A large value spreads it out.',
         action:'Move the scale radius both ways. A small value lifts the inner part of the orange line; a large value lowers it.'},
        {id:'haloNumbers',nfw:true,halo:true,cardRight:true,target:'haloMetrics',gesture:'click',title:'Read the halo’s numbers',
         text:'These cards describe the halo you have built. r₂₀₀ is the halo’s edge by convention: the radius inside which it averages 200 times the critical density of the universe. It reaches far beyond the stars. Concentration, c, is r₂₀₀ divided by the scale radius: a high c means the dark matter is packed towards the centre.',
         action:'Press each card to read what it means.'},
        {id:'fit',nfw:true,target:'applyObservedFit',gesture:'click',title:'Score your fit',
         text:'χ² (chi-squared), in the badge above the graph, scores the fit: it adds up the gaps between the orange line and the white points. Lower is better. The residual bars under the graph show each gap.',
         action:'Try to lower χ² with both halo sliders. Apply Observed Fit jumps back to the starting fit.'},
        {id:'mond',target:'modelBar',gesture:'switch',title:'Test a rival idea',
         text:'MOND has no dark matter at all. Instead, gravity is stronger than Newton’s law predicts where it is very weak, as in the outer galaxy. The halo sliders switch off.',
         action:'Switch to MOND. The orange line is now MOND’s prediction. Compare its χ² with the halo’s.'},
        {id:'mondFit',mond:true,mode:'preset',target:'mondFitWrap',gesture:'click',title:'Give MOND its best shot',
         text:'MOND has one tuning number, a₀, and the visible matter depends on the mass-to-light ratio, which is only an estimate. MOND best fit adjusts both to get the lowest χ².',
         action:'Turn on MOND best fit, at the top right of the page, and see how far χ² drops. Then switch the gravity model back to NFW.'},
        // Watching the simulation: these steps run inside fullscreen (fs), where the tour card follows.
        {id:'fsEnter',nfw:true,viewPanel:true,target:'toggleFullscreen',gesture:'click',title:'Now watch it happen',
         text:'The graph shows what the models predict. The 3D view tests it: hundreds of test stars, each pushed by the model’s gravity, step by step, for a billion years. The best details only show up in fullscreen.',
         action:'Press Fullscreen. The tour comes with you, and you can drag to spin the galaxy or scroll to zoom at any time.'},
        {id:'fsClean',fs:true,viewPanel:true,target:'toggleGas',gesture:'click',title:'Clear the stage',
         text:'Two layers sit over the stars. The blue haze is the cold gas, drawn as a still backdrop to show where it lies. The purple cloud is the dark matter halo, drawn so you can see it, though real dark matter is invisible.',
         action:'Under Show, switch off Gas and Halo, so only the test stars are left.',
         done:()=>!gasCloud.visible&&(haloHidden||mondMode||Number(massSlider.value)===0),doneText:'Clear view: everything left is a test star.'},
        {id:'fsOrbits',fs:true,viewPanel:true,target:'viewFaceOn',gesture:'click',title:'Look down on the disk',
         text:'Face-on, you look straight down on the galaxy. With Orbits on, each line is one star’s path over the whole billion years. The disk stars trace near-perfect circles: they start on circular orbits and the halo’s gravity keeps them there. The glow and black hole at the centre are decoration only.',
         action:'Press Face-on, then switch on Orbits.',
         done:()=>orbitLines.visible&&controls.view.phi<0.3,doneText:'Circles everywhere: a galaxy in balance.'},
        {id:'fsRace',fs:true,nfw:true,marks:['race'],target:'fsPlay',gesture:'click',title:'The inner star laps the outer one',
         text:'Two stars are ringed, one near the centre and one far out, with live lap counts and speeds. Their speeds are similar (that is the flat rotation curve), but the outer star has a far longer track, so it falls behind lap after lap. This is differential rotation: a galaxy does not turn like a solid wheel.',
         action:'Press Play (or Space) and watch the two labels.',
         done:()=>tutorialPlayedMyr>=400,doneText:'Similar speeds, very different lap counts.'},
        {id:'fsTilt',fs:true,viewPanel:true,marks:['tilt'],target:'viewEdgeOn',gesture:'click',title:'Turn it on its side',
         text:'Edge-on, the disk is a thin line: a typical star sits only about 0.1 kpc above or below the middle, while the disk is tens of kpc across. The ringed star’s orbit is tilted very slightly, so it rises above the middle and dips below it once every lap.',
         action:'Press Edge-on and keep it playing. Watch the height in the label change sign.',
         done:()=>Math.abs(controls.view.phi-Math.PI/2)<0.2,doneText:'Up, down, up: one bob per lap.'},
        {id:'fsEscape',fs:true,nfw:true,marks:['escaper'],ring:true,target:'fsMassSlider',gesture:'zero',title:'Pull the halo away',
         text:'The view pulls back to show the faint red circle: the ejection line, 1.5 times the galaxy’s visible edge. With no halo, the outer stars move too fast for the visible matter to hold. The ringed star is the first to cross the line, and the Ejected count in the top-left corner climbs as others follow. Ejected stars do not vanish; they keep coasting outwards.',
         action:'Drag Halo mass to 0, then press Restart.',
         done:prev=>{const k=firstEscaper();return Number(massSlider.value)===0&&k>=0&&prev<escapeTimeByStar[k]&&simTime>=escapeTimeByStar[k];},doneText:'Gone: visible matter alone cannot hold the outer stars.'},
        {id:'fsFullHalo',fs:true,nfw:true,halo:true,viewPanel:true,target:'viewFullHalo',gesture:'click',title:'Put the halo back, and step back',
         text:'With the halo back, press Restart and the same stars stay put. Now zoom right out: the bright speck in the middle is the whole visible galaxy, and the sphere is the halo’s edge, r₂₀₀. Real halos have no sharp edge; the sphere just marks it.',
         action:'Raise Halo mass, then press Full halo. Reset brings you back.'},
        {id:'galaxies',mode:'preset',target:'galaxySelect',gesture:'click',title:'Try another galaxy',
         text:'The Galaxy menu loads other real galaxies from the SPARC survey, each with its own measured speeds and starting halo.',
         action:'Choose a different galaxy. Does it show the same gap between the white points and the cyan line?'},
        {id:'ml',mode:'preset',target:'mlInfoBtn',gesture:'click',title:'How heavy are the stars?',
         text:'Telescopes measure light, not mass. The mass-to-light ratio converts one into the other and sets how high the cyan line sits. The values here are standard averages, so they are an assumption.',
         action:'Hover over or tap the ? beside Mass-to-light at the top of the page to read more.'},
        {id:'finish',target:'theoryGuide',gesture:'guide',title:'What this does and doesn’t show',
         text:'A close fit here does not prove either idea. Both models simplify a real galaxy, the measurements have uncertainties, and the mass of the stars is an estimate. One galaxy alone cannot settle which explanation is correct.',
         action:'Open the theory guide below to read more. One last thing to try: switch to MOND and watch in fullscreen. There is no halo, yet the outer stars stay put, because MOND’s boosted gravity does the halo’s job. Sandbox tour, below, lets you build a galaxy of your own.'}
      ]},
      sandbox:{label:'Sandbox tour',steps:[
        {id:'sbIntro',mode:'custom',target:'sandboxGroups',gesture:'scan',title:'Build your own galaxy',
         text:'The Custom Sandbox swaps the real galaxy for one you design. You set the mass of its stars and gas and how they are spread out; the simulator does the rest.',
         action:'There are three groups of sliders: stellar disk, bulge and gas.'},
        {id:'sbDisk',mode:'custom',target:'sbDiskMass',gesture:'adjust',title:'Make the disk heavier or lighter',
         text:'The disk is the flat, spinning part of a galaxy, where most stars are. More mass lifts the cyan line, and the halo grows too, since galaxies with more stars sit in bigger halos.',
         action:'Drag the disk mass up and down. Watch the cyan line, and the halo mass slider further down the page.'},
        {id:'sbSize',mode:'custom',target:'sbDiskScale',gesture:'adjust',title:'Spread the stars out',
         text:'The scale length sets how quickly the disk fades outwards. The truncation radius, below it, is where the visible matter ends.',
         action:'Change the scale length, then the truncation radius. Watch the peak of the cyan line move and the 3D galaxy change size.'},
        {id:'sbBulge',mode:'custom',target:'sbBulgeMass',gesture:'increase',title:'Add a bulge',
         text:'A bulge is the round ball of stars at the centre of some galaxies. It puts a lot of mass in a small space, so it speeds up the innermost orbits.',
         action:'Raise the bulge mass above zero. The inner part of the cyan line climbs, and the effective radius slider switches on so you can set the bulge’s size.'},
        {id:'sbGas',mode:'custom',target:'sbGasMass',gesture:'adjust',title:'Add the gas',
         text:'Galaxies also hold cold gas, mostly hydrogen, in a disk three times wider than the stars, so it matters most far out.',
         action:'Change the gas mass and watch the outer part of the cyan line.'},
        {id:'sbSummary',mode:'custom',target:'sandboxSummary',gesture:'scan',title:'Check the totals',
         text:'This list adds your galaxy up. You never set the dark matter directly: the simulator gives your galaxy the halo that galaxies with this mass of stars typically have, found by matching galaxy surveys to simulated halos (Moster et al. 2013).',
         action:'Move any slider above and watch the typical halo mass change with it.'},
        {id:'sbMonitor',mode:'custom',nfw:true,halo:true,target:'stabilityMonitor',gesture:'scan',title:'Meet the halo stability monitor',
         text:`Under the graph, the monitor checks that your halo has a believable shape. It reads the concentration, c = r₂₀₀ ÷ r_s: how tightly the dark matter is packed towards the centre. The shaded band on the gauge is the range this simulator treats as realistic, c from ${STABLE_C_MIN} to ${STABLE_C_MAX}. The marker is your halo.`,
         action:'Find the marker on the gauge. Is your halo inside the band?'},
        {id:'sbStability',mode:'custom',nfw:true,halo:true,target:'scaleSlider',scroll:'stabilityMonitor',gesture:'adjust',title:'Break the halo',
         text:'A tiny scale radius crams the dark matter into the centre (c too high); a huge one smears it thinly outwards (c too low). Real halos do not look like either. The monitor turns red and a warning appears on the 3D view too, but nothing is blocked: you can still run an impossible galaxy.',
         action:'Drag the scale radius to one end, then the other, and watch the marker leave the band.'},
        {id:'sbHalo',mode:'custom',nfw:true,target:'massSlider',gesture:'adjust',title:'Tune the halo yourself',
         text:'Halo mass runs from 0 to 5 times the typical halo. The stars start on orbits that balance the typical halo, so a lighter or heavier one throws them off balance.',
         action:'Try a light halo, then a heavy one. Typical Halo, under the slider, puts both halo sliders back.'},
        {id:'sbFs',mode:'custom',nfw:true,viewPanel:true,target:'toggleFullscreen',gesture:'click',title:'Watch your galaxy',
         text:'As in the main tour, the real test is watching your galaxy run. The tour comes with you into fullscreen.',
         action:'Press Fullscreen.'},
        {id:'sbOrbits',fs:true,viewPanel:true,marks:['bulge'],target:'toggleOrbits',gesture:'click',title:'Two kinds of orbit',
         text:'Disk stars trace neat circles. Bulge stars, if your galaxy has a bulge, trace tangled loops tilted every way: a bulge is held up by its stars’ random motions, not by spinning, which is why it stays round. The ringed star is one of them.',
         action:'Switch on Orbits, and switch off Gas and Halo for a clearer view. No tangle in the middle? Your galaxy has no bulge.',
         done:()=>orbitLines.visible,doneText:'Circles in the disk, a tangle in the bulge.'},
        {id:'sbFsEscape',fs:true,nfw:true,marks:['escaper'],ring:true,target:'fsMassSlider',gesture:'zero',title:'Starve it of dark matter',
         text:'The red circle is the ejection line. Take the halo away and the ringed star is the first your galaxy loses. A heavier disk holds on to more of its stars: try it after the tour.',
         action:'Drag Halo mass to 0, then press Restart.',
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
    let tutorialTour='main',tutorialSteps=TOURS.main.steps,tutorialIndex=0,tutorialTarget=null,tutorialCompleted={};
    let tutorialLastMass=0,tutorialLastScale=0,tutorialLastTime=0,tutorialDrag=null,tutorialShortcutUnlocked=false,tutorialBlockTarget=null,tutorialAimTarget=null;
    const tutorialStep=()=>tutorialSteps[tutorialIndex];
    // Watch steps finish by themselves: each render checks the step's done() test (camera angle, layers, play time).
    let tutorialPlayedMyr=0,tutorialWatchLast=0;
    function tutorialWatchTick(){
      if(tutorialCard.hidden)return;
      const prev=tutorialWatchLast;
      if(simTime>prev)tutorialPlayedMyr+=simTime-prev;
      tutorialWatchLast=simTime;
      const step=tutorialStep();
      if(step.done&&!tutorialCompleted[step.id]&&!tutorialBlocker()&&step.done(prev))completeTutorialAction(step.id,step.doneText);
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
      if(step.mode&&step.mode!==galaxyMode){
        return{target:sourceToggle.querySelector(`input[value=${step.mode}]`).parentElement,gesture:'click',
          say:step.mode==='custom'?'Switch to Custom Sandbox first. ':'Switch back to Preset Galaxies first. '};
      }
      if(step.nfw&&mondMode){
        return{target:document.querySelector('#modelBar .track'),gesture:'click',
          say:step.halo?'Switch back to NFW first; MOND does not use a dark matter halo. ':'Switch back to NFW first. '};
      }
      if(step.mond&&!mondMode)return{target:document.querySelector('#modelBar .track'),gesture:'click',say:'Switch to MOND first. '};
      // Page settings come first: fixing them needs the page, so fullscreen is asked for last.
      if(step.fs&&!fullscreen){
        const viewClosed=viewPanel.classList.contains('collapsed');
        return{target:viewClosed?viewToggle:fsBtn,gesture:'click',say:viewClosed?'Open View in the 3D panel and press Fullscreen first. ':'Press Fullscreen first. '};
      }
      if(step.halo&&Number(massSlider.value)===0)return{target:fullscreen?fsMass:massSlider,gesture:'increase',say:'Increase halo mass above 0 first so there is a halo to see. '};
      if(step.viewPanel&&viewPanel.classList.contains('collapsed'))return{target:viewToggle,gesture:'click',say:'Open View first. '};
      return null;
    }
    function tutorialInstruction(){
      const step=tutorialStep(),blocker=tutorialBlocker();
      if(blocker)return blocker.say+step.action;
      if(step.id==='removeHalo'){
        if(tutorialShortcutUnlocked)return 'The Visible Matter Only button that just appeared is a shortcut back to 0. Compare the orange line with the white points on the right of the graph.';
        if(Number(massSlider.value)===0)return 'Halo mass is already 0. Drag it up, then all the way back left to 0, and watch the orange line.';
      }
      return step.action;
    }
    function refreshTutorialGuidance(){
      if(tutorialCard.hidden)return;
      tutorialTry.textContent='Try it: '+tutorialInstruction();
      updateTutorialSpotlight();
    }
    function scrollTutorialTarget(){
      const step=tutorialStep(),blocker=tutorialBlocker();
      const target=blocker?blocker.target:(step.scroll&&innerWidth>960?document.getElementById(step.scroll):tutorialTarget);
      const wide=!blocker&&(step.gesture==='scan'||step.gesture==='drag');
      const top=target.getBoundingClientRect().top+window.scrollY-(wide?16:75);
      window.scrollTo({top:Math.max(0,top),behavior:'auto'});
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
      tutorialLaunch.textContent=galaxyMode==='custom'?'Sandbox tour':'Tutorial';
      tutorialRecheck();
    });
    function completeTutorialAction(id,message){
      if(tutorialCard.hidden||tutorialStep().id!==id||tutorialCompleted[id])return;
      tutorialCompleted[id]=message;
      tutorialFeedback.hidden=false;
      tutorialFeedback.textContent='✓ '+message;
      tutorialSpotlight.classList.add('is-complete');
    }
    function tutorialIsFullscreen(){return !!nativeFsElement()||orbitPlot.classList.contains('fullscreen-active');}
    // Keeps the highlight glued to its control when the control or the page changes size without a scroll or window
    // resize: a hover box opening inside it, a <details> section expanding above it, the sandbox panel appearing.
    const spotlightResize=window.ResizeObserver?new ResizeObserver(()=>{if(!tutorialCard.hidden)updateTutorialSpotlight();}):null;
    if(spotlightResize)spotlightResize.observe(document.body);
    function updateTutorialSpotlight(){
      if(tutorialCard.hidden||tutorialIsFullscreen()){tutorialSpotlight.hidden=true;return;}
      const step=tutorialStep(),blocker=tutorialBlocker();
      let target=blocker?blocker.target:tutorialTarget,gesture=blocker?blocker.gesture:step.gesture;
      if(gesture==='switch'){
        target=document.querySelector('#modelBar .track');gesture='click';
      }else if(gesture==='guide'){
        target=document.querySelector('#theoryGuide summary');gesture='click';
      }
      // A new control to aim at brings the pointer back; it fades again once the user's mouse reaches it.
      if(target!==tutorialAimTarget){
        if(spotlightResize){if(tutorialAimTarget)spotlightResize.unobserve(tutorialAimTarget);spotlightResize.observe(target);}
        tutorialAimTarget=target;tutorialSpotlight.classList.remove('is-reached');
      }
      const rect=target.getBoundingClientRect();
      if((!rect.width&&!rect.height)||rect.bottom<0||rect.top>innerHeight||rect.right<0||rect.left>innerWidth){tutorialSpotlight.hidden=true;return;}
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
      if(tutorialCard.hidden)return;
      if(fullscreen){
        completeTutorialAction('fsEnter','Fullscreen: the tour carries on in here. Press Next.');
        completeTutorialAction('sbFs','Fullscreen: the tour carries on in here. Press Next.');
      }
      tutorialRecheck();
    }
    function showTutorialStep(){
      const step=tutorialStep(),lastStep=tutorialIndex===tutorialSteps.length-1;
      const backHadFocus=document.activeElement===tutorialBack;
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialLastTime=Number(timeSlider.value);
      tutorialDrag=null;
      tutorialAimTarget=null;
      tutorialTarget=document.getElementById(step.target);
      const blocker=tutorialBlocker();
      tutorialBlockTarget=blocker?blocker.target:null;
      tutorialProgress.textContent=`${TOURS[tutorialTour].label} · ${tutorialIndex+1} / ${tutorialSteps.length}`;
      tutorialProgressFill.style.width=`${100*(tutorialIndex+1)/tutorialSteps.length}%`;
      tutorialTitle.textContent=step.title;
      tutorialText.textContent=step.text;
      tutorialTry.textContent='Try it: '+tutorialInstruction();
      tutorialFeedback.hidden=!tutorialCompleted[step.id];
      tutorialFeedback.textContent=tutorialCompleted[step.id]?'✓ '+tutorialCompleted[step.id]:'';
      tutorialSpotlight.classList.toggle('is-complete',!!tutorialCompleted[step.id]);
      document.getElementById('tutorialContent').scrollTop=0;
      tutorialBack.disabled=tutorialIndex===0;
      tutorialNext.textContent=lastStep?'Finish':'Next';
      tutorialMore.hidden=!(lastStep&&tutorialTour==='main');
      tutorialCard.classList.toggle('on-left',innerWidth>960&&!step.cardRight);
      setStarMarks(step.marks);
      ejectionRing.scale.setScalar(ejectionRadius);ejectionRing.visible=!!step.ring||ejectionLineShown;
      // The ejection steps pull the camera back far enough to frame the ring, so the escaping star stays in view.
      if(step.ring&&tutorialIsFullscreen())controls.flyTo({phi:0.35,radius:ejectionRadius*4,target:[0,0,0],rate:2.4});
      tutorialPlayedMyr=0;tutorialWatchLast=simTime;
      syncTutorialFullscreen();
      if(tutorialBack.disabled&&backHadFocus)tutorialNext.focus({preventScroll:true});
      if(!tutorialIsFullscreen()){
        scrollTutorialTarget();
      }
      updateTutorialSpotlight();
    }
    function startTutorial(tour){
      if(tutorialWelcome.open)tutorialWelcome.close();
      tutorialTour=TOURS[tour]?tour:'main';
      tutorialSteps=TOURS[tutorialTour].steps;
      tutorialIndex=0;
      tutorialCompleted={};
      // In the main tour the Visible Matter Only shortcut is revealed by the step that sets the halo to zero.
      visibleMatterOnly.hidden=tutorialTour==='main'&&!tutorialShortcutUnlocked;
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
      setStarMarks([]);ejectionRing.visible=ejectionLineShown;
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
      // Spot the puzzle is done once the user has looked at the right-hand half of the graph, where the gap is.
      if(tutorialStep().id==='puzzle'&&x>r.left+r.width*.5)completeTutorialAction('puzzle','Spotted: out on the right, the white points stay well above the cyan line. Press Next to find out why.');
      completeTutorialAction('sbMonitor','Found it. Next, push the marker out of the band.');
    }
    window.addEventListener('pointermove',tutorialPointerReach,{capture:true,passive:true});
    window.addEventListener('pointerdown',tutorialPointerReach,{capture:true,passive:true});
    // Capture also catches scrolling inside panels, not just the page.
    window.addEventListener('scroll',updateTutorialSpotlight,{capture:true,passive:true});
    window.addEventListener('resize',updateTutorialSpotlight);
    const sandboxInputs=SANDBOX_SLIDERS.map(slider=>slider.input);
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
        if(value>tutorialLastMass)completeTutorialAction('addHalo','Halo mass increased: more mass means stronger gravity and faster orbits.');
        if(value!==tutorialLastMass){
          completeTutorialAction('fit','Halo changed: check how χ² responded.');
          completeTutorialAction('sbHalo','Halo changed: press Play to see what it does to the orbits.');
        }
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
        completeTutorialAction('haloShape','Halo shape changed: the same mass, arranged differently.');
        completeTutorialAction('fit','Halo changed: check how χ² responded.');
        checkStabilityStep();
      }
      tutorialLastScale=value;
    });
    document.querySelectorAll('#haloMetrics button').forEach(button=>{
      button.addEventListener('click',()=>completeTutorialAction('haloNumbers','Card opened: these numbers follow the halo sliders.'));
    });
    visibleMatterOnly.addEventListener('click',()=>{
      if(!mondMode)completeTutorialAction('removeHalo','Halo removed: the orange line now shows visible matter only.');
      tutorialLastMass=Number(massSlider.value);
      tutorialRecheck();
    });
    document.getElementById('applyObservedFit').addEventListener('click',()=>{
      if(!mondMode)completeTutorialAction('fit','Starting fit restored.');
      tutorialLastMass=Number(massSlider.value);
      tutorialLastScale=Number(scaleSlider.value);
      tutorialRecheck();
    });
    timeSlider.addEventListener('input',()=>{tutorialLastTime=Number(timeSlider.value);});
    // Opening View, or a layer or camera button, can clear what a step was waiting for.
    viewToggle.addEventListener('click',tutorialRecheck);
    ['toggleGas','toggleHalo','toggleOrbits','viewFaceOn','viewEdgeOn','viewReset'].forEach(id=>document.getElementById(id).addEventListener('click',()=>{requestRender();refreshTutorialGuidance();}));
    fullHaloBtn.addEventListener('click',()=>completeTutorialAction('fsFullHalo','Zoomed out: the whole visible galaxy is that bright speck.'));
    // Restart (fullscreen corner panel): back to 0 Myr and play, so a new halo can be watched from the start.
    document.getElementById('fsRestart').addEventListener('click',()=>{onSliderInput(0);setPlaying(true);});
    document.querySelectorAll('#theoryGuide details').forEach(details=>{
      details.addEventListener('toggle',()=>{if(details.open)completeTutorialAction('finish','Theory guide opened.');});
    });
    galaxySelect.addEventListener('change',()=>{
      completeTutorialAction('galaxies','New galaxy loaded: its own data, visible matter and starting halo.');
    });
    ['pointerenter','focus','click'].forEach(type=>{
      document.getElementById('mlInfoBtn').addEventListener(type,()=>completeTutorialAction('ml','Explainer opened.'));
    });
    const SANDBOX_STEP_FEEDBACK={
      sbDiskMass:['sbDisk','Disk mass changed: the typical halo followed it.'],
      sbDiskScale:['sbSize','Disk reshaped: the same mass, spread differently.'],
      sbTruncation:['sbSize','Galaxy edge moved: the ejection radius moves with it.'],
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
    mondFitToggle.addEventListener('change',()=>{
      if(mondFitToggle.checked)completeTutorialAction('mondFit','Best fit on: a₀ and the mass-to-light ratio have been re-tuned for this galaxy.');
    });
    document.getElementById('tutorialStart').addEventListener('click',()=>startTutorial('main'));
    document.getElementById('tutorialWelcomeSkip').addEventListener('click',()=>tutorialWelcome.close());
    tutorialLaunch.addEventListener('click',()=>startTutorial(galaxyMode==='custom'?'sandbox':'main'));
    document.getElementById('tutorialSkip').addEventListener('click',endTutorial);
    tutorialBack.addEventListener('click',()=>{if(tutorialIndex>0){tutorialIndex--;showTutorialStep();}});
    tutorialNext.addEventListener('click',()=>{
      if(tutorialIndex===tutorialSteps.length-1){endTutorial();return;}
      tutorialIndex++;showTutorialStep();
    });
    // The main tour's last step offers the sandbox tour; pressing it opens the Custom Sandbox.
    tutorialMore.addEventListener('click',()=>{
      setGalaxyMode('custom');
      startTutorial('sandbox');
    });
    document.addEventListener('keydown',e=>{
      if(tutorialCard.hidden||e.repeat||e.altKey||e.ctrlKey||e.metaKey||e.shiftKey)return;
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