"""Physics checks for the Interactive Galaxy Rotation Simulator.

Run before every commit that touches the simulation:

    python tests/run_physics_checks.py

If a number changed on purpose (e.g. a deliberate physics change), accept the
new values with:

    python tests/run_physics_checks.py --update-baseline

Only uses the standard library. Needs Node.js for the checks that run the
page's JavaScript (https://nodejs.org, LTS version).
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))  # the generator imports analysis.sparc_data
GENERATOR = ROOT / "generate_interactive_html.py"
INDEX = ROOT / "index.html"
GALAXIES_JS = ROOT / "galaxies.js"
JS_CHECKS = Path(__file__).resolve().parent / "physics_checks.js"
BASELINE = Path(__file__).resolve().parent / "physics_baseline.json"

results: list[tuple[bool, str, str]] = []


def check(passed: bool, name: str, detail: str) -> None:
    results.append((bool(passed), name, detail))


def load_generator():
    spec = importlib.util.spec_from_file_location("generator_under_test", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalise(text: str) -> str:
    return text.replace("\r\n", "\n")


# ---------------------------------------------------------------- Python side
def python_checks() -> None:
    gen = load_generator()

    # 1. The deployed files really come from the generator.
    with tempfile.TemporaryDirectory() as tmp:
        gen.OUTPUT_PATH = Path(tmp) / "index.html"
        gen.GALAXIES_JS_PATH = Path(tmp) / "galaxies.js"
        with contextlib.redirect_stdout(io.StringIO()):
            gen.main()
        fresh_index = normalise(gen.OUTPUT_PATH.read_text(encoding="utf-8"))
        fresh_galaxies = normalise(gen.GALAXIES_JS_PATH.read_text(encoding="utf-8"))
    index_ok = fresh_index == normalise(INDEX.read_text(encoding="utf-8"))
    galaxies_ok = fresh_galaxies == normalise(GALAXIES_JS.read_text(encoding="utf-8"))
    check(
        index_ok and galaxies_ok,
        "Generated files match the generator",
        "index.html and galaxies.js are up to date"
        if index_ok and galaxies_ok
        else "Out of date: run `python generate_interactive_html.py`, or index.html/galaxies.js "
        "was edited by hand (edit the generator instead)",
    )

    # 2. Gas sign: SPARC's negative V_gas must be kept, never squared away.
    ngc3198 = gen.load_sparc_galaxy("NGC3198")
    negative = sum(1 for v in ngc3198["rotation"]["v_gas_kms"] if v < 0)
    squared = re.findall(r"gas\s*\*\*\s*2", GENERATOR.read_text(encoding="utf-8"))
    check(
        negative > 0 and not squared,
        "Gas sign kept (generator)",
        f"{negative} negative V_gas points kept for NGC 3198; no `gas**2` in the generator"
        if not squared
        else f"Found `gas**2` in the generator ({len(squared)}x): use gas*abs(gas)",
    )

    # 3. Star count can't run away.
    check(
        gen.BASE_STAR_COUNT == 800 and gen.MAX_STAR_COUNT <= 2400,
        "Star count limits",
        f"template {gen.BASE_STAR_COUNT} stars, cap {gen.MAX_STAR_COUNT} (must stay <= 2400)",
    )


# ------------------------------------------------------------ JavaScript side
def javascript_checks() -> dict | None:
    node = shutil.which("node")
    if node is None:
        check(False, "Node.js available",
              "Node.js not found, so the page's physics could not be tested. Install the LTS version "
              "from https://nodejs.org and reopen the terminal.")
        return None
    run = subprocess.run([node, str(JS_CHECKS), str(INDEX)], capture_output=True, text=True,
                         encoding="utf-8", cwd=ROOT)
    if run.returncode != 0:
        check(False, "Page script runs", "The page's JavaScript crashed under test:\n" + run.stderr[-1500:])
        return None
    m = json.loads(run.stdout)
    stars = m["star_count"]

    check(m["shift_kpc_halo_mass_0_vs_default"] > 1 and m["shift_kpc_halo_mass_x3"] > 0.5,
          "Halo mass slider moves the stars",
          f"removing the halo moves stars {m['shift_kpc_halo_mass_0_vs_default']:.1f} kpc on average; "
          f"tripling it moves them {m['shift_kpc_halo_mass_x3']:.1f} kpc (0 would mean fake physics)")
    check(m["shift_kpc_scale_radius"] > 0.5, "Scale radius slider moves the stars",
          f"average shift {m['shift_kpc_scale_radius']:.1f} kpc")
    check(m["escaped_visible_only"] >= 0.10 * stars and m["escaped_with_halo"] <= 0.02 * stars,
          "Without dark matter, stars escape",
          f"visible matter only: {m['escaped_visible_only']}/{stars} escape; "
          f"with the default halo: {m['escaped_with_halo']}/{stars}")
    ratio = m["median_abs_z_kpc"] / m["median_R_kpc"]
    check(m["median_abs_z_kpc"] > 0.01 and ratio < 0.1, "Disk is 3D but thin",
          f"median height {m['median_abs_z_kpc']:.2f} kpc vs median radius {m['median_R_kpc']:.1f} kpc")
    check(m["negative_gas_points"] > 0 and m["baryonic_worst_relative_mismatch"] < 1e-5,
          "Gas sign kept (page)",
          f"visible-matter curve matches V_gas*|V_gas| + 0.5 V_disk^2 + 0.7 V_bulge^2 exactly "
          f"(worst mismatch {100 * m['baryonic_worst_relative_mismatch']:.4f}%)")
    check(m["reduced_chi2_visible_only"] > 100 and m["reduced_chi2_nfw_default"] < 5,
          "NGC 3198 chi-squared makes sense",
          f"reduced chi2: visible only {m['reduced_chi2_visible_only']:.0f}, "
          f"NFW default {m['reduced_chi2_nfw_default']:.2f}, MOND {m['reduced_chi2_mond']:.2f}")
    check(m["mond_a0"] == 1.2e-10 and m["mond_halo_mass_used"] == 0
          and m["mond_shift_when_halo_sliders_move"] == 0,
          "MOND mode uses no dark matter",
          f"a0 = {m['mond_a0']:.1e} m/s^2, halo mass {m['mond_halo_mass_used']}, "
          f"halo sliders move stars {m['mond_shift_when_halo_sliders_move']:.2f} kpc")
    check(m["mond_deep_limit_relative_error"] < 1e-3 and 0 <= m["mond_strong_gravity_relative_boost"] < 0.01,
          "MOND formula limits",
          "matches Newton where gravity is strong and v^4 = G M a0 where it is weak")
    check(stars == 800 and m["star_count_for_huge_galaxy"] <= 2400, "Star count on the page",
          f"NGC 3198 uses {stars} stars; a huge sandbox galaxy is capped at {m['star_count_for_huge_galaxy']}")
    for gid, p in m["presets"].items():
        ok = "error" not in p and p["reduced_chi2_nfw_default"] < p["reduced_chi2_visible_only"]
        check(ok, f"{gid} preset loads and fits",
              p.get("error") or f"reduced chi2: visible only {p['reduced_chi2_visible_only']:.0f}, "
              f"NFW default {p['reduced_chi2_nfw_default']:.2f}")
    return m


# ----------------------------------------------------------------- Baseline
def baseline_values(m: dict) -> dict:
    keep = {k: m[k] for k in ("star_count", "default_halo", "escaped_with_halo", "escaped_visible_only",
                              "reduced_chi2_visible_only", "reduced_chi2_nfw_default", "reduced_chi2_mond")}
    keep["presets"] = m["presets"]
    return keep


def compare(old, new, path="") -> list[str]:
    changes = []
    if isinstance(old, dict):
        for k in old.keys() | new.keys():
            if k not in old or k not in new:
                changes.append(f"{path}{k}: added/removed")
            else:
                changes += compare(old[k], new[k], f"{path}{k}.")
    elif isinstance(old, (int, float)) and isinstance(new, (int, float)):
        if abs(new - old) > max(1e-9, 0.01 * abs(old)):
            changes.append(f"{path.rstrip('.')}: {old:.4g} -> {new:.4g}")
    elif old != new:
        changes.append(f"{path.rstrip('.')}: {old} -> {new}")
    return changes


def baseline_check(m: dict, update: bool) -> None:
    current = baseline_values(m)
    if update or not BASELINE.exists():
        BASELINE.write_text(json.dumps(current, indent=2) + "\n", encoding="utf-8")
        check(True, "Baseline", f"{'Updated' if update else 'Created'} {BASELINE.name}")
        return
    changes = compare(json.loads(BASELINE.read_text(encoding="utf-8")), current)
    check(not changes, "Results unchanged since the baseline",
          "same as the saved baseline" if not changes else
          "Changed: " + "; ".join(sorted(changes)) +
          ". If this was intended, run with --update-baseline and say why in the commit.")


def main() -> int:
    update = "--update-baseline" in sys.argv
    try:
        python_checks()
    except Exception as error:  # report rather than crash
        check(False, "Python checks ran", repr(error))
    measurements = javascript_checks()
    if measurements is not None:
        baseline_check(measurements, update)

    print("\nPhysics checks\n==============")
    for passed, name, detail in results:
        print(f"[{'PASS' if passed else 'FAIL'}] {name}\n       {detail}")
    failed = sum(not p for p, _, _ in results)
    print(f"\n{len(results) - failed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
