# EPQ Rotation Curves

Python models and a static GitHub Pages dashboard for exploring the rotation
curves and stellar orbits of SPARC galaxies (NGC 3198, NGC 2403 and NGC 7331).

## Set up on Windows

Install Python 3.12 (64-bit) and Git for Windows, then open this folder in
VS Code. Accept the recommended Python and Pylance extensions if prompted.
In a PowerShell terminal at the project root, create and populate the virtual
environment:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

VS Code is configured to select `.venv` automatically. If PowerShell blocks
environment activation, run the commands without activating it by replacing
`python` with `.\.venv\Scripts\python.exe`. Restart VS Code after installing
Python or Git so it can discover the new executables.

## Run the Python tools

Run commands from the project root:

```powershell
python NGC3198_Baryonic_Baseline.py
python .\Baseline_Baryonic_Model_3d
python visualize_3d_model.py
python generate_interactive_html.py
python generate_qr.py
```

`Baseline_Baryonic_Model_3d` uses the local `table2.dat` and downloads SPARC
data only if that file is missing or invalid. `NGC3198_Baryonic_Baseline.py`
fetches SPARC data from the network. The 3D Blender model script,
`NGC3198_Blender_Model.py`, must be run inside Blender 4.0 or later; Blender's
`bpy` module is not part of the project virtual environment.

## Website

The live dashboard is
<https://will-codes-astro.github.io/Interactive_Galaxy_Rotation_Simulator/>.
GitHub Pages deploys from the `main` branch using
[`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml).
The workflow publishes only `index.html`, `epq_qr_code.png`, and
`epq_qr_code.svg` from the repository root.

`generate_interactive_html.py` builds the self-contained dashboard in
`index.html` from `results/star_trajectories_baryonic.csv` and `table2.dat`.
When changing dashboard behavior, update the generator and regenerate
`index.html` so the published page includes the change. The page loads Plotly
from a CDN; Python packages are not needed to preview the existing page:

```powershell
python -m http.server 8000
```

Then open <http://localhost:8000/>. Commit the regenerated `index.html` and
`galaxies.js` when the dashboard is meant to change. Keep the checked-in data
and generated QR images in sync with any change to the page's data or
destination URL.

The Galaxy dropdown above the gravity model switch chooses which SPARC galaxy
the dashboard models. Choosing a galaxy swaps in its rotation curve, rebuilds
the baryonic curve, resets the halo sliders to that galaxy's observed NFW fit,
restarts the orbit clock and re-integrates the star orbits. The selected
gravity model (NFW or MOND) is kept. The 3D stars follow the selected galaxy
too: a rough disk scale length, outer radius, bulge and visible mass are
estimated from its SPARC mass-model speeds (`anatomyFromRotation`), and these
set the star population described under "Test stars" below.

The preset galaxies live in `galaxies.js`, a lookup keyed by SPARC id that
holds each galaxy's radii, observed speeds, error bars and gas, disk and bulge
speed contributions. The generator creates it from `table2.dat`; do not edit it
by hand. The generator also inlines the same code into `index.html`, so the
published page stays a single file and the Pages workflow needs no change. To
offer another SPARC galaxy, add a `GalaxyPreset` to `GALAXY_PRESETS` in
`generate_interactive_html.py` (its id must match the `Name` column of
`table2.dat`) and run the generator again. It reads that galaxy's rows, fits
its starting NFW halo and rewrites both files.

The Preset Galaxies / Custom Sandbox switch beside the Galaxy dropdown opens
a builder with six anatomy sliders: stellar disk mass, disk scale length,
truncation radius, bulge mass, bulge effective radius and gas mass. The bulge
radius is disabled while the bulge mass is zero. The sandbox never blocks a
configuration. Its halo stability monitor computes the concentration
c = r200 / rs of the current NFW halo from the halo mass and scale radius
sliders and shows a "⚠️ Unphysical / Unstable Galaxy Profile" badge (also over
the 3D view) whenever c falls outside 3–25 (`STABLE_CONCENTRATION` in the
generator).

In the sandbox the anatomy sliders drive the rotation curve, the χ² badge, the
halo and the 3D stars live (`PhysicsEngine` and `rebuildCustomGalaxy` in the
generator's page script). The dropdown stays visible as "Compare with": it
picks the SPARC galaxy whose measured points and radii the custom curve and χ²
are compared with, without touching the anatomy.

- Stellar disk and gas are exponential disks cut at the truncation radius, each
  holding the slider's mass inside it; the gas disk is three times more
  extended than the stars. Their speeds come from Freeman's exact Bessel
  solution for the infinite disk minus the annulus beyond the truncation
  radius (softened elliptic-integral rings); beyond it the rings alone are
  used. The bulge is a Hernquist sphere, Φ = −GM/(r + r_b), with r_b as its scale
  radius. The components add in quadrature to the baryonic curve and are
  drawn as thin dotted lines.
- The halo starts at the stellar-to-halo mass relation baseline: M200 is
  Moster et al. (2013) inverted for the stellar mass (disk + bulge) and r_s
  comes from the Dutton & Macciò (2014) mass–concentration relation.
  The existing halo mass slider becomes the fine-tuning control, spanning
  0 to 5 times the baseline (it keeps its multiple when the anatomy moves the
  baseline), and the scale radius follows the relation until it is moved.
  Apply SHMR Baseline restores both. The MOND mass-to-light fit does not apply
  to a sandbox galaxy and is hidden.

Test stars:

- The default preset (NGC 3198) and the default sandbox galaxy are the
  800-star template. Other galaxies scale it by
  (visible mass / template mass)^0.5, limited to 200–2400 stars
  (`BASE_STAR_COUNT`, `MIN_STAR_COUNT`, `MAX_STAR_COUNT` and
  `STAR_COUNT_EXPONENT` in the generator). The "Test stars" header stat, the
  3D counter and the CSV export all read the active total.
- A star is ejected once it is farther than 1.5 × R_max from the centre
  (`EJECTION_RADIUS_FACTOR`), where R_max is the sandbox truncation radius or
  the outermost SPARC radius of a preset. The 3D counter shows the live count
  and the ejection radius.
- With a bulge, a share of the stars proportional to its stellar mass sits in a
  spheroidal Hernquist cluster. Each bulge star gets an isotropic velocity
  dispersion from the Jeans equation in the reference gravity, plus a little
  rotation, so the bulge is held up by pressure and stays put; disk stars start
  on circular orbits.
- Stars start on the orbits that balance a reference galaxy: the observed
  SPARC speeds of a preset, or the baryons plus baseline halo of a sandbox
  galaxy. Moving the halo sliders, or switching to MOND, changes the gravity
  under them, so they can spiral in or be ejected (try Visible Matter Only).
- The orbit integrator uses a cubic-spline gravitational softening of
  0.5–1.0 kpc (`SOFTENING_RANGE_KPC`; 0.5 kpc for a pure disk, following the
  bulge's scale radius when there is one), exactly Newtonian beyond 2.8 ε, so
  dense centres never meet a force singularity.

The 3D view also draws the cold gas as a faint blue cloud (the View panel's
Hide gas button switches it off): an exponential disk
three times as extended as the stars, cut at R_max, whose brightness follows
the gas mass. The Reset view, Face-on and Edge-on distances scale with the
galaxy (NGC 3198 keeps its original framing), and switching galaxy or mode
re-frames the camera. The glass telemetry panel over the 3D view shows the
ejected count and ejection radius, the disk / bulge star split, the disk scale
length and edge, and the halo r200 ("no halo" in MOND or at zero halo mass).

The dashboard offers an optional hands-on tutorial when opened. Visitors can
skip it at any point and restart it using the Tutorial button beside the
gravity model switch. The tour does not reset or change simulator settings.
It covers the View controls and the scale of the full halo. On a keyboard,
the left and right arrow keys navigate the tour when a control is not being
edited. During the tutorial, the Visible Matter Only shortcut appears after
the halo mass is first moved to zero; skipping the tutorial makes it available.
In NFW mode, Visible Matter Only sets the halo mass to zero without changing
the scale radius; Apply Observed Fit restores the starting halo settings.
The View controls start collapsed. The rotation-curve guide combines the
colour key and fit explanation. Select it or a theory box to open it. Halo
values reveal their explanations on hover or when selected with a mouse, touch,
or keyboard.

## Project map

- `Baseline_Baryonic_Model_3d`: orbit simulation that writes trajectory data
  to `results/`.
- `NGC3198_Baryonic_Baseline.py`: SPARC baryonic rotation-curve analysis.
- `NGC3198_Blender_Model.py`: Blender scene and add-on for the interactive
  3D model.
- `generate_interactive_html.py`: static dashboard generator.
- `galaxies.js`: generated registry of the preset SPARC galaxies (also inlined
  into `index.html`).
- `visualize_3d_model.py`, `plot_galaxy1.py`: trajectory plots and animation.
- `results/`, `plots/`: input/output data and generated figures.
- `table2.dat`: local SPARC rotation-curve data.
