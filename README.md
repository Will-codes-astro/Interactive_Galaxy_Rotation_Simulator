# Interactive Galaxy Rotation Simulator (EPQ)

Will's EPQ, *Simulating the Invisible*, compares two explanations for why
galaxies spin too fast at large radii: a **dark matter halo (NFW)** and
**modified gravity (MOND)**. Both are tested against measured rotation curves
from the **SPARC** survey (Lelli, McGaugh & Schombert 2016), mainly NGC 3198.

- **Live site:** <https://will-codes-astro.github.io/Interactive_Galaxy_Rotation_Simulator/>
- **Main features:** rotation curve with observed SPARC speeds and gas, disk and
  bulge components; NFW halo sliders; MOND mode and best fit; χ² fit score;
  800 test stars integrated live in 3D; NGC 3198, NGC 2403 and NGC 7331 presets;
  a custom sandbox galaxy; guided tutorials; CSV export.
- **Offline analysis:** a baryonic-only baseline fit and orbit simulation for
  NGC 3198 (Python), plus a Blender 3D model.

## Documentation index

| I want to… | Read |
|---|---|
| work on the repo (rules for Claude and other assistants) | [CLAUDE.md](CLAUDE.md), [AGENTS.md](AGENTS.md) |
| understand or change a physics equation, constant or unit | [docs/PHYSICS.md](docs/PHYSICS.md) |
| know what a data column means, or where a file came from | [docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) |
| change or explain something visible on the site | [docs/WEBSITE.md](docs/WEBSITE.md) |
| see old development notes | [docs/history/](docs/history/) |

## Repository layout

```
data/sparc/table2.dat          raw SPARC data (read-only, checksum-checked)
analysis/                      offline Python
  sparc_data.py                  the only SPARC loader: validates, never writes
  baryonic_baseline.py           baryonic-only curve + χ²  -> results/, plots/
  simulate_baryonic_orbits.py    baryonic-only 2D orbits    -> results/star_trajectories_baryonic.csv
  plot_orbits_2d.py, plot_orbits_3d.py   plots of that CSV  -> plots/
blender/                       Blender scene script and saved .blend (needs Blender's Python)
generate_interactive_html.py   builds the website: index.html + galaxies.js
index.html, galaxies.js        generated website (do not edit by hand)
epq_qr_code.png/.svg, generate_qr.py   QR code to the live site
results/, plots/               generated outputs of analysis/
tests/                         physics checks (Python + Node.js)
docs/                          documentation; artefact_appendix.tex = LaTeX QR appendix
```

The main physics lives in two places. The website version is in
`generate_interactive_html.py`: the NFW fit is in the Python part, and the
NFW, MOND, χ² and orbit code are in the page JavaScript. The offline baseline is
in `analysis/`. [docs/PHYSICS.md](docs/PHYSICS.md) maps every equation to its
function.

## Set up on Windows

Install Python 3.12 (64-bit), Git for Windows and Node.js (Node is used by the
physics checks). In PowerShell at the project root:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

VS Code selects `.venv` automatically. If PowerShell blocks activation, use
`.\.venv\Scripts\python.exe` in place of `python`.

## Run

From the project root:

```powershell
python generate_interactive_html.py      # rebuild index.html and galaxies.js
python tests/run_physics_checks.py       # physics checks (add --update-baseline only for intended changes)
python -m http.server 8000               # preview at http://localhost:8000/
python analysis/baryonic_baseline.py
python analysis/simulate_baryonic_orbits.py
python analysis/plot_orbits_2d.py
python analysis/plot_orbits_3d.py
python generate_qr.py
```

`blender/ngc3198_blender_model.py` must be run inside Blender 4.0 or later.

## Data

All code reads the one SPARC table, `data/sparc/table2.dat`, through
`analysis/sparc_data.py`. That loader checks the file's checksum, header and
value ranges before any calculation runs, and it never downloads or overwrites
the file. Provenance and column definitions are in
[docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md). To add a SPARC galaxy to
the site, add a `GalaxyPreset` to `GALAXY_PRESETS` in the generator and run it again.

## Deployment

Pushing to `main` runs `.github/workflows/deploy-pages.yml`, which publishes
only `index.html`, `epq_qr_code.png` and `epq_qr_code.svg` from the repository
root. Those files must stay in the root. The page loads Plotly and Three.js
from CDNs.
