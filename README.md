# EPQ Rotation Curves

Python models and a static GitHub Pages dashboard for exploring the rotation
curve and stellar orbits of NGC 3198.

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

Then open <http://localhost:8000/>. Commit the regenerated `index.html` when
the dashboard is meant to change. Keep the checked-in data and generated QR
images in sync with any change to the page's data or destination URL.

## Project map

- `Baseline_Baryonic_Model_3d`: orbit simulation that writes trajectory data
  to `results/`.
- `NGC3198_Baryonic_Baseline.py`: SPARC baryonic rotation-curve analysis.
- `NGC3198_Blender_Model.py`: Blender scene and add-on for the interactive
  3D model.
- `generate_interactive_html.py`: static dashboard generator.
- `visualize_3d_model.py`, `plot_galaxy1.py`: trajectory plots and animation.
- `results/`, `plots/`: input/output data and generated figures.
- `table2.dat`: local SPARC rotation-curve data.
