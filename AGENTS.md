# Project instructions

- Read `README.md` before changing setup or deployment behavior.
- The GitHub Pages site is a static dashboard. Its generator is
  `generate_interactive_html.py`; regenerate and include `index.html` whenever
  generator changes are intended to affect the published site.
- `.github/workflows/deploy-pages.yml` publishes only root-level
  `index.html`, `epq_qr_code.png`, and `epq_qr_code.svg`. Keep this deploy
  contract intact unless the workflow is deliberately updated too.
- The dashboard's only data input is `data/sparc/table2.dat` (raw, never
  edit it). Prefer editing the generator over hand-editing generated HTML.
- Declare installable Python dependencies in `requirements.txt`; scripts
  should not install packages at runtime. Run project scripts with `.venv`
  activated, except Blender scripts, which require Blender's Python runtime.
- Avoid committing `.venv/`, `__pycache__/`, or Blender backup files ending
  in `.blend1`.
