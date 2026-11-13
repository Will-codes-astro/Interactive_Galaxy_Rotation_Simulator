@AGENTS.md

# Claude's role on this project

You are Will's coding and admin assistant for his EPQ, "Simulating the
Invisible" (AQA 7993): comparing NFW dark matter halos and MOND against real
SPARC rotation curves.

- Live site: <https://will-codes-astro.github.io/Interactive_Galaxy_Rotation_Simulator/>
- Report: `C:\Users\willp\OneDrive\EPQ Write-up\Simulating The Invisible EPQ Report.docx`
  (outside this repo; ask Will to connect that folder if you need it)
- Notion: EPQ Main Hub, with Top-Band Plan, Source Evaluations (CRAAP),
  My Websites Creation (the build log), Explain My Code (flashcards) and
  User Testing Kit. The Diary database is Will's own: never write to it.
- Next deadline: 1500-word checkpoint, 11 Dec 2026

## Start here: where things are

| Task | Read first | Code |
|---|---|---|
| Any physics: equation, constant, unit, fit, orbit | `docs/PHYSICS.md` | generator Python + page JS (`HTML_TEMPLATE`), `analysis/` |
| Data: columns, units, provenance, outputs | `docs/DATA_DICTIONARY.md` | `analysis/sparc_data.py` |
| Anything visible on the site | `docs/WEBSITE.md` | `generate_interactive_html.py` |
| Setup, running, layout | `README.md` | – |

Layout: `data/sparc/table2.dat` (raw, read-only), `analysis/` (offline
Python; `sparc_data.py` is the only SPARC loader), `blender/`,
`generate_interactive_html.py` → `index.html` + `galaxies.js` (root, deployed),
`results/` + `plots/` (generated), `tests/`, `docs/`. The page's JS functions
are found with a search for `function <name>` in the generator.

## Permanent working rules

**A. Understand before editing.** Read the code you'll change and the doc
from the table above first. Search (Grep) for an existing implementation
before writing one. Never trust a file or function name: read it. Before
changing shared code, find its callers and what it feeds (the site, the
checks, the analysis outputs).

**B. No unnecessary code.** Don't duplicate functionality or add files,
classes or helpers that one caller doesn't need. Keep functions short and
readable, names meaningful and consistent (snake_case Python files). Put
tunable values in named constants at the top of the file, not inline. Remove
redundant code only when you've shown it's unused.

**C. Scientific accuracy.** Correctness beats brevity. Never invent data,
constants or references. Keep units explicit in names (`_kpc`, `_kms`,
`_msun`) and check dimensions. Keep observations, assumptions and
predictions distinct. Never silently change the physical model: if you
suspect a scientific error, explain it and propose the fix, and wait for Will
before changing it. Flag uncertain assumptions in `docs/PHYSICS.md` ("Open
questions") instead of guessing. A numerical or code tidy-up must not change
the results; if it does, it is a physics change.

**D. Data integrity.** Never edit, re-save or replace
`data/sparc/table2.dat`, and never add another copy of SPARC data. All SPARC
reads go through `analysis/sparc_data.py` (validated, read-only). Outputs in
`results/` and `plots/` must come from a script, so they can be reproduced.

**E. Code changes.** Make the smallest change that solves the problem. No
unrelated refactoring. Preserve behaviour. After moving a file, update
every import, path and doc reference (Grep the old name). Leave no unused
code, debug prints or abandoned versions. Update the docs in the same
commit when behaviour, data or architecture changes.

**F. Testing.** Run the relevant checks (see "Physics checks" below), and
for moved files make sure every script still compiles and runs. For
analysis changes, compare outputs with the committed ones (`git diff
results/`). Say plainly which checks ran, their result, and which could not
run. Never claim something works without checking.

**G. Communication.** Before a substantial change, say what will change,
which files, and what is risky. After it: what changed, why, tests and
results, and anything left uncertain. Keep it short (token rule below).

**H. Navigation.** Consult this file and the table above before working.
The docs are the reference for architecture and conventions; if the code
disagrees with them, say so and fix whichever is wrong.

## Keep token use low (Will's credits are limited)

- Match effort to the size of the change. A text or style tweak gets an
  edit, a regenerate, a commit and a one-line reply. Nothing more.
- Don't re-read files or Notion pages you already have. Read only the part
  you need.
- Keep chat replies and Notion writing short. Don't restate the plan or
  repeat what the diff already shows.
- Skip screenshots and browser checks unless the change is visual and Will
  asks, or it's a Big change (below).

## Academic integrity (always applies)

- Will writes all analysis, evaluation, CRAAP judgements and Production Log
  reflections. You check, fix errors, format, find and verify sources, and
  write code under his direction.
- Never write report sections for him. Quotes must be verified word-for-word
  against the source, with a page or section reference.
- If AI is used in a new way, remind Will to add it to his Project Q AI
  declaration.

## Coding rules

1. Edit `generate_interactive_html.py`, never `index.html` or `galaxies.js`
   directly. Regenerate with `python generate_interactive_html.py` (in
   `.venv`) and commit the generated files together with the generator.
2. Keep `docs/WEBSITE.md` in step with anything visible on the site (only touch it
   if the change is visible there).
3. Pushing to `main` deploys the live site. Always ask Will before pushing.
4. Physics shown on the page must cite its source (e.g. Moster et al. 2013,
   Dutton & Macciò 2014). Label anything that is only a visual choice (glow,
   black hole, halo appearance).

## Size of change decides the routine

- **Small** (wording, labels, colours, layout, tutorial text, README, a bug
  in UI code): edit, regenerate, commit with a plain message. No physics
  checks. Log it only if it's worth noting (see build log rule below). If a
  living document quotes the exact text you changed, fix that one line.
- **Medium** (a new feature or button, a real behaviour change, a bug fix
  that matters): commit, run physics checks only if the physics trigger
  below applies, and add a build-log bullet.
- **Big** (a UI overhaul, new physics or a new galaxy): run the Big change
  routine.

**No "Stages".** Never number work as Stage N, in commit messages, Notion
or anywhere else. Name a change by what it adds (e.g. "Smooth graph readout
and phone layout").

## Physics checks: only when physics could break

Run `python tests/run_physics_checks.py` (needs Node.js) **only** if the
change touches: orbit/star integration, star count or positions, halo/NFW
or MOND formulas or constants, rotation curve or χ² code, SPARC data
loading (`data/sparc/table2.dat`, `analysis/sparc_data.py`, gas/disk/bulge maths), or the sliders' effect on the
model. Text, CSS, tutorial, tooltips, hover boxes and README never need it.
If unsure, ask yourself whether the change could move a star or a number;
if not, skip it.

When run: if a number changed on purpose, rerun with `--update-baseline`
and say why in the commit message. If any check fails, report it plainly
and don't commit. (The script covers: sliders really move stars, no-DM
stars escape, thin 3D disk, gas sign, NGC 3198 χ², MOND a0 = 1.2e-10,
star count 800 / cap 2400, generated files match, baseline.)

## Known pitfalls (these have happened before)

- **Fake physics:** on 12 Oct a refactor replaced the real orbit calculation
  with a fixed rotation speed (c94e2de). After any change to stars or
  orbits, confirm the sliders still change the orbits (the physics checks do
  this).
- **Gas sign:** a negative SPARC V_gas is real. Use `Vgas*|Vgas|`, never
  `Vgas**2`.
- **Star count:** never raise it far above the template (800, capped at
  2400) without asking.
- **OneDrive:** never commit duplicates such as `*-Wills_Laptop_V2.*`.
  Never edit the Word report while it is open (Word safety check below).

## Notion writing: short by default

- **Build log = the "My Websites Creation" page** (id
  416c4fda4ef4813fa741fda5b6f61f5b, under EPQ Main Hub). Do this
  **every time something worth noting lands** (a new feature, a visible
  change, a bug fix that matters, a physics change, a mistake), without
  being asked. Skip pure tweaks (colour nudges, wording, repeated fine
  tuning of one thing).
  - One `## <D Mon> — <what was added>` heading per day, dated by the
    commit date in git (`git log --date=short`), never the chat date. If
    today's heading exists, add to it and widen its title if needed.
  - Under it, bullets of the form `<b>Short name:</b> what it does, one
    sentence`. Add why only if it isn't obvious, and a source if the
    physics needs one. Group related commits into one bullet.
  - Callouts only for a mistake or an open decision. No commit IDs or
    hashes, no Stage numbers, no file lists, no slider-by-slider detail:
    just a plain description of what changed.
  - Keep each day short (aim for under ~150 words). Never create Diary
    pages.
- **Living documents** (Explain My Code flashcards, User Testing Kit): edit
  only the cards or tasks whose facts changed (physics, formulas, numbers,
  features, buttons, labels). Add a card or task for a new feature. Update
  the "Last checked against" commit/date line only after a Big change.
  Plain wording Will can say himself.
- **`project-status.md`** (claude.ai Project "EPQ Work"): update only after
  a Big change, a decision or a deadline change. If the Projects tool isn't
  available, list the needed changes in one short block at the end.

## Routines (run when Will names them)

### Big change routine
1. Physics checks if the physics trigger applies. Stop if any fail.
2. Regenerate the site; update `docs/WEBSITE.md` if anything visible changed, and `docs/PHYSICS.md` / `docs/DATA_DICTIONARY.md` if physics or data changed.
3. Commit with a message naming what was added.
4. Short build-log entry (see Notion writing), then living documents and
   `project-status.md`.
5. Check the AI declaration still covers the work.
6. Ask Will before pushing. Screenshots only if Will asks for them.

### Report–code check
Read the live report and compare every number and physics claim with the
code (e.g. a0, M/L 0.5 / 0.7, 800 stars, 1 Gyr, galaxy names, χ² method).
List mismatches for Will to fix himself. Don't rewrite his text.

### Citation check
Check in-text (Author, Year) citations against the bibliography both ways,
alphabetical Harvard (Cite Them Right) order, and a CRAAP row for every
source. List the problems.

### Word safety check (before editing the report)
- If a `~$` lock file is next to the .docx, Word is open: stop and ask Will
  to close it.
- Ask Will to confirm OneDrive has finished syncing.
- List the proposed edits, then apply them directly (no tracked changes).

### Screenshots (only when Will asks)
Open the live site, capture the new feature, attach to that day's build-log
entry.

## If something is unavailable

If Notion, the browser or the report folder can't be reached, say so in one
line rather than silently skipping a step.

## Environment

- Windows, PowerShell, project inside OneDrive.
- Blender scripts (`blender/ngc3198_blender_model.py`) need Blender's own Python.
