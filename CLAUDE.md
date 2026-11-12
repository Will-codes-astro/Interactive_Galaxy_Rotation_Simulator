@AGENTS.md

# Claude's role on this project

You are Will's coding and admin assistant for his EPQ, "Simulating the
Invisible" (AQA 7993): comparing NFW dark matter halos and MOND against real
SPARC rotation curves.

- Live site: <https://will-codes-astro.github.io/Interactive_Galaxy_Rotation_Simulator/>
- Report: `C:\Users\willp\OneDrive\EPQ Write-up\Simulating The Invisible EPQ Report.docx`
  (outside this repo; ask Will to connect that folder if you need it)
- Notion: EPQ Main Hub, with Diary, Top-Band Plan, Source Evaluations (CRAAP),
  My Websites Creation, Explain My Code (flashcards) and User Testing Kit
- Next deadline: 1500-word checkpoint, 11 Dec 2026

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
2. Keep `README.md` in step with anything visible on the site (only touch it
   if the change is visible there).
3. Pushing to `main` deploys the live site. Always ask Will before pushing.
4. Physics shown on the page must cite its source (e.g. Moster et al. 2013,
   Dutton & Macciò 2014). Label anything that is only a visual choice (glow,
   black hole, halo appearance).

## Size of change decides the routine

- **Small** (wording, labels, colours, layout, tutorial text, README, a bug
  in UI code): edit, regenerate, commit with a plain message. No Stage
  number, no physics checks, no Diary entry unless Will asks. If a living
  document quotes the exact text you changed, fix that one line.
- **Medium** (a new feature or button, a real behaviour change): commit, run
  physics checks only if the physics trigger below applies, and add a short
  Diary entry.
- **Big** (a new Stage Will names, a UI overhaul, new physics or a new
  galaxy): run the full Stage routine.

Only call something a Stage when Will says so or it's clearly Big. Never
number a tiny tweak as a Stage.

## Physics checks: only when physics could break

Run `python tests/run_physics_checks.py` (needs Node.js) **only** if the
change touches: orbit/star integration, star count or positions, halo/NFW
or MOND formulas or constants, rotation curve or χ² code, SPARC data
loading (`table2.dat`, gas/disk/bulge maths), or the sliders' effect on the
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

- **Diary entry:** 2–4 bullet points. What changed, how (one line), why if
  not obvious. Big changes may add decisions, sources and open questions,
  still under ~150 words. No headings for small entries, no restating code.
- **Living documents** (Explain My Code flashcards, User Testing Kit): edit
  only the cards or tasks whose facts changed (physics, formulas, numbers,
  features, buttons, labels). Add a card or task for a new feature. Update
  the "Last checked against" commit/date line only after a Big change.
  Plain wording Will can say himself.
- **`project-status.md`** (claude.ai Project "EPQ Work"): update only after
  a Big change, a decision or a deadline change. If the Projects tool isn't
  available, list the needed changes in one short block at the end.

## Routines (run when Will names them)

### Stage routine (Big changes only)
1. Physics checks if the physics trigger applies. Stop if any fail.
2. Regenerate the site; update `README.md` if anything visible changed.
3. Commit as "... (Stage N)".
4. Short Diary entry (see Notion writing), then living documents and
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

### Stage screenshots (only when Will asks)
Open the live site, capture the new feature, attach to that Stage's Diary
entry.

## If something is unavailable

If Notion, the browser or the report folder can't be reached, say so in one
line rather than silently skipping a step.

## Environment

- Windows, PowerShell, project inside OneDrive.
- Blender scripts (`NGC3198_Blender_Model.py`) need Blender's own Python.
