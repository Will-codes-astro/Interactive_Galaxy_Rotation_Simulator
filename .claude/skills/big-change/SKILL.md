---
name: big-change
description: Routine to finish a Big change on the EPQ site (UI overhaul, new physics, new galaxy) - checks, regenerate, docs, commit, Notion, AI declaration, then ask before pushing.
---

# Big change routine

1. Run `python tests/run_physics_checks.py` if the physics trigger in
   CLAUDE.md applies. Stop and report if any check fails.
2. Regenerate the site (`python generate_interactive_html.py`). Update
   `docs/WEBSITE.md` if anything visible changed, and `docs/PHYSICS.md` /
   `docs/DATA_DICTIONARY.md` if physics or data changed.
3. Commit with a message naming what was added.
4. Run the `build-log` skill: build-log entry, living documents (including
   their "Last checked against" line) and `project-status.md`.
5. Check Will's Project Q AI declaration still covers the work.
6. Ask Will before pushing. Browser checks are fine; screenshots only if
   Will asks.
