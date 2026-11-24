---
name: site-text
description: Keep Will's "Simulating the Invisible Text" editor artifact and the live site's wording in step. Use when Will says "sync my site text" (or similar), and after any commit that changes the site's wording (the hook says so).
---

# Site text editor sync

Will writes the site's wording himself in the editor artifact
<https://claude.ai/artifact/DtmC7rzNkjnfF6zV7mmVHE>. The wording's only source
is `HTML_TEMPLATE` in `generate_interactive_html.py`; every regeneration
writes it to `content/site_text.json`, which the editor shows. His edits
wait in the artifact's `edits` collection (one doc per text id:
`base` = the text he edited from, `text` = his version) until synced.
Never reword his edits: they are his own writing (academic integrity).

## Sync (Will asks, or before republishing the editor)

1. `ArtifactData list` collection `edits` (url above, `limit` 1000,
   `out_dir` = scratchpad). If empty and the editor is up to date, say so.
2. Build one JSON array `[{id, base, text}]` from the docs without a
   `problem` field (or whose `problem` is older than their `savedAt`),
   write it to the scratchpad, run `python site_text.py apply <file>`.
3. `python generate_interactive_html.py` then `python site_text.py check`
   (Node syntax check of the page scripts). If the check fails, revert the
   generator with `git checkout -- generate_interactive_html.py index.html
   content/site_text.json` (only after confirming no other uncommitted work
   in them) and tell Will which edit broke it.
4. Physics checks only if an edit touched a slider, orbit or physics
   claim's numbers (rarely: wording only).
5. Commit "Site text: Will's wording for <pages>" with the generator,
   `index.html` and `content/site_text.json`; push to `main` (Will agreed
   text syncs may go live in chunks, 24 Nov).
6. Delete the applied edit docs (`batch` delete, pinned to their versions).
   For each problem from `apply`, `update` its doc with
   `{problem: "<reason in plain words>"}` so the editor shows it.
7. Republish the editor: Artifact publish `content/site_text_editor.html`
   with `url` above and `files: {"site_text.json": "content/site_text.json"}`
   (read the artifact first if this session hasn't). Then
   `python site_text.py mark-published` and commit
   `content/site_text_published.txt`.
8. Tell Will in one or two lines: how many edits went live, any that
   couldn't and why. If his wording changes a scientific claim, flag it
   (rule C) but don't change it.

## After Claude changes the site's wording

When `content/site_text.json` changes in a commit, the hook reminds you.
Run Sync first (so pending edits aren't stranded), then step 7. An edit
whose text Claude changed meanwhile shows in the editor under "Edits whose
text changed on the site", with Will's version kept.
