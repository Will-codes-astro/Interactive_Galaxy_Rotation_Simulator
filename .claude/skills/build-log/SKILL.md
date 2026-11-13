---
name: build-log
description: Write to Will's Notion after a change lands - the "My Websites Creation" build log, the living documents (Explain My Code, User Testing Kit) and project-status.md. Use whenever a new feature, visible change, meaningful bug fix, physics change or mistake is committed, or Will asks for screenshots.
---

# Build log and living documents

Keep all Notion writing short. Never create or write to Diary pages.

## Build log: "My Websites Creation" (id 416c4fda4ef4813fa741fda5b6f61f5b, under EPQ Main Hub)

Write without being asked when something worth noting lands: a new feature,
a visible change, a bug fix that matters, a physics change, a mistake. Skip
pure tweaks (colour nudges, wording, repeated fine tuning of one thing).

- One `## <D Mon> — <what was added>` heading per day, dated by the commit
  date (`git log --date=short`), never the chat date. If that day's heading
  exists, add to it and widen its title if needed.
- Bullets: `<b>Short name:</b> what it does, one sentence`. Add why only if
  it isn't obvious, and a source if the physics needs one. Group related
  commits into one bullet.
- Callouts only for a mistake or an open decision. No commit hashes, Stage
  numbers, file lists or slider-by-slider detail.
- Aim for under ~150 words per day.

## Living documents (Explain My Code flashcards, User Testing Kit)

Edit only the cards or tasks whose facts changed (physics, formulas,
numbers, features, buttons, labels). Add a card or task for a new feature.
If a card quotes text you changed, fix that one line. Update the "Last
checked against" commit/date line only after a Big change. Plain wording
Will can say himself.

## project-status.md (claude.ai Project "EPQ Work")

Update only after a Big change, a decision or a deadline change. If the
Projects tool isn't available, list the needed changes in one short block at
the end of your reply.

## Screenshots (only when Will asks)

Open the live site, capture the new feature, attach it to that day's
build-log entry.
