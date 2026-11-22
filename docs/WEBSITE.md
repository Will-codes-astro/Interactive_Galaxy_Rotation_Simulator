# Website features

What the live site does, control by control. Code: `generate_interactive_html.py`
(the page's HTML, CSS and JavaScript are in `HTML_TEMPLATE`). The equations are
in [PHYSICS.md](PHYSICS.md). Keep this file in step with anything visible on the site.

The site is one file (`index.html`) split into five pages: The puzzle,
Simulator (Model setup and Simulation), Scorecard, Background, and Method &
sources. The top menu shows one page at a time (`showPage`); each page ends with
a "Next" button, the browser's Back button works, and a link such as
`#scorecard` or `#curvePanel` opens the right page. A hidden page keeps its
width (height 0, not `display:none`), so the 3D view and the rotation curve stay
correctly sized. If Plotly or Three.js fails to load (for example a network that blocks CDNs) or
WebGL is off, a red note at the top says which, and the page shows The puzzle,
Background and Method as one long scroll; the Simulator and Scorecard, which need
the libraries, are hidden. On phones (640 px or narrower) the menu uses short labels
("Puzzle", "Method") in one row (below 420 px the logo is hidden to make room); on
narrower screens it wraps onto a second row, so every entry stays visible. The stat tiles (test stars, timeline, mass-to-light, tilt) and the Best realistic
fits switch sit just under the Simulation heading, directly above the sliders
(Will's choice), all in one row on wide screens, so the values in use stay next to
the controls that change them. On wide screens the galaxy profile card splits into
two columns (who the galaxy is | its numbers) to keep the figures nearer the top.
The gravity model choice ("Gravity model", `#modelBar`) is the first control in
the slider bar (Will's choice), so on screens at least 961 × 700 px, where that bar
sticks to the top, the theory can be changed while watching the figures and χ².
It is two colour-coded cards (Will's choice): "Dark matter halo: extra unseen
mass" (purple, like the halo line) and "MOND: gravity changes, no extra mass"
(gold). The cards set the hidden `#mondToggle` checkbox that the rest of the page
listens to. Its hint suggests trying MOND in NFW mode and gives a₀ in MOND mode.
With MOND on, for every galaxy (preset or sandbox), the rotation curve's speed panel shades in gold every
radius where the visible matter's gravity is weaker than a₀, faint at the edge and deeper where MOND's boost is bigger ("MOND boosts
gravity here", in small gold text in the zone's empty bottom-right corner so it
never covers a curve), with a faint dotted line where it crosses a₀ (see PHYSICS.md). The
sandbox tour's "Where MOND steps in" card asks for MOND, then a heavier disk
until an unshaded centre opens up (MOND already chosen when the card opens counts
as done). On that card and "Spread the stars out", a finished first part's ring
turns green and fades before the second part and its ring appear. The halo shape
monitor card only asks the visitor to look at the monitor; clicking it (or Next)
moves on. Read-only cards (the sandbox totals, the halo
shape monitor, the Scorecard) get a box but no pointer arrow. A slider's box also
takes in its description lines and follows any change to their length. Switches
that are already in a card's end state when it opens (orbit lines on, gas or halo
off, already playing, fullscreen, the ejection line, MOND chosen) count as done.
Next from the sandbox tour's last card in the full tour leaves fullscreen before
Wrapping up opens.

The tutorial offer appears the first time the Simulator page opens in each page
load, not over the opening story. "I've been here before" closes it for the rest
of that page load; a refresh or a new tab offers it again (Will's choice:
nothing is stored between visits). Each tour step opens the page its control is on (the last two
open the Scorecard, then Background); if the visitor moves to another page mid-tour, the card asks
them to go back and points at the top-menu entry.

The puzzle page opens the site, before any control. Four numbered steps
explain what a rotation curve is, what astronomers expected (with the Solar
System as the comparison), what they measured, and the two possible answers.
Beside them a figure (`drawPuzzleFigure`) plots the default galaxy's real SPARC
speeds against the speed its visible matter predicts, with the gap between them
shaded. The speeds and the "about 5 times the pull" ratio, (V_obs / V_bar)² at
the outermost point, are read from the data, not typed in. The shading is a
visual aid only. "Take the full guided tour" runs the simulator tutorial and
then, straight on, the sandbox tutorial and a short "Wrapping up" part ("Part 1 of 3" to
"Part 3 of 3"). The simulator tutorial's two closing cards (Scorecard, Background) are held back
for that last part, so the full tour ends on them. "Open
the simulator" instead brings up the welcome offer, whose "Start simulator
tutorial" runs the simulator tutorial alone.

The Galaxy menu (Model setup) chooses which SPARC galaxy
the dashboard models. Choosing a galaxy swaps in its rotation curve, rebuilds
the baryonic curve, resets the halo sliders to that galaxy's observed NFW fit,
restarts the orbit clock and re-integrates the star orbits. The selected
gravity model (NFW or MOND) is kept. The 3D stars follow the selected galaxy
too: a rough disk scale length, outer radius, bulge and visible mass are
estimated from its SPARC mass-model speeds (`anatomyFromRotation`), and these
set the star population described under "Test stars" below.

The Galaxy menu is a custom listbox (`galaxyMenuBtn`, `galaxyMenuList`; arrow
keys, Home/End, Enter/Space, Escape), so each option shows its Hubble type badge
and the preset's one-line `summary`. Under it, the galaxy profile card
(`renderGalaxyProfile`) shows SPARC Table 1's values in plain words: the Hubble
type with its family and arms (`HUBBLE_TYPES`), whether there is a central bulge (from SPARC's mass model: any non-zero V_bulge, so NGC 3198 shows none although Sc galaxies often have a small one), distance (also in light-years),
flat rotation speed (also in km/h), how far out it was measured and how many
speeds, and a stars-vs-gas bar (stars ≈ 0.5 × L[3.6], gas = 1.33 × M_HI, the
helium factor SPARC uses). The stellar figure is a rough estimate (it ignores the
bulge's higher M/L). The tutorial's "Try another galaxy" step boxes the menu in blue and
puts a yellow watch box on the card.

The preset galaxies live in `galaxies.js`, a lookup keyed by SPARC id that
holds each galaxy's radii, observed speeds, error bars and gas, disk and bulge
speed contributions. There are four: three spirals (NGC 3198, NGC 2403,
NGC 7331) and one low-surface-brightness galaxy (UGC 128), each shown in the
galaxy picker and the scorecard. The generator creates it from `data/sparc/table2.dat` and `table1.mrt`; do not edit it
by hand. The generator also inlines the same code into `index.html`, so the
published page stays a single file and the Pages workflow needs no change. To
offer another SPARC galaxy, add a `GalaxyPreset` to `GALAXY_PRESETS` in
`generate_interactive_html.py` (its id must match the `Name` column of
`data/sparc/table2.dat`) and run the generator again. It reads that galaxy's rows, fits
its starting NFW halo and rewrites both files.

Model setup opens with two large choice cards: Preset Galaxies (real data,
the SPARC galaxies) and Custom Sandbox (build your own). The chosen card glows
blue; on phones they stack. Choosing Custom Sandbox opens
a builder with six anatomy sliders: stellar disk mass, disk scale length,
truncation radius, bulge mass, bulge Hernquist scale radius and gas mass. The bulge
radius is disabled while the bulge mass is zero. A "What you have built" card
above the sliders names the galaxy type (dwarf, gas-rich dwarf irregular, low
surface brightness, Sa/Sb/Sc spirals, bulgeless disk, lenticular, giant or
monster) from the star mass, bulge-to-total ratio, gas fraction and disk
surface density (`classifyGalaxy` in the generator), with a real example and
trait tags. The thresholds are rough rules of thumb. The sandbox never blocks a
configuration. Its halo shape monitor, under the rotation curve, computes
the concentration c = r200 / rs of the current NFW halo from the halo mass and
scale radius sliders and shows a "⚠️ Unusual halo shape"
badge (otherwise "✓ Typical halo shape") (also over the 3D view) whenever c falls outside 3–25
(`STABLE_CONCENTRATION` in the generator).

In the sandbox the anatomy sliders drive the rotation curve and the halo live,
every frame; the 3D stars and their orbits (about 0.1 s of work) are rebuilt
once the slider pauses (`PhysicsEngine`, `rebuildCustomGalaxy` and
`rebuildCustomStars` in the generator's page script). A sandbox galaxy is imaginary, so it has no measured speeds: the
Galaxy menu and profile card, the observed SPARC points, the residual bars and the χ² badge
are hidden, and the curves are drawn on an even grid of 60 radii out to
1.25 × the truncation radius, plus four radii inside the first step so every
curve starts at the centre; the radius axis starts at 0. Its CSV export leaves out the observed and
residual columns.

- Stellar disk and gas are exponential disks cut at the truncation radius, each
  holding the slider's mass inside it; the gas disk is three times more
  extended than the stars. Their speeds come from Freeman's exact Bessel
  solution for the infinite disk minus the annulus beyond the truncation
  radius (softened elliptic-integral rings); beyond it the rings alone are
  used. The bulge is a Hernquist sphere, Φ = −GM/(r + r_b), with r_b as its scale
  radius. The components add in quadrature to the baryonic curve and are
  drawn as thin dotted lines.
- The halo starts at the typical halo for the galaxy's stellar mass (the
  stellar-to-halo mass relation): M200 is
  Moster et al. (2013) inverted for the stellar mass (disk + bulge) and r_s
  comes from the Dutton & Macciò (2014) mass–concentration relation.
  The existing halo mass slider becomes the fine-tuning control, spanning
  0 to 5 times the baseline. Any change to the anatomy sliders resets both
  halo sliders to the new typical halo.
  The Typical Halo button (under the halo mass slider) restores both. The MOND mass-to-light fit does not apply
  to a sandbox galaxy and is hidden.
- The header's Mass-to-light (M/L) and Tilt tiles are hidden in the sandbox,
  because sandbox masses are set directly rather than converted from light,
  and a sandbox galaxy has no measured tilt.

Hovering over the "?" on the M/L tile fades in a short explainer (tap or
press it on touch screens and keyboards): what a mass-to-light ratio
is, why young stellar populations have low M/L and old ones high M/L, and that
the presets' 0.5 (disk) and 0.7 (bulge) at 3.6 μm are standard averages
(Lelli, McGaugh & Schombert 2016; McGaugh & Schombert 2014), not values
measured for each galaxy, and that the Best realistic fits switch lets each
model move them between about 0.6 and 1.6 times those values. It fades out when the pointer leaves, or on ×,
Escape or a tap elsewhere.

The Tilt tile beside it shows the inclination in use (SPARC's, or the best
fit's with SPARC's in brackets) and has the same kind of "?" explainer: what
inclination is, that telescopes only see line-of-sight motion so true speed =
measured ÷ sin i, that near face-on galaxies are the most sensitive, and that
Best realistic fits may move it within twice SPARC's error (Lelli, McGaugh &
Schombert 2016; Li et al. 2018). Method has a matching "Tilt (inclination)"
card with the s = sin i_SPARC / sin i factor and a worked example.

Test stars:

- The Test stars tile's "?" opens the same kind of hover explainer as the M/L
  and Tilt tiles: they are test particles (they feel the model's gravity but do
  not pull on each other), start on balanced orbits, can be ejected, and their
  number is a display choice.
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
  dense centres never meet a force singularity. Method's Simulation card names
  both simplifications (spherical pull; softening, Springel 2005) and why
  velocity-Verlet keeps the energy error from building up.

Each simulated star is drawn as a point of light whose colour and brightness
come from a mass drawn from the Kroupa (2001) IMF (L ∝ M^3.5, T ∝ M^0.475,
blackbody colour, brightness on a log scale). Bulge stars are old (≤ 1 M☉ plus
a few red giants) and the disk gets bluer outwards (de Jong 1996). Underneath,
a diffuse glow shows the unresolved starlight: an exponential disk with the
galaxy's own scale length (Freeman 1970) and a de Vaucouleurs bulge, weighted
by their masses. Colour saturation, diffraction spikes and bloom are visual
choices only. To stay smooth on large or high-DPI screens, the 3D view lowers
its render resolution in steps while frames run slow and raises it again when
there is headroom (`trackFrameTime` in the generator); 4x antialiasing is used
only on low-DPI screens at full resolution.

While a halo or sandbox slider moves, the rotation curve is redrawn every frame
through Plotly's light transition path (`updateHalo(true)`); the speed axis only
grows mid-drag, and a full redraw fits the axes once the
slider settles. Residual bars are one neutral grey (their direction shows the
sign), and the fit verdicts are green, yellow and red, so neither shares a
colour with a curve (`RESIDUAL_COLOR`).

The 3D view (the "3D galaxy simulation" panel) also draws the cold gas as a faint blue cloud (the View panel's
Gas switch turns it off): an exponential disk
three times as extended as the stars, cut at R_max, whose brightness follows
the gas mass. The Reset, Face-on and Edge-on distances scale with the
galaxy (NGC 3198 keeps its original framing), and switching galaxy or mode
re-frames the camera. In fullscreen a corner panel holds Play, Restart (back
to 0 Myr and play), the clock and the two halo sliders. A small, faint readout in the top-left corner of the 3D
view shows only the ejected count and the halo r200 ("no halo" in MOND or at
zero halo mass).

The halo, time and tracking sliders sit in one bar directly above the 3D
view and rotation curve. On a large enough screen the bar stays pinned under
the top menu while you scroll the plots, so every slider change is visible at
once, and the tutorial scrolls so the bar and the plots are both on screen.
In that bar a long readout (e.g. the orbit time) drops under its label rather than
running into the next column, and between 961 and 1280 px Play shows its icon only.
On phones (640 px or narrower) the order is 3D view, rotation curve, sliders, and
the rotation curve stays pinned under the menu while the sliders scroll beneath it,
so a slider's effect is visible while dragging (layout choice only). It is unpinned
while a tutorial runs, since the tutorial card needs the room. The curve's legend
uses shorter names there ("Observed", "Visible matter") so two fit on each row.

Next to the orbit timeline, "Track a star" takes a star number and rings that
star in the 3D view with its laps round the centre since 0 Myr, speed and
distance (also in the fullscreen corner panel, beside the round Play and
Restart buttons). It only reads the computed orbits.

The dashboard offers an optional hands-on tutorial when opened. Visitors can
skip it at any point and restart it using the Simulator tutorial button beside the
Model setup heading. The tours do not reset or change simulator settings.
While a card is open, every simulator control it does not name is locked
(inert and dimmed), and so are the top menu and the 3D camera's drag and zoom
(except from "Now watch it happen" until Face-on is done on "Look down on the
disk", and on the sandbox tour's "Watch your galaxy" and "Two kinds of orbit",
where the camera moves as normal and Reset, Face-on, Edge-on and Full halo are
free too); when something is in the way, only
the control that fixes it is free. Playback speed, the info (?) buttons and the
halo-card explainers stay free. A box ring covers an opened info box only when
the box opens right beside its control, and if an opened info box would sit
under the tutorial card, the card moves to the other side. Ending or
skipping the tour unlocks everything. Instructions with several actions show
one numbered action at a time.
The main tutorial's 19 steps start from the puzzle in the rotation curve
(measured speeds stay above the visible-matter prediction; the pointer must
cover all three stretches of the graph's right half), then remove and
restore the halo (adding it back is done once the model reaches 90% of the
measured speeds over the outer third of the points), change its shape, open all three halo cards (r200, scale
radius, concentration), then fit the halo by hand: break the fit with Halo
mass, then tune both sliders until the badge reads Close fit (Apply Observed
Fit stays locked). For a galaxy whose published halo fit is not Close, the goal
is χ²ᵣ within 1.25 × that fit's value. Then it compares MOND and its best-fit
switch. A "watch it happen" chapter then runs inside fullscreen, where the
tour card follows: hide the gas and halo, look face-on with orbit lines,
predict which of two ringed stars will finish more laps (Play stays locked
until a prediction is chosen), then watch them race (live lap counts and
speeds show differential rotation on a flat rotation curve; the lap counts start at zero when the race card
opens, so they cannot give the answer away, and start again if the clock is taken back), turn edge-on to see a ringed star bob
through the thin disk, switch on the ejection line, remove the halo and follow
the first star to cross it, then restore the halo and zoom out to the full halo.
It ends with the Galaxy menu, the mass-to-light explainer, the Scorecard
(read the At a glance table, then press Next) and the theory guide (no box or dimming, since it covers a whole section). The rings and labels are tutorial aids drawn from the
already computed trajectories; they never change the simulation. Watch steps
tick themselves off once the camera, layers or clock reach the right state. An animated pointer shows where to go on
each step; it fades out as soon as the mouse or a finger reaches the
highlighted control, and comes back only when the step points somewhere new.
Each card has one or two short sentences of explanation and an
"Instructions:" line asking for one action; a card with several actions
("1. … 2. …") shows one part at a time ("Part 1 of 2"), and each new part fades
in. Controls to use get a blue box; the numbers or graph to watch while doing
it (χ², M/L, tilt, the graph, the Ejected count) get yellow dashed boxes; a
slider's own value needs none, as the slider's blue box already holds it. On
"Give both their best shot" the χ², M/L and tilt boxes stay after Best
realistic fits is switched on, while those values change. A step, or one part of a many-part card, finished by dragging a slider
only counts as done once the slider is let go, so no highlight moves on mid-drag.
When a task is done the highlight turns green with one small ripple (every
finished highlight holds green 1.4 s, then fades over 1.2 s), the card
border tints green briefly, and the done message starts with a short word of
praise; a 10 s countdown (paused while a control is held, restarted by more
input) then fades to the next card. Green only ever marks something done on that
card: Next skips a card without a tick, a card that opens already in its end state
(because an earlier card was skipped) only ticks once its part is undone and redone,
leaving a card clears its rings at once, and going back to a finished card shows
its message without replaying the glow. Back always returns to a skipped card so it
can be done. On screens wider than 700 px every card shows all its text without
scrolling: a card too long for its space widens step by step (to 540 px at
most), and in fullscreen can also grow up to just below the view buttons. On
phones the card keeps its full-width layout and scrolls if it must. A second, 12-step Sandbox tutorial covers
the Custom Sandbox: the disk, bulge and gas sliders, the totals list, the
halo shape monitor (what it reads, then pushing a halo out of range), the
halo sliders, and a fullscreen chapter (disk circles versus tangled bulge
orbits, then the first star lost without dark matter). It follows on from the simulator
tutorial in the full guided tour, or starts from the Sandbox tutorial button on
the simulator tutorial's last step (which opens the Custom Sandbox), or from the
Simulator tutorial button, which reads "Sandbox tutorial" while the sandbox is
open. On the slider steps whose result shows far from the slider, the card
carries a copy of the one slider the current part asks for (`cardSlider` on the
step), linked to the real one, and the page scrolls to the graph, 3D view or
shape monitor it changes. A slider in the pinned slider bar (wide screens) needs no copy.
If a step needs a different mode (NFW or MOND, a halo above
zero, the other galaxy source, fullscreen or leaving it), the card says so and points at the control. Each halo and
time slider has a plain-language description underneath it, and the galaxy
picker shows a one-line summary with the number of measured speeds. On a keyboard,
the left and right arrow keys always navigate the tour while it is open, even
when the last control used (a slider, the Galaxy menu, a button) still has
focus; only a text box being typed in keeps them, and up and down still move a
focused slider. The highlight on a slider step rings the slider together with
its label row (name and value); the sandbox tutorial's first card points at
each of the three slider groups' headings in turn. Outside the tour, in fullscreen only, the arrows
step the time by 50 Myr (hold to scrub). During the tutorial, the Visible Matter Only shortcut appears after
the halo mass is first moved to zero; skipping the tutorial makes it available.
In MOND mode the halo mass and scale radius sliders (and their fullscreen
copies) are hidden, not greyed out, since MOND has no halo; on wide screens the
remaining time and Track a star controls share the bar in two columns.
Best realistic fits switch (header, beside the M/L tile; presets only): off by
default, so every model uses SPARC's standard M/L, distance and tilt and the
cyan visible-matter curve is the same in NFW and MOND mode. On, the NFW and
MOND best fits each tune those within the same realistic limits (MOND also
a₀); see PHYSICS.md. It replaced the MOND-only "MOND best fit" switch.
In NFW mode, Visible Matter Only sets the halo mass to zero without changing
the scale radius and returns to SPARC's standard M/L, distance and tilt.
Dragging the halo mass to zero does the same, and dragging it back up
restores the fit's M/L, distance and tilt until Visible Matter Only is pressed; Apply
Observed Fit restores the halo fit for the switch's state (halo only, or the
NFW best fit's halo, M/L, distance and tilt).
The M/L tile shows the values in use and, under them, where they come from
("standard SPARC values", "NFW best fit" or "MOND best fit"); the Method page's
visible-matter formula shows the same numbers and label. The fit badge's tooltip
names the fit's distance and tilt.
The View controls start collapsed. The rotation-curve guide combines the
colour key and fit explanation. Select it or a theory box to open it. Halo
values reveal their explanations on hover or when selected with a mouse, touch,
or keyboard.

Hover explainers (the halo value cards and the M/L, Tilt and Test stars "?"
boxes) open only after the mouse rests on them for 0.35 s, so sweeping past
does not open them, and fade in over 0.45 s. Opening one halo card closes
the others. In the tutorial, an explainer counts as read only once it has
stayed open for 1.5 s (`HOVER_OPEN_DELAY_MS`, `EXPLAINER_READ_MS`). While a
tutorial card is up, the page is dimmed except for every control the card
names and any explainer it has opened.

Highlight boxes: each box follows its control every frame, including while it fades out green, so a box
never hangs where a control used to be. A box whose control has scrolled under the top menu (or, for the
plots, under the pinned slider bar) is cut off at that edge rather than drawn over it. Each card scrolls
its control to just below the top menu, a slider together with its label and value, and brings a watched
box that sits just above it (the M/L and tilt tiles) or further down (the shape monitor) into view
too. The box round the Galaxy menu grows to cover its open list, and the yellow box on the profile card
hides while the list is over it. The box round a "?" button and its explainer has square-ish corners.
The pointer fades out once a step is done.

On screens 1200 px or wider the two figure headings keep one fixed height and a long fit label wraps
inside its own space, so the graph does not move when the label changes. The stat tiles are all the same
height. With no halo (MOND, or zero halo mass) the "Your NFW halo" box drops its "hover over a value" line.

Fit badge: above the rotation curve, a word ("Close fit", "Rough fit" or "Poor
fit") sits in front of χ²ᵣ; χ² itself, the point count and the values tuned are in
the number's tooltip. The grey label beside it ends with the same plain-words gap
as the Scorecard ("about 1.2 error bars off", √(χ²/N)). Beside the badge a grey label names what is being
scored and how many values were tuned (`fitLabel`): "NFW halo fit · 2 tuned",
"NFW best fit · 5 tuned",
"Your halo" (sliders moved, with or without the fit's M/L, distance and tilt),
"Visible matter only", "Plain MOND · 0 tuned" or "MOND best fit · 4 tuned". Its
tooltip says every model uses the same settings (standard, or each model's
best realistic fit). The cut-offs (χ²ᵣ below 2, 2 to 10, above 10)
are the site's own rule of thumb (`fitVerdict`), and the guide under the graph
says so. χ²ᵣ divides by the measured points minus the values tuned to the galaxy
(`tunedCount`: the halo's 2 or a₀'s 1, plus the M/L, distance and tilt values a
best fit tunes; visible matter only and plain MOND 0); see PHYSICS.md. The k
values quoted in the Method page and the Scorecard's "Tuned values" card are
written by the generator (`tuned_count_text`) from the same numbers, so they
cannot drift from the page's maths.

Scorecard: first an "At a glance" table, one row per galaxy, giving the verdict
word for visible matter only and, for NFW and MOND, standard values → best
realistic fit (a † marks a best fit that needed a value at a realistic limit; the
galaxy name opens it in the simulator). Under it, "Full scores, error bars and tuned
values" (closed by default) opens the full table, scoring every preset galaxy in two column groups.
"Standard values": visible matter only, the NFW halo fit (halo only) and MOND
at the default a₀. "Best realistic fits": the NFW and MOND best fits. The halo
and best-fit cells list how many values were tuned and what they
are (halo mass and concentration, or a₀; M/L; distance and tilt beside SPARC's). A value that stopped at a realistic limit carries a small "at limit" tag whose tooltip says which limit. Under each score it gives the root-mean-square gap in plain words
("about 1.2 error bars off", √(χ²/N)); χ² itself is in the score's tooltip.
On phones (640 px or narrower) the table becomes one stacked card per galaxy,
with each score labelled by its column, so nothing is off-screen.
It is calculated in the page when it loads (`scoreGalaxy`), with the
same functions as the fit badge, so the two always agree. "Show in simulator"
loads that galaxy (leaving the sandbox if needed). Six short cards underneath (including
"Only realistic values", which names the sources of the fit limits, and "Why MOND
scores badly here": tiny error bars, hard limits, and Li et al. 2018's acceptable
fits for most SPARC galaxies with soft limits) state
what limits the comparison; they do not say which idea wins.

Background has a third, full-width box, "Beyond rotation curves: the wider
evidence": galaxy clusters, the Bullet Cluster, the cosmic microwave background
and the radial acceleration relation, each with a tag for the idea it points
to and two short lines, "Seen:" (the observation, with its source) and "Why it
matters:". Two "Still open" cards follow: no confirmed dark matter
particle (including the single LZ candidate event of September 2026, a 2.6σ
preprint, arXiv:2609.02823) and MOND beyond galaxies. The MOND box folds a
three-step A-level derivation (g = v²/r, MOND's strong and weak limits,
v⁴ = GMa₀) under "Show the maths"; the full interpolating formula (the RAR
function, McGaugh, Lelli & Schombert 2016) appears only in Method, with a "Why this
formula?" note and a small table of plain MOND's χ²ᵣ for each galaxy with the
standard and RAR formulas (generated by `mond_formula_table`), and a line that the RAR
formula was itself fitted to SPARC data, while the halo is fitted per galaxy. The dark matter box has its own "Show the maths": v = √(GM(r)/r)
with the shell theorem, why visible matter alone gives falling speeds, why a flat
curve needs density ∝ 1/r², and the NFW density profile (1/r inside r_s, 1/r³
far outside, about 1/r² near r_s; halo speed peaks near 2.2 r_s), followed by a short "A known weakness" note on the cusp–core problem (de Blok 2010). It also has a matching schematic (not real data) in the simulator's colours:
the visible matter's share of the speed falls, the halo's grows, and together they stay level. Method's
cards flow in two columns (read down the left, then the right; one column at 900 px or narrower).
Method & sources ends with a Harvard reference list of every
source cited on the page, each linked by DOI (checked against Crossref on
14 Nov 2026).

Rotation-curve readout: hovering over the graph (or, on a touch screen, tapping
it or sliding a finger along it) shows a box with the speeds at the nearest
measured SPARC radius. The box, guide line and dots glide smoothly between radii,
but the numbers are always the values at a real data point, never interpolated.
On touch screens the readout stays until you tap outside the graph and updates
live as the sliders move. Touch screens also get larger slider thumbs, tap-to-open
"?" tips and a two-column legend on narrow screens.

"Before" line: while a halo or sandbox anatomy slider is moved, and for 4 s after
it is let go (`CURVE_GHOST_HOLD_MS`), a faint wide band labelled "before" stays
where the line it changes was: the cyan visible-matter line for the anatomy
sliders, the orange total for the halo sliders. The speed axis keeps room for it.
It fades in over 0.3 s and out over 0.9 s.
Changing galaxy, galaxy source or gravity model clears it at once. It is a visual aid only.
