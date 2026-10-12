# Module 2 Summary: GPU-Instanced Stellar Rendering

## What changed
- Stars in the 3D view are now one `THREE.InstancedMesh` (200 spheres, a single draw call) instead of `THREE.Points` sprites.
- A custom GLSL shader gives each star a 3D look: limb darkening, illumination from the galactic core, and a rim glow.
- Bulge stars (r < 3.5 kpc) are amber and brighter. Disk stars (r >= 3.5 kpc) shade from white to cyan to azure with radius.
- An illuminated 3D nucleus and soft glow replace the old core point.
- Playback now uses a `requestAnimationFrame` delta-time clock and linearly interpolates between trajectory frames. This replaces the 150 ms `setTimeout` stepping.
- Star positions are written straight into the instance matrix array, with no per-frame allocations. Slider and label updates happen only when the integer frame changes.
- The changes are in the generator, [generate_interactive_html.py](generate_interactive_html.py), and the published [index.html](index.html) was regenerated.

## Physics
- Orbits are unchanged: they come from the RK4 integration in the page, using the SPARC baryonic velocity plus the NFW halo (`V_barÂ² + V_haloÂ²`).
- Rendering only reads these positions, so the halo sliders still reshape the orbits.
- Between frames, star positions are linearly interpolated (1 Myr per frame).

## Verification
- The page was reloaded in the browser with no errors.
- Playback, the time slider, the halo mass slider (0 and the observed fit), and the face-on, edge-on and reset views all work.
- Frame rate measured about 60 FPS, with a matrix-update cost of about 0.005 ms per frame.

## Why it is faster than V1.0 (Plotly)
- One instanced draw call replaces Plotly's per-trace rebuilds.
- Updates are in-place typed-array writes, so there is no garbage collection pressure.
- The Three.js canvas only redraws when the camera, time or scene changes.

## Notes
- Bulge versus disk is split at a fixed 3.5 kpc radius. It is a visual choice, not taken from SPARC data.
- The star colours are illustrative, not computed from physical temperatures.

## V2.0 Particle Engine (implemented)

Implemented in [generate_interactive_html.py](generate_interactive_html.py) and regenerated into `index.html`. Tunables live in the `GALAXY` object. There are 1720 massless test particles (seeded PRNG), integrated in 3D with RK4 in the spherical baryonic + NFW potential from the observed SPARC speed. The mix is a visual choice, not SPARC data.

| Component | Count | Radial profile | Vertical profile | Rendering |
|---|---|---|---|---|
| Bulge | 240 | Hernquist, `a = 1 kpc` | Isotropic (spherical) | `THREE.Points`, soft radial-gradient sprite, `AdditiveBlending`, white centre to pink edge |
| Disk | 1360 | Exponential, `R_d = 4 kpc` | `sech^2`, `z_0 = 0.28 kpc` | `InstancedMesh` spheres, white to blue-white |
| Stellar halo | 120 (about 12% warm knots) | Hernquist-like, `a = 12 kpc`, 5-48 kpc | Isotropic | `THREE.Points`, faint sprites, additive |

- **Star size:** the instanced star radius is `0.32 * 0.75` (`STAR_BASE_RADIUS`).
- **Black hole:** a camera-facing shader quad at (0,0,0), radius 0.7. It draws an opaque shadow disc, a thin Doppler-beamed photon ring and a soft orange glow. A keep-out radius (`bulgeKeepOut`) keeps bulge stars off the shadow. It does not change the orbit model.
- **Purged:** grid and axes helpers, torus/ring meshes and the hard black sphere. Orbit lines are an optional "Orbit lines" toggle in the View panel, and the Reset view button gives a 3/4 view.
- **Sprites:** point size scales with the canvas height and FOV (`uScale`), and `depthWrite` is off so overlaps add into a glow.
