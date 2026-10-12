"""Create a standalone Plotly dashboard from the local NGC 3198 data."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


BASE_DIR = Path(__file__).resolve().parent
TRAJECTORY_PATH = BASE_DIR / "results" / "star_trajectories_baryonic.csv"
SPARC_PATH = BASE_DIR / "table2.dat"
OUTPUT_PATH = BASE_DIR / "index.html"

G_KPC_KMS2_PER_MSUN = 4.30091e-6
H0_KMS_PER_KPC = 0.07
UPSILON_DISK = 0.50
UPSILON_BULGE = 0.70
STAR_COUNT = 200


def load_trajectories(path: Path) -> list[dict[str, Any]]:
    by_star: dict[int, dict[str, Any]] = {}
    with path.open("r", newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            star_id = int(row["star_id"])
            star = by_star.setdefault(
                star_id,
                {
                    "id": star_id,
                    "v_rot_kms": float(row["v_obs_kms"]),
                    "initial_r_kpc": float(row["initial_r_kpc"]),
                    "frames": [],
                },
            )
            star["frames"].append(
                (
                    int(row["frame"]),
                    float(row["x_kpc"]),
                    float(row["y_kpc"]),
                    float(row["z_kpc"]),
                )
            )

    if len(by_star) < STAR_COUNT:
        raise ValueError(f"Expected at least {STAR_COUNT} stars in {path}")

    stars = []
    expected_frames: int | None = None
    for star_id in sorted(by_star):
        star = by_star[star_id]
        star["frames"].sort(key=lambda frame: frame[0])
        if expected_frames is None:
            expected_frames = len(star["frames"])
        if len(star["frames"]) != expected_frames:
            raise ValueError(f"Star {star_id} has an incomplete trajectory in {path}")
        star["positions"] = [list(frame[1:]) for frame in star.pop("frames")]
        stars.append(star)

    return stars[:STAR_COUNT]


def load_rotation_curve(path: Path) -> list[dict[str, float]]:
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as source:
        for line in source:
            columns = line.split()
            if len(columns) < 10 or columns[0].upper() != "NGC3198":
                continue
            try:
                radius, observed, error, gas, disk, bulge = map(
                    float,
                    (columns[2], columns[3], columns[4], columns[5], columns[6], columns[7]),
                )
            except ValueError:
                continue
            baryonic = math.sqrt(
                max(
                    0.0,
                    gas**2 + UPSILON_DISK * disk**2 + UPSILON_BULGE * bulge**2,
                )
            )
            rows.append(
                {
                    "radius_kpc": radius,
                    "observed_kms": abs(observed),
                    "error_kms": abs(error),
                    "baryonic_kms": baryonic,
                }
            )

    rows.sort(key=lambda row: row["radius_kpc"])
    if len(rows) < 3:
        raise ValueError(f"No usable NGC 3198 SPARC points found in {path}")
    return rows


def _nfw_shape(x: float) -> float:
    if x <= 0.0:
        return 0.0
    return math.log1p(x) - x / (1.0 + x)


def _nfw_velocity(radius_kpc: float, halo_mass_msun: float, scale_radius_kpc: float) -> float:
    critical_density = (
        3.0 * H0_KMS_PER_KPC**2 / (8.0 * math.pi * G_KPC_KMS2_PER_MSUN)
    )
    r200 = (
        3.0 * halo_mass_msun / (4.0 * math.pi * 200.0 * critical_density)
    ) ** (1.0 / 3.0)
    concentration = max(r200 / scale_radius_kpc, 1e-8)
    enclosed_mass = halo_mass_msun * _nfw_shape(radius_kpc / scale_radius_kpc)
    enclosed_mass /= max(_nfw_shape(concentration), 1e-12)
    if radius_kpc >= r200:
        enclosed_mass = halo_mass_msun
    return math.sqrt(G_KPC_KMS2_PER_MSUN * enclosed_mass / radius_kpc)


def fit_default_halo(curve: list[dict[str, float]]) -> tuple[float, float]:
    """Grid-fit a usable M200 and scale-radius starting point to SPARC."""
    critical_density = (
        3.0 * H0_KMS_PER_KPC**2 / (8.0 * math.pi * G_KPC_KMS2_PER_MSUN)
    )

    def score(mass_billion_solar: float, scale_radius: float) -> float:
        mass = mass_billion_solar * 1e9
        r200 = (3.0 * mass / (4.0 * math.pi * 200.0 * critical_density)) ** (1.0 / 3.0)
        concentration = max(r200 / scale_radius, 1e-8)
        shape_concentration = max(_nfw_shape(concentration), 1e-12)
        chi_squared = 0.0
        for point in curve:
            radius = point["radius_kpc"]
            enclosed_mass = mass * _nfw_shape(radius / scale_radius) / shape_concentration
            if radius >= r200:
                enclosed_mass = mass
            dark_velocity = math.sqrt(G_KPC_KMS2_PER_MSUN * enclosed_mass / radius)
            total_velocity = math.hypot(point["baryonic_kms"], dark_velocity)
            uncertainty = max(point["error_kms"], 1.0)
            chi_squared += ((point["observed_kms"] - total_velocity) / uncertainty) ** 2
        return chi_squared

    coarse = min(
        (
            (score(mass, scale_radius), mass, scale_radius)
            for mass in range(100, 2001, 20)
            for scale_radius in range(5, 61, 1)
        ),
        key=lambda item: item[0],
    )
    _, coarse_mass, coarse_radius = coarse
    fine = min(
        (
            (score(mass, scale_radius), mass, scale_radius)
            for mass in range(max(100, coarse_mass - 40), min(2000, coarse_mass + 40) + 1, 2)
            for scale_radius in range(max(5, coarse_radius - 3), min(60, coarse_radius + 3) + 1)
        ),
        key=lambda item: item[0],
    )
    return fine[1] * 1e9, float(fine[2])


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Interactive Galaxy Rotation Simulator for NGC 3198 orbit and SPARC rotation-curve exploration">
  <title>Interactive Galaxy Rotation Simulator</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
  <style>
    :root { color-scheme: dark; --bg:#050910; --panel:#0a111b; --line:#1b2a3a; --ink:#dbe7f4; --muted:#91a4b8; --cyan:#7be4f1; --gold:#f2c879; }
    * { box-sizing:border-box; }
    body { margin:0; background:var(--bg); color:var(--ink); font:14px/1.45 "Segoe UI",sans-serif; }
    header { display:flex; align-items:flex-end; justify-content:space-between; gap:20px; padding:22px 26px 16px; border-bottom:1px solid var(--line); background:#070d15; }
    .eyebrow { color:var(--cyan); text-transform:uppercase; letter-spacing:.12em; font-size:11px; }
    h1 { font-size:22px; line-height:1.15; margin:5px 0 0; font-weight:600; }
    .header-note { color:var(--muted); font-size:12px; text-align:right; }
    main { padding:16px; max-width:1800px; margin:auto; }
    .controls { display:grid; grid-template-columns:minmax(230px,1fr) minmax(230px,1fr) minmax(280px,1.4fr); gap:14px; padding:13px 16px; border:1px solid var(--line); background:var(--panel); border-radius:5px; margin-bottom:14px; }
    .control label { display:flex; justify-content:space-between; gap:12px; color:var(--muted); font-size:12px; margin-bottom:7px; }
    .control output { color:var(--ink); font-variant-numeric:tabular-nums; }
    .control-help { color:var(--muted); font-size:11px; margin:6px 0 0; min-height:32px; }
    .fit-button { margin-top:7px; font-size:11px; }
    input[type=range] { width:100%; accent-color:var(--cyan); }
    .panel-grid { display:grid; grid-template-columns:1.12fr 1fr; gap:14px; }
    .panel { min-width:0; background:var(--panel); border:1px solid var(--line); border-radius:5px; overflow:hidden; }
    .panel-heading { display:flex; align-items:baseline; justify-content:space-between; gap:12px; padding:13px 16px 0; }
    h2 { margin:0; font-size:14px; font-weight:600; }
    .panel-heading span { color:var(--muted); font-size:11px; }
    .plot { width:100%; height:min(68vh,700px); min-height:420px; }
    #orbitPlot { position:relative; touch-action:none; overscroll-behavior:contain; user-select:none; -webkit-user-select:none; cursor:grab; }
    #orbitPlot.dragging { cursor:grabbing; }
    #orbitPlot canvas { display:block; width:100%; height:100%; touch-action:none; }
    .view-panel { position:absolute; top:12px; right:12px; z-index:2; display:flex; flex-direction:column; gap:6px; padding:10px; min-width:132px; border-radius:12px; border:1px solid rgba(255,255,255,.14); background:rgba(14,24,38,.45); -webkit-backdrop-filter:blur(14px) saturate(140%); backdrop-filter:blur(14px) saturate(140%); box-shadow:0 8px 28px rgba(0,0,0,.35), inset 0 1px 0 rgba(255,255,255,.08); cursor:default; }
    .view-panel-title { color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; padding:0 2px 2px; }
    .view-panel button { width:100%; padding:6px 10px; font-size:12px; text-align:left; background:rgba(255,255,255,.06); border:1px solid rgba(255,255,255,.12); border-radius:8px; transition:background .15s,border-color .15s; }
    .view-panel button:hover, .view-panel button:focus-visible { background:rgba(123,228,241,.16); border-color:var(--cyan); outline:none; }
    .view-panel button[aria-pressed=true] { background:rgba(123,228,241,.22); border-color:var(--cyan); color:var(--ink); }
    .view-panel .view-row { display:grid; grid-template-columns:1fr 1fr; gap:6px; }
    .view-panel .view-row button { text-align:center; }
    .view-panel .mode-toggle { display:grid; grid-template-columns:1fr 1fr; padding:2px; margin-top:2px; border-radius:999px; background:rgba(0,0,0,.28); border:1px solid rgba(255,255,255,.1); }
    .view-panel .mode-toggle button { padding:4px 8px; font-size:11px; text-align:center; border:0; border-radius:999px; background:transparent; color:var(--muted); }
    .view-panel .mode-toggle button[aria-checked=true] { background:rgba(123,228,241,.22); color:var(--ink); box-shadow:inset 0 0 0 1px rgba(123,228,241,.5); }
    #orbitPlot.pan-mode { cursor:move; }
    #orbitPlot.pan-mode.dragging { cursor:grabbing; }
    .view-panel .view-head { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:0 2px 2px; background:transparent; border:0; border-radius:6px; color:var(--muted); font-size:10px; letter-spacing:.12em; text-transform:uppercase; }
    .view-panel .view-head:hover, .view-panel .view-head:focus-visible { background:transparent; color:var(--ink); }
    .view-panel .chev { display:inline-block; font-size:12px; transition:transform .2s; }
    .view-panel.collapsed { min-width:0; padding:7px 10px; border-radius:999px; }
    .view-panel.collapsed > :not(.view-head) { display:none; }
    .view-panel.collapsed .view-head { padding:0; }
    .view-panel.collapsed .chev { transform:rotate(90deg); }
    @media (hover:none) and (pointer:coarse) { .view-panel .mode-toggle { display:none; } .view-panel button { padding:9px 12px; font-size:13px; } }
    .timeline { display:grid; grid-template-columns:auto 1fr auto; gap:12px; align-items:center; padding:0 16px 14px; color:var(--muted); font-size:12px; }
    button { background:#102130; color:var(--ink); border:1px solid #294255; border-radius:4px; padding:7px 12px; cursor:pointer; }
    button:hover { border-color:var(--cyan); }
    footer { margin:14px 2px 0; color:var(--muted); font-size:11px; display:flex; justify-content:space-between; gap:16px; }
    @media(max-width:960px) { .panel-grid { grid-template-columns:1fr; } .plot { height:62vh; } .controls { grid-template-columns:1fr 1fr; } .control:last-child { grid-column:1/-1; } }
    @media(max-width:560px) { header { align-items:flex-start; flex-direction:column; padding:17px; } .header-note { text-align:left; } main { padding:10px; } .controls { grid-template-columns:1fr; } .control:last-child { grid-column:auto; } .plot { min-height:360px; height:56vh; } footer { flex-direction:column; } }
  </style>
</head>
<body>
  <header>
    <div><div class="eyebrow">SPARC galaxy dynamics · interactive model</div><h1>Interactive Galaxy Rotation Simulator</h1></div>
    <div class="header-note">1720 massless test particles (bulge, disk, stellar halo) · 301 frames · 1 Myr per frame<br>Local data embedded in this page · dark halo uses an NFW profile</div>
  </header>
  <main>
    <section class="controls" aria-label="Dark matter controls">
    <div class="control"><label for="massSlider"><span>Halo mass, M<sub>200</sub></span><output id="massValue"></output></label><input id="massSlider" type="range" min="0.00" max="2.00" step="0.01"><p class="control-help">Updates the rotation curve and star orbits. Set it to 0 for visible matter only.</p><button id="applyObservedFit" class="fit-button" type="button">Apply Observed Fit</button></div>
    <div class="control"><label for="scaleSlider"><span>Scale radius, r<sub>s</sub></span><output id="scaleValue"></output></label><input id="scaleSlider" type="range" min="5" max="60" step="0.5"><p class="control-help">Controls halo concentration: a smaller radius packs more dark matter near the center.</p></div>
      <div class="control"><label for="timeSlider"><span>Orbit time</span><output id="timeValue"></output></label><div style="display:flex;align-items:center;gap:12px"><input id="timeSlider" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"><button id="playButton" type="button">Play</button></div></div>
    </section>
    <section class="panel-grid">
    <article class="panel"><div class="panel-heading"><h2>Stellar orbit trajectories</h2></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional NGC 3198 star trajectories"><div class="view-panel" id="viewPanel" role="group" aria-label="3D view controls"><button type="button" class="view-head" id="viewToggle" aria-expanded="true" title="Collapse or expand the view controls"><span>View</span><span class="chev" aria-hidden="true">▾</span></button><button type="button" id="viewReset" title="Return to the default NGC 3198 view">Reset view</button><button type="button" id="toggleOrbits" aria-pressed="false" title="Show or hide the orbit trajectory lines">Orbit lines</button><div class="view-row"><button type="button" id="viewFaceOn" title="Look straight down on the galactic plane">Face-on</button><button type="button" id="viewEdgeOn" title="Look along the galactic plane">Edge-on</button></div><div class="mode-toggle" role="radiogroup" aria-label="Mouse drag mode"><button type="button" role="radio" aria-checked="true" data-mode="rotate">Rotate</button><button type="button" role="radio" aria-checked="false" data-mode="pan">Pan</button></div></div></div><div class="timeline"><span id="frameLabel">Frame 1</span><input id="timeSliderBottom" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"><span id="timeLabel">0 Myr</span></div></article>
      <article class="panel"><div class="panel-heading"><h2>Rotation curve</h2><span>SPARC data and NFW halo response</span></div><div id="curvePlot" class="plot" role="img" aria-label="NGC 3198 observed, baryonic, and fitted rotation curves"></div></article>
    </section>
    <footer><span>Baryons: V<sub>bar</sub>² = V<sub>gas</sub>² + 0.50V<sub>disk</sub>² + 0.70V<sub>bulge</sub>². Halo: NFW, M<sub>200</sub> and r<sub>s</sub> update locally.</span><span>SPARC values from the project’s local table2.dat.</span></footer>
  </main>
  <script>
    const DATA = __DATA_JSON__;
    const DEFAULTS = __DEFAULT_JSON__;
    const G = 4.30091e-6;
    const H0 = 0.07;
    const BG = '#050910';
    const GRID = '#1c2936';
    const textColor = '#bac9d8';
    const stars = DATA.stars;
    const frameCount = stars[0].positions.length;
    const maxRadius = Math.max(...stars.map(s => s.initial_r_kpc));
    const speeds = stars.map(s => s.v_rot_kms);
    const minSpeed = Math.min(...speeds), speedSpan = Math.max(...speeds) - minSpeed || 1;
        const KM_S_TO_KPC_MYR = 0.0010227121650537077;

        /* ---------- 3D galaxy morphology: Hernquist bulge + exponential / sech^2 disk ---------- */
        const GALAXY={bulgeCount:240,diskCount:1360,haloCount:120,hernquistA:1.0,bulgeMax:7,bulgeKeepOut:1.1,
            diskScaleLength:4.0,diskScaleHeight:0.28,diskMax:maxRadius,diskMin:0.8,
            haloA:12.0,haloMin:5,haloMax:48,blackHoleRadius:0.7};
        function mulberry32(seed){return()=>{seed|=0;seed=seed+0x6D2B79F5|0;let t=Math.imul(seed^seed>>>15,1|seed);
            t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
        // Hernquist (1990) inverse CDF: r = a*sqrt(u)/(1-sqrt(u)); isotropic spherical distribution gives a true 3D spherical concentration.
        function sampleBulge(rng){
            let r;
            do{const s=Math.sqrt(Math.min(Math.max(rng(),1e-6),0.999));r=GALAXY.hernquistA*s/(1-s);}
            while(r>GALAXY.bulgeMax||r<GALAXY.bulgeKeepOut);
            const cosT=2*rng()-1,sinT=Math.sqrt(Math.max(0,1-cosT*cosT)),phi=2*Math.PI*rng();
            return[r*sinT*Math.cos(phi),r*sinT*Math.sin(phi),r*cosT];
        }
        function sampleDisk(rng){
            let R;
            do{R=-GALAXY.diskScaleLength*(Math.log(1-rng())+Math.log(1-rng()));}
            while(R>GALAXY.diskMax||R<GALAXY.diskMin);
            const phi=2*Math.PI*rng();
            // sech^2 vertical profile: z = z0 * atanh(2u-1) / 2 samples rho ~ sech^2(z/z0).
            const u=Math.min(Math.max(rng(),1e-4),1-1e-4);
            const z=GALAXY.diskScaleHeight*0.5*Math.log(u/(1-u));
            return[R*Math.cos(phi),R*Math.sin(phi),z];
        }
        function sampleHalo(rng){
            let r;
            do{const s=Math.sqrt(Math.min(Math.max(rng(),1e-6),0.999));r=GALAXY.haloA*s/(1-s);}
            while(r>GALAXY.haloMax||r<GALAXY.haloMin);
            const cosT=2*rng()-1,sinT=Math.sqrt(Math.max(0,1-cosT*cosT)),phi=2*Math.PI*rng();
            return[r*sinT*Math.cos(phi),r*sinT*Math.sin(phi),r*cosT];
        }
        // Particles are ordered bulge then disk; kind selects the orbit initial conditions and colour.
        const particles=(()=>{
            const rng=mulberry32(3198),out=[];
            const randomAxis=()=>{
                const cosT=2*rng()-1,sinT=Math.sqrt(Math.max(0,1-cosT*cosT)),phi=2*Math.PI*rng();
                return[sinT*Math.cos(phi),sinT*Math.sin(phi),cosT];
            };
            for(let i=0;i<GALAXY.bulgeCount;i++){
                const p=sampleBulge(rng);
                out.push({kind:'bulge',position:p,radius:Math.hypot(p[0],p[1],p[2]),axis:randomAxis(),speedScale:1});
            }
            for(let i=0;i<GALAXY.diskCount;i++){
                const p=sampleDisk(rng);
                out.push({kind:'disk',position:p,radius:Math.hypot(p[0],p[1]),axis:randomAxis(),speedScale:1});
            }
            for(let i=0;i<GALAXY.haloCount;i++){
                const p=sampleHalo(rng);
                out.push({kind:'halo',position:p,radius:Math.hypot(p[0],p[1],p[2]),axis:randomAxis(),speedScale:1,warm:rng()<0.12});
            }
            return out;
        })();
        const indicesOf=kind=>particles.map((p,i)=>p.kind===kind?i:-1).filter(i=>i>=0);
        const bulgeIndices=indicesOf('bulge'),diskIndices=indicesOf('disk'),haloIndices=indicesOf('halo');
        const vObs=DATA.rotation.map(p=>p.observed_kms);
        function observedVelocity(radius){
            const last=radii.length-1;
            if(radius<=radii[0])return vObs[0]*radius/radii[0];
            if(radius>=radii[last])return vObs[last];
            let low=0,high=last;
            while(high-low>1){const m=(low+high)>>1;if(radii[m]<=radius)low=m;else high=m;}
            return vObs[low]+(radius-radii[low])/(radii[high]-radii[low])*(vObs[high]-vObs[low]);
        }
        const radii=DATA.rotation.map(p=>p.radius_kpc);
        const vBary=DATA.rotation.map(p=>p.baryonic_kms);
        function nfwVelocityModel(massSolar,scaleRadius) {
            if(massSolar<=0)return ()=>0;
            const rhoCrit=3*H0*H0/(8*Math.PI*G);
            const r200=Math.cbrt(3*massSolar/(4*Math.PI*200*rhoCrit));
            const shape=x=>x<=0?0:Math.log1p(x)-x/(1+x);
            const fc=Math.max(shape(r200/scaleRadius),1e-12);
            return radius=>{
                if(radius<=0)return 0;
                const enclosed=radius>=r200?massSolar:massSolar*shape(radius/scaleRadius)/fc;
                return Math.sqrt(G*enclosed/radius);
            };
        }
        function haloCurve(massT,scaleRadius) {return radii.map(nfwVelocityModel(massT*1e12,scaleRadius));}
        function baryonicVelocity(radius) {
            if(radius<=radii[0])return vBary[0]*radius/radii[0];
            const last=radii.length-1;
            if(radius>=radii[last])return vBary[last]*Math.sqrt(radii[last]/radius);
            let low=0,high=last;
            while(high-low>1){const middle=Math.floor((low+high)/2);if(radii[middle]<=radius)low=middle;else high=middle;}
            const fraction=(radius-radii[low])/(radii[high]-radii[low]);
            return vBary[low]+fraction*(vBary[high]-vBary[low]);
        }
        function orbitDerivative(state,haloVelocity) {
            const [x,y,z,vx,vy,vz]=state;
            const radius=Math.max(Math.hypot(x,y,z),1e-5);
            const baryonic=baryonicVelocity(radius);
            const halo=haloVelocity(radius);
            const k=(baryonic*baryonic+halo*halo)*KM_S_TO_KPC_MYR**2/(radius*radius);
            return [vx,vy,vz,-k*x,-k*y,-k*z];
        }
        function integrateTrajectories(massT,scaleRadius) {
            const haloVelocity=nfwVelocityModel(massT*1e12,scaleRadius),step=DATA.frame_step_myr;
            return particles.map(particle=>{
                const [x,y,z]=particle.position;
                let dx,dy,dz;
                if(particle.kind!=='disk'){
                    const a=particle.axis;
                    dx=a[1]*z-a[2]*y;dy=a[2]*x-a[0]*z;dz=a[0]*y-a[1]*x;
                    if(Math.hypot(dx,dy,dz)<1e-4){
                        const alt=[a[1],a[2],a[0]];
                        dx=alt[1]*z-alt[2]*y;dy=alt[2]*x-alt[0]*z;dz=alt[0]*y-alt[1]*x;
                    }
                }else{dx=-y;dy=x;dz=0;}
                const norm=Math.max(Math.hypot(dx,dy,dz),1e-8);
                const speed=particle.speedScale*observedVelocity(particle.radius)*KM_S_TO_KPC_MYR/norm;
                let state=[x,y,z,dx*speed,dy*speed,dz*speed];
                const positions=[[x,y,z]];
                for(let frame=1;frame<frameCount;frame++){
                      const k1=orbitDerivative(state,haloVelocity);
                      const k2=orbitDerivative(state.map((value,index)=>value+step*k1[index]/2),haloVelocity);
                      const k3=orbitDerivative(state.map((value,index)=>value+step*k2[index]/2),haloVelocity);
                      const k4=orbitDerivative(state.map((value,index)=>value+step*k3[index]),haloVelocity);
                    state=state.map((value,index)=>value+step*(k1[index]+2*k2[index]+2*k3[index]+k4[index])/6);
                    positions.push([state[0],state[1],state[2]]);
                }
                return positions;
            });
        }
        let activeTrajectories=integrateTrajectories(DEFAULTS.massT,DEFAULTS.scaleRadius);
        let currentFrame=0,orbitModelTimer=null;
    const palette = [[0,[35,145,255]],[0.5,[255,211,69]],[1,[255,69,58]]];
    function velocityColor(v, radius) {
      const t = Math.max(0,Math.min(1,(v-minSpeed)/speedSpan));
      let a=palette[0],b=palette[1];
      if(t>0.5){a=palette[1];b=palette[2];}
      const f=(t-a[0])/(b[0]-a[0]);
      const distanceLift=0.10*Math.max(0,Math.min(1,radius/maxRadius));
      const rgb=a[1].map((c,i)=>Math.round((c+f*(b[1][i]-c))*(1-distanceLift)+255*distanceLift));
      return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
    }
        const orbitPlot=document.getElementById('orbitPlot');

        /* ---------- Camera / input module: damped orbit controls (mouse, pen, touch) ---------- */
        class SmoothOrbitControls {
            // Z-up spherical orbit. Input only edits `goal`; update(dt) eases `view` toward it
            // with a frame-rate independent exponential filter, so nothing ever snaps.
            constructor(camera,element,options={}) {
                this.camera=camera;this.element=element;this.dragMode='rotate';
                this.o=Object.assign({damping:11,rotateSpeed:1,panSpeed:1,zoomSpeed:0.0016,inertia:0.14,
                    minDistance:5,maxDistance:450,minPhi:0.04,maxPhi:Math.PI-0.04,panLimit:80,
                    initial:{theta:-0.82,phi:1.01,radius:72,target:[0,0,0]},onChange:()=>{}},options);
                const i=this.o.initial;
                this.goal={theta:i.theta,phi:i.phi,radius:i.radius,target:new THREE.Vector3(...i.target)};
                this.view={theta:i.theta,phi:i.phi,radius:i.radius,target:this.goal.target.clone()};
                this.pointers=new Map();this.pinch=null;this.lastMove=0;this.velocity={theta:0,phi:0};
                this._right=new THREE.Vector3();this._up=new THREE.Vector3();
                this._bind();this._apply();
            }
            _bind() {
                const el=this.element;
                this._h={
                    down:e=>this._down(e),move:e=>this._move(e),up:e=>this._up_(e),
                    wheel:e=>this._wheel(e),menu:e=>e.preventDefault(),dbl:()=>this.reset()
                };
                el.addEventListener('pointerdown',this._h.down);
                el.addEventListener('pointermove',this._h.move);
                el.addEventListener('pointerup',this._h.up);
                el.addEventListener('pointercancel',this._h.up);
                el.addEventListener('lostpointercapture',this._h.up);
                el.addEventListener('wheel',this._h.wheel,{passive:false});
                el.addEventListener('contextmenu',this._h.menu);
                el.addEventListener('dblclick',this._h.dbl);
            }
            dispose() {
                const el=this.element,h=this._h;
                el.removeEventListener('pointerdown',h.down);el.removeEventListener('pointermove',h.move);
                el.removeEventListener('pointerup',h.up);el.removeEventListener('pointercancel',h.up);
                el.removeEventListener('lostpointercapture',h.up);el.removeEventListener('wheel',h.wheel);
                el.removeEventListener('contextmenu',h.menu);el.removeEventListener('dblclick',h.dbl);
            }
            // Eases to a preset; theta takes the shortest way round so the camera never spins the long way.
            flyTo({theta,phi,radius,target}={}) {
                const g=this.goal;
                if(theta!==undefined){
                    const turn=Math.PI*2;
                    g.theta=theta+Math.round((this.view.theta-theta)/turn)*turn;
                }
                if(phi!==undefined)g.phi=this._clampPhi(phi);
                if(radius!==undefined)g.radius=this._clampRadius(radius);
                if(target)g.target.set(...target);
                this.velocity.theta=this.velocity.phi=0;
                this.o.onChange();
            }
            reset() {this.flyTo(this.o.initial);}
            _down(e) {
                if(e.pointerType==='mouse'&&e.button>2)return;
                try{this.element.setPointerCapture(e.pointerId);}catch(_){}
                this.pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,button:e.button,shift:e.shiftKey||e.ctrlKey||e.metaKey,type:e.pointerType});
                this.velocity.theta=this.velocity.phi=0;
                if(this.pointers.size===2)this.pinch=this._pinchState();
                this.element.classList.add('dragging');
                this.o.onChange();
            }
            _pinchState() {
                const [a,b]=[...this.pointers.values()];
                return{dist:Math.hypot(a.x-b.x,a.y-b.y)||1,x:(a.x+b.x)/2,y:(a.y+b.y)/2};
            }
            _move(e) {
                const p=this.pointers.get(e.pointerId);
                if(!p)return;
                // Coalesced events keep fast swipes smooth on high-polling-rate devices.
                const events=e.getCoalescedEvents?e.getCoalescedEvents():[];
                for(const ev of (events.length?events:[e])) {
                    const dx=ev.clientX-p.x,dy=ev.clientY-p.y;
                    p.x=ev.clientX;p.y=ev.clientY;
                    if(this.pointers.size===1) {
                        const mousePan=p.type!=='touch'&&p.button===0&&this.dragMode==='pan';
                        if(p.button===2||p.button===1||p.shift||mousePan)this._pan(dx,dy);
                        else this._rotate(dx,dy,ev.timeStamp);
                    } else if(this.pointers.size===2&&this.pinch) {
                        const s=this._pinchState();
                        this.goal.radius=this._clampRadius(this.goal.radius*this.pinch.dist/s.dist);
                        this._pan(s.x-this.pinch.x,s.y-this.pinch.y);
                        this.pinch=s;
                    }
                }
            }
            _up_(e) {
                if(!this.pointers.delete(e.pointerId))return;
                try{this.element.releasePointerCapture(e.pointerId);}catch(_){}
                this.pinch=null;
                if(this.pointers.size===1)this.pinch=null;
                if(this.pointers.size===0) {
                    this.element.classList.remove('dragging');
                    // Release momentum: only if the pointer was still moving when lifted.
                    if(performance.now()-this.lastMove<80) {
                        const k=this.o.inertia*1000;
                        this.goal.theta+=this.velocity.theta*k;
                        this.goal.phi=this._clampPhi(this.goal.phi+this.velocity.phi*k);
                    }
                    this.velocity.theta=this.velocity.phi=0;
                } else if(this.pointers.size===1&&this.pinch===null) {
                    // Remaining finger continues as a rotate without a jump.
                    const rest=[...this.pointers.values()][0];rest.button=0;rest.shift=false;
                }
                this.o.onChange();
            }
            _height() {return this.element.clientHeight||1;}
            _rotate(dx,dy,t) {
                const f=2*Math.PI*this.o.rotateSpeed/this._height();
                const dTheta=-dx*f,dPhi=-dy*f;
                this.goal.theta+=dTheta;
                this.goal.phi=this._clampPhi(this.goal.phi+dPhi);
                const now=performance.now(),dt=Math.max(now-this.lastMove,1);
                this.velocity.theta=0.6*this.velocity.theta+0.4*dTheta/dt;
                this.velocity.phi=0.6*this.velocity.phi+0.4*dPhi/dt;
                this.lastMove=now;
                this.o.onChange();
            }
            _pan(dx,dy) {
                const cam=this.camera,k=2*this.goal.radius*Math.tan(THREE.MathUtils.degToRad(cam.fov)/2)/this._height()*this.o.panSpeed;
                cam.updateMatrixWorld();
                this._right.setFromMatrixColumn(cam.matrixWorld,0);
                this._up.setFromMatrixColumn(cam.matrixWorld,1);
                this.goal.target.addScaledVector(this._right,-dx*k).addScaledVector(this._up,dy*k);
                const L=this.o.panLimit;
                this.goal.target.set(...this.goal.target.toArray().map(v=>Math.max(-L,Math.min(L,v))));
                this.o.onChange();
            }
            _wheel(e) {
                e.preventDefault();
                const unit=e.deltaMode===1?16:e.deltaMode===2?100:1;
                const delta=Math.max(-400,Math.min(400,e.deltaY*unit))*(e.ctrlKey?2.5:1);
                this.goal.radius=this._clampRadius(this.goal.radius*Math.exp(delta*this.o.zoomSpeed));
                this.o.onChange();
            }
            _clampPhi(v) {return Math.max(this.o.minPhi,Math.min(this.o.maxPhi,v));}
            _clampRadius(v) {return Math.max(this.o.minDistance,Math.min(this.o.maxDistance,v));}
            // Returns true while the view is still easing toward the goal (caller keeps rendering).
            update(dt) {
                const g=this.goal,v=this.view,a=1-Math.exp(-this.o.damping*Math.min(dt,0.1));
                const dTheta=g.theta-v.theta,dPhi=g.phi-v.phi,dR=g.radius-v.radius;
                const dT=g.target.distanceTo(v.target);
                if(Math.abs(dTheta)<1e-5&&Math.abs(dPhi)<1e-5&&Math.abs(dR)<1e-3&&dT<1e-3) {
                    if(dTheta||dPhi||dR||dT){v.theta=g.theta;v.phi=g.phi;v.radius=g.radius;v.target.copy(g.target);this._apply();return true;}
                    return false;
                }
                v.theta+=dTheta*a;v.phi+=dPhi*a;v.radius+=dR*a;v.target.lerp(g.target,a);
                this._apply();
                return true;
            }
            _apply() {
                const v=this.view,s=Math.sin(v.phi);
                this.camera.up.set(0,0,1);
                this.camera.position.set(v.target.x+v.radius*s*Math.cos(v.theta),v.target.y+v.radius*s*Math.sin(v.theta),v.target.z+v.radius*Math.cos(v.phi));
                this.camera.lookAt(v.target);
            }
        }

        /* ---------- Three.js scene ---------- */
        const renderer=new THREE.WebGLRenderer({antialias:true,powerPreference:'high-performance'});
        renderer.setClearColor(BG,1);
        orbitPlot.appendChild(renderer.domElement);
        const scene=new THREE.Scene();
        const camera=new THREE.PerspectiveCamera(40,1,1,2000);
        // Three-quarter view: ~39 degrees above the disk plane, framed to fit the galaxy.
        const GALAXY_VIEW=(()=>{const e=[0.55,-1.0,0.8],r=Math.hypot(...e);return{theta:Math.atan2(e[1],e[0]),phi:Math.acos(e[2]/r),radius:r*40,target:[0,0,0]};})();
        const diskCount=diskIndices.length;

        // Disk stars are InstancedMesh spheres; bulge and halo are soft additive glow sprites.
        const STAR_BASE_RADIUS=0.32*0.75;
        const starCount=diskIndices.length;
        const starIndices=diskIndices;
        const starMesh=new THREE.InstancedMesh(new THREE.SphereGeometry(STAR_BASE_RADIUS,10,8),new THREE.MeshBasicMaterial({toneMapped:false}),starCount);
        starMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
        starMesh.frustumCulled=false;starMesh.renderOrder=1;
        (()=>{
            const m=new THREE.Matrix4(),c=new THREE.Color();
            const white=new THREE.Color('#FFFFFF'),blue=new THREE.Color('#9FC5E8');
            const smooth01=x=>{const t=Math.min(1,Math.max(0,x));return t*t*(3-2*t);};
            for(let k=0;k<starCount;k++){
                const pt=particles[starIndices[k]],p=activeTrajectories[starIndices[k]][0];
                c.copy(white).lerp(blue,smooth01((pt.radius-GALAXY.diskMin)/(maxRadius*0.6)));
                const sc=0.40+0.15*(pt.radius/maxRadius);
                m.makeScale(sc,sc,sc).setPosition(p[0],p[1],p[2]);
                starMesh.setMatrixAt(k,m);
                starMesh.setColorAt(k,c);
            }
            starMesh.instanceColor.needsUpdate=true;
        })();
        scene.add(starMesh);

        // Soft glow sprites (THREE.Points + radial-gradient texture + additive blending) for bulge and halo.
        const spriteTexture=(()=>{
            const cv=document.createElement('canvas');cv.width=cv.height=64;
            const g=cv.getContext('2d'),gr=g.createRadialGradient(32,32,0,32,32,32);
            gr.addColorStop(0,'rgba(255,255,255,1)');gr.addColorStop(0.25,'rgba(255,255,255,0.55)');
            gr.addColorStop(0.6,'rgba(255,255,255,0.12)');gr.addColorStop(1,'rgba(255,255,255,0)');
            g.fillStyle=gr;g.fillRect(0,0,64,64);
            return new THREE.CanvasTexture(cv);
        })();
        const spriteUniforms={uScale:{value:400}};
        const spriteMaterial=new THREE.ShaderMaterial({
            uniforms:{...spriteUniforms,uMap:{value:spriteTexture}},
            transparent:true,depthWrite:false,blending:THREE.AdditiveBlending,vertexColors:true,
            vertexShader:`attribute float size;varying vec3 vC;uniform float uScale;
                void main(){vC=color;vec4 mv=modelViewMatrix*vec4(position,1.0);gl_PointSize=size*uScale/-mv.z;gl_Position=projectionMatrix*mv;}`,
            fragmentShader:`uniform sampler2D uMap;varying vec3 vC;
                void main(){vec4 t=texture2D(uMap,gl_PointCoord);gl_FragColor=vec4(vC*t.a,t.a);}`
        });
        spriteUniforms.uScale=spriteMaterial.uniforms.uScale;
        const spriteSets=[];
        function makeSpritePoints(indices,colorFn,sizeFn,order){
            const n=indices.length,pos=new Float32Array(n*3),col=new Float32Array(n*3),size=new Float32Array(n);
            const c=new THREE.Color();
            for(let k=0;k<n;k++){
                const pt=particles[indices[k]];colorFn(pt,c);
                col[k*3]=c.r;col[k*3+1]=c.g;col[k*3+2]=c.b;size[k]=sizeFn(pt);
            }
            const geom=new THREE.BufferGeometry();
            const posAttr=new THREE.BufferAttribute(pos,3);posAttr.setUsage(THREE.DynamicDrawUsage);
            geom.setAttribute('position',posAttr);
            geom.setAttribute('color',new THREE.BufferAttribute(col,3));
            geom.setAttribute('size',new THREE.BufferAttribute(size,1));
            const pts=new THREE.Points(geom,spriteMaterial);
            pts.frustumCulled=false;pts.renderOrder=order;scene.add(pts);
            spriteSets.push({indices,pos,posAttr});
        }
        const bulgeIn=new THREE.Color('#FFFFFF'),bulgeOut=new THREE.Color('#FF8FA3');
        makeSpritePoints(bulgeIndices,(pt,c)=>{
            const t=Math.min(1,pt.radius/(GALAXY.bulgeMax*0.6));c.copy(bulgeIn).lerp(bulgeOut,Math.sqrt(t)).            multiplyScalar(0.42-0.22*t);
                    },pt=>3.2+1.6*Math.min(1,pt.radius/(GALAXY.bulgeMax*0.6)),0);
        makeSpritePoints(haloIndices,(pt,c)=>{
            c.set(pt.warm?'#ffc68a':'#9db8e0').multiplyScalar(pt.warm?0.7:0.35);
        },pt=>pt.warm?2.4:1.3,0);
        // Black hole: camera-facing shader quad. Opaque shadow disc, thin Doppler-beamed photon ring and a soft accretion glow.
        const BH_R=GALAXY.blackHoleRadius,BH_EXTENT=4.0;
        const blackHole=new THREE.Mesh(
            new THREE.PlaneGeometry(2*BH_R*BH_EXTENT,2*BH_R*BH_EXTENT),
            new THREE.ShaderMaterial({
                transparent:true,depthWrite:false,blending:THREE.CustomBlending,
                blendEquation:THREE.AddEquation,blendSrc:THREE.OneFactor,blendDst:THREE.OneMinusSrcAlphaFactor,
                vertexShader:`varying vec2 vP;void main(){vP=position.xy;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}`,
                fragmentShader:`
                    varying vec2 vP;
                    void main(){
                        vec2 p=vP/${BH_R.toFixed(3)};float d=length(p);
                        float shadow=1.0-smoothstep(0.96,1.0,d);
                        float beam=1.0+0.55*cos(atan(p.y,p.x)-0.6);
                        float ring=exp(-pow((d-1.07)/0.045,2.0))*1.9*beam;
                        float glow=exp(-pow((d-1.12)/0.42,2.0))*0.32;
                        float fade=1.0-smoothstep(${(BH_EXTENT-1.0).toFixed(1)},${BH_EXTENT.toFixed(1)},d);
                        vec3 col=(vec3(1.0,0.97,0.88)*ring+vec3(1.0,0.62,0.32)*glow)*fade*(1.0-shadow);
                        gl_FragColor=vec4(col,shadow);
                    }
                `
            })
        );
        blackHole.position.set(0,0,0);
        blackHole.renderOrder=3;scene.add(blackHole);

        // Optional orbit trajectory lines: full paths of an evenly spaced subset of disk stars.
        const ORBIT_LINE_STARS=90,ORBIT_LINE_STEP=3;
        const orbitLineIdx=Array.from({length:Math.min(ORBIT_LINE_STARS,diskCount)},(_,k)=>diskIndices[Math.floor(k*diskCount/Math.min(ORBIT_LINE_STARS,diskCount))]);
        const orbitSegs=Math.ceil((frameCount-1)/ORBIT_LINE_STEP);
        const orbitLineGeom=new THREE.BufferGeometry();
        orbitLineGeom.setAttribute('position',new THREE.BufferAttribute(new Float32Array(orbitLineIdx.length*orbitSegs*6),3));
        const orbitLines=new THREE.LineSegments(orbitLineGeom,new THREE.LineBasicMaterial({color:0x7be4f1,transparent:true,opacity:0.22,depthWrite:false,blending:THREE.AdditiveBlending}));
        orbitLines.frustumCulled=false;orbitLines.visible=false;orbitLines.renderOrder=0;scene.add(orbitLines);
        function fillOrbitLines(){
            const a=orbitLineGeom.attributes.position.array;let o=0;
            for(const idx of orbitLineIdx){
                const t=activeTrajectories[idx];
                for(let s=0;s<orbitSegs;s++){
                    const p=t[s*ORBIT_LINE_STEP],q=t[Math.min((s+1)*ORBIT_LINE_STEP,frameCount-1)];
                    a[o++]=p[0];a[o++]=p[1];a[o++]=p[2];a[o++]=q[0];a[o++]=q[1];a[o++]=q[2];
                }
            }
            orbitLineGeom.attributes.position.needsUpdate=true;
        }
        fillOrbitLines();

        // Sub-frame continuous orbit updates for locked 60 FPS animation
        function updateStarTransforms(timeFloat){
            const fFloat=Math.max(0,Math.min(frameCount-1,timeFloat));
            const f0=Math.floor(fFloat);
            const f1=Math.min(f0+1,frameCount-1);
            const alpha=fFloat-f0;
            const matArray=starMesh.instanceMatrix.array;
            for(let k=0;k<starCount;k++){
                const trajectory=activeTrajectories[starIndices[k]];
                const p0=trajectory[f0],p1=trajectory[f1],off=k*16;
                matArray[off+12]=p0[0]+(p1[0]-p0[0])*alpha;
                matArray[off+13]=p0[1]+(p1[1]-p0[1])*alpha;
                matArray[off+14]=p0[2]+(p1[2]-p0[2])*alpha;
            }
            starMesh.instanceMatrix.needsUpdate=true;
            for(const set of spriteSets){
                for(let k=0;k<set.indices.length;k++){
                    const trajectory=activeTrajectories[set.indices[k]];
                    const p0=trajectory[f0],p1=trajectory[f1];
                    set.pos[k*3]=p0[0]+(p1[0]-p0[0])*alpha;
                    set.pos[k*3+1]=p0[1]+(p1[1]-p0[1])*alpha;
                    set.pos[k*3+2]=p0[2]+(p1[2]-p0[2])*alpha;
                }
                set.posAttr.needsUpdate=true;
            }
        }

        // Render loop: continuous delta-time playback with throttled UI sync
        let needsRender=true,lastTime=performance.now(),simTime=0,playing=false,lastReportedFrame=-1;
        const PLAYBACK_SPEED_MYR_PER_SEC=24.0;
        const controls=new SmoothOrbitControls(camera,orbitPlot,{initial:GALAXY_VIEW,onChange:()=>{needsRender=true;}});
        function requestRender(){needsRender=true;}
        // Panel clicks must not reach the controls' drag / wheel / double-click handlers.
        const viewPanel=document.getElementById('viewPanel');
        ['pointerdown','pointermove','pointerup','wheel','dblclick','contextmenu'].forEach(type=>viewPanel.addEventListener(type,e=>e.stopPropagation(),{passive:true}));
        const viewToggle=document.getElementById('viewToggle');
        if(window.matchMedia('(hover:none) and (pointer:coarse)').matches){
            viewPanel.classList.add('collapsed');
            viewToggle.setAttribute('aria-expanded','false');
        }
        viewToggle.addEventListener('click',()=>{
            const collapsed=viewPanel.classList.toggle('collapsed');
            viewToggle.setAttribute('aria-expanded',String(!collapsed));
        });
        function resetView(){controls.reset();requestRender();}
        document.getElementById('viewReset').addEventListener('click',resetView);
        const orbitToggle=document.getElementById('toggleOrbits');
        orbitToggle.addEventListener('click',()=>{
            orbitLines.visible=!orbitLines.visible;
            orbitToggle.setAttribute('aria-pressed',String(orbitLines.visible));
            requestRender();
        });
        document.getElementById('viewFaceOn').addEventListener('click',()=>controls.flyTo({phi:0.04,target:[0,0,0],radius:GALAXY_VIEW.radius}));
        document.getElementById('viewEdgeOn').addEventListener('click',()=>controls.flyTo({phi:Math.PI/2,target:[0,0,0],radius:GALAXY_VIEW.radius}));
        // Rotate/Pan only changes what a primary mouse drag does; touch gestures are unaffected.
        const modeButtons=[...viewPanel.querySelectorAll('[data-mode]')];
        function setDragMode(mode){
            controls.dragMode=mode;
            orbitPlot.classList.toggle('pan-mode',mode==='pan');
            modeButtons.forEach(b=>b.setAttribute('aria-checked',String(b.dataset.mode===mode)));
        }
        modeButtons.forEach(b=>b.addEventListener('click',()=>setDragMode(b.dataset.mode)));
        function resizeOrbit(){
            const w=Math.max(orbitPlot.clientWidth,1),h=Math.max(orbitPlot.clientHeight,1),dpr=Math.min(window.devicePixelRatio||1,2);
            renderer.setPixelRatio(dpr);renderer.setSize(w,h,false);
            camera.aspect=w/h;camera.updateProjectionMatrix();
            spriteUniforms.uScale.value=h*dpr/(2*Math.tan(THREE.MathUtils.degToRad(camera.fov/2)));
            needsRender=true;
        }
        if(window.ResizeObserver)new ResizeObserver(resizeOrbit).observe(orbitPlot);
        window.addEventListener('resize',resizeOrbit);
        function renderLoop(now){
            const dt=Math.min((now-lastTime)/1000,0.1);
            lastTime=now;
            if(controls.update(dt))needsRender=true;
            if(playing){
                simTime+=dt*PLAYBACK_SPEED_MYR_PER_SEC;
                if(simTime>=frameCount-1){
                    simTime=frameCount-1;
                    playing=false;
                    playButton.textContent='Play';
                }
                updateStarTransforms(simTime);
                const intFrame=Math.floor(simTime);
                if(intFrame!==lastReportedFrame){
                    lastReportedFrame=intFrame;
                    syncTimelineUI(intFrame);
                }
                needsRender=true;
            }
            if(needsRender){needsRender=false;blackHole.quaternion.copy(camera.quaternion);renderer.render(scene,camera);}
            requestAnimationFrame(renderLoop);
        }
        resizeOrbit();requestAnimationFrame(renderLoop);

    const initialHalo=haloCurve(DEFAULTS.massT,DEFAULTS.scaleRadius);
    const totalCurve=vBary.map((v,i)=>Math.hypot(v,initialHalo[i]));
    const curveData=[
      {type:'scatter',mode:'markers',name:'Observed SPARC',x:radii,y:DATA.rotation.map(p=>p.observed_kms),error_y:{type:'data',array:DATA.rotation.map(p=>p.error_kms),visible:true,color:'rgba(225,235,245,0.42)',thickness:1,width:2},marker:{size:6,color:'#e8edf2',line:{color:BG,width:1}},hovertemplate:'R %{x:.2f} kpc<br>Observed %{y:.1f} km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Baryonic model',x:radii,y:vBary,line:{color:'#60bfd0',width:2.5},hovertemplate:'R %{x:.2f} kpc<br>Baryonic %{y:.1f} km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Dark matter contribution',x:radii,y:initialHalo,line:{color:'#b39ae8',width:2,dash:'dot'},hovertemplate:'R %{x:.2f} kpc<br>Halo %{y:.1f} km/s<extra></extra>'},
      {type:'scatter',mode:'lines',name:'Baryons + NFW halo',x:radii,y:totalCurve,line:{color:'#efa766',width:2.5},hovertemplate:'R %{x:.2f} kpc<br>Total %{y:.1f} km/s<extra></extra>'}
    ];
    const curveLayout={paper_bgcolor:BG,plot_bgcolor:BG,margin:{l:64,r:18,t:12,b:56},font:{color:textColor,size:11},xaxis:{title:'Radius (kpc)',gridcolor:GRID,zerolinecolor:GRID},yaxis:{title:'Orbital velocity (km/s)',range:[0,190],gridcolor:GRID,zerolinecolor:GRID},legend:{orientation:'h',y:1.12,x:0,font:{size:10}},hovermode:'x unified'};
    Plotly.newPlot('curvePlot',curveData,curveLayout,{responsive:true,displaylogo:false});

    const massSlider=document.getElementById('massSlider'),scaleSlider=document.getElementById('scaleSlider');
    massSlider.value=DEFAULTS.massT;scaleSlider.value=DEFAULTS.scaleRadius;
        function updateHalo(){
      const massT=Number(massSlider.value),rs=Number(scaleSlider.value),halo=haloCurve(massT,rs);
      const total=vBary.map((v,i)=>Math.hypot(v,halo[i]));
      document.getElementById('massValue').textContent=`${massT.toFixed(2)} × 10¹² M☉`;
      document.getElementById('scaleValue').textContent=`${rs.toFixed(1)} kpc`;
      Plotly.restyle('curvePlot',{y:[halo,total]},[2,3]);
            clearTimeout(orbitModelTimer);
            orbitModelTimer=setTimeout(()=>{
                activeTrajectories=integrateTrajectories(massT,rs);
                fillOrbitLines();
                updateStarTransforms(simTime);requestRender();
            },120);
    }
    massSlider.addEventListener('input',updateHalo);scaleSlider.addEventListener('input',updateHalo);updateHalo();
        document.getElementById('applyObservedFit').addEventListener('click',()=>{
            massSlider.value=DEFAULTS.massT;
            scaleSlider.value=DEFAULTS.scaleRadius;
            updateHalo();
        });

    const timeSlider=document.getElementById('timeSlider'),timeSliderBottom=document.getElementById('timeSliderBottom');
    const playButton=document.getElementById('playButton');
    function syncTimelineUI(intFrame){
        timeSlider.value=intFrame;timeSliderBottom.value=intFrame;
        document.getElementById('frameLabel').textContent=`Frame ${intFrame+1}`;
        document.getElementById('timeValue').textContent=`${(intFrame*DATA.frame_step_myr).toFixed(0)} Myr`;
        document.getElementById('timeLabel').textContent=`${(intFrame*DATA.frame_step_myr).toFixed(0)} Myr`;
    }
    function onSliderInput(val){
        simTime=Math.max(0,Math.min(frameCount-1,Number(val)));
        updateStarTransforms(simTime);
        const intFrame=Math.floor(simTime);
        lastReportedFrame=intFrame;
        syncTimelineUI(intFrame);
        requestRender();
    }
    timeSlider.addEventListener('input',()=>onSliderInput(timeSlider.value));
    timeSliderBottom.addEventListener('input',()=>onSliderInput(timeSliderBottom.value));
    playButton.addEventListener('click',()=>{
        if(playing){
            playing=false;playButton.textContent='Play';
        }else{
            if(simTime>=frameCount-1)simTime=0;
            playing=true;playButton.textContent='Pause';
        }
    });
    updateStarTransforms(0);syncTimelineUI(0);
  </script>
</body>
</html>
"""


def main() -> None:
    if not TRAJECTORY_PATH.is_file():
        raise FileNotFoundError(f"Missing trajectory CSV: {TRAJECTORY_PATH}")
    if not SPARC_PATH.is_file():
        raise FileNotFoundError(f"Missing SPARC table: {SPARC_PATH}")

    stars = load_trajectories(TRAJECTORY_PATH)
    rotation = load_rotation_curve(SPARC_PATH)
    mass_msun, scale_radius = fit_default_halo(rotation)
    payload = {
        "stars": stars,
        "rotation": rotation,
        "frame_step_myr": 1.0,
        "default_halo_mass_msun": mass_msun,
        "default_scale_radius_kpc": scale_radius,
    }
    defaults = {"massT": mass_msun / 1e12, "scaleRadius": scale_radius}
    html = HTML_TEMPLATE.replace(
        "__DATA_JSON__",
        json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/"),
    )
    html = html.replace("__DEFAULT_JSON__", json.dumps(defaults))
    html = html.replace("__FRAME_MAX__", str(len(stars[0]["positions"]) - 1))
    OUTPUT_PATH.write_text(html, encoding="utf-8")
    size_mb = OUTPUT_PATH.stat().st_size / (1024 * 1024)
    print(f"Generated {OUTPUT_PATH.name} ({size_mb:.2f} MiB)")
    print(f"Embedded {len(stars)} stars and {len(rotation)} NGC 3198 SPARC points")
    print(f"NFW starting fit: M200={mass_msun:.3g} Msun, rs={scale_radius:.1f} kpc")


if __name__ == "__main__":
    main()