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
    <div class="header-note">200 massless test stars · 301 frames · 1 Myr per frame<br>Local data embedded in this page · dark halo uses an NFW profile</div>
  </header>
  <main>
    <section class="controls" aria-label="Dark matter controls">
    <div class="control"><label for="massSlider"><span>Halo mass, M<sub>200</sub></span><output id="massValue"></output></label><input id="massSlider" type="range" min="0.00" max="2.00" step="0.01"><p class="control-help">Updates the rotation curve and star orbits. Set it to 0 for visible matter only.</p><button id="applyObservedFit" class="fit-button" type="button">Apply Observed Fit</button></div>
    <div class="control"><label for="scaleSlider"><span>Scale radius, r<sub>s</sub></span><output id="scaleValue"></output></label><input id="scaleSlider" type="range" min="5" max="60" step="0.5"><p class="control-help">Controls halo concentration: a smaller radius packs more dark matter near the center.</p></div>
      <div class="control"><label for="timeSlider"><span>Orbit time</span><output id="timeValue"></output></label><div style="display:flex;align-items:center;gap:12px"><input id="timeSlider" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"><button id="playButton" type="button">Play</button></div></div>
    </section>
    <section class="panel-grid">
    <article class="panel"><div class="panel-heading"><h2>Stellar orbit trajectories</h2></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional NGC 3198 star trajectories"><div class="view-panel" id="viewPanel" role="group" aria-label="3D view controls"><button type="button" class="view-head" id="viewToggle" aria-expanded="true" title="Collapse or expand the view controls"><span>View</span><span class="chev" aria-hidden="true">▾</span></button><button type="button" id="viewReset" title="Return to the default NGC 3198 view">Reset view</button><div class="view-row"><button type="button" id="viewFaceOn" title="Look straight down on the galactic plane">Face-on</button><button type="button" id="viewEdgeOn" title="Look along the galactic plane">Edge-on</button></div><button type="button" id="viewHighlight" aria-pressed="false" title="Brighten the orbit trajectory lines so they are easier to see">Highlight trajectories</button><div class="mode-toggle" role="radiogroup" aria-label="Mouse drag mode"><button type="button" role="radio" aria-checked="true" data-mode="rotate">Rotate</button><button type="button" role="radio" aria-checked="false" data-mode="pan">Pan</button></div></div></div><div class="timeline"><span id="frameLabel">Frame 1</span><input id="timeSliderBottom" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"><span id="timeLabel">0 Myr</span></div></article>
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
        const ORBIT_PATH_STEP = 6;
        const KM_S_TO_KPC_MYR = 0.0010227121650537077;
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
            const [x,y,vx,vy]=state;
            const radius=Math.max(Math.hypot(x,y),1e-5);
            const baryonic=baryonicVelocity(radius);
            const halo=haloVelocity(radius);
            const acceleration=(baryonic*baryonic+halo*halo)*KM_S_TO_KPC_MYR**2/radius;
            return [vx,vy,-acceleration*x/radius,-acceleration*y/radius];
        }
        function integrateTrajectories(massT,scaleRadius) {
            const haloVelocity=nfwVelocityModel(massT*1e12,scaleRadius),step=DATA.frame_step_myr;
            return stars.map(star=>{
                const [x,y]=star.positions[0],radius=Math.max(Math.hypot(x,y),1e-8);
                const speed=star.v_rot_kms*KM_S_TO_KPC_MYR;
                let state=[x,y,-speed*y/radius,speed*x/radius];
                const positions=[[x,y,0]];
                for(let frame=1;frame<frameCount;frame++){
                      const k1=orbitDerivative(state,haloVelocity);
                      const k2=orbitDerivative(state.map((value,index)=>value+step*k1[index]/2),haloVelocity);
                      const k3=orbitDerivative(state.map((value,index)=>value+step*k2[index]/2),haloVelocity);
                      const k4=orbitDerivative(state.map((value,index)=>value+step*k3[index]),haloVelocity);
                    state=state.map((value,index)=>value+step*(k1[index]+2*k2[index]+2*k3[index]+k4[index])/6);
                    positions.push([state[0],state[1],0]);
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
                    minDistance:12,maxDistance:420,minPhi:0.04,maxPhi:Math.PI-0.04,panLimit:80,
                    initial:{theta:-0.82,phi:1.0,radius:120,target:[0,0,0]},onChange:()=>{}},options);
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
        const GALAXY_VIEW=(()=>{const e=[1.35,-1.45,1.05],r=Math.hypot(...e);return{theta:Math.atan2(e[1],e[0]),phi:Math.acos(e[2]/r),radius:r*58,target:[0,0,0]};})();

        // Reference grid (galaxy plane) and axes
        const grid=new THREE.GridHelper(130,13,0x2a3b4c,0x1c2936);grid.rotation.x=Math.PI/2;scene.add(grid);
        const axisPoints=[0,0,0,65,0,0, 0,0,0,0,65,0, 0,0,0,0,0,65];
        const axisColors=[0.48,0.89,0.95,0.48,0.89,0.95, 0.95,0.78,0.47,0.95,0.78,0.47, 0.6,0.6,0.9,0.6,0.6,0.9];
        const axisGeometry=new THREE.BufferGeometry();
        axisGeometry.setAttribute('position',new THREE.Float32BufferAttribute(axisPoints,3));
        axisGeometry.setAttribute('color',new THREE.Float32BufferAttribute(axisColors,3));
        scene.add(new THREE.LineSegments(axisGeometry,new THREE.LineBasicMaterial({vertexColors:true,transparent:true,opacity:0.55})));

        // Orbit paths (one segment list, rewritten in place when the halo changes)
        function orbitPathSegmentCount(){return stars.length*Math.floor((frameCount-1)/ORBIT_PATH_STEP);}
        const pathPositions=new Float32Array(orbitPathSegmentCount()*6);
        const pathAttribute=new THREE.BufferAttribute(pathPositions,3);pathAttribute.setUsage(THREE.DynamicDrawUsage);
        const pathGeometry=new THREE.BufferGeometry();pathGeometry.setAttribute('position',pathAttribute);
        function fillOrbitPaths(){
            let o=0;
            for(const positions of activeTrajectories){
                for(let i=0;i+ORBIT_PATH_STEP<frameCount;i+=ORBIT_PATH_STEP){
                    const a=positions[i],b=positions[i+ORBIT_PATH_STEP];
                    pathPositions[o++]=a[0];pathPositions[o++]=a[1];pathPositions[o++]=a[2];
                    pathPositions[o++]=b[0];pathPositions[o++]=b[1];pathPositions[o++]=b[2];
                }
            }
            pathAttribute.needsUpdate=true;
        }
        fillOrbitPaths();
        const pathMaterial=new THREE.LineBasicMaterial({color:0x6e91b2,transparent:true,opacity:0.2,depthWrite:false});
        const pathLines=new THREE.LineSegments(pathGeometry,pathMaterial);
        pathLines.frustumCulled=false;scene.add(pathLines);

        // Stars: one Points draw call, round sprites via a tiny shader
        const starPositions=new Float32Array(stars.length*3),starColors=new Float32Array(stars.length*3),starSizes=new Float32Array(stars.length);
        stars.forEach((s,i)=>{
            const c=new THREE.Color(velocityColor(s.v_rot_kms,s.initial_r_kpc));
            starColors.set([c.r,c.g,c.b],i*3);
            starSizes[i]=6.5+1.5*s.initial_r_kpc/maxRadius;
        });
        const starAttribute=new THREE.BufferAttribute(starPositions,3);starAttribute.setUsage(THREE.DynamicDrawUsage);
        const starGeometry=new THREE.BufferGeometry();
        starGeometry.setAttribute('position',starAttribute);
        starGeometry.setAttribute('color',new THREE.BufferAttribute(starColors,3));
        starGeometry.setAttribute('size',new THREE.BufferAttribute(starSizes,1));
        const pointVertex='attribute float size;attribute vec3 color;varying vec3 vColor;uniform float scale;void main(){vColor=color;gl_PointSize=size*scale;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}';
        const pointFragment='varying vec3 vColor;uniform float rim;void main(){vec2 p=gl_PointCoord-0.5;float d=length(p);if(d>0.5)discard;vec3 col=mix(vColor,vec3(0.92,0.97,1.0),smoothstep(0.38,0.5,d)*rim);gl_FragColor=vec4(col,1.0);}';
        const starMaterial=new THREE.ShaderMaterial({vertexShader:pointVertex,fragmentShader:pointFragment,uniforms:{scale:{value:1},rim:{value:0.7}}});
        const starPoints=new THREE.Points(starGeometry,starMaterial);starPoints.frustumCulled=false;scene.add(starPoints);

        const coreGeometry=new THREE.BufferGeometry();
        coreGeometry.setAttribute('position',new THREE.Float32BufferAttribute([0,0,0],3));
        coreGeometry.setAttribute('color',new THREE.Float32BufferAttribute([0.94,0.84,0.61],3));
        coreGeometry.setAttribute('size',new THREE.Float32BufferAttribute([12],1));
        const coreMaterial=new THREE.ShaderMaterial({vertexShader:pointVertex,fragmentShader:pointFragment,uniforms:{scale:{value:1},rim:{value:0}}});
        scene.add(new THREE.Points(coreGeometry,coreMaterial));

        function starCoordinates(frameIndex){
            for(let i=0;i<activeTrajectories.length;i++){
                const p=activeTrajectories[i][frameIndex];
                starPositions[i*3]=p[0];starPositions[i*3+1]=p[1];starPositions[i*3+2]=p[2];
            }
            starAttribute.needsUpdate=true;
        }

        // Render loop: always cheap, draws only when the camera is easing or the scene changed.
        let needsRender=true,lastTime=performance.now();
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
        const highlightButton=document.getElementById('viewHighlight');
        highlightButton.addEventListener('click',()=>{
            const on=highlightButton.getAttribute('aria-pressed')!=='true';
            highlightButton.setAttribute('aria-pressed',String(on));
            pathMaterial.color.set(on?0x7be4f1:0x6e91b2);
            pathMaterial.opacity=on?0.85:0.2;
            pathMaterial.blending=on?THREE.AdditiveBlending:THREE.NormalBlending;
            pathMaterial.needsUpdate=true;
            requestRender();
        });
        document.getElementById('viewReset').addEventListener('click',()=>controls.reset());
        document.getElementById('viewFaceOn').addEventListener('click',()=>controls.flyTo({phi:0,target:[0,0,0],radius:GALAXY_VIEW.radius}));
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
            starMaterial.uniforms.scale.value=dpr;coreMaterial.uniforms.scale.value=dpr;
            needsRender=true;
        }
        if(window.ResizeObserver)new ResizeObserver(resizeOrbit).observe(orbitPlot);
        window.addEventListener('resize',resizeOrbit);
        function renderLoop(now){
            const dt=(now-lastTime)/1000;lastTime=now;
            if(controls.update(dt))needsRender=true;
            if(needsRender){needsRender=false;renderer.render(scene,camera);}
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
                fillOrbitPaths();starCoordinates(currentFrame);requestRender();
            },120);
    }
    massSlider.addEventListener('input',updateHalo);scaleSlider.addEventListener('input',updateHalo);updateHalo();
        document.getElementById('applyObservedFit').addEventListener('click',()=>{
            massSlider.value=DEFAULTS.massT;
            scaleSlider.value=DEFAULTS.scaleRadius;
            updateHalo();
        });

    const timeSlider=document.getElementById('timeSlider'),timeSliderBottom=document.getElementById('timeSliderBottom');
        const playButton=document.getElementById('playButton');let timer=null,playing=false;
        function setFrame(index){
      index=Math.max(0,Math.min(frameCount-1,Number(index)));
            currentFrame=index;
      timeSlider.value=index;timeSliderBottom.value=index;
      starCoordinates(index);requestRender();
      document.getElementById('frameLabel').textContent=`Frame ${index+1}`;
      document.getElementById('timeValue').textContent=`${(index*DATA.frame_step_myr).toFixed(0)} Myr`;
      document.getElementById('timeLabel').textContent=`${(index*DATA.frame_step_myr).toFixed(0)} Myr`;
            return Promise.resolve();
    }
    timeSlider.addEventListener('input',()=>setFrame(timeSlider.value));
    timeSliderBottom.addEventListener('input',()=>setFrame(timeSliderBottom.value));
        async function advanceFrame(){
            if(!playing)return;
            const next=Number(timeSlider.value)+1;
            if(next>=frameCount){playing=false;timer=null;playButton.textContent='Play';return;}
            await setFrame(next);
            if(playing)timer=setTimeout(advanceFrame,150);
        }
        playButton.addEventListener('click',async()=>{
            if(playing){playing=false;clearTimeout(timer);timer=null;playButton.textContent='Play';return;}
            if(Number(timeSlider.value)>=frameCount-1)await setFrame(0);
            playing=true;playButton.textContent='Pause';advanceFrame();
    });
    setFrame(0);
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