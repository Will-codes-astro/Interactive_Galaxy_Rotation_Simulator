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
    #orbitPlot { touch-action:none; }
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
    <article class="panel"><div class="panel-heading"><h2>Stellar orbit trajectories</h2><span>3D · kpc · drag to orbit · toolbar for Pan/Zoom · wheel to zoom</span></div><div id="orbitPlot" class="plot" role="img" aria-label="Interactive three-dimensional NGC 3198 star trajectories"></div><div class="timeline"><span id="frameLabel">Frame 1</span><input id="timeSliderBottom" type="range" min="0" max="__FRAME_MAX__" step="1" value="0"><span id="timeLabel">0 Myr</span></div></article>
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
        function starCoordinates(frameIndex) {
            return {
                x:activeTrajectories.map(positions=>positions[frameIndex][0]),
                y:activeTrajectories.map(positions=>positions[frameIndex][1]),
                z:activeTrajectories.map(positions=>positions[frameIndex][2])
            };
        }
        function starTrace(frameIndex) {
            return {
                type:'scatter3d', mode:'markers', name:'Stars',
                ...starCoordinates(frameIndex),
        customdata:stars.map(s=>[s.id,s.v_rot_kms,s.initial_r_kpc]),
        marker:{size:stars.map(s=>4.8+1.2*s.initial_r_kpc/maxRadius),color:stars.map(s=>velocityColor(s.v_rot_kms,s.initial_r_kpc)),opacity:0.98,line:{width:0.5,color:'rgba(235,248,255,0.72)'}},
        hovertemplate:'Star %{customdata[0]}<br>R = %{customdata[2]:.1f} kpc<br>v = %{customdata[1]:.1f} km/s<extra></extra>'
      };
    }
        function orbitPathCoordinates(){
            const x=[],y=[],z=[];
            for(const positions of activeTrajectories){for(let i=0;i<positions.length;i+=ORBIT_PATH_STEP){const p=positions[i];x.push(p[0]);y.push(p[1]);z.push(p[2]);}x.push(null);y.push(null);z.push(null);}
            return{x,y,z};
        }
        const orbitPaths=orbitPathCoordinates();
        const orbitPlot=document.getElementById('orbitPlot');
    const orbitData=[
    {type:'scatter3d',mode:'lines',name:'Orbit paths',x:orbitPaths.x,y:orbitPaths.y,z:orbitPaths.z,hoverinfo:'skip',line:{color:'rgba(110,145,178,0.16)',width:1}},
      starTrace(0),
      {type:'scatter3d',mode:'markers',name:'Galactic core',x:[0],y:[0],z:[0],hoverinfo:'skip',marker:{size:5,color:'#f0d69b',opacity:0.82}}
    ];
    const orbitLayout={paper_bgcolor:BG,plot_bgcolor:BG,margin:{l:0,r:0,t:8,b:0},showlegend:false,scene:{uirevision:'ngc3198-camera',dragmode:'orbit',bgcolor:BG,aspectmode:'cube',xaxis:{title:'x (kpc)',range:[-65,65],gridcolor:GRID,zerolinecolor:GRID,color:textColor},yaxis:{title:'y (kpc)',range:[-65,65],gridcolor:GRID,zerolinecolor:GRID,color:textColor},zaxis:{title:'z (kpc)',range:[-65,65],gridcolor:GRID,zerolinecolor:GRID,color:textColor},camera:{eye:{x:1.35,y:-1.45,z:1.05}}}};
    const orbitPlotReady=Plotly.newPlot(orbitPlot,orbitData,orbitLayout,{responsive:true,displaylogo:false,displayModeBar:true,scrollZoom:true});

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
                const paths=orbitPathCoordinates(),points=starCoordinates(currentFrame);
                orbitPlotReady.then(()=>Promise.all([
                    Plotly.restyle(orbitPlot,{x:[paths.x],y:[paths.y],z:[paths.z]},[0]),
                    Plotly.restyle(orbitPlot,{x:[points.x],y:[points.y],z:[points.z]},[1])
                ]));
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
        let userCamera=null,cameraRevision=0,interactionUntil=0,pointerInteracting=false,touchInteracting=false;
        orbitPlot.addEventListener('wheel',()=>{
            interactionUntil=performance.now()+350;
        },{capture:true,passive:true});
        orbitPlot.addEventListener('pointerdown',()=>{pointerInteracting=true;},true);
        orbitPlot.addEventListener('touchstart',()=>{touchInteracting=true;},{capture:true,passive:true});
        orbitPlot.addEventListener('touchmove',()=>{interactionUntil=performance.now()+350;},{capture:true,passive:true});
        window.addEventListener('pointerup',()=>{
            pointerInteracting=false;
            interactionUntil=performance.now()+250;
        },true);
        window.addEventListener('pointercancel',()=>{
            pointerInteracting=false;
            interactionUntil=performance.now()+250;
        },true);
        window.addEventListener('touchend',event=>{
            touchInteracting=event.touches.length>0;
            if(!touchInteracting)interactionUntil=performance.now()+350;
        },{capture:true,passive:true});
        window.addEventListener('touchcancel',event=>{
            touchInteracting=event.touches.length>0;
            if(!touchInteracting)interactionUntil=performance.now()+350;
        },{capture:true,passive:true});
        orbitPlot.on('plotly_relayout',event=>{
            if(Object.keys(event).some(key=>key.startsWith('scene.camera'))){
                userCamera=JSON.parse(JSON.stringify(orbitPlot.layout.scene.camera));
                cameraRevision++;
            }
        });
        function setFrame(index){
      index=Math.max(0,Math.min(frameCount-1,Number(index)));
            currentFrame=index;
      timeSlider.value=index;timeSliderBottom.value=index;
    const points=starCoordinates(index);
      document.getElementById('frameLabel').textContent=`Frame ${index+1}`;
      document.getElementById('timeValue').textContent=`${(index*DATA.frame_step_myr).toFixed(0)} Myr`;
      document.getElementById('timeLabel').textContent=`${(index*DATA.frame_step_myr).toFixed(0)} Myr`;
            const camera=JSON.parse(JSON.stringify(userCamera||orbitPlot.layout.scene.camera));
            const revision=cameraRevision;
            return Plotly.restyle(orbitPlot,{x:[points.x],y:[points.y],z:[points.z]},[1]).then(()=>{
                if(revision===cameraRevision&&camera&&JSON.stringify(orbitPlot.layout.scene.camera)!==JSON.stringify(camera)){
                    return Plotly.relayout(orbitPlot,{'scene.camera':camera});
                }
            });
    }
    timeSlider.addEventListener('input',()=>setFrame(timeSlider.value));
    timeSliderBottom.addEventListener('input',()=>setFrame(timeSliderBottom.value));
        async function advanceFrame(){
            if(!playing)return;
            if(pointerInteracting||touchInteracting||performance.now()<interactionUntil){
                timer=setTimeout(advanceFrame,50);
                return;
            }
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