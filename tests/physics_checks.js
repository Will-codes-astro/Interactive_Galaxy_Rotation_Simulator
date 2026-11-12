/*
 * Physics checks for the live simulator page (index.html).
 *
 * Runs the page's real JavaScript in Node with the graphics (Three.js, Plotly,
 * the DOM) replaced by do-nothing stand-ins, then calls the physics functions
 * directly. Prints one JSON object of measurements for run_physics_checks.py.
 *
 * Usage: node tests/physics_checks.js index.html
 */
'use strict';
const fs = require('fs');
const vm = require('vm');

const html = fs.readFileSync(process.argv[2] || 'index.html', 'utf8');
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
if (scripts.length === 0) throw new Error('No inline <script> found in index.html');

/* ---------- Stand-ins for the browser ---------- */
const stores = new WeakMap();
const store = t => { if (!stores.has(t)) stores.set(t, {}); return stores.get(t); };
// A "stub" accepts any property access, call or `new` and returns another stub.
function stub(name = 'stub') {
  const target = function () {};
  return new Proxy(target, {
    get(t, k) {
      if (k === Symbol.toPrimitive) return hint => (hint === 'string' ? '' : 0);
      if (k === 'then') return undefined;
      if (k === Symbol.iterator) return function* () { yield stub(name + '[0]'); };
      if (k === 'length') return 0;
      const s = store(t);
      if (!(k in s)) s[k] = stub(name + '.' + String(k));
      return s[k];
    },
    set(t, k, v) { store(t)[k] = v; return true; },
    apply() { return stub(name + '()'); },
    construct() { return stub('new ' + name); },
    has() { return true; },
  });
}

// Elements get their real starting attributes (slider values etc.) from the HTML.
const attrsById = {};
for (const m of html.matchAll(/<([a-z]+)\b([^>]*)>/g)) {
  const attrs = m[2];
  const id = (attrs.match(/\bid="([^"]+)"/) || [])[1];
  if (!id) continue;
  const o = { hidden: false, checked: /\bchecked\b/.test(attrs) };
  for (const key of ['value', 'min', 'max', 'step']) {
    const v = (attrs.match(new RegExp('\\b' + key + '="([^"]*)"')) || [])[1];
    if (v !== undefined) o[key] = v;
  }
  if (o.value !== undefined) o.defaultValue = o.value;
  attrsById[id] = o;
}
const elements = {};
const element = id => elements[id] || (elements[id] = Object.assign(stub('#' + id), attrsById[id] || {}));
const modeRadios = [{ value: 'preset', checked: true }, { value: 'custom', checked: false }]
  .map(o => Object.assign(stub('radio'), o));
const documentStub = stub('document');
documentStub.getElementById = id => element(id);
documentStub.querySelectorAll = sel => (/galaxyMode/.test(sel) ? modeRadios : []);
documentStub.querySelector = sel => { const m = sel.match(/^#([\w-]+)$/); return m ? element(m[1]) : stub('qs'); };

const ctx = {
  console, Math, JSON, Date, Array, Object, Number, String, Boolean, Map, Set, WeakMap, Promise, Symbol,
  Error, Float32Array, Float64Array, Uint8Array, Uint16Array, Uint32Array, Int32Array,
  isFinite, isNaN, parseFloat, parseInt, Infinity, NaN,
  setTimeout: () => 0, clearTimeout: () => 0, setInterval: () => 0, clearInterval: () => 0,
  requestAnimationFrame: () => 0, cancelAnimationFrame: () => 0, performance: { now: () => 0 },
};
for (const n of ['THREE', 'Plotly', 'navigator', 'localStorage', 'sessionStorage', 'matchMedia',
  'getComputedStyle', 'ResizeObserver', 'IntersectionObserver', 'MutationObserver', 'location',
  'history', 'Blob', 'URL', 'screen', 'devicePixelRatio', 'innerWidth', 'innerHeight',
  'addEventListener', 'removeEventListener', 'CustomEvent', 'Event', 'KeyboardEvent', 'PointerEvent',
  'Image', 'Option', 'Audio', 'HTMLElement', 'Element', 'Node', 'CSS', 'AbortController',
  'FileReader', 'TextEncoder', 'fetch', 'alert', 'confirm', 'prompt', 'queueMicrotask', 'structuredClone']) {
  ctx[n] = stub(n);
}
ctx.document = documentStub;
ctx.globalThis = ctx; ctx.window = ctx; ctx.self = ctx;
vm.createContext(ctx);

/* ---------- Run the page and expose its physics ---------- */
const exportApi = `
;globalThis.__PAGE = {
  integrate: (m, rs) => integrateTrajectories(m, rs),
  escapedCount: () => sortedEscapeTimes.length,
  activeHaloMassT: () => activeHaloMassT,
  setMond: on => { mondMode = on; },
  mondVelocity, nfwVelocityModel, baryonicVelocity,
  modelCurves, starCountFor,
  defaults: () => DEFAULTS, frameCount, MOND_A0,
  particles: () => particles,
  galaxy: () => activeGalaxy,
  obs: () => ({ v: obsV, err: obsErr, radii, vBary }),
  upsilon: () => ({ disk: DATA.upsilon_disk, bulge: DATA.upsilon_bulge }),
  selectGalaxy,
};`;
const main = scripts[scripts.length - 1];
for (const s of scripts.slice(0, -1)) vm.runInContext(s, ctx, { filename: 'page-script.js' });
vm.runInContext(main + exportApi, ctx, { filename: 'page-main.js' });
const P = ctx.__PAGE;

/* ---------- Helpers ---------- */
const last = P.frameCount - 1;
function finalPositions(trajectories) {
  // Copy now: the page reuses its trajectory buffer on every integration.
  return Array.from(trajectories, t => [t[3 * last], t[3 * last + 1], t[3 * last + 2]]);
}
function run(massT, rs) {
  const traj = P.integrate(massT, rs);
  return { final: finalPositions(traj), escaped: P.escapedCount(), haloMassT: P.activeHaloMassT() };
}
const meanShift = (a, b) => a.reduce((s, p, i) => s + Math.hypot(p[0] - b[i][0], p[1] - b[i][1], p[2] - b[i][2]), 0) / a.length;
const median = xs => { const s = [...xs].sort((x, y) => x - y); return s[Math.floor(s.length / 2)]; };
function reducedChi(total) {
  const { v, err } = P.obs();
  const chi = v.reduce((s, x, i) => s + ((x - total[i]) / err[i]) ** 2, 0);
  return chi / Math.max(1, v.length - 2);
}

/* ---------- Measurements ---------- */
const out = {};
const d = P.defaults();
const starCount = P.particles().length;
out.galaxy = P.galaxy().id;
out.star_count = starCount;
out.default_halo = { massT: d.massT, scaleRadius: d.scaleRadius };

const withHalo = run(d.massT, d.scaleRadius);
const noHalo = run(0, d.scaleRadius);
const heavyHalo = run(d.massT * 3, d.scaleRadius);
const otherRs = run(d.massT, d.scaleRadius === 40 ? 20 : 40);
out.escaped_with_halo = withHalo.escaped;
out.escaped_visible_only = noHalo.escaped;
out.shift_kpc_halo_mass_0_vs_default = meanShift(withHalo.final, noHalo.final);
out.shift_kpc_halo_mass_x3 = meanShift(withHalo.final, heavyHalo.final);
out.shift_kpc_scale_radius = meanShift(withHalo.final, otherRs.final);

// Shape of the star population at t = 0.
const ps = P.particles();
out.median_abs_z_kpc = median(ps.map(p => Math.abs(p.position[2])));
out.median_R_kpc = median(ps.map(p => Math.hypot(p.position[0], p.position[1])));

// Visible matter on the page vs. SPARC with the gas sign kept.
const g = P.galaxy().rotation, ups = P.upsilon();
let worstGasMismatch = 0, negativeGasPoints = 0;
g.radius_kpc.forEach((r, i) => {
  const gas = g.v_gas_kms[i];
  if (gas < 0) negativeGasPoints++;
  const expected = Math.sqrt(Math.max(0, gas * Math.abs(gas) + ups.disk * g.v_disk_kms[i] ** 2 + ups.bulge * g.v_bulge_kms[i] ** 2));
  if (expected > 1) worstGasMismatch = Math.max(worstGasMismatch, Math.abs(P.baryonicVelocity(r) - expected) / expected);
});
out.negative_gas_points = negativeGasPoints;
out.baryonic_worst_relative_mismatch = worstGasMismatch;

// Goodness of fit.
out.reduced_chi2_visible_only = reducedChi(P.modelCurves(0, d.scaleRadius).total);
out.reduced_chi2_nfw_default = reducedChi(P.modelCurves(d.massT, d.scaleRadius).total);

// MOND.
out.mond_a0 = P.MOND_A0;
P.setMond(true);
const mondA = run(d.massT, d.scaleRadius), mondB = run(d.massT * 3, d.scaleRadius + 10);
out.mond_halo_mass_used = mondA.haloMassT;
out.mond_shift_when_halo_sliders_move = meanShift(mondA.final, mondB.final);
out.reduced_chi2_mond = reducedChi(P.modelCurves(d.massT, d.scaleRadius).total);
P.setMond(false);
// Deep-MOND limit: g -> sqrt(gN * a0), so v^4 -> vN^2 * a0 * r.
{
  const r = 200, vN = 5; // a far, weak-gravity point (kpc, km/s)
  const rM = r * 3.0857e19, vNms = vN * 1000;
  const expected = Math.pow(vNms * vNms * P.MOND_A0 * rM, 0.25) / 1000;
  out.mond_deep_limit_relative_error = Math.abs(P.mondVelocity(vN, r) - expected) / expected;
  out.mond_strong_gravity_relative_boost = P.mondVelocity(300, 0.5) / 300 - 1;
}

// Star count guard.
out.star_count_for_huge_galaxy = P.starCountFor(1e15, 1e10);

// Other presets (best-effort: switching galaxy goes through page UI code).
out.presets = {};
for (const id of ['NGC2403', 'NGC7331']) {
  try {
    P.selectGalaxy(id);
    const dd = P.defaults();
    out.presets[id] = {
      reduced_chi2_visible_only: reducedChi(P.modelCurves(0, dd.scaleRadius).total),
      reduced_chi2_nfw_default: reducedChi(P.modelCurves(dd.massT, dd.scaleRadius).total),
    };
  } catch (e) {
    out.presets[id] = { error: String(e) };
  }
}

process.stdout.write(JSON.stringify(out, null, 2) + '\n');
