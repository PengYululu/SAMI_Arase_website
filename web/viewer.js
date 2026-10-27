// SAMI3 + Arase 3D viewer.
// Reads data/static.json (grids, coastlines, settings), data/arase.json (1-min Arase L, MLT, MLAT, Ne),
// data/manifest.json (list of frames) and one data/frames/*.bin per time step, loaded when needed.
// data/solarwind.json: IMF By, Bz and OMNI E field (1 min) for the time series below.
// data/sami_track.json: SAMI Ne at the Arase position on the same 1-min bins, for the Ne-vs-L panel of the current leg.
// Plotly.react with a fixed uirevision keeps the camera and legend choices while the time changes.
'use strict';

const DEG = Math.PI / 180;
const CACHE_MAX = 40;          // frames kept in memory
const PREFETCH = 6;            // frames loaded ahead while playing

const S = { st: null, arase: null, sami: null, leg: -1, man: null, i: 0, target: 0, playing: false, fps: 4, cache: new Map(),
            pending: null, busy: false, t0ms: 0 };
const $ = id => document.getElementById(id);

// ---------- loading ----------
async function getJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(url + ': HTTP ' + r.status);
  return r.json();
}

function decode(buf, f) {
  const plane = new Float32Array(buf, f.plane.off, f.plane.nr * f.plane.nc);
  const tec = new Float32Array(buf, f.tecarr.off, f.tecarr.nr * f.tecarr.nc);
  const meshes = f.meshes.map(m => {
    const v = new Int16Array(buf, m.voff, m.nv * 3);
    const F = m.fbytes === 2 ? new Uint16Array(buf, m.foff, m.nf * 3) : new Uint32Array(buf, m.foff, m.nf * 3);
    const x = new Float32Array(m.nv), y = new Float32Array(m.nv), z = new Float32Array(m.nv);
    for (let q = 0; q < m.nv; q++) { x[q] = v[3 * q] * f.scale; y[q] = v[3 * q + 1] * f.scale; z[q] = v[3 * q + 2] * f.scale; }
    const I = new Int32Array(m.nf), J = new Int32Array(m.nf), K = new Int32Array(m.nf);
    for (let q = 0; q < m.nf; q++) { I[q] = F[3 * q]; J[q] = F[3 * q + 1]; K[q] = F[3 * q + 2]; }
    return { x, y, z, i: I, j: J, k: K };
  });
  return { plane, tec, meshes };
}

function getFrame(i) {
  if (S.cache.has(i)) {                                   // move to the end = most recently used
    const p = S.cache.get(i); S.cache.delete(i); S.cache.set(i, p); return p;
  }
  const f = S.man.frames[i];
  const p = fetch('data/frames/' + f.file)
    .then(r => { if (!r.ok) throw new Error(f.file + ': HTTP ' + r.status); return r.arrayBuffer(); })
    .then(buf => decode(buf, f));
  p.catch(() => S.cache.delete(i));
  S.cache.set(i, p);
  while (S.cache.size > CACHE_MAX) S.cache.delete(S.cache.keys().next().value);
  return p;
}

// ---------- time helpers ----------
const hToDate = h => new Date(S.t0ms + h * 3600e3);
const fmt = d => d.toISOString().slice(0, 16).replace('T', ' ');
const fmtPlotly = d => d.toISOString().slice(0, 19).replace('T', ' ');

// Arase 1-min bin for real time h (hours since t0), as np.round(sami_h * 60) in the notebook; -1 = outside the data
function araseIndex(h) {
  const k = Math.round(h * 60);
  return (k >= 0 && k < S.arase.n) ? k : -1;
}

// dipole field line r = L cos^2(mlat); azimuth from MLT with noon at +X
function araseXYZ(k) {
  const L = S.arase.L[k], mlt = S.arase.mlt[k], mlat = S.arase.mlat[k];
  if (L === null || mlt === null || mlat === null) return [null, null, null];
  const lat = mlat * DEG, phi = (mlt - 12) * 15 * DEG, r = L * Math.cos(lat) ** 2;
  return [r * Math.cos(lat) * Math.cos(phi), r * Math.cos(lat) * Math.sin(phi), r * Math.sin(lat)];
}

// ---------- 3D figure ----------
// A grid [nr rows, nc cols] as a triangle mesh coloured by val (flat, row-major).
// Plotly mesh3d draws ~10x faster than surface traces, which is what limits the frame rate.
// Grid cells with a NaN corner are left out, like the gaps of a surface.
function gridMesh(nr, nc, xyz, val) {
  const n = nr * nc, x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n), c = new Float32Array(n);
  for (let a = 0; a < nr; a++) for (let b = 0; b < nc; b++) {
    const q = a * nc + b, p = xyz(a, b);
    x[q] = p[0]; y[q] = p[1]; z[q] = p[2]; c[q] = val[q];
  }
  const I = [], J = [], K = [];
  for (let a = 0; a < nr - 1; a++) for (let b = 0; b < nc - 1; b++) {
    const q = a * nc + b;
    if (!(isFinite(c[q]) && isFinite(c[q + 1]) && isFinite(c[q + nc]) && isFinite(c[q + nc + 1]))) continue;
    I.push(q, q + 1); J.push(q + 1, q + nc + 1); K.push(q + nc, q + nc);
  }
  for (let q = 0; q < n; q++) if (!isFinite(c[q])) c[q] = 0;    // unused vertices, keeps the colour range valid
  return { x, y, z, i: Int32Array.from(I), j: Int32Array.from(J), k: Int32Array.from(K), intensity: c };
}

function traces(i, fr) {
  const st = S.st, set = st.settings, f = S.man.frames[i];
  const phi = st.lon.map(lon => (lon - st.mlon_offset + 15 * f.ut - 180) * DEG);   // same as the notebook
  const cosp = phi.map(Math.cos), sinp = phi.map(Math.sin);
  const T = [];

  // isosurfaces
  fr.meshes.forEach((m, q) => T.push({
    type: 'mesh3d', uid: 'iso' + q, name: set.levels[q] + ' /cc', showlegend: true,
    x: m.x, y: m.y, z: m.z, i: m.i, j: m.j, k: m.k, color: set.colors[q], opacity: set.alphas[q],
    flatshading: false, lighting: { ambient: 0.5, diffuse: 0.8, specular: 0.3, roughness: 0.5 }, hoverinfo: 'name' }));

  // Earth: TEC map (or white dayside / black nightside when the TEC file is missing)
  const latE = st.lat_earth.map(v => v * DEG), nl = st.lon.length;
  const sphere = (a, b) => [Math.cos(latE[a]) * cosp[b], Math.cos(latE[a]) * sinp[b], Math.sin(latE[a])];
  // lighting of a Plotly surface trace (its defaults), so the plane and Earth look as in the notebook;
  // with flat lighting the same opacity looks more solid
  const flat = { flatshading: false, lighting: { ambient: 0.8, diffuse: 0.8, specular: 0.05, roughness: 0.5, fresnel: 0.2 } };
  if (f.tec) {
    T.push({ type: 'mesh3d', uid: 'earth', name: 'Earth (TEC)', ...gridMesh(latE.length, nl, sphere, fr.tec), ...flat,
             colorscale: 'Jet', cmin: set.tec_range[0], cmax: set.tec_range[1],
             colorbar: { title: { text: 'TEC (TECU)' }, orientation: 'h', x: 0.5, y: 0.02, yanchor: 'bottom', len: 0.5, thickness: 15 },
             hovertemplate: 'TEC = %{intensity:.1f} TECU<extra>Earth</extra>' });
  } else {
    const day = new Float32Array(latE.length * nl).map((_, q) => (cosp[q % nl] > 0 ? 1 : 0));
    T.push({ type: 'mesh3d', uid: 'earth', name: 'Earth', ...gridMesh(latE.length, nl, sphere, day), ...flat,
             colorscale: [[0, 'black'], [1, 'white']], cmin: 0, cmax: 1, showscale: false, hoverinfo: 'skip' });
  }

  // coastlines (SAMI mlat / mlon), just above the surface
  const rc = 1.005, cx = [], cy = [], cz = [];
  st.coast_mlat.forEach((la, q) => {
    const lo = st.coast_mlon[q];
    if (la === null || lo === null) { cx.push(null); cy.push(null); cz.push(null); return; }
    const p = (lo - st.mlon_offset + 15 * f.ut - 180) * DEG, l = la * DEG;
    cx.push(rc * Math.cos(l) * Math.cos(p)); cy.push(rc * Math.cos(l) * Math.sin(p)); cz.push(rc * Math.sin(l));
  });
  T.push({ type: 'scatter3d', uid: 'coast', name: 'coastlines', mode: 'lines', x: cx, y: cy, z: cz,
           line: { color: 'white', width: 2 }, connectgaps: false, hoverinfo: 'skip', visible: f.tec ? true : 'legendonly' });

  // equatorial plane ne[45, :, :] at z = 0
  const r = st.r_plane;
  const ticks = []; for (let t = set.eq_range[0]; t <= set.eq_range[1]; t++) ticks.push(t);
  T.push({ type: 'mesh3d', uid: 'plane', name: 'equatorial plane', showlegend: true,
           ...gridMesh(r.length, nl, (a, b) => [r[a] * cosp[b], r[a] * sinp[b], 0], fr.plane), ...flat,
           colorscale: 'Jet', cmin: set.eq_range[0], cmax: set.eq_range[1], opacity: set.eq_alpha,
           colorbar: { title: { text: 'Ne (/cm<sup>3</sup>)' }, tickvals: ticks, ticktext: ticks.map(t => '10<sup>' + t + '</sup>'), len: 0.6 },
           hovertemplate: 'log10 Ne = %{intensity:.2f} (/cc)<extra>equator</extra>' });

  // Arase orbit +/- track_hr and the Arase dot coloured by L3 Ne on the plane's colour scale
  const k = araseIndex(f.h);
  const tx = [], ty = [], tz = [];
  if (k >= 0 && set.track_hr > 0) {
    const n = Math.round(set.track_hr * 60);
    for (let q = Math.max(0, k - n); q <= Math.min(S.arase.n - 1, k + n); q++) {
      const [x, y, z] = araseXYZ(q); tx.push(x); ty.push(y); tz.push(z);
    }
  }
  T.push({ type: 'scatter3d', uid: 'track', name: 'Arase orbit +/- ' + set.track_hr + ' h', mode: 'lines',
           x: tx, y: ty, z: tz, line: { color: 'gray', width: 3 }, hoverinfo: 'skip', connectgaps: false });

  const [xa, ya, za] = k >= 0 ? araseXYZ(k) : [null, null, null];
  const ne = k >= 0 ? S.arase.ne[k] : null;
  const hasDot = xa !== null;
  T.push({ type: 'scatter3d', uid: 'arase', name: 'Arase', mode: 'markers',
           x: hasDot ? [xa] : [], y: hasDot ? [ya] : [], z: hasDot ? [za] : [],
           marker: { size: set.dot_size, color: ne !== null && ne > 0 ? [Math.log10(ne)] : ['gray'],
                     colorscale: 'Jet', cmin: set.eq_range[0], cmax: set.eq_range[1], showscale: false,
                     line: { color: 'black', width: 2 } },
           hovertemplate: hasDot ? ('Arase<br>L = ' + S.arase.L[k].toFixed(2) + '<br>MLT = ' + S.arase.mlt[k].toFixed(1) +
                          ' h<br>MLAT = ' + S.arase.mlat[k].toFixed(1) + ' deg<br>Ne = ' + (ne === null ? 'n/a' : ne.toPrecision(3)) +
                          ' /cc<extra></extra>') : '' });
  return T;
}

function layout3d(i) {
  const set = S.st.settings, f = S.man.frames[i];
  const ax = t => ({ title: { text: t }, range: [-set.lim, set.lim], dtick: 1, showbackground: false, gridcolor: 'lightgray', zeroline: false });
  return {
    uirevision: 'keep',                                   // keeps camera, zoom and legend clicks across time steps
    title: { text: 'SAMI3 Ne + Arase  ' + fmt(hToDate(f.h)).replace(/-/g, '/') + ' UT  (SAMI ' + f.sami.slice(11) + ', shift ' + S.st.shift + ')',
             font: { size: 15 } },
    margin: { l: 0, r: 0, t: 40, b: 0 },
    scene: { xaxis: ax('X (Re, Sun)'), yaxis: ax('Y (Re, dusk)'), zaxis: ax('Z (Re)'), aspectmode: 'cube',
             camera: set.camera, uirevision: 'keep' },
    legend: { x: 0.02, y: 0.95 },
  };
}

// ---------- time series below: IMF By/Bz, OMNI E field, Arase Ne (one figure, shared time axis) ----------
function drawTimeline() {
  const a = S.arase, w = S.sw, x = [], xs = [], ys = [];
  for (let q = 0; q < a.n; q++) x.push(fmtPlotly(hToDate((q + 0.5) / 60)));           // Arase bins: centre of each minute
  for (let q = 0; q < w.n; q++) xs.push(fmtPlotly(hToDate(q / 60)));                  // OMNI / IMF: whole minutes
  const fr = S.man.frames;
  const line = (xx, y, name, col, ax, unit) => ({ x: xx, y, type: 'scatter', mode: 'lines', name, xaxis: 'x', yaxis: ax,
    line: { color: col, width: 1.2 }, hovertemplate: '%{x}<br>' + name + ' = %{y:.3g} ' + unit + '<extra></extra>' });
  const zero = ax => ({ type: 'line', xref: 'paper', yref: ax, x0: 0, x1: 1, y0: 0, y1: 0, line: { color: 'black', width: 1, dash: 'dot' } });
  S.tlShapes = [zero('y'), zero('y2')];
  Plotly.newPlot('plot2d', [
    line(xs, w.by, 'By', '#2ca02c', 'y', 'nT'),
    line(xs, w.bz, 'Bz', '#1f4fbf', 'y', 'nT'),
    line(xs, w.ef, 'E field', '#9467bd', 'y2', 'mV/m'),
    line(x, a.ne, 'Arase L3 Ne', '#555', 'y3', '/cc'),
  ], {
    margin: { l: 70, r: 20, t: 10, b: 35 }, hovermode: 'x unified',
    yaxis:  { domain: [0.70, 1.00], title: { text: 'IMF (nT, GSM)' } },
    yaxis2: { domain: [0.36, 0.64], title: { text: 'E field (mV/m)' } },
    yaxis3: { domain: [0.00, 0.30], title: { text: 'Arase Ne (/cc)' }, type: 'log', exponentformat: 'power' },
    xaxis:  { range: [fmtPlotly(hToDate(fr[0].h)), fmtPlotly(hToDate(fr[fr.length - 1].h))], anchor: 'y3' },
    legend: { orientation: 'h', x: 1, xanchor: 'right', y: 1.0, yanchor: 'bottom' },
    shapes: S.tlShapes.concat([cursorShape(0)]),
  }, { responsive: true, displaylogo: false });
  $('plot2d').on('plotly_click', ev => {                  // click any panel -> jump to the nearest frame
    const h = (Date.parse(ev.points[0].x.replace(' ', 'T') + 'Z') - S.t0ms) / 3600e3;
    let best = 0;
    S.man.frames.forEach((f, q) => { if (Math.abs(f.h - h) < Math.abs(S.man.frames[best].h - h)) best = q; });
    pause(); request(best);
  });
}
function cursorShape(i) {
  const x = fmtPlotly(hToDate(S.man.frames[i].h));
  return { type: 'line', xref: 'x', yref: 'paper', x0: x, x1: x, y0: 0, y1: 1, line: { color: '#d62728', width: 2 } };
}

// ---------- Ne vs L along the current Arase leg: SAMI3 (solid) and Arase L3 (dotted) ----------
const C_SAMI = '#1f4fbf', C_ARASE = '#d62728';
function legOf(k) { return k < 0 ? -1 : S.arase.legs.findIndex(l => k >= l.i0 && k <= l.i1); }

function logRange(v) {
  const g = v.filter(x => x !== null && x > 0);
  if (!g.length) return [0, 5];
  return [Math.floor(Math.log10(Math.min(...g)) * 2) / 2, Math.ceil(Math.log10(Math.max(...g)) * 2) / 2];
}

function legDots(k) {
  const a = S.arase, s = S.sami.ne, L = k >= 0 ? a.L[k] : null;
  const pt = v => (L !== null && v !== null && v > 0) ? { x: [L], y: [v] } : { x: [], y: [] };
  return [pt(s[k]), pt(a.ne[k])];
}

function drawLeg(k) {
  const g = legOf(k), a = S.arase, s = S.sami.ne;
  const dot = (col, name) => ({ type: 'scatter', mode: 'markers', name, showlegend: false, x: [], y: [],
    marker: { size: 13, color: col, line: { color: 'black', width: 1.5 } },
    hovertemplate: name + '<br>L = %{x:.2f}<br>Ne = %{y:.3g} /cc<extra></extra>' });
  if (g !== S.leg) {                                       // new leg: whole SAMI and Arase lines of that leg
    S.leg = g;
    const T = [], lay = { margin: { l: 70, r: 15, t: 60, b: 50 }, uirevision: 'leg' + g, hovermode: 'closest',
      xaxis: { title: { text: 'L (90°)' } }, yaxis: { title: { text: 'Ne (/cm<sup>3</sup>)' }, type: 'log', exponentformat: 'power' },
      legend: { x: 0.99, xanchor: 'right', y: 0.99, bgcolor: 'rgba(255,255,255,0.8)' } };
    if (g < 0) {
      lay.title = { text: 'No Arase leg at this time' };
    } else {
      const l = a.legs[g], ii = [];
      for (let q = l.i0; q <= l.i1; q++) ii.push(q);
      const x = ii.map(q => a.L[q]), tt = ii.map(q => fmt(hToDate((q + 0.5) / 60)).slice(11));
      const hov = name => name + '  %{text} UT<br>L = %{x:.2f}<br>Ne = %{y:.3g} /cc<extra></extra>';
      const ys = ii.map(q => s[q]), ya = ii.map(q => a.ne[q]);
      T.push({ type: 'scatter', mode: 'lines', name: 'SAMI3', x, y: ys, text: tt, line: { color: C_SAMI, width: 2.5 }, hovertemplate: hov('SAMI3') });
      T.push({ type: 'scatter', mode: 'lines', name: 'Arase L3', x, y: ya, text: tt, line: { color: C_ARASE, width: 2.5, dash: 'dot' }, hovertemplate: hov('Arase L3') });
      const L = x.filter(v => v !== null);
      lay.xaxis.range = [Math.floor(Math.min(...L) * 2) / 2, Math.ceil(Math.max(...L) * 2) / 2];
      lay.yaxis.range = logRange(ys.concat(ya));
      const sd = l.start.slice(4, 8), ed = l.end.slice(4, 8);
      lay.title = { text: (l.dir === 'outbound' ? 'Outbound' : 'Inbound') + ' leg  ' + sd + '/' + l.start.slice(9) + ' – ' +
                          (ed === sd ? '' : ed + '/') + l.end.slice(9) + ' UT' +
                          (s.slice(l.i0, l.i1 + 1).some(v => v !== null) ? '' : '  (no SAMI data)') +
                          '<br><sub>SAMI3 (solid) vs Arase L3 (dotted); dots = current time</sub>' };
    }
    T.push(dot(C_SAMI, 'SAMI3 now'), dot(C_ARASE, 'Arase L3 now'));
    Plotly.react('plotleg', T, lay, { responsive: true, displaylogo: false });
  }
  if (g < 0) return;
  const [ds, da] = legDots(k), n = document.getElementById('plotleg').data.length;
  Plotly.restyle('plotleg', { x: [ds.x, da.x], y: [ds.y, da.y] }, [n - 2, n - 1]);
}

// ---------- drawing / navigation ----------
function labels(i) {
  const f = S.man.frames[i], k = araseIndex(f.h), a = S.arase;
  $('slider').value = i;
  $('time').innerHTML = '<b>' + fmt(hToDate(f.h)) + ' UT</b> &nbsp; SAMI file ' + f.sami + ' &nbsp; (' + (i + 1) + '/' + S.man.frames.length + ')';
  $('arase').textContent = k < 0 ? 'Arase: no data at this time' :
    'Arase: L = ' + (a.L[k] === null ? 'n/a' : a.L[k].toFixed(2)) + ', MLT = ' + (a.mlt[k] === null ? 'n/a' : a.mlt[k].toFixed(1)) +
    ' h, MLAT = ' + (a.mlat[k] === null ? 'n/a' : a.mlat[k].toFixed(1)) + ' deg, L3 Ne = ' + (a.ne[k] === null ? 'n/a' : a.ne[k].toPrecision(3)) +
    ' /cc' + (a.qf[k] > 0 ? ' (quality flag ' + a.qf[k] + ')' : '');
}

async function draw(i, fr) {
  await Plotly.react('plot3d', traces(i, fr), layout3d(i), { responsive: true, displaylogo: false });
  Plotly.relayout('plot2d', { shapes: S.tlShapes.concat([cursorShape(i)]) });
  drawLeg(araseIndex(S.man.frames[i].h));
  S.i = i;
}

// show frame i; while frames load, only the newest request is drawn
function request(i) {
  S.pending = i; S.target = i; labels(i);
  if (!S.busy) pump();
}
async function pump() {
  S.busy = true;
  while (S.pending !== null) {
    const i = S.pending; S.pending = null;
    try {
      $('status').textContent = S.cache.has(i) ? '' : 'loading…';
      const fr = await getFrame(i);
      if (S.pending !== null) continue;                     // a newer step was asked for meanwhile
      await draw(i, fr);
      $('status').textContent = '';
      for (const d of [1, -1]) if (i + d >= 0 && i + d < S.man.frames.length) getFrame(i + d);
    } catch (e) { $('status').textContent = String(e.message || e); }
  }
  S.busy = false;
}

const sleep = ms => new Promise(r => setTimeout(r, ms));
async function play() {
  if (S.playing) return;
  S.playing = true; $('play').textContent = '❚❚ Pause';
  const n = S.man.frames.length;
  while (S.playing) {
    const t = performance.now(), i = (S.i + 1) % n;
    S.target = i;
    for (let d = 1; d <= PREFETCH; d++) getFrame((i + d) % n);
    labels(i);
    try { await draw(i, await getFrame(i)); } catch (e) { $('status').textContent = String(e.message || e); pause(); break; }
    const wait = 1000 / S.fps - (performance.now() - t);
    if (wait > 0) await sleep(wait);
  }
}
function pause() { S.playing = false; $('play').textContent = '▶ Play'; }
function step(d) { pause(); request((S.target + d + S.man.frames.length) % S.man.frames.length); }   // from the last asked-for step, so fast key presses add up

// ---------- start ----------
async function main() {
  if (location.protocol === 'file:') {
    document.body.insertAdjacentHTML('afterbegin', '<div class="err">This page has to be served over http, not opened as a file. ' +
      'In a terminal: <code>cd Arase_SAMI_website/web &amp;&amp; python3 -m http.server 8000</code>, then open http://localhost:8000 (see README.txt).</div>');
    return;
  }
  [S.st, S.arase, S.man, S.sami, S.sw] = await Promise.all([getJSON('data/static.json'), getJSON('data/arase.json'),
    getJSON('data/manifest.json'), getJSON('data/sami_track.json'), getJSON('data/solarwind.json')]);
  S.t0ms = Date.parse(S.st.t0 + 'T00:00:00Z');
  S.fps = +$('fps').value;
  $('slider').max = S.man.frames.length - 1;
  $('subtitle').textContent += '  ' + S.man.frames.length + ' SAMI steps, ' + S.st.dt_sami_min + ' min apart. Keys: space = play/pause, ← → = step.';

  $('slider').addEventListener('input', e => { pause(); request(+e.target.value); });
  // blur after a click, so a later space key is not also a button press
  $('play').addEventListener('click', e => { e.currentTarget.blur(); S.playing ? pause() : play(); });
  $('prev').addEventListener('click', e => { e.currentTarget.blur(); step(-1); });
  $('next').addEventListener('click', e => { e.currentTarget.blur(); step(1); });
  $('fps').addEventListener('change', e => { S.fps = +e.target.value; });
  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT') return;
    if (e.key === ' ') { e.preventDefault(); S.playing ? pause() : play(); }
    else if (e.key === 'ArrowRight') { e.preventDefault(); step(1); }
    else if (e.key === 'ArrowLeft') { e.preventDefault(); step(-1); }
  });

  drawTimeline();
  request(0);
}
main().catch(e => { document.body.insertAdjacentHTML('afterbegin', '<div class="err">' + e + '</div>'); });
