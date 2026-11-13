## ===== Step 3 (optional): plasmasphere contour mapping and SED edge overlays, one small file per SAMI output =====
# Same as the "Plasmasphere Contour mapping to TEC" and "SED edge from the SAMI TEC gradient" sections of
# sami_arase_plasmasphere.ipynb (add_pla_mapping, sami_sed_edge, add_sed_edge, add_edge_pla), copied here.
# The page loads these files only when "add contour mapping" / "add SED edge" is ticked above the 3D view.
#   contour mapping: contours of the equatorial plane ne[nlat//2, :, MAP_K0:] at MAP_LEVELS (/cc), each contour
#                    point followed along the SAMI field line with its apex there (same mlon) to its north / south
#                    ends -> footprints on the TEC Earth (r = 1.01 Re)
#   SED edge       : |grad TEC| > SED_GRAD (TECU/km) in SED_LON x |MLAT| in SED_LAT, middle point of each run along
#                    every mlat row = edge point (on the TEC Earth, north / south); each edge point -> the SAMI field
#                    line whose north / south foot point at H_FOOT km is closest -> that line's apex on the
#                    equatorial plane (diamond), and the field line itself
# Output: web/data/<event>/overlays/contour/yymmdd_HHMM.json, .../overlays/sed/yymmdd_HHMM.json (x, y, z in Re in the
# frame of that SAMI output, +X = noon) and web/data/<event>/overlays.json (settings + which frames have files).
# Needs step 1 and step 2 (manifest.json) of the same event; the SAMI drive must be mounted.
# Run:  python3 step3_build_overlays.py --event 20241010
#       python3 step3_build_overlays.py --event 20241010 --every 20   (quick test)
import os
import sys
import json
import argparse
import datetime
import numpy as np
import contourpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
from SAMI_functions import read_info_file, read_coordinates

ap = argparse.ArgumentParser()
ap.add_argument('--event', default=DEFAULT_EVENT, choices=list(EVENTS), help='event key in config.EVENTS')
ap.add_argument('--every', type=int, default=1, help='only every N-th frame (quick test; the page shows no overlay for the others)')
args = ap.parse_args()
globals().update(load_event(args.event))
print('===== event %s: %s =====' % (EVENT, TITLE))
if not OVERLAYS.get(EVENT):
    raise SystemExit('overlays are not switched on for %s (OVERLAYS in config.py)' % EVENT)

with open(os.path.join(DATA_DIR, 'manifest.json')) as f:
    frames = json.load(f)['frames']
OUT = {k: os.path.join(DATA_DIR, 'overlays', k) for k in ('contour', 'sed')}
for d in OUT.values():
    os.makedirs(d, exist_ok=True)

# ---- SAMI field lines [nlt mlon, nf line, nz point along the line] ----
Grid, Vars = read_info_file(SAMI_GRID_DIR + 'run_information.dat')
nz, nf, nlt = Grid['Points'], Grid['nFL'], Grid['nLons']
bLon, bLat, bAlt = read_coordinates(SAMI_GRID_DIR, 'blonu.dat', 'blatu.dat', 'baltu.dat', nf, nz, nlt)
bAlt = bAlt - Re
alt_grid = np.max(bAlt, axis=2)[0] - 1                       # km, layer k of the ne cube is at alt_grid[k+1]
apex_alt = np.max(bAlt, axis=2)                              # [nlt, nf]
lon_apex = np.mod(bLon[:, 0, nz // 2], 360)                  # mlon of each line's apex
foot_a, foot_b = bLat[:, :, 0], bLat[:, :, -1]               # mlat of both ends


def to_phi(mlon, ut):
    """mlon -> MLT -> azimuth (rad) with noon at +X, same formula as sami_isosurface_fig"""
    return np.radians((np.mod((np.asarray(mlon) - MLON_OFFSET) / 15.0 + ut, 24) - 12) * 15)


def sph2xyz(r, lat_deg, phi):
    lat = np.radians(lat_deg)
    return r * np.cos(lat) * np.cos(phi), r * np.cos(lat) * np.sin(phi), r * np.sin(lat)


def L(a, nd=3):
    """array -> JSON list, rounded, NaN -> null (line breaks)"""
    return [None if not np.isfinite(v) else round(float(v), nd) for v in np.ravel(a)]


def foot_points(h_foot):
    """where every SAMI field line crosses h_foot km on its north and south halves (sami_foot_points of the notebook)
    returns {'north': (mlat, mlon), 'south': (mlat, mlon)}, arrays [nlt, nf], NaN if the apex is below h_foot"""
    fa = np.full((2, 2, nlt, nf), np.nan)
    for il in range(nlt):
        for f in range(nf):
            alt, lat = bAlt[il, f], bLat[il, f]
            lon = np.unwrap(np.radians(bLon[il, f]))
            top = np.argmax(alt)
            if alt[top] < h_foot or top == 0 or top == alt.size - 1:
                continue
            for h, idx in enumerate([np.arange(top + 1), np.arange(alt.size - 1, top - 1, -1)]):
                if alt[idx[0]] > h_foot:
                    continue
                fa[h, 0, il, f] = np.interp(h_foot, alt[idx], lat[idx])
                fa[h, 1, il, f] = np.mod(np.degrees(np.interp(h_foot, alt[idx], lon[idx])), 360)
    h_n = 0 if np.nanmean(fa[0, 0]) > np.nanmean(fa[1, 0]) else 1
    return {'north': (fa[h_n, 0], fa[h_n, 1]), 'south': (fa[1 - h_n, 0], fa[1 - h_n, 1])}


# ---- 1. plasmasphere contour -> footprints (add_pla_mapping) ----
def pla_mapping(ne, ut):
    nlat, nlon, nalt = ne.shape
    lon_c = np.linspace(0, 360, nlon + 1)[:-1] + 180. / nlon
    r_re  = (alt_grid[1:nalt + 1] + Re) / Re
    eq    = ne[nlat // 2, :, MAP_K0:]
    eq    = np.concatenate([eq, eq[:1, :]], axis=0)                       # close across 0/360 deg
    gen   = contourpy.contour_generator(np.r_[lon_c, lon_c[0] + 360], r_re[MAP_K0:], np.ma.masked_invalid(eq.T))
    out = []
    for lv, col in zip(MAP_LEVELS, MAP_COLORS):
        segs = [s for s in gen.lines(lv) if len(s) > 1]
        if not segs:
            continue
        mlon = np.concatenate([np.mod(s[:, 0], 360) for s in segs])
        r    = np.concatenate([s[:, 1] for s in segs])
        alt  = (r - 1) * Re
        ilon = np.argmin(np.abs((mlon[:, None] - lon_apex[None, :] + 180) % 360 - 180), axis=1)
        mlat_n, mlat_s = np.full(mlon.size, np.nan), np.full(mlon.size, np.nan)
        for p in range(mlon.size):
            a = apex_alt[ilon[p]]
            o = np.argsort(a)
            if not (a[o[0]] <= alt[p] <= a[o[-1]]):
                continue
            la_, lb_ = np.interp(alt[p], a[o], foot_a[ilon[p], o]), np.interp(alt[p], a[o], foot_b[ilon[p], o])
            mlat_n[p], mlat_s[p] = max(la_, lb_), min(la_, lb_)
        # contour on the equatorial plane, segments separated by NaN
        xs, ys = [], []
        for s in segs:
            ph = to_phi(s[:, 0], ut)
            xs += [s[:, 1] * np.cos(ph), [np.nan]]
            ys += [s[:, 1] * np.sin(ph), [np.nan]]
        # north + south footprints on the Earth
        ok = np.isfinite(mlat_n)
        fm, fl, fr = np.r_[mlon[ok], mlon[ok]], np.r_[mlat_n[ok], mlat_s[ok]], np.r_[r[ok], r[ok]]
        fx, fy, fz = sph2xyz(R_FOOT, fl, to_phi(fm, ut))
        out.append(dict(level=lv, color=col, x=L(np.concatenate(xs)), y=L(np.concatenate(ys)),
                        fx=L(fx), fy=L(fy), fz=L(fz), fmlat=L(fl, 1), fmlon=L(fm, 1), fr=L(fr, 2)))
    return out


# ---- 2. SED edge points of the TEC map (sami_sed_edge) ----
def sed_edge(tec):
    nlat, nlon = tec.shape
    lat_c = np.linspace(-85, 85, nlat + 1)[:-1] + 85. / nlat
    lon_c = np.linspace(0, 360, nlon + 1)[:-1] + 180. / nlon
    LAT, LON = np.meshgrid(lat_c, lon_c, indexing='ij')
    dlat_km = np.radians(lat_c[1] - lat_c[0]) * 6371.0
    dlon_km = np.radians(lon_c[1] - lon_c[0]) * 6371.0 * np.cos(np.radians(lat_c))[:, None]
    dtec_dlat, dtec_dlon = np.gradient(tec)
    grad = np.sqrt((dtec_dlat / dlat_km)**2 + (dtec_dlon / dlon_km)**2)
    region = ((LON >= SED_LON[0]) & (LON <= SED_LON[1]) & (np.abs(LAT) >= SED_LAT[0]) & (np.abs(LAT) <= SED_LAT[1]))
    with np.errstate(invalid='ignore'):
        mask = region & (grad > SED_GRAD)
    pts = []
    for ilat in range(nlat):
        d = np.diff(np.r_[0, mask[ilat].astype(int), 0])
        for start, end in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0] - 1):
            mid = (start + end) // 2
            pts.append((LON[ilat, mid], LAT[ilat, mid], grad[ilat, mid]))
    return np.array(pts, dtype=float).reshape(-1, 3)


# ---- 3. SED edge on the Earth and mapped to the equatorial plane (add_sed_edge + add_edge_pla) ----
uv = lambda lat, lon: np.stack([np.cos(np.radians(lat)) * np.cos(np.radians(lon)),
                                np.cos(np.radians(lat)) * np.sin(np.radians(lon)), np.sin(np.radians(lat))], -1)


def sed_mapping(edge, ut, foot):
    out = []
    mlon_e, mlat_e, grad_e = edge[:, 0], edge[:, 1], edge[:, 2]
    for hemi, sel, col in [('north', mlat_e > 0, EDGE_COLORS[0]), ('south', mlat_e < 0, EDGE_COLORS[1])]:
        if not sel.any():
            continue
        ex, ey, ez = sph2xyz(R_FOOT, mlat_e[sel], to_phi(mlon_e[sel], ut))
        flat, flon = foot[hemi]
        valid = np.nonzero(np.isfinite(flat.ravel()))[0]
        cosd  = uv(mlat_e[sel], mlon_e[sel]) @ uv(flat.ravel()[valid], flon.ravel()[valid]).T
        k     = valid[np.argmax(cosd, axis=1)]
        dist  = np.arccos(np.clip(cosd.max(axis=1), -1, 1)) * (Re + H_FOOT)
        il, f = np.unravel_index(k, flat.shape)
        top   = np.argmax(bAlt[il, f], axis=1)
        a_mlon, a_r = np.mod(bLon[il, f, top], 360), (bAlt[il, f, top] + Re) / Re
        ph = to_phi(a_mlon, ut)
        # field lines used (each line once), every FL_EVERY-th point plus both ends
        lx, ly, lz = [], [], []
        if EDGE_FL:
            idx = np.unique(np.r_[np.arange(0, nz, FL_EVERY), nz - 1])
            for a, b in sorted(set(zip(il.tolist(), f.tolist()))):
                x, y, z = sph2xyz((bAlt[a, b, idx] + Re) / Re, bLat[a, b, idx], to_phi(bLon[a, b, idx], ut))
                lx += [x, [np.nan]]; ly += [y, [np.nan]]; lz += [z, [np.nan]]
        cat = lambda v: np.concatenate(v) if v else np.zeros(0)
        out.append(dict(hemi=hemi, color=col,
                        ex=L(ex), ey=L(ey), ez=L(ez), emlon=L(mlon_e[sel], 1), emlat=L(mlat_e[sel], 1), egrad=L(grad_e[sel], 4),
                        ax=L(a_r * np.cos(ph)), ay=L(a_r * np.sin(ph)), ar=L(a_r, 2), amlon=L(a_mlon, 1), dist=L(dist, 0),
                        lx=L(cat(lx)), ly=L(cat(ly)), lz=L(cat(lz))))
    return out


print('SAMI field-line foot points at %g km ...' % H_FOOT)
foot = foot_points(H_FOOT)
have_c, have_s = [], []
size = {'contour': 0, 'sed': 0}
for n, fr in enumerate(frames):
    if n % args.every:
        continue
    t   = datetime.datetime.strptime(fr['sami'], '%Y-%m-%d %H:%M')
    ut  = fr['ut']
    base = fr['file'][:-4] + '.json'
    ne_file  = SAMI_NE_DIR + t.strftime('%y%m%d_%H%M') + '.npy'
    tec_file = SAMI_TEC_DIR + t.strftime('%Y%m%d_%H%M') + '_tmp_inMLON.npy'
    nc = ns = 0
    if os.path.isfile(ne_file):
        c = pla_mapping(np.load(ne_file), ut)
        s = json.dumps(dict(levels=c), separators=(',', ':'))
        with open(os.path.join(OUT['contour'], base), 'w') as fo:
            fo.write(s)
        size['contour'] += len(s); have_c.append(n); nc = sum(len(d['fx']) for d in c)
    if os.path.isfile(tec_file):
        e = sed_edge(np.load(tec_file))
        h = sed_mapping(e, ut, foot)
        s = json.dumps(dict(hemis=h), separators=(',', ':'))
        with open(os.path.join(OUT['sed'], base), 'w') as fo:
            fo.write(s)
        size['sed'] += len(s); have_s.append(n); ns = len(e)
    if n % 20 == 0 or n == len(frames) - 1:
        print('%3d/%d  %s  footprints %d, SED edge points %d' % (n + 1, len(frames), base, nc, ns))

with open(os.path.join(DATA_DIR, 'overlays.json'), 'w') as f:
    json.dump(dict(created=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
                   contour=dict(frames=have_c, levels=MAP_LEVELS, colors=MAP_COLORS, k0=MAP_K0),
                   sed=dict(frames=have_s, colors=EDGE_COLORS, grad=SED_GRAD, lon=SED_LON, lat=SED_LAT, h_foot=H_FOOT),
                   settings=dict(map_marker=MAP_MARKER, edge_marker=EDGE_MARKER)), f, separators=(',', ':'))
print('wrote overlays.json: contour %d frames (%.0f MB), SED %d frames (%.0f MB)'
      % (len(have_c), size['contour'] / 1e6, len(have_s), size['sed'] / 1e6))
write_events_json()
