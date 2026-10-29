## ===== Step 2: one binary frame file per SAMI output + manifest.json =====
# For every SAMI output (sami_t from step 1):
#   - isosurfaces of ne[:, :, K0:] at LEVELS (marching_cubes, STEP_SIZE), vertices already in x, y, z (Re)
#     in the SM-like frame of the notebook (+X = noon, +Y = dusk) at that UT
#   - equatorial plane log10(ne[45, :, :])      [nalt, nlon+1]
#   - Earth TEC map (pole rows filled)          [nlat+2, nlon+1]
# written to web/data/frames/yymmdd_HHMM.bin (SAMI file time). The byte layout of each file is in manifest.json.
# Also SAMI Ne at the Arase position on the 1-min grid (notebook SAMI-Arase step 2, sami_ne_mag) -> web/data/sami_track.json
# Run:  python3 step2_build_frames.py                 (all SAMI outputs, rewrites every frame and manifest.json)
#       python3 step2_build_frames.py --every 5       (every 5th SAMI output only, for a quick test)
import os
import sys
import json
import argparse
import datetime
import numpy as np
from skimage.measure import marching_cubes
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

ap = argparse.ArgumentParser()
ap.add_argument('--every', type=int, default=1, help='use every N-th SAMI output')
args = ap.parse_args()

st       = np.load(os.path.join(INPUT_DIR, 'sami_static.npz'))
alt_grid = st['alt_grid']
dt_sami  = float(st['dt_sami'])
sami_t   = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S.%f') for s in st['sami_t']][::args.every]
t0       = datetime.datetime.strptime(T0_STR, '%Y-%m-%d')
os.makedirs(FRAME_DIR, exist_ok=True)

# ---- Arase R, MLAT, MLT (step 1) -> SAMI Ne along the orbit, as sami_ne_mag in the notebook ----
# The SAMI cube ne[la, lo, height_pla] is on SAMI's Earth-fixed magnetic grid. Arase MLAT is la, MLT -> mlon with
# mlon = 15 (MLT - ut) + 72 at each SAMI output's own ut, altitude = (R - 1) Re. Notebook axes:
# la = calc_tec lat bin edges, lo = mlon bin centres (wrapped on both sides), height_pla = alt_grid[1:nalt+1]
trk      = np.load(os.path.join(INPUT_DIR, 'arase_track.npz'))
h_trk    = trk['hours_l3']
alt_trk  = (trk['r_1min'] - 1) * Re
mlat_trk = trk['mlat_1min'];  mlt_trk = trk['mlt_1min']
ok_trk   = np.isfinite(alt_trk) & np.isfinite(mlat_trk) & np.isfinite(mlt_trk)
ne_trk   = np.full(h_trk.size, np.nan)
prev     = None                                                                     # previous output: (h, ut, interpolator)


def phi_of(lon, ut):
    """magnetic longitude (deg) -> azimuth (rad) with noon at +X; same as the notebook's unwrapped MLT angle"""
    return np.radians(lon - MLON_OFFSET + 15. * ut - 180.)


def put(buf, arr):
    """append arr to buf, padded to 4 bytes; return its byte offset"""
    off = len(buf)
    buf += np.ascontiguousarray(arr).tobytes()
    buf += b'\0' * ((-len(buf)) % 4)
    return off


frames = []
for n, t in enumerate(sami_t):
    fname = t.strftime('%y%m%d_%H%M') + '.bin'
    ut    = t.hour + t.minute / 60.                                  # notebook uses hour + minute only
    h_real = (t - t0).total_seconds() / 3600. - SHIFT * dt_sami      # real (Arase) time, hours since T0_STR

    ne = np.load(SAMI_NE_DIR + t.strftime('%y%m%d_%H%M') + '.npy')   # [90 mlat, 180 mlon, 118 alt]
    nlat, nlon, nalt = ne.shape
    lat_c = np.linspace(-85, 85, nlat + 1)[:-1] + 85. / nlat
    lon_c = np.linspace(0, 360, nlon + 1)[:-1] + 180. / nlon
    lon_w = np.r_[lon_c, lon_c[0] + 360]
    ne_w  = np.concatenate([ne, ne[:, :1, :]], axis=1)                # closed across 0/360 deg

    # SAMI Ne along the Arase orbit: log10 Ne at this output, linear in time from the previous output
    la   = np.linspace(-85, 85, nlat + 1)[:-1]
    lo_c = np.linspace(0, 360, nlon + 1)[:-1] + 180. / nlon                        # mlon bin centres
    lg   = np.log10(np.clip(ne, 1e-3, None))
    itp  = RegularGridInterpolator((la, np.r_[lo_c[-1] - 360, lo_c, lo_c[0] + 360], alt_grid[1:nalt + 1]),
                                   np.concatenate([lg[:, -1:, :], lg, lg[:, :1, :]], axis=1),
                                   bounds_error=False, fill_value=np.nan)
    if prev is not None:
        m = (h_trk >= prev[0]) & (h_trk < h_real) & ok_trk
        w = (h_trk[m] - prev[0]) / (h_real - prev[0])
        at = lambda u, f: f(np.c_[mlat_trk[m], np.mod(15 * (mlt_trk[m] - u) + MLON_OFFSET, 360), alt_trk[m]])
        ne_trk[m] = 10.**((1 - w) * at(prev[1], prev[2]) + w * at(ut, itp))
    prev = (h_real, ut, itp)

    buf = bytearray()
    # equatorial plane, la[45] = 0
    plane = np.log10(np.clip(ne_w[45, :, :].T, 1e-3, None)).astype(np.float32)   # [nalt, nlon+1]
    plane_off = put(buf, plane)

    # Earth TEC
    tec_file = SAMI_TEC_DIR + t.strftime('%Y%m%d_%H%M') + '_tmp_inMLON.npy'
    tec_ok = os.path.isfile(tec_file)
    if tec_ok:
        tec = np.load(tec_file)
        tec = np.concatenate([tec, tec[:, :1]], axis=1)
        tec = np.vstack([np.full(nlon + 1, np.nanmean(tec[0])), tec, np.full(nlon + 1, np.nanmean(tec[-1]))])
    else:
        tec = np.full((nlat + 2, nlon + 1), np.nan)
    tec_off = put(buf, tec.astype(np.float32))

    # isosurfaces above layer K0 (same mapping as sami_isosurface_fig)
    r_iso = (alt_grid[1 + K0: 1 + nalt] + Re) / Re
    logne = np.nan_to_num(np.log10(np.clip(ne_w[:, :, K0:], 1e-3, None)), nan=-3)
    xyz, faces = [], []
    for lv in LEVELS:
        if not (logne.min() < np.log10(lv) < logne.max()):
            xyz.append(np.zeros((0, 3)));  faces.append(np.zeros((0, 3), dtype=np.int64))
            continue
        v, f, _, _ = marching_cubes(logne, np.log10(lv), step_size=STEP_SIZE)
        lat = np.radians(np.interp(v[:, 0], np.arange(nlat), lat_c))
        phi = phi_of(np.interp(v[:, 1], np.arange(nlon + 1), lon_w), ut)
        r   = np.interp(v[:, 2], np.arange(logne.shape[2]), r_iso)
        xyz.append(np.c_[r * np.cos(lat) * np.cos(phi), r * np.cos(lat) * np.sin(phi), r * np.sin(lat)])
        faces.append(f)
    vmax  = max([np.abs(a).max() for a in xyz if a.size] + [1.])
    scale = vmax / 32767.                                             # vertices stored as int16 * scale (Re)
    meshes = []
    for lv, a, f in zip(LEVELS, xyz, faces):
        fdt  = np.uint16 if a.shape[0] < 65536 else np.uint32
        voff = put(buf, np.round(a / scale).astype(np.int16))
        foff = put(buf, f.astype(fdt))
        meshes.append(dict(level=lv, nv=int(a.shape[0]), nf=int(f.shape[0]), voff=voff, foff=foff,
                           fbytes=int(np.dtype(fdt).itemsize)))

    out = os.path.join(FRAME_DIR, fname)
    with open(out, 'wb') as fo:
        fo.write(buf)
    frames.append(dict(file=fname, sami=t.strftime('%Y-%m-%d %H:%M'), ut=round(ut, 5), h=round(h_real, 6),
                       scale=scale, tec=tec_ok,
                       plane=dict(off=plane_off, nr=int(plane.shape[0]), nc=int(plane.shape[1])),
                       tecarr=dict(off=tec_off, nr=int(tec.shape[0]), nc=int(tec.shape[1])),
                       meshes=meshes, bytes=len(buf)))
    if n % 20 == 0 or n == len(sami_t) - 1:
        print('%3d/%d  %s  %.2f MB  (verts %s)' % (n + 1, len(sami_t), fname, len(buf) / 1e6, [m['nv'] for m in meshes]))

with open(os.path.join(DATA_DIR, 'sami_track.json'), 'w') as f:
    json.dump(dict(note='SAMI Ne (/cc) at the Arase GSM position, same 1-min bins as arase.json',
                   ne=[None if not np.isfinite(v) else float('%.4g' % v) for v in ne_trk]), f, separators=(',', ':'))
print('wrote sami_track.json: SAMI Ne valid at %d of %d min' % (np.isfinite(ne_trk).sum(), ne_trk.size))

with open(os.path.join(DATA_DIR, 'manifest.json'), 'w') as f:
    json.dump(dict(created=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), frames=frames), f, separators=(',', ':'))
print('wrote manifest.json: %d frames, %.0f MB in total' % (len(frames), sum(fr['bytes'] for fr in frames) / 1e6))
