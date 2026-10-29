## ===== Step 1: copy the inputs and write the static website data =====
# 1. copy the Arase orbit + HFA L3 CDFs into inputs/arase/ (download them if they are not in ARASE_SRC_DIR)
# 2. read the SAMI grid -> inputs/sami_static.npz (altitude grid, list of SAMI output times)
# 3. coastlines in SAMI magnetic coordinates (same method as coast_mag in the notebook)
# 4. web/data/static.json : grids, coastlines, plot settings
#    web/data/arase.json  : Arase L, MLT, MLAT, L3 Ne, quality flag on the 1-min grid, and the inbound / outbound
#                           legs (arase_inout of the notebook)
#    inputs/arase_track.npz: Arase R, MLAT, MLT (and GSM x, y, z) on the 1-min grid (step 2: SAMI Ne along the orbit)
#    web/data/solarwind.json: IMF By, Bz (SWMF IMF input) and OMNI E field on the 1-min grid
# 5. copy plotly.min.js (from the installed python plotly) into web/, so the page works offline
# Run:  python3 step1_prepare_inputs.py
import os
import sys
import json
import shutil
import datetime
import urllib.request
import numpy as np
import cdflib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
from SAMI_functions import read_info_file, read_time_file, read_coordinates

os.makedirs(os.path.join(INPUT_DIR, 'arase'), exist_ok=True)
os.makedirs(FRAME_DIR, exist_ok=True)

t0    = np.datetime64(T0_STR, 'ns')
dates = [(t0 + i * np.timedelta64(1, 'D')).astype('datetime64[D]').item().strftime('%Y%m%d') for i in range(NDAYS)]


def to_hours(t):
    """datetime64 (or TT2000 from a CDF) -> hours since T0_STR 00 UT"""
    t = np.asarray(t)
    if t.dtype.kind in 'iu':
        t = cdflib.cdfepoch.to_datetime(t)
    return (t.astype('datetime64[ns]') - t0) / np.timedelta64(1, 'h')


def nan_to_none(a, nd=4):
    """numpy array -> list for JSON, rounded, NaN -> null"""
    return [None if not np.isfinite(v) else round(float(v), nd) for v in np.ravel(a)]


## ===== 1. Arase CDFs -> inputs/arase =====
arase_dir = os.path.join(INPUT_DIR, 'arase')
src_url = {'orb': 'https://cdaweb.gsfc.nasa.gov/sp_phys/data/erg/orb/l2/def/%s/',
           'l3':  'https://ergsc.isee.nagoya-u.ac.jp/data/ergsc/satellite/erg/pwe/hfa/l3/%s/%s/'}
for d in dates:
    for kind, f in [('orb', 'erg_orb_l2_%s_v05.cdf' % d), ('l3', 'erg_pwe_hfa_l3_1min_%s_v06_11.cdf' % d)]:
        dst = os.path.join(arase_dir, f)
        if os.path.isfile(dst):
            continue
        if os.path.isfile(ARASE_SRC_DIR + f):
            shutil.copy2(ARASE_SRC_DIR + f, dst)
            print('copied     ' + f)
        else:
            url = (src_url[kind] % d[:4] if kind == 'orb' else src_url[kind] % (d[:4], d[4:6])) + f
            urllib.request.urlretrieve(url, dst)
            print('downloaded ' + f)

## ===== Arase orbit + L3 on the 1-min grid (same as the notebook's L3 cell and SAMI-Arase step 1) =====
tt_orb, rmlatmlt, lm, pos_gsm = [], [], [], []
tt_l3, ne_l3, qf_l3 = [], [], []
for d in dates:
    cdf = cdflib.CDF(os.path.join(arase_dir, 'erg_orb_l2_%s_v05.cdf' % d))
    tt_orb.append(cdflib.cdfepoch.to_datetime(cdf.varget('epoch')))
    rmlatmlt.append(cdf.varget('pos_rmlatmlt'))                        # [time, (R Re, MLAT deg, MLT h)]
    lm.append(cdf.varget('pos_Lm'))                                    # [time, pitch angle]
    pos_gsm.append(cdf.varget('pos_gsm'))                              # [time, 3] Re
    cdf = cdflib.CDF(os.path.join(arase_dir, 'erg_pwe_hfa_l3_1min_%s_v06_11.cdf' % d))
    tt_l3.append(cdf.varget('Epoch'))
    ne_l3.append(cdf.varget('ne_mgf'))
    qf_l3.append(cdf.varget('quality_flag'))
tt_orb    = np.concatenate(tt_orb)
orb_hours = to_hours(tt_orb)
pos_gsm   = np.concatenate(pos_gsm).astype(float)
rmlatmlt  = np.concatenate(rmlatmlt).astype(float)
lm        = np.concatenate(lm).astype(float)
lm[~(lm > 0)] = np.nan
L90       = lm[:, 0]
ne_l3     = np.concatenate(ne_l3).astype(float)
ne_l3[~(ne_l3 > 0)] = np.nan
qf_l3     = np.concatenate(qf_l3).astype(np.int64)

nbin     = NDAYS * 24 * 60
hours_l3 = (np.arange(nbin) + 0.5) / 60.
ibin     = np.floor(to_hours(np.concatenate(tt_l3)) * 60).astype(int)
ok       = (ibin >= 0) & (ibin < nbin)
ne_g     = np.full(nbin, np.nan);  ne_g[ibin[ok]] = ne_l3[ok]
qf_g     = np.full(nbin, -1);      qf_g[ibin[ok]] = qf_l3[ok]          # -1 = no data

gL        = np.isfinite(L90)
L_1min    = np.interp(hours_l3, orb_hours[gL], L90[gL])
mlt_1min  = np.mod(np.degrees(np.interp(hours_l3, orb_hours, np.unwrap(np.radians(rmlatmlt[:, 2] * 15)))) / 15, 24)
mlat_1min = np.interp(hours_l3, orb_hours, rmlatmlt[:, 1])
no_orb = (hours_l3 < orb_hours.min()) | (hours_l3 > orb_hours.max())
r_1min    = np.interp(hours_l3, orb_hours, rmlatmlt[:, 0])                                         # R (Re)
gsm_1min  = np.column_stack([np.interp(hours_l3, orb_hours, pos_gsm[:, k]) for k in range(3)])   # [n, 3] Re
L_1min[no_orb] = np.nan;  mlt_1min[no_orb] = np.nan;  mlat_1min[no_orb] = np.nan;  r_1min[no_orb] = np.nan;  gsm_1min[no_orb] = np.nan
np.savez(os.path.join(INPUT_DIR, 'arase_track.npz'), hours_l3=hours_l3, r_1min=r_1min, mlat_1min=mlat_1min,
         mlt_1min=mlt_1min, gsm_1min=gsm_1min)


## ===== inbound / outbound legs (arase_inout of the notebook, unchanged) =====
def arase_inout(t_start, t_end, tt_orb, rmlatmlt, full_only=False, min_leg=60.):
    """Split an Arase time range into inbound (R decreasing) and outbound (R increasing) legs.
    Returns inbound, outbound [n_leg, 2] of 'YYYYMMDD/HH:MM' (start, end), perigee and apogee times."""
    def tstr(ep):
        return (np.datetime64(ep, 'ns') + np.timedelta64(30, 's')).astype('datetime64[m]').item().strftime('%Y%m%d/%H:%M')

    ep = np.array([np.datetime64(datetime.datetime.strptime(s, '%Y%m%d/%H:%M'), 'ns') for s in (t_start, t_end)])
    r  = rmlatmlt[:, 0]
    ok = np.isfinite(r) & (r > 0)
    t  = tt_orb[ok].astype('datetime64[ns]');  r = r[ok]
    tm = (t - t[0]) / np.timedelta64(1, 'm')

    # sign of dR/dt over +/- 5 min
    dtm = np.median(np.diff(tm))
    k   = max(int(round(5. / dtm)), 1)
    n   = r.size
    dr  = np.zeros(n)
    dr[k:n-k] = r[2*k:] - r[:n-2*k]
    dr[:k]    = dr[k];  dr[n-k:] = dr[n-k-1]
    sgn = np.where(dr > 0, 1, -1)

    # merge legs shorter than min_leg into their neighbours
    while True:
        chg = np.nonzero(sgn[1:] != sgn[:-1])[0]
        if chg.size == 0:
            break
        b0 = np.r_[0, chg + 1];  b1 = np.r_[chg, n - 1]
        dur = tm[b1] - tm[b0]
        ishort = np.nonzero(dur < min_leg)[0]
        if ishort.size == 0:
            break
        j = ishort[np.argmin(dur[ishort])]
        sgn[b0[j]:b1[j]+1] = -sgn[b0[j]]

    # turning points: min / max of R within +/- 30 min of each sign change
    chg  = np.nonzero(sgn[1:] != sgn[:-1])[0]
    nchg = chg.size
    tp   = np.zeros(nchg, dtype=int)
    ispg = np.zeros(nchg, dtype=bool)
    for i, c in enumerate(chg):
        w = np.nonzero(np.abs(tm - tm[c]) <= 30)[0]
        if sgn[c] < 0:
            im = np.argmin(r[w]);  ispg[i] = True
        else:
            im = np.argmax(r[w])
        tp[i] = w[im]

    # legs between turning points, clipped to [t_start, t_end]
    edges   = np.concatenate([t[:1], t[tp], t[-1:]])
    full    = [False] + [True] * (nchg - 1) + [False] if nchg > 0 else [False]
    leg_out = np.r_[~ispg[0], ispg] if nchg > 0 else np.array([sgn[0] > 0])
    inbound, outbound = [], []
    for i in range(edges.size - 1):
        a = max(edges[i], ep[0]);  b = min(edges[i+1], ep[1])
        if b <= a:
            continue
        cut = (a != edges[i]) or (b != edges[i+1]) or not full[i]
        if full_only and cut:
            continue
        (outbound if leg_out[i] else inbound).append([tstr(a), tstr(b)])
    inbound  = np.array(inbound,  dtype=str).reshape(-1, 2)
    outbound = np.array(outbound, dtype=str).reshape(-1, 2)

    perigee, apogee = [], []
    for i in range(nchg):
        if t[tp[i]] < ep[0] or t[tp[i]] > ep[1]:
            continue
        (perigee if ispg[i] else apogee).append(tstr(t[tp[i]]))
    return inbound, outbound, perigee, apogee


t_last = (t0 + NDAYS * np.timedelta64(1, 'D') - np.timedelta64(1, 'm')).astype('datetime64[m]').item().strftime('%Y%m%d/%H:%M')
inbound, outbound, _, _ = arase_inout(dates[0] + '/00:00', t_last, tt_orb, rmlatmlt)
str_hours = lambda s: to_hours(np.datetime64(datetime.datetime.strptime(s, '%Y%m%d/%H:%M'), 'ns'))
legs = []
for kind, arr in [('inbound', inbound), ('outbound', outbound)]:
    for a, b in arr:
        ii = np.nonzero((hours_l3 >= str_hours(a)) & (hours_l3 <= str_hours(b)))[0]   # 1-min bins in the leg, as the notebook
        if ii.size:
            legs.append(dict(dir=kind, start=a, end=b, i0=int(ii[0]), i1=int(ii[-1])))
legs.sort(key=lambda l: l['i0'])

with open(os.path.join(DATA_DIR, 'arase.json'), 'w') as f:
    json.dump(dict(t0=T0_STR, dt_min=1, n=nbin, note='bin i is centred at t0 + (i + 0.5) min',
                   L=nan_to_none(L_1min, 3), mlt=nan_to_none(mlt_1min, 3), mlat=nan_to_none(mlat_1min, 3),
                   ne=nan_to_none(ne_g, 3), qf=[int(q) for q in qf_g], legs=legs), f, separators=(',', ':'))
print('Arase: %d min, Ne valid %d, L valid %d, %d legs' % (nbin, np.isfinite(ne_g).sum(), np.isfinite(L_1min).sum(), len(legs)))

## ===== solar wind: IMF By, Bz and OMNI E field (readers of swmf_omni_functions.py, NaN for fill values) =====
sw_dir = os.path.join(INPUT_DIR, 'solarwind')
os.makedirs(sw_dir, exist_ok=True)
for f in (IMF_FILE, OMNI_FILE):
    if not os.path.isfile(os.path.join(sw_dir, f)):
        shutil.copy2(SW_SRC_DIR + f, os.path.join(sw_dir, f))
        print('copied     ' + f)

t_imf, by, bz = [], [], []
with open(os.path.join(sw_dir, IMF_FILE)) as fi:
    for line in fi:
        if line.startswith('#START'):
            break
    for line in fi:
        a = line.split()
        if len(a) < 10:
            continue
        t_imf.append(np.datetime64('%s-%s-%sT%s:%s:%s' % tuple(a[:6])))
        by.append(float(a[8]));  bz.append(float(a[9]))
t_omni, ef = [], []
with open(os.path.join(sw_dir, OMNI_FILE)) as fi:
    data = False
    for line in fi:
        a = line.split()
        if line.startswith('YYYY'):
            data = True
            continue
        if not data or len(a) < 17 or line.startswith('<'):
            continue
        t_omni.append(np.datetime64('%s-01-01' % a[0]) + np.timedelta64(int(a[1]) - 1, 'D') +
                      np.timedelta64(int(a[2]), 'h') + np.timedelta64(int(a[3]), 'm'))
        ef.append(np.nan if a[16] == '999.99' else float(a[16]))


def on_minutes(t, v):
    """values at whole minutes t0 + i min, i = 0 .. nbin-1 (NaN where missing)"""
    out = np.full(nbin, np.nan)
    i = np.round(to_hours(np.array(t)) * 60).astype(int)
    ok = (i >= 0) & (i < nbin)
    out[i[ok]] = np.asarray(v, float)[ok]
    return out


with open(os.path.join(DATA_DIR, 'solarwind.json'), 'w') as f:
    json.dump(dict(t0=T0_STR, n=nbin, note='value i is at t0 + i min; By, Bz (nT, GSM) from ' + IMF_FILE +
                   ', E (mV/m) from ' + OMNI_FILE,
                   by=nan_to_none(on_minutes(t_imf, by), 2), bz=nan_to_none(on_minutes(t_imf, bz), 2),
                   ef=nan_to_none(on_minutes(t_omni, ef), 2)), f, separators=(',', ':'))
print('solar wind: IMF %d min, OMNI %d min' % (len(t_imf), len(t_omni)))

## ===== 2. SAMI grid and output times =====
Grid, Vars = read_info_file(SAMI_GRID_DIR + 'run_information.dat')
nz, nf, nlt = Grid['Points'], Grid['nFL'], Grid['nLons']
bLon, bLat, bAlt = read_coordinates(SAMI_GRID_DIR, 'blonu.dat', 'blatu.dat', 'baltu.dat', nf, nz, nlt)
bAlt = bAlt - Re
gLon, gLat, zAlt = read_coordinates(SAMI_GRID_DIR, 'glonu.dat', 'glatu.dat', 'zaltu.dat', nf, nz, nlt)
alt_grid = np.max(bAlt, axis=2)[0] - 1                                 # km, layer k of the ne cube is at alt_grid[k+1]

times = np.concatenate([read_time_file(f, SAMI_BASETIME) for f in SAMI_TIME_FILES])
sami_t = [t for t in times if os.path.isfile(SAMI_NE_DIR + t.strftime('%y%m%d_%H%M') + '.npy')]
dt_sami = np.median(np.diff(to_hours(np.array(times, dtype='datetime64[ns]'))))
np.savez(os.path.join(INPUT_DIR, 'sami_static.npz'), alt_grid=alt_grid, dt_sami=dt_sami,
         sami_t=np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in sami_t]))
print('SAMI outputs: %d, %s - %s, dt = %g min' % (len(sami_t), sami_t[0], sami_t[-1], dt_sami * 60))

## ===== 3. coastlines in SAMI magnetic coordinates (coast_mag of the notebook) =====
from scipy.spatial import cKDTree
import cartopy.io.shapereader as shpreader
uv = lambda lat, lon: np.c_[np.cos(np.radians(lat)) * np.cos(np.radians(lon)),
                            np.cos(np.radians(lat)) * np.sin(np.radians(lon)), np.sin(np.radians(lat))]
m    = zAlt < 300.
tree = cKDTree(uv(gLat[m], gLon[m]))
bvec = uv(bLat[m], bLon[m])
ll = []
for geom in shpreader.Reader(shpreader.natural_earth('110m', 'physical', 'coastline')).geometries():
    for line in getattr(geom, 'geoms', [geom]):
        lon, lat = np.asarray(line.coords).T
        ll += [np.c_[lat, lon], [[np.nan, np.nan]]]
ll  = np.vstack(ll)
cl  = np.full(ll.shape, np.nan)
ok  = np.isfinite(ll[:, 0])
d, i = tree.query(uv(ll[ok, 0], ll[ok, 1]), k=4)
w   = 1. / np.maximum(d, 1e-9)
v   = np.sum(bvec[i] * w[..., None], axis=1)
v  /= np.linalg.norm(v, axis=1, keepdims=True)
mag = np.c_[np.degrees(np.arcsin(v[:, 2])), np.mod(np.degrees(np.arctan2(v[:, 1], v[:, 0])), 360)]
mag[(d[:, 0] > np.radians(3)) | (np.abs(mag[:, 0]) > np.abs(bLat[m]).max() - 1)] = np.nan
cl[ok] = mag

## ===== 4. static.json: grids + settings for the viewer =====
nlat, nlon, nalt = np.load(SAMI_NE_DIR + sami_t[0].strftime('%y%m%d_%H%M') + '.npy', mmap_mode='r').shape
lat_c = np.linspace(-85, 85, nlat + 1)[:-1] + 85. / nlat               # bin centres used by calc_tec
lon_c = np.linspace(0, 360, nlon + 1)[:-1] + 180. / nlon
lon_w = np.r_[lon_c, lon_c[0] + 360]                                   # closed across 0/360 deg
static = dict(
    t0=T0_STR, shift=SHIFT, dt_sami_min=round(float(dt_sami * 60), 4), mlon_offset=MLON_OFFSET,
    lon=nan_to_none(lon_w), lat_earth=nan_to_none(np.r_[-90, lat_c, 90]),
    r_plane=nan_to_none((alt_grid[1:nalt + 1] + Re) / Re, 5),
    coast_mlat=nan_to_none(cl[:, 0], 3), coast_mlon=nan_to_none(cl[:, 1], 3),
    settings=dict(levels=LEVELS, colors=LV_COLOR, alphas=LV_ALPHA, k0=K0, step_size=STEP_SIZE, lim=LIM,
                  eq_range=EQ_RANGE, eq_alpha=EQ_ALPHA, tec_range=TEC_RANGE, dot_size=DOT_SIZE,
                  track_hr=TRACK_HR, camera=CAMERA))
with open(os.path.join(DATA_DIR, 'static.json'), 'w') as f:
    json.dump(static, f, separators=(',', ':'))
print('wrote ' + os.path.join(DATA_DIR, 'static.json'))

## ===== 5. plotly.js for the page =====
import plotly
shutil.copy2(os.path.join(os.path.dirname(plotly.__file__), 'package_data', 'plotly.min.js'), os.path.join(WEB_DIR, 'plotly.min.js'))
print('copied plotly.min.js (python plotly %s)' % plotly.__version__)
