## ===== settings shared by step1_prepare_inputs.py and step2_build_frames.py =====
# Same choices as the "SAMI-Arase 3D" cell of sami_arase_plasmasphere.ipynb.
# One entry per event in EVENTS; the scripts take the event key:
#     python3 step1_prepare_inputs.py --event 20230423
#     python3 step2_build_frames.py   --event 20230423
# After changing anything here, re-run step 1 and step 2 for that event.
import os
import json
import datetime

HERE      = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(HERE, 'inputs')            # copied Arase CDFs (inputs/arase), solar wind files (inputs/solarwind),
                                                    # per-event SAMI grid / Arase track (inputs/<event>/)
WEB_DIR   = os.path.join(HERE, 'web')               # the website: index.html, viewer.js, plotly.min.js, data/
Re          = 6370.0
ARASE_SRC_DIR = '/Users/yulupeng/Documents/AAA_research/Code/arase_pwe/'   # where the notebook keeps the CDFs (copied from here)

# ===== events =====
# SAMI_TIMES: list of SAMI output times, either ('time.dat', [files], basetime) or ('npy', file) (np.save'd datetimes)
# SAMI_NE_DIR : yymmdd_HHMM.npy, ne [nlat mlat, nlon mlon, nalt alt] from calc_tec (grid sizes are read from the files)
# SAMI_TEC_DIR: YYYYMMDD_HHMM_tmp_inMLON.npy, TEC [nlat mlat, nlon mlon] (may have its own grid)
# T0_STR, NDAYS: Arase / solar wind days loaded; all "hours" are hours since T0_STR 00 UT
# SHIFT: real (Arase) time T <-> SAMI file at T + SHIFT * dt_sami
EVENTS = {
    '20241010': dict(
        TITLE         = 'October 2024 storm',
        DESCRIPTION   = 'SAMI3 (GITM-SAMI, 180 mlon) 2024-10-10 to 10-11, 6-min outputs, shift 3 steps (+18 min)',
        SAMI_ROOT     = '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/',
        SAMI_GRID_DIR = '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/20241011_real/',
        SAMI_TIMES    = ('time.dat', ['/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/20241010_real/time.dat',
                                      '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/20241011_real/time.dat'],
                         datetime.datetime(2024, 10, 9)),
        SAMI_NE_DIR   = '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMI_equiplane/',
        SAMI_TEC_DIR  = '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/SAMIdata_npy/',
        MLON_OFFSET   = 72.,                               # MLT = (mlon - 72)/15 + UT, as in the notebook
        T0_STR        = '2024-10-10',
        NDAYS         = 3,
        SW_SRC_DIR    = '/Users/yulupeng/Documents/AAA_research/Results/2024_1010_storm/20241010_inputs/',
        IMF_FILE      = 'imf20241009_12.dat',              # SWMF IMF input: yyyy mm dd hh mm ss msc bx by bz vx vy vz n t
        OMNI_FILE     = 'omni_20241009_12.txt',            # OMNIWeb 1-min listing, column 13 = electric field (mV/m)
        SHIFT         = 3,
    ),
    '20230423': dict(
        TITLE         = 'April 2023 storm',
        DESCRIPTION   = 'SAMI3 (GITM-SAMI high-res, 720 mlon) 2023-04-23 to 04-24, 6-min outputs, no shift',
        SAMI_ROOT     = '/Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/',
        SAMI_GRID_DIR = '/Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/20230423_1/',
        SAMI_TIMES    = ('npy', '/Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/SAMI_equiplane/timesall.npy'),
        SAMI_NE_DIR   = '/Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/SAMI_equiplane/',
        SAMI_TEC_DIR  = '/Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/SAMI_tec/',
        MLON_OFFSET   = 72.,
        T0_STR        = '2023-04-23',
        NDAYS         = 3,
        SW_SRC_DIR    = '/Users/yulupeng/Documents/AAA_research/Results/2023_0423_storm/20230423_inputs/',
        IMF_FILE      = 'imf20230423.dat',
        OMNI_FILE     = 'omni_20230423.txt',
        SHIFT         = 0,
    ),
}
DEFAULT_EVENT = '20241010'

# ===== 3D figure (all events) =====
LEVELS    = [1000, 500, 50]                        # isosurface densities (/cc)
LV_COLOR  = ['red', 'yellow', 'blue']
LV_ALPHA  = [0.3, 0.3, 0.25]
K0        = 47                                     # isosurfaces only above alt layer k0 (47 = ~500 km)
STEP_SIZE = 2                                      # marching_cubes step: 1 = full mesh (~3x larger files), 2 = every 2nd grid point
MAX_NLON  = 180                                    # website mesh resolution in mlon: grids with more mlon (e.g. 720) are averaged
                                                   # to MAX_NLON bins for the isosurfaces / plane / Earth only (SAMI Ne along Arase
                                                   # uses the full grid). None = full resolution (much larger frames, slower)
LIM       = 5                                      # axis range +/- LIM Re
EQ_RANGE  = [2, 5]                                 # log10(Ne) colour range of the equatorial plane and the Arase dot
EQ_ALPHA  = 0.6                                    # opacity of the equatorial plane
TEC_RANGE = [0, 120]                               # TECU colour range of the Earth
DOT_SIZE  = 10
TRACK_HR  = 3                                      # Arase orbit +/- TRACK_HR hours drawn as a gray line; 0 = none
CAMERA    = dict(eye=dict(x=1.3, y=1.3, z=1.0))    # initial view: Sun (+X) on the left, dusk (+Y) on the right


def lon_factor(nlon):
    """number of neighbouring mlon bins averaged for the website meshes, so that at most MAX_NLON bins are left"""
    if not MAX_NLON or nlon <= MAX_NLON:
        return 1
    f = -(-nlon // MAX_NLON)
    while nlon % f:                                # use a divisor of nlon, so every website bin has the same width
        f += 1
    return f


def lon_reduce(a, axis=1):
    """mean of lon_factor neighbouring mlon bins of a along axis (NaN-aware); bin centres stay calc_tec bin centres"""
    import numpy as np
    f = lon_factor(a.shape[axis])
    if f == 1:
        return a
    a = np.moveaxis(a, axis, 0)
    a = np.nanmean(a.reshape((a.shape[0] // f, f) + a.shape[1:]), axis=1)
    return np.moveaxis(a, 0, axis)


def load_event(name):
    """settings of one event as a dict (EVENTS[name] + its folders); scripts put these into their globals"""
    if name not in EVENTS:
        raise SystemExit('unknown event %r, choose one of %s' % (name, ', '.join(EVENTS)))
    e = dict(EVENTS[name])
    e.update(EVENT=name,
             EVENT_INPUT_DIR=os.path.join(INPUT_DIR, name),             # sami_static.npz, arase_track.npz
             DATA_DIR=os.path.join(WEB_DIR, 'data', name),              # static.json, arase.json, manifest.json, ...
             FRAME_DIR=os.path.join(WEB_DIR, 'data', name, 'frames'))   # frames/*.bin
    return e


def write_events_json():
    """web/data/events.json: the events the page offers (those with a manifest.json), newest first"""
    out = []
    for name, e in EVENTS.items():
        man = os.path.join(WEB_DIR, 'data', name, 'manifest.json')
        if not os.path.isfile(man):
            continue
        with open(man) as f:
            fr = json.load(f)['frames']
        t0 = datetime.datetime.strptime(e['T0_STR'], '%Y-%m-%d')
        real = lambda h: (t0 + datetime.timedelta(hours=h)).strftime('%Y-%m-%d %H:%M')   # real (Arase) UT of a frame
        out.append(dict(id=name, title=e['TITLE'], description=e['DESCRIPTION'], nframes=len(fr),
                        first=real(fr[0]['h']), last=real(fr[-1]['h']), shift=e['SHIFT']))
    out.sort(key=lambda d: d['id'], reverse=True)
    with open(os.path.join(WEB_DIR, 'data', 'events.json'), 'w') as f:
        json.dump(dict(default=DEFAULT_EVENT, events=out), f, indent=1)
    print('wrote events.json: ' + ', '.join(d['id'] for d in out))
