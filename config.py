## ===== settings shared by step1_prepare_inputs.py and step2_build_frames.py =====
# Same choices as the "SAMI-Arase 3D" cell of sami_arase_plasmasphere.ipynb.
# After changing anything here, re-run step 1 and step 2.
import os
import datetime

HERE      = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(HERE, 'inputs')            # copied Arase CDFs + SAMI grid / coastline cache
WEB_DIR   = os.path.join(HERE, 'web')               # the website: index.html, viewer.js, plotly.min.js, data/
DATA_DIR  = os.path.join(WEB_DIR, 'data')           # static.json, arase.json, manifest.json, frames/*.bin
FRAME_DIR = os.path.join(DATA_DIR, 'frames')

# ===== SAMI3 run (read only, stays on the external drive) =====
SAMI_ROOT      = '/Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/'
SAMI_GRID_DIR  = SAMI_ROOT + 'SAMIdata/20241011_real/'                 # blonu/blatu/baltu, glonu/glatu/zaltu, run_information.dat
SAMI_TIME_FILES = [SAMI_ROOT + 'SAMIdata/20241010_real/time.dat',
                   SAMI_ROOT + 'SAMIdata/20241011_real/time.dat']
SAMI_BASETIME  = datetime.datetime(2024, 10, 9)                       # time.dat hours are counted from here
SAMI_NE_DIR    = SAMI_ROOT + 'SAMI_equiplane/'                         # yymmdd_HHMM.npy, ne [90 mlat, 180 mlon, 118 alt] (calc_tec)
SAMI_TEC_DIR   = SAMI_ROOT + 'SAMIdata/SAMIdata_npy/'                  # YYYYMMDD_HHMM_tmp_inMLON.npy, TEC [90 mlat, 180 mlon]
Re             = 6370.0
MLON_OFFSET    = 72.                                                   # MLT = (mlon - 72)/15 + UT, as in the notebook

# ===== Arase =====
ARASE_SRC_DIR = '/Users/yulupeng/Documents/AAA_research/Code/arase_pwe/'   # where the notebook keeps the CDFs (copied from here)
T0_STR        = '2024-10-10'                       # first day (00 UT); all "hours" are hours since this time
NDAYS         = 3

# ===== solar wind (as the OMNI section of swmf_omni.ipynb): By, Bz from the SWMF IMF input, E field from OMNI =====
SW_SRC_DIR = '/Users/yulupeng/Documents/AAA_research/Results/2024_1010_storm/20241010_inputs/'   # copied from here
IMF_FILE   = 'imf20241009_12.dat'                  # SWMF IMF input: yyyy mm dd hh mm ss msc bx by bz vx vy vz n t
OMNI_FILE  = 'omni_20241009_12.txt'                # OMNIWeb 1-min listing, column 13 = electric field (mV/m)

# ===== SAMI - real time shift =====
SHIFT = 3                                          # real (Arase) time T <-> SAMI file at T + SHIFT * dt_sami (6-min outputs: +18 min)

# ===== 3D figure =====
LEVELS    = [1000, 500, 50]                        # isosurface densities (/cc)
LV_COLOR  = ['red', 'yellow', 'blue']
LV_ALPHA  = [0.3, 0.3, 0.25]
K0        = 47                                     # isosurfaces only above alt layer k0 (47 = ~500 km)
STEP_SIZE = 2                                      # marching_cubes step: 1 = full mesh (~3x larger files), 2 = every 2nd grid point
LIM       = 5                                      # axis range +/- LIM Re
EQ_RANGE  = [2, 5]                                 # log10(Ne) colour range of the equatorial plane and the Arase dot
EQ_ALPHA  = 0.6                                    # opacity of the equatorial plane
TEC_RANGE = [0, 120]                               # TECU colour range of the Earth
DOT_SIZE  = 10
TRACK_HR  = 3                                      # Arase orbit +/- TRACK_HR hours drawn as a gray line; 0 = none
CAMERA    = dict(eye=dict(x=1.3, y=1.3, z=1.0))    # initial view: Sun (+X) on the left, dusk (+Y) on the right
