SAMI3 + Arase 3D website
========================

An interactive web page that steps through every SAMI3 output in time ("a video") and keeps
the 3D Plotly view interactive: drag the time bar or press Play, stop anywhere, then rotate,
zoom and hide layers as in the notebook. The camera stays where you leave it while the time
changes.

Each time step shows (same as the "SAMI-Arase 3D" cell of sami_arase_plasmasphere.ipynb):
  - SAMI3 Ne isosurfaces at 1000 / 500 / 50 /cc (above ~500 km, k0 = 47)
  - the equatorial plane ne[45, :, :] in log10(Ne), jet colormap, 10^0 - 10^5 /cc
  - the Earth coloured by the SAMI TEC map (jet, 0 - 120 TECU) with white coastlines
  - Arase as a dot at its L, MLT, MLAT (dipole: r = L cos^2(MLAT)), coloured by the Arase
    HFA L3 Ne on the same jet scale as the plane, and the Arase orbit +/- 3 h in gray
  - right of the 3D view: Ne vs L along the current Arase leg (inbound / outbound, from
    arase_inout of the notebook): SAMI3 (solid blue) and Arase L3 (dotted red) for the whole
    leg, with a dot on each line at the current time that moves along the lines as time
    advances; the panel switches to the next leg when the time enters it. SAMI3 Ne here is the
    notebook's SAMI-Arase step 2 (sami_ne_gsm): the SAMI cube taken as GSM lat / lon / alt,
    sampled at the Arase GSM position, linear in time between SAMI outputs (with SHIFT)
  - below (scroll down; the time controls stay pinned at the top of the window): three panels
    on one shared time axis, as in the OMNI section of swmf_omni.ipynb:
      IMF By and Bz (nT, GSM, from the SWMF IMF input imf20241009_12.dat)
      solar wind electric field (mV/m, OMNI 1-min, omni_20241009_12.txt; 999.99 = gap)
      the Arase L3 Ne time series with a red cursor at the current time
    (click on it to jump to that time)
Frame: +X = noon (Sun), +Y = dusk, Z up. Initial view: Sun on the left, dusk on the right.
Time: the title shows the real (Arase) UT. SAMI file at real time T is the one at T + SHIFT
steps (SHIFT = 3, 6-min outputs -> +18 min), as in the notebook.

The notebook sami_arase_plasmasphere.ipynb is not used or changed by any of this; the code
was copied from it into the scripts below.


Folder contents
---------------
  README.txt                 this file
  config.py                  all paths and plot settings (edit here)
  SAMI_functions.py          copy of ../SAMI_functions.py (SAMI grid / time readers)
  step1_prepare_inputs.py    step 1: copy Arase data, read SAMI grid, write static website data
  step2_build_frames.py      step 2: one data file per SAMI time step
  inputs/
    arase/                   copied Arase CDFs: erg_orb_l2_YYYYMMDD_v05.cdf,
                             erg_pwe_hfa_l3_1min_YYYYMMDD_v06_11.cdf (2024-10-10 .. 12)
    sami_static.npz          SAMI altitude grid, list of SAMI output times, output cadence
    arase_track.npz          Arase GSM x, y, z on the 1-min grid (for step 2)
    solarwind/               copied IMF and OMNI files (IMF_FILE, OMNI_FILE in config.py)
  web/                       the website (this is the folder that is served)
    index.html               page layout and controls
    viewer.js                loads the data and draws the figures
    plotly.min.js            copy of plotly.js from the installed python plotly (works offline)
    data/static.json         grids, coastlines in SAMI mlat/mlon, plot settings
    data/arase.json          Arase L, MLT, MLAT, L3 Ne, quality flag on a 1-min grid, legs
    data/sami_track.json     SAMI Ne at the Arase position on the same 1-min grid
    data/solarwind.json      IMF By, Bz and OMNI E field on a 1-min grid
    data/manifest.json       list of time steps and the byte layout of each frame file
    data/frames/*.bin        one file per SAMI output, named by SAMI time yymmdd_HHMM.bin
                             (480 files, ~0.55 MB each, ~270 MB in total)

Not copied (too large, ~4 GB): the SAMI cubes and TEC maps. They are read from the external
drive in step 2 only; the website itself does not need the drive.
  SAMI Ne cubes : /Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMI_equiplane/yymmdd_HHMM.npy
  SAMI TEC maps : /Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/SAMIdata_npy/YYYYMMDD_HHMM_tmp_inMLON.npy
  SAMI grid     : /Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/SAMIdata/20241011_real/


Procedure
---------
Python: the anaconda python3 (needs numpy, scipy, scikit-image, cdflib, cartopy, plotly).
The external drive ExtremePro must be mounted for steps 1 and 2.

Step 0. (optional) edit config.py
        Paths, SHIFT, isosurface levels / colours / opacity, K0, STEP_SIZE (mesh detail),
        axis range, colour ranges, plane opacity, dot size, orbit track length, initial camera.

Step 1. Prepare inputs and static data (~1 min):
            cd /Users/yulupeng/Documents/AAA_research/Code/Arase_SAMI_website
            python3 step1_prepare_inputs.py
        - copies the Arase orbit and HFA L3 CDFs from ../arase_pwe/ into inputs/arase/
          (downloads them if they are missing there)
        - puts Arase L (90 deg pitch angle), MLT, MLAT on the 1-min L3 grid, and L3 Ne,
          quality flag, and the inbound / outbound legs -> web/data/arase.json
        - Arase GSM position -> inputs/arase_track.npz
        - copies the IMF and OMNI files from SW_SRC_DIR into inputs/solarwind/ and writes
          By, Bz, E field -> web/data/solarwind.json
        - reads the SAMI grid and time.dat files -> inputs/sami_static.npz
        - maps the Natural Earth coastlines to SAMI mlat/mlon with the SAMI grid itself
          (same method as coast_mag in the notebook)
        - writes web/data/static.json and copies plotly.min.js into web/

Step 2. Build one frame per SAMI output (~1 min for all 480):
            python3 step2_build_frames.py
        For every SAMI output:
        - loads the Ne cube, closes it across 0/360 deg mlon
        - isosurfaces (skimage marching_cubes) of log10(Ne) above layer K0, vertices mapped
          to x, y, z (Re) at that UT with MLT = (mlon - 72)/15 + UT, noon at +X
        - equatorial plane log10(ne[45, :, :]) and the TEC map (poles filled)
        - SAMI Ne at the Arase GSM position for every minute between this and the previous output
        - writes web/data/frames/yymmdd_HHMM.bin and, at the end, web/data/sami_track.json and
          web/data/manifest.json
        Quick test with fewer frames:  python3 step2_build_frames.py --every 10
        (this rewrites manifest.json with only those frames; run without --every for all)

Step 3. Open the website. Browsers do not let a page opened as a file read the data files,
        so serve the web/ folder with python's built-in web server:
            cd /Users/yulupeng/Documents/AAA_research/Code/Arase_SAMI_website/web
            python3 -m http.server 8000
        then open  http://localhost:8000  in a browser. Stop the server with Ctrl-C.

Controls
--------
  Play / Pause, space bar       play through time (speed menu: 1 - 12 frames per second;
                                drawing one frame takes ~0.1 s, so ~10 fps is the real maximum)
  < >, left / right arrow keys  one SAMI step back / forward
  time bar                      drag to any time; the newest position is drawn as frames load
  click on any time series      jump to that time (drag = zoom all three, double-click = reset)
  Ne vs L panel                 drag = zoom, double-click = reset; hover a line for its UT
  mouse in the 3D view          drag = rotate, scroll = zoom, right-drag = pan
  click a legend item           hide / show that layer (stays hidden while time changes)
  camera icon (top right)       save the current view as a PNG
  home icon (top right)         back to the initial view
Frames are loaded only when needed (the next 6 are loaded ahead while playing), so the
page opens quickly even though all frames together are ~270 MB.


Sharing
-------
Copy the whole web/ folder (with data/) to another computer and run step 3 there, or put
web/ on any static web host (GitHub Pages, a university web space): no server code is needed.

GitHub Pages: this folder is the git repository github.com/PengYululu/Arase_SAMI_website.
Pages serves the main branch root; index.html there forwards to web/, so the site is
    https://pengyululu.github.io/Arase_SAMI_website/
To update it after re-running step 1 / step 2 or editing the page:
    git add -A && git commit -m "update" && git push
The site refreshes about a minute after the push. .nojekyll tells Pages to serve the files as is.


Notes
-----
  - STEP_SIZE = 2 uses every 2nd grid point for the isosurfaces (~4x fewer triangles than
    the notebook's full mesh). STEP_SIZE = 1 matches the notebook exactly but makes the
    frames ~3x larger (~800 MB) and playback slower.
  - The Arase dot is at the 1-min bin round(real hours * 60), as in the notebook. Before
    2024-10-10 00:00 (first SAMI steps after the shift) there is no Arase data and no dot.
  - The equatorial plane and the Earth are drawn as Plotly mesh3d (coloured triangles), not
    surface traces: surfaces took ~1.3 s per frame to redraw, mesh3d ~0.1 s. Hovering on the
    plane shows log10(Ne).
  - If a TEC map is missing, that frame's Earth is drawn white (dayside) / black (nightside).
  - After changing config.py, re-run step 1 and step 2, then reload the page (if the browser
    shows old data, force-reload: Cmd-Shift-R).
