SAMI3 + Arase 3D website
========================

An interactive web page that steps through every SAMI3 output in time ("a video") and keeps
the 3D Plotly view interactive. The page opens with a chooser of the events (storms) that have
been built; the menu at the top switches events later, and  ?event=<id>  links straight to one:
  20241010   October 2024 storm (SAMI 180 mlon, 6-min outputs, SHIFT = 3)
  20230423   April 2023 storm   (SAMI high-res 720 mlon, 6-min outputs, SHIFT = 0)
For each event: drag the time bar or press Play, stop anywhere, then rotate,
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
    notebook's SAMI-Arase step 2 (sami_ne_mag): Arase MLAT, MLT, R mapped to SAMI mlat, mlon, alt
    (MLT -> mlon = 15 (MLT - UT) + 72 at each SAMI output's UT), linear in time between outputs (with SHIFT)
  - below (scroll down; the time controls stay pinned at the top of the window): three panels
    on one shared time axis, as in the OMNI section of swmf_omni.ipynb:
      IMF By and Bz (nT, GSM, from the SWMF IMF input, IMF_FILE of the event)
      solar wind electric field (mV/m, OMNI 1-min, OMNI_FILE of the event; 999.99 = gap)
      the Arase L3 Ne time series with a red cursor at the current time
    (click on it to jump to that time)
  - optional, two boxes above the 3D view (October 2024 event; built by step 3; the files are
    fetched only while a box is ticked, Arase is always drawn):
      add contour mapping  the equatorial-plane Ne contours at 50 / 500 / 1000 /cc (pink / orange /
                           red, above layer K0) followed along the SAMI field line with its apex at
                           each contour point to its north and south ends -> footprints on the TEC
                           Earth (notebook: add_pla_mapping, "Plasmasphere Contour mapping to TEC")
      add SED edge         SED edge points from the SAMI TEC gradient (|grad TEC| > 0.03 TECU/km,
                           MLON 200-360, 30 <= |MLAT| <= 80; middle of each run along an mlat row)
                           on the Earth (north pink, south green); each one -> the SAMI field line
                           whose foot point at 350 km is closest -> that line's apex on the
                           equatorial plane (diamond), with the field line drawn
                           (notebook: sami_sed_edge, add_sed_edge, add_edge_pla)
    ?contour=1 and ?sed=1 in the address open the page with the box ticked (the address updates
    when you click, so the link can be shared).
  - right panel, switch at its top right: "Arase leg" (the Ne vs L panel above) or "2D mapping":
      top    : equatorial plane log10 Ne seen from the north (noon up, dawn right, a disc out to 8.6 Re, light white L-shell circles L = 2 - 8,
               Earth white day / black night), the 50 / 500 / 1000 /cc contours, and each SED edge
               point at the apex of its field line (triangles: north pink, south green)
      bottom : SAMI TEC map (MLON x MLAT, jet 0 - 120 TECU) with white coastlines, the footprints of
               the 50 / 500 / 1000 /cc contours, the SED edge points (same triangles) and the SED
               search region (dashed)
    (same as the 2D cell of the notebook, with the polar view on top). The plane and TEC map come from
    the frame files; contours, footprints and SED from the step-3 overlay files, which this panel
    loads whether or not the 3D boxes are ticked. Events without step 3 show the plane and TEC only.
    ?panel=map in the address opens the page on the 2D mapping.
Frame: +X = noon (Sun), +Y = dusk, Z up. Initial view: Sun on the left, dusk on the right.
Time: the title shows the real (Arase) UT. SAMI file at real time T is the one at T + SHIFT
steps (2024: SHIFT = 3, 6-min outputs -> +18 min; 2023: SHIFT = 0), as in the notebook.
Grids: the mlat / mlon / alt sizes of the SAMI Ne cube and of the TEC map are read from the files
(2024: 180 mlon, 2023: 720 mlon). For the website meshes (isosurfaces, plane, Earth) mlon is
averaged to at most MAX_NLON = 180 bins, so frames stay ~0.55 MB; SAMI Ne along the Arase orbit
(Ne vs L panel) uses the full grid.

The notebook sami_arase_plasmasphere.ipynb is not used or changed by any of this; the code
was copied from it into the scripts below.


Folder contents
---------------
  README.txt                 this file
  config.py                  EVENTS (paths, times, SHIFT, solar wind files per event) and plot settings
  SAMI_functions.py          copy of ../SAMI_functions.py (SAMI grid / time readers)
  step1_prepare_inputs.py    step 1: copy Arase data, read SAMI grid, write static website data
  step2_build_frames.py      step 2: one data file per SAMI time step
  step3_build_overlays.py    step 3 (optional): contour-mapping and SED-edge overlays per SAMI time step
  inputs/
    arase/                   copied Arase CDFs (all events): erg_orb_l2_YYYYMMDD_v05.cdf,
                             erg_pwe_hfa_l3_1min_YYYYMMDD_v06_11.cdf
    solarwind/               copied IMF and OMNI files (IMF_FILE, OMNI_FILE of each event)
    <event>/sami_static.npz  SAMI altitude grid, list of SAMI output times, output cadence
    <event>/arase_track.npz  Arase R, MLAT, MLT (and GSM) on the 1-min grid (for step 2)
  web/                       the website (this is the folder that is served)
    index.html               page layout, event chooser and controls
    viewer.js                loads the data and draws the figures
    plotly.min.js            copy of plotly.js from the installed python plotly (works offline)
    data/events.json         the events offered by the page (written by step 1 / step 2)
    data/<event>/static.json      grids, coastlines in SAMI mlat/mlon, plot settings
    data/<event>/arase.json       Arase L, MLT, MLAT, L3 Ne, quality flag on a 1-min grid, legs
    data/<event>/sami_track.json  SAMI Ne at the Arase position on the same 1-min grid
    data/<event>/solarwind.json   IMF By, Bz and OMNI E field on a 1-min grid
    data/<event>/manifest.json    list of time steps and the byte layout of each frame file
    data/<event>/frames/*.bin     one file per SAMI output, named by SAMI time yymmdd_HHMM.bin
                                  (~0.55 MB each; 2024: 480 files ~260 MB, 2023: 387 files ~215 MB)
    data/<event>/overlays.json    (step 3) overlay settings and which frames have overlay files
    data/<event>/overlays/contour/yymmdd_HHMM.json, overlays/sed/yymmdd_HHMM.json
                                  (step 3) x, y, z (Re) of contours, footprints, edge points, apexes and
                                  field lines in that frame (2024: ~60 + ~30 kB each, 42 MB in total)

Not copied (too large: ~4 GB for 2024, ~22 GB for 2023): the SAMI cubes and TEC maps. They are
read from the external drive in steps 1 and 2 only; the website itself does not need the drive.
  2024 (drive ExtremePro): /Volumes/ExtremePro/GITMSAMI_20241010_compass_grid/
      SAMI_equiplane/yymmdd_HHMM.npy, SAMIdata/SAMIdata_npy/YYYYMMDD_HHMM_tmp_inMLON.npy, grid SAMIdata/20241011_real/
  2023 (drive 1T): /Volumes/1T/GITMSAMI_20230422_hires/SAMIdata_v3/
      SAMI_equiplane/yymmdd_HHMM.npy (times: SAMI_equiplane/timesall.npy), SAMI_tec/YYYYMMDD_HHMM_tmp_inMLON.npy,
      grid 20230423_1/


Procedure
---------
Python: the anaconda python3 (needs numpy, scipy, scikit-image, cdflib, cartopy, plotly).
The external drive of the event must be mounted for steps 1 and 2. Every step takes the event
key (EVENTS in config.py); run steps 1 and 2 once per event.

Step 0. (optional) edit config.py
        A new event: add an entry to EVENTS (SAMI paths, SAMI_TIMES, T0_STR, NDAYS, solar wind
        files, SHIFT, TITLE, DESCRIPTION), then steps 1 and 2 with --event <key>.
        Paths, SHIFT, isosurface levels / colours / opacity, K0, STEP_SIZE (mesh detail),
        axis range, colour ranges, plane opacity, dot size, orbit track length, initial camera.

Step 1. Prepare inputs and static data (~1 min):
            cd /Users/yulupeng/Documents/AAA_research/Code/Arase_SAMI_website
            python3 step1_prepare_inputs.py --event 20230423
        (below, web/data/ means web/data/<event>/)
        - copies the Arase orbit and HFA L3 CDFs from ../arase_pwe/ into inputs/arase/
          (downloads them if they are missing there)
        - puts Arase L (90 deg pitch angle), MLT, MLAT on the 1-min L3 grid, and L3 Ne,
          quality flag, and the inbound / outbound legs -> web/data/arase.json
        - Arase R, MLAT, MLT -> inputs/<event>/arase_track.npz
        - copies the IMF and OMNI files from SW_SRC_DIR into inputs/solarwind/ and writes
          By, Bz, E field -> web/data/solarwind.json
        - reads the SAMI grid and the SAMI output times -> inputs/<event>/sami_static.npz
        - maps the Natural Earth coastlines to SAMI mlat/mlon with the SAMI grid itself
          (same method as coast_mag in the notebook)
        - writes web/data/static.json, web/data/events.json and copies plotly.min.js into web/

Step 2. Build one frame per SAMI output (~1 min for the 480 of 2024, ~3 min for the 387 of 2023):
            python3 step2_build_frames.py --event 20230423
        For every SAMI output:
        - loads the Ne cube; SAMI Ne at the Arase R, MLAT, MLT (full grid)
        - averages mlon to at most MAX_NLON bins, closes it across 0/360 deg mlon
        - isosurfaces (skimage marching_cubes) of log10(Ne) above layer K0, vertices mapped
          to x, y, z (Re) at that UT with MLT = (mlon - 72)/15 + UT, noon at +X
        - equatorial plane log10(ne[45, :, :]) and the TEC map (poles filled)
        - writes web/data/frames/yymmdd_HHMM.bin and, at the end, web/data/sami_track.json,
          web/data/manifest.json and web/data/events.json
        Quick test with fewer frames:  python3 step2_build_frames.py --event 20230423 --every 10
        (this rewrites manifest.json with only those frames; run without --every for all)

Step 3. (optional) Contour-mapping and SED-edge overlays (~1 min for the 480 of 2024):
            python3 step3_build_overlays.py --event 20241010
        Only for events in OVERLAYS (config.py); the settings there (MAP_LEVELS, SED_GRAD,
        SED_LON, SED_LAT, H_FOOT, colours, FL_EVERY) are those of the notebook cells. Needs the
        SAMI drive, steps 1 and 2 (it reads manifest.json) and contourpy. Rerun it after step 2.
        Events without overlays show the two boxes greyed out.

Step 4. Open the website. Browsers do not let a page opened as a file read the data files,
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
  add contour mapping / add SED edge (boxes above the 3D view)
                                load and draw that overlay (see above); untick to hide it again
  Arase leg / 2D mapping (top right of the right panel)
                                choose what the right panel shows; drag = zoom, double-click = reset
  camera icon (top right)       save the current view as a PNG
  home icon (top right)         back to the initial view
Frames are loaded only when needed (the next 6 are loaded ahead while playing), so the
page opens quickly even though all frames together are ~270 MB.


Sharing
-------
Copy the whole web/ folder (with data/) to another computer and run step 4 there, or put
web/ on any static web host (GitHub Pages, a university web space): no server code is needed.

GitHub Pages: this folder is the git repository github.com/PengYululu/SAMI_Arase_website.
Pages serves the main branch root; index.html there forwards to web/, so the site is
    https://pengyululu.github.io/SAMI_Arase_website/
To update it after re-running step 1 / 2 / 3 or editing the page:
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
  - After changing config.py, re-run step 1 and step 2 (and step 3), then reload the page (if the browser
    shows old data, force-reload: Cmd-Shift-R).
