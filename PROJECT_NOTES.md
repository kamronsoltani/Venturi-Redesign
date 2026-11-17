# TuriTap handpump Venturi: project notes

Working notes for the CFD study: problem, setup, results so far, and the
next stage. Read this before running anything.

## Environment

- OpenFOAM Foundation **v13** in **WSL2 Ubuntu-22.04** at `/opt/openfoam13`.
  Load with `source /opt/openfoam13/etc/bashrc`. Solver:
  `foamRun -solver incompressibleFluid`.
- Python 3 with numpy, scipy, matplotlib (WSL). gmsh 4.15 Python API
  (`pip install --user gmsh`) for STEP import, booleans and STL export.
  Ubuntu's `python3-gmsh` (4.8) is too old for reliable booleans.
- OpenFOAM 13 aborts on paths containing spaces (`fileName::stripInvalid`).
  `cfd/run_case.sh` stages each case in WSL `~/turitap_cfd/<geom>/<case>`,
  runs checkMesh + foamRun there, and copies logs, postProcessing and the
  final time back into `cfd/<geom>/<case>/`.
- From Windows, call WSL through a script file
  (`wsl -d Ubuntu-22.04 -- bash "<script.sh>"`); `$VARS` in an inline
  `bash -c` string get expanded by the Windows-side shell first. From Git
  Bash set `MSYS_NO_PATHCONV=1`.
- Meshes, time directories and logs are not tracked in git. Regenerate with
  `cfd/make_wedge_mesh.py` + `cfd/setup_case.py` + `cfd/run_case.sh`.

## The problem

- The Venturi hangs vertically, flow downward, outlet open to atmosphere.
  Upstream: rubber coupler on the pump spout, the E-Sig eccentric reducer,
  then a threaded 90 degree elbow into the Venturi inlet.
- On handpumps the current part backs water up until it spills from the
  pump's overflow holes. Average delivery is about 20-30 L/min. The pump body
  smooths the piston pulses, so the steady AVERAGE flow is the design case.
- Allowable water level above the Venturi inlet face: about 8.5 cm (Afridev),
  about 19 cm (India Mark II).
- Chlorine: a float valve holds chlorine level with the throat; a needle
  valve and Teflon tube feed a 1/4-28 port with a 1.67 mm orifice at the
  throat. The port, thread and fittings must NOT change.
- Dosing cutoff must stay at or below 8 L/min. Port suction at 6-8 L/min is
  dominated by the water column between the port and the outlet, so it only
  exists if the outlet runs full.

## Geometry

- `geometry/bore_stations.csv`: bore diameter vs z (mm from the inlet face)
  for the current part, A, B and V1-V4.
- `reference/`: current part STEP and drawing.
- `geometry/*.step`: prototypes A (14 mm throat, 19 mm exit) and B (14 mm
  throat, 26 mm exit, optional cross-vane), each with 3/4" and 1" inlets.
- `geometry/venturi_set_1in/`: V1-V4 (1" inlet, 20 degree cone starting at
  the inlet face at 31.2 mm bore, R14 fillets, 4 degree diffuser):

  | Tag | Throat | Exit | Length | Vane |
  |---|---|---|---|---|
  | V1 | 14 mm | 26 mm | 125 mm | 4 blades |
  | V2 | 16 mm | 26 mm | 125 mm | 4 blades |
  | V3 | 18 mm | 26 mm | 125 mm | 4 blades |
  | V4 | 14 mm | 30 mm | 153 mm | 6 blades |

- `geometry/prototype_C*/`: integrated long-radius elbow + Venturi concept.
- `geometry/upstream/` (for the 3-D stage):
  - `ESig_reducer_3-4in_NPT_male.step`: E-Sig reducer as used today
    (inch-based model, mm units). 44 mm bore spigot into the rubber coupler,
    eccentric flat-top S-curve transition, 21 mm bore straight section,
    3/4" NPT male end. Wall 2.54 mm. The original file name says "1 inch"
    but the end is 3/4".
  - `ESig_reducer_1in_NPT_male.step`: the same part rebuilt for 1": identical
    up to x = 14 mm past the collar, same S-curve shape, 26.64 mm bore,
    31.72 mm OD, 25 mm 1" NPT tapered blank (thread not modelled), same
    overall length 197.4 mm. Script and parameters in `esig_1in/`.
  - `LongRadius_reducer_1in_NPT_female.step`: reducer with an integrated
    large-radius bend ending in a 1" NPT female port (31.2 mm bore).
  - `Elbow_3-4in_NPT_threaded_inch_units.step`: standard 3/4" threaded
    elbow, 20.9 mm bore, about 15 mm bend centreline radius. STEP units are
    INCHES.
  - The 1" elbow is McMaster-Carr 4880K464 (standard-wall PVC, 1" NPT female
    both ends, 26.4 mm bore, 19.8 mm bend centreline radius, inch units):
    `geometry/upstream/Elbow_1in_NPT_McMaster_4880K464.step` (as downloaded
    from McMaster-Carr). This is the 1" elbow to use; no other off-the-shelf
    1" elbow is part of the study.

## Stage 1/3/3b setup (2-D, done)

Single-phase, steady, incompressible RANS on a 2-D axisymmetric wedge of the
bore only. Water at 25 C: rho 997 kg/m3, nu 8.9e-7 m2/s. Gravity ON along
the flow. k-omega SST, wall-resolved (y+ about 1, no wall functions). Inlet
fixed flow rate; outlet atmospheric with backflow-safe velocity; no-slip
walls. SIMPLEC, 4000 iterations. All 54 runs converged; every mesh passed
checkMesh (max non-orthogonality 20 deg). Mesh GCI on inlet level about 2%.

Reported per run: inlet total pressure as water level above the inlet face
(the overflow number); inlet static pressure (flag if negative); port wall
pressure; diffuser recovery vs ideal; sign of axial wall shear along the
diffuser; outlet KE flux as head and KE factor alpha; mass imbalance and
residuals. Results: `cfd/results/*.csv`, report in `report/`, viewer in
`docs/` (https://kamronsoltani.github.io/Venturi-Redesign/).

Key results, level above the Venturi inlet (cm), Venturi alone:

| Design | 20 | 25 | 30 | 40 | 50 L/min |
|---|---|---|---|---|---|
| current | 6.2 | 15.5 | 26.7 | | |
| B_1in | -3.6 | 0.4 | 5.1 | 16.7 | 31.1 |
| V1 | -5.6 | -2.2 | 1.8 | 11.7 | 23.9 |
| V2 | -8.1 | -6.0 | -3.4 | 3.0 | 10.9 |
| V3 | -9.2 | -7.5 | -5.4 | -0.3 | 6.2 |
| V4 | -9.2 | -6.2 | -2.7 | 5.9 | 16.4 |

Port suction at 30 L/min (kPa gauge): current -6.7, V1 -5.1, V4 -5.6,
V2 -3.4, V3 -2.4. The 4 degree diffusers stay attached; B's 5 degree
diffuser reverses over 30-59% of its wall. CFD level is 17-22% below the
1-D model; port pressure agrees within 3%.

Bench/field (printed B_3-4in with vane): dosing starts at 8 L/min on steady
piped flow; on a handpump it overflows later than the current part but still
overflows. Reducer + elbow with no Venturi nearly overflow on their own, so
the upstream parts are now the limit.

## Stage 5b (next): 3-D upstream + Venturi benchmark

Goal: rank complete assemblies from the pump spout to the outlet, so the
upstream parts are included and every Venturi is compared behind the same
upstream train.

| Group | Upstream | Venturis | Assemblies |
|---|---|---|---|
| 1 | 3/4" E-Sig + 3/4" threaded elbow | current, A_3-4in, B_3-4in (with vane) | 3 |
| 2 | 1" E-Sig + McMaster 4880K464 1" elbow | V1, V2, V3, V4 | 4 |
| 3 | 1" long-radius reducer | V1, V2, V3, V4 | 4 |

- Flows: 20, 25, 30, 40, 50 L/min (55 runs). Run 30 L/min for all 11 first
  to get the ranking, then the other flows (initialise from the nearest
  finished flow).
- Mesh study at 30 L/min on one assembly (three densities) before the sweep.
- Fluid volume: from the actual STEP solids (internal water path from the
  E-Sig spigot entry to the Venturi outlet), including the steps and gaps at
  each threaded joint and the cross-vanes. Assumptions unless told otherwise:
  male ends engage to standard hand-tight plus wrench make-up; Venturis
  screw directly into the female elbow/reducer with no nipple. State every
  assumption in the results.
- Add a straight inlet extension (about 5 diameters of the 44 mm spigot bore)
  so the inflow is developed, and keep the outlet open to atmosphere.
- Half domain with a symmetry plane if the assembly is symmetric (bend plane);
  full domain otherwise.
- Mesh: snappyHexMesh from STL. WSL memory (about 6 GB by default) limits the
  cell count, so wall functions are expected; report y+. Run V1 alone in 3-D
  and compare with its 2-D result to show how much the method shifts the
  numbers.
- Same physics as Stage 1 otherwise.
- Report per run: the Stage 1 quantities with the level referenced to the
  spigot entry, plus the loss of each part (reducer, elbow, Venturi), the
  velocity non-uniformity at the Venturi inlet and throat, the wall pressure
  difference around the throat (port side vs opposite), where the diffuser
  flow reverses, and the outlet KE factor. Write
  `cfd/results/assembly_3d_stage5b.csv`.
- Rough cost: about 1 day of setup, then about 30-60 min per run one at a
  time (about 40 h for 55 runs).

### Stage 5b status (2026-11-17)

First assembly done: 3/4" E-Sig (flat side down) + 3/4" threaded elbow +
current Venturi at 30 L/min. Scripts in `cfd/assembly3d/`, case in
`cfd/asm34_current/Q30`, results in `cfd/results/assembly_3d_stage5b.csv`,
3-D flow viewer in `docs/3d/`. 630k cells, wall functions (y+ mean 3-11).
Level above the inlet centre 32.1 cm (24.3 cm above the Venturi inlet, vs
26.7 cm in 2-D); losses: E-Sig 2.3 cm, elbow 1.4 cm, Venturi 14.8 cm, outlet
jet 27.9 cm. Inlet level steady to 0.17 cm over the last 500 iterations;
residuals plateaued at 2e-3 (p).

Next: the baseline sweep of the same assembly at 6, 8, 10, 20, 40, 50 and
60 L/min plus a coarse/fine mesh check at 30 L/min, run unattended with
`cfd/assembly3d/run_sweep.sh` (mesh reused across flows, each flow started
from the nearest finished one, 2000 iterations). How to run it on another
computer: `LAB_SETUP.md`.

## Later stages (not started)

- Stage 2: compare with bench measurements (inlet water level vs flow; port
  suction vs flow) when available.
- Stage 4: two-phase VoF at 4-10 L/min to see whether the outlet stays full
  (cases for the current part exist in `cfd/current_vof/`, never run).
- Stage 6: stroke-by-stroke water-level model once stroke volume and pumping
  rate are measured.

## Rules

- Never print full solver logs; read the last ~50 lines or postProcessing.
- Run checkMesh before every solve. Stop if it fails or max
  non-orthogonality exceeds 65 degrees.
- Do not tune settings to match the 1-D model. Report disagreements.
- One case folder per geometry and flow rate under `cfd/`.
- Summarise results in a CSV under `cfd/results/` and state which runs
  converged and which did not.
