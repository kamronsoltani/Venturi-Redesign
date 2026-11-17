# cfd/

OpenFOAM 13 cases, one folder per geometry and flow rate (`<geom>/Q<L/min>`).
Each case keeps its dictionaries, `case_info.txt`, `mesh_info.json`,
`results.json` and the `viz/` plots; meshes, time directories and logs are
not tracked (regenerate with the scripts below).

| Script | Purpose |
|---|---|
| `make_wedge_mesh.py` | 2-D axisymmetric wedge blockMesh from `geometry/bore_stations.csv` |
| `setup_case.py` | Writes a Stage 1 case (k-omega SST, wall-resolved, gravity on) |
| `run_case.sh` | Stages a case in a space-free WSL folder, runs checkMesh + foamRun, copies results back |
| `postprocess.py` | Level, port pressure, recovery, wall shear, KE flux, mass balance -> `results.json` |
| `make_results_csv.py` | Collects `results.json` files into `results/*.csv` |
| `viz_axisym.py` | Per-case plots and the 3-D viewer data |
| `build_gui_data.py` | Builds `results/gui/` (copied to `../docs/` for the web viewer) |
