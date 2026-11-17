# Venturi Redesign: TuriTap handpump chlorinator

CFD study of the TuriTap in-line Venturi chlorinator on India Mark II and
Afridev handpumps. The current part backs water up until the pump overflows
at the average delivery of 20-30 L/min; the goal is a Venturi that needs
less water level above its inlet while keeping chlorine dosing down to
8 L/min.

**Interactive viewer:** https://kamronsoltani.github.io/Venturi-Redesign/
(all 54 runs: 3-D views, probes, particle tracing, wall and line plots,
flow sweeps, convergence).

## Contents

| Folder | |
|---|---|
| `docs/` | Web CFD viewer (GitHub Pages) |
| `report/` | LaTeX report (`main.tex`, figures, tables, Overleaf zip) |
| `cfd/` | OpenFOAM cases, scripts, and `results/*.csv` |
| `geometry/` | Bore stations and STEP/STL files for prototypes A, B, C and V1-V4 |
| `reference/` | Current part STEP and drawing |
| `analysis/` | 1-D Bernoulli/loss model (MATLAB) and prototype drawings |

## Method

Steady, incompressible, single-phase RANS on a 2-D axisymmetric wedge of the
bore (OpenFOAM 13, `incompressibleFluid`), k-omega SST with a wall-resolved
mesh (y+ ~ 1), gravity along the flow, fixed flow-rate inlet and atmospheric
outlet. Water at 25 C. Mesh study (three densities) at 25 L/min: GCI on the
inlet level about 2%. All runs converged.

## Key result: water level above the Venturi inlet (cm)

| Design | 20 L/min | 25 | 30 | 40 | 50 |
|---|---|---|---|---|---|
| Current part | 6.2 | 15.5 | 26.7 | | |
| B (14 mm throat, 5 deg diffuser) | -3.6 | 0.4 | 5.1 | 16.7 | 31.1 |
| V1 (14 mm throat, 4 deg, 26 mm exit) | -5.6 | -2.2 | 1.8 | 11.7 | 23.9 |
| V2 (16 mm throat) | -8.1 | -6.0 | -3.4 | 3.0 | 10.9 |
| V3 (18 mm throat) | -9.2 | -7.5 | -5.4 | -0.3 | 6.2 |
| V4 (14 mm throat, 30 mm exit) | -9.2 | -6.2 | -2.7 | 5.9 | 16.4 |

The 4 degree diffusers stay attached; B's 5 degree diffuser reverses over
30-59% of its wall. Full results, derivations and the comparison with the
1-D model are in the report.

## Reproducing

Requires OpenFOAM 13 (tested in WSL Ubuntu 22.04) and Python 3 with numpy,
scipy and matplotlib. See `cfd/README.md`.
