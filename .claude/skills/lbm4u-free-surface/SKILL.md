---
name: lbm4u-free-surface
description: Set up, run, check and debug free-surface flow simulations (problemID 1) with the cLBMnet / LBM4U lattice Boltzmann solver in this repository - mold filling, dam break, sloshing, jets, droplets, inlets and outlets, surface tension and contact angle, octree refinement. Use it whenever the user wants to prepare or run a simulation case, choose the mesh size, uMax, time step or relaxation time, write parameters.txt, blocks.txt, blockBoundary or blockInit files, interpret result/output.pvd or residuum.txt, or fix a diverging or mass-losing free-surface run - also when they say it in Czech ("volná hladina", "plnění formy", "výpočet s hladinou", "připrav úlohu", "spusť simulaci"). The non-Newtonian and particle skills build on this one.
---

# Free-surface flow with cLBMnet (LBM4U)

The solver is a D3Q19 lattice Boltzmann code with a VOF / mass-tracking free surface (Körner, Thürey). Only the
liquid is simulated; the gas is empty space whose pressure enters through the interface cells. A case is a set of
text files; LBM Studio (`gui/`) writes them for CAD geometries, and `scripts/make_case.py` writes them for
geometries that can be described by boxes and expressions. The solver reads the files, runs, and writes VTK.

Read `references/parameters.md` when you need a parameter that is not in the tables below. Numerical details
(collision operators, LES, refinement, fixed bugs) are in `CHANGES.md`; the free-surface code is
`BlockFreeSurface.cpp` and the time step order is in `Mesh::timeStep` / `Mesh::advanceFreeSurfaceLevel`.

## Workflow

1. **Collect the physics.** Domain size, liquid region at t = 0, inlets or outlets, density ρ, kinematic viscosity
   ν, surface tension σ, gravity, the expected largest velocity U and the simulated time T. Ask only for what is
   missing and cannot be estimated (U can usually be estimated as `sqrt(2 g height)` or the inlet velocity).
2. **Choose the resolution and the velocity scale** (next section). Run `make_case.py --check` to see τ, Δt, the
   lattice gravity and the warnings before writing anything large.
3. **Write the case** - with LBM Studio for a STEP/STL geometry, otherwise with `make_case.py` (JSON description,
   see the docstring at the top of the script). Building by hand is fine for small changes, the formats are below.
4. **Build the solver if needed**: `build.bat` (MSVC, exe in `build-msvc/Release/`) or `build.bat mingw`.
5. **Run**: `cLBMnet --paths <geometryDir> <simulationDir>` (or `cLBMnet <geometry> <simulation>` for
   `simulations/<geometry>/<simulation>/`). Long runs go in the background.
6. **Watch** the console or `result/residuum.txt` (`iter uMax mass`). The mass of a closed case must stay constant
   (to ~1e-12 relative); uMax must stay well below `Cu` (lattice velocity < 0.3, the solver warns).
7. **Look at the results**: `result/output.pvd` (volume fields) and `result/outputFS.pvd` (free surface triangles)
   in ParaView; `fields_<iter>.bin` is read by LBM Studio.

## Units: how the physical case maps to the lattice

Everything in the input files is SI. The solver derives the lattice from three parameters:

| Quantity | Formula | Notes |
|---|---|---|
| cell size (level 0) | `h = H / ni` | `H` is the edge of a level-0 block in metres |
| velocity scale | `Cu = uMax / 0.1` | `uMax` maps to the lattice velocity 0.1 |
| time step | `Δt = h / Cu` | acoustic scaling; level L has `h/2^L`, `Δt/2^L` |
| relaxation time | `τ = 1/2 + 3 ν Δt / h²` | level L: `τ_L = 1/2 + 2^L (τ_0 − 1/2)` |
| lattice gravity | `g Δt² / h` | keep below ~1e-3 |
| physical time | `T = maxIter · Δt` | `maxIter`, `saveRate` count level-0 steps |

How to choose them:

* `uMax` must be at least the largest expected velocity (rather 1–2× it). The run is invalid when the lattice
  velocity `u/Cu` exceeds ~0.3 (compressibility); larger `uMax` costs proportionally more steps.
* The hydrostatic density change `3 g L_z / Cu²` must stay below a few percent - it is the compressibility error of
  a column of height `L_z`. Tall domains with slow flow need a larger `uMax` than the velocity alone suggests.
* Low-viscosity liquids (water: ν = 1e-6) give τ ≈ 0.5000x. That is stable only with the Smagorinsky model and a
  bounded odd relaxation time: `collisionModel = 2` (TRT) or `1` (MRT), `smagorinsky = 0.1–0.17`,
  `tauMinusMax = 1`. LBM Studio writes these values; keep them for water, metals and oils with τ < 0.51.
* τ above ~2 makes walls inaccurate and the run slow - lower `uMax` or refine.
* Resolve the smallest feature (gap, jet, film) with at least 4–6 cells; use octree refinement instead of a
  uniformly fine mesh.

Worked example (dam break, water column 10 × 20 mm, h = 1 mm): U ≈ sqrt(2·9.81·0.02) = 0.63 m/s → `uMax = 1`,
`Cu = 10 m/s`, `Δt = 1e-4 s`, τ = 0.5003 (→ LES), lattice gravity 1e-4, 3 g L_z/Cu² = 0.01. 0.5 s = 5000 steps.

## Case layout and file formats

```
<geometryDir>/mesh/blocks.txt                     # one block per line: "bi bj bk" or "bi bj bk level"
<geometryDir>/mesh/boundary/blockBoundary<tag>.txt  # cell types, one per line (missing file = all fluid)
<simulationDir>/parameters.txt                    # "name = value" lines
<simulationDir>/init/blockInit<tag>.txt           # "rho u v w flag" per cell (missing file = empty block)
<simulationDir>/result/ , out/                    # written by the solver
```

* A block has `ni × nj × nk` cells. Level-0 block files have the tag `nBi·nBj·bk + nBi·bj + bi`; refined blocks
  `L<level>_<bi>_<bj>_<bk>` with coordinates on their level. Cells are listed with **k fastest, then j, then i**.
* Cell types: 0 fluid, 1 wall (bounce-back), 2 velocity inlet (`uInlet, vInlet, wInlet`), 3 pressure outlet
  (`pOutlet` [Pa] relative), 4 moving wall (`uwall, vwall, wwall`). The border of the block grid needs no wall
  cells: a missing neighbour block acts as a bounce-back wall at the block face. Periodicity: `xPeriodicity` etc.
* `blockInit`: `rho` is the lattice density (1 at rest), `u v w` in m/s, flag 1 = liquid, 2 = gas. The interface
  layer is created by the solver between liquid and gas. `initializeWithHydroPress = 1` replaces ρ by the
  hydrostatic profile - use it for liquid initially at rest (avoids a pressure wave at t = 0).
* Octree refinement: 2:1 balance across faces, edges and corners, even `ni, nj, nk`, not D2Q9; available for
  problemID 0–4. `make_case.py` balances the tree automatically.

## Main parameters

| Parameter | Meaning |
|---|---|
| `problemID` | 0 single phase, **1 free surface**, 2 Lagrangian particles, 3 Euler FVM, 4 Euler LBM, 21 fibers |
| `model` | 319 = D3Q19 (default choice), 315 = D3Q15, 29 = D2Q9 (`nk = 1`) |
| `collisionModel` | 0 BGK, 1 MRT, 2 TRT |
| `ni nj nk`, `H` | cells per block, block edge [m] |
| `uMax`, `density`, `vKin` | velocity scale [m/s], ρ [kg/m³], ν [m²/s] |
| `gx gy gz` | gravity [m/s²]; `fluidGravity = 0` applies it to particles only |
| `maxIter`, `saveRate`, `threads` | level-0 steps, output interval, OpenMP threads |
| `smagorinsky`, `tauMinusMax`, `tauMin`, `tauMax` | LES constant, odd-τ limit, limits of the local τ |
| `surfaceTension` [N/m], `curvatureModel` | σ; curvature 0 = smoothed normals (default), 1 = height function |
| `wetAngleEq/Adv/Rec` | contact angles [deg] |
| `uInlet vInlet wInlet`, `pOutlet` | inlet velocity [m/s], outlet pressure [Pa] |

## Running, stopping, restarting

* Create a file `STOP` in the simulation directory to end the run cleanly (results are written).
* Exit code 3 = divergence (NaN). The final state is written to `out/blockInit*.txt`; copy `out/*` to `init/` to
  continue from it.
* MLUPS in the log counts liquid cells; "active blocks" counts all cells of non-empty blocks.

## When the run goes wrong

| Symptom | Likely cause and fix |
|---|---|
| NaN, exit code 3 early | τ ≈ 0.5 without LES; lattice velocity too high (raise `uMax`); inlet velocity above `uMax` |
| "lattice velocity ... is high" | raise `uMax` (smaller Δt) |
| liquid oscillates / bounces at rest | hydrostatic compressibility: raise `uMax`, set `initializeWithHydroPress = 1` |
| mass drifts in a closed box | open boundary (inlet/outlet) present, or a bug - report with the case; mass is otherwise exact, also across refinement |
| "initial condition ... not found, block is empty" | missing `blockInit<tag>.txt`: intended for empty blocks, otherwise wrong tag |
| "the mesh is not 2:1 balanced" | refine the coarse neighbour too (or let `make_case.py` do it) |
| spurious currents at the surface | surface tension with coarse cells; try `curvatureModel = 1`, refine near the surface |

Unit tests of the free surface: `build-msvc/Release/lbm_tests.exe free_surface` (and `refinement`). Run them after
changing solver code.
