# cLBMnet v2 – change log

The starting point is the first commit (copy of the branch `rotace_master_new` @ 93dd6e7).

## Building and tests

```
build.bat            # Visual Studio 2022 (MSVC, OpenMP /openmp:llvm)
build.bat mingw      # MSYS2 g++ + Ninja, statically linked exe
```

CMake builds `cLBMnet.exe` (the solver) and `lbm_tests.exe` (33 unit tests, run them in the build folder, optionally with a name filter).

Graphical interface: `gui/build_gui.bat` → `gui/bin/LBMStudio.exe`, described in [gui/README.md](gui/README.md).

The option `-DLBM_FP32=ON` stores the distribution functions in single precision (half the memory and memory traffic).

Running the solver:

```
cLBMnet                                # as before, reads args.txt
cLBMnet <geometry> <simulation>        # simulations/<geometry>/<simulation>/
cLBMnet --paths <geometryDir> <simulationDir>
```

A run can be stopped early by creating a file named `STOP` in the simulation directory (the results are saved).
Exit code 3 means divergence (NaN).

## Fixed bugs

| Area | Bug | Fix |
|---|---|---|
| Mesh step 8 | `makro`, `collision` and `newFsGeometry` in one parallel loop, neighboring blocks read half-updated `Rho/U/flagMap` | three separate passes (macroscopic values → collision + preparation of new cells → applying the changes) |
| Mass redistribution | concurrent `mass +=` into neighboring blocks, `massErr` was never returned | `omp atomic`, `massErr` is distributed among the interface cells of the block |
| Step 43 | race on `nParticles`, `maxParticlesInCell` | serial reduction per block |
| FVM particles | the flux across a block boundary was written into the wrong block, `oppLW` had wrong opposite directions, the concentration changed while neighbors were reading it | conservative flux `(F(A→B) − F(B→A))/2`, separate `computeFlux`/`applyFlux` |
| `BlockBasic` | `isBlockEmpty` not initialized; gravity 3× weaker and a `−rho0` term | initialization; Guo forcing |
| `findArtefactsFilled` | `boundary[ind]` instead of `boundary[indShift]` | fixed |
| Artefacts | simultaneous conversion of neighboring interface cells to liquid and gas created a liquid–gas contact and mass loss | a cell next to a "filled" artefact is not emptied, safety net `fixFluidGasContacts` |
| `initHydrostaticPress` | shadowed variable `lo`, reference from the top of the domain | computed from the highest liquid cell |
| Restart | `writeResultsToOutFile` used the sequential block index | correct mapping and file ids |
| `coarseLevel` | inconsistent array sizes | replaced by static refinement (see below) |
| Non-Newtonian fluids | worked only with MRT; shear rate from differences of neighbors (race) | any collision model, local shear rate from the non-equilibrium moments |
| Boundary files | `line.compare("1") == 1` | reading the integer cell type, error checks |
| File reading | `catch (int)` did not catch exceptions | `std::exception`, readable messages |
| Marching cubes | `realloc` + `delete[]`, the "surface" condition was always true | `std::vector`, min/max test of the iso value |
| `saveRate` | default `1e10` in an `int` | safe parameter reading (also the notation `1e5`) |
| `particleBarrierForce` | `omp parallel` without `for` | fixed |
| Euler-LBM particles | the new particle velocity used already overwritten components | the old values are used |
| Other | memory leaks, `delete` vs `delete[]`, `getParDoubleArray`, D2Q9 crash (`nk = 1`), wrong `tau` print, overflow of output names, `windows.h` | fixed |

## Numerics

* BGK, TRT and MRT collision operators in a common class `CollisionOperator` (parameter `collisionModel`, the old `MRT` still works).
  * TRT with the magic parameter `trtMagic` (default 3/16 – the bounce-back wall lies exactly halfway, Poiseuille flow is exact to 1e-13).
  * MRT: the odd moments relax with `s_q = 8(2 − s_ν)/(8 − s_ν)` (`mrtMagicSq = 1`), the wall position does not depend on the viscosity.
  * The MRT matrices are inverted numerically, D2Q9 has its own MRT matrix.
* `BlockBasic`: Guo forcing (the velocity includes F/2).
* Free surface: output velocity `u + F/2` (EDM), particle forces are divided by the density.
* Mass exchange between interface cells weighted by the cell types (Thürey 2007).
* Excess mass redistribution weighted by the normal.
* Curvature from a height function (`curvatureModel = 1`); the default stays the divergence of the smoothed normals, which gives smaller parasitic currents in dynamic flows.
* Boundary conditions: inlet (`uInlet, vInlet, wInlet`, cell type 2), pressure outlet (`pOutlet`, type 3, anti-bounce-back), moving wall (`uwall, vwall, wwall`, type 4, domain lid `movingLidTop`).
* Optional limits `tauMin`, `tauMax` for the variable relaxation.
* Smagorinsky LES model (`smagorinsky` = C, typically 0.1–0.17) for single phase and free surface computations, explicit relation of Hou et al. (1996).
* `tauMinusMax`: upper limit of the relaxation time of the odd moments (TRT τ⁻, MRT s_q). For τ → 0.5 (water) the magic relation gives τ⁻ ≈ 10³ and the computation diverges; with `tauMinusMax = 1` and LES, filling with water is stable. The GUI always writes this value.
* The solver warns when the lattice velocity exceeds 0.3 (high Mach number).

## Performance

* Removed the unused array `fcOld` (1/3 of the population memory).
* Neighbor tables only for the boundary cells of a block (instead of 2·nf ints per cell).
* `BlockBasic`: fused computation of moments and collision.
* Optional FP32 populations + DDF shifting (`ddfShift = 1`).
* Distance function by Dijkstra's algorithm (same result, O(N log N)).
* Output in the VTK XML format: `output_*.vtu`, free surface `outputFS_*.vtp`, `boundary.vtu`, collections `output.pvd`/`outputFS.pvd` (ParaView reads only XML files in a `.pvd`; for the older legacy `.vtk` files in the collection it reported "Could not determine the data type"). Binary data after the header (`vtkBinary = 1`, otherwise ASCII), compact `fields_<iter>.bin` for the GUI (`writeFields = 1`).
* Result compression (`compressResults = 1`, default): the VTK arrays are compressed with zlib in 1 MB blocks (`vtkZLibDataCompressor`, read directly by ParaView), `fields_<iter>.bin` has version 3 with the block data in one deflate stream. The compressor is our own (`Deflate.cpp`, no zlib dependency), blocks are compressed in parallel. Only the points used by cells are written to the `.vtu`.
* MLUPS from `omp_get_wtime` averaged over the save interval.

## Grid refinement (octree)

Static refinement driven by the geometry: the mesh is created once during generation and does not change during the computation (it does not follow the free surface).

* A block of level `L` still has `ni × nj × nk` cells of size `H / ni / 2^L`. Neighboring blocks (also across an edge and a corner) differ by at most one level (2:1), the solver rejects an unbalanced mesh.
* `mesh/blocks.txt`: a line `bi bj bk` (level 0, unchanged) or `bi bj bk L` with the block coordinates on its level. Files of level 0 blocks keep the number `nBi·nBj·bk + nBi·bj + bi`, finer blocks have the tag `L<L>_<bi>_<bj>_<bk>` (`blockBoundaryL2_5_3_1.txt`, `blockInitL2_5_3_1.txt`).
* Acoustic scaling: both `Ct` and `Ch` are halved on each level, `Cu` and the pressure are shared. A finer level makes two steps per step of the coarser one (sub-cycling), `τ_{L+1} = ½ + 2(τ_L − ½)` – fine levels are more stable. Iterations and output times are steps of level 0.
* Level interface: the volumetric method of Rohde et al. (2006), as in waLBerla. Coarse post-collision populations are distributed homogeneously into the ghost layer of the fine level, which streams them during two fine steps. A coarse cell receives the average of the populations that streamed into its 8 fine cells from the fine region.
* Ghost blocks: blocks of the same class (`Block::ghost`) surround the blocks of each level and are filled with the state of another level. That is why the existing neighbor code (`nbBlock`, `nbLocal`, `local2globalShift`) works unchanged. Streaming and balances run only in cells within distance 2 of the computed cells.
* Mass is conserved exactly (machine precision). The difference between what the coarse level hands over and what the fine level takes (and vice versa) is returned to the coarse cell (`massCorr`). The difference arises mainly where the level interface crosses a wall.
* Rescaling of the non-equilibrium part of the populations (Dupuis–Chopard) was tried and removed. In the volumetric scheme it breaks momentum conservation per population and gave an error of +21 % / −11 % in the Poiseuille test.
* Free surface across the interface: the mass exchange between interface cells of different levels is computed by the fine side and the coarse cell receives the opposite value (`massAcc`). Cell conversions (gas → interface, liquid → interface) on ghost cells are done immediately and sent to the level that owns the cell. The same applies to mass redistributed into ghost cells. A link between a liquid cell and a gas ghost cell (the conversion arrives one step later) does not transfer mass. Ghost blocks compute local quantities (interface types, levels, VOF smoothing, normals, curvature). Empty blocks are activated across ghost blocks too.
* Limitations: only `problemID` 0 and 1, not D2Q9, even `ni, nj, nk`. Where the level interface crosses a wall perpendicular to the flow, a Poiseuille channel 32 coarse cells high has a flow rate error of 1–3.5 % (a parallel interface ~1e-3). The error is in the corners of the interface and the wall, not in the mass.
* Output: VTK contains the blocks of all levels in the cell coordinates of level 0, `fields_<iter>.bin` of a refined mesh has version 2 with the block level, MLUPS counts the `2^L` steps of the fine levels.
* Tests: `refinement_*` (fluid at rest, rejection of an unbalanced mesh, mass conservation, Poiseuille flow with the interface parallel and perpendicular to the flow, free surface – dam break across the interface and a free surface at rest).

## Not done

* AA-pattern (a single population array) and SoA layout.
* Refinement for particles and fibers (`problemID` 2, 21, 3, 4).

## Format of `fields_<iter>.bin`

```
char[4] "LBMF", int32 version (1, refined mesh 2, compressed 3), int32 iteration, double time, int32 problemID,
int32 nBlock, int32 ni, nj, nk, double Ch (level 0), int32 nFields, nFields × char[16] name
for each block: int32 bi, bj, bk, [version 2, 3: int32 level]; for each field float32[ni*nj*nk] (index i + ni*j + ni*nj*k)
version 3: everything after the header (from the first block on) is one raw deflate stream (RFC 1951, no zlib header)
fields: ux, uy, uz [m/s], p [Pa], type (cell type), with a free surface also vof, flag
```
