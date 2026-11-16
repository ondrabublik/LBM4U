# cLBMnet v2 – change log

The starting point is the first commit (copy of the branch `rotace_master_new` @ 93dd6e7).

## Building and tests

```
build.bat            # Visual Studio 2022 (MSVC, OpenMP /openmp:llvm)
build.bat mingw      # MSYS2 g++ + Ninja, statically linked exe
```

CMake builds `cLBMnet.exe` (the solver) and `lbm_tests.exe` (59 unit tests, run them in the build folder, optionally with a name filter).

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
* Limitations: `problemID` 0, 1, 2, 3 and 4 (not the fibers 21), not D2Q9, even `ni, nj, nk`. Where the level interface crosses a wall perpendicular to the flow, a Poiseuille channel 32 coarse cells high has a flow rate error of 1–3.5 % (a parallel interface ~1e-3). The error is in the corners of the interface and the wall, not in the mass.
* Output: VTK contains the blocks of all levels in the cell coordinates of level 0, `fields_<iter>.bin` of a refined mesh has version 2 with the block level, MLUPS counts the `2^L` steps of the fine levels.
* Tests: `refinement_*` (fluid at rest, rejection of an unbalanced mesh, mass conservation, Poiseuille flow with the interface parallel and perpendicular to the flow, free surface – dam break across the interface and a free surface at rest).

## Particles (`problemID` 2, 3, 4)

Review of the three particle models on a simple geometry (closed or periodic box of 8³ blocks). The tests are `lagrange_*` and `euler_*` in `tests/test_particles.cpp`.

| Model | `problemID` | Input | Parameters |
|---|---|---|---|
| Lagrangian particles | 2 | `init/particles.txt`: number of particles, then lines `x y z diameter density` (SI units) | `nIntIter` (sub-steps of the particles per LBM step), `maxParticlesInCell`, `couplingWayNumber`, `particleRestitution` (0–1, default 1), `cLES` |
| Eulerian concentration, FVM | 3 | `init/concentration<block>.txt`: number of particles in each cell | `particleDiameter`, `particleDensity` |
| Eulerian concentration, LBM | 4 | the same | `particleDiameter`, `particleDensity`, `particleDiffusivity` [m²/s], `overpackDiffusivity` |

All models use the Schiller–Naumann drag with an implicit velocity update, gravity with buoyancy `(m − V ρ) g` and return the reaction of the drag to the fluid (two-way coupling).

**Gravity.** `problemID` 2 and 4 can rotate `(gx, gy)` around the z axis with `gravityRotationOmega` [rad/s] (rotating mould). The original code always rotated with 0.89 rad/s and set the z component to zero, so `gz` was ignored. Now the default is 0 (no rotation), `gz` is kept and the solver prints a line when the rotation is active. Older cases of the rotating mould have to set `gravityRotationOmega = 0.89`.

Fixed (value before → after in the tests):

| Where | Bug | Result |
|---|---|---|
| Lagrange, output | `writeResultsToVTK` filled the cell lists of particles again, the next step counted every particle twice (double repulsive force and concentration) – the result depended on `saveRate` | identical with and without output |
| Lagrange, walls | a particle closer to a wall than its radius was frozen for ever (also on vertical walls and in a tangential flow); the reflection also turned particles moving away back to the wall | the particle is moved to the distance d/2 along the wall normal, the normal velocity is reflected (`particleRestitution`), the tangential motion is kept |
| Lagrange, force spreading | the weights of the cells in walls and in the gas were dropped (21 % of the kernel for a particle 0.6 cell above a wall), `3.1415` instead of π (momentum error 8.5e-5) | weights normalised by `Σ w·vof`, momentum of fluid + particle conserved to 2e-14 |
| Lagrange, gas | the force of the last step in the liquid (uninitialised memory at the start) was applied to the fluid for ever, the collision force kept growing | no force outside of the liquid |
| Lagrange, domain | no range check of the block index: a particle outside of the domain was mapped to another block; particles in empty blocks overflowed the cell lists | `LagrangeParticles::locate`, periodic directions wrap the particle |
| Lagrange, sub-steps | only the force of the last sub-step reached the fluid | impulse summed over the sub-steps |
| Lagrange | missing `particles.txt` left the number of particles uninitialised; arrays freed with `delete` | error message, `std::vector` |
| all models | `computeGrad*` leaked `new double[nf]` in every block and step | no allocation |
| FVM, limiter | the gradient ratio was taken from the wrong side for a negative face velocity: concentration −4.8 % … +8.8 % of c0 in the settling test | flux evaluated from the upwind cell: 0 … 1.000 c0, conservation 1e-14 |
| FVM, force | `smoothForce` wrote into the neighbouring blocks from parallel threads (results differed by 1e-5 between 1 and 8 threads) | gather form, identical results |
| FVM, LBM | particle velocity only where `c > 1e-3` (`1e-4`), zero or fluid velocity elsewhere – dilute suspensions did not settle and fronts were slowed down | the velocity field of a single particle is computed in all liquid cells |
| FVM, LBM | force spread to walls and gas was lost | conservative (`normalizeForce` + `smoothForce`), sum of the force 1.000 |
| FVM, LBM | particles of a cell that lost the liquid stayed in the gas: −12 % (FVM) and −18 % (LBM) of the particles in a dam break | they are passed to the liquid cells around (`releaseDryCells`, `collectReleased`): 1e-15 and 1e-6 |
| LBM | diffusivity 1.4 times the prescribed one (`taus = ½ + 3D`, but the equilibrium has `c_s² = 1/(3θ)`) | `taus = ½ + 3θD`, measured 0.02003 for 0.02 |
| LBM | gradient of the particle velocity (convective term) 6 times smaller than it should be | exact for a linear field |
| LBM | D3Q15 lost all particles (rest weight `1/3` hard coded, negative rest population for θ = 1/1.4) | `w[0]`, θ limited from below; D3Q19, D3Q15 and D2Q9 conserve to 1e-9 |
| LBM | the fluid velocity seen by the particles was zero in the interface cells | defined in all liquid cells |
| LBM | particles in the cells next to a wall moved with the fluid, so they did not settle along vertical walls (55 % of the settling distance in a column 8 cells wide) | only the velocity towards the wall is limited, 99 % of the analytic distance |

Results of the settling tests (sphere, Schiller–Naumann terminal velocity): slip velocity of the Lagrangian particle +0.02 %, settling distance of the concentration −1.0 % (FVM) and −1.0 % (LBM).

Known properties and limits:

* Lagrange: the particle is a point with a drag law. Its own reaction force drags the fluid along (a particle of 0.8 cell settled 19 % faster than the terminal velocity in the test), so the particles should be clearly smaller than a cell.
* Lagrange: the repulsive force between particles `10 / Ca · sqrt(1.5 d − L)` does not depend on the mass of the particles, it is not a contact model; in a periodic box particles do not interact across the boundary. Inlets and outlets are walls for the particles.
* Euler: the concentration is the number of particles per cell, also in the interface cells (it is not related to the fill level). FVM has no packing limit (`PP = 0`), LBM limits it with `overpackDiffusivity` (dimensionless, in lattice units) and adds an artificial diffusion `|c'' |/1000` that depends on the units of the concentration.
* Euler LBM: `G·epsF` together with the reaction of the drag counts the buoyancy of the particles twice for a large volume fraction; `particleInterForces` and the fields `dp`, `dc`, `gVOFp` are computed or read but not used.
* The force of the particles enters the fluid as `F·vof/ρ` (interface cells get less).
* Lagrangian fibers (`problemID` 21) and the files `LagrangeParticlesVar2/3` (not built) were not reviewed.

### Particles on a refined mesh

* Lagrangian particles (`problemID` 2): one particle step per step of the level 0 (before the finer levels), the particle is coupled with the finest computed block holding it (`LagrangeParticles::locate`), interpolation and force spreading use the cells of that level only (ghost cells are skipped, the kernel is renormalised). The impulse of a particle step acts in all 2^L steps of the level L as the force `J · 4^L` per step (momentum unit 8^-L, 2^L steps), the forces are kept until the next `applyForce` instead of being cleared in the collision. The wall distance of the ghost cells is taken from their level (`Mesh::ghostWallDistance`), the wall distance of a particle is in the cells of the level 0. Particle positions are read with the cell size of the level 0 (they were scaled with the first block, which can be a finer one).
* Eulerian particles (`problemID` 3, 4): the particle step runs on every level between its propagation and collision (`Mesh::eulerLevelStep`). The concentration is the number of particles in the cell of the level, ghost cells get 1/8 of the coarser cell or the sum of the 8 finer cells. FVM: the flux through a face between the levels is computed by the finer cells, the coarser cell receives the same number of particles (`exchangeCorr`, `Mesh::forwardFluxFVM`). LBM: the scalar populations are exchanged like the fluid populations (explosion with 1/8, coalescence as the sum of the 8 finer populations) and balanced as in `computeLost`/`accountInterface` (`Mesh::scalarLost`, `scalarCoalesce`, `scalarAccount`). The number of particles is conserved to 1e-14 (FVM) and 1e-13 (LBM) when a cloud crosses the level interface.
* The force smoothing and interpolation kernels of all models skip ghost blocks, so the force of the particles is conserved on every level.
* New parameter `fluidGravity` (default 1): 0 applies the gravity only to the particles (a closed box full of liquid, where the hydrostatic pressure does not influence the flow but would vary the lattice density by 3 g L / Cu²).
* Tests: `refinement_euler_particles_conservation`, `refinement_lagrange_force_and_settling`.
### Non-spherical particles (`problemID` 22)

Rigid ellipsoids as Lagrangian point particles (`LagrangeNoSphereParticles`, block `BlockFreeSurfaceLagrangeNoSphereParticles`). Everything except the equations of motion is inherited from `problemID` 2: interpolation, force spreading, walls, cell lists, free surface, refined meshes.

* Input `init/particles.txt`: number of particles, then `x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]` (SI units). `a1..a3` are the semi-axes along the body axes, the quaternion rotates the body axes to the lab frame (it is normalised), `f1..f3` are optional factors of the drag along the body axes (default 1).
* Translation: drag tensor `K = R diag(k1, k2, k3) Rᵀ` with the exact Stokes resistance of the ellipsoid `k_i = f_i μ 16π a1a2a3 / (χ0 + a_i² α_i)` (Oberbeck; the elliptic integrals by Carlson's `R_F`, `R_D`) and the Schiller–Naumann correction with the volume-equivalent diameter. The velocity is updated implicitly (3 × 3 system), the reaction `K (v − u)` goes to the liquid as for the spheres (momentum conserved). A sphere gives exactly the result of `problemID` 2.
* Rotation: Euler equations of the rigid body in the body frame with the creeping-flow torque of Jeffery (1922), `T_i = μ KR_i (Ω_i + λ_i E_jk − ω_i)`, `λ_i = (a_j² − a_k²)/(a_j² + a_k²)`, implicit in the torque, so small particles follow the Jeffery orbits. The orientation is a quaternion advanced with the exact rotation of each step. The velocity gradient of the liquid is computed in every liquid cell (`gradU`, mirror value `2 u_wall − u` behind walls, zero normal gradient at the free surface) and interpolated to the particle with the same kernel as the velocity.
* Walls: the contact distance is the support function of the ellipsoid along the wall normal, `sqrt(nᵀ R diag(a²) Rᵀ n)`. Particle–particle repulsion uses the volume-equivalent diameter.
* Output `output_particles*.vtk`: in addition `TENSORS shape` (`R diag(a) Rᵀ`, for the Tensor Glyph filter of ParaView), `axis1..3` and `angularVelocity` [rad/s].
* Not modelled: inertial (finite Re) torque and lift, the torque of the particle on the liquid (only the force is coupled), contact torques at walls. The point-particle limits of `problemID` 2 hold, the largest semi-axis should be below the cell size.
* Gravity rotation (`gravityRotationOmega`) now works for the Lagrangian particles also on refined meshes (step 202 was missing there).
* Tests `nosphere_*` (`tests/test_nosphere.cpp`): resistance against the closed forms for a sphere and a prolate spheroid (Chwang & Wu) and against quadrature (1e-14), Jeffery period of a 3:1 spheroid (1e-5), the sphere limit identical to `problemID` 2 (difference 0), glide of a settling spheroid inclined by 45° (slip 0.06 % and 0.003° from the theory), rotation in a Couette flow (uniform mesh: 0.02 % from Jeffery; refined: the particle follows Jeffery with the interpolated gradient to 1e-4, the gradient itself is disturbed by ~1.5 % at the level interface), contact of an ellipsoid lying on a wall.

## Not done

* AA-pattern (a single population array) and SoA layout.
* Refinement for the fibers (`problemID` 21).

## Format of `fields_<iter>.bin`

```
char[4] "LBMF", int32 version (1, refined mesh 2, compressed 3), int32 iteration, double time, int32 problemID,
int32 nBlock, int32 ni, nj, nk, double Ch (level 0), int32 nFields, nFields × char[16] name
for each block: int32 bi, bj, bk, [version 2, 3: int32 level]; for each field float32[ni*nj*nk] (index i + ni*j + ni*nj*k)
version 3: everything after the header (from the first block on) is one raw deflate stream (RFC 1951, no zlib header)
fields: ux, uy, uz [m/s], p [Pa], type (cell type), with a free surface also vof, flag
```
