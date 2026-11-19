# parameters.txt - all parameters of cLBMnet

Format: one `name = value` per line, SI units. Integers accept the notation `1e5`. Required parameters are marked
with *. Source of truth: `getPar*` calls in `Block.cpp`, `Mesh.cpp`, `BlockFreeSurface*.cpp`, `LagrangeParticles.cpp`.

## Case and numerics

| Parameter | Default | Meaning |
|---|---|---|
| `problemID`* | | 0 single phase, 1 free surface, 2 Lagrangian particles, 22 Lagrangian ellipsoids (rotating), 3 Euler FVM particles, 33 Euler FVM spheroids with orientation, 4 Euler LBM particles, 42 Euler LBM spheroids with orientation, 21 Lagrangian fibers |
| `model`* | | 319 D3Q19, 315 D3Q15, 29 D2Q9 (requires `nk = 1`, no refinement) |
| `ni`*, `nj`*, `nk`* | | cells per block |
| `H`* | | edge of a level-0 block [m]; cell size `h = H / ni` |
| `uMax`* | | velocity scale [m/s]; `Cu = uMax / 0.1`, time step `h / Cu` |
| `density`* | | liquid density [kg/m³] (pressure scale `ρ Cu² / 3`) |
| `vKin`* | | kinematic viscosity [m²/s] |
| `maxIter`* | | number of level-0 steps |
| `saveRate` | 1e9 | output every N steps (VTK, residuum, free surface) |
| `threads` | 1 | OpenMP threads |
| `collisionModel` | 0 (or `MRT = 1` → 1) | 0 BGK, 1 MRT, 2 TRT |
| `trtMagic` | 3/16 | TRT magic parameter (wall exactly half way) |
| `mrtMagicSq` | 1 | MRT odd rates from the magic relation |
| `tauMin`, `tauMax` | 0.5001, 1.5 | limits of the local τ (non-Newtonian, LES) |
| `tauMinusMax` | 1e30 | upper limit of the odd relaxation time (use 1 for water with LES) |
| `smagorinsky` | 0 | Smagorinsky constant (0.1-0.17) |
| `ddfShift` | 1 | store f − w (better FP32 precision) |
| `xPeriodicity`, `yPeriodicity`, `zPeriodicity` | 0 | periodic directions |
| `gx`, `gy`, `gz` | 0 | gravity [m/s²] |
| `fluidGravity` | 1 | factor of the gravity acting on the liquid (0: particles only, closed full box) |
| `gravityRotationOmega` | 0 | rotation of (gx, gy) around z [rad/s], problemID 2, 4 and 42 |
| `initializeWithHydroPress` | 0 | hydrostatic initial density from the highest liquid cell |

## Boundary conditions

| Parameter | Default | Meaning |
|---|---|---|
| `uInlet`, `vInlet`, `wInlet` | 0 | velocity of inlet cells (type 2) [m/s] |
| `pOutlet` | 0 | relative pressure of outlet cells (type 3) [Pa] |
| `uwall`, `vwall`, `wwall` | 0 | velocity of moving wall cells (type 4) [m/s] |
| `movingLidTop` | 1 | the top face of the domain moves with (`uwall`, `vwall`, `wwall`); harmless while they are 0, set 0 otherwise |

## Free surface (problemID 1-4, 21)

| Parameter | Default | Meaning |
|---|---|---|
| `surfaceTension` | 0 | σ [N/m] |
| `curvatureModel` | 0 | 0 divergence of smoothed normals, 1 height function |
| `fsTolerance` | 0.05 | fill-level tolerance for converting interface cells |
| `wetAngleEq`, `wetAngleAdv`, `wetAngleRec` | 0 | contact angles [deg] |

## Non-Newtonian liquid

| Parameter | Default | Meaning |
|---|---|---|
| `viscosityModel` | 0 | 0 Newtonian, 1 power law, 2 Bingham/Herschel-Bulkley, 3 Cross, 4 polynomial, 5 and 6 modified power law |
| `fluidModelParameters` | | list `b0, b1, ...` of the model (kinematic viscosity [m²/s] as a function of the shear rate [1/s]) |

## Particles

| Parameter | Default | Model | Meaning |
|---|---|---|---|
| `nIntIter` | 1 | 2, 22 | particle sub-steps per level-0 step |
| `maxParticlesInCell` | 500 | 2 | capacity of the cell lists of particles |
| `couplingWayNumber` | 2 | 2 | 4: the fluid equations use the fluid volume fraction 1 − φ |
| `particleRestitution` | 1 | 2 | restitution of the wall-normal velocity (0-1) |
| `cLES` | 0 | 2 | extra viscosity proportional to concentration × shear rate |
| `particleDiameter`*, `particleDensity`* | | 3, 33, 4, 42 | d [m] (33, 42: volume-equivalent sphere), ρp [kg/m³] |
| `particleAspectRatio` | 1 | 33, 42 | r = a / b of the spheroids (> 1 prolate, < 1 oblate) |
| `orientationDiffusion` | 0 | 33, 42 | Folgar-Tucker rotary diffusion C_I (fibre suspensions ~0.001-0.01) |
| `orientationClosure` | 0 | 33, 42 | 0 hybrid, 1 quadratic closure of the orientation equation |
| `particleDiffusivity`* | | 4, 42 | diffusivity of the concentration [m²/s] (TRT τ⁻; recommended 0.001 h_L²/Δt_L, 0 allowed) |
| `frontDiffusivity` | 0 | 4, 42 | extra diffusivity at steep fronts [m²/s] × dimensionless front indicator (0-1) |
| `overpackDiffusivity` | 0 | 4, 42 | extra diffusion above the packing limit (lattice units) |
| `particleInterForces` | 0 | 4, 42 | read, not used |

## Output

| Parameter | Default | Meaning |
|---|---|---|
| `vtkBinary` | 1 | binary VTK XML (0: ASCII) |
| `compressResults` | 1 | zlib-compressed VTK arrays and `fields_<iter>.bin` v3 |
| `writeFields` | 1 | write `fields_<iter>.bin` for LBM Studio |

Fields in `output_*.vtu`: velocity, relative_pressure, VOF, FLAGMAP (1 liquid, 2 gas, 3 interface), curvature,
surfaceLevel, with a non-Newtonian liquid also shearRate and viscosity, with particles concentration and
model-specific fields. problemID 33 and 42 add orientation_tensor, orientation_direction and orientation_alignment.
`output_particles<time>us.vtk`: Lagrangian particle positions in cells of level 0 (problemID 22 also the shape tensor,
body axes and angular velocity).
