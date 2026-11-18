---
name: lbm4u-particles
description: Explain, choose, set up and evaluate particle-laden free-surface flows in the cLBMnet / LBM4U solver - Lagrangian point particles (spheres problemID 2, rotating ellipsoids / non-spherical particles problemID 22), Eulerian concentration by finite volumes (problemID 3, FVM; spheroids with an orientation tensor problemID 33) and by LBM populations (problemID 4; spheroids problemID 42) - settling, suspensions, particle clouds, orientation and Jeffery rotation, two-way coupling, octree refinement, particles.txt and concentration files. Use it whenever the user asks about particles, sediment, fillers, flakes, short fibres, concentration or volume fraction in a simulation of this repository, compares the Lagrange and Euler approaches or wants the particle methods described - also in Czech ("částice", "partikly", "nekulové částice", "elipsoidy", "usazování", "sedimentace", "koncentrace", "popis metod"). Use with lbm4u-free-surface (case files, units, running) and, for non-Newtonian carriers, lbm4u-non-newtonian.
---

# Particles in free-surface flow (cLBMnet / LBM4U)

Three models share the liquid solver (free surface, refinement, rheology) and differ in how the particles are
represented. Case files, units and running are in **lbm4u-free-surface** - read it first. The equations of all
three models are in `references/models.md` and in `docs/math/matematicky_popis.pdf`; read them when the user asks how a model works or when you need to
reason about accuracy.

This skill compares the models and helps to choose one. Once the model is chosen, use its own skill for the
parameters, input files, checks and limits:

| problemID | Model | Skill |
|---|---|---|
| 2 | Lagrange, spheres | **lbm4u-lagrange-spheres** |
| 22 | Lagrange, rotating ellipsoids | **lbm4u-lagrange-ellipsoids** |
| 3 | Euler FVM, spheres | **lbm4u-euler-fvm** |
| 33 | Euler FVM, spheroids with orientation tensor | **lbm4u-euler-fvm-spheroids** |
| 4 | Euler LBM, spheres | **lbm4u-euler-lbm** |
| 42 | Euler LBM, spheroids with orientation tensor | **lbm4u-euler-lbm-spheroids** |

## Which model

| | Lagrange (`problemID 2`) | Euler FVM (`3`) | Euler LBM (`4`) |
|---|---|---|---|
| Particles are | individual points with own d, ρp | number per cell, one size and density | number per cell, one size and density |
| Input | `init/particles.txt` | `init/concentration<tag>.txt` | `init/concentration<tag>.txt` |
| Transport | trajectories, implicit drag | conservative upwind/Lax-Wendroff fluxes, van Leer limiter | scalar LBM populations, BGK with direction-dependent τ |
| Numerical diffusion | none | small (limiter) | prescribed `particleDiffusivity`, needed for stability |
| Packing limit | repulsion between particles (not a contact model) | none | `overpackDiffusivity` above 0.74 packing |
| Cost | ∝ number of particles | ~2× liquid alone | ~2.5-3× liquid alone |
| Best for | dilute, polydisperse, trajectories, wall impacts, d ≪ h | dilute to moderate monodisperse suspensions, sharp fronts | when smooth, diffusive transport is acceptable or wanted |

Practical guidance from the tests and the settling-cloud study (21 runs, three models on a half-refined box):

* Lagrange and FVM agree to 2-8 % in cloud velocity; LBM is 18-33 % slower because its diffusivity smears the
  cloud. Prefer FVM for quantitative monodisperse work; use LBM only with a small `particleDiffusivity`.
* A Lagrangian particle "feels" part of its own disturbance; this grows as the cell shrinks (fluid velocity at the
  particle 0.26 v_t at d/h = 0.64, 0.06 v_t at d/h = 0.32). Keep d clearly below h (d/h ≲ 0.3) or accept a few %
  of extra settling speed; on a refined mesh this makes particles in the fine region settle slightly faster.
* All models conserve the particle count: Lagrange exactly, FVM to round-off, LBM to ~1e-8 (also across refinement
  levels and when liquid cells empty - particles of a drying cell are handed to liquid neighbours).

## Non-spherical particles (`problemID 22`)

Rigid ellipsoids with semi-axes a1, a2, a3, orientation and rotation; everything else (interpolation, coupling,
walls, refinement) is inherited from problemID 2. Use it for elongated or flat particles (fibres shorter than a cell,
flakes, grains) when orientation-dependent drag (gliding, sideways drift) or alignment in shear matters.

* `init/particles.txt`: number of particles, then `x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]` (SI). The
  quaternion rotates the body axes to the lab frame (`1 0 0 0` = axes along x, y, z); `f1..f3` multiply the drag along
  the body axes (default 1, use them to calibrate to measured drag). `make_case.py` writes it with `"shape"`.
* Drag: exact Stokes resistance of the ellipsoid along each body axis × Schiller-Naumann correction (volume-equivalent
  d); a 2:1 prolate spheroid has ~13 % less drag along its axis than across, a 45° tilt gives a glide of ~4°.
* Rotation: Jeffery torque (creeping flow) with the interpolated fluid velocity gradient; in shear the particles tumble
  on Jeffery orbits with period `2π/γ (r + 1/r)`. Inertial torques at finite Re (e.g. falling discs turning broadside)
  are **not** modelled, and the particle torque is not returned to the liquid.
* Keep the **largest** semi-axis below the cell size (point-particle limit). Longer fibres need another model.
* Output: `output_particles*.vtk` has `TENSORS shape` (ParaView: Tensor Glyph with a sphere source shows the
  ellipsoids), `axis1..3`, `angularVelocity`.

## Eulerian spheroids (`problemID 33` FVM, `42` LBM)

The Euler models with an orientation per cell: c particles of one shape (spheroid, aspect ratio r) and their orientation
tensor A = <p p> (trace 1; I/3 random, p p aligned). Use it instead of problemID 22 when there are too many particles
to track, or when the orientation field (alignment by the flow, flow-induced anisotropy) is the result.

* Parameters: `particleDiameter` (volume-equivalent), `particleDensity`, `particleAspectRatio`, `orientationDiffusion`
  C_I (0 for dilute, ~0.001-0.01 for interacting fibres), `orientationClosure` (0 hybrid, default). Initial A:
  `init/orientation<tag>.txt` (`A11 A22 A33 A12 A13 A23` per cell), default random; `make_case.py`:
  `"concentration": {..., "orientation": "random" | [px, py, pz] | [A11, ..., A23]}`.
* Drag: mobility averaged over the orientations (exact for Stokes drag): random prolate r = 3 settles ~10 % slower
  than the volume-equivalent sphere; aligned inclined particles glide sideways as one particle of problemID 22.
* Orientation: Folgar-Tucker equation (Jeffery rotation + rotary diffusion) with the liquid velocity gradient,
  transported with the particle flux.
* Output: `orientation_tensor` (9 components), `orientation_direction` (main axis), `orientation_alignment`
  (largest eigenvalue, 1/3 random ... 1 aligned).
* Limits: one velocity per cell (no dispersion from different orientations), one shape per run, closure error for
  partly aligned states, Jeffery (creeping flow) rotation, no particle stress in the liquid.
* problemID 42 = the same with the LBM transport of problemID 4 (`particleDiffusivity` etc. as for 4). The net
  particle flux of each link carries A of the upwind cell (as 33), so A moves and mixes only with the particles
  (also by their diffusion); r = 1 gives problemID 4. Choose 33 vs 42 as 3 vs 4: 42 smears the cloud and its orientation more.
  The fluid fraction 1 - cV of the LBM model makes the liquid velocity vary through a cloud and the spheroids turn
  slightly in that gradient (absent in 33).

## Physics common to all models

* Drag: Schiller-Naumann `Cd = 24/Re (1 + 0.15 Re^0.687)`, `Re = d |u − v| / ν`, velocity updated implicitly
  (stable for any particle response time).
* Gravity with buoyancy `(m − ρ V) g`. `gravityRotationOmega` rotates (gx, gy) about z (rotating moulds,
  problemID 2, 4 and 42). `fluidGravity = 0` puts gravity on particles only - right for a closed box completely full of
  liquid (the hydrostatic pressure only compresses the lattice), wrong whenever there is a free surface.
* Two-way coupling: the reaction of the drag is spread to the liquid (force `F · vof / ρ`, conservative kernel).
  Lagrange `couplingWayNumber = 4` additionally uses the fluid volume fraction 1 − φ in the liquid equations.
* Terminal velocity of one particle (Stokes): `v_t = (ρp − ρ) g d² / (18 μ)`; check Re_p = d v_t / ν - for Re_p > 1
  solve the Schiller-Naumann balance instead. Clouds and suspensions settle collectively faster or slower than v_t.

## Setting up

1. Start from a free-surface case (lbm4u-free-surface), set `problemID` and the model parameters.
2. Particles:
   * Lagrange: `init/particles.txt` = number of particles, then `x y z diameter density` per line (SI). Positions
     must lie in liquid cells. Optional `nIntIter` (sub-steps), `particleRestitution` (0-1, wall-normal
     restitution), `maxParticlesInCell` (cell-list capacity; raise it for dense packing).
   * Euler: `particleDiameter`, `particleDensity`; for LBM also `particleDiffusivity` [m²/s] and optionally
     `overpackDiffusivity`. `init/concentration<tag>.txt` holds the **number of particles in each cell** (one value
     per cell, block file order, k fastest). For a volume fraction φ: `c = φ · h_L³ / (π d³/6)` with `h_L` the cell
     size of the block's level. Missing file = no particles in the block.
   * `make_case.py` writes both (`"particles"` and `"concentration"` keys with a sphere, box or expression region).
3. `particleDiffusivity` for LBM: the smallest stable value is about `0.005 · h · Cu` (cell Péclet ~10); the
   diffusion length `sqrt(2 D T)` should be small compared with the cloud - report it to the user.
4. Check the run before a long computation: the console shows the particle count (`nP`, `maxP` for problemID 4 and 42,
   statistics for problemID 2). A changing count means particles leave through an inlet/outlet or a setup error.

## Refined meshes

All three models run on octree meshes. Lagrangian particles live on the finest computed block that holds them;
Euler concentrations are per cell of each level with exact exchange at the level interface. Particle-particle
repulsion (Lagrange) is not evaluated across the level interface - keep dense packings inside one level.

## Outputs and evaluation

* `output_*.vtu`: `concentration` (number per cell), Lagrange also `particle_volume_fraction`, FVM
  `particle_pressure`; volume fraction = concentration · π d³/6 / h_L³.
* `output_particles<time>us.vtk` (Lagrange): positions in **cells of level 0** (multiply by h0 for metres),
  `diameter`, `densityRatio`.
* Typical quantities: centroid height and velocity of the particle mass, spread σ, fraction deposited on walls,
  comparison with v_t. Compute them from the VTK/positions with a short Python script.

## Known limits

* Lagrange: point particles - not resolved, no rotation, no lubrication; the repulsion
  `10/Ca · sqrt(1.5 d − L)` is a numerical safeguard, not a contact model. Inlets and outlets are walls for particles.
* Euler: one particle size per run; concentration is per cell also in interface cells (not scaled by fill level).
  FVM has no packing limit. LBM: `particleInterForces` is read but not used; buoyancy counted twice at large
  volume fractions.
* Fibres (`problemID 21`) are a separate, unreviewed model without refinement support.

Unit tests: `lbm_tests.exe lagrange`, `lbm_tests.exe euler` (includes `euler_nosphere` and `euler_nosphere_lbm`), `lbm_tests.exe nosphere`, `lbm_tests.exe refinement_` (settling velocity,
momentum, conservation, wall contact) - `tests/test_particles.cpp` is also the best source of small example setups.
