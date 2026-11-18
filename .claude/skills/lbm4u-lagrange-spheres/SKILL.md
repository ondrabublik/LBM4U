---
name: lbm4u-lagrange-spheres
description: Set up, run, check and evaluate the Lagrangian spherical point-particle model (problemID 2) of the cLBMnet / LBM4U free-surface solver - init/particles.txt, polydisperse spheres with their own diameter and density, implicit Schiller-Naumann drag, two-way coupling (couplingWayNumber), nIntIter sub-steps, wall contact and particleRestitution, particle repulsion and maxParticlesInCell, refined meshes, output_particles VTK. Use it whenever the user works with problemID 2, tracks individual particles or trajectories, wall impacts or deposition of spheres, or asks how the Lagrangian spheres work or how accurate they are - also in Czech ("Lagrangeovy částice", "kulové částice", "trajektorie částic", "problemID 2"). Use with lbm4u-particles (choice between models) and lbm4u-free-surface (case files, units, running).
---

# Lagrangian spheres (problemID 2)

Every particle is a point with its own position, velocity, diameter d and density ρp, moved through the free-surface
LBM liquid. It is the reference model for dilute and polydisperse suspensions and the only one that gives
trajectories. Case files, units and running are in **lbm4u-free-surface**; the choice between the particle models is in
**lbm4u-particles**. Code: `LagrangeParticles.cpp` (particles), `BlockFreeSurfaceLagrangeParticles.cpp` (liquid side),
`Mesh::timeStep`. Equations: `.claude/skills/lbm4u-particles/references/models.md`, sections 1-3 and 6, and the
mathematical description `docs/math/matematicky_popis.pdf`.

## When to use it

* Dilute suspensions (volume fraction up to a few %), particles of different sizes or densities in one run.
* Trajectories, residence times, wall impacts, where particles end up in a mould.
* d clearly smaller than the cell, d/h ≲ 0.3 (point-particle limit, see Accuracy).
* Not for: many millions of particles (cost ∝ N, use problemID 3 or 4), dense packings (no contact model),
  non-spherical particles (problemID 22).

## Parameters (`parameters.txt`)

| Parameter | Default | Meaning |
|---|---|---|
| `problemID` | – | `2` |
| `nIntIter` | 1 | particle sub-steps per LBM step of level 0 (raise for fast, light particles or strong repulsion) |
| `couplingWayNumber` | 2 | 4 = two-way coupling + fluid fraction 1 − φ in the liquid equilibrium; any other value = two-way coupling only (there is no one-way option) |
| `particleRestitution` | 1.0 | wall-normal restitution 0-1 (0 = particles stop at the wall, 1 = elastic) |
| `maxParticlesInCell` | 500 | capacity of the cell lists for the repulsion; raise it for dense packings |
| `fluidGravity` | 1 | factor of gravity acting on the liquid; 0 = gravity on the particles only (closed full box) |
| `gravityRotationOmega` | 0 | rotation of (gx, gy) about z [rad/s] (rotating moulds) |

Gravity `gx, gy, gz`, `density`, `vKin` are the liquid parameters of the free-surface case.

## Particle file `init/particles.txt`

```
N
x y z d rho          (one line per particle, SI units, positions in metres in the domain coordinates)
```

Positions must lie in liquid cells (a particle in gas falls freely until it reaches the liquid). `make_case.py`
writes it: `"particles": {"sphere": [x, y, z, R] | "box": [x0, y0, z0, x1, y1, z1], "d": 1e-4, "rho": 2500,
"phi": 0.005, "seed": 1}` or `{"list": [[x, y, z, d, rho], ...]}`. Report the number of particles to the user - it
sets the cost.

## Physics and numerics (summary)

* Drag: Schiller-Naumann `Cd = 24/Re (1 + 0.15 Re^0.687)` (0.44 above Re = 1000), velocity updated implicitly -
  stable for any particle response time τp = ρp d² / (18 ρ ν), also τp ≪ Δt.
* Gravity with buoyancy (m − ρV) g. Fluid velocity at the particle: 4-point cosine (Peskin) kernel over liquid cells.
* Two-way coupling: the reaction is spread with the same kernel, normalised over the liquid cells it reaches - momentum
  of liquid + particles is conserved to round-off.
* Walls: a particle closer than d/2 is pushed back along the wall normal, normal velocity × `particleRestitution`.
* Particle-particle: repulsion `10/Ca · sqrt(1.5 d − L)` for centre distance L < 1.5 d - a safeguard against
  overlap, not a collision model.
* Inlets and outlets act as walls for particles; inlets do not bring particles in.

## Accuracy and checks

* A point particle feels part of its own disturbance: the fluid velocity at the particle is ~0.26 v_t at d/h = 0.64
  and ~0.06 v_t at d/h = 0.32. Keep d/h ≲ 0.3 or accept a few % faster settling; on refined meshes particles in
  the fine region settle slightly faster.
* Terminal velocity of one particle (Stokes) `v_t = (ρp − ρ) g d² / (18 μ)`; check Re_p = d v_t / ν, for Re_p > 1
  solve the Schiller-Naumann balance.
* Studies: a settling cloud falls 2-8 % faster than with Euler FVM (problemID 3); on an equidistant mesh sphere cloud
  7.17 mm/s vs FVM 6.68 mm/s, mesh h/2 changes it by 3.5 % (`docs/studie_*.pdf`, local only).
* The particle count is exact; the console prints particle statistics every output.

## Refined meshes

Each particle is coupled to the finest computed block that contains it; the impulse of one level-0 particle step acts
in all 2^L steps of level L. Repulsion is not evaluated across a level interface - keep dense packings in one level.

## Output

* `output_particles<time>us.vtk`: positions in **cells of level 0** (multiply by h0 for metres), `diameter`,
  `densityRatio`, `particleDensity`, vectors `rychlost` (particle velocity), `gradP`, `force`.
* `output_*.vtu`: `concentration` (particles per cell), `particle_volume_fraction`.
* Typical evaluation (short Python script): centroid height and velocity, spread σ, fraction deposited on walls,
  comparison with v_t.

## Limits

No rotation, no lubrication, no contact model; point particles only (d < h); inlets do not inject particles;
gravity on particles and liquid is the same vector (except `fluidGravity`).

Tests: `lbm_tests.exe lagrange`, `refinement_lagrange_force_and_settling`, `particles_gravity_rotation`
(`tests/test_particles.cpp` has small example setups).
