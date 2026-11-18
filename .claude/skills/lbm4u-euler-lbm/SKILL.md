---
name: lbm4u-euler-lbm
description: Set up, run, check and evaluate the Eulerian lattice-Boltzmann particle model (problemID 4) of the cLBMnet / LBM4U free-surface solver - particle concentration carried by a second set of LBM populations, particleDiameter, particleDensity, particleDiffusivity and its stability limit, overpackDiffusivity above close packing, direction-dependent relaxation, init/concentration files, two-way coupling, refined meshes, rotating moulds. Use it whenever the user works with problemID 4, wants diffusive particle transport or packing limits, chooses particleDiffusivity, or asks how the LBM particle transport works or why it smears a cloud - also in Czech ("Eulerovy částice LBM", "koncentrace", "difuzivita částic", "problemID 4"). Use with lbm4u-particles (choice between models) and lbm4u-free-surface (case files, running).
---

# Eulerian particles with LBM populations (problemID 4)

The particle concentration c (number per cell, one size and density) is carried by a second set of populations g_i on
the liquid lattice - an advection-diffusion LBM with the particle velocity v. It is simple, local and robust and has a
packing limit, but needs a minimum diffusivity for stability, which smears clouds and fronts. Case files and running:
**lbm4u-free-surface**; choice of model: **lbm4u-particles**. Code: `BlockFreeSurfaceEulerParticlesLBM.cpp`,
`Mesh::eulerLevelStep`. Equations: `.claude/skills/lbm4u-particles/references/models.md` sections 1, 2, 5, 6 and
`docs/math/matematicky_popis.pdf`.

## When to use it

* Suspensions where particle diffusion is physical or acceptable (dispersion, mixing), packing above 0.74 must be
  limited, or a smooth concentration field is wanted.
* Not for: quantitative settling of sharp clouds or fronts (use problemID 3), polydisperse mixtures, trajectories,
  non-spherical particles (problemID 42).

## Parameters (`parameters.txt`)

| Parameter | Default | Meaning |
|---|---|---|
| `problemID` | – | `4` |
| `particleDiameter` | required | d [m] |
| `particleDensity` | required | ρp [kg/m³] |
| `particleDiffusivity` | required | D [m²/s], sets the relaxation time of g |
| `overpackDiffusivity` | 0 | extra relaxation above close packing (c/c_max − 1), spreads overpacked regions |
| `particleInterForces` | 0 | read but not used |
| `fluidGravity` | 1 | factor of gravity on the liquid; 0 = gravity on the particles only |
| `gravityRotationOmega` | 0 | rotation of (gx, gy) about z [rad/s] (rotating moulds) |

### Choosing `particleDiffusivity`

* The smallest stable value is about `D = 0.005 · h · Cu` (cell Péclet ~10; D_lat = 0.005 in lattice units), h the
  cell size of level 0 and Cu = h/Δt. Below it the solution oscillates or diverges.
* The diffusion length `sqrt(2 D T)` over the simulated time T must be small compared with the cloud or the features
  of interest - compute it and report it to the user before a long run.
* D is fixed in m²/s: on a finer mesh the same D is a larger lattice value; to keep the same numerical smearing keep
  D_lat, i.e. D ∝ h. The result converges towards FVM as the mesh is refined.

## Concentration file `init/concentration<tag>.txt`

Same as problemID 3: number of particles per cell, block file order (k fastest), missing file = no particles,
`c = φ · h_L³ / (π d³/6)`. `make_case.py`: `"concentration": {"sphere" | "box": ..., "d": ..., "phi": ...}`.

## Physics and numerics (summary)

* `g_i(x + e_i, t + 1) = g_i − (g_i − g_i^eq)/τ_i`, `g_i^eq = w_i c (1 + 3θ e_i·v)`, reproducing
  `∂c/∂t + ∇·(c v) = ∇·(D ∇c)`, `D = (τ − ½)/(3θ)`.
* τ per direction pair: base from `particleDiffusivity` + artificial diffusion at steep fronts + overpacking term;
  c is reconstructed so the collision conserves Σg exactly.
* Particle velocity: single-particle equation with Schiller-Naumann drag (implicit) and an explicit convective term;
  next to walls only the velocity towards the wall is limited.
* Force on the liquid `c k (v − u)`, conservative kernel. The liquid sees the fluid fraction 1 − cV.

## Accuracy and checks

* Studies (local only): a settling cloud falls 18-33 % slower than with Lagrange or FVM because of the smearing
  (equidistant mesh: 5.29 mm/s vs FVM 6.68 mm/s, σ_z 0.34 cm vs 0.20 cm); on mesh h/2 with the same D_lat the
  gap to FVM shrinks from 21 % to 12 %. In advection-dominated flows (mould filling) it agrees closely with FVM.
* Particle count conserved to ~1e-8 (also across refinement levels and drying cells). The console prints `nP` and
  `maxP`; a changing nP means a setup error or particles leaving through an outlet.

## Refined meshes

g is exchanged at the level interface like the liquid populations (explode with 1/8, coalesce as the sum) with a
balance correction that keeps the count exact.

## Output

`output_*.vtu`: `concentration`, `particle_velocity`; volume fraction = concentration · π d³/6 / h_L³.

## Limits

Minimum diffusivity (smearing), one size per run, one velocity per cell, `particleInterForces` unused, buoyancy is
counted twice at large volume fractions, inlets do not bring particles in.

Tests: `lbm_tests.exe euler_lbm` (settling, diffusion, lattices, velocity gradient, interface cells),
`refinement_euler_particles_conservation`.
