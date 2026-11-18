---
name: lbm4u-euler-fvm
description: Set up, run, check and evaluate the Eulerian finite-volume particle model (problemID 3) of the cLBMnet / LBM4U free-surface solver - particle concentration per cell, particleDiameter and particleDensity, init/concentration files, conservative upwind / Lax-Wendroff fluxes with the van Leer limiter, one particle velocity field from the drag balance, two-way coupling, refined meshes, sharp settling fronts. Use it whenever the user works with problemID 3, simulates sedimentation, settling clouds or monodisperse suspensions with a concentration field, or asks how the FVM particle transport works or how accurate it is - also in Czech ("Eulerovy částice", "koncentrace částic", "FVM", "konečné objemy", "sedimentace", "usazování", "problemID 3"). Use with lbm4u-particles (choice between models) and lbm4u-free-surface (case files, running).
---

# Eulerian particles with finite volumes (problemID 3)

The particles are a concentration c = number of particles per cell, all of one diameter d and density ρp. One
particle velocity v is computed in every liquid cell from the drag-gravity balance of a single particle, and c is
transported by conservative finite-volume fluxes. It is the most accurate model for monodisperse settling
(sharp fronts, little numerical diffusion) and much cheaper than Lagrange for many particles. Case files and running:
**lbm4u-free-surface**; choice of model: **lbm4u-particles**. Code: `BlockFreeSurfaceEulerParticlesFVM.cpp`,
`Mesh::eulerLevelStep`. Equations: `.claude/skills/lbm4u-particles/references/models.md` sections 1, 2, 4, 6 and
`docs/math/matematicky_popis.pdf`.

## When to use it

* Dilute to moderate monodisperse suspensions, settling clouds, sedimentation fronts, particle fields in a whole part.
* Quantitative work where the front sharpness matters (prefer it over problemID 4).
* Not for: polydisperse mixtures (one size per run - run several cases or use Lagrange), trajectories, dense packings
  (no packing limit), non-spherical particles (problemID 33).

## Parameters (`parameters.txt`)

| Parameter | Default | Meaning |
|---|---|---|
| `problemID` | – | `3` |
| `particleDiameter` | required | d [m] |
| `particleDensity` | required | ρp [kg/m³] |
| `fluidGravity` | 1 | factor of gravity on the liquid; 0 = gravity on the particles only (closed full box) |
| `gravityRotationOmega` | 0 | rotation of (gx, gy) about z [rad/s] |

## Concentration file `init/concentration<tag>.txt`

The **number of particles in each cell**, one value per cell, in the block file order (k fastest), one file per block
(`<tag>` as in the other block files). Missing file = no particles in the block. For a volume fraction φ:
`c = φ · h_L³ / (π d³/6)` with h_L the cell size of the block's level. `make_case.py`:
`"concentration": {"sphere": [x, y, z, R] | "box": [...], "d": 1e-4, "phi": 0.005}` (d must match
`particleDiameter`). Cells must be liquid; gas cells hold no particles.

## Physics and numerics (summary)

* Particle velocity per cell: the single-particle equation with Schiller-Naumann drag, implicit (stable for any τp),
  local form without the convective term, with the fluid velocity of the cell.
* Transport `∂c/∂t + ∇·(c v) = 0`: face velocity = mean of the two cells; upwind flux blended with Lax-Wendroff by the
  van Leer limiter; first order next to walls and gas. Each face flux is applied with opposite signs to both cells -
  exactly conservative.
* Force on the liquid `c k (v − u)` spread by a conservative kernel; momentum conserved.
* Free surface: particles of a cell that stops being liquid are handed to liquid neighbours. c is per cell also in
  interface cells (not scaled by the fill level).

## Accuracy and checks

* CFL: the particle moves less than a cell per step in practice (|v| ≲ uMax); check the console for growing c.
* Studies (local only): a settling cloud falls 2-8 % slower than with Lagrange (point particles feel part of their
  own disturbance and fall slightly faster), 18-33 % faster than with Euler LBM; on an equidistant mesh sphere cloud 6.68 mm/s, mesh h/2 +4.5 %.
* The particle count is conserved to round-off (also across refinement levels and drying cells).

## Refined meshes

c is the number in a cell of its level; a fine ghost cell gets 1/8 of the coarse cell, a coarse ghost the sum of its 8
fine cells. The flux through a level-interface face is computed by the fine cells and the coarse cell receives the
same number - exact conservation.

## Output

`output_*.vtu`: `concentration` (number per cell), `particle_velocity`, `particle_pressure`,
`particle_pressure_gradient`;
volume fraction = concentration · π d³/6 / h_L³. Evaluate centroid height and velocity, spread σ, deposited fraction
with a short Python script.

## Limits

One size and density per run, one velocity per cell (no crossing trajectories), no packing limit (concentration can
exceed close packing at the bottom), no Brownian diffusion, inlets do not bring particles in.

Tests: `lbm_tests.exe euler_fvm`, `euler_force_spreading_is_conservative`, `euler_particles_follow_free_surface`,
`refinement_euler_particles_conservation`.
