---
name: lbm4u-euler-lbm
description: Set up, run, check and evaluate the Eulerian lattice-Boltzmann particle model (problemID 4) of the cLBMnet / LBM4U free-surface solver - particle concentration carried by a second set of LBM populations, particleDiameter, particleDensity, particleDiffusivity, frontDiffusivity, overpackDiffusivity above close packing, TRT collision of the scalar, init/concentration files, two-way coupling, refined meshes, rotating moulds. Use it whenever the user works with problemID 4, wants diffusive particle transport or packing limits, chooses particleDiffusivity, or asks how the LBM particle transport works or why it smears a cloud - also in Czech ("Eulerovy částice LBM", "koncentrace", "difuzivita částic", "problemID 4"). Use with lbm4u-particles (choice between models) and lbm4u-free-surface (case files, running).
---

# Eulerian particles with LBM populations (problemID 4)

The particle concentration c (number per cell, one size and density) is carried by a second set of populations g_i on
the liquid lattice - an advection-diffusion LBM (TRT) with the particle velocity v. It is simple, local and robust and
has a packing limit; its accuracy is set by `particleDiffusivity` (smearing) and by the undershoots of the scheme at
sharp fronts. Case files and running:
**lbm4u-free-surface**; choice of model: **lbm4u-particles**. Code: `BlockFreeSurfaceEulerParticlesLBM.cpp`,
`Mesh::eulerLevelStep`. Equations: `.claude/skills/lbm4u-particles/references/models.md` sections 1, 2, 5, 6 and
`docs/math/matematicky_popis.pdf`.

## When to use it

* Suspensions where particle diffusion is physical or acceptable (dispersion, mixing), packing above 0.74 must be
  limited, or a smooth concentration field is wanted.
* Sharp clouds and fronts: usable with a small D (settling cloud within 5-10 % of FVM), but FVM (problemID 3) stays
  sharper and has no undershoots. Not for polydisperse mixtures, trajectories, non-spherical particles (problemID 42).

## Parameters (`parameters.txt`)

| Parameter | Default | Meaning |
|---|---|---|
| `problemID` | – | `4` |
| `particleDiameter` | required | d [m] |
| `particleDensity` | required | ρp [kg/m³] |
| `particleDiffusivity` | required | D [m²/s], sets τ⁻ of the antisymmetric (flux) part of g; 0 is allowed |
| `frontDiffusivity` | 0 | extra D [m²/s] at steep fronts, × dimensionless indicator `|c₊ − 2c + c₋| / (|c₊| + 2|c| + |c₋|)` (1/3 at the edge of a step) |
| `overpackDiffusivity` | 0 | extra relaxation above close packing (c/c_max − 1), spreads overpacked regions |
| `particleInterForces` | 0 | read but not used |
| `fluidGravity` | 1 | factor of gravity on the liquid; 0 = gravity on the particles only |
| `gravityRotationOmega` | 0 | rotation of (gx, gy) about z [rad/s] (rotating moulds) |

### Choosing `particleDiffusivity`

* The TRT scheme is stable down to D = 0; the lower limit is the accuracy: below `D_lat = D Δt_L / h_L² ≈ 0.001`
  (on the level that holds the particles) sharp fronts get undershoots (c < 0, up to 13 % of c at D = 0) and
  overshoots of 13-17 %. Recommended: `D ≈ 0.001 · h_L² / Δt_L` (cloud benchmark: 2.5e-7 m²/s at h = 0.5 mm,
  Δt = 1 ms; undershoot 1.4 % of c at the start, 0.05 % later). A 4× smaller D with `frontDiffusivity ≈ 3 ·
  D_recommended` smears less (cloud 5 % instead of 10 % slower than FVM) at the price of an 11 % overshoot.
* Before the TRT change (BGK) D_lat ≳ 0.005 was needed; old cases with that value still run, but smear 5× more
  than necessary - lower D.
* The diffusion length `sqrt(2 D T)` over the simulated time T must be small compared with the cloud or the features
  of interest - compute it and report it to the user before a long run.
* D is fixed in m²/s: on a finer mesh the same D is a larger lattice value; to keep the same numerical smearing keep
  D_lat, i.e. D ∝ h. The result converges towards FVM as the mesh is refined.

## Concentration file `init/concentration<tag>.txt`

Same as problemID 3: number of particles per cell, block file order (k fastest), missing file = no particles,
`c = φ · h_L³ / (π d³/6)`. `make_case.py`: `"concentration": {"sphere" | "box": ..., "d": ..., "phi": ...}`.

## Physics and numerics (summary)

* TRT collision: `g± = (g_i ± g_ī)/2`, the symmetric part relaxes with τ⁺ = 1 (reset to the equilibrium every step),
  the antisymmetric part with τ⁻ of the pair; `g_i^eq = w_i c (1 + 3θ e_i·v)`, reproducing
  `∂c/∂t + ∇·(c v) = ∇·(D ∇c)`, `D = (τ⁻ − ½)/(3θ)`.
* τ⁻ per direction pair: `particleDiffusivity` + `frontDiffusivity` × front indicator + overpacking term. Σg is
  conserved by the collision for any τ⁻, `c = Σg/θ`.
* Why τ⁺ = 1: BGK (τ⁺ = τ⁻ → ½) over-relaxes the symmetric modes (oscillations, minimum D); a large τ⁺ from a fixed
  "magic" Λ = 1/4 keeps the symmetric non-equilibrium created by bounce-back at walls and slows the transport along
  walls (−8 % diffusion, −9 % settling in a column 8 cells wide).
* Particle velocity: single-particle equation with Schiller-Naumann drag (implicit) and an explicit convective term;
  next to walls only the velocity towards the wall is limited.
* Force on the liquid `c k (v − u)`, conservative kernel. The liquid sees the fluid fraction 1 − cV.

## Accuracy and checks

* Settling cloud benchmark (Nitsche & Batchelor 1997, `simulations/cloud/`, R = 8 cells, refined box): U/U_theory
  FVM 0.764, LBM TRT D = 2.5e-7 0.689 (−10 %), TRT D = 6.25e-8 + frontDiffusivity 7.5e-7 0.729 (−5 %); the old BGK
  with D = 1.25e-6 gave 0.471 (−38 %). σ_z at the end: FVM 2.10, TRT 2.50, TRT + front 2.29, BGK 3.62 mm.
* Older studies (BGK, local only): a cloud fell 18-33 % slower than with Lagrange or FVM. In advection-dominated
  flows (mould filling) it agrees closely with FVM.
* Particle count conserved to ~1e-8 (also across refinement levels and drying cells). The console prints `nP` and
  `maxP`; a changing nP means a setup error or particles leaving through an outlet.

## Refined meshes

g is exchanged at the level interface like the liquid populations (explode with 1/8, coalesce as the sum) with a
balance correction that keeps the count exact.

## Output

`output_*.vtu`: `concentration`, `particle_velocity`; volume fraction = concentration · π d³/6 / h_L³.

## Limits

Smearing by D, small undershoots (c < 0) at sharp fronts, one size per run, one velocity per cell, `particleInterForces` unused, buoyancy is
counted twice at large volume fractions, inlets do not bring particles in.

Tests: `lbm_tests.exe euler_lbm` (settling, diffusion, lattices, velocity gradient, interface cells),
`refinement_euler_particles_conservation`.
