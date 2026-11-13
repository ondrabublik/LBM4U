# The particle models - equations and numerics

Lattice units inside the solver (h = Δt = ρ = 1 on level 0); the formulas below are written in physical form. The
liquid is the free-surface LBM (D3Q19, BGK/TRT/MRT); u is its velocity, ν its viscosity, ρ its density.
Files: `LagrangeParticles.cpp`, `BlockFreeSurfaceLagrangeParticles.cpp` (problemID 2),
`BlockFreeSurfaceEulerParticlesFVM.cpp` (3), `BlockFreeSurfaceEulerParticlesLBM.cpp` (4), time stepping in
`Mesh::timeStep`, `Mesh::eulerLevelStep`.

## Contents

1. Drag and particle equation of motion (all models)
2. Coupling to the liquid (all models)
3. Lagrangian point particles
4. Euler FVM
5. Euler LBM
6. Refinement and free surface

## 1. Drag and the equation of motion

For a sphere of diameter d, volume V = πd³/6, mass m = ρp V, cross-section A = πd²/4 and slip velocity
w = v − u:

```
m dv/dt = −½ ρ Cd A |w| w + (m − ρ V) g          Cd = 24/Re (1 + 0.15 Re^0.687) for Re < 1000, else 0.44
Re = d |w| / ν
```

The drag coefficient is evaluated with the old slip, the velocity is updated implicitly (backward Euler):

```
v^{n+1} = ( m v^n / Δt + k u + (m − ρV) g + F_other ) / ( m / Δt + k ),     k = ½ ρ Cd A |w^n|
```

so the update is stable for any particle response time `τp = ρp d² / (18 ρ ν)`, also τp ≪ Δt (the particle then
moves with u + v_t). Gravity can be rotated about z (`gravityRotationOmega`); `fluidGravity` scales only the part
acting on the liquid.

## 2. Coupling to the liquid (two-way)

The reaction `−k (v − u)` (summed over sub-steps as an impulse) is the force density F given to the liquid. The
liquid uses the exact difference method: `acc = fluidGravity · g + F · vof / ρ_cell`, applied in liquid and interface
cells only. The spreading kernels are normalised over the liquid cells they reach (walls, gas and ghost blocks are
skipped and the weights renormalised), so the total force equals the total reaction - momentum of liquid + particles
is conserved to round-off. Lagrange `couplingWayNumber = 4` additionally weights the equilibrium by the fluid
fraction ε = 1 − φ (volume-averaged equations).

## 3. Lagrangian point particles (problemID 2)

* State per particle: position, velocity, d, ρp. Sub-steps per LBM step: `nIntIter`.
* Fluid velocity at the particle: interpolation with the 4-point cosine kernel
  `δ(r) = (1 + cos(π r / 2)) / 4` for |r| < 2 (Peskin), product over x, y, z, liquid cells only, normalised.
* The reaction is spread with the same kernel (minimum kernel sum 0.1 guards a particle surrounded by walls).
* Walls: a particle closer than d/2 to a wall (signed distance function from the wall cells) is moved to d/2 along
  the wall normal, the normal velocity is reflected with `particleRestitution`, tangential motion kept. Near a wall
  the interpolated fluid velocity has no component into the wall.
* Particle-particle: repulsive force `10/Ca · sqrt(1.5 d − L)` for centre distance L < 1.5 d, found through cell
  lists (`maxParticlesInCell`). It prevents overlap; it is not a contact or collision model.
* A particle in gas falls freely (no drag, no force on the liquid).
* Periodic directions wrap positions; inlets and outlets behave like walls for particles.
* Output fields: `concentration` (particles per cell), `particle_volume_fraction`.

Resolution: the particle is a point; its own reaction disturbs the interpolated u. The disturbance scales roughly
like d/h, so keep d < h (d/h ≲ 0.3 for a few % error in the settling velocity).

## 4. Euler FVM (problemID 3)

Unknown: c = number of particles per cell (one particle size and density per run). One particle velocity field
v(x) is computed in every liquid cell from the equation of motion of a single particle (section 1, local form
without the convective term), using the fluid velocity of the cell. The concentration obeys

```
∂c/∂t + ∇·(c v) = 0
```

discretised with fluxes through the 6 faces of each cell: face velocity = mean of the two cells' particle
velocities; flux = upwind flux blended with the Lax-Wendroff flux by the van Leer limiter
`ψ(r) = (r + |r|)/(1 + |r|)` evaluated from the upwind side; first order next to walls and gas. Every face flux is
computed once and applied with opposite signs to both cells (exactly conservative). The force on the liquid is
`c · k (v − u)`, smoothed by a conservative kernel. There is no packing limit (`PP = 0`).

## 5. Euler LBM (problemID 4)

The concentration is carried by a second set of populations g_i on the same lattice:

```
g_i(x + e_i, t + 1) = g_i − (g_i − g_i^eq) / τ_i
g_i^eq = w_i c (1 + 3 θ e_i·v),  g_0^eq = θ c − Σ_{i≠0} g_i^eq          (θ ≈ 1/1.4, limited for D3Q15)
```

which reproduces `∂c/∂t + ∇·(c v) = ∇·(D ∇c)` with `D = (τ − ½) / (3θ)` in lattice units; τ is set from
`particleDiffusivity`. The relaxation time is per direction pair:
`τ_i = τ_s + |c_{+} − 2c + c_{−}| / 1000 + max(c/c_max − 1, 0) · overpackDiffusivity`, where `c_max` is the 0.74
packing; the second term is an artificial diffusion at steep fronts, the third spreads overpacked regions. Because τ
differs between directions, c is reconstructed so that the collision conserves Σg exactly
(`c = Σ (g_i/τ_i) / Σ (w_i'/τ_i)`).

The particle velocity v solves the single-particle equation of motion with an explicit convective term
`(v·∇)v`; next to walls only the velocity towards the wall is limited to the liquid velocity. Force on the liquid as
in FVM. D below ~`0.005 h Cu` tends to be unstable; larger D smears the solution (diffusion length `sqrt(2 D T)`).

## 6. Refinement and free surface

* Lagrange: each particle is coupled to the finest computed block that contains it; the impulse of one level-0
  particle step acts in all 2^L steps of level L as the force `J · 4^L` per step. Kernel cells in ghost blocks are
  skipped. Particle-particle search does not cross the level interface.
* Euler: the particle step runs on every level between propagation and collision. c is the number in a cell of the
  level; a ghost cell gets 1/8 of the coarser cell, a coarse ghost the sum of its 8 fine cells. FVM: the flux through
  a level-interface face is computed by the fine cells and the coarse cell receives exactly the same number. LBM:
  g is exchanged like the liquid populations (explode with 1/8, coalesce as the sum) with a balance correction.
* Free surface: when a cell stops being liquid, its particles are handed to liquid neighbours (FVM, LBM); gas
  cells hold no concentration. c is per cell also in interface cells (not scaled by the fill level).
