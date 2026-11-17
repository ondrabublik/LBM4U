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
7. Non-spherical particles (problemID 22)
8. Eulerian spheroids (problemID 33)

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

## 7. Non-spherical particles (problemID 22)

Rigid ellipsoids, semi-axes a_i along the body axes, rotation matrix R (body → lab) from a unit quaternion q.

```
chi0 = a1a2a3 ∫0^∞ dλ/Δ = 2 a1a2a3 R_F(a1², a2², a3²)
alpha_i = a1a2a3 ∫0^∞ dλ/((a_i²+λ)Δ) = 2/3 a1a2a3 R_D(a_j², a_k², a_i²),   Δ = sqrt((a1²+λ)(a2²+λ)(a3²+λ))
KT_i = 16π a1a2a3 / (chi0 + a_i² alpha_i)                         (translation, sphere: 6πa)
KR_i = 16π a1a2a3 (a_j² + a_k²) / (3 (a_j² alpha_j + a_k² alpha_k))   (rotation, sphere: 8πa³)
```

Translation (implicit, lab frame): `(m/Δt I + K) v^{n+1} = m v^n/Δt + K u + (m − ρV) g + F_other`,
`K = R diag(f_i μ KT_i g(Re)) Rᵀ`, `g(Re) = 1 + 0.15 Re^0.687` (Re with the volume-equivalent diameter, 0.44 Re/24
above Re = 1000). The reaction `K (v − u)` is spread to the liquid as in section 3.

Rotation (body frame, implicit in the torque): `I_i dω_i/dt + (I_k − I_j) ω_j ω_k = μ KR_i (ω^J_i − ω_i)` with the
Jeffery angular velocity `ω^J_i = Ω_i + (a_j² − a_k²)/(a_j² + a_k²) E_jk` (Ω half vorticity, E strain rate, both
transformed to the body frame), `I_1 = m (a2² + a3²)/5`. Then `q ← q ⊗ (cos(|ω|Δt/2), ω/|ω| sin(|ω|Δt/2))`.
For small particles the relaxation time `I/(μ KR)` is far below Δt and ω = ω^J (Jeffery orbits).

Velocity gradient per liquid cell: `∂u_a/∂x_b = 3 Σ_i w_i e_ib u_a(x + e_i)`; behind a wall the mirror value
`2 u_wall − u`, behind gas or an outlet the cell value. Interpolated with the velocity kernel (scaled by 2^L to cells of
level 0 on refined meshes). Wall contact uses the support function `sqrt(nᵀ R diag(a²) Rᵀ n)`.

## 8. Eulerian spheroids (problemID 33)

Per cell: c (number of particles) and A = <p p> (p the symmetry axis, trace A = 1). Semi-axes a = r b along p,
b across, with a b² = (d/2)³.

```
<M> = I / KT_perp + (1/KT_par - 1/KT_perp) A          (exact average of the Stokes mobility of a spheroid)
(m I + K) v_new = m v + K u + (m - ρV) g,   K = <M>^-1 μ g(Re)       force on the liquid c K (v - u)
```

Transport of q = c A with the flux of c (limited upwind/Lax-Wendroff): flux_q = flux_c · A_upwind; then
A = q / trace(q). Orientation (Folgar-Tucker), RK4 in each liquid cell, L = ∇u, D, W its symmetric and antisymmetric
parts, γ̇ = sqrt(2 D:D), λ = (r² - 1)/(r² + 1):

```
dA/dt = W A - A W + λ (D A + A D - 2 A4:D) + 2 C_I γ̇ (I - 3A)
A4:D quadratic = A (A:D);  linear = -(δ trD + 2D)/35 + (A trD + 2(AD + DA) + δ A:D)/7
hybrid = f quadratic + (1 - f) linear,  f = 1 - 27 det A
```

The hybrid closure is exact for aligned (f = 1) and random (f = 0) particles. After each step A is renormalised to trace
1 and negative eigenvalues are cut. Refinement: q is exchanged with c at the level interface (exchangeOrient), coarse
ghosts get the c-weighted mean of A.
