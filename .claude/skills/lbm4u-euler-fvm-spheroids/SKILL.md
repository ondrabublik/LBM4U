---
name: lbm4u-euler-fvm-spheroids
description: Set up, run, check and evaluate the Eulerian non-spherical particle model with finite volumes (problemID 33) of the cLBMnet / LBM4U free-surface solver - concentration of spheroids (aspect ratio particleAspectRatio) with an orientation tensor A per cell, orientation-averaged Stokes mobility, Folgar-Tucker equation with rotary diffusion orientationDiffusion and the hybrid closure orientationClosure, init/orientation files, orientation_tensor / orientation_direction / orientation_alignment output, fibre orientation in mould filling. Use it whenever the user works with problemID 33, wants the orientation field of short fibres, flakes or elongated particles in a whole part, flow-induced alignment or anisotropy, or asks about the closure error - also in Czech ("nekulové částice Euler", "orientační tenzor", "orientace vláken", "Folgar-Tucker", "problemID 33"). Use with lbm4u-euler-fvm (transport and common parameters), lbm4u-particles (choice between models) and lbm4u-free-surface.
---

# Eulerian spheroids with finite volumes (problemID 33)

problemID 3 (**lbm4u-euler-fvm**) with particles of one spheroidal shape: every cell holds the concentration c and the
orientation tensor A = ⟨p p⟩ of the particle axes p (trace 1; I/3 random, p p fully aligned). The drag uses the mobility
averaged over the orientations and A evolves by the Folgar-Tucker equation. Transport, concentration files, coupling,
refinement and the free surface are those of problemID 3 - read that skill first. Code:
`BlockFreeSurfaceEulerNoSphereParticlesFVM.cpp`, `SpheroidOrientation.cpp` (shared with problemID 42). Equations:
`.claude/skills/lbm4u-particles/references/models.md` section 8 and `docs/math/matematicky_popis.pdf`.

## When to use it

* The orientation field is the result: short fibres or flakes in mould filling, alignment by shear and extension,
  orientation-dependent settling of many particles.
* Too many particles for problemID 22, or a smooth field per cell is wanted.
* Not for: individual trajectories or tumbling of single particles (problemID 22), several shapes in one run,
  quantitatively exact partly aligned states (closure error, see Accuracy).

## Parameters (`parameters.txt`)

All parameters of problemID 3 (`particleDiameter` = **volume-equivalent** diameter, `particleDensity`, `fluidGravity`,
...) plus:

| Parameter | Default | Meaning |
|---|---|---|
| `problemID` | – | `33` |
| `particleAspectRatio` | 1 | r = a/b of the spheroid; r > 1 prolate (fibres), r < 1 oblate (flakes), 1 = problemID 3 |
| `orientationDiffusion` | 0 | Folgar-Tucker C_I; 0 for dilute, ~0.001-0.01 for interacting fibres |
| `orientationClosure` | 0 | closure of A4:D, 0 = hybrid (recommended), other values = quadratic / linear |

Semi-axes: a = r b along p, b across, a b² = (d/2)³.

## Initial orientation `init/orientation<tag>.txt`

Per cell `A11 A22 A33 A12 A13 A23`, block file order (k fastest); missing file = random (A = I/3). `make_case.py`:
`"concentration": {..., "orientation": "random" | [px, py, pz] | [A11, A22, A33, A12, A13, A23]}`.

## Physics and numerics (summary)

* Mobility averaged exactly over the orientations (Stokes): `⟨M⟩ = I/KT⊥ + (1/KT∥ − 1/KT⊥) A`; drag
  `K = ⟨M⟩⁻¹ μ g(Re)`, implicit velocity update. Random prolate r = 3 settles ~10 % slower than the volume-equivalent
  sphere; an aligned inclined suspension glides sideways like one particle of problemID 22.
* Transport of q = c A with the limited particle flux of c (flux_q = flux_c · A_upwind), then A = q / trace q -
  orientation moves only with the particles.
* Folgar-Tucker (RK4 per liquid cell): `dA/dt = W A − A W + λ (D A + A D − 2 A4:D) + 2 C_I γ̇ (I − 3A)`,
  λ = (r² − 1)/(r² + 1); hybrid closure exact for aligned and random states; A renormalised and kept positive.
* Refinement: q is exchanged with c at the level interface; coarse ghosts get the c-weighted mean of A.

## Accuracy and checks

* The closure overpredicts alignment of **partly aligned** prolate particles (r ≈ 2-5): in shear and in mould filling
  ⟨p_x²⟩ is 0.05-0.2 higher than with Lagrange (problemID 22), e.g. 0.71-0.76 vs 0.60-0.67 in a gate for r = 3.
  Nearly aligned particles (r = 10 in a gate) agree (0.92-0.94). Eulerian clouds keep their alignment longer.
* r = 1 reproduces problemID 3 (`euler_nosphere_sphere_limit`); Jeffery orbits of an aligned suspension follow the
  analytical period (`euler_nosphere_jeffery_orbit`); a uniform A is kept by the transport.
* Check trace A = 1 and 1/3 ≤ `orientation_alignment` ≤ 1 in the output; values outside mean a bug.

## Output

`output_*.vtu`: the fields of problemID 3 plus `orientation_tensor` (9 components), `orientation_direction` (main
eigenvector) and `orientation_alignment` (largest eigenvalue, 1/3 random ... 1 aligned). ParaView: glyphs of
`orientation_direction` scaled by `orientation_alignment`; for statistics average A weighted by c over a region.

## Limits

One shape and size per run, one velocity per cell (no dispersion by different orientations), closure error for
partly aligned states, Jeffery (creeping flow) rotation, no particle stress in the liquid (no fibre-induced
viscosity anisotropy), inlets do not bring particles in.

Tests: `lbm_tests.exe euler_nosphere_` (excluding `_lbm` ones: orientation equations, sphere limit, isotropic and
aligned inclined settling, Jeffery orbit, rotary diffusion, uniform orientation, free surface) and
`refinement_euler_nosphere_conservation`.
