---
name: lbm4u-lagrange-ellipsoids
description: Set up, run, check and evaluate the Lagrangian non-spherical particle model (problemID 22) of the cLBMnet / LBM4U free-surface solver - rigid rotating ellipsoids and spheroids (short fibres, flakes, grains) with semi-axes a1 a2 a3, quaternion orientation, exact Stokes resistance tensors, Jeffery rotation in the velocity gradient, drag factors f1-f3, wall contact with the support function, init/particles.txt with shape and orientation, tensor glyphs in ParaView. Use it whenever the user works with problemID 22, tracks individual non-spherical particles, their orientation, tumbling, Jeffery orbits, gliding or sideways drift - also in Czech ("nekulové částice", "elipsoidy", "rotující částice", "vlákna", "destičky", "Jefferyho orbity", "problemID 22"). Use with lbm4u-particles (choice between models), lbm4u-lagrange-spheres (common Lagrangian machinery) and lbm4u-free-surface (case files, running).
---

# Lagrangian ellipsoids (problemID 22)

Rigid ellipsoids with semi-axes a1, a2, a3 along their body axes, a position, velocity, orientation (unit quaternion)
and angular velocity. Everything that is not about shape - interpolation, two-way coupling, walls, refinement,
sub-steps, the particle file workflow - is inherited from problemID 2 (**lbm4u-lagrange-spheres**); read that skill for
those parameters. Code: `LagrangeNoSphereParticles.cpp`, `BlockFreeSurfaceLagrangeNoSphereParticles.cpp`. Equations:
`.claude/skills/lbm4u-particles/references/models.md` section 7, `docs/math/matematicky_popis.pdf`.

## When to use it

* Elongated or flat particles shorter than a cell (short fibres, flakes, grains) when orientation matters:
  orientation-dependent drag (gliding, sideways drift), alignment in shear or extension, tumbling.
* Reference solution for the Eulerian spheroid models (problemID 33, 42), especially for partly aligned states where
  their closure is inaccurate.
* Not for: fibres longer than a cell (the largest semi-axis must stay below h), particles at finite Re where inertial
  torques matter (falling discs turning broadside are not modelled), very many particles (use 33 or 42).

## Particle file `init/particles.txt`

```
N
x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]        (SI units)
```

* `a1..a3` semi-axes [m] along the body axes; the quaternion rotates body axes to the lab frame
  (`1 0 0 0` = axes along x, y, z). A prolate spheroid of aspect ratio r: a1 = r b, a2 = a3 = b.
* `f1..f3` optional factors multiplying the drag along the body axes (default 1) - for calibration to measured drag.
* `make_case.py`: `"particles": {..., "shape": [a1, a2, a3], "orientation": "random" | [qw, qx, qy, qz],
  "drag_factors": [f1, f2, f3]}` instead of `"d"`.

Parameters in `parameters.txt` are those of problemID 2 (`nIntIter`, `couplingWayNumber`, `particleRestitution`,
`maxParticlesInCell`, `fluidGravity`, `gravityRotationOmega`) with `problemID 22`.

## Physics and numerics (summary)

* Translation: exact Stokes resistance KT_i along each body axis (elliptic integrals R_F, R_D), multiplied by the
  Schiller-Naumann correction g(Re) with the volume-equivalent diameter; implicit update with the full tensor
  `K = R diag(f_i μ KT_i g) Rᵀ`. A 2:1 prolate spheroid has ~13 % less drag along its axis than across; a 45° tilt
  glides ~4° sideways.
* Rotation: Euler equations in the body frame with the Jeffery torque `μ KR_i (ω^J_i − ω_i)`, ω^J from the
  interpolated vorticity and strain rate; quaternion update. Small particles follow Jeffery orbits with period
  `2π/γ̇ (r + 1/r)`.
* Velocity gradient per cell from the lattice neighbours (mirror value behind walls), interpolated with the kernel.
* Wall contact with the support function `sqrt(nᵀ R diag(a²) Rᵀ n)` - long particles stay further from walls.
* The particle torque is not returned to the liquid (only the force).

## Accuracy and checks

* Keep the **largest** semi-axis below the cell size (point-particle limit); d/h ≲ 0.3 as for spheres with the
  volume-equivalent d.
* r = 1 reproduces problemID 2 (test `nosphere_sphere_limit_equals_problem2`).
* Check Jeffery orbits in simple shear (test `nosphere_jeffery_period`) when changing the code.
* Studies (local only): falling cloud on an equidistant mesh - shape changes the cloud velocity by at most 2 %
  (vertical prolate particles −4 to −6 %), Euler models agree within 1-2 percentage points; in mould filling the
  Eulerian models overpredict alignment of partly aligned r = 3 particles (⟨p_x²⟩ 0.71-0.76 vs 0.60-0.67), r = 10
  agrees (0.92-0.94). The Lagrangian result has statistical noise ∝ 1/sqrt(particles per cell).

## Output

`output_particles<time>us.vtk` additionally holds `TENSORS shape` (ParaView: Tensor Glyph with a sphere source shows
the ellipsoids), `axis1..3`, `angularVelocity`. For orientation statistics compute ⟨p p⟩ of the main axis per region
(or per cell for comparison with the Eulerian `orientation_tensor`).

## Limits

Creeping-flow (Jeffery) rotation only, no inertial lift or torque, no rotational Brownian motion, no particle stress
or torque in the liquid, point particles (shape below the cell), contacts only by the repulsion of problemID 2.

Tests: `lbm_tests.exe nosphere` (resistance of ellipsoids, Jeffery period, sphere limit, inclined settling, rotation
in shear, wall contact).
