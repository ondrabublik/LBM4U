---
name: lbm4u-euler-lbm-spheroids
description: Set up, run, check and evaluate the Eulerian non-spherical particle model with LBM transport (problemID 42) of the cLBMnet / LBM4U free-surface solver - spheroid concentration carried by LBM populations with particleDiffusivity, orientation tensor A carried by the net link flux of the upwind cell, Folgar-Tucker equation with particleAspectRatio, orientationDiffusion and orientationClosure, init/orientation files, refined meshes, orientation output. Use it whenever the user works with problemID 42, wants a diffusive or packing-limited suspension of fibres or flakes with an orientation field, compares 33 and 42, or asks how the orientation is transported by the populations - also in Czech ("nekulové částice LBM", "orientační tenzor", "orientace vláken", "problemID 42"). Use with lbm4u-euler-lbm (transport and particleDiffusivity), lbm4u-euler-fvm-spheroids (orientation physics), lbm4u-particles and lbm4u-free-surface.
---

# Eulerian spheroids with LBM transport (problemID 42)

The concentration c is carried by the LBM populations of problemID 4 (**lbm4u-euler-lbm**: `particleDiffusivity`,
stability limit, packing) and every cell holds the orientation tensor A with the drag and the Folgar-Tucker equation
of problemID 33 (**lbm4u-euler-fvm-spheroids**: `particleAspectRatio`, `orientationDiffusion`, `orientationClosure`,
`init/orientation<tag>.txt`, output fields). Read those two skills for the parameters; this one adds what is specific
to the combination. Code: `BlockFreeSurfaceEulerNoSphereParticlesLBM.cpp`, `SpheroidOrientation.cpp`, the refinement
exchange in `Mesh.cpp` (`scalarCoalesce`, `scalarAccount`, step 44). Equations:
`.claude/skills/lbm4u-particles/references/models.md` section 9 and `docs/math/matematicky_popis.pdf`, section 7.5.

## When to use it

* As problemID 33 when the particle transport should be diffusive (dispersion, mixing) or packing-limited, or when the
  LBM transport is preferred for simplicity.
* Choose 33 vs 42 like 3 vs 4: 42 smears clouds and fronts somewhat more (`particleDiffusivity`) and their orientation with them;
  in advection-dominated flows (mould filling) both give nearly the same orientation field.

## Parameters

`problemID 42`, the parameters of problemID 4 (`particleDiameter` volume-equivalent, `particleDensity`,
`particleDiffusivity` (recommended `0.001 · h_L² / Δt_L`), `frontDiffusivity`, `overpackDiffusivity`) and of problemID 33 (`particleAspectRatio`,
`orientationDiffusion`, `orientationClosure`). Concentration and orientation files as for 33. `gravityRotationOmega`
works as in problemID 4.

## Orientation transport (specific to 42)

The net particle flux of each link (out = g_i leaving x, in = g_opp arriving from x + e_i) carries A of the upwind
cell; bounce-back links carry nothing:

```
S A(x, t+1) = S A(x, t) + Σ_links max(in − out, 0) (A(x + e_i, t) − A(x, t)),     S = Σ_i g_i(x, t+1)
```

* Σ_x S A is conserved like Σ g, A stays a convex combination (trace 1, positive), a uniform A is kept, and A mixes
  only with the particles (their diffusion plus the first-order upwind error) - no transport without a particle flux.
* On refined meshes the net fluxes of the exchange balance carry A of their upwind cell (exchangeOrient), applied in
  step 44 before the Folgar-Tucker step.
* Do **not** let every population carry A of its source cell: the gross population exchange is much larger than the
  net flux and A would diffuse ~100× faster than the particles (diffusivity ~1/(6θ) lattice units). This was the bug
  of the first version (commit a54cbdb, fixed in 28fed0c); it made orientation fields nearly random.

## Physics specific to the LBM velocity

The particle velocity includes the explicit convective term `(v·∇)v`; the liquid sees the fluid fraction 1 − cV,
so the liquid velocity varies through a cloud and the spheroids turn slightly in that gradient (absent in 33).

## Accuracy and checks

* Equidistant-mesh study (local only, BGK transport before the TRT change, D_lat = 0.005): sphere cloud 5.29 mm/s vs FVM 6.68 mm/s and Lagrange 7.17 mm/s (smearing,
  σ_z 0.34 vs 0.20 cm); mesh h/2 +16 %. Shape effect on the cloud velocity agrees with the other models within 1-2
  percentage points. Orientation at the end of the cloud fall: ⟨p_z²⟩ 0.53 (r = 3) vs 0.50 FVM, 0.42 Lagrange.
* Mould filling: same orientation structure as 33 and 22; slightly lower extremes than 33 (⟨p_x²⟩ behind the gate
  0.59 vs 0.68 for r = 3) because the particle diffusion also mixes the orientation.
* Particle count conserved to ~1e-8; Σ S A conserved likewise (`euler_nosphere_lbm_orientation_measure_is_conserved`).
* An orientation step in a liquid at rest must stay sharp (`euler_nosphere_lbm_orientation_needs_particle_flux`).

## Limits

Those of problemID 4 (smearing by D, small undershoots at sharp fronts, one size per run) and of problemID 33 (closure error for partly aligned
states, Jeffery rotation, no particle stress in the liquid); first-order upwind smearing of A at orientation fronts;
inlets do not bring particles in.

Tests: `lbm_tests.exe euler_nosphere_lbm` and `lbm_tests.exe refinement_euler_nosphere_lbm`.
