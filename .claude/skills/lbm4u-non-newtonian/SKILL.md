---
name: lbm4u-non-newtonian
description: Set up and check free-surface simulations of non-Newtonian liquids (power law, Bingham / Herschel-Bulkley, Cross, polynomial, modified power law) with the cLBMnet / LBM4U solver - choosing viscosityModel and fluidModelParameters, converting rheometer data (dynamic viscosity, consistency K, flow index n, yield stress) to the solver's kinematic form, and making sure the viscosity range fits between tauMin and tauMax. Use it whenever the user mentions shear-thinning or shear-thickening fluids, pastes, slurries, concrete, polymer melts, blood, ceramic or food masses, yield stress, rheology or viscosity curves in a simulation of this repository, also in Czech ("nenewtonovská kapalina", "reologie", "mez kluzu", "mocninový model", "pseudoplastická"). Use together with lbm4u-free-surface, which covers the case files and running.
---

# Non-Newtonian liquids with free surface (cLBMnet / LBM4U)

The liquid viscosity is a function of the local shear rate. Each step, every liquid cell computes the shear rate γ̇
from the non-equilibrium part of its populations (local, no neighbours), evaluates the viscosity model, and sets its
own relaxation time `τ = 1/2 + 3 ν(γ̇) / Cvis`, clipped to `[tauMin, tauMax]`. It works with BGK, TRT and MRT, with
LES, on refined meshes and with all particle models (problemID 1-4). Case setup, units and running are in the
**lbm4u-free-surface** skill - read it first; this skill adds what changes for the rheology.

## The models

Set `viscosityModel` and `fluidModelParameters = b0, b1, ...` in `parameters.txt`. The model gives the
**kinematic** viscosity ν [m²/s] from γ̇ [1/s]:

| `viscosityModel` | ν(γ̇) | Typical use |
|---|---|---|
| 1 power law | `b0 · γ̇^(b1−1)` | `b0 = K/ρ`, `b1 = n` |
| 2 Bingham | `b0 + b1 · γ̇^(b2−1)` | Bingham: `b0 = μp/ρ`, `b1 = τ0/ρ`, `b2 = 0` |
| 3 Cross | `b0 + (b1 − b0) / (1 + b2 · γ̇^b3)` | `b0 = ν∞`, `b1 = ν0`, `b2 = λ^m`, `b3 = m` |
| 4 polynomial | `b0 + b1 γ̇ + b2 γ̇² + b3 γ̇³` | fitted data in a limited range |
| 5 modified power law | `b0 · γ̇^(b1−1 + b2 ln γ̇)` | γ̇ is limited from below at the minimum of the curve |
| 6 modified power law 2 | `b0 · γ̇^(b1−1 + b2 ln γ̇ + b3 ln² γ̇)` | as 5 with a quadratic term |

Converting rheometer data: the solver wants kinematic quantities, so divide every dynamic quantity (Pa·s, Pa·sⁿ,
Pa) by the density. Herschel-Bulkley `μ = τ0/γ̇ + K γ̇^(n−1)` has two power terms and no model matches it exactly.
Either use Bingham (model 2) with the plastic viscosity `K γ̇c^(n−1)` at the characteristic shear rate, or fit
model 5, 6 or 4 to the H-B curve in the γ̇ range of the flow - and tell the user which approximation you made and
how large its error is in that range.

`vKin` is still required: it sets the base τ printed at the start and used for cells that have no shear rate yet.
Set it to the viscosity at the characteristic shear rate `γ̇c ≈ U/L`.

## The key constraint: the representable viscosity range

The lattice viscosity is `ν / Cvis` with `Cvis = h² / Δt = h · Cu` (level 0, halves on every finer level). Because
τ is clipped, only viscosities in

```
ν_min = (tauMin − 1/2)/3 · Cvis     ν_max = (tauMax − 1/2)/3 · Cvis      (defaults tauMin 0.5001, tauMax 1.5)
```

are represented; outside of it the liquid silently has the clipped value. This matters most for shear-thinning
liquids and yield-stress liquids, whose ν → ∞ at γ̇ → 0: the unyielded or slow regions get `ν_max`, which acts as a
regularisation. Before running:

1. Estimate the γ̇ range of the flow (`U/L` for the bulk, `U/(few cells)` near walls and the front, ~0 in dead
   zones) and evaluate ν(γ̇) of the model over it.
2. Run `python .claude/skills/lbm4u-free-surface/scripts/make_case.py case.json --check`: it prints
   `ν_min .. ν_max` for every refinement level.
3. If the relevant viscosities exceed `ν_max`: raise `uMax` (Cvis grows with Cu, at the cost of more steps), use a
   coarser h, or raise `tauMax` (TRT handles τ up to ~3-5 with less accurate walls). Finer refinement levels make it
   worse: `τ_L − 1/2 = 2^L (τ_0 − 1/2)`.
4. If viscosities are below `ν_min` (τ → 0.5): add `smagorinsky = 0.1-0.17` and `tauMinusMax = 1` as for water.
5. For a yield-stress liquid, report what `ν_max` corresponds to: the effective "plug" viscosity. A plug that
   should not move will creep with velocity ~ `τ0 L / (ρ ν_max)`; if that matters, push `ν_max` up.

At start-up the solver prints `tau(gamma)` for γ̇ = 0, 0.1, 1, 10, 1e4 and the τ limits - read these lines and check
that the values in the relevant range are inside the limits.

## Checking the results

* `output_*.vtu` contains `shearRate` [1/s] and `viscosity` [m²/s] in liquid cells. Look at their range: cells sitting
  at ν_max or ν_min show where the clipping is active.
* Compare with a simple analytic case before trusting a new model setup: a gravity-driven film or channel
  (`tests/test_free_surface.cpp`, the power-law test checks the shear rate against the analytic profile to 3 %).
* A Newtonian run with `vKin` = ν(γ̇c) is a useful reference: differences should be explainable by the rheology.

## Example

Ceramic slurry, ρ = 1800 kg/m³, power law K = 2.5 Pa·sⁿ, n = 0.45, filling at U ≈ 0.5 m/s, h = 1 mm:

```
viscosityModel = 1
fluidModelParameters = 1.389e-3, 0.45      # K/ρ, n
vKin = 2.4e-4                              # ν at γ̇c = U/L = 0.5/0.02 = 25 1/s
uMax = 1                                   # Cu = 10 m/s, Cvis = h·Cu = 0.01 m²/s
tauMax = 3                                 # ν_max = 2.5/3 · 0.01 = 8.3e-3 m²/s  → γ̇ ≥ ~0.04 1/s represented
collisionModel = 2
```
