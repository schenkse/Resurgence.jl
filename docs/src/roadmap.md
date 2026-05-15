# Roadmap

Resurgence.jl ships a working classical-resummation toolkit: Shanks, Wynn-ε, Richardson, Cesàro/Abel, Levin u/t/v, Weniger δ, Brezinski θ/ρ, Aitken–Steffensen, and Sidi S sequence accelerators; Padé and Borel/Borel–Le Roy/(single-singularity) conformal Borel–Padé; Meijer-G via the Slater-collapse onto `pFq`; optimal truncation; lateral/median/discontinuity Borel–Padé with Stokes/large-order diagnostics on top; and first-class trans-series (sectors, arithmetic, and per-sector resummation through the unified `resum` surface).

The remaining gaps are on the resurgence-specific side: hyperasymptotic remainders, multi-singularity conformal maps, and order-dependent mapping.
This page is a working list of those items, sized and sketched enough to start a focused implementation task.

## Status at a glance

| ID  | Technique                                | Section            | Size | Status  |
| --- | ---------------------------------------- | ------------------ | ---- | ------- |
| B1  | Hyperasymptotic remainders               | Borel/resurgence | M    | Planned |
| B2  | Multi-singularity conformal map          | Borel/resurgence | M    | Planned |
| B3  | Order-dependent mapping (ODM)            | Borel/resurgence | M    | Planned |

Sizes: **S** ≤ ~50 LOC, **M** ≤ ~300 LOC, **L** = needs design.

Per-item entries below use the same fixed layout — **Goal**, **Why now**, **Sketch**, **Reuses**, **API tag** — so they can be skimmed in any order.

## B. Borel side/resurgence-specific (core)

### B1 — Hyperasymptotic remainders (Berry–Howls/Costin)  · *Medium*

**Goal.** Extend [src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl) `superasymptotic_remainder` by adding terminant-function corrections that incorporate the leading instanton contribution: error scales like `e^{−2|S|/g}` instead of `e^{−|S|/g}`.
Iteratively, hyperasymptotic level `k` gives error `e^{−(k+1)|S|/g}`.

**Why now.** The natural next step beyond superasymptotic, and the *quantitative* face of resurgence: precision improves geometrically by re-summing the tail of the original series as a new asymptotic series weighted by terminants.

**Sketch.** New `hyperasymptotic(a, x; level=1, action=nothing)` in [src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl).
Terminant `T_p(σ)` via the incomplete-gamma representation `T_p(σ) = e^{iπp} Γ(p) Γ(1−p, σ)/(2πi)`.
If `action` not provided, derive it from the existing Stokes-fit utilities.

**Reuses.** [src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl), `SpecialFunctions.gamma_inc`, existing Stokes diagnostics for `S` extraction.
**API tag.** `Hyperasymptotic(x; level, action)`.

### B2 — Multi-singularity conformal map  · *Medium*

**Goal.** Generalise [src/conformal.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/conformal.jl) `conformal_map` from a single negative-real singularity at `t = −sing` to a finite set of complex singularities `{tⱼ}`.
Two flavours: (i) symmetric pair `±i·sing` (common in problems with complex-conjugate singularity pairs — quartic anharmonic oscillator, etc.); (ii) generic uniformizing map for arbitrary `{tⱼ}` via Schwarz–Christoffel-style construction.

**Why now.** Many physics problems have Borel singularities at `±i` (PT-symmetric), at multiple instanton actions, or as a branch cut.
The current single-pole assumption silently mis-treats these.

**Sketch.**

- `conformal_map_pair(t; a)` — `w(t) = t/√(1 + (t/a)²)` style map for `±i·a`.
- `conformal_map_set(t, sings)` — generic; harder, defer to second pass.
- Add `conformal_borel_pade_pair(a; n, m, x, sing, kwargs...)`.

**Reuses.** [src/conformal.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/conformal.jl) `conformal_reseries` (the re-expansion bookkeeping is reusable; only the map changes).
**API tag.** `ConformalBorelPadePair(...)`.

### B3 — Order-dependent mapping (ODM)/variational Bender–Boettcher  · *Medium*

**Goal.** Treat the conformal-map exponent (or Le Roy `b`) as a variational parameter and choose it order-by-order to minimise sensitivity (`d/db = 0`).

**Why now.** Standard physics trick (Bender–Boettcher 1969 onwards) for accelerating Borel–Le Roy and conformal Borel.

**Sketch.** `borel_leroy_pade_odm(a; n, m, x, kwargs...)` — sweep `b` over a grid, locate the stationary point of the result vs. `b`, return the value there.
Same construction for the conformal exponent.

**Reuses.** [src/borel_pade.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/borel_pade.jl) `borel_leroy_pade` and `conformal_borel_pade`.
**API tag.** `BorelLeroyPadeODM(...)`, `ConformalBorelPadeODM(...)`.
