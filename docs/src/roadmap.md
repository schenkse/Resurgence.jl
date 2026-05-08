# Roadmap

Resurgence.jl ships a working classical-resummation toolkit: Shanks and
Richardson acceleration, Padé and Borel / Borel–Le Roy / (single-singularity)
conformal Borel–Padé, Meijer-G via the Slater-collapse onto `pFq`, optimal
truncation, and lateral / median / discontinuity Borel–Padé with Stokes /
large-order diagnostics on top.

The remaining gaps are mostly on the resurgence-specific side — first-class
trans-series, hyperasymptotic remainders, multi-singularity conformal maps,
order-dependent mapping — and a handful of standard sequence-acceleration
tools (Levin, Weniger, Wynn ε, Brezinski θ/ρ, Sidi S) that users expect to
find. This page is a working list of those items, sized and sketched enough
to start a focused implementation task.

## Status at a glance

| ID  | Technique                                | Section            | Size | Status  |
| --- | ---------------------------------------- | ------------------ | ---- | ------- |
| A1  | Wynn ε-algorithm                         | Acceleration       | S    | Planned |
| A2  | Levin u / t / v transforms               | Acceleration       | M    | Planned |
| A3  | Weniger δ-transformation                 | Acceleration       | S    | Planned |
| A4  | Brezinski θ-algorithm                    | Acceleration       | S    | Planned |
| A5  | Brezinski ρ-algorithm                    | Acceleration       | S    | Planned |
| A6  | Aitken–Steffensen for fixed points       | Acceleration       | S    | Planned |
| A7  | Sidi S-transformation (W-algorithm)      | Acceleration       | S    | Planned |
| A8  | Cesàro / Abel summation                  | Acceleration       | S    | Planned |
| B1  | Trans-series type and arithmetic         | Borel / resurgence | L    | Planned |
| B2  | Per-sector trans-series resummation      | Borel / resurgence | S    | Planned |
| B3  | Hyperasymptotic remainders               | Borel / resurgence | M    | Planned |
| B4  | Multi-singularity conformal map          | Borel / resurgence | M    | Planned |
| B5  | Order-dependent mapping (ODM)            | Borel / resurgence | M    | Planned |

Sizes: **S** ≤ ~50 LOC, **M** ≤ ~300 LOC, **L** = needs design.

Per-item entries below use the same fixed layout — **Goal**, **Why now**,
**Sketch**, **Reuses**, **API tag** — so they can be skimmed in any order.

## A. Sequence acceleration (general-purpose utilities)

### A1 — Wynn ε-algorithm  · *Small*

**Goal.** Tableau form of iterated Shanks via
`ε^{(n)}_{k+1} = ε^{(n+1)}_{k-1} + 1/(ε^{(n+1)}_k − ε^{(n)}_k)`.

**Why now.** Cleaner and numerically more stable than the partial-sum
recursion in [src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl)
`_shanks_iter` for high `depth`. Standard in the literature; users expect it.

**Sketch.** New `wynn_eps(a, n; depth)` in
[src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl). Operates on the
tableau directly; even-column entries are the meaningful approximants.

**Reuses.** Same partial-sum prelude as `_shanks_iter`.
**API tag.** `WynnEps <: AbstractResummation` in [src/api.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/api.jl).

### A2 — Levin transforms (u, t, v)  · *Medium*

**Goal.** `L_k^{(n)} = Σⱼ (−1)ʲ C(k,j) (n+j+β)^{k−1} Aₙ₊ⱼ / ωₙ₊ⱼ
            / Σⱼ (−1)ʲ C(k,j) (n+j+β)^{k−1} / ωₙ₊ⱼ`
with three remainder estimators ωₙ: `aₙ` (t), `(n+β) aₙ` (u),
`aₙ aₙ₊₁/(aₙ₊₁−aₙ)` (v).

**Why now.** Frequently outperform Shanks on monotone divergent / asymptotic
series — exactly the regime resurgence cares about. Not currently available.

**Sketch.** New `src/levin.jl` exposing `levin(a, n; variant=:u, β=1)`.
Pure-Julia, element-type generic.

**Reuses.** `inv_factorials` from [src/utils.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/utils.jl) for binomial
coefficients (or compute `C(k,j)` iteratively).
**API tag.** `Levin(n; variant, β)`.

### A3 — Weniger δ-transformation  · *Small* (after A2)

**Goal.** Levin-style transform with ωₙ = `(β+n)(β+n+1) … aₙ` factorial
weighting; ratio of two finite sums.

**Why now.** Often the single best black-box accelerator for factorially
divergent series — the textbook Stieltjes-series accelerator.

**Sketch.** Same module as Levin (`src/levin.jl`); near-trivial once A2 lands.

**Reuses.** A2 infrastructure.
**API tag.** `Weniger(n; β)`.

### A4 — Brezinski θ-algorithm  · *Small*

**Goal.** ε-algorithm variant that avoids the breakdown of Wynn's ε on
irregular sequences. Recursion mixes 4 tableau entries.

**Why now.** Useful when Wynn ε divides by zero; complements A1.

**Sketch.** New `theta_brezinski(a, n; depth)` in
[src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl).

**Reuses.** Same partial-sum prelude as A1.
**API tag.** `BrezinskiTheta(n; depth)`.

### A5 — Brezinski ρ-algorithm  · *Small*

**Goal.** Like ε but with `(k+1)` numerator, tuned for logarithmically
convergent sequences.

**Why now.** Companion to Richardson when partial sums decay like `1/n^p`
with non-integer p.

**Sketch.** Sibling function to A4 in
[src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl).

**Reuses.** Same scaffolding as A4.
**API tag.** `BrezinskiRho(n; depth)`.

### A6 — Aitken–Steffensen for fixed points  · *Small*

**Goal.** `x_{k+1} = g(x_k) − (g(x_k)−x_k)² / (g(g(x_k)) − 2g(x_k) + x_k)`.

**Why now.** Useful when a series can be written as a fixed-point iteration;
complements Shanks on the sum side.

**Sketch.** New `aitken_steffensen(g, x0; tol, maxiter)` in
[src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl).

**Reuses.** Standalone — no shared scaffolding.
**API tag.** None (function-style entry point only).

### A7 — Sidi S-transformation (W-algorithm)  · *Small*

**Goal.** Levin-style transformation for oscillatory integrands; W-algorithm
gives an efficient recursion.

**Why now.** Useful when partial sums oscillate and Levin u/t struggle.

**Sketch.** Append to `src/levin.jl`.

**Reuses.** A2 infrastructure.
**API tag.** `SidiS(n; …)`.

### A8 — Cesàro / Abel summation  · *Small*

**Goal.** `cesaro(a, n)` and `abel(a; x)` (radial limit
`lim_{x→1⁻} Σ aₖ xᵏ`).

**Why now.** Mainly for completeness and as baselines; very cheap to
implement and useful for sanity checks against more sophisticated methods.

**Sketch.** Append to
[src/series_acceleration.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/series_acceleration.jl).

**Reuses.** Existing partial-sum machinery.
**API tag.** `Cesaro`, `Abel(x)`.

## B. Borel side / resurgence-specific (core)

### B1 — Trans-series type and arithmetic  · *Large*

**Goal.** A `TransSeries{T}` representing
`Σⱼ e^{−Sⱼ/g} g^{βⱼ} (Σₖ aⱼₖ gᵏ)`. Stores a vector of `Sector{T}` records
`(S, β, a::Vector{T})`. Arithmetic: addition (concat sectors, merge equal
`(S,β)`), scalar multiplication, multiplication (Cauchy product per sector +
cross-sector convolution), and a `transseries_exp` for the canonical
`exp(small)` build-up.

**Why now.** Resurgence is fundamentally trans-series, not series. Even
without alien calculus, having a first-class type lets users define and
resum non-perturbative completions (instanton + perturbative + …).

**Sketch.** New `src/transseries.jl`. Keep ops minimal (add, scale, mul);
avoid premature alien-derivative scope.

**Reuses.** Nothing internally, but `resum(::TransSeries, g)` (B2) builds on
the existing per-sector resummers.
**API tag.** `TransSeries`, `Sector`, `+`, `*`, `transseries_exp`.

### B2 — Per-sector trans-series resummation  · *Small* (after B1)

**Goal.** `resum_transseries(ts::TransSeries, g; method=BorelPade(...))`
resums each sector's perturbative series with `method`, multiplies by
`e^{−Sⱼ/g} g^{βⱼ}`, and sums.

**Why now.** The point of the type. Enables one-line "instanton-completed"
answers.

**Sketch.** Appended to `src/transseries.jl`. Loop over sectors, dispatch
through the existing `resum(::AbstractResummation, ...)`.

**Reuses.** [src/api.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/api.jl) entire dispatch surface.
**API tag.** Standalone function; no new struct.

### B3 — Hyperasymptotic remainders (Berry–Howls / Costin)  · *Medium*

**Goal.** Extend [src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl)
`superasymptotic_remainder` by adding terminant-function corrections that
incorporate the leading instanton contribution: error scales like
`e^{−2|S|/g}` instead of `e^{−|S|/g}`. Iteratively, hyperasymptotic level
`k` gives error `e^{−(k+1)|S|/g}`.

**Why now.** The natural next step beyond superasymptotic, and the
*quantitative* face of resurgence: precision improves geometrically by
re-summing the tail of the original series as a new asymptotic series
weighted by terminants.

**Sketch.** New `hyperasymptotic(a, x; level=1, action=nothing)` in
[src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl). Terminant `T_p(σ)` via the
incomplete-gamma representation
`T_p(σ) = e^{iπp} Γ(p) Γ(1−p, σ) / (2πi)`. If `action` not provided, derive
it from the existing Stokes-fit utilities.

**Reuses.** [src/truncation.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/truncation.jl), `SpecialFunctions.gamma_inc`,
existing Stokes diagnostics for `S` extraction.
**API tag.** `Hyperasymptotic(x; level, action)`.

### B4 — Multi-singularity conformal map  · *Medium*

**Goal.** Generalise [src/conformal.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/conformal.jl) `conformal_map`
from a single negative-real singularity at `t = −sing` to a finite set of
complex singularities `{tⱼ}`. Two flavours: (i) symmetric pair `±i·sing`
(common in problems with complex-conjugate singularity pairs — quartic
anharmonic oscillator, etc.); (ii) generic uniformizing map for arbitrary
`{tⱼ}` via Schwarz–Christoffel-style construction.

**Why now.** Many physics problems have Borel singularities at `±i`
(PT-symmetric), at multiple instanton actions, or as a branch cut. The
current single-pole assumption silently mis-treats these.

**Sketch.**

- `conformal_map_pair(t; a)` — `w(t) = t / √(1 + (t/a)²)` style map for
  `±i·a`.
- `conformal_map_set(t, sings)` — generic; harder, defer to second pass.
- Add `conformal_borel_pade_pair(a; n, m, x, sing, kwargs...)`.

**Reuses.** [src/conformal.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/conformal.jl) `conformal_reseries` (the
re-expansion bookkeeping is reusable; only the map changes).
**API tag.** `ConformalBorelPadePair(...)`.

### B5 — Order-dependent mapping (ODM) / variational Bender–Boettcher  · *Medium*

**Goal.** Treat the conformal-map exponent (or Le Roy `b`) as a variational
parameter and choose it order-by-order to minimise sensitivity (`d/db = 0`).

**Why now.** Standard physics trick (Bender–Boettcher 1969 onwards) for
accelerating Borel–Le Roy and conformal Borel.

**Sketch.** `borel_leroy_pade_odm(a; n, m, x, kwargs...)` — sweep `b` over a
grid, locate the stationary point of the result vs. `b`, return the value
there. Same construction for the conformal exponent.

**Reuses.** [src/borel_pade.jl](https://github.com/schenkse/Resurgence.jl/blob/main/src/borel_pade.jl) `borel_leroy_pade` and
`conformal_borel_pade`.
**API tag.** `BorelLeroyPadeODM(...)`, `ConformalBorelPadeODM(...)`.
