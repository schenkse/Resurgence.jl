# Trans-series

A trans-series

```math
F(g) \;=\; \sum_j e^{-S_j/g} \, g^{\beta_j} \, \sum_k a_{jk} \, g^k
```

is the natural object of resurgence theory: the perturbative answer (`j = 0`) sitting alongside instanton sectors (`j ≥ 1`) that perturbation theory can't see on its own.
This page builds one by hand, evaluates it numerically with [`resum_transseries`](@ref), and uses the action-additive arithmetic on top.

Prerequisite: [Lateral and median sums](lateral_sums.md).
The worked example below is the Stokes-completion of the non-Borel-summable driver from that page, so it helps to have seen the lateral/median story once first.

## Driver: the all-positive series `aₖ = k!`

```@example transseries
using Resurgence
using Printf

a = ComplexF64[ComplexF64(factorial(big(k))) for k in 0:24]
nothing # hide
```

`ComplexF64` rather than `Float64` because we'll combine the perturbative series with a sector carrying an `iπ` prefactor in a moment — keeping everything in one eltype lets the `TransSeries` constructor stay type-stable.

The bare resummations from [Lateral and median sums](lateral_sums.md):

```@example transseries
v_median  = borel_pade_median(a;  n = 10, m = 10, x = 1)
v_lateral = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = +1)
@printf("median      %.6f + %.6fi\n", real(v_median),  imag(v_median))
@printf("lateral L+  %.6f + %.6fi\n", real(v_lateral), imag(v_lateral))
```

The lateral sum picks up an imaginary part `≈ -1.1557 = -π/e` at `g = 1` that the median doesn't have.
Resurgence theory says that imaginary piece *is* an explicit one-instanton sector: action `S = 1` (the Borel-plane pole sits at `t = 1`), exponent `β = -1` (from the `1/g` in the residue at `t = 1/g`), and leading coefficient `-iπ`.
The sign of the `iπ` here is glued to the package's `side = +1` convention; using `side = -1` would flip both signs together — see the convention note in [Lateral and median sums](lateral_sums.md).

## Building the trans-series

```@example transseries
# 1-instanton sector: e^{-S/g} g^β · (-iπ + 0·g + 0·g² + …).
# Pad with zeros to keep the same length as the perturbative tail.
a_inst = ComplexF64[-im * π; zeros(24)]

ts = TransSeries([
    Sector(0.0,  0.0, a),       # perturbative sector
    Sector(1.0, -1.0, a_inst),  # one-instanton sector
])
length(ts.sectors)
```

`TransSeries` is just a `Vector{Sector{T}}`; each `Sector(S, β, a)` records one term in `Σⱼ e^{-Sⱼ/g} g^{βⱼ} (Σₖ aⱼₖ gᵏ)`.
The two sectors above auto-promote to a common element type, so `ts::TransSeries{ComplexF64}`.

## Numerical evaluation

```@example transseries
g = 1.0
v = resum_transseries(ts, g; method = BorelPadeMedian(10, 10))
@printf("trans-series at g = 1:  %.6f + %.6fi\n", real(v), imag(v))
@printf("lateral L+ (for ref):   %.6f + %.6fi\n", real(v_lateral), imag(v_lateral))
@printf("difference:             %.2e\n", abs(v - v_lateral))
```

[`resum_transseries`](@ref) walks the sector list, calls `resum(BorelPadeMedian(10, 10; x = g), sec.a)` on each perturbative tail (the method's evaluation point is rebound from `1` to `g`), multiplies by the prefactor `e^{-S/g} g^β`, and sums.
For the perturbative sector that gives the median; for the trivial instanton sector `[-iπ, 0, 0, …]` it gives `-iπ`; the prefactor `e^{-1/g} g^{-1}` at `g = 1` is `1/e`; total: `median - iπ/e ≈ L⁺`.

The `(S, β)` labels carry the full `g`-dependence — change `g` and the prefactor follows along automatically:

```@example transseries
g = 0.5
v_g  = resum_transseries(ts, g; method = BorelPadeMedian(10, 10))
v_lp = borel_pade_lateral(a; n = 10, m = 10, x = g, side = +1)
@printf("at g = 0.5: trans-series = %.6f + %.6fi\n", real(v_g),  imag(v_g))
@printf("            lateral L+   = %.6f + %.6fi\n", real(v_lp), imag(v_lp))
```

The `β = -1` matters here: the imaginary part now scales as `-2π·e^{-2}` rather than just `-π·e^{-2}`.

## Action-additive arithmetic

The interesting algebraic fact: multiplying two trans-series adds their actions.

```@example transseries
seed = TransSeries([Sector(1.0, 0.0, [1.0])])   # bare 1-instanton seed
[(s.S, s.β, s.a) for s in (seed * seed).sectors]
```

`seed * seed` produces a single sector at `S = 2, β = 0` with Cauchy product `[1] ⋆ [1] = [1]` as its perturbative tail — the standard two-instanton action.
Addition does the dual operation: sectors with the same `(S, β)` get merged, and their perturbative tails are summed (padded to the longer length).

## Dilute instanton gas via `transseries_exp`

`transseries_exp(seed; order)` builds the truncated Taylor series `1 + seed + seed²/2! + ⋯ + seed^order / order!`.
With a single-instanton sector as the seed, this is exactly the dilute-instanton-gas expansion:

```@example transseries
gas = transseries_exp(seed; order = 4)
[(s.S, s.a) for s in gas.sectors]
```

The `k`-instanton sector has action `k · S` and Cauchy-product coefficient `1/k!` — what you'd get from formally evaluating `exp(c · e^{-S/g})` and reading off the multi-instanton coefficients.

## Where this fits

- **You have an explicit instanton prediction** — from a saddle-point analysis, an `e^{-S/g}` correction in a textbook, the Stokes constant of an existing series, etc. — and want to add it to a bare perturbative series.
  Build a two-sector `TransSeries` and call `resum_transseries`.
- **You want to study how a non-perturbative correction interacts with the perturbative one across orders in `g`**.
  Keep both sectors as full series and let the action-additive arithmetic combine them.
- **You're working with a dilute-instanton-gas approximation** (e.g., a partition function of the form `Σₙ (c·e^{-S/g})ⁿ / n!`).
  [`transseries_exp`](@ref) is the direct construction.
- **You want to verify a hypothesised trans-series structure numerically**.
  Drop the `(S, β, A)` tuple from [`stokes_fit`](@ref) straight into a `Sector` constructor and check.

API reference: [`TransSeries`](@ref), [`Sector`](@ref), [`resum_transseries`](@ref), [`transseries_exp`](@ref).
