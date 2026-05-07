# Resurgence.jl

Resummation techniques for divergent / asymptotic series, in pure Julia.

The package targets resurgence and perturbative-physics workflows: take a
formal power series with factorially-growing coefficients, recover a
Borel-summable result. All methods are generic over the input element type —
`Float64`, `BigFloat`, and complex variants are all first-class.

[![Julia ≥ 1.10](https://img.shields.io/badge/Julia-≥1.10-9558B2?logo=julia)](https://julialang.org)
[![Polynomials v4](https://img.shields.io/badge/Polynomials-v4-blue)](https://github.com/JuliaMath/Polynomials.jl)
[![PolynomialRoots v1](https://img.shields.io/badge/PolynomialRoots-v1-blue)](https://github.com/giordano/PolynomialRoots.jl)
[![QuadGK v2](https://img.shields.io/badge/QuadGK-v2-blue)](https://github.com/JuliaMath/QuadGK.jl)
[![SpecialFunctions v2](https://img.shields.io/badge/SpecialFunctions-v2-blue)](https://github.com/JuliaMath/SpecialFunctions.jl)
[![GenericLinearAlgebra v0.4](https://img.shields.io/badge/GenericLinearAlgebra-v0.4-blue)](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl)

> **Built with LLMs:** This project was developed with the help of AI coding tools, primarily [Claude Code](https://claude.ai/code) by Anthropic. All code has been reviewed and is maintained by the author.

## Key features

- Element-type generic: `Float64`, `BigFloat`, and complex variants are all first-class. (Generic SVD for rank-deficient `BigFloat` Padé fits is provided by [GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl).)
- Multiple resummation families: Shanks, Richardson, Padé (linear, continued-fraction qd, and Hermite/quadratic for algebraic branch points), Borel / Borel–Le Roy / Borel–Padé / conformal Borel–Padé, Meijer-G, and optimal-truncation.
- Both a functional API and a unified `resum(::AbstractResummation, a)` dispatch layer.

## References

- C. M. Bender & S. A. Orszag, *Advanced Mathematical Methods for Scientists and Engineers* (1978) — asymptotic series and optimal truncation.
- O. Costin, *Asymptotics and Borel Summability* (2008) — modern treatment of Borel summation and resurgence.

## Setup

> Resurgence.jl is not yet registered in Julia's General registry, so
> `Pkg.add("Resurgence")` will not work — install from GitHub for now.

### Install from GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/schenkse/Resurgence.jl")
using Resurgence
```

### Clone for development

```bash
git clone https://github.com/schenkse/Resurgence.jl
cd Resurgence.jl
```

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
using Resurgence
```

## Methods

| Method | Functional API | Tag (for `resum`) |
|---|---|---|
| Shanks transformation (with iterated depth) | `shanks(a, n; depth)` | `Shanks(n; depth)` |
| Richardson extrapolation (with iterated depth) | `richardson(a, n; depth)` | `Richardson(n; depth)` |
| Padé approximant `[n/m]` | `pade(a, n, m)`, `pade_value(a, n, m, x)` | `Pade(n, m; x)` |
| Padé via qd continued fraction (`n == m` or `n + 1 == m`) | `pade_cf(a, n, m)`, `pade_cf_value(a, n, m, x)` | `PadeCF(n, m; x)` |
| Hermite/quadratic Padé (algebraic branch points) | `hermite_pade(a, n, m, l)`, `hermite_pade_value(a, n, m, l, x; branch)` | `HermitePade(n, m, l; x, branch)` |
| Borel transform | `borel_transform(a)` | — |
| Borel–Le Roy transform | `borel_leroy_transform(a, b)` | — |
| Borel–Padé resummation | `borel_pade(a; n, m, x, ...)` | `BorelPade(n, m; x, ...)` |
| Lateral Borel–Padé sum (`side = ±1`) | `borel_pade_lateral(a; n, m, x, side, ...)` | `BorelPadeLateral(n, m; x, side, ...)` |
| Median Borel–Padé sum `(L⁺+L⁻)/2` | `borel_pade_median(a; n, m, x, ...)` | `BorelPadeMedian(n, m; x, ...)` |
| Borel–Padé Stokes discontinuity `(L⁺−L⁻)/(2i)` | `borel_pade_discontinuity(a; n, m, x, ...)` | — |
| Borel–Le Roy–Padé resummation | `borel_leroy_pade(a; n, m, b, x, ...)` | `BorelLeRoyPade(n, m; b, x, ...)` |
| Lateral / median / discontinuity Borel–Le Roy–Padé | `borel_leroy_pade_lateral`, `_median`, `_discontinuity` | — |
| Conformal-Borel–Padé resummation | `conformal_borel_pade(a; n, m, x, sing, ...)` | `ConformalBorelPade(n, m; x, sing, ...)` |
| Meijer-G resummation (Slater-collapsed onto pFq) | `borel_meijerg(a; n, x, ...)` | `MeijerG(n; x, ...)` |
| Optimal-truncation / superasymptotics | `optimal_truncation(a)`, `superasymptotic_remainder(a)` | — |
| Stokes / large-order extraction (`S`, `β`, `A`, subleading `cⱼ`) | `stokes_action`, `stokes_exponent`, `stokes_constant`, `stokes_fit` | — |

`quadgk` keyword arguments (`rtol`, `atol`, `order`) flow through to the
Laplace integration step in every Borel-based method.

## Quick start: Stieltjes / Euler's series

The asymptotic series `S(z) = Σ (-1)^k k! z^k` is divergent for every nonzero
`z`, but is Borel summable to `(1/z) e^{1/z} E₁(1/z)`. At `z = 1` the exact
value is `e · E₁(1) ≈ 0.5963473623`.

```julia
using Resurgence

# 25 coefficients  a[k+1] = (-1)^k k!
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]

borel_pade(a; n = 10, m = 10, x = 1)             # ≈ 0.59634736
conformal_borel_pade(a; n = 10, m = 10, x = 1)   # ≈ 0.59634736
resum(BorelPade(10, 10), a)                       # same, via the unified API

# Higher precision: BigFloat propagates end-to-end
ab = BigFloat[(-1)^k * factorial(big(k)) for k in 0:40]
borel_pade(ab; n = 20, m = 20, x = BigFloat(1))  # ≈ 0.59634736…
```

See [`examples/stieltjes.jl`](examples/stieltjes.jl) for a side-by-side
comparison of all methods.

## Non-Borel-summable series and lateral sums

When the Borel transform has a singularity on the *positive* real axis (i.e.
the series is not Borel-summable in the conventional sense), the Laplace
integral is ill-defined; deforming the contour above (`L⁺`) or below (`L⁻`)
the singularity gives two distinct lateral sums, each complex-valued. Their
average is the ambiguity-free median sum, and their difference exposes the
Stokes discontinuity.

Driver: `aₖ = k!`. Borel transform `1/(1−t)` has a pole at `t = 1`. The
median sum recovers the Cauchy principal value `Ei(1)/e ≈ 0.6972`, and the
Stokes discontinuity has magnitude `π/e ≈ 1.1557`.

```julia
using Resurgence

a = Float64[Float64(factorial(big(k))) for k in 0:24]   # a_k = k!

borel_pade_median(a; n = 10, m = 10, x = 1)         # ≈  0.69717 + 0im
borel_pade_discontinuity(a; n = 10, m = 10, x = 1)  # ≈ -1.15573 + 0im

# The two lateral sums are complex conjugates of each other for real `a`.
borel_pade_lateral(a; n = 10, m = 10, x = 1, side = +1)   # ≈ 0.69717 - 1.15573im
borel_pade_lateral(a; n = 10, m = 10, x = 1, side = -1)   # ≈ 0.69717 + 1.15573im
```

Convention: `side = +1` shifts the offending Padé poles by `+iε` (upper
half-plane), so the integration contour passes *below* them — matching the
lateral sum traditionally written `S_θ⁻` in the resurgence literature.
The discontinuity is computed as `(L_{+1} − L_{−1}) / (2i) = (S_θ⁻ − S_θ⁺) / (2i)`,
which has magnitude equal to the textbook `(S_θ⁺ − S_θ⁻) / (2i)` but the
opposite sign. Negate if you want the textbook convention.

## Large-order / Stokes diagnostics

A divergent perturbative series with leading large-order behaviour

    a[k+1] ∼ A · Γ(k + β) / S^{k+β} · (1 + c₁/k + c₂/k² + …)

encodes its instanton action `S`, exponent `β`, prefactor `A`, and
subleading `cⱼ` directly in the coefficient ratios. `stokes_fit` reads them
off via Richardson-accelerated tail extrapolation; `stokes_action`,
`stokes_exponent`, and `stokes_constant` are the underlying single-quantity
functions.

```julia
using Resurgence

# a_k = (k+1)!  ⇒  exact (S, β, A) = (1, 2, 1).
a = BigFloat[BigFloat(factorial(big(k + 1))) for k in 0:49]

fit = stokes_fit(a)
fit.S, fit.β, fit.A   # ≈ (1, 2, 1)

# Recover the leading 1/k correction too:
fit2 = stokes_fit(a; subleading = 2)
fit2.c                # ≈ [c₁, c₂]
```

The fit is element-type generic (`Float64`, `BigFloat`, `Complex{T}`); for
PT-symmetric problems where the action is complex (`S = ±i|S|`), a
`Complex{T}` coefficient vector returns complex `S`, `β`, `A`.

## Tests

```julia
using Pkg; Pkg.test()
```
