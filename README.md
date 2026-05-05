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

> **Built with LLMs:** This project was developed with the help of AI coding tools, primarily [Claude Code](https://claude.ai/code) by Anthropic. All code has been reviewed and is maintained by the author.

## Key features

- Element-type generic: `Float64`, `BigFloat`, and complex variants are all first-class.
- Multiple resummation families: Shanks, Richardson, Padé, Borel / Borel–Le Roy / Borel–Padé / conformal Borel–Padé, and optimal-truncation.
- Both a functional API and a unified `resum(::AbstractResummation, a)` dispatch layer.

## References

- C. M. Bender & S. A. Orszag, *Advanced Mathematical Methods for Scientists and Engineers* (1978) — asymptotic series and optimal truncation.
- O. Costin, *Asymptotics and Borel Summability* (2008) — modern treatment of Borel summation and resurgence.

## Setup

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
| Padé approximant `[n/m]` | `pade(a, m, n)`, `pade_value(a, m, n, x)` | `Pade(m, n; x)` |
| Borel transform | `borel_transform(a)` | — |
| Borel–Le Roy transform | `borel_leroy_transform(a, b)` | — |
| Borel–Padé resummation | `borel_pade(a; n, m, x, ...)` | `BorelPade(n, m; x, ...)` |
| Borel–Le Roy–Padé resummation | `borel_leroy_pade(a; n, m, b, x, ...)` | `BorelLeRoyPade(n, m; b, x, ...)` |
| Conformal-Borel–Padé resummation | `conformal_borel_pade(a; n, m, x, sing, ...)` | `ConformalBorelPade(n, m; x, sing, ...)` |
| Optimal-truncation / superasymptotics | `optimal_truncation(a)`, `superasymptotic_remainder(a)` | — |

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

## Tests

```julia
using Pkg; Pkg.test()
```
