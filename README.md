# Resurgence.jl

Resummation techniques for divergent / asymptotic series, in pure Julia.

[![Test and Release](https://github.com/schenkse/Resurgence.jl/actions/workflows/test-release.yml/badge.svg)](https://github.com/schenkse/Resurgence.jl/actions/workflows/test-release.yml)
[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://schenkse.github.io/Resurgence.jl/stable/)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://schenkse.github.io/Resurgence.jl/dev/)
[![Julia ≥ 1.10](https://img.shields.io/badge/Julia-≥1.10-9558B2?logo=julia)](https://julialang.org)

The package targets resurgence and perturbative-physics workflows: take a
formal power series with factorially-growing coefficients, recover a
Borel-summable result. All methods are generic over the input element type —
`Float64`, `BigFloat`, and complex variants are all first-class.

> **Built with LLMs:** This project was developed with the help of AI coding
> tools, primarily [Claude Code](https://claude.ai/code) by Anthropic. All
> code has been reviewed and is maintained by the author.

## Install

> Resurgence.jl is not yet registered in Julia's General registry, so
> `Pkg.add("Resurgence")` will not work — install from GitHub for now.

```julia
using Pkg
Pkg.add(url = "https://github.com/schenkse/Resurgence.jl")
using Resurgence
```

## Quick start

The Stieltjes / Euler series `S(z) = Σ (-1)^k k! z^k` is divergent for every
nonzero `z`, but Borel summable to `(1/z) e^{1/z} E₁(1/z)`. At `z = 1` the
exact value is `e · E₁(1) ≈ 0.5963473623`.

```julia
using Resurgence

# 25 coefficients  a[k+1] = (-1)^k k!
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]

borel_pade(a; n = 10, m = 10, x = 1)             # ≈ 0.59634736
resum(BorelPade(10, 10), a)                       # same, via the unified API
```

## Documentation

Full docs, tutorials, methods guide, and API reference:
**[schenkse.github.io/Resurgence.jl](https://schenkse.github.io/Resurgence.jl/dev/)**

Tutorials cover the [Stieltjes / Euler series](https://schenkse.github.io/Resurgence.jl/dev/tutorials/stieltjes/),
[lateral and median sums](https://schenkse.github.io/Resurgence.jl/dev/tutorials/lateral_sums/) for
non-Borel-summable inputs, and
[Stokes / large-order diagnostics](https://schenkse.github.io/Resurgence.jl/dev/tutorials/stokes_diagnostics/)
for reading `S`, `β`, `A` off the coefficients.

## Tests

```julia
using Pkg; Pkg.test()
```
