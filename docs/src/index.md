# Resurgence.jl

Resummation techniques for divergent and asymptotic series, in pure Julia.

The package targets resurgence and perturbative-physics workflows: take a formal power series with factorially-growing coefficients, recover a Borel-summable result, and read off the large-order data (instanton action, exponent, prefactor) when you need it.
Every method is generic over the input element type — `Float64`, `BigFloat`, `Complex{T}` — and there are no Python dependencies.

## What's included

| Family                  | Methods |
| ----------------------- | ------- |
| Sequence acceleration   | Shanks (with iterated depth), Richardson |
| Padé approximants       | Linear, continued-fraction (qd), Hermite / quadratic |
| Borel side              | Borel, Borel–Le Roy, Borel–Padé, conformal Borel–Padé |
| Lateral / median sums   | Lateral, median, and discontinuity for Borel–Padé and Borel–Le Roy–Padé |
| Closed-form             | Meijer-G via Slater collapse onto `pFq` |
| Truncation              | Optimal truncation, superasymptotic remainder |
| Diagnostics             | Stokes action, exponent, constant, joint fit |

For deciding which method to use on which series, see the [Methods guide](methods_guide.md).

## Five-line install + use

```julia
using Pkg
Pkg.add(url = "https://github.com/schenkse/Resurgence.jl")

using Resurgence

a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
borel_pade(a; n = 10, m = 10, x = 1)   # ≈ 0.5963473623
```

The series `aₖ = (−1)ᵏ k!` is divergent for every nonzero `z` but Borel summable to `(1/z) e^{1/z} E₁(1/z)`.
At `z = 1` the exact value is `e · E₁(1) ≈ 0.5963473623`.
See [Getting started](getting_started.md) for the full first run.

## Two APIs over the same engine

A functional API:

```julia
borel_pade(a; n = 10, m = 10, x = 1)
borel_leroy_pade(a; n = 10, m = 10, x = 1)
conformal_borel_pade(a; n = 10, m = 10, x = 1, sing = 1)
```

…and a unified dispatch API on top:

```julia
resum(BorelPade(10, 10), a)
resum(BorelLeRoyPade(10, 10), a)
resum(ConformalBorelPade(10, 10; sing = 1), a)
```

Per-method functions are canonical; the `resum(::AbstractResummation, …)` layer is a thin dispatcher that lets you swap methods without changing caller code.

## Where to start

- New to the package?
  Read [Getting started](getting_started.md).
- Have a divergent series and want to know which method fits?
  Jump to the [Methods guide](methods_guide.md).
- Want worked physics-style examples?
  See the [Stieltjes tutorial](tutorials/stieltjes.md), [lateral and median sums](tutorials/lateral_sums.md), or [Stokes / large-order diagnostics](tutorials/stokes_diagnostics.md).
- Looking up a specific function?
  The API pages (sidebar) auto-generate from source docstrings.

## License and citation

Resurgence.jl is released under the MIT license.
If it's useful in published work, please cite the repository directly along with the relevant classical references — see the [References](references.md) page.
