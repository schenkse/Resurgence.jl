# Getting started

This page walks through installing the package, computing one Borel sum,
and pointing you at the next read.

## Install

Resurgence.jl is not yet registered in Julia's General registry, so install
from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/schenkse/Resurgence.jl")
```

For development, clone and `Pkg.activate`:

```bash
git clone https://github.com/schenkse/Resurgence.jl
cd Resurgence.jl
```

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
using Resurgence
```

The package targets Julia ≥ 1.10 (LTS).

## A first end-to-end run

The Stieltjes / Euler series `S(z) = Σ (−1)ᵏ k! zᵏ` is the canonical
divergent-but-Borel-summable example. At `z = 1` it sums to
`e · E₁(1) ≈ 0.5963473623`. Here's the full pipeline:

```@example getting_started
using Resurgence

# 25 coefficients in Float64. Use factorial(big(k)) so we don't overflow at k = 21.
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]

# Direct partial sum — diverges; the more terms you add, the worse it gets.
sum(a)
```

```@example getting_started
# Optimal truncation: stop at the smallest term. Superasymptotic accuracy.
Nstar, partial_opt, εN = optimal_truncation(a)
(Nstar, partial_opt, εN)
```

```@example getting_started
# Borel–Padé: the standard tool for factorially divergent, alternating series.
borel_pade(a; n = 10, m = 10, x = 1)
```

The reference value is `0.5963473623…`, so Borel–Padé already gives ~8 digits
from 25 input coefficients.

## Higher precision

Every method propagates the input element type. Hand it `BigFloat` and the
intermediate Padé fit and Laplace integral are `BigFloat` too.

```@example getting_started
using Resurgence  # hide

ab = BigFloat[(-1)^k * factorial(big(k)) for k in 0:40]
borel_pade(ab; n = 20, m = 20, x = BigFloat(1))
```

!!! note "QuadGK and BigFloat"
    `quadgk` evaluates integration nodes in `Float64` by default, so even
    with `BigFloat` input the quadrature error sits around `1e-14` unless
    you ask for tighter tolerances. Pass `rtol = BigFloat("1e-30")` (and a
    matching `atol`) explicitly to push further:

    ```julia
    borel_pade(ab; n = 20, m = 20, x = BigFloat(1),
               rtol = BigFloat("1e-30"), atol = BigFloat("1e-30"))
    ```

## Functional vs. unified API

Both of these compute the same number:

```@example getting_started
using Resurgence  # hide
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]  # hide

borel_pade(a; n = 10, m = 10, x = 1)
```

```@example getting_started
resum(BorelPade(10, 10), a)
```

Use whichever fits your code. Per-method functions are canonical; the
`resum(::AbstractResummation, …)` layer is a thin dispatch shim that lets
callers swap methods by changing one struct.

## Next steps

- Side-by-side comparison of every method on the Stieltjes series:
  [Stieltjes tutorial](tutorials/stieltjes.md).
- Series whose Borel transform has a singularity on the *positive* real axis
  (so the standard Laplace integral is ill-defined):
  [Lateral and median sums](tutorials/lateral_sums.md).
- Reading the instanton action `S`, exponent `β`, and prefactor `A` off the
  large-order behaviour of the coefficients themselves:
  [Stokes / large-order diagnostics](tutorials/stokes_diagnostics.md).
- Picking a method for your problem:
  [Methods guide](methods_guide.md).
