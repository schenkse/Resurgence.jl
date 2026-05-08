# Stieltjes/Euler series

The Stieltjes/Euler series

```math
S(z) \;=\; \sum_{k \ge 0} (-1)^k \, k! \, z^k
```

is divergent for every nonzero `z`, but Borel summable to

```math
S(z) \;=\; \frac{1}{z} \, e^{1/z} \, E_1\!\left(\tfrac{1}{z}\right)
\quad\text{for }z > 0.
```

At `z = 1` the exact value is `e · E₁(1) ≈ 0.5963473623`.
This page runs every method in the package against that reference value, side by side.

## Setup

```@example stieltjes
using Resurgence
using Printf

const STIELTJES_AT_1 = 0.5963473623231940743410784993

# 25 coefficients in Float64. factorial(big(k)) avoids overflow at k = 21.
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
nothing # hide
```

## Direct partial sum (asymptotic)

```@example stieltjes
partial = sum(a)
@printf("partial sum (all 25 terms): %+.6e   (err %.2e)\n",
        partial, partial - STIELTJES_AT_1)
```

The asymptotic series gets *worse* the more terms you add — Stieltjes is the textbook divergent example.

## Optimal truncation (superasymptotic)

```@example stieltjes
Nstar, partial_opt, εN = optimal_truncation(a)
@printf("optimal truncation (N* = %d): %+.6e   (err %.2e, est %.2e)\n",
        Nstar, partial_opt, partial_opt - STIELTJES_AT_1, εN)
```

`N*` is the index of the smallest term; the superasymptotic estimate of the remainder is `εN`.
Useful as a free lower bound — you usually want a real resummation method on top.

## Borel–Padé

```@example stieltjes
v_borel_pade = borel_pade(a; n = 10, m = 10, x = 1)
@printf("Borel–Padé[10/10]: %+.6e   (err %.2e)\n",
        v_borel_pade, v_borel_pade - STIELTJES_AT_1)
```

The standard tool for factorially divergent, alternating series with a single Borel-plane singularity on the negative real axis.

## Borel–Le Roy–Padé

```@example stieltjes
v_borel_lr = borel_leroy_pade(a; n = 10, m = 10, x = 1)   # b = -1/2 by default
@printf("Borel–Le Roy–Padé[10/10] (b = -1/2): %+.6e   (err %.2e)\n",
        v_borel_lr, v_borel_lr - STIELTJES_AT_1)
```

The Le Roy parameter `b` is a knob: setting `b = -1/2` is a common physics choice that often helps Padé convergence, but the optimal value is problem-dependent.
See [`borel_leroy_pade`](@ref) for the default and keyword arguments.

## Conformal Borel–Padé

```@example stieltjes
v_conformal = conformal_borel_pade(a; n = 10, m = 10, x = 1, sing = 1)
@printf("conformal-Borel–Padé[10/10]: %+.6e   (err %.2e)\n",
        v_conformal, v_conformal - STIELTJES_AT_1)
```

The conformal map `t → w(t)` accumulates the singularity at `t = -sing` to the boundary of the unit disk in `w`, then Padé works on the well-behaved re-expansion.
Best when you know where the leading Borel singularity is.

## Same problem at `BigFloat` precision

```@example stieltjes
ab = BigFloat[(-1)^k * factorial(big(k)) for k in 0:40]
vb = borel_pade(ab; n = 20, m = 20, x = BigFloat(1))
abs(vb - BigFloat(STIELTJES_AT_1))
```

`Float64` → `BigFloat` is a one-character change at the input.
The `pade` fast path uses LU; for rank-deficient fits it falls back to `pinv` via [GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl).

!!! note "Pushing the integration tolerance"
    `quadgk`'s integration nodes are computed in `Float64` by default.
    To push past ~`1e-14` you have to ask:

    ```julia
    borel_pade(ab; n = 20, m = 20, x = BigFloat(1),
               rtol = BigFloat("1e-30"), atol = BigFloat("1e-30"))
    ```

## Meijer-G

The Stieltjes series itself is a degenerate driver for Meijer-G (its Borel transform `1/(1+t)` is degree-1 rational, so the rational fit `P(k)/Q(k)` the method uses is rank-deficient at every order).
The shifted series `aₖ = (−1)ᵏ (k+1)!` has Borel transform `1/(1+t)²` and is the canonical non-degenerate driver.
Its Borel sum at `x = 1` is `1 − e · E₁(1) ≈ 0.40365`.

```@example stieltjes
const SHIFTED_REF = 1.0 - STIELTJES_AT_1
a_shifted = Float64[(-1.0)^k * Float64(factorial(big(k + 1))) for k in 0:24]

v_meijerg = borel_meijerg(a_shifted; n = 3, x = 1)
@printf("borel_meijerg (n = 3): %+.6e   (err %.2e)\n",
        v_meijerg, v_meijerg - SHIFTED_REF)
```

`borel_meijerg` collapses the resulting G-function onto a single `HypergeometricFunctions.pFq` call via the Slater identity (the doubled `b_h = 1` from Borel–Laplace cancels two of the numerator parameters), so no Meijer-G is ever evaluated directly — and there's no Python dependency.

## Side by side

| Method                      | Result        | Error vs. reference |
| --------------------------- | ------------- | ------------------- |
| Partial sum (25 terms)      | diverges      | grows with `k`      |
| Optimal truncation (N* = 6) | ~0.5963       | ~`1e-3`             |
| Borel–Padé[10/10]           | ~0.5963474    | ~`1e-7`             |
| Borel–Le Roy–Padé[10/10]    | ~0.5963474    | ~`1e-7`             |
| Conformal Borel–Padé[10/10] | ~0.5963474    | ~`1e-8`             |
| BigFloat Borel–Padé[20/20]  | full precision| `≲ 1e-30` w/ rtol  |

The point isn't that one method "wins" — it's that for a benign single-singularity problem like Stieltjes, all the Borel-side methods agree to many digits.
The differences show up on harder problems: see [Lateral and median sums](lateral_sums.md) for non-Borel-summable inputs.
