# Stokes / large-order diagnostics

A divergent perturbative series with leading large-order behaviour

```math
a_{k+1} \;\sim\; A \cdot \frac{\Gamma(k + \beta)}{S^{\,k+\beta}}
                  \cdot \left( 1 + \frac{c_1}{k} + \frac{c_2}{k^2} + \cdots \right)
```

encodes its instanton action `S`, exponent `β`, prefactor `A`, and
subleading corrections `cⱼ` directly in the coefficient ratios.
Resurgence.jl reads them off via Richardson-accelerated tail
extrapolation.

## Driver: `aₖ = (k+1)!`

This driver has exact `(S, β, A) = (1, 2, 1)` and is the canonical
non-degenerate test case (see [`stokes_fit`](@ref) and
[CLAUDE.md](https://github.com/schenkse/Resurgence.jl/blob/main/CLAUDE.md)
for why the plain Stieltjes series `aₖ = (−1)ᵏ k!` is too degenerate).

```@example stokes
using Resurgence

# 50 BigFloat coefficients gives Richardson plenty of headroom.
a = BigFloat[BigFloat(factorial(big(k + 1))) for k in 0:49]
nothing # hide
```

## All-in-one fit

```@example stokes
fit = stokes_fit(a)
(S = fit.S, β = fit.β, A = fit.A)
```

Returns a named tuple with `S`, `β`, `A` — and (optionally) the subleading
`cⱼ` if you ask for them.

## Subleading corrections

```@example stokes
fit2 = stokes_fit(a; subleading = 2)
fit2.c
```

Order-by-order Richardson extrapolation reads off `[c₁, c₂]` directly.
Higher `subleading` values cost more coefficients (and more precision) for
diminishing accuracy returns.

## Single-quantity helpers

The underlying functions are also exposed individually. They take what
they need from earlier extractions: `stokes_exponent` needs `S`,
`stokes_constant` needs both `S` and `β`.

```@example stokes
S = stokes_action(a)
β = stokes_exponent(a, S)
A = stokes_constant(a, S, β)
(S, β, A)
```

Useful when you want the action `S` alone to feed into a hyperasymptotic
remainder estimate, or when you're sweeping a parameter and plotting how
`S` varies.

## Complex `S` for PT-symmetric problems

The fit is element-type generic. For problems where the dominant Borel
singularity is a complex-conjugate pair (`S = ±i|S|`), pass a `Complex{T}`
coefficient vector:

```julia
ac = ComplexF64.(a)
stokes_fit(ac).S    # complex S
```

`Float64` and `BigFloat` real input return real `S`, `β`, `A`; `Complex{T}`
input returns the corresponding complex tuple. No silent promotion.

## How the fit works (sketch)

The recipe (Bender & Wu 1969 / 1973 [BenderWu1969](@cite),
[BenderWu1973](@cite); Suslov 2005 [Suslov2005](@cite)):

1. From `aₖ`, form the ratio sequence `rₖ = aₖ₊₁ / aₖ`. To leading order
   `rₖ ≈ k / S`, so `1 / S` is the slope of the line `rₖ vs. k` at large
   `k`.
2. The exponent `β` shows up as the next-to-leading correction
   `rₖ ≈ k / S + (β − 1) / S + O(1/k)` — extract it from
   `S · rₖ − k`.
3. The prefactor `A` is read off from the asymptotic relation
   `aₖ ≈ A · Γ(k + β) / S^{k+β}` once `S` and `β` are known.
4. Richardson extrapolation (with iterated depth) accelerates the
   convergence of each tail-extracted quantity, which is critical when only
   30–60 coefficients are available.

For the underlying mathematical literature see the
[References](../references.md) page.

## When this won't work

- **Subleading singularities**. If the Borel transform has more than one
  singularity at comparable distance, the leading-order extraction picks up
  whichever wins; subleading `cⱼ` get contaminated. See the planned
  multi-singularity conformal map (B4) on the [roadmap](../roadmap.md).
- **Sign-alternating with sub-leading complex pair**. The leading
  large-order behaviour is real, but the next contribution can be complex
  — a `Complex{T}` fit picks this up.
- **Insufficient precision**. `Float64` runs out of digits for
  `Γ(k + β) / S^{k+β}` quickly. Use `BigFloat` whenever you suspect digit
  starvation.
