# Methods guide

Which method should you use for which problem?
Resurgence.jl ships several families that overlap in coverage but disagree on cost and on what information they need from you.
This page is a decision tree.

## Quick lookup

| Your series looks like…                                               | Try first                                          | Then maybe                                                                 |
| --------------------------------------------------------------------- | -------------------------------------------------- | -------------------------------------------------------------------------- |
| Convergent but slow                                                   | [`shanks`](@ref) / [`richardson`](@ref)            | [`pade`](@ref) for analytic-continuation needs                             |
| Factorially divergent, alternating sign (Borel-summable)              | [`borel_pade`](@ref)                               | [`conformal_borel_pade`](@ref) if you know where the singularity is        |
| Factorially divergent, all-positive (non-Borel-summable)              | [`borel_pade_median`](@ref)                        | [`borel_pade_lateral`](@ref) + [`borel_pade_discontinuity`](@ref)          |
| Borel transform has known branch point at `t = -sing`                 | [`conformal_borel_pade`](@ref)                     | [`borel_leroy_pade`](@ref) with tuned `b`                                  |
| Borel transform has algebraic branch point in `x`                     | [`hermite_pade`](@ref)                             | [`hermite_pade_value`](@ref) with branch selector                          |
| You suspect closed form / hypergeometric                              | [`borel_meijerg`](@ref)                            | The Stieltjes tutorial covers degenerate cases                             |
| Just need a quick error bound on a partial sum                        | [`optimal_truncation`](@ref)                       | [`superasymptotic_remainder`](@ref) for the tail estimate alone            |
| Want `S`, `β`, `A` from large-order behaviour                         | [`stokes_fit`](@ref)                               | [`stokes_action`](@ref) / [`stokes_exponent`](@ref) one at a time          |

## Decision tree

```
                ┌──────────────────────────────────────┐
                │ Does the partial sum settle as you   │
                │ add terms?                           │
                └─┬────────────────────────────────┬───┘
                  │ Yes (convergent)               │ No (divergent)
                  ▼                                ▼
             ┌────────┐                     ┌──────────────────────────┐
             │ shanks │                     │ Are aₖ alternating in    │
             │ pade   │                     │ sign?                    │
             └────────┘                     └─┬──────────────────────┬─┘
                                              │ Yes                  │ No
                                              ▼                      ▼
                                       ┌──────────────┐       ┌────────────────────┐
                                       │ borel_pade   │       │ borel_pade_median  │
                                       │ borel_leroy  │       │ + _lateral / _disc │
                                       │ conformal_   │       │ for ambiguity      │
                                       │   borel_pade │       └────────────────────┘
                                       └──────────────┘
```

If you don't know whether the series is alternating, look at [`stokes_fit`](@ref): real `S > 0` means the leading Borel singularity is on the negative real axis (alternating-sign / Borel-summable case); real `S < 0` puts it on the positive real axis (non-Borel-summable — reach for the lateral / median methods).

## Sequence acceleration: Shanks vs. Richardson

Both [`shanks`](@ref) and [`richardson`](@ref) take a coefficient vector and a starting offset `n`, and both accept iterated `depth` (apply the transform `depth` times).

- **Shanks** (Aitken-Δ²) [Shanks1955](@cite) cancels a geometric tail `Aρⁿ` exactly.
  Best when the partial sums look like a geometric approach to the limit.
- **Richardson** [BenderOrszag1978](@cite) cancels a polynomial tail `c/n + c'/n² + …`.
  Best when the partial sums approach the limit like `1/nᵖ`.

If you don't know which tail you have, try both and pick the one whose output stabilises faster.

## Padé family

Three flavours, increasing in specialty:

- [`pade`](@ref) — linear `[n/m]` Padé.
  Solves an `(n+m+1)×(n+m+1)` linear system; falls back to `pinv` (via [GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl)) for rank-deficient inputs.
  The default Padé you want.
- [`pade_cf`](@ref) — Padé via the qd algorithm, returning a continued fraction.
  Cheaper to evaluate at many `x` values once built, but only works for `n == m` or `n + 1 == m`.
- [`hermite_pade`](@ref) — quadratic / Hermite Padé.
  Use this when the function you're approximating has a *branch point* in `x` rather than a pole — algebraic branch points need quadratic, not linear, denominators.

## Borel-side methods

All Borel methods share the same three-step pipeline:

1. Construct the Borel (or Borel–Le Roy) transform of the coefficient sequence in `t`.
2. Build a Padé approximant of the transform to extend it past the radius of convergence in `t`.
3. Laplace-integrate against `e^{-t/x}` to recover the resummed value at `x`.

The variants differ in how they handle the Padé step:

- [`borel_pade`](@ref) — vanilla Padé in `t`.
  Works when the Borel transform has a single isolated singularity on the negative real axis.
- [`borel_leroy_pade`](@ref) — Borel–Le Roy with parameter `b`: the factorial damping is `Γ(k + 1 + b)` instead of `k!`.
  The `b = -1/2` default is a common physics choice; sweeping `b` by hand often improves convergence on a particular series.
- [`conformal_borel_pade`](@ref) — first apply a conformal map that pushes the known Borel singularity to the boundary of the unit disk, then Padé.
  Best when you know `sing` (or have estimated it via [`stokes_action`](@ref)).
  Subsumes `borel_pade` for benign cases and beats it when the singularity is known.

For a series whose Borel transform has a singularity on the *positive* real axis (`S < 0`), the standard Laplace integral is ill-defined.
Use the lateral / median variants:

- [`borel_pade_lateral`](@ref) `(side = ±1)` — contour-deformed integral.
- [`borel_pade_median`](@ref) — `(L⁺ + L⁻) / 2`, the ambiguity-free median.
- [`borel_pade_discontinuity`](@ref) — `(L⁺ − L⁻) / (2i)`, the imaginary Stokes jump.

The Le Roy variants (`borel_leroy_pade_lateral`, `_median`, `_discontinuity`) parallel these one-for-one.

## Closed-form: Meijer-G

[`borel_meijerg`](@ref) tries to identify the Borel transform with a hypergeometric `pFq`, which when Laplace-integrated gives a Meijer-G function.
The implementation collapses the resulting G onto a single `HypergeometricFunctions.pFq` call via the Slater identity (the doubled `b_h = 1` from Borel–Laplace cancels two of the numerator parameters), so no Meijer-G is ever evaluated directly.

This wins when it works — closed-form-quality answers from a finite series — but the rational fit `P(k) / Q(k)` it does internally is rank-deficient on degenerate drivers (the unshifted Stieltjes series is one such trap).
See the [Stieltjes tutorial](tutorials/stieltjes.md) for details.

## Truncation diagnostics

[`optimal_truncation`](@ref) returns `(N*, partial_sum, εN)`: stop summing at the smallest term, take the partial sum, and an estimate of the error.
This is "free" — it doesn't compute anything beyond the partial sums you already have — and is a useful sanity check on top of any of the resummation methods above.

[`superasymptotic_remainder`](@ref) is the underlying tail-size estimate on its own.

## Stokes diagnostics

[`stokes_fit`](@ref) reads `S`, `β`, `A` (and optionally subleading `cⱼ`) off the coefficient ratios `aₖ₊₁ / aₖ` via Richardson-accelerated tail extrapolation.
Single-quantity helpers exist: [`stokes_action`](@ref), [`stokes_exponent`](@ref), [`stokes_constant`](@ref).

Useful when:

- You want to know where the leading Borel singularity is (`sing = S`) before reaching for [`conformal_borel_pade`](@ref).
- You're studying the resurgent structure of a problem and want to see whether `S` matches a known instanton action.
- You want to compare the discontinuity from [`borel_pade_discontinuity`](@ref) against the analytic prediction `2π · A / Γ(β) · e^{-S/x}`.

See [Stokes / large-order diagnostics](tutorials/stokes_diagnostics.md) for a full worked example.

## Unified API

Every method has a struct counterpart for the `resum(::AbstractResummation, a)` dispatcher:

```julia
resum(BorelPade(10, 10),               a; x = 1)
resum(BorelLeRoyPade(10, 10; b = -1/2), a; x = 1)
resum(ConformalBorelPade(10, 10),      a; x = 1, sing = 1)
resum(MeijerG(3),                       a; x = 1)
```

Use the per-method functions when the call is one-off; use the structs when you want to pass a method around as a value (sweeping parameters, configuring a benchmark, etc.).
The struct layer is intentionally thin; all the math lives in the per-method functions.
