# Lateral and median sums

When the Borel transform of a series has a singularity on the *positive* real axis, the standard Laplace integral

```math
S(x) \;=\; \int_0^\infty e^{-t/x} \, B(t) \, dt
```

runs straight through that singularity and is ill-defined.
The classical fix is to deform the contour above (`L⁺`) or below (`L⁻`) the singularity.
That gives two distinct lateral sums, each complex-valued.
Their average is the ambiguity-free **median** sum and their difference exposes the **Stokes discontinuity** — the imaginary, non-perturbative content the plain perturbation series can't see on its own.

## Driver: `aₖ = k!`

The series `S(z) = Σ k! zᵏ` has Borel transform `1/(1−t)`, with a pole at `t = 1`.
At `z = 1` the median sum recovers the Cauchy principal value `Ei(1)/e ≈ 0.6972`, and the Stokes discontinuity has magnitude `π/e ≈ 1.1557`.

```@example lateral
using Resurgence
using Printf

a = Float64[Float64(factorial(big(k))) for k in 0:24]   # a_k = k!
nothing # hide
```

## Median sum

```@example lateral
v_med = borel_pade_median(a; n = 10, m = 10, x = 1)
@printf("median:        %.6f + %.6fim\n", real(v_med), imag(v_med))
```

The imaginary part should be numerically zero — `borel_pade_median` averages `L⁺` and `L⁻`, which are complex conjugates of each other for real input.

## Stokes discontinuity

```@example lateral
v_disc = borel_pade_discontinuity(a; n = 10, m = 10, x = 1)
@printf("discontinuity: %.6f + %.6fim\n", real(v_disc), imag(v_disc))
```

Magnitude `≈ π/e ≈ 1.1557`, real-valued.
This is the imaginary jump that encodes the leading instanton contribution in resurgence-style problems.

## The two lateral sums

```@example lateral
v_plus  = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = +1)
v_minus = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = -1)
@printf("L+ (side = +1): %.6f + %.6fim\n", real(v_plus),  imag(v_plus))
@printf("L- (side = -1): %.6f + %.6fim\n", real(v_minus), imag(v_minus))
```

For real input these are exact complex conjugates of each other.

!!! warning "Sign convention"
    `side = +1` shifts the offending Padé poles by `+iε` (upper half-plane), so the integration contour passes *below* them — matching the lateral sum traditionally written `S_θ⁻` in the resurgence literature.

    `borel_pade_discontinuity` is computed as

    ```math
    (L_{+1} - L_{-1}) \big/ (2i) \;=\; (S_\theta^- - S_\theta^+) \big/ (2i),
    ```

    which has magnitude equal to the textbook `(S_θ⁺ − S_θ⁻)/(2i)` but the *opposite sign*.
    Negate the return value if you want the textbook convention.

## Le Roy variants

The same lateral/median/discontinuity trio exists for Borel–Le Roy:

```julia
borel_leroy_pade_lateral(a; n, m, b, x, side, ...)
borel_leroy_pade_median(a;  n, m, b, x, ...)
borel_leroy_pade_discontinuity(a; n, m, b, x, ...)
```

Useful when the Le Roy parameter `b` improves Padé conditioning on the particular series you have in hand.

## When to reach for these

- **Median sum**: the right "value" to report when the standard Borel sum doesn't exist.
  Real-valued for real input series with complex-conjugate singularity pairs or with isolated positive-real-axis singularities.
- **Discontinuity**: the numerical face of resurgence — proportional to the leading non-perturbative contribution `e^{−S/g}` weighted by the Stokes constant.
  To extract `S` and `β` directly from the coefficients, see [Stokes/large-order diagnostics](stokes_diagnostics.md).
- **Lateral sums individually**: useful if downstream code expects a particular contour choice rather than the median.
