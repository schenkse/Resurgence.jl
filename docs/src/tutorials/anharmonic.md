# Quartic anharmonic oscillator

The quartic anharmonic oscillator
```math
H \;=\; \tfrac{1}{2} p^2 + \tfrac{1}{2} x^2 + g \, x^4
```
is the textbook driver of resurgence: Bender and Wu showed in 1969 [BenderWu1969](@cite) that the ground-state perturbation series
```math
E_0(g) \;=\; \frac{1}{2} + \frac{3}{4} g - \frac{21}{8} g^2 + \frac{333}{16} g^3 - \frac{30\,885}{128} g^4 + \cdots
```
is divergent — its coefficients grow factorially with alternating sign — yet *Borel summable* to the exact ground-state energy.
The Borel transform has its nearest singularity at `t = -1/3` on the negative real axis (the instanton sets the scale `S = 1/3`), so the canonical Borel–Laplace integral is well defined and tools like [`borel_pade`](@ref) and [`conformal_borel_pade`](@ref) work out of the box.

This page runs every applicable resummation method on the Bender-Wu series at `g = 0.1` and compares against the high-precision reference `E_0(0.1) ≈ 0.5591463272…` from numerical diagonalisation.

## Computing the Bender-Wu coefficients

The recurrence below is Rayleigh-Schrödinger perturbation in the ansatz `ψ_n(x) = e^{-x²/2} · ∑_{k=0}^{2n} A_{n,k} x^{2k}`, applied to `H` order-by-order in `g`:

```math
A_{n,k} \;=\; \frac{(k+1)(2k+1)\, A_{n,k+1} \;-\; A_{n-1,\,k-2} \;+\; \sum_{m=1}^{n-1} E_0^{(m)}\, A_{n-m,\,k}}{2k},
```

with `A_{n,0} = 0` (intermediate normalisation) and `E_0^{(n)} = -A_{n,1}`.

```@example anharmonic
using Resurgence
using Printf

function bender_wu_aho_coefficients(N::Integer)
    R = Rational{BigInt}
    E = R[1//2]
    A = Vector{Vector{R}}()
    push!(A, R[1])
    for n in 1:N
        An = zeros(R, 2n + 1)
        for k in 2n:-1:1
            Akp1     = (k + 1 ≤ 2n)               ? An[k + 2]   : R(0)
            Anm1km2  = (0 ≤ k - 2 ≤ 2(n - 1))     ? A[n][k - 1] : R(0)
            s = R(0)
            for m in 1:n-1
                if k ≤ 2(n - m)
                    s += E[m + 1] * A[n - m + 1][k + 1]
                end
            end
            An[k + 1] = (R(k + 1) * R(2k + 1) * Akp1 - Anm1km2 + s) // R(2k)
        end
        push!(A, An)
        push!(E, -An[2])
    end
    return E
end

BW = bender_wu_aho_coefficients(20)
a  = Float64.(BW)
BW[1:8]
```

## Direct partial sum (asymptotic)

The series is asymptotic: summing more terms eventually makes the result *worse*.

```@example anharmonic
const G      = 0.1
const E0_REF = 0.55914633237109099
partial = sum(a[k+1] * G^k for k in 0:length(a)-1)
@sprintf("partial sum (all 21 terms): %+.6e   (err %.2e)", partial, partial - E0_REF)
```

## Optimal truncation (superasymptotic)

```@example anharmonic
ak_g = [a[k+1] * G^k for k in 0:length(a)-1]
Nstar, partial_opt, εN = optimal_truncation(ak_g)
@sprintf("optimal truncation (N* = %d): %+.10f   (err %.2e, est %.2e)",
         Nstar, partial_opt, partial_opt - E0_REF, εN)
```

`N* = 4` for `g = 0.1` matches the textbook smallest-term estimate `N* ≈ 1/(3g)`.
The superasymptotic error sits at `1e-2` — useful as a free lower bound, but a real Borel method does dramatically better.

## Borel-Padé

```@example anharmonic
v_bp = borel_pade(a; n = 10, m = 10, x = G)
@sprintf("Borel-Padé[10/10]: %+.10f   (err %.2e)", v_bp, v_bp - E0_REF)
```

Nine decimals of agreement from a divergent series of 21 terms — this is what Borel-Laplace buys you on an alternating Borel-summable input.

## Borel–Le Roy–Padé

```@example anharmonic
v_lr = borel_leroy_pade(a; n = 10, m = 10, x = G)   # b = -1/2 default
@sprintf("Borel-Le Roy-Padé[10/10] (b = -1/2): %+.10f   (err %.2e)", v_lr, v_lr - E0_REF)
```

The Le Roy parameter `b = -1/2` is the conventional physics choice; for the anharmonic oscillator it ties the Borel and Borel–Le Roy answers to roughly the same precision.

## Conformal Borel-Padé

Plugging in the known singularity location `t = -sing = -1/3` lets the conformal map push the Borel singularity to the unit-disk boundary before the Padé fit, which often improves convergence:

```@example anharmonic
v_cb = conformal_borel_pade(a; n = 10, m = 10, x = G, sing = 1//3)
@sprintf("conformal-Borel-Padé[10/10] (sing = 1/3): %+.10f   (err %.2e)", v_cb, v_cb - E0_REF)
```

## Diagnostics

Confirm the Bender-Wu large-order structure from the coefficients themselves.
For an alternating series the extracted action `S` is negative real, so we pass `ComplexF64.(a)` to keep [`stokes_fit`](@ref) happy when it raises `S^(k+β)`:

```@example anharmonic
d = diagnose(ComplexF64.(a))
(growth = d.growth, alternating = d.alternating, S = d.S, β = d.β, A = d.A,
 recommended = d.recommended)
```

`S ≈ -0.333` recovers the Bender-Wu prediction `S = -1/3` — the Borel-plane singularity sits at `t = -1/3` exactly as expected.

For raw numerical inspection (no Stokes fit, no method routing):

```@example anharmonic
cd = coefficient_diagnostics(a)
(growth = cd.growth, ratio_growth = cd.ratio_growth,
 alternation_score = cd.alternation_score)
```

## Side-by-side via `compare`

```@example anharmonic
rows = compare([
    BorelPade(10, 10; x = G),
    BorelLeRoyPade(10, 10; x = G),
    ConformalBorelPade(10, 10; x = G, sing = 1//3),
    Pade(10, 10; x = G),
], a; reference = E0_REF)
```

## High precision via BigFloat

The same code in `BigFloat` arithmetic, requesting tighter `quadgk` tolerances, pushes the answer well past Float64 precision:

```@example anharmonic
a_big = BigFloat.(BW)
setprecision(BigFloat, 256) do
    v_big = borel_pade(a_big; n = 10, m = 10, x = big"0.1",
                        rtol = BigFloat("1e-30"), atol = BigFloat("1e-30"))
    v_big
end
```

The residual `~1e-14` is dominated by `QuadGK`'s Float64-default integration nodes — push that further by raising `order` and tightening `atol`/`rtol`.
