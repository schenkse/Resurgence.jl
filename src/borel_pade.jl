"""
    borel_pade(a; n, m, x = 1, regularize_poles = false, ε = nothing,
               side = +1, return_error = false, quad_kwargs...)

Borel–Padé resummation of the formal power series with coefficients `a`,
evaluated at argument `x`:

1. Borel-transform `a` to obtain `B[k] = a[k] / k!`.
2. Replace `B(t)` by its Padé approximant `[n/m]` `Bₚ(t) / B_q(t)`.
3. Laplace-transform `Bₚ(t) / B_q(t)` along the positive real axis:

        ∫₀^∞ (Bₚ(x t) / B_q(x t)) e^{-t} dt.

When `regularize_poles = true`, any poles of the Padé denominator on the
positive real axis are first shifted by `side · iε` (default
`ε = 100·eps(real(T))`, default `side = +1`) to push them off the contour.
This branch is only meaningful for `x > 0` because the integration runs along
positive `t`; for `x ≤ 0` an `ArgumentError` is raised. Selecting `side = +1`
or `side = −1` produces the two lateral Borel sums; see
[`borel_pade_lateral`](@ref).

`quad_kwargs` are forwarded to `QuadGK.quadgk` (e.g. `rtol`, `atol`, `order`).
With `return_error = true`, returns a `(value, abserr)` tuple; otherwise just
the value.

Requires `length(a) ≥ n + m + 1`.
"""
function borel_pade(a::AbstractVector{T};
                    n::Integer, m::Integer, x = 1,
                    regularize_poles::Bool = false,
                    ε::Union{Real,Nothing} = nothing,
                    side::Integer = +1,
                    return_error::Bool = false,
                    quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    if regularize_poles && !(x isa Real && x > 0)
        throw(ArgumentError("regularize_poles=true requires positive real x (got $x)"))
    end
    Ba = borel_transform(a[1:n+m+1])
    Bp, Bq = pade(Ba, n, m)

    integrand = if regularize_poles
        ε_eff = ε === nothing ? 100 * eps(real(T)) : ε
        Bqr = poles_regularized(Polynomials.coeffs(Bq), ε_eff; side = side)
        # default to order=24 for the regularized contour unless caller overrides
        if !haskey(quad_kwargs, :order)
            quad_kwargs = (; quad_kwargs..., order = 24)
        end
        t -> Bp(x * t) / Bqr(x * t) * exp(-t)
    else
        t -> Bp(x * t) / Bq(x * t) * exp(-t)
    end

    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    borel_leroy_pade(a; n, m, b = -1//2, x = 1, regularize_poles = false,
                     ε = nothing, side = +1, return_error = false, quad_kwargs...)

Borel–Le Roy–Padé resummation:

1. Borel–Le Roy transform `a` with parameter `b` to obtain `B[k] = a[k] / Γ(k+1+b)`.
2. Replace `B(t)` by its Padé approximant `[n/m]`.
3. Laplace-transform with the Le Roy weight:

        ∫₀^∞ t^b (Bₚ(x t) / B_q(x t)) e^{-t} dt.

`b` may be any real number (including non-integer values); the default `-1//2`
matches the original code. Pole regularization (`regularize_poles`, `ε`,
`side`) works as in [`borel_pade`](@ref); selecting opposite sides produces the
two lateral Le Roy–Borel sums. `quad_kwargs` are forwarded to `QuadGK.quadgk`.
"""
function borel_leroy_pade(a::AbstractVector{T};
                          n::Integer, m::Integer,
                          b::Real = -1//2, x = 1,
                          regularize_poles::Bool = false,
                          ε::Union{Real,Nothing} = nothing,
                          side::Integer = +1,
                          return_error::Bool = false,
                          quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_leroy_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    if regularize_poles && !(x isa Real && x > 0)
        throw(ArgumentError("regularize_poles=true requires positive real x (got $x)"))
    end
    Ba = borel_leroy_transform(a[1:n+m+1], b)
    Bp, Bq = pade(Ba, n, m)

    integrand = if regularize_poles
        ε_eff = ε === nothing ? 100 * eps(real(T)) : ε
        Bqr = poles_regularized(Polynomials.coeffs(Bq), ε_eff; side = side)
        if !haskey(quad_kwargs, :order)
            quad_kwargs = (; quad_kwargs..., order = 24)
        end
        t -> t^b * Bp(x * t) / Bqr(x * t) * exp(-t)
    else
        t -> t^b * Bp(x * t) / Bq(x * t) * exp(-t)
    end

    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    conformal_borel_pade(a; n, m, x = 1, sing = 1, return_error = false, quad_kwargs...)

Conformal-Borel–Padé resummation:

1. Borel-transform `a`.
2. Re-expand the Borel transform under the conformal map
   `t = 4·sing·w / (1-w)²` (see [`conformal_reseries`](@ref)).
3. Padé-approximate the conformal series in `w`.
4. Map back via `w = (√(1+t/sing) - 1)/(√(1+t/sing) + 1)` and Laplace-integrate
   along the positive real `t` axis.

`sing > 0` is the assumed location of the nearest singularity of the Borel
transform on the negative real axis (i.e. the singularity sits at `t = -sing`).
"""
function conformal_borel_pade(a::AbstractVector{T};
                              n::Integer, m::Integer, x = 1,
                              sing::Real = 1,
                              return_error::Bool = false,
                              quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("conformal_borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_transform(a[1:n+m+1])
    Bw = conformal_reseries(Ba, sing, n + m)
    Bp, Bq = pade(Bw, n, m)
    integrand = function (t)
        w = conformal_map(x * t; a = sing)
        return Bp(w) / Bq(w) * exp(-t)
    end
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    conformal_borel_pade_pair(a; n, m, x = 1, sing = 1, return_error = false, quad_kwargs...)

Conformal-Borel–Padé resummation for a complex-conjugate Borel-plane
singularity pair at `t = ±i·sing`:

1. Borel-transform `a`.
2. Re-expand the Borel transform under the singularity-pair conformal map
   `t = 2·sing·v / (1 - v²)` (see [`conformal_reseries_pair`](@ref)).
3. Padé-approximate the conformal series in `v`.
4. Map back via `v = (√(t² + sing²) − sing) / t` (numerically as
   `v = t / (√(t² + sing²) + sing)`) and Laplace-integrate along the
   positive real `t` axis.

`sing > 0` is the assumed magnitude of the imaginary singularity pair
(singularities at `t = ±i·sing`). PT-symmetric problems and the quartic
anharmonic oscillator are typical use cases.
"""
function conformal_borel_pade_pair(a::AbstractVector{T};
                                   n::Integer, m::Integer, x = 1,
                                   sing::Real = 1,
                                   return_error::Bool = false,
                                   quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("conformal_borel_pade_pair needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_transform(a[1:n+m+1])
    Bv = conformal_reseries_pair(Ba, sing, n + m)
    Bp, Bq = pade(Bv, n, m)
    integrand = function (t)
        v = conformal_map_pair(x * t; a = sing)
        return Bp(v) / Bq(v) * exp(-t)
    end
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    borel_pade_lateral(a; n, m, x = 1, side = +1, ε = nothing, kwargs...)

Lateral Borel–Padé sum: shorthand for
`borel_pade(a; n, m, x, regularize_poles = true, side, ε, kwargs...)`. The
positive-real-axis poles of the Padé denominator are pushed off the contour
by `side · iε`, giving one of the two lateral Borel sums. For real `a` and
real `x > 0` the result is generally complex.

Convention: `side = +1` shifts poles to `+iε` (upper half-plane), so the
real-axis integration contour passes below them — this matches the lateral
sum traditionally written `S_θ⁻` in the resurgence literature; `side = -1`
gives `S_θ⁺`. Requires `x > 0`.

Defaults `ε` to `√eps(real(eltype(a)))`. The lateral integrand has a
near-singularity of height `∼1/ε` at the shifted pole, so the very tight
default of `borel_pade` (`100·eps`) makes `quadgk` underflow; `√eps` keeps
both the pole shift small and the integrand resolvable.
"""
function borel_pade_lateral(a::AbstractVector{T};
                            n::Integer, m::Integer, x = 1,
                            side::Integer = +1,
                            ε::Union{Real,Nothing} = nothing,
                            kwargs...) where {T<:Number}
    ε_use = ε === nothing ? sqrt(eps(real(T))) : ε
    return borel_pade(a; n = n, m = m, x = x,
                      regularize_poles = true, side = side, ε = ε_use, kwargs...)
end

"""
    borel_pade_median(a; n, m, x = 1, ε = nothing, kwargs...)

Median Borel–Padé sum `(L⁺ + L⁻)/2`, the average of the two lateral sums. For
real `a` the result is purely real (up to roundoff and `quadgk` tolerance);
this is the ambiguity-free real-valued resummation of a non-Borel-summable
series.
"""
function borel_pade_median(a::AbstractVector;
                           n::Integer, m::Integer, x = 1,
                           ε::Union{Real,Nothing} = nothing,
                           kwargs...)
    Lp = borel_pade_lateral(a; n = n, m = m, x = x, side = +1, ε = ε, kwargs...)
    Lm = borel_pade_lateral(a; n = n, m = m, x = x, side = -1, ε = ε, kwargs...)
    return (Lp + Lm) / 2
end

"""
    borel_pade_discontinuity(a; n, m, x = 1, ε = nothing, kwargs...)

Stokes discontinuity `(L⁺ − L⁻) / (2i)` of the Borel–Padé sum. Encodes the
ambiguity of the lateral sums and is proportional to the leading instanton
contribution `Stokes_constant · e^{−|S|/x}`. For real `a` the result is real.
"""
function borel_pade_discontinuity(a::AbstractVector;
                                  n::Integer, m::Integer, x = 1,
                                  ε::Union{Real,Nothing} = nothing,
                                  kwargs...)
    Lp = borel_pade_lateral(a; n = n, m = m, x = x, side = +1, ε = ε, kwargs...)
    Lm = borel_pade_lateral(a; n = n, m = m, x = x, side = -1, ε = ε, kwargs...)
    return (Lp - Lm) / (2 * im)
end

"""
    borel_leroy_pade_lateral(a; n, m, b = -1//2, x = 1, side = +1, ε = nothing, kwargs...)

Lateral Borel–Le Roy–Padé sum. Same lateral conventions as
[`borel_pade_lateral`](@ref); the only difference is the Le Roy weight `t^b`
in the Laplace integrand. `ε` defaults to `√eps(real(eltype(a)))`, large
enough to keep the integrand resolvable across the shifted pole.
"""
function borel_leroy_pade_lateral(a::AbstractVector{T};
                                  n::Integer, m::Integer,
                                  b::Real = -1//2, x = 1,
                                  side::Integer = +1,
                                  ε::Union{Real,Nothing} = nothing,
                                  kwargs...) where {T<:Number}
    ε_use = ε === nothing ? sqrt(eps(real(T))) : ε
    return borel_leroy_pade(a; n = n, m = m, b = b, x = x,
                            regularize_poles = true, side = side,
                            ε = ε_use, kwargs...)
end

"""
    borel_leroy_pade_median(a; n, m, b = -1//2, x = 1, ε = nothing, kwargs...)

Median Borel–Le Roy–Padé sum `(L⁺ + L⁻) / 2`.
"""
function borel_leroy_pade_median(a::AbstractVector;
                                 n::Integer, m::Integer,
                                 b::Real = -1//2, x = 1,
                                 ε::Union{Real,Nothing} = nothing,
                                 kwargs...)
    Lp = borel_leroy_pade_lateral(a; n = n, m = m, b = b, x = x,
                                  side = +1, ε = ε, kwargs...)
    Lm = borel_leroy_pade_lateral(a; n = n, m = m, b = b, x = x,
                                  side = -1, ε = ε, kwargs...)
    return (Lp + Lm) / 2
end

"""
    borel_leroy_pade_discontinuity(a; n, m, b = -1//2, x = 1, ε = nothing, kwargs...)

Stokes discontinuity `(L⁺ − L⁻) / (2i)` of the Borel–Le Roy–Padé sum.
"""
function borel_leroy_pade_discontinuity(a::AbstractVector;
                                        n::Integer, m::Integer,
                                        b::Real = -1//2, x = 1,
                                        ε::Union{Real,Nothing} = nothing,
                                        kwargs...)
    Lp = borel_leroy_pade_lateral(a; n = n, m = m, b = b, x = x,
                                  side = +1, ε = ε, kwargs...)
    Lm = borel_leroy_pade_lateral(a; n = n, m = m, b = b, x = x,
                                  side = -1, ε = ε, kwargs...)
    return (Lp - Lm) / (2 * im)
end
