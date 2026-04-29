"""
    borel_pade(a; n, m, x = 1, regularize_poles = false, ε = nothing,
               return_error = false, quad_kwargs...)

Borel–Padé resummation of the formal power series with coefficients `a`,
evaluated at argument `x`:

1. Borel-transform `a` to obtain `B[k] = a[k] / k!`.
2. Replace `B(t)` by its Padé approximant `[n/m]` `Bₚ(t) / B_q(t)`.
3. Laplace-transform `Bₚ(t) / B_q(t)` along the positive real axis:

        ∫₀^∞ (Bₚ(x t) / B_q(x t)) e^{-t} dt.

When `regularize_poles = true`, any poles of the Padé denominator on the
positive real axis are first shifted by `+iε` (default `ε = 100·eps(real(T))`)
to push them off the contour. This branch is only meaningful for `x > 0`
because the integration runs along positive `t`; for `x ≤ 0` an
`ArgumentError` is raised.

`quad_kwargs` are forwarded to `QuadGK.quadgk` (e.g. `rtol`, `atol`, `order`).
With `return_error = true`, returns a `(value, abserr)` tuple; otherwise just
the value.

Requires `length(a) ≥ n + m + 1`.
"""
function borel_pade(a::AbstractVector{T};
                    n::Integer, m::Integer, x::Real = 1,
                    regularize_poles::Bool = false,
                    ε::Union{Real,Nothing} = nothing,
                    return_error::Bool = false,
                    quad_kwargs...) where {T<:Real}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    if regularize_poles && x ≤ 0
        throw(ArgumentError("regularize_poles=true requires x > 0 (got $x)"))
    end
    Ba = borel_transform(a[1:n+m+1])
    Bp, Bq = pade(Ba, m, n)

    integrand = if regularize_poles
        ε_eff = ε === nothing ? 100 * eps(real(T)) : ε
        Bqr = poles_regularized(Polynomials.coeffs(Bq), ε_eff)
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
    borel_leroy_pade(a; n, m, b = -1//2, x = 1, return_error = false, quad_kwargs...)

Borel–Le Roy–Padé resummation:

1. Borel–Le Roy transform `a` with parameter `b` to obtain `B[k] = a[k] / Γ(k+1+b)`.
2. Replace `B(t)` by its Padé approximant `[n/m]`.
3. Laplace-transform with the Le Roy weight:

        ∫₀^∞ t^b (Bₚ(x t) / B_q(x t)) e^{-t} dt.

`b` may be any real number (including non-integer values); the default `-1//2`
matches the original code. `quad_kwargs` are forwarded to `QuadGK.quadgk`.
"""
function borel_leroy_pade(a::AbstractVector{T};
                          n::Integer, m::Integer,
                          b::Real = -1//2, x::Real = 1,
                          return_error::Bool = false,
                          quad_kwargs...) where {T<:Real}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_leroy_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_leroy_transform(a[1:n+m+1], b)
    Bp, Bq = pade(Ba, m, n)
    integrand = t -> t^b * Bp(x * t) / Bq(x * t) * exp(-t)
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
                              n::Integer, m::Integer, x::Real = 1,
                              sing::Real = 1,
                              return_error::Bool = false,
                              quad_kwargs...) where {T<:Real}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("conformal_borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_transform(a[1:n+m+1])
    Bw = conformal_reseries(Ba, sing, n + m)
    Bp, Bq = pade(Bw, m, n)
    integrand = function (t)
        w = conformal_map(x * t; a = sing)
        return Bp(w) / Bq(w) * exp(-t)
    end
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end
